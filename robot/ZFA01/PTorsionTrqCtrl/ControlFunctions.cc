//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2022/01/13
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

// 基本のインクルードファイル
#include <cmath>
#include "ControlFunctions.hh"
#include "ARCSprint.hh"
#include "ARCSassert.hh"
#include "ScreenParams.hh"
#include "InterfaceFunctions.hh"
#include "GraphPlot.hh"
#include "DataMemory.hh"

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "../robot/ZFA01/addon/ZFA01.hh"
#include "DisturbanceObsrv.hh"

using namespace ARCS;

//! @brief スレッド間通信用グローバル変数の無名名前空間
namespace {
	// スレッド間で共有したい変数をここに記述
	std::array<double, ConstParams::ACTUATOR_NUM> PositionRes = {0};	//!< [rad] 位置応答
	std::array<double, ConstParams::ACTUATOR_NUM> VelocityRes = {0};	//!< [rad/s] 速度応答
	std::array<double, ConstParams::ACTUATOR_NUM> TorqueRes = {0};		//!< [Nm] トルク応答
	std::array<double, ConstParams::ACTUATOR_NUM> CurrentRef = {0};		//!< [Nm]  電流指令
}

//! @brief 制御用周期実行関数1
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction1(double t, double Tact, double Tcmp){
	// 制御用定数設定
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[0]*1e-9;	// [s]	制御周期
	static constexpr Matrix<1,ZFA01::N_AX> Rg  = ZFA01::Rg;	// [-] 減速比
	static constexpr Matrix<1,ZFA01::N_AX> Kpt = {5, 4, 2};	// [-] ねじれトルク制御 比例ゲイン
	
	// 制御用変数宣言
	static Matrix<1,ZFA01::N_AX> thm;	// [rad] モータ側位置ベクトル
	static Matrix<1,ZFA01::N_AX> thl;	// [rad] 負荷側位置ベクトル
	static Matrix<1,ZFA01::N_AX> wm;	// [rad/s] モータ側速度ベクトル
	static Matrix<1,ZFA01::N_AX> wl;	// [rad/s] 負荷側速度ベクトル
	static Matrix<1,ZFA01::N_AX> taus;	// [Nm] ねじれトルクベクトル
	static Matrix<1,ZFA01::N_AX> iqref;	// [A] q軸電流指令ベクトル
	static Matrix<1,3> p1;				// [m] 作業空間位置ベクトル[px,py,pz]^T 1S軸
	static Matrix<1,3> p2;				// [m] 作業空間位置ベクトル[px,py,pz]^T 2L軸
	static Matrix<1,3> p3;				// [m] 作業空間位置ベクトル[px,py,pz]^T 3U軸
	static Matrix<1,3> pG;				// [m] 作業空間位置ベクトル[px,py,pz]^T 先端G(グラインド点)
	static Matrix<1,3> vG;				// [m/s] 作業空間速度ベクトル[vx,vy,vz]^T 先端G(グラインド点)
	static Matrix<1,ZFA01::N_AX> tsref;	// [Nm] ねじれトルク指令ベクトル
	
	static ZFA01 Robot;					// ZFA01ロボット制御クラス
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
		Initializing = true;				// 初期化中ランプ点灯
		Screen.InitOnlineSetVar();			// オンライン設定変数の初期値の設定
		Interface.ServoON();				// サーボON指令の送出
		Initializing = false;				// 初期化中ランプ消灯
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[0] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		Interface.GetPositionAndVelocity(PositionRes, VelocityRes);	// [rad,rad/s] 位置と速度応答の取得
		Interface.GetTorque(TorqueRes);		// [Nm] ねじれトルク応答の取得
		Screen.GetOnlineSetVar();			// オンライン設定変数の読み込み
		
		// センサ系データを縦ベクトルとして読み込んで整理
		thm.LoadArray(PositionRes);	// [rad] モータ側位置
		wm.LoadArray(VelocityRes);	// [rad/s] モータ側速度
		taus.LoadArray(TorqueRes);	// [Nm] ねじれトルク
		thl = thm % Rg;				// [rad] 負荷側位置(簡易的換算)
		wl = wm % Rg;				// [rad/s] 負荷側速度(簡易的換算)
		
		// 関節空間から作業空間への変換
		Robot.GetWorkspacePosition(thl, p1, p2, p3, pG);// [m] 作業空間位置の計算
		Robot.UpdateJacobi(thl);	// Jacobi行列の更新
		vG = Robot.Jaco(wl);		// [m/s] 作業空間速度の計算
		
		// 比例ねじれトルク制御
		tsref = Matrix<1,3>::zeros();		// [Nm] ねじれトルク指令
		iqref = (tsref - taus) & Kpt % Rg;	// [A]  電流指令
		// ↑やってみた感想：ゲインが全然上がらない。
		
		// 安全系と電流指令設定
		Robot.CheckSafety(thl, pG);		// 安全監視（範囲外になると緊急停止）
		Robot.CurrentLimiter(iqref);	// [A] q軸電流リミッタ
		iqref.StoreArray(CurrentRef);	// [A] 最終的な全軸の電流指令
		
		Interface.SetCurrent(CurrentRef);	// [A] 電流指令の出力
		Screen.SetVarIndicator(pG[1], pG[2], pG[3], 0, 0, 0, 0, 0, 0, 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);											// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, iqref[1], iqref[2], iqref[3], 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, thm[1], thm[2], thm[3], 0, 0, 0, 0, 0);		// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, wm[1], wm[2], wm[3], 0, 0, 0, 0, 0);			// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, taus[1], taus[2], taus[3], 0, 0, 0, 0, 0);		// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(4, vG[1], vG[2], vG[3], 0, 0, 0, 0, 0);			// グラフプロット4 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetWorkspace(p2, p3, pG);									// 作業空間グラフプロット
		Memory.SetData(Tact, t, tsref[1], taus[1], tsref[2], taus[2], tsref[3], taus[3], 0, 0, 0);	// CSVデータ保存変数 (周期, A列, B列, ..., J列)
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
		Interface.SetZeroCurrent();	// 電流指令を零に設定
		Interface.ServoOFF();		// サーボOFF信号の送出
	}
	return true;	// クロックオーバーライドフラグ(falseにすると次の周期時刻を待たずにスレッドが即刻動作する)
}

//! @brief 制御用周期実行関数2
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction2(double t, double Tact, double Tcmp){
	// 制御用定数宣言
	[[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[1]*1e-9;	// [s]	制御周期
	
	// 制御用変数宣言
	
	// 制御器等々の宣言
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[1] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	return true;	// クロックオーバーライドフラグ(falseにすると次の周期時刻を待たずにスレッドが即刻動作する)
}

//! @brief 制御用周期実行関数3
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction3(double t, double Tact, double Tcmp){
	// 制御用定数宣言
	[[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[2]*1e-9;	// [s]	制御周期
	
	// 制御用変数宣言
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[2] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	return true;	// クロックオーバーライドフラグ(falseにすると次の周期時刻を待たずにスレッドが即刻動作する)
}

//! @brief 制御用変数値を更新する関数
void ControlFunctions::UpdateControlValue(void){
	// ARCS画面パラメータに値を書き込む
	Screen.SetNetworkLink(NetworkLink);						// ネットワークリンクフラグを書き込む
	Screen.SetInitializing(Initializing);					// ロボット初期化フラグを書き込む
	Screen.SetCurrentAndPosition(CurrentRef, PositionRes);	// 電流指令と位置応答を書き込む
}

