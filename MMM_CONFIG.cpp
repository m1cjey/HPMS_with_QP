#include "stdafx.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"MMM_CONFIG.h"	//class CON定義

MMM_config::MMM_config()
{
	///////解析条件
	Hf_type=1;		//外部磁場　0:一様磁場 1:磁石
	Hf_H = 0.03 / (4 * PI * 1e-07);		//外部磁場Hの強さ
	isWLSM = 0;
	eForce = 0;             //電磁力(0:並進力,1:Kelvin力)
	force_t = 1;            //電磁力計算手法(0:応力のみ,1:体積力のみ)

	P_interval = 20;
}