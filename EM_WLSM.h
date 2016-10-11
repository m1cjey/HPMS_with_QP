#include "stdafx.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納

class WLSM_point3D//節点クラス
{
public:
	double r[DIMENSION];
	int material; //材質
	int boundary_condition; //境界条件 0=普通　1,2=固定境界 -1 自然境界条件(勾配がゼロ)
	int particleID;			//対応する粒子番号 存在しないときは-1を格納
	int remesh;				//リメッシュ領域に属するか、しないか
	int surface;			//物体の表面であるか、ないか。1=ON 0=OFF, 境界条件とは別
	int depth;				//界面からの深さ　界面を0とする
	int calced;				//calced=ONなら、ラプラシアンのための周辺粒子の寄与率は計算済みなので、そのまま使用すればよい。OFFなら再計算して求める
	double rp;				//比誘電率もしくは比透磁率

	int index;//格納されている格子の番号
	int N2;   //re2内に存在する周辺粒子数
	vector<int> NEI2;
	double L;				//平均粒子間距離
	double R;				//影響半径
	vector<double> fai;		//ラプラシアンを離散化するのに必要な周辺粒子の重み

	double E[3];			//電界または磁界ベクトル
	double potential;		
	double F[3];			//電磁力
	double Fs;				//電磁応力[N/m^2]
	double normal[3];		//外向き法線ベクトル

	int root;				//空気計算点に関して、水計算点の法線より生成されたならroot=FLUID,という風に、根っこの材質を格納する。解析領域なら空気。　空気計算点以外のrootは参照しない(-1を格納)
	
	int flag;				//適当な使い回し用変数
};

class INDEX//格子
{
public:
	vector<int> node;		//格納する節点番号
};