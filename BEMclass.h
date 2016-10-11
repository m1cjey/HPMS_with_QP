#include "stdafx.h"	//主要なヘッダーファイルはまとめてこのなか。
//#include"define.h"	//#define 格納


/////BEM2D用のクラス

class point2D//節点クラス
{
public:
	double r[2];
	int boundary_condition; //境界条件 0=Diric 1=Neumn
	double potential;		//ポテンシャル
	double slop1;			//法線方向の勾配
	double slop2;			//法線方向の勾配
	double C;				//内角
	double L;				//一定要素の場合のみ使用。一定要素の場合、NODE[n].Lで節点nを中点とする一定要素の要素の長さが得られる
	int particle;			//対応する粒子番号格納 対応するものがない場合は－1をダミーとして格納
};

class element2D//要素クラス
{
public:
    int node[2];//要素を構築するnodeの番号 
	double r[2];		//中点の座標
	double L;//長さ
	double direct[2];		//法線ベクトル
	int boundary_condition; //境界条件 0=Diric 1=Neumn 一定要素の時に使用
	int map;//マッピング
	int material;//材質
};

class REGION//領域クラス
{
public:
	int start;	//startする計算点番号
	int end;
};

////////////


/////BEM3D用のクラス/////////
class BEMpoint3D//節点クラス
{
public:
	double r[3];
	int boundary_condition; //境界条件 0=Diric 1=Neumn
	double potential;		//ポテンシャル
	double slop1;			//法線方向の勾配
	double slop2;			//法線方向の勾配
	double C;				//内角
	int particle;			//対応する粒子番号格納 対応するものがない場合は－1をダミーとして格納
};

class BEMelement3D//要素クラス
{
public:
    int node[3];//要素を構築するnodeの番号 
	double r[3];		//中点の座標
	double S;		//面積
	double direct[3];		//法線ベクトル
	int boundary_condition; //境界条件 0=Diric 1=Neumn 一定要素の時に使用
	int map;//マッピング
	int material;//材質
};
