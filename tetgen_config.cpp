#include "stdafx.h"

tetgen_config::tetgen_config()
{
	//静電霧化装置の寸法
	radius_column=0.0005;		//円柱電極半径
	length_column=0.02;			//円柱電極長さ
	height_plate=0.0075;		//平板電極高さ
	length_plate=0.02;			//平板電極一辺長さ
	thickness_plate=0.0005;		//平板電極厚さ
	length_base=0.02;			//土台一辺長さ
	thickness_base=0.005;		//土台厚さ
	
	//メッシュの粗さ(境界の1要素の辺の長さをどれぐらいにするか)
	fine_air=0.01;				//空気領域境界の粗さ		0.01
	fine_plate_t=0.00025;		//平板電極厚み方向の粗さ	0.00025
	fine_plate_L=0.00050;		//平板電極平面方向の粗さ	0.0005
	fine_column_L=0.00015;		//円柱電極の長さ方向の粗さ	0.00015
	fine_base=0.0010;			//土台表面の粗さ			0.0010

	//水滴表面のメッシュ層の設定
	num_layer_out=1;	//流体外側メッシュ層数
	num_layer_in=0;		//流体内側メッシュ層数
	thick_layer=0.3;	//境界メッシュ1層の厚さ(leの何倍か)

	//長い流体要素を削除する閾値 
	del_length=3.0;		//leの何倍以上の辺を持つ要素を消すか2.0
}

tetgen_config::tetgen_config(mpsconfig &CON)
{
	//磁石寸法 cf. MPSTOFEM_MRE()<-MPS_TO_FEM3D.cpp
	magnet_height=CON.get_magnet_H();
	magnet_radius=CON.get_magnet_r();

	//メッシュの粗さ(境界の1要素の辺の長さをどれぐらいにするか)
	fine_air=0.01;				//空気領域境界の粗さ		0.01
	fine_plate_t=0.00025;		//平板電極厚み方向の粗さ	0.00025
	fine_plate_L=0.00050;		//平板電極平面方向の粗さ	0.0005
	fine_column_L=0.00015;		//円柱電極の長さ方向の粗さ	0.00015
	fine_base=0.0010;			//土台表面の粗さ			0.0010

	//水滴表面のメッシュ層の設定
	num_layer_out=1;	//流体外側メッシュ層数
	num_layer_in=0;		//流体内側メッシュ層数
	thick_layer=0.3;	//境界メッシュ1層の厚さ(leの何倍か)

	//長い流体要素を削除する閾値 
	del_length=2.0;		//leの何倍以上の辺を持つ要素を消すか

}
