#include "stdafx.h"	
#include"define.h"	//#define 格納
#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include"PART.h"		//class PART定義
//#include"BEMclass.h"	//BEM2D関係のclass 定義
//#include"FEM3Dclass.h"	//FEM3D関係のclass 定義
#include"EM_WLSM.h"
//#include"FEM2Dclass.h"	//FEM2D class 定義
#include<omp.h>
#include<vector>
#include"function.h"
#include"MMM_CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい

//PartMM2 




void Cal_Mag_Density(double Ki[][3], double *p_vector, double src_size){
	double ri[3];
	double ri_norm;
	int weight;
	for(int i = 0; i < 2; i ++){
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){
				ri[0] = p_vector[0] - pow(-1.0, i+1) * src_size;
				ri[1] = p_vector[1] - pow(-1.0, j+1) * src_size;
				ri[2] = p_vector[2] - pow(-1.0, k+1) * src_size;
				ri_norm = sqrt(ri[0] * ri[0] + ri[1] * ri[1] + ri[2] * ri[2]);
				//if(ri[0] == 0)Ki[0][0] += pow(-1.0, i + j + k + 3) * (-PI / 2);
				//else 
				Ki[0][0] += pow(-1.0, i + j + k + 3) * (-atan2(ri[1] * ri[2], ri[0] * ri_norm));
				Ki[0][1] += pow(-1.0, i + j + k + 3) * log((ri_norm + ri[2]));
				Ki[0][2] += pow(-1.0, i + j + k + 3) * log((ri_norm + ri[1]));
				//if(ri[1] == 0)Ki[1][1] += pow(-1.0, i + j + k + 3) * (-PI / 2);
				//else 
				Ki[1][1] += pow(-1.0, i + j + k + 3) * (-atan2(ri[0] * ri[2], ri[1] * ri_norm));
				Ki[1][2] += pow(-1.0, i + j + k + 3) * log((ri_norm + ri[0]));
				//if(ri[2] == 0)Ki[2][2] += pow(-1.0, i + j + k + 3) * (-PI / 2);
				//else 
				Ki[2][2] += pow(-1.0, i + j + k + 3) * (-atan2(ri[0] * ri[1], ri[2] * ri_norm));

			}
		}
	}
	Ki[1][0] = Ki[0][1];
	Ki[2][1] = Ki[1][2];
	Ki[2][0] = Ki[0][2];
}

void Cal_Mag_Density(double Ki[][3], double *p_vector, double src_size_x, double src_size_y, double src_size_z){
	double ri[3];
	double ri_norm;
	int weight;
	for(int i = 0; i < 2; i ++){
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){
				ri[0] = p_vector[0] - pow(-1.0, i+1) * src_size_x;
				ri[1] = p_vector[1] - pow(-1.0, j+1) * src_size_y;
				ri[2] = p_vector[2] - pow(-1.0, k+1) * src_size_z;
				ri_norm = sqrt(ri[0] * ri[0] + ri[1] * ri[1] + ri[2] * ri[2]);

				//if(ri[0] == 0)Ki[0][0] += pow(-1.0, i + j + k + 3) * (-PI / 2);
				//else 
				Ki[0][0] += pow(-1.0, i + j + k + 3) * (-atan2(ri[1] * ri[2], ri[0] * ri_norm));
				Ki[0][1] += pow(-1.0, i + j + k + 3) * log((ri_norm + ri[2]));
				Ki[0][2] += pow(-1.0, i + j + k + 3) * log((ri_norm + ri[1]));
				//if(ri[1] == 0)Ki[1][1] += pow(-1.0, i + j + k + 3) * (-PI / 2);
				//else 
				Ki[1][1] += pow(-1.0, i + j + k + 3) * (-atan2(ri[0] * ri[2], ri[1] * ri_norm));
				Ki[1][2] += pow(-1.0, i + j + k + 3) * log((ri_norm + ri[0]));
				//if(ri[2] == 0)Ki[2][2] += pow(-1.0, i + j + k + 3) * (-PI / 2);
				//else 
				Ki[2][2] += pow(-1.0, i + j + k + 3) * (-atan2(ri[0] * ri[1], ri[2] * ri_norm));

			}
		}
	}
	Ki[1][0] = Ki[0][1];
	Ki[2][1] = Ki[1][2];
	Ki[2][0] = Ki[0][2];
}

//20140904:行列生成時に繰り返し数を1/4にするための計算
void Cal_Mag_Density_To_Matrix(double K1[][3], double K2[][3], double *p_vector, double src_size){
	double r1[3];
	double r2[3];
	double r1_norm;
	double r2_norm;
	for(int i = 0; i < 2; i ++){
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){
				r1[0] = p_vector[0] - pow(-1.0, i+1) * src_size;
				r1[1] = p_vector[1] - pow(-1.0, j+1) * src_size;
				r1[2] = p_vector[2] - pow(-1.0, k+1) * src_size;
				r1_norm = sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
				r2[0] = -p_vector[0] - pow(-1.0, i+1) * src_size;
				r2[1] = -p_vector[1] - pow(-1.0, j+1) * src_size;
				r2[2] = -p_vector[2] - pow(-1.0, k+1) * src_size;
				r2_norm = sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]);

				K1[0][0] += pow(-1.0, i + j + k + 3) * (-atan2(r1[1] * r1[2], r1[0] * r1_norm));
				K1[0][1] += pow(-1.0, i + j + k + 3) * log((r1_norm + r1[2]));
				K1[0][2] += pow(-1.0, i + j + k + 3) * log((r1_norm + r1[1]));
				K1[1][1] += pow(-1.0, i + j + k + 3) * (-atan2(r1[0] * r1[2], r1[1] * r1_norm));
				K1[1][2] += pow(-1.0, i + j + k + 3) * log((r1_norm + r1[0]));
				K1[2][2] += pow(-1.0, i + j + k + 3) * (-atan2(r1[0] * r1[1], r1[2] * r1_norm));

				K2[0][0] += pow(-1.0, i + j + k + 3) * (-atan2(r2[1] * r2[2], r2[0] * r2_norm));
				K2[0][1] += pow(-1.0, i + j + k + 3) * log((r2_norm + r2[2]));
				K2[0][2] += pow(-1.0, i + j + k + 3) * log((r2_norm + r2[1]));
				K2[1][1] += pow(-1.0, i + j + k + 3) * (-atan2(r2[0] * r2[2], r2[1] * r2_norm));
				K2[1][2] += pow(-1.0, i + j + k + 3) * log((r2_norm + r2[0]));
				K2[2][2] += pow(-1.0, i + j + k + 3) * (-atan2(r2[0] * r2[1], r2[2] * r2_norm));

			}
		}
	}
	K1[1][0] = K1[0][1];
	K1[2][1] = K1[1][2];
	K1[2][0] = K1[0][2];

	K2[1][0] = K2[0][1];
	K2[2][1] = K2[1][2];
	K2[2][0] = K2[0][2];
}

void GaussSeidel(vector<vector<double> >&A, int np, vector<double>&b, int dim){
	double tmp;
	double e = 0.0;
	vector<double> x(np);
	for(int k = 0; k < 10000; k++){
		e = 0.0;
		for(int i = 0; i < x.size(); i++){
			tmp = x[i];
			x[i] = b[i];
			for(int j = 0; j < x.size(); j++){
				x[i] -= (j != i ? A[i][j] * x[j] : 0.0);
			}
			x[i] /= A[i][i];
			e+=fabs(tmp - x[i]);
		}
		if(e <= 1e-6)break;
	}
	double temp[3];
	int part = 0;
	for(int i = 0; i < x.size(); i++){
		b[i] = x[i];
	}
}

void GaussSeidelv2(vector<vector<double> >&A, int np, vector<double>&b, int dim, vector<vector<int> > &index, vector<int>& number){
	double tmp;
	double e = 0.0;
	int diag;
	vector<double> x(np * dim);
	for(int k = 0; k < 10000; k++){
		e = 0.0;
		for(int i = 0; i < np; i++){
			for(int j = 0; j < dim; j++){
				tmp = x[i * dim + j];
				x[i * dim + j] = b[i * dim + j];
				for(int k = 0; k < number[i]; k++){
					if(i == index[i][k])diag = k;
					if(index[i][k] * dim + j != i * dim + j)x[i * dim + j]-=A[i * dim + j][k * dim + j] * x[index[i][k] * dim + j]; 
				}
				x[i * dim + j] /= A[i * dim + j][diag * dim +j];
				e+=fabs(tmp - x[i * dim + j]);
			}
		}
		cout<<"error="<<e<<endl;
		if(e <= 1e-5){
			cout<<"Solver converged.\n";
			break;
		}
	}
	double temp[3];
	int part = 0;
	for(int i = 0; i < x.size(); i++){
		b[i] = x[i];
	}
}

///H勾配計算関数ver.1
void H_gradient1(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,double **Hgrad,double **H)
{
	double le=CON.get_distancebp();//初期粒子間距離
	double r=CON.get_re()*le;
	int d=CON.get_dimension();

	double *HH=new double[fluid_number];//各粒子位置での磁場強度H格納
	for(int i=0;i<fluid_number;i++) HH[i]=H[A_X][i] + H[A_Y][i] + H[A_Z][i];


	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++) Hgrad[D][i]=0;//初期化

		double W=0;//粒子数密度　OUTを除いたりするのでPND[i]は微妙
		for(int k=0;k<PART[i].N;k++)
		{       
			int j=PART[i].NEI[k];
			if(PART[j].type==FLUID || PART[j].type==HYPERELAST)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);

				double w=kernel(r,dis);
				W+=w;

				Hgrad[A_X][i]+=(HH[j]-HH[i])*X*w/(dis*dis);
				Hgrad[A_Y][i]+=(HH[j]-HH[i])*Y*w/(dis*dis);
				Hgrad[A_Z][i]+=(HH[j]-HH[i])*Z*w/(dis*dis);
			}

		}
		for(int D=0;D<DIMENSION;D++) if(W!=0) Hgrad[D][i]= Hgrad[D][i]*d/W;
	}///////////////Pgrad[D][i]計算終了



	delete [] HH;
}

//ガウスの消去法 解は最終的にBのなかへ
void gauss_s(double *matrix,double *B,int N)
{
	for(int k=0;k<N;k++)
	{
		double akk=matrix[k*N+k];

		for(int i=0;i<N;i++)
		{
			if(i!=k)
			{
				double A=matrix[i*N+k]/akk;
				//for(int j=0;j<N;j++)
				for(int j=k;j<N;j++)
				{					
					matrix[i*N+j]-=A*matrix[k*N+j];
				}
				B[i]-=A*B[k];
			}
		}
	}
	for(int k=0;k<N;k++) B[k]/=matrix[k*N+k];

}

//任意の点における磁化を推定する
void get_M(double *cx, double *cy, double *cz, double x, double y, double z, double bx, double by, double bz, double *M){
	double dx = x - bx;
	double dy = y - by;
	double dz = z - bz;
	M[0] = cx[0] * dx + cx[1] * dy + cx[2] * dz + cx[3] * dx * dx + cx[4] * dy * dy + cx[5] * dz * dz + cx[6] * dx * dy + cx[7] * dy * dz + cx[8] * dx * dz;
	M[1] = cy[0] * dx + cy[1] * dy + cy[2] * dz + cy[3] * dx * dx + cy[4] * dy * dy + cy[5] * dz * dz + cy[6] * dx * dy + cy[7] * dy * dz + cy[8] * dx * dz;
	M[2] = cz[0] * dx + cz[1] * dy + cz[2] * dz + cz[3] * dx * dx + cz[4] * dy * dy + cz[5] * dz * dz + cz[6] * dx * dy + cz[7] * dy * dz + cz[8] * dx * dz;
}

void get_M2(double *cx, double *cy, double *cz, double x, double y, double z, double bx, double by, double bz, double *M){
	double dx = x - bx;
	double dy = y - by;
	double dz = z - bz;
	M[0] = cx[0] * dx + cx[1] * dy + cx[2] * dz + cx[3] * dx * dx + cx[4] * dy * dy + cx[5] * dz * dz + cx[6] * dx * dy + cx[7] * dy * dz + cx[8] * dx * dz;
	M[1] = cy[0] * dx + cy[1] * dy + cy[2] * dz + cy[3] * dx * dx + cy[4] * dy * dy + cy[5] * dz * dz + cy[6] * dx * dy + cy[7] * dy * dz + cy[8] * dx * dz;
	M[2] = cz[0] * dx + cz[1] * dy + cz[2] * dz + cz[3] * dx * dx + cz[4] * dy * dy + cz[5] * dz * dz + cz[6] * dx * dy + cz[7] * dy * dz + cz[8] * dx * dz;
}
//最小二乗法
void cal_WLSM_for_M(vector<vector<double>>&mesh, vector<vector<int>>&mesh_n, vector<double> &mesh_M, int mesh_id, double mesh_width, double *coeff_x, double *coeff_y, double *coeff_z){
	double matrix[100];
	double matrix_save[100];
	int N = 10;
	for(int k = 0; k < mesh_n[mesh_id].size(); k++){
		int j = mesh_n[mesh_id][k];
		double X = mesh[j][0] - mesh[mesh_id][0];
		double Y = mesh[j][1] - mesh[mesh_id][1];
		double Z = mesh[j][2] - mesh[mesh_id][2];
		double U = (mesh_M[j * 3] - mesh_M[mesh_id * 3]);
		double V = (mesh_M[j * 3 + 1] - mesh_M[mesh_id * 3 + 1]);
		double W = (mesh_M[j * 3 + 2] - mesh_M[mesh_id * 3 + 2]);

		double dis = sqrt(X * X + Y * Y + Z * Z);
		double w = 10;
		if(dis>mesh_width) w*=pow((mesh_width / dis),4);
		matrix[0]+=X*X*w;		//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*X*X*w;		//ΣXj^3wj
		matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
		matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
		matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
		matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj
		matrix[9]+=X*w;		//ΣXj^2Zjwj

		matrix[11]+=Y*Y*w;		
		matrix[12]+=Y*Z*w;		
		matrix[13]+=X*X*Y*w;		
		matrix[14]+=Y*Y*Y*w;		
		matrix[15]+=Y*Z*Z*w;
		matrix[16]+=X*Y*Y*w;
		matrix[17]+=Y*Y*Z*w;
		matrix[19]+=Y*w;

		matrix[22]+=Z*Z*w;			
		matrix[23]+=X*X*Z*w;
		matrix[24]+=Y*Y*Z*w;
		matrix[25]+=Y*Y*Y*w;
		matrix[29]+=Z*w;

		matrix[33]+=X*X*X*X*w;
		matrix[34]+=X*X*Y*Y*w;
		matrix[35]+=X*X*Z*Z*w;	
		matrix[36]+=X*X*X*Y*w;	
		matrix[37]+=X*X*Y*Z*w;	
		matrix[38]+=X*X*X*Z*w;	

		matrix[44]+=Y*Y*Y*Y*w;
		matrix[45]+=Y*Y*Z*Z*w;
		matrix[46]+=X*Y*Y*Y*w;
		matrix[47]+=Y*Y*Y*Z*w;
		matrix[48]+=X*Y*Y*Z*w;

		matrix[55]+=Z*Z*Z*Z*w;		//6行目
		matrix[56]+=X*Y*Z*Z*w;
		matrix[57]+=Y*Z*Z*Z*w;
		matrix[58]+=X*Z*Z*Z*w;

		matrix[99]+=w;

		coeff_x[0]+=U*X*w;		//a
		coeff_x[1]+=U*Y*w;		//b
		coeff_x[2]+=U*Z*w;		//c
		coeff_x[3]+=U*X*X*w;		//d
		coeff_x[4]+=U*Y*Y*w;		//e
		coeff_x[5]+=U*Z*Z*w;		//f
		coeff_x[6]+=U*X*Y*w;		//g
		coeff_x[7]+=U*Y*Z*w;		//h
		coeff_x[8]+=U*X*Z*w;		//i
		coeff_x[9]+=U*w;

		coeff_y[0]+=V*X*w;		//a
		coeff_y[1]+=V*Y*w;		//b
		coeff_y[2]+=V*Z*w;		//c
		coeff_y[3]+=V*X*X*w;		//d
		coeff_y[4]+=V*Y*Y*w;		//e
		coeff_y[5]+=V*Z*Z*w;		//f
		coeff_y[6]+=V*X*Y*w;		//g
		coeff_y[7]+=V*Y*Z*w;		//h
		coeff_y[8]+=V*X*Z*w;		//i
		coeff_y[9]+=V*w;	

		coeff_z[0]+=W*X*w;		//a
		coeff_z[1]+=W*Y*w;		//b
		coeff_z[2]+=W*Z*w;		//c
		coeff_z[3]+=W*X*X*w;		//d
		coeff_z[4]+=W*Y*Y*w;		//e
		coeff_z[5]+=W*Z*Z*w;		//f
		coeff_z[6]+=W*X*Y*w;		//g
		coeff_z[7]+=W*Y*Z*w;		//h
		coeff_z[8]+=W*X*Z*w;		//i
		coeff_z[9]+=W*w;
	}
	matrix[10]=matrix[1];
	matrix[18]=matrix[7];

	matrix[20]=matrix[2];
	matrix[21]=matrix[12];
	matrix[24]=matrix[16];
	matrix[26]=matrix[7];
	matrix[27]=matrix[15];
	matrix[28]=matrix[5];

	for(int k=0;k<=2;k++) matrix[30+k]=matrix[3+10*k];//30～32要素
	matrix[39]=matrix[0];

	for(int k=0;k<=3;k++) matrix[40+k]=matrix[4+10*k];//40～43要素
	matrix[49]=matrix[11];

	for(int k=0;k<=4;k++) matrix[50+k]=matrix[5+10*k];//50～54要素
	matrix[59]=matrix[22];

	for(int k=0;k<=5;k++) matrix[60+k]=matrix[6+10*k];//60～65要素
	matrix[66]=matrix[34];
	matrix[67]=matrix[48];
	matrix[68]=matrix[37];
	matrix[69]=matrix[1];

	for(int k=0;k<=6;k++) matrix[70+k]=matrix[7+10*k];//70～76要素
	matrix[77]=matrix[54];
	matrix[78]=matrix[56];
	matrix[79]=matrix[12];

	for(int k=0;k<=7;k++) matrix[80+k]=matrix[8+10*k];//80～87要素
	matrix[88]=matrix[35];
	matrix[89]=matrix[20];

	for(int k=0;k<=8;k++) matrix[90+k]=matrix[9+10*k];//90～98要素

	matrix[99]+=20;//自身
	coeff_x[9]+=mesh_M[mesh_id * 3] * 10;
	coeff_y[9]+=mesh_M[mesh_id * 3 + 1] * 10;
	coeff_z[9]+=mesh_M[mesh_id * 3 + 2] * 10;
	for(int L=0;L<100;L++) matrix_save[L]=matrix[L];//行列の値を保存

	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;
	int Yflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Zflag=OFF;
	for(int k=0;k<9;k++)
	{
		if(coeff_x[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(coeff_y[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(coeff_z[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}
	if(Xflag==ON)
	{
		gauss(matrix,coeff_x,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_save[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) coeff_x[k]=0;} //flagがOFFならどのみちdudxはゼロ。

	if(Yflag==ON) 
	{
		gauss(matrix,coeff_y,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_save[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) coeff_y[k]=0;}	//flagがOFFならどのみちdvdyはゼロ。

	if(Zflag==ON) gauss(matrix,coeff_z,N);
	else {for(int k=0;k<N;k++) coeff_z[k]=0;}	//flagがOFFならどのみちdwdzはゼロ。
}
void cal_WLSM_for_grad_2D2(vector<mpselastic>&PART, int part_id, double R, double *M[3], double *H[3], double *grad){
	//この関数は，電磁力計算を行うためのもの. 解く式は∇(M・H)
	double matrix[36] = {.0};
	double coeff[6] = {.0};
	int N = 6;
	for(int k = 0; k < PART[part_id].N3; k++){
		int j = PART[part_id].NEI3[k];
		if(PART[j].type != FLUID && PART[j].type!=HYPERELAST) continue; //流体じゃないと磁束密度と磁化を持たない．加えてvectorも確保していない．
		double X = PART[j].r[0] - PART[part_id].r[0];
		double Y = PART[j].r[1] - PART[part_id].r[1];
		double Z = M[0][j] * H[0][j] + M[1][j] * H[1][j] + M[2][j] * H[2][j]- (M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id] + M[2][part_id] * H[2][part_id]);


		double dis = sqrt(X * X + Y * Y);
		double w = 1;
		//if(dis < 1.2 * mesh_width) w = 1.;
		w = (1. - 6. * pow((dis / R),2) +8. * pow((dis / R),3) - 3. * pow((dis /R),4));
		//w = (exp(-pow(dis / 4,2)) - exp(-pow(R / 4,2))) / (1 - exp(-pow(R / 4,2)));
		//else w*=pow((mesh_width / dis), 4.0);
		matrix[0]+=X*X*X*X*w;			
		matrix[1]+=X*X*Y*Y*w;	
		matrix[2]+=X*X*X*Y*w;
		matrix[3]+=X*X*X*w;
		matrix[4]+=X*X*Y*w;
		matrix[5]+=X*X*w;

		matrix[7]+=Y*Y*Y*Y*w;
		matrix[8]+=X*Y*Y*Y*w;
		matrix[9]+=X*Y*Y*w;
		matrix[10]+=Y*Y*Y*w;
		matrix[11]+=Y*Y*w;

		matrix[17]+=X*Y*w;

		matrix[23]+=X*w;

		matrix[29]+=Y*w;

		matrix[35]+=w;

		coeff[0]+=X*X*Z*w;
		coeff[1]+=Y*Y*Z*w;
		coeff[2]+=X*Y*Z*w;
		coeff[3]+=X*Z*w;
		coeff[4]+=Y*Z*w;
		coeff[5]+=Z*w;

	}
	matrix[6]=matrix[1];

	matrix[12]=matrix[2];
	matrix[13]=matrix[8];
	matrix[14]=matrix[1];
	matrix[15]=matrix[4];
	matrix[16]=matrix[9];

	matrix[18]=matrix[3];
	matrix[19]=matrix[9];
	matrix[20]=matrix[15];
	matrix[21]=matrix[5];
	matrix[22]=matrix[17];

	matrix[24]=matrix[4];
	matrix[25]=matrix[10];
	matrix[26]=matrix[16];
	matrix[27]=matrix[22];
	matrix[28]=matrix[11];

	matrix[30]=matrix[5];
	matrix[31]=matrix[11];
	matrix[32]=matrix[17];
	matrix[33]=matrix[23];
	matrix[34]=matrix[29];

	matrix[35]+=1;//自分自身
	coeff[5]+=M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id] + M[2][part_id] * H[2][part_id];


	//行列をガウスの消去法で解く　解はBに格納される
	bool flag = OFF;
	for(int k=0;k<6;k++)
	{
		if(coeff[k]!=0) flag=ON;
	}
	if(flag==ON)gauss(matrix,coeff,N);
	//各係数を用いて勾配の計算
	grad[0] = coeff[3];
	grad[1] = coeff[4];
}
void cal_WLSM_for_grad2(vector<mpselastic>&PART, int part_id, double R, double *M[3], double *H[3], double *grad){
	double matrix[100]={.0};
	double coeff[10]={.0};
	int N = 10;
	for(int k = 0; k < PART[part_id].N3; k++){
		int j = PART[part_id].NEI3[k];
		if(PART[j].type != FLUID&&PART[j].type != HYPERELAST)continue; //流体じゃないと磁束密度と磁化を持たない．加えてvectorも確保していない．
		double X = PART[j].r[0] - PART[part_id].r[0];
		double Y = PART[j].r[1] - PART[part_id].r[1];
		double Z = PART[j].r[2] - PART[part_id].r[2];
		double K = M[0][j] * H[0][j] + M[1][j] * H[1][j] + M[2][j] * H[2][j] - (M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id] + M[2][part_id] * H[2][part_id]);

		double dis = sqrt(X * X + Y * Y + Z * Z);
		double w = (1. - 6. * pow((dis / R),2) +8. * pow((dis / R),3) - 3. * pow((dis /R),4));
		matrix[0]+=X*X*w;		//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*X*X*w;		//ΣXj^3wj
		matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
		matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
		matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
		matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj
		matrix[9]+=X*w;		//ΣXj^2Zjwj

		matrix[11]+=Y*Y*w;		
		matrix[12]+=Y*Z*w;		
		matrix[13]+=X*X*Y*w;		
		matrix[14]+=Y*Y*Y*w;		
		matrix[15]+=Y*Z*Z*w;
		matrix[16]+=X*Y*Y*w;
		matrix[17]+=Y*Y*Z*w;
		matrix[19]+=Y*w;

		matrix[22]+=Z*Z*w;			
		matrix[23]+=X*X*Z*w;
		matrix[24]+=Y*Y*Z*w;
		matrix[25]+=Y*Y*Y*w;
		matrix[29]+=Z*w;

		matrix[33]+=X*X*X*X*w;
		matrix[34]+=X*X*Y*Y*w;
		matrix[35]+=X*X*Z*Z*w;	
		matrix[36]+=X*X*X*Y*w;	
		matrix[37]+=X*X*Y*Z*w;	
		matrix[38]+=X*X*X*Z*w;	

		matrix[44]+=Y*Y*Y*Y*w;
		matrix[45]+=Y*Y*Z*Z*w;
		matrix[46]+=X*Y*Y*Y*w;
		matrix[47]+=Y*Y*Y*Z*w;
		matrix[48]+=X*Y*Y*Z*w;

		matrix[55]+=Z*Z*Z*Z*w;	//6行目
		matrix[56]+=X*Y*Z*Z*w;
		matrix[57]+=Y*Z*Z*Z*w;
		matrix[58]+=X*Z*Z*Z*w;

		matrix[99]+=w;

		coeff[0]+=K*X*w;		//a
		coeff[1]+=K*Y*w;		//b
		coeff[2]+=K*Z*w;		//c
		coeff[3]+=K*X*X*w;		//d
		coeff[4]+=K*Y*Y*w;		//e
		coeff[5]+=K*Z*Z*w;		//f
		coeff[6]+=K*X*Y*w;		//g
		coeff[7]+=K*Y*Z*w;		//h
		coeff[8]+=K*X*Z*w;		//i
		coeff[9]+=K*w;
	}
	matrix[10]=matrix[1];
	matrix[18]=matrix[7];

	matrix[20]=matrix[2];
	matrix[21]=matrix[12];
	matrix[24]=matrix[16];
	matrix[26]=matrix[7];
	matrix[27]=matrix[15];
	matrix[28]=matrix[5];

	for(int k=0;k<=2;k++) matrix[30+k]=matrix[3+10*k];//30～32要素
	matrix[39]=matrix[0];

	for(int k=0;k<=3;k++) matrix[40+k]=matrix[4+10*k];//40～43要素
	matrix[49]=matrix[11];

	for(int k=0;k<=4;k++) matrix[50+k]=matrix[5+10*k];//50～54要素
	matrix[59]=matrix[22];

	for(int k=0;k<=5;k++) matrix[60+k]=matrix[6+10*k];//60～65要素
	matrix[66]=matrix[34];
	matrix[67]=matrix[48];
	matrix[68]=matrix[37];
	matrix[69]=matrix[1];

	for(int k=0;k<=6;k++) matrix[70+k]=matrix[7+10*k];//70～76要素
	matrix[77]=matrix[54];
	matrix[78]=matrix[56];
	matrix[79]=matrix[12];

	for(int k=0;k<=7;k++) matrix[80+k]=matrix[8+10*k];//80～87要素
	matrix[88]=matrix[35];
	matrix[89]=matrix[20];

	for(int k=0;k<=8;k++) matrix[90+k]=matrix[9+10*k];//90～98要素

	matrix[99]+=1;//自身
	coeff[9]+=M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id] + M[2][part_id] * H[2][part_id];

	//行列をガウスの消去法で解く　解はBに格納される
	int flag=OFF;
	for(int k=0;k<N;k++)
	{
		if(coeff[k]!=0) flag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}
	if(flag==ON)gauss(matrix,coeff,N);
	else {for(int k=0;k<N;k++) coeff[k]=0;} //flagがOFFならどのみちdudxはゼロ。
	grad[0] = coeff[0];
	grad[1] = coeff[1];
	grad[2] = coeff[2];
}
//最小二乗法(自身を通らない)
void cal_WLSM_for_M2(vector<vector<double>>&mesh, vector<vector<int>>&mesh_n, vector<double> &mesh_M, int mesh_id, double mesh_width, double R, double *coeff_x, double *coeff_y, double *coeff_z){
	double matrix[81] = {.0};
	double matrix_save[81] = {.0};
	int N = 9;
	for(int k = 0; k < mesh_n[mesh_id].size(); k++){
		int j = mesh_n[mesh_id][k];
		double X = mesh[j][0] - mesh[mesh_id][0];
		double Y = mesh[j][1] - mesh[mesh_id][1];
		double Z = mesh[j][2] - mesh[mesh_id][2];
		double U = (mesh_M[j * 3] - mesh_M[mesh_id * 3]);
		double V = (mesh_M[j * 3 + 1] - mesh_M[mesh_id * 3 + 1]);
		double W = (mesh_M[j * 3 + 2] - mesh_M[mesh_id * 3 + 2]);

		double dis = sqrt(X * X + Y * Y + Z * Z);
		double w = 1;
		//if(dis < 1.2 * mesh_width) w = 1.;
		//w = (1. - 6. * pow((dis / R),2) +8. * pow((dis / R),3) - 3. * pow((dis /R),4));
		//w = (exp(-pow(dis / 4,2)) - exp(-pow(R / 4,2))) / (1 - exp(-pow(R / 4,2)));
		//else w*=pow((mesh_width / dis), 4.0);
		matrix[0]+=X*X*w;		//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*X*X*w;		//ΣXj^3wj
		matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
		matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
		matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
		matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj

		matrix[10]+=Y*Y*w;		
		matrix[11]+=Y*Z*w;		
		matrix[12]+=X*X*Y*w;		
		matrix[13]+=Y*Y*Y*w;		
		matrix[14]+=Y*Z*Z*w;
		matrix[15]+=X*Y*Y*w;
		matrix[16]+=Y*Y*Z*w;

		matrix[20]+=Z*Z*w;			
		matrix[23]+=Z*Z*Z*w;		

		matrix[30]+=X*X*X*X*w;
		matrix[31]+=X*X*Y*Y*w;
		matrix[32]+=X*X*Z*Z*w;	
		matrix[33]+=X*X*X*Y*w;	
		matrix[34]+=X*X*Y*Z*w;	
		matrix[35]+=X*X*X*Z*w;	

		matrix[40]+=Y*Y*Y*Y*w;
		matrix[41]+=Y*Y*Z*Z*w;
		matrix[42]+=X*Y*Y*Y*w;
		matrix[43]+=Y*Y*Y*Z*w;
		matrix[44]+=X*Y*Y*Z*w;

		matrix[50]+=Z*Z*Z*Z*w;	//6行目
		matrix[51]+=X*Y*Z*Z*w;
		matrix[52]+=Y*Z*Z*Z*w;
		matrix[53]+=X*Z*Z*Z*w;

		//matrix[99]+=w;

		coeff_x[0]+=U*X*w;		//a
		coeff_x[1]+=U*Y*w;		//b
		coeff_x[2]+=U*Z*w;		//c
		coeff_x[3]+=U*X*X*w;		//d
		coeff_x[4]+=U*Y*Y*w;		//e
		coeff_x[5]+=U*Z*Z*w;		//f
		coeff_x[6]+=U*X*Y*w;		//g
		coeff_x[7]+=U*Y*Z*w;		//h
		coeff_x[8]+=U*X*Z*w;		//i

		coeff_y[0]+=V*X*w;		//a
		coeff_y[1]+=V*Y*w;		//b
		coeff_y[2]+=V*Z*w;		//c
		coeff_y[3]+=V*X*X*w;		//d
		coeff_y[4]+=V*Y*Y*w;		//e
		coeff_y[5]+=V*Z*Z*w;		//f
		coeff_y[6]+=V*X*Y*w;		//g
		coeff_y[7]+=V*Y*Z*w;		//h
		coeff_y[8]+=V*X*Z*w;		//i

		coeff_z[0]+=W*X*w;		//a
		coeff_z[1]+=W*Y*w;		//b
		coeff_z[2]+=W*Z*w;		//c
		coeff_z[3]+=W*X*X*w;		//d
		coeff_z[4]+=W*Y*Y*w;		//e
		coeff_z[5]+=W*Z*Z*w;		//f
		coeff_z[6]+=W*X*Y*w;		//g
		coeff_z[7]+=W*Y*Z*w;		//h
		coeff_z[8]+=W*X*Z*w;		//i
	}
	matrix[9]=matrix[1];		//ΣXjYjwj
	matrix[17]=matrix[7];

	matrix[18]=matrix[2];
	matrix[19]=matrix[11];
	matrix[21]=matrix[8];
	matrix[22]=matrix[16];
	matrix[24]=matrix[7];
	matrix[25]=matrix[14];
	matrix[26]=matrix[5];

	matrix[27]=matrix[3];
	matrix[28]=matrix[12];
	matrix[29]=matrix[21];

	matrix[36]=matrix[4];
	matrix[37]=matrix[13];
	matrix[38]=matrix[22];
	matrix[39]=matrix[31];

	matrix[45]=matrix[5];
	matrix[46]=matrix[14];
	matrix[47]=matrix[23];
	matrix[48]=matrix[32];
	matrix[49]=matrix[41];

	matrix[54]=matrix[6];
	matrix[55]=matrix[15];
	matrix[56]=matrix[24];
	matrix[57]=matrix[33];
	matrix[58]=matrix[42];
	matrix[59]=matrix[51];
	matrix[60]=matrix[31];
	matrix[61]=matrix[44];
	matrix[62]=matrix[34];

	matrix[63]=matrix[7];
	matrix[64]=matrix[16];
	matrix[65]=matrix[25];
	matrix[66]=matrix[34];
	matrix[67]=matrix[43];
	matrix[68]=matrix[52];
	matrix[69]=matrix[61];
	matrix[70]=matrix[41];
	matrix[71]=matrix[51];

	matrix[72]=matrix[8];
	matrix[73]=matrix[17];
	matrix[74]=matrix[26];
	matrix[75]=matrix[35];
	matrix[76]=matrix[44];
	matrix[77]=matrix[53];
	matrix[78]=matrix[62];
	matrix[79]=matrix[71];
	matrix[80]=matrix[32];
	for(int L=0;L<81;L++) matrix_save[L]=matrix[L];//行列の値を保存

	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;
	int Yflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Zflag=OFF;
	for(int k=0;k<9;k++)
	{
		if(coeff_x[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(coeff_y[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(coeff_z[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}
	if(Xflag==ON)
	{
		gauss(matrix,coeff_x,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_save[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) coeff_x[k]=0;} //flagがOFFならどのみちdudxはゼロ。

	if(Yflag==ON) 
	{
		gauss(matrix,coeff_y,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_save[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) coeff_y[k]=0;}	//flagがOFFならどのみちdvdyはゼロ。

	if(Zflag==ON) gauss(matrix,coeff_z,N);
	else {for(int k=0;k<N;k++) coeff_z[k]=0;}	//flagがOFFならどのみちdwdzはゼロ。
}
void cal_WLSM_for_grad(vector<mpselastic>&PART, int part_id, double R, double *M[3], double *H[3], double *grad){
	//この関数は，電磁力計算を行うためのもの. 解く式は∇(M・H)
	double matrix[81] = {.0};
	double coeff[9] = {.0};
	int N = 9;
	for(int k = 0; k < PART[part_id].N; k++)
	{
		int j = PART[part_id].NEI[k];
		if(PART[j].type != FLUID&&PART[j].type != HYPERELAST)continue; //流体じゃないと磁束密度と磁化を持たない．加えてvectorも確保していない．
		double X = PART[j].r[0] - PART[part_id].r[0];
		double Y = PART[j].r[1] - PART[part_id].r[1];
		double Z = PART[j].r[2] - PART[part_id].r[2];
		double K = M[0][j] * H[0][j] + M[1][j] * H[1][j] + M[2][j] * H[2][j] - (M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id] + M[2][part_id] * H[2][part_id]);


		double dis = sqrt(X * X + Y * Y + Z * Z);
		double w = 1;
		//if(dis < 1.2 * mesh_width) w = 1.;
		w = (1. - 6. * pow((dis / R),2) +8. * pow((dis / R),3) - 3. * pow((dis /R),4));
		//w=1.;
		//w = (exp(-pow(dis / 4,2)) - exp(-pow(R / 4,2))) / (1 - exp(-pow(R / 4,2)));
		//else w*=pow((mesh_width / dis), 4.0);
		matrix[0]+=X*X*w;		//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*X*X*w;		//ΣXj^3wj
		matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
		matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
		matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
		matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj

		matrix[10]+=Y*Y*w;		
		matrix[11]+=Y*Z*w;		
		matrix[12]+=X*X*Y*w;		
		matrix[13]+=Y*Y*Y*w;		
		matrix[14]+=Y*Z*Z*w;
		matrix[15]+=X*Y*Y*w;
		matrix[16]+=Y*Y*Z*w;

		matrix[20]+=Z*Z*w;			
		matrix[23]+=Z*Z*Z*w;		

		matrix[30]+=X*X*X*X*w;
		matrix[31]+=X*X*Y*Y*w;
		matrix[32]+=X*X*Z*Z*w;	
		matrix[33]+=X*X*X*Y*w;	
		matrix[34]+=X*X*Y*Z*w;	
		matrix[35]+=X*X*X*Z*w;	

		matrix[40]+=Y*Y*Y*Y*w;
		matrix[41]+=Y*Y*Z*Z*w;
		matrix[42]+=X*Y*Y*Y*w;
		matrix[43]+=Y*Y*Y*Z*w;
		matrix[44]+=X*Y*Y*Z*w;

		matrix[50]+=Z*Z*Z*Z*w;	//6行目
		matrix[51]+=X*Y*Z*Z*w;
		matrix[52]+=Y*Z*Z*Z*w;
		matrix[53]+=X*Z*Z*Z*w;

		//解行列
		coeff[0]+=K*X*w;		//a
		coeff[1]+=K*Y*w;		//b
		coeff[2]+=K*Z*w;		//c
		coeff[3]+=K*X*X*w;		//d
		coeff[4]+=K*Y*Y*w;		//e
		coeff[5]+=K*Z*Z*w;		//f
		coeff[6]+=K*X*Y*w;		//g
		coeff[7]+=K*Y*Z*w;		//h
		coeff[8]+=K*X*Z*w;		//i
	}
	//cout<<"coeff"<<coeff[0]<<","<<coeff[1]<<","<<coeff[2]<<","<<coeff[3]<<","<<coeff[4]<<","<<coeff[5]<<","<<coeff[6]<<","<<coeff[7]<<","<<coeff[8]<<endl;
	matrix[9]=matrix[1];		//ΣXjYjwj
	matrix[17]=matrix[7];

	matrix[18]=matrix[2];
	matrix[19]=matrix[11];
	matrix[21]=matrix[8];
	matrix[22]=matrix[16];
	matrix[24]=matrix[7];
	matrix[25]=matrix[14];
	matrix[26]=matrix[5];

	matrix[27]=matrix[3];
	matrix[28]=matrix[12];
	matrix[29]=matrix[21];

	matrix[36]=matrix[4];
	matrix[37]=matrix[13];
	matrix[38]=matrix[22];
	matrix[39]=matrix[31];

	matrix[45]=matrix[5];
	matrix[46]=matrix[14];
	matrix[47]=matrix[23];
	matrix[48]=matrix[32];
	matrix[49]=matrix[41];

	matrix[54]=matrix[6];
	matrix[55]=matrix[15];
	matrix[56]=matrix[24];
	matrix[57]=matrix[33];
	matrix[58]=matrix[42];
	matrix[59]=matrix[51];
	matrix[60]=matrix[31];
	matrix[61]=matrix[44];
	matrix[62]=matrix[34];

	matrix[63]=matrix[7];
	matrix[64]=matrix[16];
	matrix[65]=matrix[25];
	matrix[66]=matrix[34];
	matrix[67]=matrix[43];
	matrix[68]=matrix[52];
	matrix[69]=matrix[61];
	matrix[70]=matrix[41];
	matrix[71]=matrix[51];

	matrix[72]=matrix[8];
	matrix[73]=matrix[17];
	matrix[74]=matrix[26];
	matrix[75]=matrix[35];
	matrix[76]=matrix[44];
	matrix[77]=matrix[53];
	matrix[78]=matrix[62];
	matrix[79]=matrix[71];
	matrix[80]=matrix[32];

	//行列をガウスの消去法で解く　解はBに格納される
	bool flag = OFF;
	for(int k=0;k<9;k++)
	{
		if(coeff[k]!=0) flag=ON;
	}
	if(flag==ON)	gauss_s(matrix,coeff,N);
	//各係数を用いて勾配の計算
	grad[0] = coeff[0];
	grad[1] = coeff[1];
	grad[2] = coeff[2];
	//cout<<"grad="<<grad[0]<<","<<grad[1]<<","<<grad[2]<<endl;

}
void cal_WLSM_for_grad_2D(vector<mpselastic>&PART, int part_id, double R, double *M[3], double *H[3], double *grad){
	//この関数は，電磁力計算を行うためのもの. 解く式は∇(M・H)
	double matrix[25] = {.0};
	double coeff[5] = {.0};
	int N = 5;
	for(int k = 0; k < PART[part_id].N3; k++){
		int j = PART[part_id].NEI3[k];
		if(PART[j].type != FLUID && PART[j].type != HYPERELAST)continue; //流体じゃないと磁束密度と磁化を持たない．加えてvectorも確保していない．
		double X = PART[j].r[0] - PART[part_id].r[0];
		double Y = PART[j].r[1] - PART[part_id].r[1];
		double Z = M[0][j] * H[0][j] + M[1][j] * H[1][j] - (M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id]);


		double dis = sqrt(X * X + Y * Y);
		double w = 1;
		//if(dis < 1.2 * mesh_width) w = 1.;
		w = (1. - 6. * pow((dis / R),2) +8. * pow((dis / R),3) - 3. * pow((dis /R),4));
		//w = (exp(-pow(dis / 4,2)) - exp(-pow(R / 4,2))) / (1 - exp(-pow(R / 4,2)));
		//else w*=pow((mesh_width / dis), 4.0);
		matrix[0]+=X*X*w;			//ΣXjwj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*X*X*w;			//ΣXj^3wj
		matrix[3]+=X*X*Y*w;			//ΣXj^2Yjwj
		matrix[4]+=X*Y*Y*w;			//ΣXjYj^2wj

		matrix[6]+=Y*Y*w;			//ΣYj^2wj
		matrix[9]+=Y*Y*Y*w;			//ΣYj^3wj

		matrix[12]+=X*X*X*X*w;			//ΣXj^4wj
		matrix[13]+=X*X*X*Y*w;			//ΣXj^3Yjwj
		matrix[14]+=X*X*Y*Y*w;			//ΣXj^2Yj^2wj

		matrix[19]+=X*Y*Y*Y*w;			//ΣXjYj^3wj

		matrix[24]+=Y*Y*Y*Y*w;

		coeff[0]+=Z*X*w;//ΣdPjXjwj
		coeff[1]+=Z*Y*w;//ΣdPjYjwj
		coeff[2]+=Z*X*X*w;//ΣdPjXj^2wj
		coeff[3]+=Z*X*Y*w;//ΣdPjXjYjwj
		coeff[4]+=Z*Y*Y*w;//ΣdPjYj^2wj

	}
	matrix[5]=matrix[1];		//ΣXjYjwj
	matrix[7]=matrix[3];		//ΣXj^2Yjwj
	matrix[8]=matrix[4];		//ΣXjYj^2wj
	matrix[10]=matrix[2];		//ΣXj^3Yjwj
	matrix[11]=matrix[3];		//ΣXj^2Yjwj
	matrix[15]=matrix[3];		//ΣXj^2Yjwj
	matrix[16]=matrix[4];		//ΣXjYj^2wj
	matrix[17]=matrix[13];		//ΣXj^3Yjwj
	matrix[18]=matrix[14];		//ΣXj^2Yj^2wj
	matrix[20]=matrix[4];		//ΣXjYj^2wj
	matrix[21]=matrix[9];		//ΣYj^3wj
	matrix[22]=matrix[14];		//ΣXj^2Yj^2wj
	matrix[23]=matrix[19];		//ΣXjYj^3wj


	//行列をガウスの消去法で解く　解はBに格納される
	bool flag = OFF;
	for(int k=0;k<5;k++)
	{
		if(coeff[k]!=0) flag=ON;
	}
	if(flag==ON)gauss(matrix,coeff,N);
	//各係数を用いて勾配の計算
	grad[0] = coeff[0];
	grad[1] = coeff[1];
}

void calc_MPS(mpsconfig &CON, vector<hyperelastic>HYPER,vector<mpselastic>PART,int part_id,double *M[3], double *H[3], double *grad)
{
	double IP=0.;
	double r=CON.get_h_dis();
	int h_num=HYPER.size();
	int d=3;
	int i=part_id;

	double pnd0=HYPER[i].pnd0;
	double qij[DIMENSION]={0,0,0};
	double dis=0.;
	double d_dis=0.;
	double w=0.;
	double p_grad[DIMENSION]={0,0,0};

	for(int j=0;j<h_num;j++)
	{
		if(j!=i)
		{
			qij[A_X]=PART[j].r[A_X]-PART[i].r[A_X];
			qij[A_Y]=PART[j].r[A_Y]-PART[i].r[A_Y];
			qij[A_Z]=PART[j].r[A_Z]-PART[i].r[A_Z];

			dis=sqrt((PART[j].r[A_X]-PART[i].r[A_X])*(PART[j].r[A_X]-PART[i].r[A_X])+(PART[j].r[A_Y]-PART[i].r[A_Y])*(PART[j].r[A_Y]-PART[i].r[A_Y])+(PART[j].r[A_Z]-PART[i].r[A_Z])*(PART[j].r[A_Z]-PART[i].r[A_Z]));
			d_dis=((PART[j].r[A_X]-PART[i].r[A_X])*(PART[j].r[A_X]-PART[i].r[A_X])+(PART[j].r[A_Y]-PART[i].r[A_Y])*(PART[j].r[A_Y]-PART[i].r[A_Y])+(PART[j].r[A_Z]-PART[i].r[A_Z])*(PART[j].r[A_Z]-PART[i].r[A_Z]));

			if(dis<r)
			{
				IP=M[0][j] * H[0][j] + M[1][j] * H[1][j] - (M[0][part_id] * H[0][part_id] + M[1][part_id] * H[1][part_id]);
				w=kernel4(r,dis);
				p_grad[A_X]+=d/pnd0*IP*qij[A_X]*w/d_dis;
				p_grad[A_Y]+=d/pnd0*IP*qij[A_Y]*w/d_dis;
				p_grad[A_Z]+=d/pnd0*IP*qij[A_Z]*w/d_dis;
			}
		}
	}
	grad[A_X]=p_grad[A_X];
	grad[A_Y]=p_grad[A_Y];
	grad[A_Z]=p_grad[A_Z];

}


//発散
double diver(mpsconfig &CON,vector<mpselastic> &PART,int i,double n0, double **M)
{
	double W=0;										//粒子数密度
	double R=CON.get_distancebp()*CON.get_re();	//影響半径
	double div=0;									//発散の値

	for(int k=0;k<PART[i].N;k++)
	{    
		int j=PART[i].NEI[k]; 
		double X=PART[j].r[A_X]-PART[i].r[A_X];
		double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
		double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
		double dis=sqrt(X*X+Y*Y+Z*Z);

		double w=kernel(R,dis);

		div+=(M[0][j] - M[0][i])*X*w/(dis*dis);
		div+=(M[1][j] - M[1][i])*Y*w/(dis*dis);
		div+=(M[2][j] - M[2][i])*Z*w/(dis*dis);
		W+=w;
	}
	if(W!=0)
	{
		div*=CON.get_dimension()/W;
	}
	return div;
}

///磁気モーメント法計算関数
///磁気モーメント法計算関数
void Magnetic_Moment_Method(mpsconfig &CON,vector<mpselastic>&PART,double **F,double n0,double lamda,int fluid_number,int particle_number, double current_time, int t)
{	
	double u0=4*PI*1e-7;
	MMM_config MCON;

	cout<<"磁気モーメント法---";

	double le=CON.get_distancebp();
	double R=CON.get_re()*le;	//発散用影響半径
	//int d=CON.get_dimention();						//次元
	int d = 3;
	unsigned int timeA=GetTickCount();				//計算開始時刻
	int count=0;
	int pn=fluid_number*d;							//未知数:粒子数×速度D(次元)成分 
	double co=1/(4*PI*u0);							//計算によく現れる係数
	double RP0=CON.get_RP();						//比透磁率
	double kai=(RP0-1.0);								//磁気感受率
	double l=pow(CON.get_particle_volume(), double(1/3));
	double S=le;
	double mt = 30000;

	//double mag_center[3]={0,0,-2.1e-05};
	double mag_center[3]={0,0,-0.3};
	double mag_M[3] = {0,0,1.3};
	//double size[3] = {0.003001, 0.0005001, CON.get_distancebp() / 2.};
	double size[3] = {0.5, 0.5, 0.5};

	if(d==3 && CON.get_model_set_way()==0) S=le*le;
	if(d==3 && CON.get_model_set_way()==1) S=sqrt(3.0)/4*le*le;//断面積

	//行列確保 //未知数の順番はMx0,My0,Mz0,Mx1,My1,Mz1,Mx2,My2,Mz2・・・
	vector<vector<double> > matrix(pn, vector<double>(pn));
	vector<double> B(pn);
	vector<double> Bo(pn);
	double *direct[DIMENSION];
	for(int D=0;D<DIMENSION;D++) direct[D]=new double [particle_number];//外向き法線ベクトル

	double max_z = 0.;
	for(int i = 0; i < fluid_number; i++){
		if(max_z < PART[i].r[2]){
			if(max_z < PART[i].r[2] && sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1])<0.002){
				max_z = PART[i].r[2];
			}
		}
	}


	if(t==1)
	{
		ofstream t1("height.dat");		//横軸ステップ数、縦軸dtのグラフ
		t1.close();
	}


	ofstream t1("height.dat",ios :: app);
	t1<<current_time<<" "<<max_z<<endl;
	t1.close();

	double z_max = 0;
	double z_min = 1;
	for(int i = 0; i < fluid_number; i++){
		if(z_max < PART[i].r[d-1])z_max = PART[i].r[d-1];
		if(z_min > PART[i].r[d-1])z_min = PART[i].r[d-1];
	}
	ofstream fh;
	cout<<"Diameter="<<z_max-z_min<<endl;
	if(t==1){
		ofstream fh("./freq.dat");
		fh<<current_time<<" "<<z_max<<endl;
		fh.close();
	}else{
		ofstream fh("./freq.dat" , iostream::app);
		fh<<current_time<<" "<<z_max<<endl;
		fh.close();
	}

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++) direct[D][i]=0;				//初期化
		if(PART[i].surface==ON)
		{
			direct_f(CON,PART,i,direct);
			for(int D=0;D<d;D++) direct[D][i]*=-1;//外向きが欲しいから反転する
		}
	}

	for(int i=0;i<pn;i++)
	{
		B[i]=0;
		for(int j=0;j<pn;j++) matrix[i][j]=0;		//初期化
	}
	//////////////////////////////////////////////////////////////////

	//解行列B[]作成
	if(MCON.get_Hf_type()==0)
	{
		cout<<"\n一様磁場\n";
		double sign[3]={0,0,0};
		sign[d-1]=1;
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<d;D++) B[i*d+D]=MCON.get_Hf_H()*sign[D] * (RP0 - 1) * MU_0;//2次元ならA_Y,3次元ならA_Zの方向だけ値が入る。
			//for(int D=0;D<d;D++) Bo[i*d+D]=MCON.get_Hf_H()*sign[D];
		}
		cout<<"飽和磁化の考慮"<<endl;
		for(int i=0;i<fluid_number;i++)
		{
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(B[i * d + D],2.0);
			}
			nB = sqrt(nB);
			double ratio;
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
			}
		}
	}
	else if(MCON.get_Hf_type()==1)
	{
		cout<<"\n磁石\n";
		for(int i=0;i<fluid_number;i++)
		{
			double ki[3][3] = {.0};
			//double p[3] = {PART[i].r[A_X]-mag_center[A_X], PART[i].r[A_Y]-mag_center[A_Y], 0};
			double p[3] = {PART[i].r[A_X]-mag_center[A_X], PART[i].r[A_Y]-mag_center[A_Y], PART[i].r[A_Z]-mag_center[A_Z]};
			Cal_Mag_Density(ki, p, size[0], size[1], size[2]);
			for(int j = 0; j < d; j++){
				for(int k = 0; k < d; k++){
					//Bo[i * d + j] += ki[k][j] * mag_M[k] / (-4.0 * PI);
					B[i * d + j] += ki[k][j] * mag_M[k] / (-4.0 * PI) * (RP0 - 1);
				}
			}
		}
		cout<<"飽和磁化の考慮"<<endl;
		for(int i=0;i<fluid_number;i++)
		{
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(B[i * d + D],2.0);
			}
			nB = sqrt(nB);
			double ratio;
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
			}
		}
	}
	else cout<<"外部磁場のタイプが未実装"<<endl;
	/////
	//磁気モーメント法で，要素の面が厳密にくっついてしまうと，誤差が大きくなる．
	//なので，計算上粒子のサイズを1%縮める．こっちの方が数段まし．
	/////////////////////matrixに値を格納
	for(int dst=0;dst<fluid_number;dst++)
	{
		for(int src=0;src<fluid_number;src++)
		{
			double ki[3][3] = {.0};
			//double k1[3][3] = {.0};
			//double k2[3][3] = {.0};
			//double rA[3]={PART[dst].r[A_X]-PART[src].r[A_X],PART[dst].r[A_Y]-PART[src].r[A_Y],0};
			double rA[3]={PART[dst].r[A_X]-PART[src].r[A_X],PART[dst].r[A_Y]-PART[src].r[A_Y],PART[dst].r[A_Z]-PART[src].r[A_Z]};
			//if(sqrt(rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2]) > CON.get_distancebp() * 7)continue;
			Cal_Mag_Density(ki, rA, CON.get_distancebp()/2, CON.get_distancebp()/2 , CON.get_distancebp()/2);
			//変更
			//Cal_Mag_Density_To_Matrix(k1, k2, rA, CON.get_distancebp() / 4.0);
			for(int i = 0; i < d; i++){
				for(int j = 0; j < d; j++){
					//matrix[dst * d + i][src * d + j] += k1[i][j] / (-4 * PI) * (RP0 - 1);
					//matrix[src * d + i][dst * d + j] += k2[i][j] / (-4 * PI) * (RP0 - 1);
					//磁束密度ベースの式を用いていることに注意！
					if(src == dst && i == j){
						matrix[dst * d + i][src * d + j] += (ki[j][i] + (4 * PI)) / (4 * PI) * (RP0 - 1) + 1;
					}else matrix[dst * d + i][src * d + j] += ki[j][i] / (4 * PI) * (RP0 - 1);
				}
			}

		}
	}
	cout<<"行列作成"<<endl;

	//行列値出力
	//output_matrix(matrix,B, pn);

	cout<<"未知数:"<<pn<<"ガウスザイデル法--";
	unsigned int timeB=GetTickCount();
	GaussSeidel(matrix, pn, B, d);
	//jacobi(matrix,B,pn);//答えはBに格納
	cout<<"ok time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

	double *M[DIMENSION];					
	for(int D=0;D<DIMENSION;D++) M[D]=new double [fluid_number];
	double *H[DIMENSION];					
	for(int D=0;D<DIMENSION;D++) H[D]=new double [fluid_number];//粒子位置での磁場H

	cout<<"飽和磁化の考慮"<<endl;
	for(int i=0;i<fluid_number;i++)
	{
		double nB = 0.0;
		for(int D=0;D<d;D++){
			nB += pow(B[i * d + D],2.0);
		}
		nB = sqrt(nB);
		double ratio;
		nB > mt ? ratio = mt / nB : ratio = 1.0;
		if(MCON.get_Hf_type() == 0){
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
				M[D][i]=B[i*d+D];
			}
		}else if(MCON.get_Hf_type() == 1){
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
				M[D][i]=B[i*d+D];

			}
		}
	}

	//粒子位置での磁場H
	for(int i = 0; i < fluid_number; i++)for(int D = 0; D < d; D++)H[D][i] = 0;
	cout<<"外部磁場計算---";
	if(MCON.get_Hf_type() == 0){
		double sign[3]={0,0,0};
		sign[d-1]=1;
		for(int i=0; i < fluid_number; i++){
			for(int D=0;D<d;D++) H[D][i]=MCON.get_Hf_H()*sign[D];
		}
	}

	//ある点における磁束密度を出す(丸め誤差が気色悪いので，1e+7をかける)
	//外部磁場については，計算のしやすさとかもあるので，ソースの磁場だけにしておく．
	else if(MCON.get_Hf_type() == 1){
		for(int i=0;i<fluid_number;i++)
		{
			double ki[3][3] = {.0};
			double p[3] = {PART[i].r[A_X]-mag_center[A_X], PART[i].r[A_Y]-mag_center[A_Y], PART[i].r[A_Z]-mag_center[A_Z]};
			Cal_Mag_Density(ki, p, size[0], size[1], size[2]);
			for(int j = 0; j < d; j++){
				for(int k = 0; k < d; k++){
					H[j][i] += ki[k][j] * mag_M[k] / (-16 * PI * PI) * 10000000;
				}
			}
		}
	}

	/*
	for(int dst=0;dst<fluid_number;dst++)
	{
	for(int src=0;src<fluid_number;src++)
	{
	if(dst == src) continue;
	double ki[3][3] = {.0};
	double rA[3]={PART[dst].r[A_X]-PART[src].r[A_X],PART[dst].r[A_Y]-PART[src].r[A_Y],0};
	//if(sqrt(rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2]) > CON.get_distancebp() * 7)continue;
	Cal_Mag_Density(ki, rA, CON.get_distancebp()/2, CON.get_distancebp()/2, 0.0005);
	for(int i = 0; i < d; i++){
	for(int j = 0; j < d; j++){
	//if(dst == src && i == j)H[i][dst] += (ki[j][i] + 4 * PI) * M[j][src] / (-4.0 * PI);
	//else H[i][dst] += ki[j][i] * M[j][src] / (-4.0 * PI);
	H[i][dst] += ki[j][i] * M[j][src] / (-16 * PI * PI) * 10000000;
	}
	}
	}
	}*/

	ofstream fp("M.dat");
	ofstream fq("H.dat");
	double times=1e-11;
	if(d==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<M[A_X][i]<<" "<<M[A_Y][i]<<endl;
			fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i])<<endl;
		}
	}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]<<" "<<M[A_Y][i]<<" "<<M[A_Z][i]<<endl;
			fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
		}
	}
	fp.close();
	fq.close();

	//力を求める
	double *Fs[DIMENSION];
	for(int D=0;D<DIMENSION;D++) Fs[D]=new double [particle_number];//単位面積あたりの力
	double *Fv[DIMENSION];
	for(int D=0;D<DIMENSION;D++) Fv[D]=new double [particle_number];//単位体積あたりの力
	double *Hgrad[DIMENSION];			//H勾配ﾍﾞｸﾄﾙ格納
	double *Hgrad1[DIMENSION];
	double *Hgrad2[DIMENSION];
	double *Gaiseki[DIMENSION];
	double *naiseki = new double[fluid_number];
	for(int D=0;D<DIMENSION;D++){
		Hgrad[D]=new double [fluid_number];
		Hgrad1[D]=new double [fluid_number];
		Hgrad2[D]=new double [fluid_number];
		Gaiseki[D]=new double [fluid_number];
	}

	if(MCON.get_force_t() == 0){
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENSION;D++) Fs[D][i]=0;		//初期化
			if(PART[i].surface==ON)
			{
				double H_n = H[0][i] * direct[0][i] + H[1][i] * direct[1][i] + H[2][i] * direct[2][i];
				//double Hr = sqrt(H_n[0] * H_n[0] + H_n[1] * H_n[1] + H_n[2] * H_n[2]);
				//double H_n = M[0][i] / u0 * 2 * direct[0][i] + M[1][i]  / u0 * 2 * direct[1][i] + M[2][i] / u0 * 2 * direct[2][i];			
				for(int D=0;D<DIMENSION;D++) Fs[D][i] = 0.5 * u0 * H_n * H_n * direct[D][i];
				for(int D=0;D<DIMENSION;D++) F[D][i]=Fs[D][i]*CON.get_distancebp();
			}
		}
	}


	//体積力
	else{
		double Fx = 0, Fy = 0;
		if(MCON.get_eForce() == 0){
			/*
			for(int i = 0; i < fluid_number; i++)for(int D = 0; D < d; D++)H[D][i] *= M[D][i];
			H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENSION;D++)F[D][i] = Hgrad[D][i];
			*/

			for(int i = 0; i < fluid_number; i++){

				double f[3] = {0.};
				//WLSMによる勾配の計算
				//今回はre2を使う

				cal_WLSM_for_grad(PART, i, CON.get_re() * le, M, H, f);
				for(int D=0;D<DIMENSION;D++){
					F[D][i] = f[D];
				}
			}

			//各方向の電磁力を出力
		}else if(MCON.get_eForce() == 1){
			H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
			//H_gradient2(CON,PART, fluid_number,Hgrad,H[0]);//∇Hx計算
			//H_gradient2(CON,PART, fluid_number,Hgrad1,H[1]);//∇Hy計算
			//H_gradient2(CON,PART, fluid_number,Hgrad2,H[2]);//∇Hz計算
			for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++)Fv[D][i] = Hgrad[D][i] * sqrt(M[0][i] * M[0][i] + M[1][i] * M[1][i] + M[2][i]* M[2][i]);
			//for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
			//	Fv[0][i] = M[0][i] * Hgrad[0][i] + M[1][i] * Hgrad[1][i] + M[2][i] * Hgrad[2][i];
			//	Fv[1][i] = M[0][i] * Hgrad1[0][i] + M[1][i] * Hgrad1[1][i] + M[2][i] * Hgrad1[2][i];
			//	Fv[2][i] = M[0][i] * Hgrad2[0][i] + M[1][i] * Hgrad2[1][i] + M[2][i] * Hgrad2[2][i];
			//}

			//偶力
			/*
			for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
			Gaiseki[0][i] = (M[1][i]*H[2][i] - M[2][i]*H[1][i]);
			Gaiseki[1][i] = (M[2][i]*H[0][i] - M[0][i]*H[2][i]);
			Gaiseki[2][i] = (M[0][i]*H[1][i] - M[1][i]*H[0][i]);
			}
			H_gradient2(CON,PART, fluid_number,Hgrad,Gaiseki[0]);//∇(M×H)x計算
			H_gradient2(CON,PART, fluid_number,Hgrad1,Gaiseki[1]);//∇(M×H)y計算
			H_gradient2(CON,PART, fluid_number,Hgrad2,Gaiseki[2]);//∇(M×H)z計算
			for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
			Fv[0][i] += .5*(Hgrad1[2][i] - Hgrad2[1][i]);
			Fv[1][i] += .5*(Hgrad2[0][i] - Hgrad[2][i]);
			Fv[2][i] += .5*(Hgrad[1][i] - Hgrad1[0][i]);
			}*/
		}
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENSION;D++){
			//F[D][i] += Fv[D][i] * CON.get_distancebp() * CON.get_distancebp() * CON.get_distancebp();
			PART[i].F[D] = F[D][i];
		}
	}

	//test計算領域全体の磁束密度の計算
	int x = (CON.get_maxX() - CON.get_minX()) / (2*le);
	int y = (CON.get_maxY() - CON.get_minY()) / (2*le);

	ofstream bb("B_Area.dat");
	for(int i = 0; i < x; i++){
		for(int j = 0; j < y; j++){
			double Barea[3] = {.0};
			for(int k = 0; k < fluid_number; k++){
				double ki[3][3] = {0.};
				double rA[3] = {(i * 2 * le) + CON.get_minX() + le / 2 - PART[k].r[0], (j * 2 * le) + CON.get_minY() + le / 2 - PART[k].r[1], 0};
				Cal_Mag_Density(ki, rA, CON.get_distancebp()/2, CON.get_distancebp()/2, 0.0005);
				for(int d = 0; d < 3; d ++){
					for(int dd = 0; dd < 3; dd++){
						Barea[d] += ki[d][dd] * M[dd][k] / (-4 * PI);
					}
				}
			}
			double ki[3][3] = {0.};
			double rA[3] = {(i * 2 * le) + CON.get_minX() + le / 2 - mag_center[0], (j * 2 * le) + CON.get_minY() + le / 2 - mag_center[1], 0};
			Cal_Mag_Density(ki, rA, size[0], size[1], size[2]);
			for(int d = 0; d < 3; d ++){
				for(int dd = 0; dd < 3; dd++){
					Barea[d] += ki[d][dd] * mag_M[dd] / (-4 * PI);
				}
			}
			bb<<(i * 2 * le) + CON.get_minX() + le / 2 <<" "<<(j * 2 * le) + CON.get_minY() + le / 2<<" "<<Barea[0]<<" "<<Barea[1]<<endl;
		}
	}
	bb.close();

	////ｽﾑｰｼﾞﾝｸﾞ
	//smoothingF3D(CON,PART,fluid_number,F);

	ofstream fr("Fs.dat");
	times=5e-2/MCON.get_Hf_H();
	if(d==2){ for(int i=0;i<fluid_number;i++) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Y][i]*times<<endl;}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Z][i]*times<<endl;
	}
	fr.close();
	ofstream ft("Fv.dat");
	times=1e+3;
	if(d==3) for(int i=0;i<fluid_number;i++) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Y]*times<<endl;
	ft.close();

	if((CON.get_dir_for_P()==2 ||CON.get_dir_for_P()==3) && MCON.get_force_t() == 0 )
	{
		//ofstream bb("electromagnetic_P.dat");
		for(int i=0;i<fluid_number;i++)
		{
			double fs=0;//表面力
			if(PART[i].surface==ON)
			{
				fs=sqrt(Fs[A_X][i]*Fs[A_X][i]+Fs[A_Y][i]*Fs[A_Y][i]+Fs[A_Z][i]*Fs[A_Z][i]);
				for(int D=0;D<DIMENSION;D++) F[D][i]=0;
			}
			PART[i].dir_Pem=-fs;
			//bb<<-fs<<endl;
		}
		//bb.close();
	}

	for(int D=0;D<DIMENSION;D++) delete [] M[D];
	for(int D=0;D<DIMENSION;D++) delete [] H[D];
	for(int D=0;D<DIMENSION;D++) delete [] direct[D];
	for(int D=0;D<DIMENSION;D++) delete [] Fs[D];
	for(int D=0;D<DIMENSION;D++) delete [] Fv[D];
	for(int D=0;D<DIMENSION;D++) delete [] Hgrad[D];


}

///20140904磁気モーメント法メッシュベース・・・高速化とメモリ節約はコレ
void Magnetic_Moment_Methodv2(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER, double **F,double n0,double lamda,int fluid_number,int particle_number, double current_time, int t)
{	
	double u0=4*PI*1e-7;
	MMM_config MCON;

	cout<<"磁気モーメント法---メッシュベース"<<endl;

	double le=CON.get_distancebp();
	double R=CON.get_re()*le;	//発散用影響半径
	int d=CON.get_dimension();						//次元
	unsigned int timeA=GetTickCount();				//計算開始時刻
	int count=0;
	int pn=fluid_number*d;							//未知数:粒子数×速度D(次元)成分 
	double co=1/(4*PI*u0);							//計算によく現れる係数
	double RP0=CON.get_RP();						//比透磁率
	double kai=(RP0-1.0);								//磁気感受率
	//double l=pow(CON.get_particle_volume(), double(1/3));
	double S=le;
	double mt = 1;                                //流体に対す磁性体の占める体積によって変化させる
	double ratio;
	//double mag_center[3]={0,0,-2.1e-05};
	double mag_center_Z_init =0.02452; //磁石の重心の初期地
	double mag_center_Z_goal =0.02452; //磁石の重心の目標地．初期地と同じにすれば動かない．
	double mag_center_dx = 0.2;       //磁石の移動速度(m/s)
	double mag_center_Z_curr = mag_center_Z_init + current_time * mag_center_dx;         //磁石の現在地の更新
	double mag_M[3] = {0,0,0.1};    //[T](交流の場合は振幅)


	int mesh_size;
	double size[3] = {0.005, 0.005, 0.005};	//各辺の半分の寸法


	/*//------交流磁場
	double mag_Hz = 7.0;
	//double duty = (mag_Hz * 0.15) / 7.25;
	double duty = 0.15;
	double ct = current_time;
	while(ct > (1.0 / mag_Hz)){
	ct -= (1.0 / mag_Hz);
	}
	if(ct > (1.0 / mag_Hz) * duty)mag_M[2] = 0;
	//impulse
	if(t > 500)mag_M[2]=0;else mag_M[2]=0.5;
	//------交流磁場F.dat
	*/



	if(mag_center_Z_curr > mag_center_Z_goal)mag_center_Z_curr = mag_center_Z_goal;
	double mag_center[3] = {0,0,mag_center_Z_curr};

	cout<<"磁石の位置={0,0,"<<mag_center_Z_curr<<"}"<<endl;
	cout<<"磁石の強さ={0,0,"<<mag_M[2]<<"}"<<endl;

	if(d==3 && CON.get_model_set_way()==0) S=le*le;
	if(d==3 && CON.get_model_set_way()==1) S=sqrt(3.0)/4*le*le;//断面積
	double mesh_scale =2.;//2.72;                                 //磁気モーメント法用メッシュの大きさ　磁石の寸法によっては変える必要あり　2
	double rm = 3;                                      //メッシュの中心から見込んだメッシュの数(1.9:26近傍,1:6近傍)．
	double width=CON.get_distancebp()*mesh_scale;		//格子幅
	int mesh_x = ((CON.get_maxX() - CON.get_minX()) / width + 0.5);
	int mesh_y = ((CON.get_maxY() - CON.get_minY()) / width + 0.5);
	int mesh_z = ((CON.get_maxZ() - CON.get_minZ()) / width + 0.5);
	mesh_size = mesh_x * mesh_y * mesh_z;

	vector<vector<double>>MESH;              //粒子の入っているメッシュの中心座標
	vector<int>mesh_index(fluid_number);     //粒子の入っているメッシュの番号(粒子の入っているメッシュのみでスケーリング)
	vector<int>mesh_p_num;	                 //メッシュ内に入っている粒子の数
	vector<vector<int>>p_index;              //メッシュ内にある粒子の番号

	//メッシュ内にどれだけの流体粒子があるか
	vector<int> mesh_i;
	vector<vector<int>> MESH_ALL(mesh_size, vector<int>(5));
	vector<vector<int>> MESH_N(mesh_size, vector<int>(0));
	vector<vector<int>> part_on_mesh(mesh_size, vector<int>(0));
	mesh_i.reserve(fluid_number);
	part_on_mesh.reserve(fluid_number);
	for(int i = 0; i < mesh_size; i++){
		MESH_ALL[i][3] = 0;
	}

	//電磁力関係
	double *Fs[DIMENSION];
	for(int D=0;D<DIMENSION;D++) Fs[D]=new double [particle_number];//単位面積あたりの力
	double *Fv[DIMENSION];
	for(int D=0;D<DIMENSION;D++) Fv[D]=new double [particle_number];//単位体積あたりの力
	double *Hgrad[DIMENSION];			//H勾配ﾍﾞｸﾄﾙ格納
	double *naiseki = new double[fluid_number];
	for(int D=0;D<DIMENSION;D++){
		Hgrad[D]=new double [fluid_number];
	}

	//磁化と磁界強度
	double *M[DIMENSION];					
	for(int D=0;D<DIMENSION;D++) M[D]=new double [fluid_number];
	double *H[DIMENSION];					
	for(int D=0;D<DIMENSION;D++) H[D]=new double [fluid_number];
	for(int i = 0; i < fluid_number; i++)for(int D = 0; D < DIMENSION; D++)H[D][i] = 0;



	double max_z = 0.;
	for(int i = 0; i < fluid_number; i++){
		if(max_z < PART[i].r[2]){
			if(max_z < PART[i].r[2] && sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1])<0.001){
				max_z = PART[i].r[2];
			}
		}
	}
	double z_max = 0;
	double z_min = 1;
	for(int i = 0; i < fluid_number; i++){
		if(z_max < PART[i].r[d-1]&&sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1]) < 2*le)z_max = PART[i].r[d-1];
		if(z_min > PART[i].r[d-1]&&sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1]) < 2*le)z_min = PART[i].r[d-1];
	}
	double x_max = 0;
	double x_min = 1;
	for(int i = 0; i < fluid_number; i++){
		if(x_max < PART[i].r[0])x_max = PART[i].r[0];
		if(x_min > PART[i].r[0])x_min = PART[i].r[0];
	}
	ofstream fh;
	cout<<"Diameter="<<x_max-x_min<<endl;
	cout<<"  Height="<<z_max-z_min<<endl;

	if(t==1){
		ofstream fh("./freq.dat");
		fh<<current_time<<" "<<z_max-z_min<<endl;
		fh.close();
	}else{
		ofstream fh("./freq.dat" , iostream::app);
		fh<<current_time<<" "<<z_max-z_min<<endl;
		fh.close();
	}

	//=======================================================================
	//	出力
	/////////////////////////////////////////////////////////////////////////
	if(t==1)
	{
		ofstream t1("height.dat");		//横軸時間、縦軸流体高さのグラフ
		t1.close();
	}
	ofstream t1("height.dat",ios :: app);
	t1<<current_time<<" "<<max_z-0.005<<endl;
	t1.close();

	//////////////////////////////////////////////////////////////////////////

	//バックグラウンドセルの番号を設定
	for(int i = 0; i < mesh_x; i++){
		for(int j = 0; j < mesh_y; j++){
			for(int k = 0; k < mesh_z; k++){
				MESH_ALL[k * mesh_x * mesh_y + j * mesh_x + i][0] = i;
				MESH_ALL[k * mesh_x * mesh_y + j * mesh_x + i][1] = j;
				MESH_ALL[k * mesh_x * mesh_y + j * mesh_x + i][2] = k;
			}
		}
	}

	//流体粒子がどのセルに属しているのか計算
	for(int i=0;i<fluid_number;i++)
	{
		int X=(int)((PART[i].r[A_X]-CON.get_minX())/width + 0.000001);//X方向に何個目の格子か 
		int Y=(int)((PART[i].r[A_Y]-CON.get_minY())/width + 0.000001);//Y方向に何個目の格子か
		int Z=(int)((PART[i].r[A_Z]- le -CON.get_minZ())/width + 0.000001);//Z方向に何個目の格子か
		int number=Z*mesh_x*mesh_y+Y*mesh_x+X;//粒子iを含む格子の番号
		MESH_ALL[number][3]++;
		mesh_i.emplace_back(number);
		part_on_mesh[number].emplace_back(i);
	}

	cout<<"メッシュの排除"<<endl;
	//計算の都合上，いらないメッシュを排除
	//ここでいういらないメッシュとは，流体が入っていないセルのこと
	for(int i = 0; i < mesh_size; i++){
		vector<double>work(3);
		vector<int>p_work;
		//セルiの中の粒子の数が0でないか検査
		if(MESH_ALL[i][3] != 0){
			work[0] = (((double)MESH_ALL[i][0] + 0.5) * width) + CON.get_minX();
			work[1] = (((double)MESH_ALL[i][1] + 0.5) * width) + CON.get_minY();
			work[2] = (((double)MESH_ALL[i][2]) * width) + CON.get_minZ() + le;
			//セルiの情報を計算に関与するセル情報に追加
			MESH.push_back(work);
			//追加したセルの中に入っている粒子の数をコピー
			mesh_p_num.push_back(MESH_ALL[i][3]);
			//流体粒子の入っているセルを改めてナンバリング
			for(int j = 0; j < MESH_ALL[i][3]; j++){
				p_work.push_back(part_on_mesh[i][j]);
				mesh_index[part_on_mesh[i][j]] = i;
			}
			//ナンバリングしたセルと粒子の番号を関連づけ
			p_index.push_back(p_work);
		}
	}

	cout<<"隣接メッシュ"<<endl;
	//あるメッシュに対して，隣接するメッシュがいくつあるのか計算する
	for(int i = 0; i < MESH.size(); i++){
		for(int j = 0; j < MESH.size(); j++){
			if(i == j)continue; //自分自身は計算に入れない
			double xm = (MESH[i][0] - MESH[j][0]);
			double ym = (MESH[i][1] - MESH[j][1]);
			double zm = (MESH[i][2] - MESH[j][2]);
			double r = sqrt(xm * xm + ym * ym + zm * zm);
			if(r <= rm * width){
				MESH_N[i].push_back(j);
			}
		}
	}


	cout<<"対象となるメッシュ:"<<MESH.size()<<endl;

	//行列確保 //未知数の順番はMx0,My0,Mz0,Mx1,My1,Mz1,Mx2,My2,Mz2・・・
	//メッシュごとの磁化を求めようとしている．
	vector<vector<double> > matrix(MESH.size() * 3, vector<double>(MESH.size() * 3));
	vector<double> B(MESH.size() * 3);
	double *direct[DIMENSION];
	for(int D=0;D<DIMENSION;D++) direct[D]=new double [particle_number];//外向き法線ベクトル

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++) direct[D][i]=0;				//初期化
		if(PART[i].surface==ON)
		{
			direct_f(CON,PART,i,direct);
			for(int D=0;D<d;D++) direct[D][i]*=-1;//外向きが欲しいから反転する
		}
	}

	//この段階では，磁気モーメント法での要素数はナンバリングされたセルの数になっている
	for(int i=0;i<MESH.size();i++)
	{
		B[i]=0;
		for(int j=0;j<MESH.size();j++) matrix[i][j]=0;		//初期化
	}
	//////////////////////////////////////////////////////////////////

	//解行列B[]作成

	//各セル内で磁化を求める．
	cout<<"解行列の作成-----"<<endl;
	if(MCON.get_Hf_type()==0)
	{
		cout<<"\t一様磁場---"<<endl;
		double sign[3]={0,0,0};
		sign[d-1]=1;
		for(int i=0;i<MESH.size();i++)
		{
			for(int D=0;D<d;D++) B[i*d+D]=MCON.get_Hf_H()*sign[D] * (RP0 - 1) * MU_0;//2次元ならA_Y,3次元ならA_Zの方向だけ値が入る。
		}
		cout<<"OK"<<endl;
		cout<<"\t飽和磁化の考慮---";
		for(int i=0;i<MESH.size();i++)
		{
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(B[i * d + D],2.0);
			}
			nB = sqrt(nB);
			double ratio;
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
			}
		}
		cout<<"OK"<<endl;
	}
	else if(MCON.get_Hf_type()==1)
	{
		cout<<"\t磁石からの影響計算---";
		for(int i=0;i<MESH.size();i++)
		{
			double ki[3][3] = {.0};
			double p[3] = {MESH[i][0] - mag_center[0], MESH[i][1] - mag_center[1], MESH[i][2] - mag_center[2]};
			Cal_Mag_Density(ki, p, size[0], size[1], size[2]);
			for(int j = 0; j < d; j++){
				for(int k = 0; k < d; k++){
					B[i * d + j] += ki[k][j] * mag_M[k] / (-4.0 * PI) * (RP0 - 1);
				}
			}
		}
		cout<<"OK"<<endl;
		cout<<"\t飽和磁化の考慮---";
		for(int i=0;i<MESH.size();i++)
		{
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(B[i * d + D],2.0);
			}
			nB = sqrt(nB);
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
			}
		}
		cout<<"OK"<<endl;
	}else cout<<"外部磁場のタイプが未実装"<<endl;
	cout<<"-----完了"<<endl;
	/////
	if(mag_M[2]!=0){
		/////////////////////matrixに値を格納
		for(int dst=0;dst<MESH.size();dst++)
		{
			for(int src=0;src<MESH.size();src++)
			{
				double px = MESH[dst][0] - MESH[src][0];
				double py = MESH[dst][1] - MESH[src][1];
				double pz = MESH[dst][2] - MESH[src][2];
				double ki[3][3] = {.0};
				double rA[3]={px, py, pz};
				//Cal_Mag_Density(ki, rA, (width / 2.), (width / 2.), 0);
				Cal_Mag_Density(ki, rA, (width / 2.));
				for(int i = 0; i < d; i++){
					for(int j = 0; j < d; j++){
						if(dst == src && i == j){
							matrix[dst * d + i][src * d + j] += (ki[i][j] + 4 * PI) / (4 * PI) * (RP0 - 1) + 1;
						}else matrix[dst * d + i][src * d + j] += ki[i][j] / (4 * PI) * (RP0 - 1);
					}
				}
			}
		}
		cout<<"行列作成"<<endl;

		//行列値出力
		//output_matrix(matrix,B, pn);

		cout<<"magcenter="<<mag_center[0]<<","<<mag_center[1]<<","<<mag_center[2]<<endl;

		cout<<"未知数:"<<MESH.size() * 3<<".ガウスザイデル法--";
		unsigned int timeB=GetTickCount();
		//行列を解く（ガウスザイデル法）
		GaussSeidel(matrix, MESH.size() * 3, B, d);
		//jacobi(matrix,B,pn);//答えはBに格納
		cout<<"OK time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

		//粒子位置での磁場H
		if(MCON.get_isWLSM() == 0){
			cout<<"重みつき最小二乗法による推定"<<endl;
			for(int i=0;i<MESH.size();i++)
			{
				//関数近似後の各係数
				double cx[9] = {.0};
				double cy[9] = {.0};
				double cz[9] = {.0};
				//ある集合体内で最小二乗法を行う
				cal_WLSM_for_M2(MESH, MESH_N, B, i, width, width * rm, cx, cy, cz);
				for(int j = 0; j <p_index[i].size(); j++){
					double Mb[3] = {0};
					//係数から変化量を求める．この中で関数を再現して同一セル内にある粒子の座標に応じた磁化を付与する．
					get_M2(cx, cy, cz, PART[p_index[i][j]].r[A_X], PART[p_index[i][j]].r[A_Y], PART[p_index[i][j]].r[A_Z], MESH[i][0], MESH[i][1], MESH[i][2], Mb);
					for(int D=0;D<d;D++)M[D][p_index[i][j]]=Mb[D] + B[i * 3 + D];
				}
			}
		}
		//最小二乗法を用いた磁化補間ができた

		cout<<"OK"<<endl;
		//一応飽和磁化の考慮をしておく．なぜなら，近似関数の関係で飽和磁化より大きい値を持つ可能性があるから．
		cout<<"飽和磁化の考慮---";
		for(int i = 0; i < fluid_number; i++){
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(M[D][i],2.0);
			}
			nB = sqrt(nB);
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				M[D][i] *= ratio;
			}
		}
		cout<<"OK"<<endl;


		cout<<"外部磁場計算---";
		//磁石による磁性流体に対する外部磁場の計算
		if(MCON.get_Hf_type()!=0){
			cout<<"磁石"<<endl;
			for(int i=0;i<fluid_number;i++)
			{
				double ki[3][3] = {.0};
				double p[3] = {PART[i].r[A_X]-mag_center[A_X], PART[i].r[A_Y]-mag_center[A_Y], PART[i].r[A_Z]-mag_center[A_Z]};
				Cal_Mag_Density(ki, p, size[0], size[1], size[2]);
				for(int j = 0; j < d; j++){
					for(int k = 0; k < d; k++){
						H[j][i] = ki[j][k] * mag_M[k] / (-16. * PI * PI) * 10000000;
					}
				}
			}
		}

		/*
		//本来は，全粒子からの影響を計算すべきではあるものの，計算量が膨大となるので，
		//ここでは，影響半径内に入った粒子からの磁場の寄与のみを考える．
		for(int dst=0;dst<fluid_number;dst++)
		{
		for(int src=0;src<PART[dst].N;src++)
		{
		if(dst == src)continue; //自己場はいらない．
		double k1[3][3] = {.0};
		double rA[3]={PART[dst].r[A_X]-PART[PART[dst].NEI[src]].r[A_X],PART[dst].r[A_Y]-PART[PART[dst].NEI[src]].r[A_Y],PART[dst].r[A_Z]-PART[PART[dst].NEI[src]].r[A_Z]};
		double r = sqrt(rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2]);
		double l = CON.get_distancebp();
		Cal_Mag_Density(k1, rA, l / 2.);
		for(int i = 0; i < d; i++){
		for(int j = 0; j < d; j++){
		H[i][dst] += k1[i][j] * M[j][src] / (-16. * PI * PI) * 10000000;
		}
		}
		}
		}*/
		cout<<"OK"<<endl;

		//応力を計算する場合
		if(MCON.get_force_t() == 0){
			cout<<"応力"<<endl;
			for(int i=0;i<fluid_number;i++)
			{
				for(int D=0;D<DIMENSION;D++) Fs[D][i]=0;		//初期化
				if(PART[i].surface==ON)
				{
					for(int D=0;D<DIMENSION;D++)H[D][i] = M[D][i] / (u0 * (RP0 - 1.));
					//cout<<"a"<<endl;
					//double Mn=(H[A_X][i]*direct[A_X][i]+H[A_Y][i]*direct[A_Y][i]+H[A_Z][i]*direct[A_Z][i]);
					//double Mn=(M[A_X][i]*direct[A_X][i]+M[A_Y][i]*direct[A_Y][i]+M[A_Z][i]*direct[A_Z][i]) / (RP0 - 1.);
					//double val=0.5 * Mn * Mn / u0;		//応力値
					//double val = 0.5 * Mn * Mn * u0;		//応力値
					//for(int D=0;D<DIMENSION;D++) Fs[D][i]=val*direct[D][i];
					//for(int D=0;D<DIMENSION;D++) F[D][i]=Fs[D][i]*CON.get_distancebp();
					//一項目
					double H_n = H[0][i] * direct[0][i] + H[1][i] * direct[1][i] + H[2][i] * direct[2][i];
					//double Mn = M[0][i] * direct[0][i] + M[1][i]
					//double Hr = H[0][i] * H[0][i] + H[1][i] * H[1][i] + H[2][i] * H[2][i];
					//double H_n = M[0][i] / u0 * 2 * direct[0][i] + M[1][i]  / u0 * 2 * direct[1][i] + M[2][i] / u0 * 2 * direct[2][i];			
					for(int D=0;D<DIMENSION;D++) Fs[D][i] = 0.5 * H_n * H_n * u0 * direct[D][i];
					//二項目
					//double H_H = H[0][i] * H[0][i] + H[1][i] * H[1][i] + H[2][i] * H[2][i];
					//合計
					//for(int D=0;D<DIMENSION;D++) Fs[D][i] -= (.5 * u0 * (H_H * direct[D][i]));
					//for(int D=0;D<DIMENSION;D++) F[D][i] = Fs[D][i] * CON.get_distancebp() * CON . get_distancebp();
				}
			}
			if(CON.get_dir_for_P()==2 ||CON.get_dir_for_P()==3 )
			{
				//ofstream bb("electromagnetic_P.dat");
				for(int i=0;i<fluid_number;i++)
				{
					double fs=0;//表面力
					if(PART[i].surface==ON)
					{
						fs=sqrt(Fs[A_X][i]*Fs[A_X][i]+Fs[A_Y][i]*Fs[A_Y][i]+Fs[A_Z][i]*Fs[A_Z][i]);
						for(int D=0;D<DIMENSION;D++) F[D][i]=0;
					}
					PART[i].dir_Pem=-fs;
					//bb<<-fs<<endl;
				}
				//bb.close();
			}
		}
		//体積力を計算する場合
		else if(MCON.get_force_t() == 1){
			if(MCON.get_eForce() == 0){
				cout<<"体積力:並進力ver.2"<<endl;
				for(int i = 0; i < fluid_number; i++){
					double f[3]={0.};
					cal_WLSM_for_grad(PART, i, CON.get_re3() * CON.get_distancebp(), M, H ,f);
					//calc_MPS(CON,HYPER,PART,i,M,H,f);
					for(int D = 0; D < 3; D++)F[D][i] = f[D];
				}
				//for(int i = 0; i < fluid_number; i++)for(int D = 0; D < d; D++)H[D][i] *= M[D][i];
				//H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
				//for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENSION;D++) F[D][i] = Hgrad[D][i];

				//単純な勾配計算だが，MPSよりもWLSMの方が圧力が安定する
				/*
				for(int i = 0; i < fluid_number; i++){
				double f[3];
				//WLSMによる勾配の計算
				//今回はre3を使う
				cal_WLSM_for_grad(PART, i, CON.get_re3(), M, H, f);
				for(int D=0;D<DIMENSION;D++)F[D][i] = f[D];
				}
				*/

			}else if(MCON.get_eForce() == 1){
				//*ケルビン力は妙な結果になるので，放置．
				cout<<"体積力:Kelvin力";
				/*H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
				for(int i = 0; i < fluid_number; i++){
				double Mn = sqrt(M[0][i] * M[0][i] + M[1][i] * M[1][i] + M[2][i] * M[2][i]);
				for(int D = 0; D < d; D++)Fv[D][i] = Mn * Hgrad[D][i];
				}

				H_gradient2(CON,PART, fluid_number,Hgrad,H[0]);//∇Hx計算
				H_gradient2(CON,PART, fluid_number,Hgrad1,H[1]);//∇Hy計算
				H_gradient2(CON,PART, fluid_number,Hgrad2,H[2]);//∇Hz計算
				for(int i = 0; i < fluid_number; i++){
				Fv[0][i] = M[0][i] * Hgrad[0][i] + M[1][i] * Hgrad[1][i] + M[2][i] * Hgrad[2][i];
				Fv[1][i] = M[0][i] * Hgrad1[0][i] + M[1][i] * Hgrad1[1][i] + M[2][i] * Hgrad1[2][i];
				Fv[2][i] = M[0][i] * Hgrad2[0][i] + M[1][i] * Hgrad2[1][i] + M[2][i] * Hgrad2[2][i];
				}

				//偶力
				/*
				cout<<"+偶力モーメント"<<endl;
				for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
				Gaiseki[0][i] = (M[1][i]*H[2][i] - M[2][i]*H[1][i]);
				Gaiseki[1][i] = (M[2][i]*H[0][i] - M[0][i]*H[2][i]);
				Gaiseki[2][i] = (M[0][i]*H[1][i] - M[1][i]*H[0][i]);
				}
				H_gradient2(CON,PART, fluid_number,Hgrad,Gaiseki[0]);//∇(M×H)x計算
				H_gradient2(CON,PART, fluid_number,Hgrad1,Gaiseki[1]);//∇(M×H)y計算
				H_gradient2(CON,PART, fluid_number,Hgrad2,Gaiseki[2]);//∇(M×H)z計算
				for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
				Fv[0][i] += .5*(Hgrad1[2][i] - Hgrad2[1][i]);
				Fv[1][i] += .5*(Hgrad2[0][i] - Hgrad[2][i]);
				Fv[2][i] += .5*(Hgrad[1][i] - Hgrad1[0][i]);
				}
				*/

			}
			//
			//for(int i = 0; i < fluid_number; i++){
			//	for(int D = 0; D < d; D++)F[D][i] = Fv[D][i] * CON.get_distancebp() * CON.get_distancebp() * CON.get_distancebp();
			//}
			cout<<endl;
		}
	}else{
		for(int i = 0; i < fluid_number; i++){
			for(int D = 0; D < DIMENSION; D++){
				PART[i].F[D] = 0;
			}
		}
	}




	if(t==1)
	{
		ofstream init0("M.dat",ios::trunc);
		ofstream init1("H.dat",ios::trunc);
		ofstream init2("F.dat",ios::trunc);
		ofstream init3("E_f.csv",ios::trunc);

		init0.close();
		init1.close();
		init2.close();
		init3.close();
	}

	ofstream fp("M.dat", ios::app);
	ofstream fq("H.dat", ios::app);
	ofstream ff("F.dat", ios::app);
	ofstream fe("E_f.csv",ios::app);

	ofstream ffal("F_all.dat");
	ofstream fn("N.dat");
	double times=1e-11;
	double sum_fe=0;
	double V=get_volume(&CON);
	fp<<"t="<<t<<endl;
	fq<<"t="<<t<<endl;
	ff<<"t="<<t<<endl;

	for(int i=0;i<fluid_number;i++)
	{
		fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]<<" "<<M[A_Y][i]<<" "<<M[A_Z][i]<<endl;
		fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
		if(PART[i].r[1] > -CON.get_distancebp() / 2. && PART[i].r[1] < CON.get_distancebp() / 2.0){
			ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]<<" "<<F[A_Z][i]<<endl;
		}
		//cout<<i<<"=F"<<F[A_X][i]<<","<<F[A_Y][i]<<","<<F[A_Z][i]<<", M"<<M[A_X][i]<<", "<<M[A_Y][i]<<","<<M[A_Z][i]<<", H"<<H[A_X][i]<<", "<<H[A_Y][i]<<","<<H[A_Z][i]<<endl;
		HYPER[i].Ef=V*(M[A_X][i]*H[A_X][i]+M[A_Y][i]*H[A_Y][i]+M[A_Z][i]*H[A_Z][i]);
		fe<<V*(M[A_X][i]*H[A_X][i]+M[A_Y][i]*H[A_Y][i]+M[A_Z][i]*H[A_Z][i])<<",";
		sum_fe+=V*(M[A_X][i]*H[A_X][i]+M[A_Y][i]*H[A_Y][i]+M[A_Z][i]*H[A_Z][i]);
		ffal<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
		fn<<PART[i].N<<endl;
	}

	fe<<sum_fe<<endl;
	ffal.close();
	fp.close();
	fq.close();
	ff.close();
	fe.close();
	fn.close();


	/*
	ofstream fp("M.dat", ios::app);
	ofstream fq("H.dat");
	ofstream ff("F.dat");
	double times=1e-11;
	if(d==2)
	{
	for(int i=0;i<fluid_number;i++)
	{
	fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<M[A_X][i]*times<<" "<<M[A_Y][i]*times<<endl;
	fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i])<<endl;
	ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<endl;
	}
	}
	else if(d==3)
	{
	for(int i=0;i<fluid_number;i++)
	{
	fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]<<" "<<M[A_Y][i]<<" "<<M[A_Z][i]<<endl;
	fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
	if(PART[i].r[1] > -CON.get_distancebp() / 2. && PART[i].r[1] < CON.get_distancebp() / 2.0){
	ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]<<" "<<F[A_Z][i]<<endl;
	}
	}
	}
	fp.close();
	fq.close();
	ff.close();*/

	////ｽﾑｰｼﾞﾝｸﾞ
	//smoothingF3D(CON,PART,fluid_number,F);

	ofstream fr("Fs.dat");
	times=5e-2/MCON.get_Hf_H();
	if(d==2){ for(int i=0;i<fluid_number;i++) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Y][i]*times<<endl;}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Z][i]*times<<endl;
	}
	fr.close();
	ofstream ft("Fv.dat");
	times=1e-6;
	if(d==2) for(int i=0;i<fluid_number;i++) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fv[A_X][i]*times<<" "<<Fv[A_Y][i]*times<<endl;
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fv[A_X][i]<<" "<<Fv[A_Z][i]<<endl;
	}
	ft.close();





	for(int D=0;D<DIMENSION;D++) delete [] M[D];
	for(int D=0;D<DIMENSION;D++) delete [] H[D];
	for(int D=0;D<DIMENSION;D++) delete [] direct[D];
	for(int D=0;D<DIMENSION;D++) delete [] Fs[D];
	for(int D=0;D<DIMENSION;D++) delete [] Fv[D];
	for(int D=0;D<DIMENSION;D++) delete [] Hgrad[D];
	delete  naiseki;
}

void Magnetic_Moment_Methodv2_mag2(mpsconfig &CON,vector<mpselastic> &PART,double **F,double n0,double lamda,int fluid_number,int particle_number, double current_time, int t)
{	
	double u0=4*PI*1e-7;
	MMM_config MCON;

	cout<<"磁気モーメント法---メッシュベース"<<endl;

	double le=CON.get_distancebp();
	double R=CON.get_re()*le;	//発散用影響半径
	int d=CON.get_dimension();						//次元
	unsigned int timeA=GetTickCount();				//計算開始時刻
	int count=0;
	int pn=fluid_number*d;							//未知数:粒子数×速度D(次元)成分 
	double co=1/(4*PI*u0);							//計算によく現れる係数
	double RP0=CON.get_RP();						//比透磁率
	double kai=(RP0-1.0);								//磁気感受率
	//double l=pow(CON.get_particle_volume(), double(1/3));
	double S=le;
	double mt = 1;                                //流体に対す磁性体の占める体積によって変化させる
	double ratio;
	//double mag_center[3]={0,0,-2.1e-05};

	int mesh_size;

	double Rm=CON.get_R1();

	//磁石１の寸法
	double size1[3] = {4*Rm, 4*Rm, 1.0*le};	//各辺の半分の寸法

	double mag1_center_Z_init =4*Rm+size1[2]; //磁石の重心の初期地
	double mag1_center_Z_goal =4*Rm+size1[2]; //磁石の重心の目標地．初期地と同じにすれば動かない．
	double mag1_center_dx = 0.2;        //磁石の移動速度(m/s)
	double mag1_center_Z_curr = mag1_center_Z_init + current_time * mag1_center_dx;         //磁石の現在地の更新
	double mag1_M[3] = {0,0,2.0};    //[T](交流の場合は振幅)


	//磁石２の寸法
	double size2[3] = {4*Rm, 4*Rm, 1.0*le};	//各辺の半分の寸法

	double mag2_center_Z_init =-1*(4*Rm+size2[2]); //磁石の重心の初期地
	double mag2_center_Z_goal =-1*(4*Rm+size2[2]); //磁石の重心の目標地．初期地と同じにすれば動かない．
	double mag2_center_dx = 0.2;        //磁石の移動速度(m/s)
	double mag2_center_Z_curr = mag2_center_Z_init + current_time * mag2_center_dx;         //磁石の現在地の更新
	double mag2_M[3] = {0,0,2.0};    //[T](交流の場合は振幅)





	/*//------交流磁場
	double mag_Hz = 7.0;
	//double duty = (mag_Hz * 0.15) / 7.25;
	double duty = 0.15;
	double ct = current_time;
	while(ct > (1.0 / mag_Hz)){
	ct -= (1.0 / mag_Hz);
	}
	if(ct > (1.0 / mag_Hz) * duty)mag_M[2] = 0;
	//impulse
	if(t > 500)mag_M[2]=0;else mag_M[2]=0.5;
	//------交流磁場F.dat
	*/



	if(mag1_center_Z_curr > mag1_center_Z_goal)mag1_center_Z_curr = mag1_center_Z_goal;
	if(mag2_center_Z_curr > mag2_center_Z_goal)mag2_center_Z_curr = mag2_center_Z_goal;

	double mag1_center[3] = {0,0,mag1_center_Z_curr};
	double mag2_center[3] = {0,0,mag2_center_Z_curr};

	cout<<"磁石1の位置={0,0,"<<mag1_center_Z_curr<<"}"<<endl;
	cout<<"磁石1の強さ={0,0,"<<mag1_M[2]<<"}"<<endl;
	cout<<"磁石2の位置={0,0,"<<mag2_center_Z_curr<<"}"<<endl;
	cout<<"磁石2の強さ={0,0,"<<mag2_M[2]<<"}"<<endl;


	if(d==3 && CON.get_model_set_way()==0) S=le*le;
	if(d==3 && CON.get_model_set_way()==1) S=sqrt(3.0)/4*le*le;//断面積
	double mesh_scale = 2.1;                                 //磁気モーメント法用メッシュの大きさ　磁石の寸法によっては変える必要あり　2
	double rm = 3;                                      //メッシュの中心から見込んだメッシュの数(1.9:26近傍,1:6近傍)．
	double width=CON.get_distancebp()*mesh_scale;		//格子幅
	int mesh_x = ((CON.get_maxX() - CON.get_minX()) / width + 0.5);
	int mesh_y = ((CON.get_maxY() - CON.get_minY()) / width + 0.5);
	int mesh_z = ((CON.get_maxZ() - CON.get_minZ()) / width + 0.5);
	mesh_size = mesh_x * mesh_y * mesh_z;

	vector<vector<double>>MESH;              //粒子の入っているメッシュの中心座標
	vector<int>mesh_index(fluid_number);     //粒子の入っているメッシュの番号(粒子の入っているメッシュのみでスケーリング)
	vector<int>mesh_p_num;	                 //メッシュ内に入っている粒子の数
	vector<vector<int>>p_index;              //メッシュ内にある粒子の番号

	//メッシュ内にどれだけの流体粒子があるか
	vector<int> mesh_i;
	vector<vector<int>> MESH_ALL(mesh_size, vector<int>(5));
	vector<vector<int>> MESH_N(mesh_size, vector<int>(0));
	vector<vector<int>> part_on_mesh(mesh_size, vector<int>(0));
	mesh_i.reserve(fluid_number);
	part_on_mesh.reserve(fluid_number);
	for(int i = 0; i < mesh_size; i++){
		MESH_ALL[i][3] = 0;
	}

	//電磁力関係
	double *Fs[DIMENSION];
	for(int D=0;D<DIMENSION;D++) Fs[D]=new double [particle_number];//単位面積あたりの力
	double *Fv[DIMENSION];
	for(int D=0;D<DIMENSION;D++) Fv[D]=new double [particle_number];//単位体積あたりの力
	double *Hgrad[DIMENSION];			//H勾配ﾍﾞｸﾄﾙ格納
	double *naiseki = new double[fluid_number];
	for(int D=0;D<DIMENSION;D++){
		Hgrad[D]=new double [fluid_number];
	}

	//磁化と磁界強度
	double *M[DIMENSION];					
	for(int D=0;D<DIMENSION;D++) M[D]=new double [fluid_number];
	double *H[DIMENSION];					
	for(int D=0;D<DIMENSION;D++) H[D]=new double [fluid_number];
	for(int i = 0; i < fluid_number; i++)for(int D = 0; D < DIMENSION; D++)H[D][i] = 0;



	double max_z = 0.;
	for(int i = 0; i < fluid_number; i++){
		if(max_z < PART[i].r[2]){
			if(max_z < PART[i].r[2] && sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1])<0.001){
				max_z = PART[i].r[2];
			}
		}
	}
	double z_max = 0;
	double z_min = 1;
	for(int i = 0; i < fluid_number; i++){
		if(z_max < PART[i].r[d-1]&&sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1]) < 2*le)z_max = PART[i].r[d-1];
		if(z_min > PART[i].r[d-1]&&sqrt(PART[i].r[0] * PART[i].r[0] + PART[i].r[1] * PART[i].r[1]) < 2*le)z_min = PART[i].r[d-1];
	}
	double x_max = 0;
	double x_min = 1;
	for(int i = 0; i < fluid_number; i++){
		if(x_max < PART[i].r[0])x_max = PART[i].r[0];
		if(x_min > PART[i].r[0])x_min = PART[i].r[0];
	}
	ofstream fh;
	cout<<"Diameter="<<x_max-x_min<<endl;
	cout<<"  Height="<<z_max-z_min<<endl;

	if(t==1){
		ofstream fh("./freq.dat");
		fh<<current_time<<" "<<z_max-z_min<<endl;
		fh.close();
	}else{
		ofstream fh("./freq.dat" , iostream::app);
		fh<<current_time<<" "<<z_max-z_min<<endl;
		fh.close();
	}

	//=======================================================================
	//	出力
	/////////////////////////////////////////////////////////////////////////
	if(t==1)
	{
		ofstream t1("height.dat");		//横軸時間、縦軸流体高さのグラフ
		t1.close();
	}
	ofstream t1("height.dat",ios :: app);
	t1<<current_time<<" "<<max_z-0.005<<endl;
	t1.close();

	//////////////////////////////////////////////////////////////////////////

	//バックグラウンドセルの番号を設定
	for(int i = 0; i < mesh_x; i++){
		for(int j = 0; j < mesh_y; j++){
			for(int k = 0; k < mesh_z; k++){
				MESH_ALL[k * mesh_x * mesh_y + j * mesh_x + i][0] = i;
				MESH_ALL[k * mesh_x * mesh_y + j * mesh_x + i][1] = j;
				MESH_ALL[k * mesh_x * mesh_y + j * mesh_x + i][2] = k;
			}
		}
	}

	//流体粒子がどのセルに属しているのか計算
	for(int i=0;i<fluid_number;i++)
	{
		int X=(int)((PART[i].r[A_X]-CON.get_minX())/width + 0.000001);//X方向に何個目の格子か 
		int Y=(int)((PART[i].r[A_Y]-CON.get_minY())/width + 0.000001);//Y方向に何個目の格子か
		int Z=(int)((PART[i].r[A_Z]- le -CON.get_minZ())/width + 0.000001);//Z方向に何個目の格子か
		int number=Z*mesh_x*mesh_y+Y*mesh_x+X;//粒子iを含む格子の番号
		MESH_ALL[number][3]++;
		mesh_i.emplace_back(number);
		part_on_mesh[number].emplace_back(i);
	}

	cout<<"メッシュの排除"<<endl;
	//計算の都合上，いらないメッシュを排除
	//ここでいういらないメッシュとは，流体が入っていないセルのこと
	for(int i = 0; i < mesh_size; i++){
		vector<double>work(3);
		vector<int>p_work;
		//セルiの中の粒子の数が0でないか検査
		if(MESH_ALL[i][3] != 0){
			work[0] = (((double)MESH_ALL[i][0] + 0.5) * width) + CON.get_minX();
			work[1] = (((double)MESH_ALL[i][1] + 0.5) * width) + CON.get_minY();
			work[2] = (((double)MESH_ALL[i][2]) * width) + CON.get_minZ() + le;
			//セルiの情報を計算に関与するセル情報に追加
			MESH.push_back(work);
			//追加したセルの中に入っている粒子の数をコピー
			mesh_p_num.push_back(MESH_ALL[i][3]);
			//流体粒子の入っているセルを改めてナンバリング
			for(int j = 0; j < MESH_ALL[i][3]; j++){
				p_work.push_back(part_on_mesh[i][j]);
				mesh_index[part_on_mesh[i][j]] = i;
			}
			//ナンバリングしたセルと粒子の番号を関連づけ
			p_index.push_back(p_work);
		}
	}

	cout<<"隣接メッシュ"<<endl;
	//あるメッシュに対して，隣接するメッシュがいくつあるのか計算する
	for(int i = 0; i < MESH.size(); i++){
		for(int j = 0; j < MESH.size(); j++){
			if(i == j)continue; //自分自身は計算に入れない
			double xm = (MESH[i][0] - MESH[j][0]);
			double ym = (MESH[i][1] - MESH[j][1]);
			double zm = (MESH[i][2] - MESH[j][2]);
			double r = sqrt(xm * xm + ym * ym + zm * zm);
			if(r <= rm * width){
				MESH_N[i].push_back(j);
			}
		}
	}


	cout<<"対象となるメッシュ:"<<MESH.size()<<endl;

	//行列確保 //未知数の順番はMx0,My0,Mz0,Mx1,My1,Mz1,Mx2,My2,Mz2・・・
	//メッシュごとの磁化を求めようとしている．
	vector<vector<double> > matrix(MESH.size() * 3, vector<double>(MESH.size() * 3));
	vector<double> B(MESH.size() * 3);
	double *direct[DIMENSION];
	for(int D=0;D<DIMENSION;D++) direct[D]=new double [particle_number];//外向き法線ベクトル

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++) direct[D][i]=0;				//初期化
		if(PART[i].surface==ON)
		{
			direct_f(CON,PART,i,direct);
			for(int D=0;D<d;D++) direct[D][i]*=-1;//外向きが欲しいから反転する
		}
	}

	//この段階では，磁気モーメント法での要素数はナンバリングされたセルの数になっている
	for(int i=0;i<MESH.size();i++)
	{
		B[i]=0;
		for(int j=0;j<MESH.size();j++) matrix[i][j]=0;		//初期化
	}
	//////////////////////////////////////////////////////////////////

	//解行列B[]作成

	//各セル内で磁化を求める．
	cout<<"解行列の作成-----"<<endl;
	if(MCON.get_Hf_type()==0)
	{
		cout<<"\t一様磁場---"<<endl;
		double sign[3]={0,0,0};
		sign[d-1]=1;
		for(int i=0;i<MESH.size();i++)
		{
			for(int D=0;D<d;D++) B[i*d+D]=MCON.get_Hf_H()*sign[D] * (RP0 - 1) * MU_0;//2次元ならA_Y,3次元ならA_Zの方向だけ値が入る。
		}
		cout<<"OK"<<endl;
		cout<<"\t飽和磁化の考慮---";
		for(int i=0;i<MESH.size();i++)
		{
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(B[i * d + D],2.0);
			}
			nB = sqrt(nB);
			double ratio;
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
			}
		}
		cout<<"OK"<<endl;
	}
	else if(MCON.get_Hf_type()==1)
	{
		cout<<"\t磁石からの影響計算---";
		for(int i=0;i<MESH.size();i++)
		{
			//磁石１
			double ki1[3][3] = {.0};
			double p1[3] = {MESH[i][0] - mag1_center[0], MESH[i][1] - mag1_center[1], MESH[i][2] - mag1_center[2]};
			Cal_Mag_Density(ki1, p1, size1[0], size1[1], size1[2]);
			for(int j = 0; j < d; j++){
				for(int k = 0; k < d; k++){
					B[i * d + j] += ki1[k][j] * mag1_M[k] / (-4.0 * PI) * (RP0 - 1);
				}
			}

			//磁石２
			double ki2[3][3] = {.0};
			double p2[3] = {MESH[i][0] - mag2_center[0], MESH[i][1] - mag2_center[1], MESH[i][2] - mag2_center[2]};
			Cal_Mag_Density(ki2, p2, size2[0], size2[1], size2[2]);
			for(int j = 0; j < d; j++){
				for(int k = 0; k < d; k++){
					B[i * d + j] += ki2[k][j] * mag2_M[k] / (-4.0 * PI) * (RP0 - 1);
				}
			}

		}
		cout<<"OK"<<endl;

		cout<<"\t飽和磁化の考慮---";
		for(int i=0;i<MESH.size();i++)
		{
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(B[i * d + D],2.0);
			}
			nB = sqrt(nB);
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				B[i * d + D] *= ratio;
			}
		}
		cout<<"OK"<<endl;
	}else cout<<"外部磁場のタイプが未実装"<<endl;
	cout<<"-----完了"<<endl;
	/////
	if(mag1_M[2]!=0 || mag2_M[2]!=0)
	{
		/////////////////////matrixに値を格納
		for(int dst=0;dst<MESH.size();dst++)
		{
			for(int src=0;src<MESH.size();src++)
			{
				double px = MESH[dst][0] - MESH[src][0];
				double py = MESH[dst][1] - MESH[src][1];
				double pz = MESH[dst][2] - MESH[src][2];
				double ki[3][3] = {.0};
				double rA[3]={px, py, pz};
				//Cal_Mag_Density(ki, rA, (width / 2.), (width / 2.), 0);
				Cal_Mag_Density(ki, rA, (width / 2.));
				for(int i = 0; i < d; i++){
					for(int j = 0; j < d; j++){
						if(dst == src && i == j){
							matrix[dst * d + i][src * d + j] += (ki[i][j] + 4 * PI) / (4 * PI) * (RP0 - 1) + 1;
						}else matrix[dst * d + i][src * d + j] += ki[i][j] / (4 * PI) * (RP0 - 1);
					}
				}
			}
		}
		cout<<"行列作成"<<endl;

		//行列値出力
		//output_matrix(matrix,B, pn);

		cout<<"未知数:"<<MESH.size() * 3<<".ガウスザイデル法--";
		unsigned int timeB=GetTickCount();
		//行列を解く（ガウスザイデル法）
		GaussSeidel(matrix, MESH.size() * 3, B, d);
		//jacobi(matrix,B,pn);//答えはBに格納
		cout<<"OK time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

		//粒子位置での磁場H
		if(MCON.get_isWLSM() == 0){
			cout<<"重みつき最小二乗法による推定"<<endl;
			for(int i=0;i<MESH.size();i++)
			{
				//関数近似後の各係数
				double cx[9] = {.0};
				double cy[9] = {.0};
				double cz[9] = {.0};
				//ある集合体内で最小二乗法を行う
				cal_WLSM_for_M2(MESH, MESH_N, B, i, width, width * rm, cx, cy, cz);
				for(int j = 0; j <p_index[i].size(); j++){
					double Mb[3] = {0};
					//係数から変化量を求める．この中で関数を再現して同一セル内にある粒子の座標に応じた磁化を付与する．
					get_M2(cx, cy, cz, PART[p_index[i][j]].r[A_X], PART[p_index[i][j]].r[A_Y], PART[p_index[i][j]].r[A_Z], MESH[i][0], MESH[i][1], MESH[i][2], Mb);
					for(int D=0;D<d;D++)M[D][p_index[i][j]]=Mb[D] + B[i * 3 + D];
				}
			}
		}
		//最小二乗法を用いた磁化補間ができた

		cout<<"OK"<<endl;
		//一応飽和磁化の考慮をしておく．なぜなら，近似関数の関係で飽和磁化より大きい値を持つ可能性があるから．
		cout<<"飽和磁化の考慮---";
		for(int i = 0; i < fluid_number; i++){
			double nB = 0.0;
			for(int D=0;D<d;D++){
				nB += pow(M[D][i],2.0);
			}
			nB = sqrt(nB);
			nB > mt ? ratio = mt / nB : ratio = 1.0;
			for(int D=0;D<d;D++){
				M[D][i] *= ratio;
			}
		}
		cout<<"OK"<<endl;


		cout<<"外部磁場計算---";
		//磁石による磁性流体に対する外部磁場の計算
		if(MCON.get_Hf_type()!=0){
			cout<<"磁石"<<endl;
			for(int i=0;i<fluid_number;i++)
			{
				double ki1[3][3] = {.0};
				double p1[3] = {PART[i].r[A_X]-mag1_center[A_X], PART[i].r[A_Y]-mag1_center[A_Y], PART[i].r[A_Z]-mag1_center[A_Z]};
				Cal_Mag_Density(ki1, p1, size1[0], size1[1], size1[2]);

				double ki2[3][3] = {.0};
				double p2[3] = {PART[i].r[A_X]-mag2_center[A_X], PART[i].r[A_Y]-mag2_center[A_Y], PART[i].r[A_Z]-mag2_center[A_Z]};
				Cal_Mag_Density(ki2, p2, size2[0], size2[1], size2[2]);

				for(int j = 0; j < d; j++){
					for(int k = 0; k < d; k++){
						H[j][i] = (ki1[j][k] * mag1_M[k]+ki2[j][k] * mag2_M[k]) / (-16. * PI * PI) * 10000000;
					}
				}
			}
		}

		/*
		//本来は，全粒子からの影響を計算すべきではあるものの，計算量が膨大となるので，
		//ここでは，影響半径内に入った粒子からの磁場の寄与のみを考える．
		for(int dst=0;dst<fluid_number;dst++)
		{
		for(int src=0;src<PART[dst].N;src++)
		{
		if(dst == src)continue; //自己場はいらない．
		double k1[3][3] = {.0};
		double rA[3]={PART[dst].r[A_X]-PART[PART[dst].NEI[src]].r[A_X],PART[dst].r[A_Y]-PART[PART[dst].NEI[src]].r[A_Y],PART[dst].r[A_Z]-PART[PART[dst].NEI[src]].r[A_Z]};
		double r = sqrt(rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2]);
		double l = CON.get_distancebp();
		Cal_Mag_Density(k1, rA, l / 2.);
		for(int i = 0; i < d; i++){
		for(int j = 0; j < d; j++){
		H[i][dst] += k1[i][j] * M[j][src] / (-16. * PI * PI) * 10000000;
		}
		}
		}
		}*/
		cout<<"OK"<<endl;

		//応力を計算する場合
		if(MCON.get_force_t() == 0){
			cout<<"応力"<<endl;
			for(int i=0;i<fluid_number;i++)
			{
				for(int D=0;D<DIMENSION;D++) Fs[D][i]=0;		//初期化
				if(PART[i].surface==ON)
				{
					for(int D=0;D<DIMENSION;D++)H[D][i] = M[D][i] / (u0 * (RP0 - 1.));
					//cout<<"a"<<endl;
					//double Mn=(H[A_X][i]*direct[A_X][i]+H[A_Y][i]*direct[A_Y][i]+H[A_Z][i]*direct[A_Z][i]);
					//double Mn=(M[A_X][i]*direct[A_X][i]+M[A_Y][i]*direct[A_Y][i]+M[A_Z][i]*direct[A_Z][i]) / (RP0 - 1.);
					//double val=0.5 * Mn * Mn / u0;		//応力値
					//double val = 0.5 * Mn * Mn * u0;		//応力値
					//for(int D=0;D<DIMENSION;D++) Fs[D][i]=val*direct[D][i];
					//for(int D=0;D<DIMENSION;D++) F[D][i]=Fs[D][i]*CON.get_distancebp();
					//一項目
					double H_n = H[0][i] * direct[0][i] + H[1][i] * direct[1][i] + H[2][i] * direct[2][i];
					//double Mn = M[0][i] * direct[0][i] + M[1][i]
					//double Hr = H[0][i] * H[0][i] + H[1][i] * H[1][i] + H[2][i] * H[2][i];
					//double H_n = M[0][i] / u0 * 2 * direct[0][i] + M[1][i]  / u0 * 2 * direct[1][i] + M[2][i] / u0 * 2 * direct[2][i];			
					for(int D=0;D<DIMENSION;D++) Fs[D][i] = 0.5 * H_n * H_n * u0 * direct[D][i];
					//二項目
					//double H_H = H[0][i] * H[0][i] + H[1][i] * H[1][i] + H[2][i] * H[2][i];
					//合計
					//for(int D=0;D<DIMENSION;D++) Fs[D][i] -= (.5 * u0 * (H_H * direct[D][i]));
					//for(int D=0;D<DIMENSION;D++) F[D][i] = Fs[D][i] * CON.get_distancebp() * CON . get_distancebp();
				}
			}
			if(CON.get_dir_for_P()==2 ||CON.get_dir_for_P()==3 )
			{
				//ofstream bb("electromagnetic_P.dat");
				for(int i=0;i<fluid_number;i++)
				{
					double fs=0;//表面力
					if(PART[i].surface==ON)
					{
						fs=sqrt(Fs[A_X][i]*Fs[A_X][i]+Fs[A_Y][i]*Fs[A_Y][i]+Fs[A_Z][i]*Fs[A_Z][i]);
						for(int D=0;D<DIMENSION;D++) F[D][i]=0;
					}
					PART[i].dir_Pem=-fs;
					//bb<<-fs<<endl;
				}
				//bb.close();
			}
		}
		//体積力を計算する場合
		else if(MCON.get_force_t() == 1){
			if(MCON.get_eForce() == 0){
				cout<<"体積力:並進力ver.2"<<endl;
				for(int i = 0; i < fluid_number; i++){
					double f[3]={0};
					cal_WLSM_for_grad(PART, i, CON.get_re3() * CON.get_distancebp(), M, H ,f);
					for(int D = 0; D < 3; D++)F[D][i] = f[D];
				}
				//for(int i = 0; i < fluid_number; i++)for(int D = 0; D < d; D++)H[D][i] *= M[D][i];
				//H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
				//for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENSION;D++) F[D][i] = Hgrad[D][i];

				//単純な勾配計算だが，MPSよりもWLSMの方が圧力が安定する
				/*
				for(int i = 0; i < fluid_number; i++){
				double f[3];
				//WLSMによる勾配の計算
				//今回はre3を使う
				cal_WLSM_for_grad(PART, i, CON.get_re3(), M, H, f);
				for(int D=0;D<DIMENSION;D++)F[D][i] = f[D];
				}
				*/

			}else if(MCON.get_eForce() == 1){
				//*ケルビン力は妙な結果になるので，放置．
				cout<<"体積力:Kelvin力";
				/*H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
				for(int i = 0; i < fluid_number; i++){
				double Mn = sqrt(M[0][i] * M[0][i] + M[1][i] * M[1][i] + M[2][i] * M[2][i]);
				for(int D = 0; D < d; D++)Fv[D][i] = Mn * Hgrad[D][i];
				}

				H_gradient2(CON,PART, fluid_number,Hgrad,H[0]);//∇Hx計算
				H_gradient2(CON,PART, fluid_number,Hgrad1,H[1]);//∇Hy計算
				H_gradient2(CON,PART, fluid_number,Hgrad2,H[2]);//∇Hz計算
				for(int i = 0; i < fluid_number; i++){
				Fv[0][i] = M[0][i] * Hgrad[0][i] + M[1][i] * Hgrad[1][i] + M[2][i] * Hgrad[2][i];
				Fv[1][i] = M[0][i] * Hgrad1[0][i] + M[1][i] * Hgrad1[1][i] + M[2][i] * Hgrad1[2][i];
				Fv[2][i] = M[0][i] * Hgrad2[0][i] + M[1][i] * Hgrad2[1][i] + M[2][i] * Hgrad2[2][i];
				}

				//偶力
				/*
				cout<<"+偶力モーメント"<<endl;
				for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
				Gaiseki[0][i] = (M[1][i]*H[2][i] - M[2][i]*H[1][i]);
				Gaiseki[1][i] = (M[2][i]*H[0][i] - M[0][i]*H[2][i]);
				Gaiseki[2][i] = (M[0][i]*H[1][i] - M[1][i]*H[0][i]);
				}
				H_gradient2(CON,PART, fluid_number,Hgrad,Gaiseki[0]);//∇(M×H)x計算
				H_gradient2(CON,PART, fluid_number,Hgrad1,Gaiseki[1]);//∇(M×H)y計算
				H_gradient2(CON,PART, fluid_number,Hgrad2,Gaiseki[2]);//∇(M×H)z計算
				for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++){
				Fv[0][i] += .5*(Hgrad1[2][i] - Hgrad2[1][i]);
				Fv[1][i] += .5*(Hgrad2[0][i] - Hgrad[2][i]);
				Fv[2][i] += .5*(Hgrad[1][i] - Hgrad1[0][i]);
				}
				*/

			}
			//
			//for(int i = 0; i < fluid_number; i++){
			//	for(int D = 0; D < d; D++)F[D][i] = Fv[D][i] * CON.get_distancebp() * CON.get_distancebp() * CON.get_distancebp();
			//}
			cout<<endl;
		}
	}
	else
	{
		for(int i = 0; i < fluid_number; i++){
			for(int D = 0; D < DIMENSION; D++){
				PART[i].F[D] = 0;
			}
		}
	}




	if(t==1)
	{
		ofstream init0("M.dat",ios::trunc);
		ofstream init1("H.dat",ios::trunc);
		ofstream init2("F.dat",ios::trunc);
		ofstream init3("E_f.csv",ios::trunc);

		init0.close();
		init1.close();
		init2.close();
		init3.close();
	}

	ofstream fp("M.dat", ios::app);
	ofstream fq("H.dat", ios::app);
	ofstream ff("F.dat", ios::app);
	ofstream fe("E_f.csv",ios::app);

	double times=1e-11;
	double sum_fe=0;
	double V=get_volume(&CON);
	fp<<"t="<<t<<endl;
	fq<<"t="<<t<<endl;
	ff<<"t="<<t<<endl;

	for(int i=0;i<fluid_number;i++)
	{
		fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]<<" "<<M[A_Y][i]<<" "<<M[A_Z][i]<<endl;
		fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
		if(PART[i].r[1] > -CON.get_distancebp() / 2. && PART[i].r[1] < CON.get_distancebp() / 2.0){
			ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]<<" "<<F[A_Z][i]<<endl;
		}
		fe<<V*(M[A_X][i]*H[A_X][i]+M[A_Y][i]*H[A_Y][i]+M[A_Z][i]*H[A_Z][i])<<",";
		sum_fe+=V*(M[A_X][i]*H[A_X][i]+M[A_Y][i]*H[A_Y][i]+M[A_Z][i]*H[A_Z][i]);
	}
	fe<<sum_fe<<endl;

	fp.close();
	fq.close();
	ff.close();
	fe.close();
	/*
	ofstream fp("M.dat", ios::app);
	ofstream fq("H.dat");
	ofstream ff("F.dat");
	double times=1e-11;
	if(d==2)
	{
	for(int i=0;i<fluid_number;i++)
	{
	fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<M[A_X][i]*times<<" "<<M[A_Y][i]*times<<endl;
	fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i])<<endl;
	ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<endl;
	}
	}
	else if(d==3)
	{
	for(int i=0;i<fluid_number;i++)
	{
	fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]<<" "<<M[A_Y][i]<<" "<<M[A_Z][i]<<endl;
	fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
	if(PART[i].r[1] > -CON.get_distancebp() / 2. && PART[i].r[1] < CON.get_distancebp() / 2.0){
	ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]<<" "<<F[A_Z][i]<<endl;
	}
	}
	}
	fp.close();
	fq.close();
	ff.close();*/

	////ｽﾑｰｼﾞﾝｸﾞ
	//smoothingF3D(CON,PART,fluid_number,F);

	ofstream fr("Fs.dat");
	times=5e-2/MCON.get_Hf_H();
	if(d==2){ for(int i=0;i<fluid_number;i++) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Y][i]*times<<endl;}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Z][i]*times<<endl;
	}
	fr.close();
	ofstream ft("Fv.dat");
	times=1e-6;
	if(d==2) for(int i=0;i<fluid_number;i++) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fv[A_X][i]*times<<" "<<Fv[A_Y][i]*times<<endl;
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fv[A_X][i]<<" "<<Fv[A_Z][i]<<endl;
	}
	ft.close();





	for(int D=0;D<DIMENSION;D++) delete [] M[D];
	for(int D=0;D<DIMENSION;D++) delete [] H[D];
	for(int D=0;D<DIMENSION;D++) delete [] direct[D];
	for(int D=0;D<DIMENSION;D++) delete [] Fs[D];
	for(int D=0;D<DIMENSION;D++) delete [] Fv[D];
	for(int D=0;D<DIMENSION;D++) delete [] Hgrad[D];
	delete  naiseki;
}



//2階テンソル用
void H_gradient2(mpsconfig &CON,vector<mpsparticle> &PART,int fluid_number,double **Hgrad,double *H)
{
	double le=CON.get_distancebp();//初期粒子間距離
	double r=CON.get_re()*le;
	int d=CON.get_dimension();

	double *HH=new double[fluid_number];//各粒子位置での磁場強度H格納


	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++) Hgrad[D][i]=0;//初期化

		double W=0;//粒子数密度　OUTを除いたりするのでPND[i]は微妙

		for(int k=0;k<PART[i].N;k++)
		{       
			int j=PART[i].NEI[k];
			if(PART[j].type==FLUID || PART[j].type==HYPERELAST)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);

				double w=kernel(r,dis);
				W+=w;

				Hgrad[A_X][i]+=(H[j]-H[i])*X*w/(dis*dis);
				Hgrad[A_Y][i]+=(H[j]-H[i])*Y*w/(dis*dis);
				Hgrad[A_Z][i]+=(H[j]-H[i])*Z*w/(dis*dis);
			}

		}
		for(int D=0;D<DIMENSION;D++) if(W!=0) Hgrad[D][i]= Hgrad[D][i]*d/W;
	}///////////////Pgrad[D][i]計算終了



	delete [] HH;
}








//行列値出力関数(管理用)
void output_matrix(double **matrix,double *Bmatrix,int node_num)
{
	//node_num=30;					//左上だけ表示したい時
	int *DDN=new int[node_num];		//対角優位　diagonally dominant matrix
	ofstream fp2("matrixMMM.dat");
	for(int n=0;n<node_num;n++)
	{
		double val=0;
		DDN[n]=OFF;
		for(int m=0;m<node_num;m++)
		{
			fp2<<matrix[n][m]<<"\t";
			if(n!=m) val+=sqrt(matrix[n][m]*matrix[n][m]);
		}
		fp2<<endl;
		if(val<sqrt(matrix[n][n]*matrix[n][n])) DDN[n]=ON;//対角優位
	}
	fp2.close();
	ofstream fp4("BmatrixMMM.dat");
	for(int n=0;n<node_num;n++) fp4<<Bmatrix[n]<<endl;
	fp4.close();

	//for(int n=0;n<node_num;n++) if(DDN[n]==ON) cout<<n<<"は対角優位"<<endl;

	delete [] DDN;
}

//jacobiの反復法 解は最終的にBのなかへ
void jacobi(double **matrix,double *B,int N)
{
	double ep=1e-8;//収束判定
	double E=10;		//誤差
	int count=0;

	double *X=new double[N];//解
	for(int k=0;k<N;k++) X[k]=0;//初期値
	while(E>ep)
	{
		E=0;
		for(int i=0;i<N;i++)
		{
			double L=0;
			double U=0;
			for(int j=0;j<i;j++) L+=matrix[i][j]*X[j]; 
			for(int j=i+1;j<N;j++) U+=matrix[i][j]*X[j];
			double Xnew=(B[i]-L-U)/matrix[i][i];
			E+=fabs(Xnew-X[i]);
			X[i]=Xnew;
		}
		E/=N;
		count++;
		cout<<count<<" E="<<E<<endl;
	}

	for(int k=0;k<N;k++) B[k]=X[k];//Bに答えを格納
	delete [] X;
}