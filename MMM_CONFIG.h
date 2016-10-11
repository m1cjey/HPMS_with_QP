////新しい粒子を追加したときいじらないといけないのは①u_laplacian_fのﾉﾝｽﾘｯﾌﾟ②P_gradient③P_gradient2④calc_Temperatureのk[i]
////⑤ freeon
#ifndef mmmconfig
#define mmmconfig
class MMM_config
{       
	///////解析条件

	
	
	int Hf_type;	//外部磁場　0:一様磁場 1:磁石
	double Hf_H;			//外部磁場Hの強さ
	int isWLSM;             //重みつき最小二乗法=0, 算術平均=1
	int eForce;             //電磁力計算法(並進力=0，kelvin力=1)
	int force_t;			//電磁力タイプ(応力=0, 体積力=1)

	int P_interval;         //圧力コンター出力間隔

	
	public:
	MMM_config();
	int		get_Hf_type()	{return Hf_type;}
	int		get_Hf_H()		{return Hf_H;}
	int		get_isWLSM()	{return isWLSM;}
	int		get_eForce()	{return eForce;}
	int		get_force_t()	{return force_t;}

	int		get_P_interval(){return P_interval;}

};


#endif
