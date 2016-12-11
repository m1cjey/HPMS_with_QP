#include "stdafx.h"	

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num,double V, double mi, double Dt, double nG[DIMENSION], double E0);
void p_QP(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num,double V, double mi, double Dt, double nG[DIMENSION], double E0);
void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,double nG[DIMENSION]);
void q_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT,double *rL,double *dL,int h_num,int count,int count_min,int t,double E);

//Hyperelstic2からの引継ぎ
void calc_W(mpsconfig &CON,vector<hyperelastic> &HYPER,int h_num);
void calc_pi(mpsconfig &CON,vector<hyperelastic> &HYPER);
void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void previous_strage(vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num);

void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F);
void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1, int t);
void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t);


void calc_HYPER_QP_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{

	ofstream time("time_log.dat",ios::app);
	clock_t	start_t=clock();
	time<<start_t*CLOCKS_PER_SEC<<"	";

	////////////定義///////////////
	int h_num=HYPER.size();

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();

	cout<<"h_num="<<h_num<<endl;
	cout<<"Hypercalculation starts."<<endl;


	///////////初期計算	派生クラスとして計算できたら理想
	if(t==1)
	{
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=PART[i].r[D];
		calc_constant(CON,PART,HYPER,HYPER1);
		calc_pi(CON,HYPER);
		calc_W(CON,HYPER,h_num);
	}
	
	///////////データ出力
	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,t);
		momentum_movie_AVS(CON,t,PART,HYPER,F);
		output_energy(CON,PART,HYPER,t);
	}

	///////////前ステップ計算結果保管
	previous_strage(PART,HYPER,HYPER1,h_num);





	double nG[DIMENSION]={0,0,1};

	double E0=0;

	cout<<"QP start-------------------------"<<endl;
	for(int i=0;i<h_num;i++)
	{
		E0+=0.5/mi*(HYPER[i].p_n[A_X]*HYPER[i].p_n[A_X]+HYPER[i].p_n[A_Y]*HYPER[i].p_n[A_Y]+HYPER[i].p_n[A_Z]*HYPER[i].p_n[A_Z])+V*HYPER[i].W_n;
	}
	cout<<"E0="<<E0<<endl<<endl;
	ofstream fq0("q0_QP.csv");
	ofstream fqn0("qn0_QP.csv");
	ofstream fpn0("pn0_QP.csv");
	ofstream fphn0("phn0_QP.csv");
	for(int i=0;i<h_num;i++)
	{
		fq0<<HYPER[i].q_n[A_X]<<","<<HYPER[i].q_n[A_Y]<<","<<HYPER[i].q_n[A_Z]<<endl;
		fqn0<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		fpn0<<HYPER[i].p_n[A_X]<<","<<HYPER[i].p_n[A_Y]<<","<<HYPER[i].p_n[A_Z]<<endl;
		fphn0<<HYPER[i].ph_n[A_X]<<","<<HYPER[i].ph_n[A_Y]<<","<<HYPER[i].ph_n[A_Z]<<endl;
	}
	fq0.close();
	fqn0.close();
	fpn0.close();
	fphn0.close();

	////////////QP計算///////////////		
	q_QP(CON,PART,HYPER,HYPER1,h_num,V,mi,Dt,nG,E0);

	ofstream fqn("qn_QP.csv");
	ofstream fphn("phn_QP.csv");
	for(int i=0;i<h_num;i++)
	{
		fqn<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		fphn<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
	}
	fqn.close();
	fphn.close();


	cout<<"--------------------------OK"<<endl;






	cout<<"Hypercalculation ends."<<endl;

	clock_t end_t=clock();
	time<<end_t*CLOCKS_PER_SEC<<"	";
	time.close();

}




void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num, double V, double mi, double Dt, double nG[DIMENSION], double E0)
{
	////////////定義///////////////
	int count=0;
	int count_min=0;
	int c_max=10000;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1e-5;
	double ep_min=1e-5;

	double r=1.0E10;

	double En=0;
	double T=0;

	double *dE=new double [h_num];
	double *dT=new double [h_num];	

	double *rE=new double [h_num*h_num];
	double *rT=new double [h_num*h_num];

	double *g=new double [h_num];
	double **dg=new double *[h_num];

	for(int i=0;i<h_num;i++)
	{
		dg[i]=new double [h_num];
	}

	double *th_g=new double [h_num];
		
	////////////初期化算///////////////
	for(int i=0;i<h_num;i++)
	{
		dE[i]=0;
		dT[i]=0;
		g[i]=0;
		th_g[i]=0;
		for(int j=0;j<h_num;j++)
		{
			dg[i][j]=0;
			rE[i*h_num+j]=0;
			rT[i*h_num+j]=0;
		}
	}


	while(E_min>ep_min)
	{
		count_min++;
	
		for(int i=0;i<h_num;i++)
		{
			HYPER[i].old_lam=HYPER[i].lam;
		}
		
		E=1;
		count=0;
		while(E>ep)
		{
			count++;

			q_variables(CON,PART,HYPER,HYPER1);
			q_nab_lap(PART,HYPER,HYPER1,dE,rE,dg,Dt,V,mi,nG);

			En=0;
			for(int i=0;i<h_num;i++)	En+=0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])+V*HYPER[i].W;

			T=(En-E0)*(En-E0);
			//cout<<"T="<<T<<", En"<<En<<endl;
			for(int k=0;k<h_num;k++)
			{	
				dT[k]=2.*dE[k]*(E-E0);
				for(int l=0;l<h_num;l++)	rT[l*h_num+k]=2.*rE[l*h_num+k]*(E-E0)+2.*dE[k]*dE[l];
			}
			for(int k=0;k<h_num;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			for(int i=0;i<h_num;i++)
			{	
				g[i]=V*(1-HYPER[i].J);
				
				T+=0.5*r*(g[i]*g[i]+th_g[i])*(g[i]*g[i]+th_g[i]);
			}

			for(int k=0;k<h_num;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					dT[k]+=2.*r*dg[i][k]*g[i]*(g[i]*g[i]+th_g[i]);
					//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
				}
				for(int l=0;l<h_num;l++)
				{
					for(int i=0;i<h_num;i++)
					{	
						rT[l*h_num+k]+=2.*r*g[i]*dg[i][k]*dg[i][l]*(3.*g[i]*g[i]+th_g[i]);
					}
				}
			}
			for(int k=0;k<h_num;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			//if(count==1||count%100==0)
			{

				stringstream sg;
				sg<<"./g/g"<<count_min<<"_"<<count<<".csv";
				ofstream fg(sg.str());
				stringstream sdg;
				sdg<<"./dg/dg"<<count_min<<"_"<<count<<".csv";
				ofstream fdg(sdg.str());

				for(int i=0;i<h_num;i++)
				{
					fg<<g[i]<<","<<g[i]*g[i]<<","<<th_g[i]<<endl;
					for(int j=0;j<h_num;j++)
					{
						fdg<<dg[i][j]<<",";
					}
					fdg<<endl;
				}
				fdg.close();

				stringstream std;
				std<<"./dT/dT"<<count_min<<"_"<<count<<".csv";
				ofstream fd(std.str());
				stringstream sTr;
				sTr<<"./rT/rT"<<count_min<<"_"<<count<<".csv";
				ofstream fr(sTr.str());

				stringstream sed;
				sed<<"./dE/dE"<<count_min<<"_"<<count<<".csv";
				ofstream fde(sed.str());
				stringstream str;
				str<<"./rE/rE"<<count_min<<"_"<<count<<".csv";
				ofstream fre(str.str());

				for(int i=0;i<h_num;i++)
				{
					fd<<dT[i]<<endl;
					fde<<dE[i]<<endl;
					for(int j=0;j<h_num;j++)
					{
						fr<<rT[i*h_num+j]<<",";
						fre<<rE[i*h_num+j]<<",";
					}
					fr<<endl;
					fre<<endl;
				}
				fd.close();
				fr.close();
				fde.close();
				fre.close();

			}
			gauss(rT,dT,h_num);
			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam+=-dT[i];
				E_sum+=dT[i]*dT[i];
			}
			E=sqrt(E_sum);

			/////////////pn1_2計算
			cout<<"E"<<count<<"="<<E<<" ,En="<<En<<endl;
			
			//if(count==1||count%100==0)
			{
				stringstream stpnh;
				stpnh<<"./p/phn_QP"<<count_min<<"_"<<count<<".csv";
				ofstream fp(stpnh.str());
				stringstream stq;
				stq<<"./q/qn_QP"<<count_min<<"_"<<count<<".csv";
				ofstream fq(stq.str());
				stringstream sd;
				sd<<"./lam/d"<<count_min<<"_"<<count<<".csv";
				ofstream fd(sd.str());
				stringstream sj;
				sj<<"./lam/J"<<count_min<<"_"<<count<<".csv";
				ofstream fj(sj.str());
				stringstream stl;
				stl<<"./lam/lam"<<count_min<<"_"<<count<<".csv";
				ofstream fl(stl.str());
				stringstream stpn;
				stpn<<"./n0/pin0"<<count_min<<"_"<<count<<".csv";
				ofstream fpn(stpn.str());
				stringstream stldg;
				stldg<<"./DgDq/ldg"<<count_min<<"_"<<count<<".csv";
				ofstream fldg(stldg.str());


				for(int i=0;i<h_num;i++)
				{
					fp<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
					fq<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
					fd<<dT[i]<<endl;
					fj<<HYPER[i].J<<endl;
					fl<<HYPER[i].lam<<endl;

					double pin0[DIMENSION]={0,0,0};
					double ldg[DIMENSION]={0,0,0};
					int Ni=HYPER[i].N;
					for(int j=0;j<Ni;j++)
					{		
						int k=HYPER[i].NEI[j];
						pin0[A_X]+=HYPER[k].pi_n[0][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[0][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[0][2]*HYPER1[i*h_num+k].n0ij[2];
						pin0[A_Y]+=HYPER[k].pi_n[1][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[1][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[1][2]*HYPER1[i*h_num+k].n0ij[2];
						pin0[A_Z]+=HYPER[k].pi_n[2][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[2][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[2][2]*HYPER1[i*h_num+k].n0ij[2];
						ldg[A_X]+=HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[0];
						ldg[A_Y]+=HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[1];
						ldg[A_Z]+=HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[2];
					}//jに関するfor文の終わり	
					pin0[A_X]+=HYPER[i].pi_n[0][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[0][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[0][2]*HYPER1[i*h_num+i].n0ij[2];
					pin0[A_Y]+=HYPER[i].pi_n[1][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[1][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[1][2]*HYPER1[i*h_num+i].n0ij[2];
					pin0[A_Z]+=HYPER[i].pi_n[2][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[2][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[2][2]*HYPER1[i*h_num+i].n0ij[2];
					ldg[A_X]+=HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[0];
					ldg[A_Y]+=HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[1];
					ldg[A_Z]+=HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[2];
					fpn<<pin0[A_X]<<","<<pin0[A_Y]<<","<<pin0[A_Z]<<endl;
					fldg<<ldg[A_X]<<","<<ldg[A_Y]<<","<<ldg[A_Z]<<endl;

				}
				fp.close();
				fq.close();
				fd.close();
				fj.close();
				fl.close();
				fpn.close();
				fldg.close();
			}
			if(count>100)	exit(0);
		}

		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=(HYPER[i].old_lam-HYPER[i].lam)*(HYPER[i].old_lam-HYPER[i].lam);
		E_min=sqrt(E_sum);
		if(E_min<ep_min*1000)	r*=4;
		cout<<"Emin"<<count_min<<"="<<E_min<<endl;

		for(int i=0;i<h_num;i++)
		{
			th_g[i]+=g[i]*g[i];
		}
	}

	delete[]	dE;
	delete[]	dT;
	delete[]	rE;
	delete[]	rT;
	delete[]	g;
	for(int i=0;i<h_num;i++)
	{
		delete[]	dg[i];
	}
	delete[]	dg;


	/////////////pn1_2計算
	for(int i=0;i<h_num;i++)
	{
		double p_half_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{		
			int k=HYPER[i].NEI[j];
			p_half_p[A_X]+=HYPER[k].pi_n[0][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[0][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[0][2]*HYPER1[i*h_num+k].n0ij[2]-HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[0];
			p_half_p[A_Y]+=HYPER[k].pi_n[1][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[1][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[1][2]*HYPER1[i*h_num+k].n0ij[2]-HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[1];
			p_half_p[A_Z]+=HYPER[k].pi_n[2][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[2][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[2][2]*HYPER1[i*h_num+k].n0ij[2]-HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[2];
		}//jに関するfor文の終わり	
		p_half_p[A_X]+=HYPER[i].pi_n[0][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[0][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[0][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[0];
		p_half_p[A_Y]+=HYPER[i].pi_n[1][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[1][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[1][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[1];
		p_half_p[A_Z]+=HYPER[i].pi_n[2][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[2][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[2][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[2];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_half_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_half_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_half_p[A_Z];
	
		/////////////q計算
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}

}

void q_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double nG[DIMENSION]={0,0,1};


	/////////////pn1_2計算
	for(int i=0;i<h_num;i++)
	{
		double p_half_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{		
			int k=HYPER[i].NEI[j];
			p_half_p[A_X]+=HYPER[k].pi_n[0][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[0][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[0][2]*HYPER1[i*h_num+k].n0ij[2]-HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[0];
			p_half_p[A_Y]+=HYPER[k].pi_n[1][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[1][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[1][2]*HYPER1[i*h_num+k].n0ij[2]-HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[1];
			p_half_p[A_Z]+=HYPER[k].pi_n[2][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].pi_n[2][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].pi_n[2][2]*HYPER1[i*h_num+k].n0ij[2]-HYPER[k].lam*HYPER1[k*h_num+i].DgDq_n[2];
		}//jに関するfor文の終わり	
		p_half_p[A_X]+=HYPER[i].pi_n[0][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[0][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[0][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[0];
		p_half_p[A_Y]+=HYPER[i].pi_n[1][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[1][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[1][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[1];
		p_half_p[A_Z]+=HYPER[i].pi_n[2][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[2][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[2][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[2];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_half_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_half_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_half_p[A_Z];
	
		/////////////q計算
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}


	/////////////F, J, t_inverseの計算
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		int Ni=HYPER[i].N;	

		for(int jn=0;jn<Ni;jn++)
		{
			int jnn=HYPER[i].NEI[jn];
			double w=HYPER1[i*h_num+jnn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+jnn].aiin[A_X],	HYPER1[i*h_num+jnn].aiin[A_Y],	HYPER1[i*h_num+jnn].aiin[A_Z]};
			
			fi[0][0]+=w*(PART[jnn].r[A_X]-PART[i].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[jnn].r[A_X]-PART[i].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[jnn].r[A_X]-PART[i].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[jnn].r[A_Y]-PART[i].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[jnn].r[A_Y]-PART[i].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[jnn].r[A_Y]-PART[i].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[jnn].r[A_Z]-PART[i].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[jnn].r[A_Z]-PART[i].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[jnn].r[A_Z]-PART[i].r[A_Z])*a[A_Z];
		}

		p_Fi[0][0]=fi[0][0]*HYPER[i].inverse_Ai[0][0]+fi[0][1]*HYPER[i].inverse_Ai[1][0]+fi[0][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[0][1]=fi[0][0]*HYPER[i].inverse_Ai[0][1]+fi[0][1]*HYPER[i].inverse_Ai[1][1]+fi[0][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[0][2]=fi[0][0]*HYPER[i].inverse_Ai[0][2]+fi[0][1]*HYPER[i].inverse_Ai[1][2]+fi[0][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[1][0]=fi[1][0]*HYPER[i].inverse_Ai[0][0]+fi[1][1]*HYPER[i].inverse_Ai[1][0]+fi[1][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[1][1]=fi[1][0]*HYPER[i].inverse_Ai[0][1]+fi[1][1]*HYPER[i].inverse_Ai[1][1]+fi[1][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[1][2]=fi[1][0]*HYPER[i].inverse_Ai[0][2]+fi[1][1]*HYPER[i].inverse_Ai[1][2]+fi[1][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[2][0]=fi[2][0]*HYPER[i].inverse_Ai[0][0]+fi[2][1]*HYPER[i].inverse_Ai[1][0]+fi[2][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[2][1]=fi[2][0]*HYPER[i].inverse_Ai[0][1]+fi[2][1]*HYPER[i].inverse_Ai[1][1]+fi[2][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[2][2]=fi[2][0]*HYPER[i].inverse_Ai[0][2]+fi[2][1]*HYPER[i].inverse_Ai[1][2]+fi[2][2]*HYPER[i].inverse_Ai[2][2];
		
		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	

		double J=calc_det3(p_Fi);
		HYPER[i].J=J;	

		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];
	}	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;


	/////////////DgDq, Stress, W, S, dSdc計算
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double **in_Ci=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	in_Ci[D]=new double[DIMENSION];
	double in_Ci2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double trace_dC=0;
	double trace_dC2=0;

	double Ic=0;
	double IIc=0;
	
	double S[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dSdc[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};


	double F_dSdC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double F_dSdC_F[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<h_num;j++)
	{	

		////////W計算
		double J=HYPER[j].J;	
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[j].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[j].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[0][0]*dC[0][0]+dC[0][1]*dC[1][0]+dC[0][2]*dC[2][0];
		dC2[0][1]=dC[0][0]*dC[0][1]+dC[0][1]*dC[1][1]+dC[0][2]*dC[2][1];
		dC2[0][2]=dC[0][0]*dC[0][2]+dC[0][1]*dC[1][2]+dC[0][2]*dC[2][2];
		dC2[1][0]=dC[1][0]*dC[0][0]+dC[1][1]*dC[1][0]+dC[1][2]*dC[2][0];
		dC2[1][1]=dC[1][0]*dC[0][1]+dC[1][1]*dC[1][1]+dC[1][2]*dC[2][1];
		dC2[1][2]=dC[1][0]*dC[0][2]+dC[1][1]*dC[1][2]+dC[1][2]*dC[2][2];
		dC2[2][0]=dC[2][0]*dC[0][0]+dC[2][1]*dC[1][0]+dC[2][2]*dC[2][0];
		dC2[2][1]=dC[2][0]*dC[0][1]+dC[2][1]*dC[1][1]+dC[2][2]*dC[2][1];
		dC2[2][2]=dC[2][0]*dC[0][2]+dC[2][1]*dC[1][2]+dC[2][2]*dC[2][2];

		trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		Ic=trace_dC;
		IIc=0.50*(trace_dC*trace_dC-trace_dC2);
	
		HYPER[j].W=c10*(Ic-3)+c01*(IIc-3);

	
		
		////////S, dSdc計算
		in_Ci[0][0]=dC[0][0];	in_Ci[0][1]=dC[0][1];	in_Ci[0][2]=dC[0][2];
		in_Ci[1][0]=dC[1][0];	in_Ci[1][1]=dC[1][1];	in_Ci[1][2]=dC[1][2];
		in_Ci[2][0]=dC[2][0];	in_Ci[2][1]=dC[2][1];	in_Ci[2][2]=dC[2][2];

		inverse(in_Ci,DIMENSION);

		in_Ci2[0][0]=in_Ci[0][0]*in_Ci[0][0]+in_Ci[0][1]*in_Ci[1][0]+in_Ci[0][2]*in_Ci[2][0];
		in_Ci2[0][1]=in_Ci[0][0]*in_Ci[0][1]+in_Ci[0][1]*in_Ci[1][1]+in_Ci[0][2]*in_Ci[2][1];
		in_Ci2[0][2]=in_Ci[0][0]*in_Ci[0][2]+in_Ci[0][1]*in_Ci[1][2]+in_Ci[0][2]*in_Ci[2][2];
		in_Ci2[1][0]=in_Ci[1][0]*in_Ci[0][0]+in_Ci[1][1]*in_Ci[1][0]+in_Ci[1][2]*in_Ci[2][0];
		in_Ci2[1][1]=in_Ci[1][0]*in_Ci[0][1]+in_Ci[1][1]*in_Ci[1][1]+in_Ci[1][2]*in_Ci[2][1];
		in_Ci2[1][2]=in_Ci[1][0]*in_Ci[0][2]+in_Ci[1][1]*in_Ci[1][2]+in_Ci[1][2]*in_Ci[2][2];
		in_Ci2[2][0]=in_Ci[2][0]*in_Ci[0][0]+in_Ci[2][1]*in_Ci[1][0]+in_Ci[2][2]*in_Ci[2][0];
		in_Ci2[2][1]=in_Ci[2][0]*in_Ci[0][1]+in_Ci[2][1]*in_Ci[1][1]+in_Ci[2][2]*in_Ci[2][1];
		in_Ci2[2][2]=in_Ci[2][0]*in_Ci[0][2]+in_Ci[2][1]*in_Ci[1][2]+in_Ci[2][2]*in_Ci[2][2];

		if(J<0)
		{
			S[0][0]=-2.*1/pow(-J,2/3.)*( c10*(1.-1/3.*Ic*in_Ci[0][0]) + c01*(Ic-dC[0][0]-2./3.*IIc*in_Ci[0][0]) );
			S[0][1]=-2.*1/pow(-J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[0][1]) + c01*(  -dC[0][1]-2./3.*IIc*in_Ci[0][1]) );
			S[0][2]=-2.*1/pow(-J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[0][2]) + c01*(  -dC[0][2]-2./3.*IIc*in_Ci[0][2]) );
			S[1][0]=-2.*1/pow(-J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[1][0]) + c01*(  -dC[1][0]-2./3.*IIc*in_Ci[1][0]) );
			S[1][1]=-2.*1/pow(-J,2/3.)*( c10*(1.-1/3.*Ic*in_Ci[1][1]) + c01*(Ic-dC[1][1]-2./3.*IIc*in_Ci[1][1]) );
			S[1][2]=-2.*1/pow(-J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[1][2]) + c01*(  -dC[1][2]-2./3.*IIc*in_Ci[1][2]) );
			S[2][0]=-2.*1/pow(-J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[2][0]) + c01*(  -dC[2][0]-2./3.*IIc*in_Ci[2][0]) );
			S[2][1]=-2.*1/pow(-J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[2][1]) + c01*(  -dC[2][1]-2./3.*IIc*in_Ci[2][1]) );	
			S[2][2]=-2.*1/pow(-J,2/3.)*( c10*(1.-1/3.*Ic*in_Ci[2][2]) + c01*(Ic-dC[2][2]-2./3.*IIc*in_Ci[2][2]) );

			dSdc[0][0]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[0][0]-in_Ci[0][0])+c01*(-2.*Ic*in_Ci[0][0]+5./3.*(1+IIc*in_Ci2[0][0])) );
			dSdc[0][1]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[0][1]-in_Ci[0][1])+c01*(-2.*Ic*in_Ci[0][1]+5./3.*( +IIc*in_Ci2[0][1])) );
			dSdc[0][2]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[0][2]-in_Ci[0][2])+c01*(-2.*Ic*in_Ci[0][2]+5./3.*( +IIc*in_Ci2[0][2])) );
			dSdc[1][0]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[1][0]-in_Ci[1][0])+c01*(-2.*Ic*in_Ci[1][0]+5./3.*( +IIc*in_Ci2[1][0])) );
			dSdc[1][1]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[1][1]-in_Ci[1][1])+c01*(-2.*Ic*in_Ci[1][1]+5./3.*(1+IIc*in_Ci2[1][1])) );
			dSdc[1][2]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[1][2]-in_Ci[1][2])+c01*(-2.*Ic*in_Ci[1][2]+5./3.*( +IIc*in_Ci2[1][2])) );
			dSdc[2][0]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[2][0]-in_Ci[2][0])+c01*(-2.*Ic*in_Ci[2][0]+5./3.*( +IIc*in_Ci2[2][0])) );
			dSdc[2][1]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[2][1]-in_Ci[2][1])+c01*(-2.*Ic*in_Ci[2][1]+5./3.*( +IIc*in_Ci2[2][1])) );
			dSdc[2][2]=-4./3.*1/pow(-J,4/3.)*( c10*(2./3.*Ic*in_Ci2[2][2]-in_Ci[2][2])+c01*(-2.*Ic*in_Ci[2][2]+5./3.*(1+IIc*in_Ci2[2][2])) );
		}
		else
		{
			S[0][0]=2.*1/pow(J,2/3.)*( c10*(1.-1/3.*Ic*in_Ci[0][0]) + c01*(Ic-dC[0][0]-2./3.*IIc*in_Ci[0][0]) );
			S[0][1]=2.*1/pow(J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[0][1]) + c01*(  -dC[0][1]-2./3.*IIc*in_Ci[0][1]) );
			S[0][2]=2.*1/pow(J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[0][2]) + c01*(  -dC[0][2]-2./3.*IIc*in_Ci[0][2]) );
			S[1][0]=2.*1/pow(J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[1][0]) + c01*(  -dC[1][0]-2./3.*IIc*in_Ci[1][0]) );
			S[1][1]=2.*1/pow(J,2/3.)*( c10*(1.-1/3.*Ic*in_Ci[1][1]) + c01*(Ic-dC[1][1]-2./3.*IIc*in_Ci[1][1]) );
			S[1][2]=2.*1/pow(J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[1][2]) + c01*(  -dC[1][2]-2./3.*IIc*in_Ci[1][2]) );
			S[2][0]=2.*1/pow(J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[2][0]) + c01*(  -dC[2][0]-2./3.*IIc*in_Ci[2][0]) );
			S[2][1]=2.*1/pow(J,2/3.)*( c10*(  -1/3.*Ic*in_Ci[2][1]) + c01*(  -dC[2][1]-2./3.*IIc*in_Ci[2][1]) );	
			S[2][2]=2.*1/pow(J,2/3.)*( c10*(1.-1/3.*Ic*in_Ci[2][2]) + c01*(Ic-dC[2][2]-2./3.*IIc*in_Ci[2][2]) );


			dSdc[0][0]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[0][0]-in_Ci[0][0])+c01*(-2.*Ic*in_Ci[0][0]+5./3.*(1+IIc*in_Ci2[0][0])) );
			dSdc[0][1]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[0][1]-in_Ci[0][1])+c01*(-2.*Ic*in_Ci[0][1]+5./3.*( +IIc*in_Ci2[0][1])) );
			dSdc[0][2]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[0][2]-in_Ci[0][2])+c01*(-2.*Ic*in_Ci[0][2]+5./3.*( +IIc*in_Ci2[0][2])) );
			dSdc[1][0]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[1][0]-in_Ci[1][0])+c01*(-2.*Ic*in_Ci[1][0]+5./3.*( +IIc*in_Ci2[1][0])) );
			dSdc[1][1]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[1][1]-in_Ci[1][1])+c01*(-2.*Ic*in_Ci[1][1]+5./3.*(1+IIc*in_Ci2[1][1])) );
			dSdc[1][2]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[1][2]-in_Ci[1][2])+c01*(-2.*Ic*in_Ci[1][2]+5./3.*( +IIc*in_Ci2[1][2])) );
			dSdc[2][0]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[2][0]-in_Ci[2][0])+c01*(-2.*Ic*in_Ci[2][0]+5./3.*( +IIc*in_Ci2[2][0])) );
			dSdc[2][1]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[2][1]-in_Ci[2][1])+c01*(-2.*Ic*in_Ci[2][1]+5./3.*( +IIc*in_Ci2[2][1])) );
			dSdc[2][2]=4./3.*1/pow(J,4/3.)*( c10*(2./3.*Ic*in_Ci2[2][2]-in_Ci[2][2])+c01*(-2.*Ic*in_Ci[2][2]+5./3.*(1+IIc*in_Ci2[2][2])) );
		}
		F_dSdC[0][0]=HYPER[j].Fi[0][0]*dSdc[0][0]+HYPER[j].Fi[0][1]*dSdc[1][0]+HYPER[j].Fi[0][2]*dSdc[2][0];
		F_dSdC[0][1]=HYPER[j].Fi[0][0]*dSdc[0][1]+HYPER[j].Fi[0][1]*dSdc[1][1]+HYPER[j].Fi[0][2]*dSdc[2][1];
		F_dSdC[0][2]=HYPER[j].Fi[0][0]*dSdc[0][2]+HYPER[j].Fi[0][1]*dSdc[1][2]+HYPER[j].Fi[0][2]*dSdc[2][2];
		F_dSdC[1][0]=HYPER[j].Fi[1][0]*dSdc[0][0]+HYPER[j].Fi[1][1]*dSdc[1][0]+HYPER[j].Fi[1][2]*dSdc[2][0];
		F_dSdC[1][1]=HYPER[j].Fi[1][0]*dSdc[0][1]+HYPER[j].Fi[1][1]*dSdc[1][1]+HYPER[j].Fi[1][2]*dSdc[2][1];
		F_dSdC[1][2]=HYPER[j].Fi[1][0]*dSdc[0][2]+HYPER[j].Fi[1][1]*dSdc[1][2]+HYPER[j].Fi[1][2]*dSdc[2][2];
		F_dSdC[2][0]=HYPER[j].Fi[2][0]*dSdc[0][0]+HYPER[j].Fi[2][1]*dSdc[1][0]+HYPER[j].Fi[2][2]*dSdc[2][0];
		F_dSdC[2][1]=HYPER[j].Fi[2][0]*dSdc[0][1]+HYPER[j].Fi[2][1]*dSdc[1][1]+HYPER[j].Fi[2][2]*dSdc[2][1];
		F_dSdC[2][2]=HYPER[j].Fi[2][0]*dSdc[0][2]+HYPER[j].Fi[2][1]*dSdc[1][2]+HYPER[j].Fi[2][2]*dSdc[2][2];

		F_dSdC_F[0][0]=F_dSdC[0][0]*HYPER[j].Fi[0][0]+F_dSdC[0][1]*HYPER[j].Fi[0][1]+F_dSdC[0][2]*HYPER[j].Fi[0][2];
		F_dSdC_F[0][1]=F_dSdC[0][0]*HYPER[j].Fi[1][0]+F_dSdC[0][1]*HYPER[j].Fi[1][1]+F_dSdC[0][2]*HYPER[j].Fi[1][2];
		F_dSdC_F[0][2]=F_dSdC[0][0]*HYPER[j].Fi[2][0]+F_dSdC[0][1]*HYPER[j].Fi[2][1]+F_dSdC[0][2]*HYPER[j].Fi[2][2];
		F_dSdC_F[1][0]=F_dSdC[1][0]*HYPER[j].Fi[0][0]+F_dSdC[1][1]*HYPER[j].Fi[0][1]+F_dSdC[1][2]*HYPER[j].Fi[0][2];
		F_dSdC_F[1][1]=F_dSdC[1][0]*HYPER[j].Fi[1][0]+F_dSdC[1][1]*HYPER[j].Fi[1][1]+F_dSdC[1][2]*HYPER[j].Fi[1][2];
		F_dSdC_F[1][2]=F_dSdC[1][0]*HYPER[j].Fi[2][0]+F_dSdC[1][1]*HYPER[j].Fi[2][1]+F_dSdC[1][2]*HYPER[j].Fi[2][2];
		F_dSdC_F[2][0]=F_dSdC[2][0]*HYPER[j].Fi[0][0]+F_dSdC[2][1]*HYPER[j].Fi[0][1]+F_dSdC[2][2]*HYPER[j].Fi[0][2];
		F_dSdC_F[2][1]=F_dSdC[2][0]*HYPER[j].Fi[1][0]+F_dSdC[2][1]*HYPER[j].Fi[1][1]+F_dSdC[2][2]*HYPER[j].Fi[1][2];
		F_dSdC_F[2][2]=F_dSdC[2][0]*HYPER[j].Fi[2][0]+F_dSdC[2][1]*HYPER[j].Fi[1][1]+F_dSdC[2][2]*HYPER[j].Fi[2][2];

		HYPER[j].pi[0][0]=HYPER[j].Fi[0][0]*S[0][0]+HYPER[j].Fi[0][1]*S[1][0]+HYPER[j].Fi[0][2]*S[2][0];
		HYPER[j].pi[0][1]=HYPER[j].Fi[0][0]*S[0][1]+HYPER[j].Fi[0][1]*S[1][1]+HYPER[j].Fi[0][2]*S[2][1];
		HYPER[j].pi[0][2]=HYPER[j].Fi[0][0]*S[0][2]+HYPER[j].Fi[0][1]*S[1][2]+HYPER[j].Fi[0][2]*S[2][2];
		HYPER[j].pi[1][0]=HYPER[j].Fi[1][0]*S[0][0]+HYPER[j].Fi[1][1]*S[1][0]+HYPER[j].Fi[1][2]*S[2][0];
		HYPER[j].pi[1][1]=HYPER[j].Fi[1][0]*S[0][1]+HYPER[j].Fi[1][1]*S[1][1]+HYPER[j].Fi[1][2]*S[2][1];
		HYPER[j].pi[1][2]=HYPER[j].Fi[1][0]*S[0][2]+HYPER[j].Fi[1][1]*S[1][2]+HYPER[j].Fi[1][2]*S[2][2];
		HYPER[j].pi[2][0]=HYPER[j].Fi[2][0]*S[0][0]+HYPER[j].Fi[2][1]*S[1][0]+HYPER[j].Fi[2][2]*S[2][0];
		HYPER[j].pi[2][1]=HYPER[j].Fi[2][0]*S[0][1]+HYPER[j].Fi[2][1]*S[1][1]+HYPER[j].Fi[2][2]*S[2][1];
		HYPER[j].pi[2][2]=HYPER[j].Fi[2][0]*S[0][2]+HYPER[j].Fi[2][1]*S[1][2]+HYPER[j].Fi[2][2]*S[2][2];

		HYPER[j].dpidF[0][0]=S[0][0]+F_dSdC_F[0][0];
		HYPER[j].dpidF[0][1]=S[0][1]+F_dSdC_F[0][1];
		HYPER[j].dpidF[0][2]=S[0][2]+F_dSdC_F[0][2];
		HYPER[j].dpidF[1][0]=S[1][0]+F_dSdC_F[1][0];
		HYPER[j].dpidF[1][1]=S[1][1]+F_dSdC_F[1][1];
		HYPER[j].dpidF[1][2]=S[1][2]+F_dSdC_F[1][2];
		HYPER[j].dpidF[2][0]=S[2][0]+F_dSdC_F[2][0];
		HYPER[j].dpidF[2][1]=S[2][1]+F_dSdC_F[2][1];
		HYPER[j].dpidF[2][2]=S[2][2]+F_dSdC_F[2][2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	in_Ci[D];
	delete[]	in_Ci;
}

void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,double nG[DIMENSION])
{

	int h_num=HYPER.size();
	int Nx=2*h_num;

	/////////////////dE計算///////////////////
	for(int k=0;k<h_num;k++)
	{
		dE[k]=0;

		for(int i=0;i<h_num;i++)
		{
			dg[i][k]=0;
			for(int m=0;m<h_num;m++)
			{
				dE[k]+=0.5*Dt*Dt/mi*((HYPER[i].pi[A_X][0]*HYPER1[m*h_num+i].n0ij[0]+HYPER[i].pi[A_X][1]*HYPER1[m*h_num+i].n0ij[1]+HYPER[i].pi[A_X][2]*HYPER1[m*h_num+i].n0ij[2])*HYPER1[k*h_num+m].DgDq_n[A_X]
									 +(HYPER[i].pi[A_Y][0]*HYPER1[m*h_num+i].n0ij[0]+HYPER[i].pi[A_Y][1]*HYPER1[m*h_num+i].n0ij[1]+HYPER[i].pi[A_Y][2]*HYPER1[m*h_num+i].n0ij[2])*HYPER1[k*h_num+m].DgDq_n[A_Y]
									 +(HYPER[i].pi[A_Z][0]*HYPER1[m*h_num+i].n0ij[0]+HYPER[i].pi[A_Z][1]*HYPER1[m*h_num+i].n0ij[1]+HYPER[i].pi[A_Z][2]*HYPER1[m*h_num+i].n0ij[2])*HYPER1[k*h_num+m].DgDq_n[A_Z]);			

				dg[i][k]+=-0.5*Dt*Dt/mi*HYPER[i].J*((HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[m*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[m*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[m*h_num+i].n0ij[2])*HYPER1[k*h_num+m].DgDq_n[A_X]
											   	   +(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[m*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[m*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[m*h_num+i].n0ij[2])*HYPER1[k*h_num+m].DgDq_n[A_Y]
												   +(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[m*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[m*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[m*h_num+i].n0ij[2])*HYPER1[k*h_num+m].DgDq_n[A_Z]);
			}
			dE[k]+=-0.5*Dt/mi*(HYPER1[k*h_num+i].DgDq_n[A_X]*HYPER[i].half_p[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*HYPER[i].half_p[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*HYPER[i].half_p[A_Z]);
		}
	}


	/////////////////rE計算///////////////////
	double llE=0;
	double mlE=0;
	double lmE=0;
	double mmE=0;
	for(int l=0;l<h_num;l++)
	{
		for(int k=0;k<h_num;k++)
		{
			rE[l*h_num+k]=0;
			for(int i=0;i<h_num;i++)
			{
				int Ni=HYPER[i].N;
				for(int m=0;m<Ni;m++)
				{				
					int mn=HYPER[i].NEI[m];
					for(int s=0;s<Ni;s++)
					{
						int sn=HYPER[i].NEI[s];
						rE[l*h_num+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
							((HYPER[i].dpidF[A_X][A_X]*HYPER1[sn*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[sn*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[sn*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+sn].DgDq_n[A_X]
						    +(HYPER[i].dpidF[A_Y][A_X]*HYPER1[sn*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[sn*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[sn*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+sn].DgDq_n[A_Y]
							+(HYPER[i].dpidF[A_Z][A_X]*HYPER1[sn*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[sn*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[sn*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+sn].DgDq_n[A_Z])
							*(HYPER1[mn*h_num+i].n0ij[A_X]*HYPER1[k*h_num+mn].DgDq_n[A_X]+HYPER1[mn*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+mn].DgDq_n[A_Y]+HYPER1[mn*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+mn].DgDq_n[A_Z]);
					}
					rE[l*h_num+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
						((HYPER[i].dpidF[A_X][A_X]*HYPER1[i*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[i*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[i*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+i].DgDq_n[A_X]
						+(HYPER[i].dpidF[A_Y][A_X]*HYPER1[i*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[i*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[i*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+i].DgDq_n[A_Y]
						+(HYPER[i].dpidF[A_Z][A_X]*HYPER1[i*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[i*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[i*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+i].DgDq_n[A_Z])
						*(HYPER1[mn*h_num+i].n0ij[A_X]*HYPER1[k*h_num+mn].DgDq_n[A_X]+HYPER1[mn*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+mn].DgDq_n[A_Y]+HYPER1[mn*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+mn].DgDq_n[A_Z]);
					rE[l*h_num+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
						((HYPER[i].dpidF[A_X][A_X]*HYPER1[mn*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[mn*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[mn*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+mn].DgDq_n[A_X]
						+(HYPER[i].dpidF[A_Y][A_X]*HYPER1[mn*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[mn*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[mn*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+mn].DgDq_n[A_Y]
						+(HYPER[i].dpidF[A_Z][A_X]*HYPER1[mn*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[mn*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[mn*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+mn].DgDq_n[A_Z])
						*(HYPER1[i*h_num+i].n0ij[A_X]*HYPER1[k*h_num+i].DgDq_n[A_X]+HYPER1[i*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+i].DgDq_n[A_Y]+HYPER1[i*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+i].DgDq_n[A_Z]);
				}
				rE[l*h_num+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
					((HYPER[i].dpidF[A_X][A_X]*HYPER1[i*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[i*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[i*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+i].DgDq_n[A_X]
					+(HYPER[i].dpidF[A_Y][A_X]*HYPER1[i*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[i*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[i*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+i].DgDq_n[A_Y]
					+(HYPER[i].dpidF[A_Z][A_X]*HYPER1[i*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[i*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[i*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+i].DgDq_n[A_Z])
					*(HYPER1[i*h_num+i].n0ij[A_X]*HYPER1[k*h_num+i].DgDq_n[A_X]+HYPER1[i*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+i].DgDq_n[A_Y]+HYPER1[i*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+i].DgDq_n[A_Z]);
				rE[l*h_num+k]+=0.25*Dt*Dt/mi*(HYPER1[k*h_num+i].DgDq_n[A_X]*HYPER1[l*h_num+i].DgDq_n[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*HYPER1[l*h_num+i].DgDq_n[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*HYPER1[l*h_num+i].DgDq_n[A_Z]);
			}
		}
	}
}



void p_QP(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0)
{
	int count=1;
	int count_min=1;
	int c_max=10000;

	double E=1;
	double E_min=1;
	double ep=1e-5;
	double ep_min=1e-5;

	double r=0;
	double r0=1000;

			r=r0*0.01;
			//
			//theta_dh=0;
			//mu=90;
			//E_min=1;
			//count_min=0;
			//while(E_min>ep_min)
			//{
			//	count_min++;

			//	old_mu=mu;

			//	E=1;
			//	count=0;
			//	while(E>ep)
			//	{
			//		count++;

			//		pn1[A_X]=pn1_2[A_X]-0.5*Dt*mi*G*nG[A_X]+0.5*Dt*n[A_X]*mu;	
			//		pn1[A_Y]=pn1_2[A_Y]-0.5*Dt*mi*G*nG[A_Y]+0.5*Dt*n[A_Y]*mu;	
			//		pn1[A_Z]=pn1_2[A_Z]-0.5*Dt*mi*G*nG[A_Z]+0.5*Dt*n[A_Z]*mu;	

			//		En=0.5/mi*(pn1[A_X]*pn1[A_X]+pn1[A_Y]*pn1[A_Y]+pn1[A_Z]*pn1[A_Z])+mi*G*(qn1[A_X]*nG[A_X]+qn1[A_Y]*nG[A_Y]+qn1[A_Z]*nG[A_Z]);
	
			//		Tr=(En-E0)*(En-E0);
			//		
			//		dE=0.5/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);
			//		dTr=2*dE*(En-E0);

			//		rE=0.25*Dt*Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
			//		rTr=2*rE*(En-E0)+2*dE*dE;

			//		dhdt=-1/mi*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

			//		d_dhdt=-0.5*Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
			//		if(dhdt+theta_dh>0)
			//		{
			//			Tr+=0.5*r*(dhdt+theta_dh)*(dhdt+theta_dh);
			//			dTr+=r*d_dhdt*(dhdt+theta_dh);
			//			rTr+=r*d_dhdt*d_dhdt;
			//		}
			//		else if(dhdt+theta_dh==0)
			//		{
			//			dTr+=r*d_dhdt*(dhdt+theta_dh);
			//			rTr+=r*d_dhdt*d_dhdt;
			//		}

			//		mu+=-1*dTr/rTr;
			//		E=sqrt(dTr/rTr*dTr/rTr);



			//}

			//	E_min=sqrt((old_mu-mu)*(old_mu-mu));	
			//	if(E_min<ep_min*1000)	r*=4;
			//	if(dhdt+theta_dh>0)	theta_dh+=dhdt;
	

			//	pn1[A_X]=pn1_2[A_X]-0.5*Dt*mi*G*nG[A_X]+0.5*Dt*n[A_X]*mu;	
			//	pn1[A_Y]=pn1_2[A_Y]-0.5*Dt*mi*G*nG[A_Y]+0.5*Dt*n[A_Y]*mu;	
			//	pn1[A_Z]=pn1_2[A_Z]-0.5*Dt*mi*G*nG[A_Z]+0.5*Dt*n[A_Z]*mu;	

			//	output_QP(t,count_min,E_min,1);
			//	if(count_min>c_max)	break;
			//}


}

void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT,double *rL,double *dL,int h_num,int count,int count_min,int t,double E)
{

	int Nx=4*h_num;

	stringstream ss0;
	ss0<<"./rT/rT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss1;
	ss1<<"./dT/dT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./h_lam/h_lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./lam/lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss4;
	ss4<<"./h_mu/h_mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss5;
	ss5<<"./mu/mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./rL/rL_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss7;
	ss7<<"./dL/dL_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./p/p_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./q/q_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./w/w_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./E/E_"<<t<<"_"<<count_min<<".csv";

	ofstream f_rt(ss0.str(), ios::trunc);
	ofstream f_dt(ss1.str(), ios::trunc);
	ofstream f_h_lam(ss2.str(), ios::trunc);
	ofstream f_lam(ss3.str(), ios::trunc);
	ofstream f_h_mu(ss4.str(), ios::trunc);
	ofstream f_mu(ss5.str(), ios::trunc);
	ofstream f_rL(ss6.str(), ios::trunc);
	ofstream f_dL(ss7.str(), ios::trunc);
	ofstream f_p(ss8.str(), ios::trunc);
	ofstream f_q(ss9.str(), ios::trunc);
	ofstream f_w(ss10.str(), ios::trunc);
	if(count==1)
	{
		ofstream f_E(ss11.str(), ios::trunc);
		f_E<<E<<endl;
		f_E.close();
	}	
	else
	{
		ofstream f_E(ss11.str(), ios::app);
		f_E<<E<<endl;
		f_E.close();
	}

	for(int i=0;i<Nx;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_rt<<rT[i*Nx+j]<<",";
			f_rL<<rL[i*Nx+j]<<",";
		}
		f_rt<<endl;
		f_rL<<endl;

		f_dt<<dT[i]<<endl;
		f_dL<<dL[i]<<endl;
	}
	for(int i=0;i<h_num;i++)
	{
		f_h_lam<<HYPER[i].h_lam<<endl;
		f_lam<<HYPER[i].lam<<endl;
		f_h_mu<<HYPER[i].h_mu<<endl;
		f_mu<<HYPER[i].mu<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		f_w<<HYPER[i].W<<endl;
	}

	f_rt.close();
	f_dt.close();
	f_h_lam.close();
	f_lam.close();
	f_h_mu.close();
	f_mu.close();
	f_rL.close();
	f_dL.close();
	f_p.close();
	f_q.close();
	f_w.close();
}


////////Hyperelastic2からの引継ぎ
void calc_W(mpsconfig &CON,vector<hyperelastic> &HYPER,int h_num)
{
	//	cout<<"弾性ポテンシャル計算";

	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();

	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	for(int i=0;i<h_num;i++)
	{
		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		HYPER[i].W=c10*(Ic-3)+c01*(IIc-3);
	}
//	cout<<"----------OK"<<endl;
}


void calc_pi(mpsconfig &CON,vector<hyperelastic> &HYPER)
{
	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double nG[DIMENSION]={0,0,1};

	/////////////DgDq, Stress, W, S, dSdc計算
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double trace_dC=0;
	double trace_dC2=0;
	double Ic=0;
	double IIc=0;

	double **in_Ci=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	in_Ci[D]=new double[DIMENSION];

	double first_term[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double S[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};


	for(int j=0;j<h_num;j++)
	{	

		double J=HYPER[j].J;	
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3.)*HYPER[j].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3.)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3.)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3.)*HYPER[j].Fi[2][2];
		}
		
		//cout<<"dF"<<j<<"="<<d_Fi[0][0]<<","<<d_Fi[0][1]<<","<<d_Fi[0][2]<<endl;
		//cout<<d_Fi[1][0]<<","<<d_Fi[1][1]<<","<<d_Fi[1][2]<<endl;
		//cout<<d_Fi[2][0]<<","<<d_Fi[2][1]<<","<<d_Fi[2][2]<<endl;

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[0][0]*dC[0][0]+dC[0][1]*dC[1][0]+dC[0][2]*dC[2][0];
		dC2[0][1]=dC[0][0]*dC[0][1]+dC[0][1]*dC[1][1]+dC[0][2]*dC[2][1];
		dC2[0][2]=dC[0][0]*dC[0][2]+dC[0][1]*dC[1][2]+dC[0][2]*dC[2][2];
		dC2[1][0]=dC[1][0]*dC[0][0]+dC[1][1]*dC[1][0]+dC[1][2]*dC[2][0];
		dC2[1][1]=dC[1][0]*dC[0][1]+dC[1][1]*dC[1][1]+dC[1][2]*dC[2][1];
		dC2[1][2]=dC[1][0]*dC[0][2]+dC[1][1]*dC[1][2]+dC[1][2]*dC[2][2];
		dC2[2][0]=dC[2][0]*dC[0][0]+dC[2][1]*dC[1][0]+dC[2][2]*dC[2][0];
		dC2[2][1]=dC[2][0]*dC[0][1]+dC[2][1]*dC[1][1]+dC[2][2]*dC[2][1];
		dC2[2][2]=dC[2][0]*dC[0][2]+dC[2][1]*dC[1][2]+dC[2][2]*dC[2][2];

		//cout<<"dC"<<j<<"="<<dC[0][0]<<","<<dC[0][1]<<","<<dC[0][2]<<endl;
		//cout<<dC[1][0]<<","<<dC[1][1]<<","<<dC[1][2]<<endl;
		//cout<<dC[2][0]<<","<<dC[2][1]<<","<<dC[2][2]<<endl;
		//
		//cout<<"dC2"<<j<<"="<<dC2[0][0]<<","<<dC2[0][1]<<","<<dC2[0][2]<<endl;
		//cout<<dC2[1][0]<<","<<dC2[1][1]<<","<<dC2[1][2]<<endl;
		//cout<<dC2[2][0]<<","<<dC2[2][1]<<","<<dC2[2][2]<<endl;


		trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		//cout<<"trace_dC="<<trace_dC<<", trace_dC2="<<trace_dC2;

		Ic=trace_dC;
		IIc=0.50*(trace_dC*trace_dC-trace_dC2);

		//cout<<" ,Ic="<<Ic<<", IIc="<<IIc<<endl;
	
		
		////////S計算
		in_Ci[0][0]=dC[0][0];	in_Ci[0][1]=dC[0][1];	in_Ci[0][2]=dC[0][2];
		in_Ci[1][0]=dC[1][0];	in_Ci[1][1]=dC[1][1];	in_Ci[1][2]=dC[1][2];
		in_Ci[2][0]=dC[2][0];	in_Ci[2][1]=dC[2][1];	in_Ci[2][2]=dC[2][2];

		inverse(in_Ci,DIMENSION);


		//cout<<"in_Ci"<<j<<"="<<in_Ci[0][0]<<","<<in_Ci[0][1]<<","<<in_Ci[0][2]<<endl;
		//cout<<in_Ci[1][0]<<","<<in_Ci[1][1]<<","<<in_Ci[1][2]<<endl;
		//cout<<in_Ci[2][0]<<","<<in_Ci[2][1]<<","<<in_Ci[2][2]<<endl;

		first_term[0][0]=1-1/3.*Ic*in_Ci[0][0];	first_term[0][1]=-1/3.*Ic*in_Ci[0][1];	first_term[0][2]=-1/3.*Ic*in_Ci[0][2];
		first_term[1][0]=-1/3.*Ic*in_Ci[1][0];	first_term[1][1]=1-1/3.*Ic*in_Ci[1][1];	first_term[1][2]=-1/3.*Ic*in_Ci[1][0];
		first_term[2][0]=-1/3.*Ic*in_Ci[2][0];	first_term[2][1]=-1/3.*Ic*in_Ci[2][1];	first_term[2][2]=1-1/3.*Ic*in_Ci[2][2];
		
		//cout<<"1/3.*Ic*in_Ci[0][0]="<<1/3.*Ic*in_Ci[0][0]<<endl;
		//cout<<"in_Ci[0][0]="<<in_Ci[0][0]<<endl;
		//cout<<"Ic="<<Ic<<endl;
		//cout<<"Ic*in_Ci[0][0]="<<Ic*in_Ci[0][0]<<endl;
		//cout<<"first_term"<<j<<"="<<first_term[0][0]<<","<<first_term[0][1]<<","<<first_term[0][2]<<endl;
		//cout<<first_term[1][0]<<","<<first_term[1][1]<<","<<first_term[1][2]<<endl;
		//cout<<first_term[2][0]<<","<<first_term[2][1]<<","<<first_term[2][2]<<endl;

		if(J<0)
		{
			S[0][0]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3.*trace_dC2*in_Ci[0][0]) );
			S[0][1]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3.*trace_dC2*in_Ci[0][1]) );
			S[0][2]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3.*trace_dC2*in_Ci[0][2]) );
			S[1][0]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3.*trace_dC2*in_Ci[1][0]) );
			S[1][1]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3.*trace_dC2*in_Ci[1][1]) );
			S[1][2]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3.*trace_dC2*in_Ci[1][2]) );
			S[2][0]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3.*trace_dC2*in_Ci[2][0]) );
			S[2][1]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3.*trace_dC2*in_Ci[2][1]) );		
			S[2][2]=-2*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3.*trace_dC2*in_Ci[2][2]) );
		}
		else
		{
			S[0][0]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3.*trace_dC2*in_Ci[0][0]) );
			S[0][1]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3.*trace_dC2*in_Ci[0][1]) );
			S[0][2]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3.*trace_dC2*in_Ci[0][2]) );
			S[1][0]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3.*trace_dC2*in_Ci[1][0]) );
			S[1][1]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3.*trace_dC2*in_Ci[1][1]) );
			S[1][2]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3.*trace_dC2*in_Ci[1][2]) );
			S[2][0]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3.*trace_dC2*in_Ci[2][0]) );
			S[2][1]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3.*trace_dC2*in_Ci[2][1]) );		
			S[2][2]=2*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3.*trace_dC2*in_Ci[2][2]) );
		}
		//cout<<"S"<<j<<"="<<S[0][0]<<","<<S[0][1]<<","<<S[0][2]<<endl;
		//cout<<S[1][0]<<","<<S[1][1]<<","<<S[1][2]<<endl;
		//cout<<S[2][0]<<","<<S[2][1]<<","<<S[2][2]<<endl;



		HYPER[j].pi[0][0]=HYPER[j].Fi[0][0]*S[0][0]+HYPER[j].Fi[0][1]*S[1][0]+HYPER[j].Fi[0][2]*S[2][0];
		HYPER[j].pi[0][1]=HYPER[j].Fi[0][0]*S[0][1]+HYPER[j].Fi[0][1]*S[1][1]+HYPER[j].Fi[0][2]*S[2][1];
		HYPER[j].pi[0][2]=HYPER[j].Fi[0][0]*S[0][2]+HYPER[j].Fi[0][1]*S[1][2]+HYPER[j].Fi[0][2]*S[2][2];
		HYPER[j].pi[1][0]=HYPER[j].Fi[1][0]*S[0][0]+HYPER[j].Fi[1][1]*S[1][0]+HYPER[j].Fi[1][2]*S[2][0];
		HYPER[j].pi[1][1]=HYPER[j].Fi[1][0]*S[0][1]+HYPER[j].Fi[1][1]*S[1][1]+HYPER[j].Fi[1][2]*S[2][1];
		HYPER[j].pi[1][2]=HYPER[j].Fi[1][0]*S[0][2]+HYPER[j].Fi[1][1]*S[1][2]+HYPER[j].Fi[1][2]*S[2][2];
		HYPER[j].pi[2][0]=HYPER[j].Fi[2][0]*S[0][0]+HYPER[j].Fi[2][1]*S[1][0]+HYPER[j].Fi[2][2]*S[2][0];
		HYPER[j].pi[2][1]=HYPER[j].Fi[2][0]*S[0][1]+HYPER[j].Fi[2][1]*S[1][1]+HYPER[j].Fi[2][2]*S[2][1];
		HYPER[j].pi[2][2]=HYPER[j].Fi[2][0]*S[0][2]+HYPER[j].Fi[2][1]*S[1][2]+HYPER[j].Fi[2][2]*S[2][2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	in_Ci[D];
	delete[]	in_Ci;
}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
//	cout<<"初期値計算";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//考慮が必要かもしれない
	double mi=V*CON.get_hyper_density();
	int h_num=HYPER.size();
	int model=CON.get_model_number();
//	cout<<"V"<<V<<endl;

	ofstream fc("constant.dat");
	fc<<"le"<<","<<le<<endl;
	fc<<"r"<<","<<r<<endl;
	fc<<"Dt"<<","<<Dt<<endl;

	//垂直降下
	//for(int i=0;i<h_num;i++)	HYPER[i].p[A_Z]=-1.0*mi;
	
	//曲げねじり
	/*if(model==21)
	{		
		int t=30,b=2;
		double min=PART[0].q0[A_Z];
		double max=PART[0].q0[A_Z];

		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double Y=PART[i].q0[A_Y];
			double X=PART[i].q0[A_X];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*(-t*part_p*part_p*part_p*Y+b*(3*part_p*part_p-1));
			HYPER[i].p[A_Y]=mi*t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}//*/

	//曲げ
	if(model==21)
	{
		int b=10;
		double max=0,min=0;

		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*b*(3*part_p*part_p-1);
		}
	}//*/

	//ねじり
	/*if(model==21)
	{
		int t=30;
		double max=0,min=0;

		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double Y=PART[i].q0[A_Y];
			double X=PART[i].q0[A_X];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=-mi*t*part_p*part_p*part_p*Y;
			HYPER[i].p[A_Y]=mi*t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}//*/


	//回転
	if(model==22)
	{
		for(int i=0;i<h_num;i++)
		{
			PART[i].p[A_X]=mi*0.4*(PART[i].q0[A_Z]-PART[i].q0[A_Y]);
			PART[i].p[A_Y]=mi*0.4*(PART[i].q0[A_X]-PART[i].q0[A_Z]);
			PART[i].p[A_Z]=mi*0.4*(PART[i].q0[A_Y]-PART[i].q0[A_X]);
		}
	}

	////角運動量計算
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
		HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
		HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
	}


	////近傍粒子の記憶とaiin,wiin,Aiの計算
	for(int i=0;i<h_num;i++)
	{
		int N=0;
		double dis=0;
		double wiin=0;
		double aiin[DIMENSION];
		for(int j=0;j<h_num;j++)
		{
			wiin=0;
			aiin[A_X]=PART[j].q0[A_X]-PART[i].q0[A_X];	aiin[A_Y]=PART[j].q0[A_Y]-PART[i].q0[A_Y];	aiin[A_Z]=PART[j].q0[A_Z]-PART[i].q0[A_Z];
			
			HYPER1[i*h_num+j].aiin[A_X]=aiin[A_X];
			HYPER1[i*h_num+j].aiin[A_Y]=aiin[A_Y];
			HYPER1[i*h_num+j].aiin[A_Z]=aiin[A_Z];

			dis=sqrt(aiin[A_X]*aiin[A_X]+aiin[A_Y]*aiin[A_Y]+aiin[A_Z]*aiin[A_Z]);
			if(dis<r && j!=i)
			{	
				wiin=kernel4(r,dis);
				HYPER[i].NEI[N]=j;
				HYPER[i].pnd+=wiin;
				N++;
			}
			else	wiin=0;

			HYPER1[i*h_num+j].wiin=wiin;
		}
		HYPER[i].N=N;
		HYPER[i].pnd0=N;
	}
	

	////Ai, Fi関連の計算
	double **p_Ai=new double *[DIMENSION];
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Ai[D]=new double [DIMENSION];
		p_Fi[D]=new double[DIMENSION];
	}
	for(int i=0;i<h_num;i++)
	{
		//Aiの計算
		double Ai[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		Ai[0][0]=0;	Ai[0][1]=0;	Ai[0][2]=0;
		Ai[1][0]=0;	Ai[1][1]=0;	Ai[1][2]=0;
		Ai[2][0]=0;	Ai[2][1]=0;	Ai[2][2]=0;

		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)	
		{
			int k=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+k].wiin;		
			double a[DIMENSION]={HYPER1[i*h_num+k].aiin[A_X],	a[A_Y]=HYPER1[i*h_num+k].aiin[A_Y], a[A_Z]=HYPER1[i*h_num+k].aiin[A_Z]};
			
			Ai[0][0]+=w*a[0]*a[0];		Ai[0][1]+=w*a[0]*a[1];		Ai[0][2]+=w*a[0]*a[2];
			Ai[1][0]+=w*a[1]*a[0];		Ai[1][1]+=w*a[1]*a[1];		Ai[1][2]+=w*a[1]*a[2];
			Ai[2][0]+=w*a[2]*a[0];		Ai[2][1]+=w*a[2]*a[1];		Ai[2][2]+=w*a[2]*a[2];
		}
		HYPER[i].Ai[0][0]=Ai[0][0];	HYPER[i].Ai[0][1]=Ai[0][1];	HYPER[i].Ai[0][2]=Ai[0][2];
		HYPER[i].Ai[1][0]=Ai[1][0];	HYPER[i].Ai[1][1]=Ai[1][1];	HYPER[i].Ai[1][2]=Ai[1][2];
		HYPER[i].Ai[2][0]=Ai[2][0];	HYPER[i].Ai[2][1]=Ai[2][1];	HYPER[i].Ai[2][2]=Ai[2][2];
		
		//inverse_Ai,t_inverse_Aiの計算
		p_Ai[0][0]=Ai[0][0];	p_Ai[0][1]=Ai[0][1];	p_Ai[0][2]=Ai[0][2];
		p_Ai[1][0]=Ai[1][0];	p_Ai[1][1]=Ai[1][1];	p_Ai[1][2]=Ai[1][2];
		p_Ai[2][0]=Ai[2][0];	p_Ai[2][1]=Ai[2][1];	p_Ai[2][2]=Ai[2][2];
		inverse(p_Ai,DIMENSION);
		HYPER[i].inverse_Ai[0][0]=p_Ai[0][0];		HYPER[i].inverse_Ai[0][1]=p_Ai[0][1];		HYPER[i].inverse_Ai[0][2]=p_Ai[0][2];
		HYPER[i].inverse_Ai[1][0]=p_Ai[1][0];		HYPER[i].inverse_Ai[1][1]=p_Ai[1][1];		HYPER[i].inverse_Ai[1][2]=p_Ai[1][2];
		HYPER[i].inverse_Ai[2][0]=p_Ai[2][0];		HYPER[i].inverse_Ai[2][1]=p_Ai[2][1];		HYPER[i].inverse_Ai[2][2]=p_Ai[2][2];

		HYPER[i].t_inverse_Ai[0][0]=p_Ai[0][0];		HYPER[i].t_inverse_Ai[0][1]=p_Ai[1][0];		HYPER[i].t_inverse_Ai[0][2]=p_Ai[2][0];
		HYPER[i].t_inverse_Ai[1][0]=p_Ai[0][1];		HYPER[i].t_inverse_Ai[1][1]=p_Ai[1][1];		HYPER[i].t_inverse_Ai[1][2]=p_Ai[2][1];
		HYPER[i].t_inverse_Ai[2][0]=p_Ai[0][2];		HYPER[i].t_inverse_Ai[2][1]=p_Ai[1][2];		HYPER[i].t_inverse_Ai[2][2]=p_Ai[2][2];
		
		//Fiの計算
		p_Fi[0][0]=Ai[0][0]*p_Ai[0][0]+Ai[0][1]*p_Ai[1][0]+Ai[0][2]*p_Ai[2][0];	p_Fi[0][1]=Ai[0][0]*p_Ai[0][1]+Ai[0][1]*p_Ai[1][1]+Ai[0][2]*p_Ai[2][1];	p_Fi[0][2]=Ai[0][0]*p_Ai[0][2]+Ai[0][1]*p_Ai[1][2]+Ai[0][2]*p_Ai[2][2];
		p_Fi[1][0]=Ai[1][0]*p_Ai[0][0]+Ai[1][1]*p_Ai[1][0]+Ai[1][2]*p_Ai[2][0];	p_Fi[1][1]=Ai[1][0]*p_Ai[0][1]+Ai[1][1]*p_Ai[1][1]+Ai[1][2]*p_Ai[2][1];	p_Fi[1][2]=Ai[1][0]*p_Ai[0][2]+Ai[1][1]*p_Ai[1][2]+Ai[1][2]*p_Ai[2][2];
		p_Fi[2][0]=Ai[2][0]*p_Ai[0][0]+Ai[2][1]*p_Ai[1][0]+Ai[2][2]*p_Ai[2][0];	p_Fi[2][1]=Ai[2][0]*p_Ai[0][1]+Ai[2][1]*p_Ai[1][1]+Ai[2][2]*p_Ai[2][1];	p_Fi[2][2]=Ai[2][0]*p_Ai[0][2]+Ai[2][1]*p_Ai[1][2]+Ai[2][2]*p_Ai[2][2];		

		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	
		//Jの計算
		double J=calc_det3(p_Fi);
		HYPER[i].J=J;
//		cout<<"J["<<i<<"]="<<HYPER[i].J<<endl;

		//t_inverse_Fiの計算
		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];		HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];		HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];		HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];		HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];		HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];		HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];
	}
	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Ai[D];
		delete[]	p_Fi[D];
	}
	delete[]	p_Ai;
	delete[]	p_Fi;


	////n0ijの計算
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;

		double p_n0ij[DIMENSION]={0,0,0};
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];			
			p_n0ij[A_X]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_X];
			p_n0ij[A_Y]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_Y];
			p_n0ij[A_Z]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_Z];

			HYPER1[i*h_num+k].n0ij[A_X]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_X][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_X][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_X][2]*HYPER1[i*h_num+k].aiin[2]);
			HYPER1[i*h_num+k].n0ij[A_Y]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_Y][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_Y][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_Y][2]*HYPER1[i*h_num+k].aiin[2]);
			HYPER1[i*h_num+k].n0ij[A_Z]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_Z][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_Z][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_Z][2]*HYPER1[i*h_num+k].aiin[2]);
		}
		HYPER1[i*h_num+i].n0ij[A_X]=V*(HYPER[i].t_inverse_Ai[A_X][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_X][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_X][2]*p_n0ij[2]);
		HYPER1[i*h_num+i].n0ij[A_Y]=V*(HYPER[i].t_inverse_Ai[A_Y][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_Y][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_Y][2]*p_n0ij[2]);
		HYPER1[i*h_num+i].n0ij[A_Z]=V*(HYPER[i].t_inverse_Ai[A_Z][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_Z][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_Z][2]*p_n0ij[2]);
	}
		
	////DgDqの計算
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);		
//			cout<<"DgDq["<<k<<","<<i<<"]="<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
//		cout<<"DgDq["<<i<<","<<i<<"]="<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	}	
//	cout<<"----------OK"<<endl;


 }

void previous_strage(vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num)
{
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].q_n[A_X]=PART[i].r[A_X];	HYPER[i].q_n[A_Y]=PART[i].r[A_Y];	HYPER[i].q_n[A_Z]=PART[i].r[A_Z];
		HYPER[i].p_n[A_X]=HYPER[i].p[A_X];	HYPER[i].p_n[A_Y]=HYPER[i].p[A_Y];	HYPER[i].p_n[A_Z]=HYPER[i].p[A_Z];
		HYPER[i].ph_n[A_X]=HYPER[i].half_p[A_X];	HYPER[i].ph_n[A_Y]=HYPER[i].half_p[A_Y];	HYPER[i].ph_n[A_Z]=HYPER[i].half_p[A_Z];
		HYPER[i].W_n=HYPER[i].W;
		HYPER[i].stress_n[A_X][A_X]=HYPER[i].stress[A_X][A_X];	HYPER[i].stress_n[A_X][A_Y]=HYPER[i].stress[A_X][A_Y];	HYPER[i].stress_n[A_X][A_Z]=HYPER[i].stress[A_X][A_Z];
		HYPER[i].stress_n[A_Y][A_X]=HYPER[i].stress[A_Y][A_X];	HYPER[i].stress_n[A_Y][A_Y]=HYPER[i].stress[A_Y][A_Y];	HYPER[i].stress_n[A_Y][A_Z]=HYPER[i].stress[A_Y][A_Z];
		HYPER[i].stress_n[A_Z][A_X]=HYPER[i].stress[A_Z][A_X];	HYPER[i].stress_n[A_Z][A_Y]=HYPER[i].stress[A_Z][A_Y];	HYPER[i].stress_n[A_Z][A_Z]=HYPER[i].stress[A_Z][A_Z];
		HYPER[i].pi_n[A_X][A_X]=HYPER[i].pi[A_X][A_X];	HYPER[i].pi_n[A_X][A_Y]=HYPER[i].pi[A_X][A_Y];	HYPER[i].pi_n[A_X][A_Z]=HYPER[i].pi[A_X][A_Z];
		HYPER[i].pi_n[A_Y][A_X]=HYPER[i].pi[A_Y][A_X];	HYPER[i].pi_n[A_Y][A_Y]=HYPER[i].pi[A_Y][A_Y];	HYPER[i].pi_n[A_Y][A_Z]=HYPER[i].pi[A_Y][A_Z];
		HYPER[i].pi_n[A_Z][A_X]=HYPER[i].pi[A_Z][A_X];	HYPER[i].pi_n[A_Z][A_Y]=HYPER[i].pi[A_Z][A_Y];	HYPER[i].pi_n[A_Z][A_Z]=HYPER[i].pi[A_Z][A_Z];
		
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq_n[A_X]=HYPER1[k*h_num+i].DgDq[A_X];	HYPER1[k*h_num+i].DgDq_n[A_Y]=HYPER1[k*h_num+i].DgDq[A_Y];	HYPER1[k*h_num+i].DgDq_n[A_Z]=HYPER1[k*h_num+i].DgDq[A_Z];
		}
		HYPER1[i*h_num+i].DgDq_n[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	HYPER1[i*h_num+i].DgDq_n[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	HYPER1[i*h_num+i].DgDq_n[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
	}
}

//圧力など粒子の持つ情報をコンター図で表示する関数
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	int h_num=HYPER.size();
	double TIME=CON.get_step()*CON.get_dt();
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_interval()+1;		//出力する総ステップ数
	int step;
	int n=3;	//出力次元
	double mi=get_volume(&CON)*CON.get_hyper_density();

	if(CON.get_interval()==1)	step=t/CON.get_interval();
	else step=t/CON.get_interval()+1;
	

	if(t==1) 
	{
		ofstream fp("momentum.inp", ios::trunc);			
		fp<<STEP<<endl;//総ステップ数
		fp<<"data_geom"<<endl;
		fp.close();
	}

	//mainファイル書き込み
	ofstream fp("momentum.inp",ios :: app);
	fp<<"step"<<step<<" TIME="<<TIME<<endl;

	//fp<<step<<endl;
	
	//fp<<"data_geom"<<endl;
	//fp<<"step1"<<endl;
	//fp<<"step"<<t/CON->get_interval()+1<<" TIME="<<TIME<<endl;
	
	if(n==3)
	{
		fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力
	
		//節点番号とその座標の出力 
		for(int i=0;i<h_num;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
		//要素番号と要素形状の種類、そして要素を構成する節点番号出力
		for(int i=0;i<h_num;i++)	fp<<i<<"  0 pt "<<i<<endl;

		//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
		fp<<"10 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
		fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		fp<<"p_x,"<<endl;
		fp<<"p_y,"<<endl;
		fp<<"p_z,"<<endl;
		fp<<"lambda,"<<endl;
		fp<<"F_x,"<<endl;
		fp<<"F_y,"<<endl;
		fp<<"F_z,"<<endl;
		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;
		//fp<<"P,N/m^2"<<endl;
		//fp<<"value1,??"<<endl;

		//各節点の情報値入力
		for(int i=0;i<h_num;i++)
		{
			fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<" "<<HYPER[i].vis_force[A_X]<<" "<<HYPER[i].vis_force[A_Y]<<" "<<HYPER[i].vis_force[A_Z]<<endl;
			//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
		}
	}
	else if(n==2)
	{
		int num=0;
		int *Nxz=new int [h_num];
		fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力
	
		//節点番号とその座標の出力 
		for(int i=0;i<h_num;i++)
		{
			if(PART[i].q0[A_Y]<1/2*le&&PART[i].q0[A_Y]>-1/2*le)
			{
				fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				Nxz[num]=i;
				num++;
			}
		}
		//要素番号と要素形状の種類、そして要素を構成する節点番号出力
		for(int i=0;i<num;i++)	fp<<i<<"  0 pt "<<i<<endl;

		//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
//		fp<<"10 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
//		fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		fp<<"7 0"<<endl;
		fp<<"7 1 1 1 1 1 1 1"<<endl;
		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;
		fp<<"lambda,"<<endl;
		fp<<"F_x,"<<endl;
		fp<<"F_y,"<<endl;
		fp<<"F_z,"<<endl;
/*		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;*/
		//fp<<"P,N/m^2"<<endl;
		//fp<<"value1,??"<<endl;

		//各節点の情報値入力
		for(int i=0;i<num;i++)
		{
			int j=Nxz[i];
			fp<<j<<" "<<1/mi*HYPER[j].p[A_X]<<" "<<1/mi*HYPER[j].p[A_Y]<<" "<<1/mi*HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<endl;
			//fp<<j<<" "<<HYPER[j].p[A_X]<<" "<<HYPER[j].p[A_Y]<<" "<<HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<" "<<HYPER[j].vjs_force[A_X]<<" "<<HYPER[j].vjs_force[A_Y]<<" "<<HYPER[j].vjs_force[A_Z]<<endl;
			//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
		}
		delete[] Nxz;
	}
	fp.close();
}

void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int t)
{
	int h_num=HYPER.size();


	ofstream j("J.csv", ios::app);
	ofstream lam("lambda.csv", ios::app);
//	ofstream p("P.csv", ios::app);
	ofstream p_an("P_ave_norm.csv", ios::app);
	ofstream h("model_height.csv", ios::app);

//	if(t==1)	p<<"h_p,,,d_p,,,p\n";
//	p<<"t"<<t<<endl;
	double sum_lam=0;
	double sum_j=0;
	double p_sum=0,d_p_sum=0,h_p_sum=0;
	double h_min=PART[0].r[A_Z];
	double h_max=PART[0].r[A_Z];
	for(int i=0;i<h_num;i++)
	{
		lam<<HYPER[i].lambda<<",";
		j<<HYPER[i].J<<",";
		sum_lam+=HYPER[i].lambda;
		sum_j+=HYPER[i].J;
//		p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		p_sum+=sqrt(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z]);
		d_p_sum+=sqrt(HYPER[i].differential_p[A_X]*HYPER[i].differential_p[A_X]+HYPER[i].differential_p[A_Y]*HYPER[i].differential_p[A_Y]+HYPER[i].differential_p[A_Z]*HYPER[i].differential_p[A_Z]);
		h_p_sum+=sqrt(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z]);
		if(h_min>PART[i].r[A_Z])	h_min=PART[i].r[A_Z];
		if(h_max<PART[i].r[A_Z])	h_max=PART[i].r[A_Z];
	}
	
	lam<<sum_lam/h_num<<endl;
	j<<sum_j/h_num<<endl;
//	p<<endl;
	if(t==1)
	{
		p_an<<"p_sum,d_p_sum,h_p_sum\n";
		h<<"h_min,h_max\n";
	}
	p_an<<p_sum/h_num<<","<<d_p_sum/h_num<<","<<h_p_sum/h_num<<endl;
	h<<h_min<<","<<h_max<<endl;

	lam.close();
	j.close();
//	p.close();
	p_an.close();
	h.close();

	/*if(t==1)
	{
		////計算した各定数の出力
		ofstream ai("Ai.csv");
		ofstream inai("inverse_Ai.csv");
		ofstream aiin("aiin.csv");
		ofstream n0("noij.csv");
		ofstream wiin("wiin.csv");

		aiin<<"i"<<","<<"j"<<","<<"x"<<","<<"y"<<","<<"z"<<endl;
		n0<<"i"<<","<<"j"<<","<<"x"<<","<<"y"<<","<<"z"<<endl;

		for(int j=0;j<h_num;j++)	wiin<<j<<",";
		wiin<<endl;

		for(int i=0;i<h_num;i++)
		{
			inai<<i<<",";
			inai<<HYPER[i].inverse_Ai[0][A_X]<<","<<HYPER[i].inverse_Ai[0][A_Y]<<","<<HYPER[i].inverse_Ai[0][A_Z]<<","<<endl;		
			inai<<","<<HYPER[i].inverse_Ai[1][A_X]<<","<<HYPER[i].inverse_Ai[1][A_Y]<<","<<HYPER[i].inverse_Ai[1][A_Z]<<","<<endl;	
			inai<<","<<HYPER[i].inverse_Ai[2][A_X]<<","<<HYPER[i].inverse_Ai[2][A_Y]<<","<<HYPER[i].inverse_Ai[2][A_Z]<<","<<endl;	

			ai<<i<<",";
			ai<<HYPER[i].Ai[0][A_X]<<","<<HYPER[i].Ai[0][A_Y]<<","<<HYPER[i].Ai[0][A_Z]<<","<<endl;
			ai<<","<<HYPER[i].Ai[1][A_X]<<","<<HYPER[i].Ai[1][A_Y]<<","<<HYPER[i].Ai[1][A_Z]<<","<<endl;
			ai<<","<<HYPER[i].Ai[2][A_X]<<","<<HYPER[i].Ai[2][A_Y]<<","<<HYPER[i].Ai[2][A_Z]<<","<<endl;
			
			for(int j=0;j<h_num;j++)
			{
				aiin<<i<<","<<j<<","<<HYPER1[i*h_num+j].aiin[A_X]<<","<<HYPER1[i*h_num+j].aiin[A_Y]<<","<<HYPER1[i*h_num+j].aiin[A_Z]<<","<<endl;
				n0<<i<<","<<j<<","<<HYPER1[i*h_num+j].n0ij[A_X]<<","<<HYPER1[i*h_num+j].n0ij[A_Y]<<","<<HYPER1[i*h_num+j].n0ij[A_Z]<<","<<endl;
				wiin<<i<<","<<HYPER1[i*h_num+j].wiin<<",";
			}
			wiin<<endl;
		}
		ai.close();
		inai.close();
		n0.close();
		wiin.close();
	}*/

/*	ofstream d_p("d_P.csv", ios::app);
	ofstream h_p("h_P.csv", ios::app);*/
	stringstream ss_dgdq;
	ss_dgdq<<"./DgDq/DgDq"<<t<<".csv";
	ofstream dg(ss_dgdq.str());
	stringstream ss_n0;
	ss_n0<<"./n0/n0"<<t<<".csv";
	ofstream fn0(ss_n0.str());

	dg<<"t,"<<t<<endl;
	fn0<<"t,"<<t<<endl;
	for(int i=0;i<h_num;i++)
	{
		int N=HYPER[i].N;
		dg<<",,"<<i<<endl;
		fn0<<",,"<<i<<endl;
		for(int j=0;j<N;j++)
		{
			int nei=HYPER[i].NEI[j];
			dg<<",,,"<<nei<<","<<HYPER1[i*h_num+nei].DgDq[A_X]<<","<<HYPER1[i*h_num+nei].DgDq[A_Y]<<","<<HYPER1[i*h_num+nei].DgDq[A_Z]<<endl;
			fn0<<",,,"<<nei<<","<<HYPER1[nei*h_num+i].n0ij[A_X]<<","<<HYPER1[nei*h_num+i].n0ij[A_Y]<<","<<HYPER1[nei*h_num+i].n0ij[A_Z]<<endl;
		}
		dg<<",,,"<<i<<","<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
		fn0<<",,,"<<i<<","<<HYPER1[i*h_num+i].n0ij[A_X]<<","<<HYPER1[i*h_num+i].n0ij[A_Y]<<","<<HYPER1[i*h_num+i].n0ij[A_Z]<<endl;
	}
	dg.close();
	fn0.close();

	if(t==1)
	{
		ofstream fsdg("stress_DgDq.csv",ios::trunc);
		ofstream fpi("pi.csv",ios::trunc);
		ofstream fpn("pi_n0.csv",ios::trunc);
		ofstream fst("stress.csv",ios::trunc);
		fsdg.close();
		fpn.close();
		fpi.close();
		fst.close();
	}

	ofstream fsdg("stress_DgDq.csv",ios::app);
	ofstream fpi("pi.csv",ios::app);
	ofstream fst("stress.csv",ios::app);
	ofstream fpn("pi_n0.csv",ios::app);


	fsdg<<"t,"<<t<<endl;
	fpn<<"t,"<<t<<endl;
	fpi<<"t,"<<t<<endl;
	fst<<"t,"<<t<<endl;
	double sdgdq[DIMENSION]={0,0,0};
	double pin0[DIMENSION]={0,0,0};
	for(int i=0;i<h_num;i++)
	{
		int N=HYPER[i].N;
		fsdg<<",,"<<i<<endl;
		fpn<<",,"<<i<<endl;
		for(int j=0;j<N;j++)
		{
			int jn=HYPER[i].NEI[j];
			sdgdq[A_X]=HYPER[jn].stress[A_X][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[A_X][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[A_X][2]*HYPER1[jn*h_num+i].DgDq[2];
			sdgdq[A_Y]=HYPER[jn].stress[A_Y][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[A_Y][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[A_Y][2]*HYPER1[jn*h_num+i].DgDq[2];
			sdgdq[A_Z]=HYPER[jn].stress[A_Z][0]*HYPER1[jn*h_num+i].DgDq[0]+HYPER[jn].stress[A_Z][1]*HYPER1[jn*h_num+i].DgDq[1]+HYPER[jn].stress[A_Z][2]*HYPER1[jn*h_num+i].DgDq[2];

			fsdg<<",,,"<<jn<<","<<sdgdq[A_X]<<","<<sdgdq[A_Y]<<","<<sdgdq[A_Z]<<endl;

			pin0[A_X]=HYPER[jn].pi[A_X][0]*HYPER1[i*h_num+jn].n0ij[0]+HYPER[jn].pi[A_X][1]*HYPER1[i*h_num+jn].n0ij[1]+HYPER[jn].pi[A_X][2]*HYPER1[i*h_num+jn].n0ij[2];
			pin0[A_Y]=HYPER[jn].pi[A_Y][0]*HYPER1[i*h_num+jn].n0ij[0]+HYPER[jn].pi[A_Y][1]*HYPER1[i*h_num+jn].n0ij[1]+HYPER[jn].pi[A_Y][2]*HYPER1[i*h_num+jn].n0ij[2];
			pin0[A_Z]=HYPER[jn].pi[A_Z][0]*HYPER1[i*h_num+jn].n0ij[0]+HYPER[jn].pi[A_Z][1]*HYPER1[i*h_num+jn].n0ij[1]+HYPER[jn].pi[A_Z][2]*HYPER1[i*h_num+jn].n0ij[2];
			fpn<<",,,"<<jn<<","<<pin0[A_X]<<","<<pin0[A_Y]<<","<<pin0[A_Z]<<endl;
		}
		sdgdq[A_X]=HYPER[i].stress[A_X][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+i].DgDq[2];
		sdgdq[A_Y]=HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+i].DgDq[2];
		sdgdq[A_Z]=HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+i].DgDq[2];
		fsdg<<",,,"<<i<<","<<sdgdq[A_X]<<","<<sdgdq[A_Y]<<","<<sdgdq[A_Z]<<endl;

		pin0[A_X]=HYPER[i].pi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi[A_X][2]*HYPER1[i*h_num+i].n0ij[2];
		pin0[A_Y]=HYPER[i].pi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2];
		pin0[A_Z]=HYPER[i].pi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2];
		fpn<<",,,"<<i<<","<<pin0[A_X]<<","<<pin0[A_Y]<<","<<pin0[A_Z]<<endl;
	
		fpi<<",,"<<i<<","<<HYPER[i].pi[A_X][A_X]<<","<<HYPER[i].pi[A_X][A_Y]<<","<<HYPER[i].pi[A_X][A_Z]<<endl;
		fpi<<",,,"<<HYPER[i].pi[A_Y][A_X]<<","<<HYPER[i].pi[A_Y][A_Y]<<","<<HYPER[i].pi[A_Y][A_Z]<<endl;
		fpi<<",,,"<<HYPER[i].pi[A_Z][A_X]<<","<<HYPER[i].pi[A_Z][A_Y]<<","<<HYPER[i].pi[A_Z][A_Z]<<endl;
		fst<<",,"<<i<<","<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<endl;
		fst<<",,,"<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<endl;
		fst<<",,,"<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<endl;
	}
	fsdg.close();
	fpi.close();
	fpn.close();
	fst.close();

//
//	ofstream r("r_CG.csv", ios::app);
//	double r_cg_x=0;
//	double r_cg_y=0;
//	double r_cg_z=0;
//
//	for(int i=0;i<h_num;i++)
//	{
//		r_cg_x+=PART[i].r[A_X];
//		r_cg_y+=PART[i].r[A_Y];
//		r_cg_z+=PART[i].r[A_Z];
//	}
//
//	r<<r_cg_x/h_num<<","<<r_cg_y/h_num<<","<<r_cg_z/h_num<<endl;
//	r.close();

/*	ofstream J("J.csv", ios::app);
	ofstream ti_Fi("ti_Fi.csv", ios::app);
	ofstream Fi("Fi.csv", ios::app);*/
	/*
	ofstream lam("lambda.csv", ios::app);
	lam<<HYPER[182].lambda<<","<<HYPER[911].lambda<<","<<HYPER[1640].lambda<<","<<HYPER[2369].lambda<<endl;
	lam.close();

	ofstream j("J.csv", ios::app);
	j<<HYPER[182].J<<","<<HYPER[911].J<<","<<HYPER[1640].J<<","<<HYPER[2369].J<<endl;
	j.close();


	ofstream p("P.csv", ios::app);
	p<<sqrt(HYPER[182].p[A_X]*HYPER[182].p[A_X]+HYPER[182].p[A_Z]*HYPER[182].p[A_Z])<<",";
	p<<sqrt(HYPER[911].p[A_X]*HYPER[911].p[A_X]+HYPER[911].p[A_Z]*HYPER[911].p[A_Z])<<",";
	p<<sqrt(HYPER[1640].p[A_X]*HYPER[1640].p[A_X]+HYPER[1640].p[A_Z]*HYPER[1640].p[A_Z])<<",";	
	p<<sqrt(HYPER[2369].p[A_X]*HYPER[2369].p[A_X]+HYPER[2369].p[A_Z]*HYPER[2369].p[A_Z])<<endl;
	p.close();

	/*
	if(t==1)
	{
		p<<"t"<<",";		
		d_p<<"t"<<",";
		h_p<<"t"<<",";/
		lam<<"t"<<",";
		stress<<"t"<<",";
		ti_Fi<<"t"<<",";
		Fi<<"t"<<",";
		J<<"t"<<",";//

		for(int i=0;i<h_num;i++)
		{
			p<<i<<","<<","<<",";
			lam<<i<<",";
			d_p<<i<<","<<","<<",";
			h_p<<i<<","<<","<<",";
			stress<<i<<","<<","<<",";
			ti_Fi<<i<<","<<","<<",";
			Fi<<i<<","<<","<<",";
			J<<i<<",";
		}
		p<<endl;
		lam<<endl;
		d_p<<endl;
		h_p<<endl;
		stress<<endl;
		ti_Fi<<endl;
		Fi<<endl;
		J<<endl;
	}

	p<<t<<",";
	lam<<t<<",";
	d_p<<t<<",";
	h_p<<t<<",";
	stress<<t<<",";
	ti_Fi<<t<<",";
	Fi<<t<<",";
	J<<t<<",";
	for(int i=0;i<h_num;i++)
	{
		
		p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<",";
		lam<<HYPER[i].lambda<<",";
/*		d_p<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<endl;
		h_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<endl;

		stress<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_X][A_X]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Z]<<",";	
		Fi<<HYPER[i].Fi[0][0]<<","<<HYPER[i].Fi[0][1]<<","<<HYPER[i].Fi[0][2]<<",";
		J<<HYPER[i].J<<",";
		for(int j=0;j<h_num;j++)	dg<<i<<","<<j<<","<<HYPER1[j*h_num+i].DgDq[A_X]<<","<<HYPER1[j*h_num+i].DgDq[A_Y]<<","<<HYPER1[j*h_num+i].DgDq[A_Z]<<","<<endl;
		
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num;i++)
	{
		stress<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Y][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Z]<<",";
		Fi<<HYPER[i].Fi[1][0]<<","<<HYPER[i].Fi[1][1]<<","<<HYPER[i].Fi[1][2]<<",";
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num;i++)
	{
		stress<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Z][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Z]<<",";
		Fi<<HYPER[i].Fi[2][0]<<","<<HYPER[i].Fi[2][1]<<","<<HYPER[i].Fi[2][2]<<",";
	}

	p<<endl;
	lam<<endl;
	d_p<<endl;
	h_p<<endl;
	stress<<endl;
	ti_Fi<<endl;
	Fi<<endl;
	J<<endl;

	dg.close();
	stress.close();
	ti_Fi.close();
	Fi.close();
	d_p.close();
	h_p.close();
	J.close();
	p.close();
	lam.close();*/
}

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t)
{
//	cout<<"弾性ポテンシャル計算";
	int h_num=HYPER.size();
	int p_num=PART.size();
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	vector<double>	W;
	W.reserve(h_num);
	for(int i=0;i<h_num;i++)	W.emplace_back(0);
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double G=9.8;
	/*
	for(int j=0;j<4;j++)
	{
		if(j==0)	i=182;
		else if(j==1) i=911;
		else if(j==2) i=1640;
		else if(j==3) i=2369;

		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		W[i]=c10*(Ic-3)+c01*(IIc-3);
	}*/

	
	for(int i=0;i<h_num;i++)
	{
		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		W[i]=c10*(Ic-3)+c01*(IIc-3);
	}//*/
//	cout<<"----------OK"<<endl;

/*	ofstream e("E.csv", ios::app);
	ofstream e_T("E_T.csv", ios::app);
	ofstream e_g("E_g.csv", ios::app);
	ofstream e_W("E_W.csv", ios::app);
	ofstream e_lam("E_lam.csv", ios::app);//*/
	ofstream e_sum("E_sum.csv", ios::app);

	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
/*	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int j=0;j<4;j++)
	{
		if(j==0)	i=182;
		else if(j==1) i=911;
		else if(j==2) i=1640;
		else if(j==3) i=2369;

		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		//energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*9.8*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";
	}*/
	
/*	if(t==1)
	{
		e<<"t"<<",";
		e_T<<"t"<<",";
		e_g<<"t"<<",";
		e_W<<"t"<<",";
		e_lam<<"t"<<",";
		for(int i=0;i<h_num;i++)
		{
			e<<i<<",";
			e_T<<i<<",";
			e_g<<i<<",";
			e_W<<i<<",";
			e_lam<<i<<",";
		}
		e<<endl;
		e_T<<endl;
		e_g<<endl;
		e_W<<endl;
		e_lam<<endl;
	}

	e<<t<<",";
	e_T<<t<<",";
	e_g<<t<<",";
	e_W<<t<<",";
	e_lam<<t<<",";*/

	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int i=0;i<h_num;i++)
	{
		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		//energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		if(CON.get_flag_G()==ON)	energy=0.5/mi*vv+mi*G*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		if(CON.get_flag_G()==OFF)	energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_T+=0.5/mi*vv;
		sum_e_g+=mi*G*PART[i].r[A_Z];
		sum_e_lam+=HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_W+=W[i]*V;
		sum_e+=energy;
	/*	e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*G*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";//*/
	}
/*	e<<endl;
	e_T<<endl;
	e_g<<endl;
	e_W<<endl;
	e_lam<<endl;*/

	
	if(t==1)	e_sum<<"E,E_T,E_g,E_W,E_lam\n";
	e_sum<<sum_e<<","<<sum_e_T<<","<<sum_e_g<<","<<sum_e_W<<","<<sum_e_lam<<endl;//*/

/*	e.close();
	e_T.close();
	e_g.close();
	e_W.close();
	e_lam.close();//*/
	e_sum.close();
}

double calc_det3(double **M)
{
	double det=0;
	det+=M[0][0]*M[1][1]*M[2][2]+M[0][1]*M[1][2]*M[2][0]+M[0][2]*M[1][0]*M[2][1];
	det-=M[0][2]*M[1][1]*M[2][0]+M[0][1]*M[1][0]*M[2][2]+M[0][0]*M[1][2]*M[2][1];
	return det;
}

void inverse(double **a,int N)
{
	double d=0;
	double *col=new double [N];
	double **y=new double *[N];
	for(int i=0;i<N;i++)	y[i]=new double [N];
	int *index=new int [N];

	//int p_revision=0,m_revision=0;
	int fig=0,x=0;
	int	flag=OFF;

	for(int D=0;D<N;D++)
	{
		col[D]=0;
		index[D]=0;
		for(int D2=0;D2<N;D2++)	y[D][D2]=0;
	}

	for(int D=0;D<N;D++)	for(int D2=0;D2<N;D2++)	if(fabs(a[D][D2])<DBL_EPSILON)	a[D][D2]=0;

	for(int D=0;D<N;D++)
	{
		for(int D2=0;D2<N;D2++)
		{
			if(fabs(a[D][D2])>pow(10.0,5.0))
			{
				flag=ON;
				x=log10(fabs(a[D][D2]));
				if(x>fig)	fig=x;
			}
		}
	}	
	
	if(flag==ON)	for(int D=0;D<N;D++) for(int D2=0;D2<N;D2++) a[D][D2]/=pow(10.0,fig);

	ludcmp(a,N,index,&d);
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)	col[i]=0.0;
		col[j]=1.0;
		lubksb(a,N,index,col);
		for(int i=0;i<N;i++)	y[i][j]=col[i];
	}

	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	a[i][j]=y[i][j];
	if(flag==ON)	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	a[i][j]/=pow(10.0,fig);

	for(int i=0;i<N;i++)	delete[]	y[i];
	delete[]	y;
	delete[]	col;
	delete[]	index;
}

void ludcmp(double **a,int N,int *index,double *d)
{
	int imax=0;
	int count=0;
	double max=0,temp=0,sum=0,dum=0;
	*d=1.0;
	vector<float> buf(N);

	for(int i=0;i<N;i++)	buf[i]=0;

	for(int i=0;i<N;i++)
	{
		max=0.0;
		for(int j=0;j<N;j++)
		{
			temp=fabs(a[i][j]);
			if(temp>max)	max=temp;
		}
		if(max==0.0&&count==0)
		{
			count++;
			cout<<"特異行列である"<<endl;
		}
		buf[i]=1.0/max;
	}

	for(int j=0;j<N;j++)
	{
		for(int i=0;i<j;i++)
		{
			sum=a[i][j];
			for(int k=0;k<i;k++)	sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		max=0.0;
		for(int i=j;i<N;i++)
		{
			sum=a[i][j];
			for(int k=0;k<j;k++)	sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
			dum=buf[i]*fabs(sum);
			if(dum>=max)
			{
				max=dum;
				imax=i;
			}
		}

		if(j!=imax)
		{
			for(int k=0;k<N;k++)
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d=-(*d);
			buf[imax]=buf[j];
		}
		index[j]=imax;
		if(fabs(a[j][j])<DBL_EPSILON)
		{
			cout<<"error:a["<<j<<"]["<<j<<"]=0"<<endl;
			a[j][j]=0;
		}
		if(j!=N-1)
		{
			dum=1.0/(a[j][j]);
			for(int i=j+1;i<N;i++)a[i][j]*=dum;
		}
	}
	buf.clear();
}

void lubksb(double **a,int N,int *index,double b[])
{
	int ip=0,ii=0,count=0;
	double sum=0;

	for(int i=0;i<N;i++)
	{
		ip=index[i];
		sum=b[ip];
		b[ip]=b[i];
		if(count)	for(int j=ii;j<i;j++)	sum-=a[i][j]*b[j];
		else if(sum)
		{
			count++;
			ii=i;
		}
		b[i]=sum;
	}

	for(int i=N-1;i>=0;i--)
	{
		sum=b[i];
		for(int j=i+1;j<N;j++)	sum-=a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
