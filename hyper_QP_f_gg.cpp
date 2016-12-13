#include "stdafx.h"	

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0);
void q_nab_lap(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, double *dg, double Dt, double V, double mi,double nG[DIMENSION]);
void q_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_n(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void p_QP(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION]);
void p_variables(vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int h_num,double Dt,double mi);
void p_nab(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *df, double *dh, double Dt, double V, double mi,double nG[DIMENSION]);
void output_data_p(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double *dG, double *G,double *f, double *dh, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E);
//void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT, double *dN,double *g, double **dg, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E);
void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double *df, double *f,double *g, double *dg, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E);


void calc_HYPER_QP_gh(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{
	////////////íËã`///////////////
	int h_num=HYPER.size();
	int Nx=h_num*2;

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
	double nG[DIMENSION]={0,0,1.};

	double E0=0;

	cout<<"QP start-------------------------"<<endl;
	ofstream fq0("q0_QP.csv");
	ofstream fqn0("qn0_QP.csv");
	ofstream fpn0("pn0_QP.csv");
	ofstream fhpn0("hpn0_QP.csv");
	for(int i=0;i<h_num;i++)
	{
		fq0<<HYPER[i].q_n[A_X]<<","<<HYPER[i].q_n[A_Y]<<","<<HYPER[i].q_n[A_Z]<<endl;
		fqn0<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		fpn0<<HYPER[i].p_n[A_X]<<","<<HYPER[i].p_n[A_Y]<<","<<HYPER[i].p_n[A_Z]<<endl;
		fhpn0<<HYPER[i].ph_n[A_X]<<","<<HYPER[i].ph_n[A_Y]<<","<<HYPER[i].ph_n[A_Z]<<endl;
	}
	fq0.close();
	fqn0.close();
	fpn0.close();
	fhpn0.close();

	////////////QPåvéZ///////////////		
	q_QP(CON,PART,HYPER,HYPER1,t,h_num,Nx,V,mi,Dt,nG,E0);
	for(int i=0;i<h_num;i++)	HYPER[i].lambda=HYPER[i].lam;
	calc_n(CON,PART,HYPER,HYPER1);

	p_QP(HYPER,HYPER1,t,h_num,Nx,V,mi,Dt,nG);

	cout<<"--------------------------OK"<<endl;


}

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=10000;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1e-5;
	double ep_min=1e-5;
	double d_ep=1e-20;

	double p_p[DIMENSION]={0,0,0};

	double r=1000.;

	double En=0;

	double *f=new double [h_num];	
	double *df=new double [h_num*h_num];

	double *fx=new double [h_num];	
	double *dfx=new double [h_num*h_num];

	double *g=new double [h_num];
	double *h=new double [h_num];

	double *dg=new double [h_num*h_num];
	double *dh=new double [h_num*h_num];

	double *th_h=new double [h_num];

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		f[i]=0.;
		fx[i]=0.;
		g[i]=0.;
		h[i]=0.;
		th_h[i]=0.;
		for(int j=0;j<h_num;j++)
		{
			dg[i*h_num+j]=0.;
			dh[i*h_num+j]=0.;
			df[i*h_num+j]=0.;
			dfx[i*h_num+j]=0.;
		}
	}

	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			dh[in*h_num+k]=(double)0.5*Dt*Dt/mi*(HYPER1[k*h_num+in].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+in].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+in].DgDq_n[A_Z]*nG[A_Z]);
		}
		dh[k*h_num+k]=(double)0.5*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+k].DgDq_n[A_Z]*nG[A_Z]);
	}

	stringstream ss0;
	ss0<<"dh"<<t<<".csv";
	ofstream f_dh(ss0.str());
	for(int i=0;i<h_num;i++)
	{
		for(int k=0;k<h_num;k++)	f_dh<<dh[i*h_num+k]<<",";
		f_dh<<endl;
	}
	f_dh.close();


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
			q_nab_lap(PART,HYPER,HYPER1,dg,Dt,V,mi,nG);


			for(int i=0;i<h_num;i++)
			{
				g[i]=(double)V*(1-HYPER[i].J);
				
				//cout<<"f"<<count<<"="<<f[i]<<endl;
			}
			for(int i=0;i<h_num;i++)	f[i]=(double)g[i];
			//cout<<"T="<<T<<", En"<<En<<endl;
			for(int i=0;i<h_num;i++)
			{
				for(int j=0;j<h_num;j++)
				{
					df[i*h_num+j]=dg[i*h_num+j];
				}				
				//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			}
			////for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;
	
			for(int i=0;i<h_num;i++)
			{	
				h[i]=-1*(PART[i].r[A_X]*nG[A_X]+PART[i].r[A_Y]*nG[A_Y]+PART[i].r[A_Z]*nG[A_Z]);
				//cout<<"h"<<count<<"="<<h[i]<<endl;

				if(h[i]+th_h[i]>0)
				{
					f[i]+=0.5*r*(h[i]+th_h[i])*(h[i]+th_h[i]);
					cout<<"	h"<<count<<"="<<h[i]<<endl;
				}
				//cout<<"h="<<h[i]<<endl;
				for(int j=0;j<h_num;j++)
				{
					//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
					if(h[i]+th_h[i]>0)	df[i*h_num+j]+=r*dh[i*h_num+j]*(h[i]+th_h[i]);
				}
				//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			}
			//for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			//stringstream ssf;
			//ssf<<"f"<<count_min<<"_"<<count<<".csv";
			//ofstream ff(ssf.str());
			//stringstream ssdf;
			//ssdf<<"df"<<count_min<<"_"<<count<<".csv";
			//ofstream fdf(ssdf.str());
			//for(int i=0;i<h_num;i++)
			//{
			//	ff<<f[i]<<","<<HYPER[i].N_f<<endl;
			//	for(int j=0;j<h_num;j++)
			//	{
			//		fdf<<df[i*h_num+j]<<",";
			//	}		
			//	fdf<<endl;
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
			//fdf<<endl;
			//for(int i=0;i<h_num;i++)
			//{
			//	for(int j=0;j<h_num;j++)
			//	{
			//		fdf<<HYPER1[i*h_num+j].N_Df<<",";
			//	}		
			//	fdf<<endl;
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
			//ff.close();
			//fdf.close();


			//for(int i=0;i<h_num;i++)	fx[i]=f[i];
			////cout<<"T="<<T<<", En"<<En<<endl;
			//for(int i=0;i<h_num;i++)
			//{
			//	for(int j=0;j<h_num;j++)
			//	{
			//		dfx[i*h_num+j]=df[i*h_num+j];
			//	}				
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
	
			//for(int i=0;i<h_num;i++)
			//{
			//	fx[i]=g[i];
			//	for(int j=0;j<h_num;j++)
			//	{
			//		dfx[i*h_num+j]=HYPER1[i*h_num+j].N_Df;
			//	}		
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}

			gauss(df,f,h_num);
			//gauss(dfx,fx,h_num);

			//for(int i=0;i<h_num;i++)
			//{
			//	//cout<<"f"<<count<<"="<<f[i]<<", fx="<<fx[i]<<endl;
			//	cout<<"f"<<count<<"="<<f[i]<<endl;
			//}
			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=f[i];
				E_sum+=fabs(f[i]);
			}
			E=E_sum;
			///////////////pn1_2åvéZ
			//for(int i=0;i<h_num;i++)
			//{
			//	p_p[A_X]=0;	p_p[A_Y]=0;	p_p[A_Z]=0;	
			//	int Ni=HYPER[i].N;
			//	for(int j=0;j<Ni;j++)
			//	{
			//		int jn=HYPER[i].NEI[j];
			//		double DgDq_n[DIMENSION]={HYPER1[jn*h_num+i].DgDq_n[A_X], HYPER1[jn*h_num+i].DgDq_n[A_Y], HYPER1[jn*h_num+i].DgDq_n[A_Z]};
			//		p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*DgDq_n[A_Z];
			//		p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*DgDq_n[A_Z];
			//		p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*DgDq_n[A_Z];
			//	}
			//	double DgDq_nii[DIMENSION]={HYPER1[i*h_num+i].DgDq_n[A_X], HYPER1[i*h_num+i].DgDq_n[A_Y], HYPER1[i*h_num+i].DgDq_n[A_Z]};
			//	p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].lam)*DgDq_nii[A_X]+HYPER[i].stress_n[A_X][A_Y]*DgDq_nii[A_Y]+HYPER[i].stress_n[A_X][A_Z]*DgDq_nii[A_Z];
			//	p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*DgDq_nii[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].lam)*DgDq_nii[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*DgDq_nii[A_Z];
			//	p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*DgDq_nii[A_X]+HYPER[i].stress_n[A_Z][A_Y]*DgDq_nii[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].lam)*DgDq_nii[A_Z];
			//	HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
			//	HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
			//	HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];

			//	/////////////qåvéZ
			//	PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
			//	PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
			//	PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
			//}

			if(count==1||count%500==0)cout<<"E"<<count<<"="<<E<<endl;
			output_data(PART,HYPER,HYPER1,df,f,g,dg,h,th_h,Nx,h_num,count,count_min,t,E);
		}

		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=(HYPER[i].old_lam-HYPER[i].lam)*(HYPER[i].old_lam-HYPER[i].lam);
		E_min=sqrt(E_sum);
		if(E_min<ep_min*1000)	r*=4;
		cout<<"Emin"<<count_min<<"="<<E_min<<endl;

		for(int i=0;i<h_num;i++)
		{
			if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		}
	}




	///////////////pn1_2åvéZ
	//for(int i=0;i<h_num;i++)
	//{
	//	p_p[A_X]=0;	p_p[A_Y]=0;	p_p[A_Z]=0;	
	//	int Ni=HYPER[i].N;
	//	for(int j=0;j<Ni;j++)
	//	{
	//		int jn=HYPER[i].NEI[j];
	//		double DgDq_n2[DIMENSION]={HYPER1[jn*h_num+i].DgDq_n[A_X], HYPER1[jn*h_num+i].DgDq_n[A_Y], HYPER1[jn*h_num+i].DgDq_n[A_Z]};
	//		p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*DgDq_n2[A_X]+HYPER[jn].stress_n[A_X][A_Y]*DgDq_n2[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*DgDq_n2[A_Z];
	//		p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*DgDq_n2[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*DgDq_n2[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*DgDq_n2[A_Z];
	//		p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*DgDq_n2[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*DgDq_n2[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*DgDq_n2[A_Z];
	//	}
	//	double DgDq_nii2[DIMENSION]={HYPER1[i*h_num+i].DgDq_n[A_X], HYPER1[i*h_num+i].DgDq_n[A_Y], HYPER1[i*h_num+i].DgDq_n[A_Z]};
	//	p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].lam)*DgDq_nii2[A_X]+HYPER[i].stress_n[A_X][A_Y]*DgDq_nii2[A_Y]+HYPER[i].stress_n[A_X][A_Z]*DgDq_nii2[A_Z];
	//	p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*DgDq_nii2[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].lam)*DgDq_nii2[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*DgDq_nii2[A_Z];
	//	p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*DgDq_nii2[A_X]+HYPER[i].stress_n[A_Z][A_Y]*DgDq_nii2[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].lam)*DgDq_nii2[A_Z];
	//	HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
	//	HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
	//	HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];

	//	/////////////qåvéZ
	//	PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
	//	PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
	//	PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	//}

	delete[]	f;
	delete[]	df;
	delete[]	g;
	delete[]	h;
	delete[]	dg;
	delete[]	dh;
	delete[]	th_h;
	delete[]	dfx;
	delete[]	fx;

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


	/////////////pn1_2åvéZ
	for(int i=0;i<h_num;i++)
	{
		double p_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double DgDq_n[DIMENSION]={HYPER1[jn*h_num+i].DgDq_n[A_X], HYPER1[jn*h_num+i].DgDq_n[A_Y], HYPER1[jn*h_num+i].DgDq_n[A_Z]};
			p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*DgDq_n[A_Z];
			p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*DgDq_n[A_Z];
			p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*DgDq_n[A_Z];
		}
		double DgDq_nii[DIMENSION]={HYPER1[i*h_num+i].DgDq_n[A_X], HYPER1[i*h_num+i].DgDq_n[A_Y], HYPER1[i*h_num+i].DgDq_n[A_Z]};
		p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].lam)*DgDq_nii[A_X]+HYPER[i].stress_n[A_X][A_Y]				 *DgDq_nii[A_Y]+HYPER[i].stress_n[A_X][A_Z]				  *DgDq_nii[A_Z];
		p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]				*DgDq_nii[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].lam)*DgDq_nii[A_Y]+HYPER[i].stress_n[A_Y][A_Z]				  *DgDq_nii[A_Z];
		p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]				*DgDq_nii[A_X]+HYPER[i].stress_n[A_Z][A_Y]				 *DgDq_nii[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].lam)*DgDq_nii[A_Z];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];

		/////////////qåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}


	/////////////F, J, t_inverseÇÃåvéZ
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	double J=0.;
	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		int Ni=HYPER[i].N;	

		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+jn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+jn].aiin[A_X],	HYPER1[i*h_num+jn].aiin[A_Y],	HYPER1[i*h_num+jn].aiin[A_Z]};

			fi[0][0]+=w*(PART[jn].r[A_X]-PART[i].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[jn].r[A_X]-PART[i].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[jn].r[A_X]-PART[i].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[jn].r[A_Y]-PART[i].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[jn].r[A_Y]-PART[i].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[jn].r[A_Y]-PART[i].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[jn].r[A_Z]-PART[i].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[jn].r[A_Z]-PART[i].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[jn].r[A_Z]-PART[i].r[A_Z])*a[A_Z];
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

		J=calc_det3(p_Fi);
		HYPER[i].J=J;	

		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];

		for(int j=0;j<Ni;j++)
		{			
			int k=HYPER[i].NEI[j];
			HYPER1[i*h_num+k].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			HYPER1[i*h_num+k].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			HYPER1[i*h_num+k].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
			//cout<<"i"<<i<<"j"<<k<<"	"<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);
	}	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;


//cout<<"----------OK"<<endl;
}


void q_nab_lap(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dg, double Dt, double V, double mi,double nG[DIMENSION])
{

	int h_num=HYPER.size();
	int Nx=2*h_num;
	double dm=0.;
	/////////////////dEåvéZ///////////////////
	for(int k=0;k<h_num;k++)
	{
		for(int i=0;i<h_num;i++)	dg[i*h_num+k]=0.;
	}

	//for(int k=0;k<h_num;k++)
	//{
	//	int Nk=HYPER[k].N;
	//	for(int i=0;i<Nk;i++)
	//	{
	//		int in=HYPER[k].NEI[i];

	//		for(int j=0;j<Nk;j++)
	//		{
	//			int jn=HYPER[k].NEI[j];

	//			dg[in*h_num+jn]-=0.5*Dt*Dt/mi*(HYPER1[in*h_num+k].DgDq[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+HYPER1[in*h_num+k].DgDq[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+HYPER1[in*h_num+k].DgDq[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
	//		}
	//		dg[in*h_num+k]-=0.5*Dt*Dt/mi*(HYPER1[in*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]+HYPER1[in*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]+HYPER1[in*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);
	//		dg[k*h_num+in]-=0.5*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[in*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[in*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[in*h_num+k].DgDq_n[A_Z]);

	//	}
	//	dg[k*h_num+k]-=0.5*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);
	//}
	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		double DgDq_kk[DIMENSION]={HYPER1[k*h_num+k].DgDq[A_X], HYPER1[k*h_num+k].DgDq[A_Y], HYPER1[k*h_num+k].DgDq[A_Z]};
		double DgDq_nkk[DIMENSION]={HYPER1[k*h_num+k].DgDq_n[A_X], HYPER1[k*h_num+k].DgDq_n[A_Y], HYPER1[k*h_num+k].DgDq_n[A_Z]};
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];

			double DgDq_ik[DIMENSION]={HYPER1[in*h_num+k].DgDq[A_X], HYPER1[in*h_num+k].DgDq[A_Y], HYPER1[in*h_num+k].DgDq[A_Z]};
			double DgDq_nik[DIMENSION]={HYPER1[k*h_num+in].DgDq_n[A_X], HYPER1[k*h_num+in].DgDq_n[A_Y], HYPER1[k*h_num+in].DgDq_n[A_Z]};

			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				dg[in*h_num+jn]-=0.5*Dt*Dt/mi*(DgDq_ik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+DgDq_ik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+DgDq_ik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
			}
			dg[in*h_num+k]-=0.5*Dt*Dt/mi*(DgDq_ik[A_X]*DgDq_nkk[A_X]+DgDq_ik[A_Y]*DgDq_nkk[A_Y]+DgDq_ik[A_Z]*DgDq_nkk[A_Z]);
			dg[k*h_num+in]-=0.5*Dt*Dt/mi*(DgDq_kk[A_X]*DgDq_nik[A_X]+DgDq_kk[A_Y]*DgDq_nik[A_Y]+DgDq_kk[A_Z]*DgDq_nik[A_Z]);

		}
		dg[k*h_num+k]-=0.5*Dt*Dt/mi*(DgDq_kk[A_X]*DgDq_nkk[A_X]+DgDq_kk[A_Y]*DgDq_nkk[A_Y]+DgDq_kk[A_Z]*DgDq_nkk[A_Z]);
	}


}

void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double *df, double *f,double *g, double *dg, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E)
{
	stringstream ss1;
	ss1<<"./dT/df_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./lam/lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./g/g_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss7;
	ss7<<"./dg/dg_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./p/hp_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./q/q_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./dT/d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./E/E_"<<t<<"_"<<count_min<<".csv";
	stringstream ss12;
	ss12<<"./h/h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss13;
	ss13<<"./g/J_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss14;
	ss14<<"./g/DgDq_"<<t<<"_"<<count_min<<"_"<<count<<".csv";

	ofstream f_dt(ss1.str());
	ofstream f_lam(ss3.str());
	ofstream f_g(ss6.str());
	ofstream f_dg(ss7.str());
	ofstream f_p(ss8.str());
	ofstream f_q(ss9.str());
	ofstream f_d(ss10.str());
	ofstream f_h(ss12.str());
	ofstream f_J(ss13.str());
	ofstream f_dgdq(ss14.str());

	if(count==1)
	{
		ofstream f_E(ss11.str(), ios::trunc);
		f_E<<count<<","<<E<<endl;
		f_E.close();
	}	
	else
	{
		ofstream f_E(ss11.str(), ios::app);
		f_E<<count<<","<<E<<endl;
		f_E.close();
	}



	for(int i=0;i<h_num;i++)
	{				
		for(int j=0;j<h_num;j++)
		{
			f_dt<<df[i*h_num+j]<<",";
		}
		f_dt<<endl;
		f_d<<f[i]<<endl;
	}
	for(int i=0;i<h_num;i++)
	{
		f_lam<<HYPER[i].lam<<endl;
		f_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
		f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		f_dgdq<<i<<endl;

		for(int k=0;k<h_num;k++)
		{
			f_dg<<dg[i*h_num+k]<<",";
		}
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			f_dgdq<<","<<jn<<","<<HYPER1[i*h_num+jn].DgDq[A_X]<<","<<HYPER1[i*h_num+jn].DgDq[A_Y]<<","<<HYPER1[i*h_num+jn].DgDq[A_Z]<<endl;
		}

		f_dg<<endl;
		f_g<<g[i]<<endl;
		f_h<<h[i]<<","<<th_h[i]<<endl;
		f_J<<HYPER[i].J<<endl;
	}
	f_dt.close();
	f_lam.close();
	f_dg.close();
	f_g.close();
	f_p.close();
	f_q.close();
	f_d.close();
	f_h.close();
	f_J.close();
	f_dgdq.close();


	//stringstream ss0;
	//ss0<<"./rT/rT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss1;
	//ss1<<"./dT/dT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss3;
	//ss3<<"./lam/lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss5;
	//ss5<<"./mu/mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss6;
	//ss6<<"./g/g_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss7;
	//ss7<<"./dg/dg_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss8;
	//ss8<<"./p/hp_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss9;
	//ss9<<"./q/q_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss10;
	//ss10<<"./dT/d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss11;
	//ss11<<"./E/E_"<<t<<"_"<<count_min<<".csv";
	//stringstream ss12;
	//ss12<<"./h/h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss13;
	//ss13<<"./g/J_"<<t<<"_"<<count_min<<"_"<<count<<".csv";

	//ofstream f_rt(ss0.str());
	//ofstream f_dt(ss1.str());
	//ofstream f_lam(ss3.str());
	//ofstream f_mu(ss5.str());
	//ofstream f_g(ss6.str());
	//ofstream f_dg(ss7.str());
	//ofstream f_p(ss8.str());
	//ofstream f_q(ss9.str());
	//ofstream f_d(ss10.str());
	//ofstream f_h(ss12.str());
	//ofstream f_J(ss13.str());

	//if(count==1)
	//{
	//	ofstream f_E(ss11.str(), ios::trunc);
	//	f_E<<count<<","<<E<<endl;
	//	f_E.close();
	//}	
	//else
	//{
	//	ofstream f_E(ss11.str(), ios::app);
	//	f_E<<count<<","<<E<<endl;
	//	f_E.close();
	//}



	//for(int i=0;i<Nx;i++)
	//{				
	//	for(int j=0;j<Nx;j++)
	//	{
	//		f_rt<<rT[i*Nx+j]<<",";
	//	}
	//	f_rt<<endl;
	//	f_dt<<dT[i]<<endl;
	//	f_d<<dN[i]<<endl;
	//}
	//for(int i=0;i<h_num;i++)
	//{
	//	f_lam<<HYPER[i].lam<<endl;
	//	f_mu<<HYPER[i].mu<<endl;
	//	f_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
	//	f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
	//	for(int k=0;k<Nx;k++)
	//	{
	//		f_dg<<dg[i][k]<<",";
	//	}
	//	f_dg<<endl;
	//	f_g<<g[i]<<endl;
	//	f_h<<h[i]<<","<<th_h[i]<<endl;
	//	f_J<<HYPER[i].J<<endl;
	//}
	//f_rt.close();
	//f_dt.close();
	//f_lam.close();
	//f_mu.close();
	//f_dg.close();
	//f_g.close();
	//f_p.close();
	//f_q.close();
	//f_d.close();
	//f_h.close();
	//f_J.close();
}


void calc_n(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double nG[DIMENSION]={0,0,1};

	/////////////pn1_2åvéZ
	for(int i=0;i<h_num;i++)
	{
		double p_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double DgDq_n[DIMENSION]={HYPER1[jn*h_num+i].DgDq_n[A_X], HYPER1[jn*h_num+i].DgDq_n[A_Y], HYPER1[jn*h_num+i].DgDq_n[A_Z]};
			p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*DgDq_n[A_Z];
			p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*DgDq_n[A_Z];
			p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*DgDq_n[A_Z];
		}
		double DgDq_nii[DIMENSION]={HYPER1[i*h_num+i].DgDq_n[A_X], HYPER1[i*h_num+i].DgDq_n[A_Y], HYPER1[i*h_num+i].DgDq_n[A_Z]};
		p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].lam)*DgDq_nii[A_X]+HYPER[i].stress_n[A_X][A_Y]				 *DgDq_nii[A_Y]+HYPER[i].stress_n[A_X][A_Z]				  *DgDq_nii[A_Z];
		p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]				*DgDq_nii[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].lam)*DgDq_nii[A_Y]+HYPER[i].stress_n[A_Y][A_Z]				  *DgDq_nii[A_Z];
		p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]				*DgDq_nii[A_X]+HYPER[i].stress_n[A_Z][A_Y]				 *DgDq_nii[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].lam)*DgDq_nii[A_Z];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];

		/////////////qåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}


	double **d_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	d_Fi[D]=new double [DIMENSION];

	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

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

		b[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[0][1]*d_Fi[0][1]+d_Fi[0][2]*d_Fi[0][2];
		b[0][1]=d_Fi[0][0]*d_Fi[1][0]+d_Fi[0][1]*d_Fi[1][1]+d_Fi[0][2]*d_Fi[1][2];
		b[0][2]=d_Fi[0][0]*d_Fi[2][0]+d_Fi[0][1]*d_Fi[2][1]+d_Fi[0][2]*d_Fi[2][2];
		b[1][0]=d_Fi[1][0]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[0][2];
		b[1][1]=d_Fi[1][0]*d_Fi[1][0]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[1][2]*d_Fi[1][2];
		b[1][2]=d_Fi[1][0]*d_Fi[2][0]+d_Fi[1][1]*d_Fi[2][1]+d_Fi[1][2]*d_Fi[2][2];
		b[2][0]=d_Fi[2][0]*d_Fi[0][0]+d_Fi[2][1]*d_Fi[0][1]+d_Fi[2][2]*d_Fi[0][2];
		b[2][1]=d_Fi[2][0]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[1][2];
		b[2][2]=d_Fi[2][0]*d_Fi[2][0]+d_Fi[2][1]*d_Fi[2][1]+d_Fi[2][2]*d_Fi[2][2];

		bb[0][0]=b[0][0]*b[0][0]+b[0][1]*b[1][0]+b[0][2]*b[2][0];
		bb[0][1]=b[0][0]*b[0][1]+b[0][1]*b[1][1]+b[0][2]*b[2][1];
		bb[0][2]=b[0][0]*b[0][2]+b[0][1]*b[1][2]+b[0][2]*b[2][2];
		bb[1][0]=b[1][0]*b[0][0]+b[1][1]*b[1][0]+b[1][2]*b[2][0];
		bb[1][1]=b[1][0]*b[0][1]+b[1][1]*b[1][1]+b[1][2]*b[2][1];
		bb[1][2]=b[1][0]*b[0][2]+b[1][1]*b[1][2]+b[1][2]*b[2][2];
		bb[2][0]=b[2][0]*b[0][0]+b[2][1]*b[1][0]+b[2][2]*b[2][0];
		bb[2][1]=b[2][0]*b[0][1]+b[2][1]*b[1][1]+b[2][2]*b[2][1];
		bb[2][2]=b[2][0]*b[0][2]+b[2][1]*b[1][2]+b[2][2]*b[2][2];

		double trace_b=b[0][0]+b[1][1]+b[2][2];
		double trace_bb=bb[0][0]+bb[1][1]+bb[2][2];
	
		HYPER[j].stress[0][0]=2./J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[j].stress[0][1]=2./J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[j].stress[0][2]=2./J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[j].stress[1][0]=2./J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[j].stress[1][1]=2./J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[j].stress[1][2]=2./J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[j].stress[2][0]=2./J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[j].stress[2][1]=2./J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[j].stress[2][2]=2./J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));
	}

	for(int D=0;D>DIMENSION;D++)	delete[]	d_Fi[D];
	delete[]	d_Fi;



	/////////////F, J, t_inverseÇÃåvéZ
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	double J=0.;
	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		int Ni=HYPER[i].N;	

		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+jn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+jn].aiin[A_X],	HYPER1[i*h_num+jn].aiin[A_Y],	HYPER1[i*h_num+jn].aiin[A_Z]};

			fi[0][0]+=w*(PART[jn].r[A_X]-PART[i].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[jn].r[A_X]-PART[i].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[jn].r[A_X]-PART[i].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[jn].r[A_Y]-PART[i].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[jn].r[A_Y]-PART[i].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[jn].r[A_Y]-PART[i].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[jn].r[A_Z]-PART[i].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[jn].r[A_Z]-PART[i].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[jn].r[A_Z]-PART[i].r[A_Z])*a[A_Z];
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

		J=calc_det3(p_Fi);
		HYPER[i].J=J;	

		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];

		for(int j=0;j<Ni;j++)
		{			
			int k=HYPER[i].NEI[j];
			HYPER1[i*h_num+k].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			HYPER1[i*h_num+k].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			HYPER1[i*h_num+k].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
			//cout<<"i"<<i<<"j"<<k<<"	"<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);
	}	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;


//cout<<"----------OK"<<endl;
}


void p_QP(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION])
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=10000;

	double E=1.;
	double E_min=1.;
	double E_sum=0.;
	double ep=1e-5;
	double ep_min=1e-5;
	double d_ep=1e-20;

	double p_p[DIMENSION]={0,0,0};

	double r=1000.;

	double *f=new double [h_num];	
	double *df=new double [h_num*h_num];

	double *G=new double [h_num];	
	double *dG=new double [h_num*h_num];

	double *h=new double [h_num];

	double *dh=new double [h_num*h_num];
	double *th_h=new double [h_num];

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		f[i]=0.;
		G[i]=0.;
		h[i]=0.;
		th_h[i]=0.;
		for(int j=0;j<h_num;j++)
		{
			dh[i*h_num+j]=0.;
			df[i*h_num+j]=0.;
			dG[i*h_num+j]=0.;
		}
	}


	p_nab(HYPER,HYPER1,dG,dh,Dt,V,mi,nG);
	
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

			/////////////påvéZ
			double DgDq_ii[DIMENSION]={0,0,0};
			double DgDq_ij[DIMENSION]={0,0,0};
			double DgDq_ji[DIMENSION]={0,0,0};
			double p_p[DIMENSION]={0,0,0};
			for(int i=0;i<h_num;i++)
			{
				p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;
				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{
					int jn=HYPER[i].NEI[j];
					DgDq_ji[A_X]=HYPER1[jn*h_num+i].DgDq[A_X];	DgDq_ji[A_Y]=HYPER1[jn*h_num+i].DgDq[A_Y];	DgDq_ji[A_Z]=HYPER1[jn*h_num+i].DgDq[A_Z];
					p_p[A_X]+=(HYPER[jn].stress[A_X][A_X]-HYPER[jn].lam)*DgDq_ji[A_X]+ HYPER[jn].stress[A_X][A_Y]				*DgDq_ji[A_Y]+ HYPER[jn].stress[A_X][A_Z]				*DgDq_ji[A_Z];
					p_p[A_Y]+= HYPER[jn].stress[A_Y][A_X]				*DgDq_ji[A_X]+(HYPER[jn].stress[A_Y][A_Y]-HYPER[jn].lam)*DgDq_ji[A_Y]+ HYPER[jn].stress[A_Y][A_Z]				*DgDq_ji[A_Z];
					p_p[A_Z]+= HYPER[jn].stress[A_Z][A_X]				*DgDq_ji[A_X]+ HYPER[jn].stress[A_Z][A_Y]				*DgDq_ji[A_Y]+(HYPER[jn].stress[A_Z][A_Z]-HYPER[jn].lam)*DgDq_ji[A_Z];
				}
				
				DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
				p_p[A_X]+=(HYPER[i].stress[A_X][A_X]-HYPER[i].lam)*DgDq_ii[A_X]+ HYPER[i].stress[A_X][A_Y]				*DgDq_ii[A_Y]+ HYPER[i].stress[A_X][A_Z]			  *DgDq_ii[A_Z];
				p_p[A_Y]+= HYPER[i].stress[A_Y][A_X]			  *DgDq_ii[A_X]+(HYPER[i].stress[A_Y][A_Y]-HYPER[i].lam)*DgDq_ii[A_Y]+ HYPER[i].stress[A_Y][A_Z]			  *DgDq_ii[A_Z];
				p_p[A_Z]+= HYPER[i].stress[A_Z][A_X]			  *DgDq_ii[A_X]+ HYPER[i].stress[A_Z][A_Y]				*DgDq_ii[A_Y]+(HYPER[i].stress[A_Z][A_Z]-HYPER[i].lam)*DgDq_ii[A_Z];
				HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+0.5*Dt*p_p[A_X];
				HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+0.5*Dt*p_p[A_Y];
				HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+0.5*Dt*p_p[A_Z];
			}

			//stringstream ss0;
			//ss0<<"fi"<<count_min<<"_"<<count<<".csv";
			//ofstream f_fi(ss0.str());
			//double fi=0.;
			for(int i=0;i<h_num;i++)
			{
				G[i]=0.;
				//fi=0.;
				int Ni=HYPER[i].N;
				DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
		
				for(int j=0;j<Ni;j++)
				{
					int jn=HYPER[i].NEI[j];
					DgDq_ij[A_X]=HYPER1[i*h_num+jn].DgDq[A_X];	DgDq_ij[A_Y]=HYPER1[i*h_num+jn].DgDq[A_Y];	DgDq_ij[A_Z]=HYPER1[i*h_num+jn].DgDq[A_Z];
					//fi+=(double)(1./mi)*(DgDq_ij[A_X]*HYPER[jn].p[A_X]+DgDq_ij[A_Y]*HYPER[jn].p[A_Y]+DgDq_ij[A_Z]*HYPER[jn].p[A_Z]);			
					G[i]+=1./mi*(DgDq_ij[A_X]*HYPER[jn].p[A_X]+DgDq_ij[A_Y]*HYPER[jn].p[A_Y]+DgDq_ij[A_Z]*HYPER[jn].p[A_Z]);
					//f_fi<<",,"<<jn<<","<<1./mi*(DgDq_ij[A_X]*HYPER[jn].p[A_X]+DgDq_ij[A_Y]*HYPER[jn].p[A_Y]+DgDq_ij[A_Z]*HYPER[jn].p[A_Z])<<","<<f[i]<<endl;
					//cout<<"	f"<<i<<"_"<<jn<<"="<<1./mi*(DgDq_ij[A_X]*HYPER[jn].p[A_X]+DgDq_ij[A_Y]*HYPER[jn].p[A_Y]+DgDq_ij[A_Z]*HYPER[jn].p[A_Z])<<endl;
				}
				G[i]+=1./mi*(DgDq_ii[A_X]*HYPER[i].p[A_X]+DgDq_ii[A_Y]*HYPER[i].p[A_Y]+DgDq_ii[A_Z]*HYPER[i].p[A_Z]);
				//fi+=1./mi*(DgDq_ii[A_X]*HYPER[i].p[A_X]+DgDq_ii[A_Y]*HYPER[i].p[A_Y]+DgDq_ii[A_Z]*HYPER[i].p[A_Z]);
				//f_fi<<",,"<<i<<","<<1./mi*(DgDq_ii[A_X]*HYPER[i].p[A_X]+DgDq_ii[A_Y]*HYPER[i].p[A_Y]+DgDq_ii[A_Z]*HYPER[i].p[A_Z])<<","<<f[i]<<endl;
				//cout<<"	f"<<i<<"_"<<i<<"="<<1./mi*(DgDq_ij[A_X]*HYPER[i].p[A_X]+DgDq_ij[A_Y]*HYPER[i].p[A_Y]+DgDq_ij[A_Z]*HYPER[i].p[A_Z])<<endl;
				h[i]=-1./mi*(nG[A_X]*HYPER[i].p[A_X]+nG[A_Y]*HYPER[i].p[A_Y]+nG[A_Z]*HYPER[i].p[A_Z]);
				//f_fi<<i<<","<<f[i]<<endl;
				//cout<<"f"<<count<<"="<<f[i]<<",	"<<fi<<endl<<endl;
			}
			for(int i=0;i<h_num;i++)
			{
				f[i]=G[i]*G[i];
				if(h[i]+th_h[i]>0)	f[i]+=0.5*r*(h[i]+th_h[i])*(h[i]+th_h[i]);
				for(int j=0;j<h_num;j++)
				{
					df[i*h_num+j]=2*dG[i*h_num+j]*G[i];
					if(h[i]+th_h[i]>0)	df[i*h_num+j]+=r*dh[i*h_num+j]*(h[i]+th_h[i]);
				}
			}
			//f_fi.close();
			//for(int i=0;i<h_num;i++)
			//{	
			//	//cout<<"h"<<count<<"="<<h[i]<<endl;

			//	if(h[i]+th_h[i]>0)
			//	{
			//		f[i]+=0.5*r*(h[i]+th_h[i])*(h[i]+th_h[i]);
			//		cout<<"	h"<<count<<"="<<h[i]<<endl;
			//	}
			//	//cout<<"h="<<h[i]<<endl;
			//	for(int j=0;j<h_num;j++)
			//	{
			//		//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
			//		if(h[i]+th_h[i]>0)	df[i*h_num+j]+=r*dh[i*h_num+j]*(h[i]+th_h[i]);
			//	}
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
			////for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			//stringstream ssf;
			//ssf<<"f"<<count_min<<"_"<<count<<".csv";
			//ofstream ff(ssf.str());
			//stringstream ssdf;
			//ssdf<<"df"<<count_min<<"_"<<count<<".csv";
			//ofstream fdf(ssdf.str());
			//for(int i=0;i<h_num;i++)
			//{
			//	ff<<f[i]<<","<<HYPER[i].N_f<<endl;
			//	for(int j=0;j<h_num;j++)
			//	{
			//		fdf<<df[i*h_num+j]<<",";
			//	}		
			//	fdf<<endl;
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
			//fdf<<endl;
			//for(int i=0;i<h_num;i++)
			//{
			//	for(int j=0;j<h_num;j++)
			//	{
			//		fdf<<HYPER1[i*h_num+j].N_Df<<",";
			//	}		
			//	fdf<<endl;
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
			//ff.close();
			//fdf.close();


			//for(int i=0;i<h_num;i++)	G[i]=f[i];
			////cout<<"T="<<T<<", En"<<En<<endl;
			//for(int i=0;i<h_num;i++)
			//{
			//	for(int j=0;j<h_num;j++)
			//	{
			//		dG[i*h_num+j]=df[i*h_num+j];
			//	}				
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}
	
			//for(int i=0;i<h_num;i++)
			//{
			//	G[i]=g[i];
			//	for(int j=0;j<h_num;j++)
			//	{
			//		dG[i*h_num+j]=HYPER1[i*h_num+j].N_Df;
			//	}		
			//	//cout<<"df"<<count<<"="<<df[i*h_num+h_num-1]<<endl;
			//}

			gauss(df,f,h_num);
			//gauss(dG,G,h_num);

			//for(int i=0;i<h_num;i++)
			//{
			//	//cout<<"f"<<count<<"="<<f[i]<<", G="<<G[i]<<endl;
			//	cout<<"f"<<count<<"="<<f[i]<<endl;
			//}
			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=f[i];
				E_sum+=fabs(f[i]);
			}
			E=E_sum;
			///////////////pn1_2åvéZ
			//for(int i=0;i<h_num;i++)
			//{
			//	p_p[A_X]=0;	p_p[A_Y]=0;	p_p[A_Z]=0;	
			//	int Ni=HYPER[i].N;
			//	for(int j=0;j<Ni;j++)
			//	{
			//		int jn=HYPER[i].NEI[j];
			//		double DgDq_n[DIMENSION]={HYPER1[jn*h_num+i].DgDq_n[A_X], HYPER1[jn*h_num+i].DgDq_n[A_Y], HYPER1[jn*h_num+i].DgDq_n[A_Z]};
			//		p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*DgDq_n[A_Z];
			//		p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*DgDq_n[A_Z];
			//		p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*DgDq_n[A_Z];
			//	}
			//	double DgDq_nii[DIMENSION]={HYPER1[i*h_num+i].DgDq_n[A_X], HYPER1[i*h_num+i].DgDq_n[A_Y], HYPER1[i*h_num+i].DgDq_n[A_Z]};
			//	p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].lam)*DgDq_nii[A_X]+HYPER[i].stress_n[A_X][A_Y]*DgDq_nii[A_Y]+HYPER[i].stress_n[A_X][A_Z]*DgDq_nii[A_Z];
			//	p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*DgDq_nii[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].lam)*DgDq_nii[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*DgDq_nii[A_Z];
			//	p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*DgDq_nii[A_X]+HYPER[i].stress_n[A_Z][A_Y]*DgDq_nii[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].lam)*DgDq_nii[A_Z];
			//	HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
			//	HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
			//	HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];

			//	/////////////qåvéZ
			//	PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
			//	PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
			//	PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
			//}

			cout<<"E"<<count<<"="<<E<<endl;
			output_data_p(HYPER,HYPER1,dG,G,f,dh,h,th_h,Nx,h_num,count,count_min,t,E);
		}

		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=(HYPER[i].old_lam-HYPER[i].lam)*(HYPER[i].old_lam-HYPER[i].lam);
		E_min=sqrt(E_sum);
		if(E_min<ep_min*1000)	r*=4;
		cout<<"Emin"<<count_min<<"="<<E_min<<endl;

		for(int i=0;i<h_num;i++)
		{
			if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		}
	}




	///////////////pn1_2åvéZ
	//for(int i=0;i<h_num;i++)
	//{
	//	p_p[A_X]=0;	p_p[A_Y]=0;	p_p[A_Z]=0;	
	//	int Ni=HYPER[i].N;
	//	for(int j=0;j<Ni;j++)
	//	{
	//		int jn=HYPER[i].NEI[j];
	//		double DgDq_n2[DIMENSION]={HYPER1[jn*h_num+i].DgDq_n[A_X], HYPER1[jn*h_num+i].DgDq_n[A_Y], HYPER1[jn*h_num+i].DgDq_n[A_Z]};
	//		p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*DgDq_n2[A_X]+HYPER[jn].stress_n[A_X][A_Y]*DgDq_n2[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*DgDq_n2[A_Z];
	//		p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*DgDq_n2[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*DgDq_n2[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*DgDq_n2[A_Z];
	//		p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*DgDq_n2[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*DgDq_n2[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*DgDq_n2[A_Z];
	//	}
	//	double DgDq_nii2[DIMENSION]={HYPER1[i*h_num+i].DgDq_n[A_X], HYPER1[i*h_num+i].DgDq_n[A_Y], HYPER1[i*h_num+i].DgDq_n[A_Z]};
	//	p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].lam)*DgDq_nii2[A_X]+HYPER[i].stress_n[A_X][A_Y]*DgDq_nii2[A_Y]+HYPER[i].stress_n[A_X][A_Z]*DgDq_nii2[A_Z];
	//	p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*DgDq_nii2[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].lam)*DgDq_nii2[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*DgDq_nii2[A_Z];
	//	p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*DgDq_nii2[A_X]+HYPER[i].stress_n[A_Z][A_Y]*DgDq_nii2[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].lam)*DgDq_nii2[A_Z];
	//	HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
	//	HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
	//	HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];

	//	/////////////qåvéZ
	//	PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
	//	PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
	//	PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	//}

	delete[]	f;
	delete[]	df;
	delete[]	h;
	delete[]	dh;
	delete[]	th_h;
	delete[]	dG;
	delete[]	G;

}

void p_nab(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dG, double *dh, double Dt, double V, double mi,double nG[DIMENSION])
{

	int h_num=HYPER.size();
	double dm=0.;
	/////////////////dEåvéZ///////////////////
	for(int k=0;k<h_num;k++)
	{
		for(int i=0;i<h_num;i++)	dG[i*h_num+k]=0.;
	}

	//for(int k=0;k<h_num;k++)
	//{
	//	int Nk=HYPER[k].N;
	//	for(int i=0;i<Nk;i++)
	//	{
	//		int in=HYPER[k].NEI[i];

	//		for(int j=0;j<Nk;j++)
	//		{
	//			int jn=HYPER[k].NEI[j];

	//			dg[in*h_num+jn]-=0.5*Dt*Dt/mi*(HYPER1[in*h_num+k].DgDq[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+HYPER1[in*h_num+k].DgDq[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+HYPER1[in*h_num+k].DgDq[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
	//		}
	//		dg[in*h_num+k]-=0.5*Dt*Dt/mi*(HYPER1[in*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]+HYPER1[in*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]+HYPER1[in*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);
	//		dg[k*h_num+in]-=0.5*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[in*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[in*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[in*h_num+k].DgDq_n[A_Z]);

	//	}
	//	dg[k*h_num+k]-=0.5*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);
	//}
	double DgDq_kk[DIMENSION]={0,0,0};
	double DgDq_ik[DIMENSION]={0,0,0};
	double DgDq_ki[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		DgDq_kk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	DgDq_kk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	DgDq_kk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];

		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];

			DgDq_ik[A_X]=HYPER1[in*h_num+k].DgDq[A_X];	DgDq_ik[A_Y]=HYPER1[in*h_num+k].DgDq[A_Y];	DgDq_ik[A_Z]=HYPER1[in*h_num+k].DgDq[A_Z];
			DgDq_ki[A_X]=HYPER1[k*h_num+in].DgDq[A_X];	DgDq_ki[A_Y]=HYPER1[k*h_num+in].DgDq[A_Y];	DgDq_ki[A_Z]=HYPER1[k*h_num+in].DgDq[A_Z];

			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				dG[in*h_num+jn]-=0.5*Dt/mi*(DgDq_ik[A_X]*HYPER1[jn*h_num+k].DgDq[A_X]+DgDq_ik[A_Y]*HYPER1[jn*h_num+k].DgDq[A_Y]+DgDq_ik[A_Z]*HYPER1[jn*h_num+k].DgDq[A_Z]);
			}
			dG[in*h_num+k]-=0.5*Dt/mi*(DgDq_ik[A_X]*DgDq_kk[A_X]+DgDq_ik[A_Y]*DgDq_kk[A_Y]+DgDq_ik[A_Z]*DgDq_kk[A_Z]);
			dG[k*h_num+in]-=0.5*Dt/mi*(DgDq_kk[A_X]*DgDq_ik[A_X]+DgDq_kk[A_Y]*DgDq_ik[A_Y]+DgDq_kk[A_Z]*DgDq_ik[A_Z]);

			dh[in*h_num+k]=0.5*Dt/mi*(nG[A_X]*DgDq_ki[A_X]+nG[A_Y]*DgDq_ki[A_Y]+nG[A_Z]*DgDq_ki[A_Z]);
		}
		dh[k*h_num+k]=0.5*Dt/mi*(nG[A_X]*DgDq_kk[A_X]+nG[A_Y]*DgDq_kk[A_Y]+nG[A_Z]*DgDq_kk[A_Z]);
		dG[k*h_num+k]-=0.5*Dt/mi*(DgDq_kk[A_X]*DgDq_kk[A_X]+DgDq_kk[A_Y]*DgDq_kk[A_Y]+DgDq_kk[A_Z]*DgDq_kk[A_Z]);
	}

}


void p_variables(vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int h_num,double Dt,double mi)
{
	/////////////påvéZ
	for(int i=0;i<h_num;i++)
	{
		double p_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double DgDq_ji[DIMENSION]={HYPER1[jn*h_num+i].DgDq[A_X], HYPER1[jn*h_num+i].DgDq[A_Y], HYPER1[jn*h_num+i].DgDq[A_Z]};
			p_p[A_X]+=(HYPER[jn].stress[A_X][A_X]-HYPER[jn].lam)*DgDq_ji[A_X]+HYPER[jn].stress[A_X][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_X][A_Z]*DgDq_ji[A_Z];
			p_p[A_Y]+=HYPER[jn].stress[A_Y][A_X]*DgDq_ji[A_X]+(HYPER[jn].stress[A_Y][A_Y]-HYPER[jn].lam)*DgDq_ji[A_Y]+HYPER[jn].stress[A_Y][A_Z]*DgDq_ji[A_Z];
			p_p[A_Z]+=HYPER[jn].stress[A_Z][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_Z][A_Y]*DgDq_ji[A_Y]+(HYPER[jn].stress[A_Z][A_Z]-HYPER[jn].lam)*DgDq_ji[A_Z];
		}
		double DgDqii[DIMENSION]={HYPER1[i*h_num+i].DgDq[A_X], HYPER1[i*h_num+i].DgDq[A_Y], HYPER1[i*h_num+i].DgDq[A_Z]};
		p_p[A_X]+=(HYPER[i].stress[A_X][A_X]-HYPER[i].lam)*DgDqii[A_X]+HYPER[i].stress[A_X][A_Y]				 *DgDqii[A_Y]+HYPER[i].stress[A_X][A_Z]				  *DgDqii[A_Z];
		p_p[A_Y]+=HYPER[i].stress[A_Y][A_X]				  *DgDqii[A_X]+(HYPER[i].stress[A_Y][A_Y]-HYPER[i].lam)  *DgDqii[A_Y]+HYPER[i].stress[A_Y][A_Z]				  *DgDqii[A_Z];
		p_p[A_Z]+=HYPER[i].stress[A_Z][A_X]			      *DgDqii[A_X]+HYPER[i].stress[A_Z][A_Y]				 *DgDqii[A_Y]+(HYPER[i].stress[A_Z][A_Z]-HYPER[i].lam)*DgDqii[A_Z];
		HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+0.5*Dt*p_p[A_X];
		HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+0.5*Dt*p_p[A_Y];
		HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+0.5*Dt*p_p[A_Z];
	}
}



void output_data_p(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double *dG, double *G,double *f, double *dh, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E)
{
	if(count_min==1&&count==1)
	{
		stringstream ss1;
		ss1<<"./dT/p_dG_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
		stringstream ss13;
		ss13<<"./dh/p_dh_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
		stringstream ss14;
		ss14<<"stress_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
		stringstream ss15;
		ss15<<"DgDq_"<<t<<"_"<<count_min<<"_"<<count<<".csv";


		ofstream f_dt(ss1.str());
		ofstream f_dh(ss13.str());
		ofstream f_s(ss14.str());
		ofstream f_dg(ss15.str());
		for(int i=0;i<h_num;i++)
		{				
			for(int j=0;j<h_num;j++)
			{
				f_dt<<dG[i*h_num+j]<<",";
				f_dh<<dh[i*h_num+j]<<",";
			}
			f_dt<<endl;
			f_dh<<endl;
			f_dg<<i<<endl;
			f_s<<i<<","<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<endl;
			f_s<<","<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<endl;
			f_s<<","<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<endl;
			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{
				int jn=HYPER[i].NEI[j];
				f_dg<<","<<jn<<","<<HYPER1[i*h_num+jn].DgDq[A_X]<<","<<HYPER1[i*h_num+jn].DgDq[A_Y]<<","<<HYPER1[i*h_num+jn].DgDq[A_Z]<<endl;
			}
			
		}
		f_s.close();
		f_dg.close();
		f_dt.close();
		f_dh.close();
	}

	




	stringstream ss2;
	ss2<<"./dT/p_G_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./lam/p_lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./p/p_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./dT/p_d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./E/p_E_"<<t<<"_"<<count_min<<".csv";
	stringstream ss12;
	ss12<<"./h/p_h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	ofstream f_fp(ss2.str());
	ofstream f_lam(ss3.str());
	ofstream f_p(ss8.str());
	ofstream f_d(ss10.str());
	ofstream f_h(ss12.str());

	if(count==1)
	{
		ofstream f_E(ss11.str(), ios::trunc);
		f_E<<count<<","<<E<<endl;
		f_E.close();
	}	
	else
	{
		ofstream f_E(ss11.str(), ios::app);
		f_E<<count<<","<<E<<endl;
		f_E.close();
	}

	for(int i=0;i<h_num;i++)
	{				
		f_fp<<G[i]<<endl;
		f_d<<f[i]<<endl;
	}
	for(int i=0;i<h_num;i++)
	{
		f_lam<<HYPER[i].lam<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;

		f_h<<h[i]<<","<<th_h[i]<<endl;
	}
	f_fp.close();
	f_lam.close();
	f_p.close();
	f_d.close();
	f_h.close();

}










