#include "stdafx.h"	

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, double V, double mi, double Dt,vector<double > NEIw);
void q_nab_lap(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double **dg, double Dt, double V, double mi, vector<double > NEIw);
void q_variables(mpsconfig &CON, vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, vector<double > NEIw);
void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double **dg, double *h,double **dh, double *d, int h_num,int count,int t,double E,double mi, double V,vector<double > NEIw);

void calc_n(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void p_QP(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt,vector<double > NEIw);
void p_variables(vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int h_num,double Dt,double mi);
void p_nab(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, double *G,double *H, double Dt, double V, double mi);
void p_lap(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double **dH, double **dG, double Dt, double V, double mi);
void output_data_p(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double **dG, double *G,double **dH,double *H, double *rT, double *dT, double T, double *d, int Nx, int h_num,int count,int t,double E);
//void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT, double *dN,double *g, double **dg, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E);


void calc_HYPER_QP_gh(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F,vector<double > NEIw)
{
	////////////íËã`///////////////
	int h_num=HYPER.size();
	int Nx=h_num*2;

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
	double nG[DIMENSION]={0,0,1.};

	double E0=0;	
	if(t==1)
	{
		HYPER[0].E=0.;
		ofstream fe0("E0.csv");
		for(int i=0;i<h_num;i++)
		{
			HYPER[0].E+=0.5/mi*(HYPER[i].p_n[A_X]*HYPER[i].p_n[A_X]+HYPER[i].p_n[A_Y]*HYPER[i].p_n[A_Y]+HYPER[i].p_n[A_Z]*HYPER[i].p_n[A_Z])+V*HYPER[i].W_n;
			fe0<<0.5/mi*(HYPER[i].p_n[A_X]*HYPER[i].p_n[A_X]+HYPER[i].p_n[A_Y]*HYPER[i].p_n[A_Y]+HYPER[i].p_n[A_Z]*HYPER[i].p_n[A_Z])+V*HYPER[i].W_n<<","<<0.5/mi*(HYPER[i].p_n[A_X]*HYPER[i].p_n[A_X]+HYPER[i].p_n[A_Y]*HYPER[i].p_n[A_Y]+HYPER[i].p_n[A_Z]*HYPER[i].p_n[A_Z])<<","<<V*HYPER[i].W_n<<endl;
		}
		fe0.close();
	}
	E0=HYPER[0].E;
	cout<<endl;
	cout<<"E0="<<E0<<endl;


	cout<<"QP start-------------------------"<<endl;
	ofstream fq0("q0_QP.csv");
	ofstream fp0("hp0_QP.csv");
	ofstream fhp0("hp0_QP.csv");

	ofstream fqn("qn_QP.csv");
	ofstream fpn("pn_QP.csv");
	ofstream fhpn("hpn_QP.csv");
	for(int i=0;i<h_num;i++)
	{
		fq0<<HYPER[i].q_n[A_X]<<","<<HYPER[i].q_n[A_Y]<<","<<HYPER[i].q_n[A_Z]<<endl;
		fp0<<HYPER[i].p_n[A_X]<<","<<HYPER[i].p_n[A_Y]<<","<<HYPER[i].p_n[A_Z]<<endl;
		fhp0<<HYPER[i].ph_n[A_X]<<","<<HYPER[i].ph_n[A_Y]<<","<<HYPER[i].ph_n[A_Z]<<endl;
		fqn<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		fpn<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		fhpn<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
	}
	fq0.close();
	fp0.close();
	fhp0.close();
	fqn.close();
	fpn.close();
	fhpn.close();

	////////////QPåvéZ///////////////		
	q_QP(CON,PART,HYPER,HYPER1,t,h_num,V,mi,Dt,NEIw);
	for(int i=0;i<h_num;i++)	HYPER[i].lambda=HYPER[i].lam;

	calc_n(CON,PART,HYPER,HYPER1);

	p_QP(HYPER,HYPER1,t,h_num,Nx,V,mi,Dt,NEIw);

	cout<<"--------------------------OK"<<endl;


}

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, double V, double mi, double Dt,vector<double > NEIw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;
	double E_old=1.;
	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1.e-5;
	double ep_min=1.e-5;
	double d_ep=1.e-20;

	double hp[DIMENSION]={0,0,0};
	double W=0.;
	double Dgki_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};
	double p_p[DIMENSION]={0,0,0};
	double Dgji_n[DIMENSION]={0,0,0};
	int fw=0;
	double lam=0.;
	double mu=0.;

	double r=10;

	int Nw=NEIw.size();
	int Nx=Nw+h_num;
	double T=0.;
	double *dT=new double [Nx];	
	double *rT=new double [Nx*Nx];

	double *d=new double [Nx];
	double *B=new double [Nx*Nx];

	double *g=new double [h_num];
	double **dg=new double *[Nx];
	double *h=new double [Nw];
	double **dh=new double *[Nw];

	for(int i=0;i<h_num;i++)
	{
		dg[i]=new double [Nx];
		dh[i]=new double [Nx];
	}
	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		g[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			dg[i][j]=0.;
		}
	}
	for(int i=0;i<Nw;i++)
	{
		h[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			dh[i][j]=0.;
		}
	}


	for(int i=0;i<Nx;i++)
	{
		dT[i]=0.;
		d[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			rT[i*Nx+j]=0.;
			B[i*Nx+j]=0.;
		}
	}

	for(int i=0;i<Nw;i++)
	{
		int in=NEIw[i];
		int Ni=HYPER[in].N;
		Dgii_n[A_Z]=HYPER1[in*h_num+in].DgDq_n[A_Z];
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[in].NEI[k];
			Dgki_n[A_Z]=HYPER1[kn*h_num+in].DgDq_n[A_Z];
			dh[i][kn]+=0.5*Dt*Dt/mi*Dgki_n[A_Z];
		}
		dh[i][in]+=0.5*Dt*Dt/mi*Dgii_n[A_Z];
	}
	for(int k=0;k<Nw;k++)
	{
		dh[k][k+h_num]=-Dt*Dt*0.5/mi;
	}


	E=1;
	count=0;
	while(E>ep)
	{
		count++;
		q_variables(CON,PART,HYPER,HYPER1,NEIw);
		q_nab_lap(PART,HYPER,HYPER1,dg,Dt,V,mi,NEIw);

		T=0;
		for(int i=0;i<h_num;i++)
		{	
			g[i]=V*(1-HYPER[i].J);
			T+=g[i]*g[i];
		}
		for(int i=0;i<Nw;i++)
		{
			int in=NEIw[i];
			h[i]=-1.*PART[in].r[A_Z];
			if(h[i]>0)	T+=h[i]*h[i];
		}

		for(int k=0;k<Nx;k++)
		{	
			dT[k]=0.;
			for(int i=0;i<h_num;i++)
			{	
				dT[k]+=2.*dg[i][k]*g[i];
			}
			for(int i=0;i<Nw;i++)
			{
				if(h[i]>0)	dT[k]+=2.*dh[i][k]*h[i];
			}

			for(int l=0;l<Nx;l++)
			{
				rT[l*Nx+k]=0.;
				for(int i=0;i<h_num;i++)
				{	
					rT[l*Nx+k]+=2.*dg[i][k]*dg[i][l];
				}
				for(int i=0;i<Nw;i++)
				{
					if(h[i]>0)	rT[l*Nx+k]+=2.*dh[i][k]*dh[i][l];
					else if(h[i]>-1e-20 && dh[i][k]*dh[i][l]>0)	rT[l*Nx+k]+=2.*dh[i][k]*dh[i][l];
				}
			}
		}



		for(int i=0;i<Nx;i++)
		{
			d[i]=dT[i];
			for(int j=0;j<Nx;j++)
			{
				B[i*Nx+j]=rT[i*Nx+j];
			}
		}
		gauss(B,d,Nx);

		E_sum=0;
		for(int i=0;i<h_num;i++)
		{
			HYPER[i].lam-=d[i];
			E_sum+=fabs(d[i]);
		}
		for(int i=0;i<Nw;i++)
		{
			int in=NEIw[i];
			HYPER[in].mu-=d[i+h_num];
			E_sum+=fabs(d[i+h_num]);
		}
		E_old=E;
		E=E_sum;
		//if(E>E_old&&count!=1)	break;
		/////////////pn1_2åvéZ
		for(int i=0;i<h_num;i++)
		{
			p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{
				int jn=HYPER[i].NEI[j];
				Dgji_n[A_X]=HYPER1[jn*h_num+i].DgDq_n[A_X];	Dgji_n[A_Y]=HYPER1[jn*h_num+i].DgDq_n[A_Y];	Dgji_n[A_Z]=HYPER1[jn*h_num+i].DgDq_n[A_Z];
				lam=HYPER[jn].lam;


				p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-lam)*Dgji_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*Dgji_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*Dgji_n[A_Z];
				p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*Dgji_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-lam)*Dgji_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*Dgji_n[A_Z];
				p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*Dgji_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*Dgji_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-lam)*Dgji_n[A_Z];
			}
			Dgii_n[A_X]=HYPER1[i*h_num+i].DgDq_n[A_X];	Dgii_n[A_Y]=HYPER1[i*h_num+i].DgDq_n[A_Y];	Dgii_n[A_Z]=HYPER1[i*h_num+i].DgDq_n[A_Z];
			lam=HYPER[i].lam;
			p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-lam)*Dgii_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*Dgii_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*Dgii_n[A_Z];
			p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*Dgii_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-lam)*Dgii_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*Dgii_n[A_Z];
			p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*Dgii_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*Dgii_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-lam)*Dgii_n[A_Z];
			HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
			HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
			HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];
			fw=HYPER[i].fw;
			if(fw==ON)
			{
				mu=HYPER[i].mu;
				HYPER[i].half_p[A_Z]+=Dt*0.5*mu;
			}
			//HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X]+0.5*Dt*nG[A_X]*HYPER[i].mu;
			//HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y]+0.5*Dt*nG[A_Y]*HYPER[i].mu;
			//HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*nG[A_Z]*HYPER[i].mu;

			/////////////qåvéZ
			PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
			PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
			PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
		}	

		cout<<"	E,T"<<count<<"="<<E<<","<<T<<endl;
		output_data(PART,HYPER,HYPER1,T,dT,rT,g,dg,h,dh,d,h_num,count,t,E,mi,V,NEIw);
	}

	for(int i=0;i<h_num;i++)
	{
		delete[]	dg[i];
		delete[]	dh[i];
	}

	delete[]	g;
	delete[]	dg;
	delete[]	h;
	delete[]	dh;
	delete[]	dT;
	delete[]	rT;
	delete[]	d;
	delete[]	B;

}

void q_variables(mpsconfig &CON, vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, vector<double > NEIw)
{
	int h_num=HYPER.size();
	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();

	/////////////pn1_2åvéZ
	double p_p[DIMENSION]={0,0,0};
	double Dgji_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};
	int fw=0;
	double lam=0.;
	double mu=0.;
	for(int i=0;i<h_num;i++)
	{
		p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;

		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			Dgji_n[A_X]=HYPER1[jn*h_num+i].DgDq_n[A_X];	Dgji_n[A_Y]=HYPER1[jn*h_num+i].DgDq_n[A_Y];	Dgji_n[A_Z]=HYPER1[jn*h_num+i].DgDq_n[A_Z];
			lam=HYPER[jn].lam;


			p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-lam)*Dgji_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*Dgji_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*Dgji_n[A_Z];
			p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*Dgji_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-lam)*Dgji_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*Dgji_n[A_Z];
			p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*Dgji_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*Dgji_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-lam)*Dgji_n[A_Z];
		}
		Dgii_n[A_X]=HYPER1[i*h_num+i].DgDq_n[A_X];	Dgii_n[A_Y]=HYPER1[i*h_num+i].DgDq_n[A_Y];	Dgii_n[A_Z]=HYPER1[i*h_num+i].DgDq_n[A_Z];
		lam=HYPER[i].lam;
		p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-lam)*Dgii_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*Dgii_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*Dgii_n[A_Z];
		p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*Dgii_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-lam)*Dgii_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*Dgii_n[A_Z];
		p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*Dgii_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*Dgii_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-lam)*Dgii_n[A_Z];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z];
		fw=HYPER[i].fw;
		if(fw==ON)
		{
			mu=HYPER[i].mu;
			HYPER[i].half_p[A_Z]+=Dt*0.5*mu;
		}
		//HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X]+0.5*Dt*nG[A_X]*HYPER[i].mu;
		//HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y]+0.5*Dt*nG[A_Y]*HYPER[i].mu;
		//HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*nG[A_Z]*HYPER[i].mu;

		/////////////qåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}	


	/////////////F, J, t_inverseÇÃåvéZ
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	double J=0.;
	double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double a[DIMENSION]={0,0,0};
	for(int i=0;i<h_num;i++)
	{
		fi[A_X][0]=0.;	fi[A_X][1]=0.;	fi[A_X][2]=0.;
		fi[A_Y][0]=0.;	fi[A_Y][1]=0.;	fi[A_Y][2]=0.;
		fi[A_Z][0]=0.;	fi[A_Z][1]=0.;	fi[A_Z][2]=0.;

		int Ni=HYPER[i].N;	

		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+jn].wiin;
			a[A_X]=HYPER1[i*h_num+jn].aiin[A_X];	a[A_Y]=HYPER1[i*h_num+jn].aiin[A_Y];	a[A_Z]=HYPER1[i*h_num+jn].aiin[A_Z];

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


	/////////////DgDq, Stress, W, S, dSdcåvéZ
	double c10=CON.get_c10();
	double c01=CON.get_c01();
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

	double trace_b=0., trace_bb=0.;
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double S[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dSdc[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double F_dSdC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double F_dSdC_F[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double n0kj[DIMENSION]={0,0,0};
	double n0jj[DIMENSION]={0,0,0};

	for(int j=0;j<h_num;j++)
	{	

		////////WåvéZ
		double J=HYPER[j].J;	
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1./3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1./3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1./3.)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1./3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1./3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1./3.)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1./3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1./3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1./3.)*HYPER[j].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1./3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1./3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1./3.)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1./3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1./3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1./3.)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1./3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1./3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1./3.)*HYPER[j].Fi[2][2];
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

		HYPER[j].W=c10*(Ic-3.)+c01*(IIc-3.);

		////////StressåvéZ	
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

		trace_b=b[0][0]+b[1][1]+b[2][2];
		trace_bb=bb[0][0]+bb[1][1]+bb[2][2];

		HYPER[j].stress[0][0]=2./J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[j].stress[0][1]=2./J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[j].stress[0][2]=2./J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[j].stress[1][0]=2./J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[j].stress[1][1]=2./J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[j].stress[1][2]=2./J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[j].stress[2][0]=2./J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[j].stress[2][1]=2./J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[j].stress[2][2]=2./J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));


		////////S, dSdcåvéZ
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
			S[0][0]=-2.*1/pow(-J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[0][0]) + c01*(Ic-dC[0][0]-2./3.*IIc*in_Ci[0][0]) );
			S[0][1]=-2.*1/pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][1]) + c01*(  -dC[0][1]-2./3.*IIc*in_Ci[0][1]) );
			S[0][2]=-2.*1/pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][2]) + c01*(  -dC[0][2]-2./3.*IIc*in_Ci[0][2]) );
			S[1][0]=-2.*1/pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][0]) + c01*(  -dC[1][0]-2./3.*IIc*in_Ci[1][0]) );
			S[1][1]=-2.*1/pow(-J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[1][1]) + c01*(Ic-dC[1][1]-2./3.*IIc*in_Ci[1][1]) );
			S[1][2]=-2.*1/pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][2]) + c01*(  -dC[1][2]-2./3.*IIc*in_Ci[1][2]) );
			S[2][0]=-2.*1/pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][0]) + c01*(  -dC[2][0]-2./3.*IIc*in_Ci[2][0]) );
			S[2][1]=-2.*1/pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][1]) + c01*(  -dC[2][1]-2./3.*IIc*in_Ci[2][1]) );	
			S[2][2]=-2.*1/pow(-J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[2][2]) + c01*(Ic-dC[2][2]-2./3.*IIc*in_Ci[2][2]) );

			dSdc[0][0]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][0]-in_Ci[0][0])+c01*(-2.*Ic*in_Ci[0][0]+5./3.*(1.+IIc*in_Ci2[0][0])) );
			dSdc[0][1]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][1]-in_Ci[0][1])+c01*(-2.*Ic*in_Ci[0][1]+5./3.*( +IIc*in_Ci2[0][1])) );
			dSdc[0][2]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][2]-in_Ci[0][2])+c01*(-2.*Ic*in_Ci[0][2]+5./3.*( +IIc*in_Ci2[0][2])) );
			dSdc[1][0]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][0]-in_Ci[1][0])+c01*(-2.*Ic*in_Ci[1][0]+5./3.*( +IIc*in_Ci2[1][0])) );
			dSdc[1][1]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][1]-in_Ci[1][1])+c01*(-2.*Ic*in_Ci[1][1]+5./3.*(1.+IIc*in_Ci2[1][1])) );
			dSdc[1][2]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][2]-in_Ci[1][2])+c01*(-2.*Ic*in_Ci[1][2]+5./3.*( +IIc*in_Ci2[1][2])) );
			dSdc[2][0]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][0]-in_Ci[2][0])+c01*(-2.*Ic*in_Ci[2][0]+5./3.*( +IIc*in_Ci2[2][0])) );
			dSdc[2][1]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][1]-in_Ci[2][1])+c01*(-2.*Ic*in_Ci[2][1]+5./3.*( +IIc*in_Ci2[2][1])) );
			dSdc[2][2]=-4./3.*1/pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][2]-in_Ci[2][2])+c01*(-2.*Ic*in_Ci[2][2]+5./3.*(1.+IIc*in_Ci2[2][2])) );
		}
		else
		{
			S[0][0]=2.*1/pow(J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[0][0]) + c01*(Ic-dC[0][0]-2./3.*IIc*in_Ci[0][0]) );
			S[0][1]=2.*1/pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][1]) + c01*(  -dC[0][1]-2./3.*IIc*in_Ci[0][1]) );
			S[0][2]=2.*1/pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][2]) + c01*(  -dC[0][2]-2./3.*IIc*in_Ci[0][2]) );
			S[1][0]=2.*1/pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][0]) + c01*(  -dC[1][0]-2./3.*IIc*in_Ci[1][0]) );
			S[1][1]=2.*1/pow(J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[1][1]) + c01*(Ic-dC[1][1]-2./3.*IIc*in_Ci[1][1]) );
			S[1][2]=2.*1/pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][2]) + c01*(  -dC[1][2]-2./3.*IIc*in_Ci[1][2]) );
			S[2][0]=2.*1/pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][0]) + c01*(  -dC[2][0]-2./3.*IIc*in_Ci[2][0]) );
			S[2][1]=2.*1/pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][1]) + c01*(  -dC[2][1]-2./3.*IIc*in_Ci[2][1]) );	
			S[2][2]=2.*1/pow(J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[2][2]) + c01*(Ic-dC[2][2]-2./3.*IIc*in_Ci[2][2]) );


			dSdc[0][0]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][0]-in_Ci[0][0])+c01*(-2.*Ic*in_Ci[0][0]+5./3.*(1.+IIc*in_Ci2[0][0])) );
			dSdc[0][1]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][1]-in_Ci[0][1])+c01*(-2.*Ic*in_Ci[0][1]+5./3.*( +IIc*in_Ci2[0][1])) );
			dSdc[0][2]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][2]-in_Ci[0][2])+c01*(-2.*Ic*in_Ci[0][2]+5./3.*( +IIc*in_Ci2[0][2])) );
			dSdc[1][0]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][0]-in_Ci[1][0])+c01*(-2.*Ic*in_Ci[1][0]+5./3.*( +IIc*in_Ci2[1][0])) );
			dSdc[1][1]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][1]-in_Ci[1][1])+c01*(-2.*Ic*in_Ci[1][1]+5./3.*(1.+IIc*in_Ci2[1][1])) );
			dSdc[1][2]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][2]-in_Ci[1][2])+c01*(-2.*Ic*in_Ci[1][2]+5./3.*( +IIc*in_Ci2[1][2])) );
			dSdc[2][0]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][0]-in_Ci[2][0])+c01*(-2.*Ic*in_Ci[2][0]+5./3.*( +IIc*in_Ci2[2][0])) );
			dSdc[2][1]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][1]-in_Ci[2][1])+c01*(-2.*Ic*in_Ci[2][1]+5./3.*( +IIc*in_Ci2[2][1])) );
			dSdc[2][2]=4./3.*1/pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][2]-in_Ci[2][2])+c01*(-2.*Ic*in_Ci[2][2]+5./3.*(1.+IIc*in_Ci2[2][2])) );
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

		int Nj=HYPER[j].N;
		for(int k=0;k<Nj;k++)
		{
			int kn=HYPER[j].NEI[k];
			n0kj[A_X]=HYPER1[kn*h_num+j].n0ij[A_X];	n0kj[A_Y]=HYPER1[kn*h_num+j].n0ij[A_Y];	n0kj[A_Z]=HYPER1[kn*h_num+j].n0ij[A_Z];

			HYPER1[j*h_num+kn].DpiDq[A_X]=(S[A_X][0]+F_dSdC_F[A_X][0])*n0kj[0]+(S[A_X][1]+F_dSdC_F[A_X][1])*n0kj[1]+(S[A_X][2]+F_dSdC_F[A_X][2])*n0kj[2];
			HYPER1[j*h_num+kn].DpiDq[A_Y]=(S[A_Y][0]+F_dSdC_F[A_Y][0])*n0kj[0]+(S[A_Y][1]+F_dSdC_F[A_Y][1])*n0kj[1]+(S[A_Y][2]+F_dSdC_F[A_Y][2])*n0kj[2];
			HYPER1[j*h_num+kn].DpiDq[A_Z]=(S[A_Z][0]+F_dSdC_F[A_Z][0])*n0kj[0]+(S[A_Z][1]+F_dSdC_F[A_Z][1])*n0kj[1]+(S[A_Z][2]+F_dSdC_F[A_Z][2])*n0kj[2];
		}
		n0jj[A_X]=HYPER1[j*h_num+j].n0ij[A_X];	n0jj[A_Y]=HYPER1[j*h_num+j].n0ij[A_Y];	n0jj[A_Z]=HYPER1[j*h_num+j].n0ij[A_Z];
		HYPER1[j*h_num+j].DpiDq[A_X]=(S[A_X][0]+F_dSdC_F[A_X][0])*n0jj[0]+(S[A_X][1]+F_dSdC_F[A_X][1])*n0jj[1]+(S[A_X][2]+F_dSdC_F[A_X][2])*n0jj[2];
		HYPER1[j*h_num+j].DpiDq[A_Y]=(S[A_Y][0]+F_dSdC_F[A_Y][0])*n0jj[0]+(S[A_Y][1]+F_dSdC_F[A_Y][1])*n0jj[1]+(S[A_Y][2]+F_dSdC_F[A_Y][2])*n0jj[2];
		HYPER1[j*h_num+j].DpiDq[A_Z]=(S[A_Z][0]+F_dSdC_F[A_Z][0])*n0jj[0]+(S[A_Z][1]+F_dSdC_F[A_Z][1])*n0jj[1]+(S[A_Z][2]+F_dSdC_F[A_Z][2])*n0jj[2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	in_Ci[D];
	delete[]	in_Ci;
}


void q_nab_lap(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double **dg, double Dt, double V, double mi, vector<double > NEIw)
{
	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;
	/////////////////èâä˙âª///////////////////
	for(int k=0;k<Nx;k++)
	{
		for(int i=0;i<h_num;i++)
		{
			dg[i][k]=0.;
		}
	}


	/////////////////dE, dgåvéZ///////////////////
	//////////íËã`
	double Dgkk_n[DIMENSION]={0,0,0};
	double Dgki_n[DIMENSION]={0,0,0};

	double Dgik[DIMENSION]={0,0,0};
	double Dgik_n[DIMENSION]={0,0,0};
	double Dgkk[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		/////////////dlam E1
		Dgkk_n[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];	Dgkk_n[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];	Dgkk_n[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];
		/////////////dlam E1
		Dgkk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	Dgkk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];
		/////////////dlam E3
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			Dgki_n[A_X]=HYPER1[k*h_num+in].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[k*h_num+in].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[k*h_num+in].DgDq_n[A_Z];			/////////////dlam E2
			Dgik[A_X]=HYPER1[in*h_num+k].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+k].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+k].DgDq[A_Z];
			Dgik_n[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];	Dgik_n[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];	Dgik_n[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				dg[in][jn]-=0.5*Dt*Dt/mi*(Dgik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+Dgik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+Dgik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);

			}
			dg[in][k]-=0.5*Dt*Dt/mi*(Dgik[A_X]*Dgkk_n[A_X]+Dgik[A_Y]*Dgkk_n[A_Y]+Dgik[A_Z]*Dgkk_n[A_Z]);
			dg[k][in]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgik_n[A_X]+Dgkk[A_Y]*Dgik_n[A_Y]+Dgkk[A_Z]*Dgik_n[A_Z]);
		}
		dg[k][k]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgkk_n[A_X]+Dgkk[A_Y]*Dgkk_n[A_Y]+Dgkk[A_Z]*Dgkk_n[A_Z]);
	}

	for(int k=0;k<Nw;k++)
	{
		int kn=NEIw[k];
		/////////////dlam E1
		Dgkk[A_X]=HYPER1[kn*h_num+kn].DgDq[A_X];	Dgkk[A_Y]=HYPER1[kn*h_num+kn].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[kn*h_num+kn].DgDq[A_Z];
		/////////////dlam E3
		int Nk=HYPER[kn].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[kn].NEI[i];
			Dgik[A_X]=HYPER1[in*h_num+kn].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+kn].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+kn].DgDq[A_Z];
			dg[in][k+h_num]+=0.5*Dt*Dt/mi*Dgik[A_Z];
		}
		dg[kn][k+h_num]+=0.5*Dt*Dt/mi*Dgkk[A_Z];
	}

}

void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double **dg, double *h,double **dh, double *d, int h_num,int count,int t,double E,double mi, double V,vector<double > NEIw)
{

	stringstream ss0;
	ss0<<"./T/T_"<<t<<".csv";
	stringstream ss1;
	ss1<<"./T/dT_"<<t<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./T/rT_"<<t<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./lam/lam_"<<t<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./g/g_"<<t<<"_"<<count<<".csv";
	stringstream ss7;
	ss7<<"./g/dg_"<<t<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./p/hp_"<<t<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./q/q_"<<t<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./T/d_"<<t<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./E/E_"<<t<<".csv";
	//stringstream ss12;
	//ss12<<"./h/h_"<<t<<"_"<<count<<".csv";
	stringstream ss13;
	ss13<<"./g/J_"<<t<<"_"<<count<<".csv";
	stringstream ss14;
	ss14<<"./g/DgDq_"<<t<<"_"<<count<<".csv";
	stringstream ss15;
	ss15<<"./g/stress_"<<t<<"_"<<count<<".csv";
	stringstream ss16;
	ss16<<"./h/h_"<<t<<"_"<<count<<".csv";
	stringstream ss17;
	ss17<<"./h/dh_"<<t<<"_"<<count<<".csv";
	stringstream ss18;
	ss18<<"./mu/mu_"<<t<<"_"<<count<<".csv";


	ofstream f_dt(ss1.str());
	ofstream f_rt(ss2.str());
	ofstream f_lam(ss3.str());
	ofstream f_g(ss6.str());
	ofstream f_dg(ss7.str());
	ofstream f_p(ss8.str());
	ofstream f_q(ss9.str());
	ofstream f_d(ss10.str());
	//ofstream f_h(ss12.str());
	ofstream f_J(ss13.str());
	ofstream f_dgdq(ss14.str());
	ofstream f_s(ss15.str());
	ofstream f_h(ss16.str());
	ofstream f_dh(ss17.str());
	ofstream f_mu(ss18.str());

	if(count==1)
	{
		ofstream f_T(ss0.str(), ios::trunc);
		f_T<<count<<","<<T<<endl;
		f_T.close();
		ofstream f_E(ss11.str(), ios::trunc);
		f_E<<E<<endl;
		f_E.close();
	}	
	else
	{
		ofstream f_T(ss0.str(), ios::app);
		f_T<<count<<","<<T<<endl;
		f_T.close();
		ofstream f_E(ss11.str(), ios::app);
		f_E<<E<<endl;
		f_E.close();
	}

	int Nw=NEIw.size();
	int Nx=Nw+h_num;
	for(int i=0;i<Nw;i++)
	{
		for(int k=0;k<Nx;k++)
		{
			f_dh<<dh[i][k]<<",";
		}
		f_dh<<endl;
		f_h<<h[i]<<endl;
	}

	for(int i=0;i<h_num;i++)
	{
		f_lam<<HYPER[i].lam<<endl;
		f_mu<<HYPER[i].mu<<endl;
		for(int k=0;k<Nx;k++)
		{
			f_dg<<dg[i][k]<<",";
		}
		f_dg<<endl;

		f_g<<g[i]<<endl;
		f_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
		f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		f_d<<d[i]<<endl;
		f_J<<HYPER[i].J<<endl;

		f_dgdq<<i<<endl;
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			f_dgdq<<","<<jn<<","<<HYPER1[i*h_num+jn].DgDq[A_X]<<","<<HYPER1[i*h_num+jn].DgDq[A_Y]<<","<<HYPER1[i*h_num+jn].DgDq[A_Z]<<endl;
		}

		//f_h<<h[i]<<","<<th_h[i]<<endl;
		f_s<<i<<","<<HYPER[i].stress[A_X][0]<<","<<HYPER[i].stress[A_X][1]<<","<<HYPER[i].stress[A_X][2]<<endl;
		f_s<<","<<HYPER[i].stress[A_Y][0]<<","<<HYPER[i].stress[A_Y][1]<<","<<HYPER[i].stress[A_Y][2]<<endl;
		f_s<<","<<HYPER[i].stress[A_Z][0]<<","<<HYPER[i].stress[A_Z][1]<<","<<HYPER[i].stress[A_Z][2]<<endl;
	}
	for(int i=0;i<Nx;i++)
	{
		f_dt<<dT[i]<<endl;

		for(int k=0;k<Nx;k++)
		{
			f_rt<<rT[i*Nx+k]<<",";
		}
		f_rt<<endl;

		f_d<<d[i]<<endl;

	}
	f_dt.close();
	f_rt.close();
	f_lam.close();
	f_dg.close();
	f_g.close();
	f_p.close();
	f_q.close();
	f_d.close();
	//f_h.close();
	f_J.close();
	f_dgdq.close();
	f_s.close();
	f_h.close();
	f_dh.close();
	f_mu.close();


}


void calc_n(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	int h_num=HYPER.size();

	/////////////F, J, t_inverseÇÃåvéZ
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	double J=0.;
	double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double a[DIMENSION]={0,0,0};
	for(int i=0;i<h_num;i++)
	{
		fi[A_X][0]=0.;	fi[A_X][1]=0.;	fi[A_X][2]=0.;
		fi[A_Y][0]=0.;	fi[A_Y][1]=0.;	fi[A_Y][2]=0.;
		fi[A_Z][0]=0.;	fi[A_Z][1]=0.;	fi[A_Z][2]=0.;

		int Ni=HYPER[i].N;	

		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+jn].wiin;
			a[A_X]=HYPER1[i*h_num+jn].aiin[A_X];	a[A_Y]=HYPER1[i*h_num+jn].aiin[A_Y];	a[A_Z]=HYPER1[i*h_num+jn].aiin[A_Z];

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


	/////////////DgDq, Stress, W, S, dSdcåvéZ
	/////////////DgDq, Stress, W, S, dSdcåvéZ
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};


	double trace_dC=0;
	double trace_dC2=0;

	double Ic=0;
	double IIc=0;

	double trace_b=0., trace_bb=0.;
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<h_num;j++)
	{	

		////////WåvéZ
		double J=HYPER[j].J;	
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1./3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1./3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1./3.)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1./3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1./3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1./3.)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1./3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1./3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1./3.)*HYPER[j].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1./3.)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1./3.)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1./3.)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1./3.)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1./3.)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1./3.)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1./3.)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1./3.)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1./3.)*HYPER[j].Fi[2][2];
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

		////////StressåvéZ	
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

		trace_b=b[0][0]+b[1][1]+b[2][2];
		trace_bb=bb[0][0]+bb[1][1]+bb[2][2];

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
	//cout<<"----------OK"<<endl;
}


void p_QP(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt,vector<double > NEIw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1.e-6;
	double ep_min=1.e-6;
	double d_ep=1.e-20;

	double p_p[DIMENSION]={0,0,0};
	double p[DIMENSION]={0,0,0};
	double W=0.;
	double g=0.;
	double r=0.01;

	double T=0.;
	double *dT=new double [Nx];	
	double *rT=new double [Nx*Nx];

	double *d=new double [Nx];
	double *B=new double [Nx*Nx];

	double *G=new double [h_num];
	double **dG=new double *[h_num];
	double *H=new double [h_num];
	double **dH=new double *[h_num];

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		dG[i]=new double[Nx];
		dH[i]=new double[Nx];
	}


	for(int i=0;i<h_num;i++)
	{
		dT[i]=0.;
		G[i]=0.;
		H[i]=0.;
		d[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			dG[i][j]=0.;
			dH[i][j]=0.;
		}
	}
	for(int i=0;i<Nx;i++)
	{
		dT[i]=0.;
		d[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			rT[i*Nx+j]=0.;
			B[i*Nx+j]=0.;
		}
	}

	p_lap(HYPER,HYPER1,dH,dG,Dt,V,mi);

	E=1;
	count=0;
	while(E>ep)
	{
		count++;
		p_variables(HYPER,HYPER1,h_num,Dt,mi);
		p_nab(HYPER,HYPER1,G,H,Dt,V,mi);
		T=0.;
		for(int i=0;i<h_num;i++)
		{	

			T+=G[i]*G[i];
			if(H[i]>0)	T+=H[i]*H[i];
		}

		for(int k=0;k<Nx;k++)
		{	
			dT[k]=0.;
			for(int i=0;i<h_num;i++)
			{	
				dT[k]+=2.*dG[i][k]*G[i];
				if(H[i]>0)	dT[k]+=2.*dH[i][k]*H[i];
			}

			for(int l=0;l<Nx;l++)
			{
				rT[l*Nx+k]=0.;
				for(int i=0;i<h_num;i++)
				{	
					//rT[l*h_num+k]+=r*2.*dG[i*h_num+k]*dG[i*h_num+l]*(3*G[i]*G[i]+th_G[i]);
					rT[l*Nx+k]+=2.*dG[i][k]*dG[i][l];
					if(H[i]>0)	rT[l*Nx+k]+=2.*dH[i][k]*dH[i][l];
					if(H[i]>-1.e-20&&dH[i][k]*dH[i][l]>0)	rT[l*Nx+k]+=2.*dH[i][k]*dH[i][l];
					//rT[l*h_num+k]+=r*dG[i*h_num+k]*dG[i*h_num+l];
				}
			}
		}

		for(int i=0;i<Nx;i++)
		{
			d[i]=dT[i];
			for(int j=0;j<Nx;j++)	B[i*Nx+j]=rT[i*Nx+j];
		}
		gauss(B,d,Nx);
		E_sum=0;
		for(int i=0;i<h_num;i++)
		{
			HYPER[i].lam-=d[i];
			HYPER[i].mu-=d[i+h_num];
			E_sum+=fabs(d[i])+fabs(d[i+h_num]);
		}
		E=E_sum;
		cout<<"E"<<count<<"="<<E<<endl;


		/////////////påvéZ
		double p_p[DIMENSION]={0,0,0};
		double DgDq_ji[DIMENSION]={0,0,0};
		double DgDq_ii[DIMENSION]={0,0,0};
		double lam=0.;
		double mu=0.;
		for(int i=0;i<h_num;i++)
		{
			p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;
			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{
				int jn=HYPER[i].NEI[j];
				lam=HYPER[jn].lam;
				DgDq_ji[A_X]=HYPER1[jn*h_num+i].DgDq[A_X];	DgDq_ji[A_Y]=HYPER1[jn*h_num+i].DgDq[A_Y];	DgDq_ji[A_Z]=HYPER1[jn*h_num+i].DgDq[A_Z];
				p_p[A_X]+=(HYPER[jn].stress[A_X][A_X]-lam)*DgDq_ji[A_X]+HYPER[jn].stress[A_X][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_X][A_Z]*DgDq_ji[A_Z];
				p_p[A_Y]+=HYPER[jn].stress[A_Y][A_X]*DgDq_ji[A_X]+(HYPER[jn].stress[A_Y][A_Y]-lam)*DgDq_ji[A_Y]+HYPER[jn].stress[A_Y][A_Z]*DgDq_ji[A_Z];
				p_p[A_Z]+=HYPER[jn].stress[A_Z][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_Z][A_Y]*DgDq_ji[A_Y]+(HYPER[jn].stress[A_Z][A_Z]-lam)*DgDq_ji[A_Z];
			}
			DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
			lam=HYPER[i].lam;
			mu=HYPER[i].mu;
			p_p[A_X]+=(HYPER[i].stress[A_X][A_X]-lam)*DgDq_ii[A_X]+HYPER[i].stress[A_X][A_Y]	  *DgDq_ii[A_Y]+HYPER[i].stress[A_X][A_Z]			    *DgDq_ii[A_Z];
			p_p[A_Y]+=HYPER[i].stress[A_Y][A_X]		 *DgDq_ii[A_X]+(HYPER[i].stress[A_Y][A_Y]-lam)*DgDq_ii[A_Y]+HYPER[i].stress[A_Y][A_Z]			    *DgDq_ii[A_Z];
			p_p[A_Z]+=HYPER[i].stress[A_Z][A_X]		 *DgDq_ii[A_X]+HYPER[i].stress[A_Z][A_Y]	  *DgDq_ii[A_Y]+(HYPER[i].stress[A_Z][A_Z]-HYPER[i].lam)*DgDq_ii[A_Z];
			HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+0.5*Dt*p_p[A_X];
			HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+0.5*Dt*p_p[A_Y];
			HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*mu;
		}
		//if(count==1||count%100==0)
		{
			//cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
			output_data_p(HYPER,HYPER1,dG,G,dH,H,rT,dT,T,d,Nx,h_num,count,t,E);
		}

	}



	for(int i=0;i<h_num;i++)
	{
		delete[]	dG[i];
		delete[]	dH[i];
	}



	delete[]	G;
	delete[]	dG;
	delete[]	dT;
	delete[]	rT;
	delete[]	H;
	delete[]	dH;
	delete[]	d;
	delete[]	B;
}

void p_nab(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, double *G,double *H, double Dt, double V, double mi)
{

	int h_num=HYPER.size();
	/////////////////dEåvéZ///////////////////
	for(int k=0;k<h_num;k++)
	{
		G[k]=0.;
		H[k]=0.;
	}


	double Dgkk[DIMENSION]={0,0,0};
	double Dgki[DIMENSION]={0,0,0};
	double pk[DIMENSION]={0,0,0};
	double pi[DIMENSION]={0,0,0};

	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		Dgkk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	Dgkk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];
		pk[A_X]=HYPER[k].p[A_X];	pk[A_Y]=HYPER[k].p[A_Y];	pk[A_Z]=HYPER[k].p[A_Z];

		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];

			Dgki[A_X]=HYPER1[k*h_num+in].DgDq[A_X];	Dgki[A_Y]=HYPER1[k*h_num+in].DgDq[A_Y];	Dgki[A_Z]=HYPER1[k*h_num+in].DgDq[A_Z];
			pi[A_X]=HYPER[in].p[A_X];	pi[A_Y]=HYPER[in].p[A_Y];	pi[A_Z]=HYPER[in].p[A_Z];

			G[k]+=1./mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
		}
		G[k]+=1./mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
		H[k]=-1./mi*HYPER[k].p[A_Z];
	}

}


void p_lap(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double **dH, double **dG, double Dt, double V, double mi)
{

	int h_num=HYPER.size();
	/////////////////dEåvéZ///////////////////
	for(int k=0;k<2*h_num;k++)
	{
		for(int i=0;i<h_num;i++)
		{
			dG[i][k]=0.;
			dH[i][k]=0.;
		}
	}

	double Dgjj[DIMENSION]={0,0,0};
	double Dgkj[DIMENSION]={0,0,0};
	double Dgjk[DIMENSION]={0,0,0};
	double Dglj[DIMENSION]={0,0,0};
	for(int j=0;j<h_num;j++)
	{
		int Nj=HYPER[j].N;
		Dgjj[A_X]=HYPER1[j*h_num+j].DgDq[A_X];	Dgjj[A_Y]=HYPER1[j*h_num+j].DgDq[A_Y];	Dgjj[A_Z]=HYPER1[j*h_num+j].DgDq[A_Z];

		for(int k=0;k<Nj;k++)
		{
			int kn=HYPER[j].NEI[k];

			Dgkj[A_X]=HYPER1[kn*h_num+j].DgDq[A_X];	Dgkj[A_Y]=HYPER1[kn*h_num+j].DgDq[A_Y];	Dgkj[A_Z]=HYPER1[kn*h_num+j].DgDq[A_Z];
			Dgjk[A_Z]=HYPER1[j*h_num+kn].DgDq[A_Z];

			for(int l=0;l<Nj;l++)
			{
				int ln=HYPER[j].NEI[l];
				Dglj[A_X]=HYPER1[ln*h_num+j].DgDq[A_X];	Dglj[A_Y]=HYPER1[ln*h_num+j].DgDq[A_Y];	Dglj[A_Z]=HYPER1[ln*h_num+j].DgDq[A_Z];

				dG[kn][ln]-=0.5*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
			}

			dG[kn][j]-=0.5*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			dG[j][kn]-=0.5*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

			dG[j][kn+h_num]+=0.5*Dt/mi*Dgjk[A_Z];
			dH[j][kn]+=0.5*Dt/mi*Dgkj[A_Z];
		}

		dG[j][j]-=0.5*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);
		dG[j][j+h_num]+=0.5*Dt/mi*Dgjj[A_Z];
		dH[j][j]+=0.5*Dt/mi*Dgjj[A_Z];
		dH[j][j+h_num]=-0.5*Dt/mi;
	}

}


void p_variables(vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int h_num,double Dt,double mi)
{
	/////////////påvéZ
	double p_p[DIMENSION]={0,0,0};
	double DgDq_ji[DIMENSION]={0,0,0};
	double DgDq_ii[DIMENSION]={0,0,0};
	double mu=0.;
	double lam=0.;
	for(int i=0;i<h_num;i++)
	{
		p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			lam=HYPER[jn].lam;
			DgDq_ji[A_X]=HYPER1[jn*h_num+i].DgDq[A_X];	DgDq_ji[A_Y]=HYPER1[jn*h_num+i].DgDq[A_Y];	DgDq_ji[A_Z]=HYPER1[jn*h_num+i].DgDq[A_Z];
			p_p[A_X]+=(HYPER[jn].stress[A_X][A_X]-lam)*DgDq_ji[A_X]+HYPER[jn].stress[A_X][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_X][A_Z]*DgDq_ji[A_Z];
			p_p[A_Y]+=HYPER[jn].stress[A_Y][A_X]*DgDq_ji[A_X]+(HYPER[jn].stress[A_Y][A_Y]-lam)*DgDq_ji[A_Y]+HYPER[jn].stress[A_Y][A_Z]*DgDq_ji[A_Z];
			p_p[A_Z]+=HYPER[jn].stress[A_Z][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_Z][A_Y]*DgDq_ji[A_Y]+(HYPER[jn].stress[A_Z][A_Z]-lam)*DgDq_ji[A_Z];
		}
		DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
		lam=HYPER[i].lam;
		mu=HYPER[i].mu;
		p_p[A_X]+=(HYPER[i].stress[A_X][A_X]-lam)*DgDq_ii[A_X]+HYPER[i].stress[A_X][A_Y]	  *DgDq_ii[A_Y]+HYPER[i].stress[A_X][A_Z]			    *DgDq_ii[A_Z];
		p_p[A_Y]+=HYPER[i].stress[A_Y][A_X]		 *DgDq_ii[A_X]+(HYPER[i].stress[A_Y][A_Y]-lam)*DgDq_ii[A_Y]+HYPER[i].stress[A_Y][A_Z]			    *DgDq_ii[A_Z];
		p_p[A_Z]+=HYPER[i].stress[A_Z][A_X]		 *DgDq_ii[A_X]+HYPER[i].stress[A_Z][A_Y]	  *DgDq_ii[A_Y]+(HYPER[i].stress[A_Z][A_Z]-HYPER[i].lam)*DgDq_ii[A_Z];
		HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+0.5*Dt*p_p[A_X];
		HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+0.5*Dt*p_p[A_Y];
		HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*mu;
	}
}



void output_data_p(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double **dG, double *G,double **dH,double *H, double *rT, double *dT, double T, double *d, int Nx, int h_num,int count,int t,double E)
{
	stringstream ss0;
	ss0<<"./g/dG_"<<t<<"_"<<count<<".csv";
	stringstream ss1;
	ss1<<"./g/G_"<<t<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./T/p_rT_"<<t<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./T/p_dT_"<<t<<"_"<<count<<".csv";
	stringstream ss4;
	ss4<<"./T/p_T_"<<t<<".csv";
	stringstream ss5;
	ss5<<"./p/p_"<<t<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./E/p_E_"<<t<<".csv";
	stringstream ss7;
	ss7<<"./T/p_d_"<<t<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./h/dH_"<<t<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./h/H_"<<t<<"_"<<count<<".csv";


	if(count==1)
	{
		ofstream f_t(ss4.str(),ios::trunc);
		f_t<<T<<endl;
		f_t.close();
		ofstream f_e(ss6.str(),ios::trunc);
		f_e<<E<<endl;
		f_e.close();
	}
	else
	{
		ofstream f_t(ss4.str(),ios::app);
		f_t<<T<<endl;
		f_t.close();
		ofstream f_e(ss6.str(),ios::app);
		f_e<<E<<endl;
		f_e.close();
	}

	ofstream f_dG(ss0.str());
	ofstream f_G(ss1.str());
	ofstream f_rt(ss2.str());
	ofstream f_dt(ss3.str());
	ofstream f_p(ss5.str());
	ofstream f_d(ss7.str());

	ofstream f_dH(ss8.str());
	ofstream f_H(ss9.str());


	for(int i=0;i<h_num;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_dG<<dG[i][j]<<",";
			f_dH<<dH[i][j]<<",";
		}
		f_dG<<endl;
		f_dH<<endl;
		f_G<<G[i]<<endl;
		f_H<<H[i]<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
	}
	for(int i=0;i<Nx;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_rt<<rT[i*Nx+j]<<",";
		}
		f_rt<<endl;
		f_dt<<dT[i]<<endl;
		f_d<<d[i]<<endl;
	}
	f_dG.close();
	f_G.close();
	f_rt.close();
	f_dt.close();
	f_p.close();
	f_d.close();
	f_dH.close();
	f_H.close();

}









