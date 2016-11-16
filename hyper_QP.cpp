#include "stdafx.h"	
void calc_nabla_laplacian(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dL,double *rL,double **dg,double **dh,double **d_dgdt,double **d_dhdt,double **rg,double **rh,double **r_dgdt,double **r_dhdt);


void calc_HYPER_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{

	////定義
	int h_num=HYPER.size();
	int Nx=h_num*4;

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();

	double L=0;
	double L_n=0;
	double f=0;
	double T=0;

	double *df=new double [Nx];
	double *dT=new double [Nx];
	double *dL=new double [Nx];

	double *rf=new double [Nx*Nx];
	double *rT=new double [Nx*Nx];
	double *rL=new double [Nx*Nx];

	double *g=new double [h_num];
	double *h=new double [h_num];
	double *dgdt=new double [h_num];
	double *dhdt=new double [h_num];

	double **dg=new double *[h_num];
	double **dh=new double *[h_num];
	double **d_dgdt=new double *[h_num];
	double **d_dhdt=new double *[h_num];

	double **rg=new double *[h_num];
	double **rh=new double *[h_num];
	double **r_dgdt=new double *[h_num];
	double **r_dhdt=new double *[h_num];


	for(int i=0;i<h_num;i++)
	{
		dg[i]=new double [Nx];
		dh[i]=new double [Nx];
		d_dgdt[i]=new double [Nx];
		d_dhdt[i]=new double [Nx];

		rg[i]=new double [Nx*Nx];
		rh[i]=new double [Nx*Nx];
		r_dgdt[i]=new double [Nx*Nx];
		r_dhdt[i]=new double [Nx*Nx];
	}

	double r=1;
	double ep=1e-10;
	double *th_g=new double [h_num];
	double *th_h=new double [h_num];
	double *th_dgdt=new double [h_num];
	double *th_dhdt=new double [h_num];

	//初期化
	for(int i=0;i<h_num;i++)
	{
		L_n+=0.5/mi*(HYPER[i].p_n[A_X]*HYPER[i].p_n[A_X]+HYPER[i].p_n[A_Y]*HYPER[i].p_n[A_Y]+HYPER[i].p_n[A_Z]*HYPER[i].p_n[A_Z])-V*HYPER[i].W_n;

		df[i]=0;
		dT[i]=0;
		dL[i]=0;

		g[i]=0;
		h[i]=0;
		dgdt[i]=0;
		dhdt[i]=0;

		th_g[i]=0;
		th_h[i]=0;
		th_dgdt[i]=0;
		th_dhdt[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dg[i][j]=0;
			dh[i][j]=0;
			d_dgdt[i][j]=0;
			d_dhdt[i][j]=0;
			rf[i*Nx+j]=0;
			rT[i*Nx+j]=0;
			rL[i*Nx+j]=0;
			for(int k=0;k<Nx;k++)
			{
				rg[i][j*Nx+k]=0;
				rh[i][j*Nx+k]=0;
				r_dgdt[i][j*Nx+k]=0;
				r_dhdt[i][j*Nx+k]=0;
			}
		}
	}
	for(int i=h_num;i<Nx;i++)
	{
		df[i]=0;
		dT[i]=0;
		dL[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dg[i][j]=0;
			dh[i][j]=0;
			d_dgdt[i][j]=0;
			d_dhdt[i][j]=0;
			rf[i*Nx+j]=0;
			rT[i*Nx+j]=0;

			rL[i*Nx+j]=0;
			for(int k=0;k<Nx;k++)
			{
				rg[i][j*Nx+k]=0;
				rh[i][j*Nx+k]=0;
				r_dgdt[i][j*Nx+k]=0;
				r_dhdt[i][j*Nx+k]=0;
			}
		}
	}








	////t=n　変数計算結果は既値とする。下記での計算はなし。
	cout<<"QP START----------------------";
	double E_min=1;
	int count_min=0;
	while(E_min>ep)
	{
		count_min++;

		//初期化
		for(int i=0;i<h_num;i++)
		{
			HYPER[i].old_h_lam=HYPER[i].h_lam;
			HYPER[i].old_lam=HYPER[i].lam;
			HYPER[i].old_h_mu=HYPER[i].h_mu;
			HYPER[i].old_mu=HYPER[i].mu;
		}


		double E=1;
		int count=0;
		while(E>ep)
		{
			count++;

			//t=n+1 変数計算





			calc_nabla_laplacian(CON,PART,HYPER,HYPER1,dL,rL,dg,dh,d_dgdt,d_dhdt,rg,rh,r_dgdt,r_dhdt);


			//目的関数計算
			T=0;
			for(int i=0;i<h_num;i++)
			{
				L+=0.5/mi*(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z])-V*HYPER[i].W;
				T+=0.5*r*( (g[i]+th_g[i])*(g[i]+th_g[i]) + (dgdt[i]+th_dgdt[i])*(dgdt[i]+th_dgdt[i]) );
				if(h[i]+th_h[i]>0)	T+=0.5*r*(h[i]+th_h[i])*(h[i]+th_h[i]);
				if(dhdt[i]+th_dhdt[i]>0)	T+=0.5*r*(dhdt[i]+th_dhdt[i])*(dhdt[i]+th_dhdt[i]);
			}			
			f=(L-L_n)*(L-L_n);
			T=+f;

			for(int i=0;i<Nx;i++)
			{
				dT[i]=0;

				df[i]=2*dL[i]*(L-L_n);
				dT[i]+=df[i];

				for(int j=0;j<h_num;j++)
				{
					dT[i]+=r*( dg[j][i]*(g[j]+th_g[j]) + d_dgdt[j][i]*(dgdt[j]+th_dgdt[j]) );
					if(h[j]+th_h[j]>0)	dT[i]+=r*dh[j][i]*(h[j]+th_h[j]);
					if(dhdt[j]+th_dhdt[j]>0)	dT[i]+=r*d_dhdt[j][i]*(dhdt[j]+th_dhdt[j]);
				}

				for(int j=0;j<Nx;j++)
				{
					rT[i*Nx+j]=0;
					
					rf[i*Nx+j]=2*(rL[i*Nx+j]*(L-L_n)+dL[i]*dL[j]);
					rT[i*Nx+j]+=rf[i*Nx+j];
					for(int k=0;k<h_num;k++)
					{
						rT[i*Nx+j]+=r*( rg[k][i*Nx+j]*(g[k]+th_g[k])+dg[k][i]*dg[k][j] + r_dgdt[k][i*Nx+j]*(dgdt[k]+th_dgdt[k])+d_dgdt[k][i]*d_dgdt[k][j] );
					
						double p_rh=rh[k][i*Nx+j]*(h[k]+th_h[k])+dh[k][i]*dh[k][j];
						if(h[j]+th_h[j]>0)	rT[i*Nx+j]+=r*p_rh;
						else if(p_rh>0 && (h[j]+th_h[j]==0))	rT[i*Nx+j]+=r*p_rh;

						double p_rdhdt=r_dhdt[k][i*Nx+j]*(dhdt[k]+th_dhdt[k])+d_dhdt[k][i]*d_dhdt[k][j];
						if(dhdt[j]+th_dhdt[j]>0)	rT[i*Nx+j]+=r*p_rdhdt;
						else if(p_rdhdt>0 && (dhdt[j]+th_dhdt[j]==0))	rT[i*Nx+j]+=r*p_rdhdt;
					}
				}
			}

			gauss(rT,dT,Nx);

			for(int i=0;i<h_num;i++)
			{
				HYPER[i].h_lam-=dT[i];
				HYPER[i].lam-=dT[i+h_num];
				HYPER[i].h_mu-=dT[i+2*h_num];
				HYPER[i].mu-=dT[i+3*h_num];
			}	

			double p_E=0;
			for(int i=0;i<Nx;i++)	p_E+=dT[i]*dT[i];
			E=sqrt(p_E);
		}

		double p_E_min=0;
		for(int i=0;i<h_num;i++)
		{
			p_E_min+=(HYPER[i].old_h_lam-HYPER[i].h_lam)*(HYPER[i].old_h_lam-HYPER[i].h_lam)+(HYPER[i].old_lam-HYPER[i].lam)*(HYPER[i].old_lam-HYPER[i].lam)
					 +(HYPER[i].old_h_mu-HYPER[i].h_mu)*(HYPER[i].old_h_mu-HYPER[i].h_mu)+(HYPER[i].old_mu-HYPER[i].mu)*(HYPER[i].old_mu-HYPER[i].mu);
		}
		E_min=sqrt(p_E_min);

		if(E_min<ep*1000)	r*=4;
		for(int i=0;i<h_num;i++)
		{
			th_g[i]+=g[i];
			th_dgdt[i]+=dgdt[i];
			if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
			if(dhdt[i]+th_dhdt[i]>0)	th_dhdt[i]+=dhdt[i];
		}
	}
	cout<<"OK"<<endl;


	
	//メモリ解放
	delete[]	g;
	delete[]	h;
	delete[]	dgdt;
	delete[]	dhdt;

	delete[]	df;
	delete[]	dT;
	delete[]	dL;
	delete[]	rf;
	delete[]	rT;
	delete[]	rL;

	for(int i=0;i<h_num;i++)
	{
		delete[]	dg[i];
		delete[]	dh[i];
		delete[]	d_dgdt[i];
		delete[]	d_dhdt[i];
		delete[]	rg[i];
		delete[]	rh[i];
		delete[]	r_dgdt[i];
		delete[]	r_dhdt[i];
	}
	delete[]	dg;
	delete[]	dh;
	delete[]	d_dgdt;
	delete[]	d_dhdt;
	delete[]	rg;
	delete[]	rh;
	delete[]	r_dgdt;
	delete[]	r_dhdt;

	delete[]	th_g;
	delete[]	th_h;
	delete[]	th_dgdt;
	delete[]	th_dhdt;


}

void calc_nabla_laplacian(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dL,double *rL,double **dg,double **dh,double **d_dgdt,double **d_dhdt,double **rg,double **rh,double **r_dgdt,double **r_dhdt)
{

}
