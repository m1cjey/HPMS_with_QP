#include "stdafx.h"	
void calc_nabla_laplacian(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dL,double *rL,double **dg,double **dh,double **d_dgdt,double **d_dhdt);
void calc_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void calc_HYPER_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{

	////////////定義///////////////
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

	////////////初期化算///////////////
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








	////////////QP計算///////////////		
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

			////////////////////////////////t=n+1 変数計算//////////////////////////////	/////t=n　変数計算結果は既値とする。下記での計算はなし。






			calc_nabla_laplacian(CON,PART,HYPER,HYPER1,dL,rL,dg,dh,d_dgdt,d_dhdt);


			////////////////////////////////目的関数計算/////////////////////////////////
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
				//勾配計算
				dT[i]=0;
				df[i]=2*dL[i]*(L-L_n);
				dT[i]+=df[i];
				for(int j=0;j<h_num;j++)
				{
					dT[i]+=r*( dg[j][i]*(g[j]+th_g[j]) + d_dgdt[j][i]*(dgdt[j]+th_dgdt[j]) );
					if(h[j]+th_h[j]>0)	dT[i]+=r*dh[j][i]*(h[j]+th_h[j]);
					if(dhdt[j]+th_dhdt[j]>0)	dT[i]+=r*d_dhdt[j][i]*(dhdt[j]+th_dhdt[j]);
				}

				//発散計算
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


	
	///////////メモリ解放///////////////		
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

void calc_nabla_laplacian(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dL,double *rL,double **dg,double **dh,double **d_dgdt,double **d_dhdt)
{

	int h_num=HYPER.size();
	int Nx=4*h_num;

	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double nG[DIMENSION]={0,0,1};


	/////////////////dL計算///////////////////
	for(int k=0;k<h_num;k++)
	{
		dL[k]=0;
		dL[k+h_num]=0;
		dL[k+2*h_num]=0;
		dL[k+3*h_num]=0;

		double dL_k=0, dL_p=0;
		double dL2_k=0;
		double dL3_p=0;


		for(int i=0;i<h_num;i++)
		{			
			dL_k+=HYPER[i].p[A_X]*HYPER1[k*h_num+i].DgDq_n[A_X]+HYPER[i].p[A_Y]*HYPER1[k*h_num+i].DgDq_n[A_Y]+HYPER[i].p[A_Z]*HYPER1[k*h_num+i].DgDq_n[A_Z];
			dL2_k+=HYPER[i].p[A_X]*HYPER1[k*h_num+i].DgDq[A_X]+HYPER[i].p[A_Y]*HYPER1[k*h_num+i].DgDq[A_Y]+HYPER[i].p[A_Z]*HYPER1[k*h_num+i].DgDq[A_Z];

			dg[i][k]=0;
			dg[i][k+2*h_num]=0;
			d_dgdt[i][k]=0;
			d_dgdt[i][k+h_num]=0;
			d_dgdt[i][k+2*h_num]=0;
			d_dgdt[i][k+3*h_num]=0;
			dh[i][k]=0;
			dh[i][k+2*h_num]=0;
			d_dhdt[i][k]=0;
			d_dhdt[i][k+h_num]=0;
			d_dhdt[i][k+2*h_num]=0;
			d_dhdt[i][k+3*h_num]=0;
			for(int j=0;j<h_num;j++)
			{
				dL_p+=Dt*( (HYPER[i].stress[A_X][0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+j].DgDq[2])*HYPER1[k*h_num+j].DgDq_n[A_X]
				+(HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+j].DgDq[2])*HYPER1[k*h_num+j].DgDq_n[A_Y]
				+(HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+j].DgDq[2])*HYPER1[k*h_num+j].DgDq_n[A_Z] );
				dL3_p+=Dt*( (HYPER[i].stress[A_X][0]*HYPER1[i*h_num+k].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+k].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+k].DgDq[2])*nG[A_X]
				+(HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+k].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+k].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+k].DgDq[2])*nG[A_Y]
				+(HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+k].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+k].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+k].DgDq[2])*nG[A_Z] );
	
				dg[i][k]+=-0.5*Dt*Dt/mi*(HYPER1[i*h_num+j].DgDq[A_X]*HYPER1[k*h_num+j].DgDq_n[A_X]+HYPER1[i*h_num+j].DgDq[A_Y]*HYPER1[k*h_num+j].DgDq_n[A_Y]+HYPER1[i*h_num+j].DgDq[A_Z]*HYPER1[k*h_num+j].DgDq_n[A_Z]);				
				d_dgdt[i][k]+=-0.5*Dt/mi*(HYPER1[i*h_num+j].DgDq[A_X]*HYPER1[k*h_num+j].DgDq_n[A_X]+HYPER1[i*h_num+j].DgDq[A_Y]*HYPER1[k*h_num+j].DgDq_n[A_Y]+HYPER1[i*h_num+j].DgDq[A_Z]*HYPER1[k*h_num+j].DgDq_n[A_Z]);
				d_dgdt[i][k+h_num]+=-0.5*Dt/mi*(HYPER1[i*h_num+j].DgDq[A_X]*HYPER1[k*h_num+j].DgDq[A_X]+HYPER1[i*h_num+j].DgDq[A_Y]*HYPER1[k*h_num+j].DgDq[A_Y]+HYPER1[i*h_num+j].DgDq[A_Z]*HYPER1[k*h_num+j].DgDq[A_Z]);				
			}	
			dg[i][k+2*h_num]=0.5*Dt*Dt/mi*(HYPER1[i*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[i*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[i*h_num+k].DgDq[A_Z]*nG[A_Z]);
			d_dgdt[i][k+2*h_num]=0.5*Dt/mi*(HYPER1[i*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[i*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[i*h_num+k].DgDq[A_Z]*nG[A_Z]);
			d_dgdt[i][k+3*h_num]=0.5*Dt/mi*(HYPER1[i*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[i*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[i*h_num+k].DgDq[A_Z]*nG[A_Z]);
			dh[i][k]=0.5*Dt*Dt/mi*(HYPER1[k*h_num+i].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*nG[A_Z]);
			dh[i][k+2*h_num]=-0.5*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
			d_dhdt[i][k]=0.5*Dt/mi*(HYPER1[k*h_num+i].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*nG[A_Z]);
			d_dhdt[i][k+h_num]=0.5*Dt/mi*(HYPER1[k*h_num+i].DgDq[A_X]*nG[A_X]+HYPER1[k*h_num+i].DgDq[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].DgDq[A_Z]*nG[A_Z]);
			d_dhdt[i][k+2*h_num]=-0.5*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
			d_dhdt[i][k+3*h_num]=-0.5*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
		}
		dL[k]=-0.5*Dt/mi*(dL_k+dL_p);
		dL[k+h_num]=-0.5*Dt/mi*dL2_k;
		dL[k+2*h_num]=0.5*Dt/mi*(HYPER[k].p[A_X]*nG[A_X]+HYPER[k].p[A_Y]*nG[A_Y]+HYPER[k].p[A_Z]*nG[A_Z]+dL3_p);
		dL[k+3*h_num]=0.5*Dt/mi*(HYPER[k].p[A_X]*nG[A_X]+HYPER[k].p[A_Y]*nG[A_Y]+HYPER[k].p[A_Z]*nG[A_Z]);
	}





	/////////////////rL計算///////////////////
	for(int l=0;l<h_num;l++)
	{
		for(int k=0;k<h_num;k++)
		{
			rL[l*Nx+k]=0;	rL[(l+h_num)*Nx+k]=0;	rL[(l+2*h_num)*Nx+k]=0;	rL[(l+3*h_num)*Nx+k]=0;
			rL[l*Nx+k+h_num]=0;	rL[(l+h_num)*Nx+k+h_num]=0;	rL[(l+2*h_num)*Nx+k+h_num]=0;	rL[(l+3*h_num)*Nx+k+h_num]=0;
			
			double rL_k=0, rL_p=0;
			double rL10_k=0;
			double rL20_p=0;

			double rL01_k=0, rL01_p=0;
			double rL11_k=0;
			double rL21_p=0;

			double rL02_p=0;
			double rL22_p=0;
			for(int i=0;i<h_num;i++)
			{
				////////rL_h_lam_h_lam
				rL_k+=HYPER1[k*h_num+i].DgDq_n[A_X]*HYPER1[l*h_num+i].DgDq_n[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*HYPER1[l*h_num+i].DgDq_n[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*HYPER1[l*h_num+i].DgDq_n[A_Z];
				////////rL_lam_h_lam
				rL10_k+=HYPER1[k*h_num+i].DgDq_n[A_X]*HYPER1[l*h_num+i].DgDq[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*HYPER1[l*h_num+i].DgDq[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*HYPER1[l*h_num+i].DgDq[A_Z];
		
				////////rL_h_lam_lam
				rL01_k+=HYPER1[l*h_num+i].DgDq_n[A_X]*HYPER1[k*h_num+i].DgDq[A_X]+HYPER1[l*h_num+i].DgDq_n[A_Y]*HYPER1[k*h_num+i].DgDq[A_Y]+HYPER1[l*h_num+i].DgDq_n[A_Z]*HYPER1[k*h_num+i].DgDq[A_Z];
				////////rL_lam_lam
				rL11_k+=HYPER1[l*h_num+i].DgDq[A_X]*HYPER1[k*h_num+i].DgDq[A_X]+HYPER1[l*h_num+i].DgDq[A_Y]*HYPER1[k*h_num+i].DgDq[A_Y]+HYPER1[l*h_num+i].DgDq[A_Z]*HYPER1[k*h_num+i].DgDq[A_Z];
				////////rL_h_mu_h_mu
				double p_rL22_p=( (HYPER[i].S[A_X][A_X]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[l*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_X][A_Y]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[l*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_X][A_Z]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[l*h_num+i].n0ij[A_Z] )
						*nG[A_X]+
			
						( (HYPER[i].S[A_Y][A_X]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[l*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Y][A_Y]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[l*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Y][A_Z]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[l*h_num+i].n0ij[A_Z] )
						*nG[A_Y]+
							
						( (HYPER[i].S[A_Z][A_X]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[l*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Z][A_Y]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[l*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Z][A_Z]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[l*h_num+i].n0ij[A_Z] )
						*nG[A_Z];
				rL22_p+=p_rL22_p*(HYPER1[k*h_num+i].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+i].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].n0ij[A_Z]*nG[A_Z]);

				for(int j=0;j<h_num;j++)
				{
					////////rL_h_lam_h_lam
					double p_rL_p=0;
					double p_rg=0;
					for(int m=0;m<h_num;m++)
					{
						p_rL_p+=( (HYPER[i].S[A_X][A_X]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[m*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_X][A_Y]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[m*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_X][A_Z]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[m*h_num+i].n0ij[A_Z] )
						*HYPER1[l*h_num+m].DgDq_n[A_X]
					
						+( (HYPER[i].S[A_Y][A_X]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[m*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Y][A_Y]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[m*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Y][A_Z]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[m*h_num+i].n0ij[A_Z] )
						*HYPER1[l*h_num+m].DgDq_n[A_Y]
						
						+( (HYPER[i].S[A_Z][A_X]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[m*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Z][A_Y]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[m*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Z][A_Z]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[m*h_num+i].n0ij[A_Z] )
						*HYPER1[l*h_num+m].DgDq_n[A_Z];

					}
					rL_p+=p_rL_p*(HYPER1[j*h_num+i].n0ij[A_X]*HYPER1[k*h_num+j].DgDq_n[A_X]+HYPER1[j*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+j].DgDq_n[A_Y]+HYPER1[j*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+j].DgDq_n[A_Z]);
				
					////////rL_h_mu_h_lam
					double p_rL20_p=( (HYPER[i].S[A_X][A_X]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[l*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_X][A_Y]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[l*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_X][A_Z]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[l*h_num+i].n0ij[A_Z] )
						*nG[A_X]+
			
						( (HYPER[i].S[A_Y][A_X]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[l*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Y][A_Y]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[l*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Y][A_Z]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[l*h_num+i].n0ij[A_Z] )
						*nG[A_Y]+
							
						( (HYPER[i].S[A_Z][A_X]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[l*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Z][A_Y]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[l*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Z][A_Z]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[l*h_num+i].n0ij[A_Z] )
						*nG[A_Z];
						rL20_p+=p_rL20_p*(HYPER1[j*h_num+i].n0ij[A_X]*HYPER1[k*h_num+j].DgDq_n[A_X]+HYPER1[j*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+j].DgDq_n[A_Y]+HYPER1[j*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+j].DgDq_n[A_Z]);

					////////rL_h_lam_h_mu
					double p_rL02_p=( (HYPER[i].S[A_X][A_X]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[j*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_X][A_Y]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[j*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_X][A_Z]+2*(HYPER[i].dSdc[A_X][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_X][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_X][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[j*h_num+i].n0ij[A_Z] )
						*HYPER1[l*h_num+j].DgDq_n[A_X]+
			
						( (HYPER[i].S[A_Y][A_X]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[j*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Y][A_Y]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[j*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Y][A_Z]+2*(HYPER[i].dSdc[A_Y][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Y][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Y][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[j*h_num+i].n0ij[A_Z] )
						*HYPER1[l*h_num+j].DgDq_n[A_Y]+
							
						( (HYPER[i].S[A_Z][A_X]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_X]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_X]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_X]))*HYPER1[j*h_num+i].n0ij[A_X]
						+(HYPER[i].S[A_Z][A_Y]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Y]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Y]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Y]))*HYPER1[j*h_num+i].n0ij[A_Y]
						+(HYPER[i].S[A_Z][A_Z]+2*(HYPER[i].dSdc[A_Z][A_X]*HYPER[i].Fi[A_X][A_Z]+HYPER[i].dSdc[A_Z][A_Y]*HYPER[i].Fi[A_Y][A_Z]+HYPER[i].dSdc[A_Z][A_Z]*HYPER[i].Fi[A_Z][A_Z]))*HYPER1[j*h_num+i].n0ij[A_Z] )
						*HYPER1[l*h_num+j].DgDq_n[A_Z];
					rL02_p+=p_rL02_p*(HYPER1[k*h_num+i].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+i].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].n0ij[A_Z]*nG[A_Z]);
				}
			}
			rL[l*Nx+k]=0.25*Dt*Dt/mi*(rL_k-rL_p); //rL_h_lam_h_lam
			rL[(l+h_num)*Nx+k]=0.25*Dt*Dt/mi*rL10_k; //rL_lam_h_Lam
			rL[(l+2*h_num)*Nx+k]=-0.25*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq_n[A_Z]*nG[A_Z]-Dt*Dt/mi/V*rL20_p); //rL_h_mu_h_lam
			rL[(l+3*h_num)*Nx+k]=-0.25*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq_n[A_Z]*nG[A_Z]); //rL_mu_h_lam

			rL[l*Nx+k+h_num]=0.25*Dt*Dt/mi*rL10_k; //rL_h_lam_lam
			rL[(l+h_num)*Nx+k+h_num]=0.25*Dt*Dt/mi*rL11_k; //rL_lam_lam
			rL[(l+2*h_num)*Nx+k+h_num]=-0.25*Dt*Dt/mi*(nG[A_X]*HYPER1[k*h_num+l].DgDq[A_X]+nG[A_Y]*HYPER1[k*h_num+l].DgDq[A_Y]+nG[A_Z]*HYPER1[k*h_num+l].DgDq[A_Z]); //rL_h_mu_lam
			rL[(l+3*h_num)*Nx+k+h_num]=-0.25*Dt*Dt/mi*(nG[A_X]*HYPER1[k*h_num+l].DgDq[A_X]+nG[A_Y]*HYPER1[k*h_num+l].DgDq[A_Y]+nG[A_Z]*HYPER1[k*h_num+l].DgDq[A_Z]); //rL_mu_lam

			rL[l*Nx+k+2*h_num]=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*nG[A_Z]-Dt*Dt/mi/V*rL02_p); //rL_h_lam_h_mu
			rL[(l+h_num)*Nx+k+2*h_num]=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq[A_Z]*nG[A_Z]); //rL_lam_h_mu
			rL[(l+2*h_num)*Nx+k+2*h_num]=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]-Dt*Dt/mi/V*rL22_p); //rL_h_mu_h_mu
			rL[(l+3*h_num)*Nx+k+2*h_num]=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]); //rL_mu_h_mu

			rL[l*Nx+k+3*h_num]=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*nG[A_Z]); //rL_h_lam_mu
			rL[(l+h_num)*Nx+k+3*h_num]=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq[A_Z]*nG[A_Z]); //rL_lam_mu
			rL[(l+2*h_num)*Nx+k+3*h_num]=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]); //rL_h_mu_mu
			rL[(l+3*h_num)*Nx+k+3*h_num]=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]); //rL_mu_mu
		}
	}
}



void calc_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double nG[DIMENSION])
{
	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

	//stressとDgDqを計算

	for(int i=0;i<h_num;i++)
	{
		double p_r[DIMENSION]={0,0,0};
		double p_p[DIMENSION]={0,0,0};
		for(int j=0;j<h_num;j++)
		{
			p_r[A_X]=(HYPER[j].stress_n[A_X][A_X]-HYPER[j].h_lam)*HYPER1[j*h_num+i].DgDq_n[A_X]+HYPER[j].stress_n[A_X][A_Y]*HYPER1[j*h_num+i].DgDq_n[A_Y]+HYPER[j].stress_n[A_X][A_Z]*HYPER1[j*h_num+i].DgDq_n[A_Z];
			p_r[A_Y]=HYPER[j].stress_n[A_Y][A_X]*HYPER1[j*h_num+i].DgDq_n[A_X]+(HYPER[j].stress_n[A_Y][A_Y]-HYPER[j].h_lam)*HYPER1[j*h_num+i].DgDq_n[A_Y]+HYPER[j].stress_n[A_Y][A_Z]*HYPER1[j*h_num+i].DgDq_n[A_Z];
			p_r[A_Z]=HYPER[j].stress_n[A_Z][A_X]*HYPER1[j*h_num+i].DgDq_n[A_X]+HYPER[j].stress_n[A_Z][A_Y]*HYPER1[j*h_num+i].DgDq_n[A_Y]+(HYPER[j].stress_n[A_Z][A_Z]-HYPER[j].h_lam)*HYPER1[j*h_num+i].DgDq_n[A_Z];
				
			p_p[A_X]=(HYPER[j].stress[A_X][A_X]-HYPER[j].lam)*HYPER1[j*h_num+i].DgDq_n[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*h_num+i].DgDq_n[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*h_num+i].DgDq[A_Z];
			p_p[A_Y]=HYPER[j].stress[A_Y][A_X]*HYPER1[j*h_num+i].DgDq_n[A_X]+(HYPER[j].stress[A_Y][A_Y]-HYPER[j].lam)*HYPER1[j*h_num+i].DgDq_n[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*h_num+i].DgDq[A_Z];
			p_p[A_Z]=HYPER[j].stress[A_Z][A_X]*HYPER1[j*h_num+i].DgDq_n[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*h_num+i].DgDq_n[A_Y]+(HYPER[j].stress[A_Z][A_Z]-HYPER[j].lam)*HYPER1[j*h_num+i].DgDq[A_Z];
		}
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*( HYPER[i].p_n[A_X]+0.5*Dt*(p_r[A_X]+nG[A_X]*HYPER[i].h_mu) );
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*( HYPER[i].p_n[A_Y]+0.5*Dt*(p_r[A_Y]+nG[A_Y]*HYPER[i].h_mu) );
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*( HYPER[i].p_n[A_Z]+0.5*Dt*(p_r[A_Z]+nG[A_Z]*HYPER[i].h_mu) );
		HYPER[i].p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*( p_r[A_X]+p_p[A_X]+nG[A_X]*(HYPER[i].h_mu+HYPER[i].mu) );
		HYPER[i].p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*( p_r[A_Y]+p_p[A_Y]+nG[A_Y]*(HYPER[i].h_mu+HYPER[i].mu) );
		HYPER[i].p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*( p_r[A_Z]+p_p[A_Z]+nG[A_Z]*(HYPER[i].h_mu+HYPER[i].mu) );
	}


}


















}