#include "stdafx.h"
int MM_method();
void QP2(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
void calc_wg(mpsconfig &CON, vector<mpselastic> PART, vector<hyperelastic> &HYPER, vector<hyperelastic2> &HYPER1, double **dq, double *w, double *g, double *J, double **F, double **ti_F, double **s);

void QP(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	int h_num=HYPER.size();

	double **dp=new double *[h_num];
	double **old_dp=new double *[h_num];

	double **dq=new double *[h_num];
	double **old_dq=new double *[h_num];

	double Tr=0;
	double **dpTr=new double *[h_num];
	double **dqTr=new double *[h_num];

	double *J=new double [h_num];
	double **F=new double *[h_num];
	double **ti_F=new double *[h_num];
	double **s=new double *[h_num];

	for(int i=0;i<h_num;i++)
	{
		dp[i]=new double [DIMENSION];
		old_dp[i]=new double [DIMENSION];
	
		dq[i]=new double [DIMENSION];
		old_dq[i]=new double [DIMENSION];

		dpTr[i]=new double [DIMENSION];
		dqTr[i]=new double [DIMENSION];

		F[i]=new double [DIMENSION*DIMENSION];
		ti_F[i]=new double [DIMENSION*DIMENSION];

		s[i]=new double [DIMENSION*DIMENSION];
	}

	double V=get_volume(&CON);
	double mi=CON.get_h_dis()*V;

	double *w=new double [h_num];
	double *g=new double [h_num];
	double *h=new double [h_num];

	double *seta_g=new double [h_num];
	double *seta_h=new double [h_num];

	double a[DIMENSION]={0,0,0};
	double n[DIMENSION]={0,0,1};

	double *B=new double [36*h_num*h_num];
	double *d=new double [6*h_num];	
	double *Nr=new double [6*h_num];


	//èâä˙âª
	for(int i=0; i<h_num; i++)
	{
		dp[i][A_X]=0;	dp[i][A_Y]=0;	dp[i][A_Z]=0;
		dq[i][A_X]=0;	dq[i][A_Y]=0;	dq[i][A_Z]=0;
		
		old_dp[i][A_X]=0;	old_dp[i][A_Y]=0;	old_dp[i][A_Z]=0;
		old_dq[i][A_X]=0;	old_dq[i][A_Y]=0;	old_dq[i][A_Z]=0;
		
		dpTr[i][A_X]=0;	dpTr[i][A_Y]=0;	dpTr[i][A_Z]=0;
		dqTr[i][A_X]=0;	dqTr[i][A_Y]=0;	dqTr[i][A_Z]=0;

		F[i][A_X*DIMENSION+A_X]=0;	F[i][A_X*DIMENSION+A_Y]=0;	F[i][A_X*DIMENSION+A_Z]=0;	
		F[i][A_Y*DIMENSION+A_X]=0;	F[i][A_Y*DIMENSION+A_Y]=0;	F[i][A_Y*DIMENSION+A_Z]=0;	
		F[i][A_Z*DIMENSION+A_X]=0;	F[i][A_Z*DIMENSION+A_Y]=0;	F[i][A_Z*DIMENSION+A_Z]=0;	

		ti_F[i][A_X*DIMENSION+A_X]=0;	ti_F[i][A_X*DIMENSION+A_Y]=0;	ti_F[i][A_X*DIMENSION+A_Z]=0;	
		ti_F[i][A_Y*DIMENSION+A_X]=0;	ti_F[i][A_Y*DIMENSION+A_Y]=0;	ti_F[i][A_Y*DIMENSION+A_Z]=0;	
		ti_F[i][A_Z*DIMENSION+A_X]=0;	ti_F[i][A_Z*DIMENSION+A_Y]=0;	ti_F[i][A_Z*DIMENSION+A_Z]=0;	

		s[i][A_X*DIMENSION+A_X]=0;	s[i][A_X*DIMENSION+A_Y]=0;	s[i][A_X*DIMENSION+A_Z]=0;	
		s[i][A_Y*DIMENSION+A_X]=0;	s[i][A_Y*DIMENSION+A_Y]=0;	s[i][A_Y*DIMENSION+A_Z]=0;	
		s[i][A_Z*DIMENSION+A_X]=0;	s[i][A_Z*DIMENSION+A_Y]=0;	s[i][A_Z*DIMENSION+A_Z]=0;	


		w[i]=0;
		g[i]=0;
		h[i]=0;
		seta_g[i]=0;
		seta_h[i]=0;
		J[i]=0;
	}
	for(int i=0;i<h_num*6;i++)
	{
		for(int j=0;j<h_num*6;j++)
		{
			if(j==i)	B[i*h_num*6+j]=1;
			else
			{
				B[i*h_num*6+j]=0;
			}			
		}
		d[i]=0;
		Nr[i]=0;
	}




	double r=0;
	double ep=1e-10;

	double E_min=1;
	int count_min=0;

	while(E_min>ep)
	{
		count_min++;

		for(int i=0;i<h_num;i++)
		{
			old_dp[i][A_X]=dp[i][A_X];	old_dp[i][A_Y]=dp[i][A_Y];	old_dp[i][A_Z]=dp[i][A_Z];
			old_dq[i][A_X]=dq[i][A_X];	old_dq[i][A_Y]=dq[i][A_Y];	old_dq[i][A_Z]=dq[i][A_Z];
		}
		for(int i=0;i<h_num*6;i++)
		{
			for(int j=0;j<h_num*6;j++)
			{
				if(j==i)	B[i*h_num*6+j]=1;
				else
				{
					B[i*h_num*6+j]=0;
				}			
			}
			d[i]=0;
			Nr[i]=0;
		}

		double E=1;
		int count=0;

		calc_wg(CON, PART, HYPER, HYPER1, dq, w, g, J, F, ti_F, s);

		Tr=0;
		for(int i=0;i<h_num;i++)
		{
			Tr+=0.5/mi*( (HYPER[i].p[A_X]+dp[i][A_X])*(HYPER[i].p[A_X]+dp[i][A_X]) + (HYPER[i].p[A_Y]+dp[i][A_Y])*(HYPER[i].p[A_Y]+dp[i][A_Y]) + (HYPER[i].p[A_Z]+dp[i][A_Z])*(HYPER[i].p[A_Z]+dp[i][A_Z]) )
				-V*w[i]+0.5*r*(g[i]+seta_g[i])*(g[i]+seta_g[i]);			

			h[i]=-V*( (PART[i].r[A_X]+dp[i][A_X] - a[A_X])*n[A_X] + (PART[i].r[A_Y]+dp[i][A_Y] - a[A_Y])*n[A_Y] +(PART[i].r[A_Z]+dp[i][A_Z] - a[A_Z])*n[A_Z] );	
			if(h[i]+seta_h[i]>0)	Tr+=0.5*r*(h[i]+seta_h[i])*(h[i]+seta_h[i]);
		}
		for(int i=0;i<h_num;i++)
		{

			dpTr[i][A_X]=1/mi*(HYPER[i].p[A_X]+dp[i][A_X]);
			dpTr[i][A_Y]=1/mi*(HYPER[i].p[A_X]+dp[i][A_Y]);
			dpTr[i][A_Z]=1/mi*(HYPER[i].p[A_X]+dp[i][A_Z]);

			dqTr[i][A_X]=0;	dqTr[i][A_Y]=0;	dqTr[i][A_Z]=0;

			for(int j=0;j<h_num;j++)	dqTr[j][A_X]=-r*V*J[j]*


			for(int j=0;j<Nx;j++)
			{
				dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_X]*(-g[j]+seta[j]);

				dTr[i+Nx]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_Y]*(-g[j]+seta[j]);
				dTr[i+2*Nx]+=-V*(HYPER[j].stress[A_Z][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Z][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_Z]*(-g[j]+seta[j]);			
			}

			h[i]=-V*( (qx[i]-ax)*nx + (qy[i]-ay)*ny +(qz[i]-az)*nz );	
			if(h[i]+seta[i+Nc/2]>0)
			{
				Tr+=0.5*r*(h[i]+seta[i+Nc/2])*(h[i]+seta[i+Nc/2]);
				dTr[i]+=r*(-V*nx+seta[i+Nc/2]);
				dTr[i+Nx]+=r*(-V*ny+seta[i+Nc/2]);
				dTr[i+2*Nx]+=r*(-V*nz+seta[i+Nc/2]);
			}

			dTr[i+Nx*3]=1/mi*px[i];
			dTr[i+Nx*4]=1/mi*py[i];
			dTr[i+Nx*5]=1/mi*pz[i];

			E+=sqrt(dTr[i]*dTr[i] + dTr[i+Nx]*dTr[i+Nx] + dTr[i+Nx*2]*dTr[i+Nx*2] + dTr[i+Nx*3]*dTr[i+Nx*3] + dTr[i+Nx*4]*dTr[i+Nx*4] + dTr[i+Nx*5]*dTr[i+Nx*5]);
		}
		cout<<"E0="<<E<<endl;


		if(E<ep)	break;
		else
		{
			double *qx_k=new double [Nx];
			double *qy_k=new double [Nx];
			double *qz_k=new double [Nx];
			double *px_k=new double [Nx];
			double *py_k=new double [Nx];
			double *pz_k=new double [Nx];
			double *dTr_k=new double [Nx*6];

			double *qx_a=new double [Nx];
			double *qy_a=new double [Nx];
			double *qz_a=new double [Nx];
			double *px_a=new double [Nx];
			double *py_a=new double [Nx];
			double *pz_a=new double [Nx];
			double *w_a=new double [Nx];
			double *g_a=new double [Nx];

			double *s=new double [Nx*6];
			double *y=new double [Nx*6];
			double *sB=new double [Nx*6];
			double *Bs=new double [Nx*6];
			double *BssB=new double [Nx*Nx*36];


			while(E>ep)
			{
				count++;

				for(int i=0;i<Nx;i++)
				{
					qx_k[i]=qx[i];
					qy_k[i]=qy[i];
					qz_k[i]=qz[i];
					px_k[i]=px[i];
					py_k[i]=py[i];
					pz_k[i]=pz[i];
					dTr_k[i]=dTr[i]; dTr_k[i+Nx]=dTr[i+Nx]; dTr_k[i+2*Nx]=dTr[i+2*Nx];	dTr_k[i+3*Nx]=dTr[i+3*Nx];	dTr_k[i+4*Nx]=dTr[i+4*Nx];	dTr_k[i+5*Nx]=dTr[i+5*Nx];	
					
					s[i]=0;	s[i+Nx]=0; s[i+2*Nx]=0; s[i+3*Nx]=0; s[i+4*Nx]=0; s[i+5*Nx]=0;
					y[i]=0;	y[i+Nx]=0; y[i+2*Nx]=0;	y[i+3*Nx]=0; y[i+4*Nx]=0; y[i+5*Nx]=0;
					sB[i]=0; sB[i+Nx]=0; sB[i+2*Nx]=0; sB[i+3*Nx]=0; sB[i+4*Nx]=0; sB[i+5*Nx]=0;
					Bs[i]=0; Bs[i+Nx]=0; Bs[i+2*Nx]=0; Bs[i+3*Nx]=0; Bs[i+4*Nx]=0; Bs[i+5*Nx]=0;

					for(int j=0;j<Nx*6;j++)
					{
						BssB[i*(Nx*6)+j]=0;
						BssB[(i+Nx)*(Nx*6)+j+Nx]=0;
						BssB[(i+2*Nx)*(Nx*6)+j+2*Nx]=0;
						BssB[(i+3*Nx)*(Nx*6)+j+3*Nx]=0;
						BssB[(i+4*Nx)*(Nx*6)+j+4*Nx]=0;
						BssB[(i+5*Nx)*(Nx*6)+j+5*Nx]=0;
					}
				}

				double Tr_min=Tr;
				double a_min=1e-3;
				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;

					double Tr_a=0;
					for(int i=0;i<Nx;i++)
					{
						qx_a[i]=qx[i]+d[i]*alpha;
						qy_a[i]=qy[i]+d[i+Nx]*alpha;
						qz_a[i]=qz[i]+d[i+2*Nx]*alpha;
						px_a[i]=px[i]+d[i+3*Nx]*alpha;
						py_a[i]=py[i]+d[i+4*Nx]*alpha;
						pz_a[i]=pz[i]+d[i+5*Nx]*alpha;

						Tr_a+=0.5/mi*(px_a[i]*px_a[i]+py_a[i]*py_a[i]+pz_a[i]*pz_a[i])-V*w_a[i]+0.5*r*(g_a[i]+seta[i])*(g_a[i]+seta[i]);
			
						double hi_a=-V*( (qx_a[i]-ax)*nx + (qy_a[i]-ay)*ny +(qz_a[i]-az)*nz );
						if(hi_a+seta[i+Nc/2]>0)	Tr_a+=0.5*r*(hi_a+seta[i+Nc/2])*(hi_a+seta[i+Nc/2]);
					}

					if(Tr_a<Tr_min)
					{
						Tr_min=Tr_a;
						a_min=alpha;
					}
				}
				cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
				for(int i=0;i<Nx;i++)
				{
					qx[i]+=d[i]*=a_min;
					qy[i]+=d[i+Nx]*a_min;
					qz[i]+=d[i+2*Nx]*a_min;
					px[i]+=d[i+3*Nx]*a_min;
					py[i]+=d[i+4*Nx]*a_min;
					pz[i]+=d[i+5*Nx]*a_min;
				}

				Tr=0;
				for(int i=0;i<Nx;i++)
				{
					Tr+=0.5/mi*(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-V*w[i]+0.5*r*(g[i]+seta[i])*(g[i]+seta[i]);			

					dTr[i]=0;
					dTr[i+Nx]=0;
					dTr[i+2*Nx]=0;
					for(int j=0;j<Nx;j++)
					{
						dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_X]*(-g[j]+seta[j]);
						dTr[i+Nx]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_Y]*(-g[j]+seta[j]);
						dTr[i+2*Nx]+=-V*(HYPER[j].stress[A_Z][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Z][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_Z]*(-g[j]+seta[j]);			
					}

					h[i]=-V*( (qx[i]-ax)*nx + (qy[i]-ay)*ny +(qz[i]-az)*nz );	
					if(h[i]+seta[i+Nc/2]>0)
					{
						Tr+=0.5*r*(h[i]+seta[i+Nc/2])*(h[i]+seta[i+Nc/2]);
						dTr[i]+=r*(-V*nx+seta[i+Nc/2]);
						dTr[i+Nx]+=r*(-V*ny+seta[i+Nc/2]);
						dTr[i+2*Nx]+=r*(-V*nz+seta[i+Nc/2]);
					}

					dTr[i+Nx*3]=1/mi*px[i];
					dTr[i+Nx*4]=1/mi*py[i];
					dTr[i+Nx*5]=1/mi*pz[i];

					E+=sqrt(dTr[i]*dTr[i] + dTr[i+Nx]*dTr[i+Nx] + dTr[i+Nx*2]*dTr[i+Nx*2] + dTr[i+Nx*3]*dTr[i+Nx*3] + dTr[i+Nx*4]*dTr[i+Nx*4] + dTr[i+Nx*5]*dTr[i+Nx*5]);
				}
				if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
				if(E<ep)	break;

				double beta=0;
				double sigma=0;

				for(int i=0;i<Nx;i++)
				{
					s[i]=qx[i]-qx_k[i];	s[i+Nx]=qy[i]-qy_k[i]; s[i+2*Nx]=qz[i]-qz_k[i]; s[i+3*Nx]=px[i]-px_k[i]; s[i+4*Nx]=py[i]-py_k[i]; s[i+5*Nx]=pz[i]-pz_k[i];
	
					y[i]=dTr[i]-dTr_k[i]; y[i+Nx]=dTr[i+Nx]-dTr_k[i+Nx]; y[i+2*Nx]=dTr[i+2*Nx]-dTr_k[i+2*Nx]; y[i+3*Nx]=dTr[i+3*Nx]-dTr_k[i+3*Nx]; y[i+4*Nx]=dTr[i+4*Nx]-dTr_k[i+4*Nx]; y[i+5*Nx]=dTr[i+5*Nx]-dTr_k[i+5*Nx];
					beta+=s[i]*y[i]+s[i+Nx]*y[i+Nx]+s[i+2*Nx]*y[i+2*Nx]+s[i+3*Nx]*y[i+3*Nx]+s[i+4*Nx]*y[i+4*Nx]+s[i+5*Nx]*y[i+5*Nx];
	
					sB[i]=0; sB[i+Nx]=0; sB[i+2*Nx]=0; sB[i+3*Nx]=0; sB[i+4*Nx]=0; sB[i+5*Nx]=0;
					for(int j=0;j<Nx;j++)
					{
						sB[i]+=s[j]*B[j*Nx*6+i]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i];
						sB[i+Nx]+=s[j]*B[j*Nx*6+i+Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+Nx];
						sB[i+2*Nx]+=s[j]*B[j*Nx*6+i+2*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+2*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+2*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+2*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+2*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+2*Nx];
						sB[i+3*Nx]+=s[j]*B[j*Nx*6+i+3*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+3*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+3*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+3*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+3*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+3*Nx];
						sB[i+4*Nx]+=s[j]*B[j*Nx*6+i+4*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+4*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+4*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+4*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+4*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+4*Nx];
						sB[i+5*Nx]+=s[j]*B[j*Nx*6+i+5*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+5*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+5*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+5*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+5*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+5*Nx];
					}
					sigma+=sB[i]*s[i]+sB[i+Nx]*s[i+Nx]+sB[i+2*Nx]*s[i+2*Nx]+sB[i+3*Nx]*s[i+3*Nx]+sB[i+4*Nx]*s[i+4*Nx]+sB[i+5*Nx]*s[i+5*Nx];
				}
				cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

				if(beta>0)//(beta>=0.2*sigma)
				{
					for(int i=0;i<Nx;i++)
					{
						Bs[i]=0; Bs[i+Nx]=0; Bs[i+2*Nx]=0; Bs[i+3*Nx]=0; Bs[i+4*Nx]=0; Bs[i+5*Nx]=0;
						for(int j=0;j<Nx;j++)
						{
							Bs[i]+=B[i*Nx*6+j]*s[j]+B[i*6*Nx+j+Nx]*s[j+Nx]+B[i*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+Nx]+=B[i+Nx*Nx*6+j]*s[j]+B[i+Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+2*Nx]+=B[i+2*Nx*Nx*6+j]*s[j]+B[i+2*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+2*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+2*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+2*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+2*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+3*Nx]+=B[i+3*Nx*Nx*6+j]*s[j]+B[i+3*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+3*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+3*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+3*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+3*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+4*Nx]+=B[i+4*Nx*Nx*6+j]*s[j]+B[i+4*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+4*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+4*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+4*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+4*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+5*Nx]+=B[i+5*Nx*Nx*6+j]*s[j]+B[i+5*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+5*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+5*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+5*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+5*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
						}
					}
					for(int i=0;i<6*Nx;i++)
					{
						for(int j=0;j<6*Nx;j++)
						{
							BssB[i*Nx*6+j]=0;
							for(int k=0;k<6*Nx;k++)	BssB[i*Nx*6+j]+=Bs[i]*s[k]*B[k*Nx*6+j];
							B[i*Nx*6+j]+=1/beta*y[i]*y[j]-1/sigma*BssB[i*Nx*6+j];
						}
					}
				}
			}
			delete[]	qx_k;
			delete[]	qy_k;
			delete[]	qz_k;
			delete[]	px_k;
			delete[]	py_k;
			delete[]	pz_k;
			delete[]	dTr_k;

			delete[]	qx_a;
			delete[]	qy_a;
			delete[]	qz_a;
			delete[]	px_a;
			delete[]	py_a;
			delete[]	pz_a;
			delete[]	w_a;
			delete[]	g_a;

			delete[]	s;
			delete[]	y;
			delete[]	sB;
			delete[]	Bs;
			delete[]	BssB;		

		}
		for(int i=0;i<Nx;i++)
		{
			E_min+=sqrt((old_qx[i]-qx[i])*(old_qx[i]-qx[i]) + (old_qy[i]-qy[i])*(old_qy[i]-qy[i]) + (old_qz[i]-qz[i])*(old_qz[i]-qz[i])
				+ (old_px[i]-px[i])*(old_px[i]-px[i]) + (old_py[i]-py[i])*(old_py[i]-py[i]) + (old_pz[i]-pz[i])*(old_pz[i]-pz[i]));

			seta[i]+=g[i];
			seta[i+Nx]+=h[i];
		}
		if(E_min<ep*1000)	r*=4;
	}

	double *fc=new double [Nx];
	int *Nv_n=new int [Nx];

	for(int i=0;i<Nx;i++)
	{
		fc[i]=-g[i];
		fc[i+Nx]=h[i];
		HYPER[i].lambda=0;
		HYPER[i].mu=0;
	}

	int nv=Nx;
	for(int i=0;i<Nx;i++)
	{
		if(fc[i+Nx]>ep)	HYPER[i].mu=0;
		else
		{
			Nv_n[nv]=i;
			nv++;
		}
	}

	double *Nrv=new double [nv];
	double *Nlv=new double [nv*nv];

	for(int i=0;i<Nx;i++)
	{
		dTr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z]);
		}
		Nrv[i]=-dTr[i];
	}

	for(int i=Nx;i<nv;i++)
	{
		dTr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dTr[i]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z]);
		}
		Nrv[i]=-dTr[i];
	}
	for(int i=0;i<nv;i++)
	{

	}


	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<nv;j++)
		{
			if(j<Nx)		Nlv[i*nv+j]=HYPER1[j*Nx+i].DgDq[A_X];
			else
			{
				int jv=Nv_n[j];
				Nlv[i*nv+jv]=-V*nx;
			}
			int jv=Nv_n[j];

			if(jv==0)	Nlv[i*nv+j]=dc0[i];
			else if(jv==1)	Nlv[i*nv+j]=dc1[i];
			else if(jv==2)	Nlv[i*nv+j]=dc2[i];
			cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
		}
		//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
		Nrv[i]=-dTr[i];
		cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
	}


	delete[]	fc;
	delete[]	Nv_n;


	int *Nv_n=new int [3];
	Nv_n[0]=0;	Nv_n[1]=0;	Nv_n[2]=0;

	//v0Å@åvéZ
	int nv=0;
	cout<<"fc="<<c0<<", "<<c1<<", "<<c2<<endl;
	if(c0>ep)	v[0]=0;
	else
	{
		Nv_n[nv]=0;
		nv++;
	}
	if(c1>ep)	v[1]=0;
	else
	{
		Nv_n[nv]=1;
		nv++;
	}
	if(c2>ep)	v[2]=0;
	else
	{
		Nv_n[nv]=2;
		nv++;
	}
	cout<<"nv="<<nv<<endl;
	if(nv>1)
	{
		double  *Nrv=new double [nv];
		double  *Nlv=new double [nv*nv];

		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nlv[i*nv+j]=dc0[i];
				else if(jv==1)	Nlv[i*nv+j]=dc1[i];
				else if(jv==2)	Nlv[i*nv+j]=dc2[i];
				cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
			}
			//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
			Nrv[i]=-dfx[i];
			cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
		}
		gauss(Nlv,Nrv,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nrv[i];
			else if(iv==1)	v[1]=Nrv[i];
			else if(iv==2)	v[2]=Nrv[i];
		}
		delete[]	Nrv;
		delete[]	Nlv;
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	}
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;

	delete[]	v;
	delete[]	Nv_n;


	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"óùò_â"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;







	delete[]	qx;
	delete[]	qy;
	delete[]	qz;
	delete[]	old_qx;
	delete[]	old_qy;
	delete[]	old_qz;
	delete[]	px;
	delete[]	py;
	delete[]	pz;
	delete[]	old_px;
	delete[]	old_py;
	delete[]	old_pz;
	delete[]	dTr;
	delete[]	w;
	delete[]	g;
	delete[]	h;
	delete[]	d;
	delete[]	B;
	delete[]	Nr;
	delete[]	seta;

}





void QP(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	int Nx=HYPER.size();
	double *qx=new double [Nx];
	double *qy=new double [Nx];
	double *qz=new double [Nx];
	double *old_qx=new double [Nx];
	double *old_qy=new double [Nx];
	double *old_qz=new double [Nx];

	double *px=new double [Nx];
	double *py=new double [Nx];
	double *pz=new double [Nx];
	double *old_px=new double [Nx];
	double *old_py=new double [Nx];
	double *old_pz=new double [Nx];


	double Tr=0;
	double *dTr=new double [Nx*6];
	double *w=new double [Nx];
	double *g=new double [Nx];
	double *h=new double [Nx];

	double V=get_volume(&CON);
	double mi=CON.get_h_dis()*V;

	double ax=0;
	double ay=0;
	double az=0;
	double nx=0;
	double ny=0;
	double nz=1;

	double *B=new double[36*Nx*Nx];
	double *d=new double [Nx*6];	
	double *Nr=new double [Nx*6];

	for(int i=0; i<Nx; i++)
	{
		qx[i]=PART[i].r[A_X];
		qy[i]=PART[i].r[A_Y];
		qz[i]=PART[i].r[A_Z];
		px[i]=HYPER[i].p[A_X];
		py[i]=HYPER[i].p[A_Y];
		pz[i]=HYPER[i].p[A_Z];
		dTr[i]=0;
		w[i]=0;
		g[i]=0;
		h[i]=0;
	}
	for(int i=Nx; i<Nx*6; i++)	dTr[i]=0;

	int Nc=2*HYPER.size();	//çSë©èåèÇÃå¬êî

	double *seta=new double [Nc];


	for(int i=0;i<Nc;i++)			seta[i]=0;

	double r=0;
	double ep=1e-10;

	double E_min=1;
	int count_min=0;

	while(E_min>ep)
	{
		count_min++;

		for(int i=0;i<Nx;i++)
		{
			old_qx[i]=qx[i];
			old_qy[i]=qy[i];
			old_qz[i]=qz[i];
			old_px[i]=px[i];
			old_py[i]=py[i];
			old_pz[i]=pz[i];
		}
		for(int i=0;i<Nx*6;i++)
		{
			d[i]=0;
			Nr[i]=0;
			for(int j=0;j<Nx*6;j++)
			{
				if(j==i)	B[i*Nx+j]=1;
				else
				{
					B[i*Nx+j]=0;
				}
			}
		}

		double E=1;
		int count=0;

		calc_wg(CON,HYPER,HYPER1,qx,qy,qz,w,g);
		Tr=0;
		for(int i=0;i<Nx;i++)
		{
			Tr+=0.5/mi*(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-V*w[i]+0.5*r*(-g[i]+seta[i])*(-g[i]+seta[i]);			

			dTr[i]=0;
			dTr[i+Nx]=0;
			dTr[i+2*Nx]=0;
			for(int j=0;j<Nx;j++)
			{
				dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_X]*(-g[j]+seta[j]);

				dTr[i+Nx]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_Y]*(-g[j]+seta[j]);
				dTr[i+2*Nx]+=-V*(HYPER[j].stress[A_Z][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Z][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_Z]*(-g[j]+seta[j]);			
			}

			h[i]=-V*( (qx[i]-ax)*nx + (qy[i]-ay)*ny +(qz[i]-az)*nz );	
			if(h[i]+seta[i+Nc/2]>0)
			{
				Tr+=0.5*r*(h[i]+seta[i+Nc/2])*(h[i]+seta[i+Nc/2]);
				dTr[i]+=r*(-V*nx+seta[i+Nc/2]);
				dTr[i+Nx]+=r*(-V*ny+seta[i+Nc/2]);
				dTr[i+2*Nx]+=r*(-V*nz+seta[i+Nc/2]);
			}

			dTr[i+Nx*3]=1/mi*px[i];
			dTr[i+Nx*4]=1/mi*py[i];
			dTr[i+Nx*5]=1/mi*pz[i];

			E+=sqrt(dTr[i]*dTr[i] + dTr[i+Nx]*dTr[i+Nx] + dTr[i+Nx*2]*dTr[i+Nx*2] + dTr[i+Nx*3]*dTr[i+Nx*3] + dTr[i+Nx*4]*dTr[i+Nx*4] + dTr[i+Nx*5]*dTr[i+Nx*5]);
		}
		cout<<"E0="<<E<<endl;


		if(E<ep)	break;
		else
		{
			double *qx_k=new double [Nx];
			double *qy_k=new double [Nx];
			double *qz_k=new double [Nx];
			double *px_k=new double [Nx];
			double *py_k=new double [Nx];
			double *pz_k=new double [Nx];
			double *dTr_k=new double [Nx*6];

			double *qx_a=new double [Nx];
			double *qy_a=new double [Nx];
			double *qz_a=new double [Nx];
			double *px_a=new double [Nx];
			double *py_a=new double [Nx];
			double *pz_a=new double [Nx];
			double *w_a=new double [Nx];
			double *g_a=new double [Nx];

			double *s=new double [Nx*6];
			double *y=new double [Nx*6];
			double *sB=new double [Nx*6];
			double *Bs=new double [Nx*6];
			double *BssB=new double [Nx*Nx*36];


			while(E>ep)
			{
				count++;

				for(int i=0;i<Nx;i++)
				{
					qx_k[i]=qx[i];
					qy_k[i]=qy[i];
					qz_k[i]=qz[i];
					px_k[i]=px[i];
					py_k[i]=py[i];
					pz_k[i]=pz[i];
					dTr_k[i]=dTr[i]; dTr_k[i+Nx]=dTr[i+Nx]; dTr_k[i+2*Nx]=dTr[i+2*Nx];	dTr_k[i+3*Nx]=dTr[i+3*Nx];	dTr_k[i+4*Nx]=dTr[i+4*Nx];	dTr_k[i+5*Nx]=dTr[i+5*Nx];	
					
					s[i]=0;	s[i+Nx]=0; s[i+2*Nx]=0; s[i+3*Nx]=0; s[i+4*Nx]=0; s[i+5*Nx]=0;
					y[i]=0;	y[i+Nx]=0; y[i+2*Nx]=0;	y[i+3*Nx]=0; y[i+4*Nx]=0; y[i+5*Nx]=0;
					sB[i]=0; sB[i+Nx]=0; sB[i+2*Nx]=0; sB[i+3*Nx]=0; sB[i+4*Nx]=0; sB[i+5*Nx]=0;
					Bs[i]=0; Bs[i+Nx]=0; Bs[i+2*Nx]=0; Bs[i+3*Nx]=0; Bs[i+4*Nx]=0; Bs[i+5*Nx]=0;

					for(int j=0;j<Nx*6;j++)
					{
						BssB[i*(Nx*6)+j]=0;
						BssB[(i+Nx)*(Nx*6)+j+Nx]=0;
						BssB[(i+2*Nx)*(Nx*6)+j+2*Nx]=0;
						BssB[(i+3*Nx)*(Nx*6)+j+3*Nx]=0;
						BssB[(i+4*Nx)*(Nx*6)+j+4*Nx]=0;
						BssB[(i+5*Nx)*(Nx*6)+j+5*Nx]=0;
					}
				}

				double Tr_min=Tr;
				double a_min=1e-3;
				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;

					double Tr_a=0;
					for(int i=0;i<Nx;i++)
					{
						qx_a[i]=qx[i]+d[i]*alpha;
						qy_a[i]=qy[i]+d[i+Nx]*alpha;
						qz_a[i]=qz[i]+d[i+2*Nx]*alpha;
						px_a[i]=px[i]+d[i+3*Nx]*alpha;
						py_a[i]=py[i]+d[i+4*Nx]*alpha;
						pz_a[i]=pz[i]+d[i+5*Nx]*alpha;

						Tr_a+=0.5/mi*(px_a[i]*px_a[i]+py_a[i]*py_a[i]+pz_a[i]*pz_a[i])-V*w_a[i]+0.5*r*(g_a[i]+seta[i])*(g_a[i]+seta[i]);
			
						double hi_a=-V*( (qx_a[i]-ax)*nx + (qy_a[i]-ay)*ny +(qz_a[i]-az)*nz );
						if(hi_a+seta[i+Nc/2]>0)	Tr_a+=0.5*r*(hi_a+seta[i+Nc/2])*(hi_a+seta[i+Nc/2]);
					}

					if(Tr_a<Tr_min)
					{
						Tr_min=Tr_a;
						a_min=alpha;
					}
				}
				cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
				for(int i=0;i<Nx;i++)
				{
					qx[i]+=d[i]*=a_min;
					qy[i]+=d[i+Nx]*a_min;
					qz[i]+=d[i+2*Nx]*a_min;
					px[i]+=d[i+3*Nx]*a_min;
					py[i]+=d[i+4*Nx]*a_min;
					pz[i]+=d[i+5*Nx]*a_min;
				}

				Tr=0;
				for(int i=0;i<Nx;i++)
				{
					Tr+=0.5/mi*(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-V*w[i]+0.5*r*(g[i]+seta[i])*(g[i]+seta[i]);			

					dTr[i]=0;
					dTr[i+Nx]=0;
					dTr[i+2*Nx]=0;
					for(int j=0;j<Nx;j++)
					{
						dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_X]*(-g[j]+seta[j]);
						dTr[i+Nx]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_Y]*(-g[j]+seta[j]);
						dTr[i+2*Nx]+=-V*(HYPER[j].stress[A_Z][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Z][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_Z]*(-g[j]+seta[j]);			
					}

					h[i]=-V*( (qx[i]-ax)*nx + (qy[i]-ay)*ny +(qz[i]-az)*nz );	
					if(h[i]+seta[i+Nc/2]>0)
					{
						Tr+=0.5*r*(h[i]+seta[i+Nc/2])*(h[i]+seta[i+Nc/2]);
						dTr[i]+=r*(-V*nx+seta[i+Nc/2]);
						dTr[i+Nx]+=r*(-V*ny+seta[i+Nc/2]);
						dTr[i+2*Nx]+=r*(-V*nz+seta[i+Nc/2]);
					}

					dTr[i+Nx*3]=1/mi*px[i];
					dTr[i+Nx*4]=1/mi*py[i];
					dTr[i+Nx*5]=1/mi*pz[i];

					E+=sqrt(dTr[i]*dTr[i] + dTr[i+Nx]*dTr[i+Nx] + dTr[i+Nx*2]*dTr[i+Nx*2] + dTr[i+Nx*3]*dTr[i+Nx*3] + dTr[i+Nx*4]*dTr[i+Nx*4] + dTr[i+Nx*5]*dTr[i+Nx*5]);
				}
				if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
				if(E<ep)	break;

				double beta=0;
				double sigma=0;

				for(int i=0;i<Nx;i++)
				{
					s[i]=qx[i]-qx_k[i];	s[i+Nx]=qy[i]-qy_k[i]; s[i+2*Nx]=qz[i]-qz_k[i]; s[i+3*Nx]=px[i]-px_k[i]; s[i+4*Nx]=py[i]-py_k[i]; s[i+5*Nx]=pz[i]-pz_k[i];
	
					y[i]=dTr[i]-dTr_k[i]; y[i+Nx]=dTr[i+Nx]-dTr_k[i+Nx]; y[i+2*Nx]=dTr[i+2*Nx]-dTr_k[i+2*Nx]; y[i+3*Nx]=dTr[i+3*Nx]-dTr_k[i+3*Nx]; y[i+4*Nx]=dTr[i+4*Nx]-dTr_k[i+4*Nx]; y[i+5*Nx]=dTr[i+5*Nx]-dTr_k[i+5*Nx];
					beta+=s[i]*y[i]+s[i+Nx]*y[i+Nx]+s[i+2*Nx]*y[i+2*Nx]+s[i+3*Nx]*y[i+3*Nx]+s[i+4*Nx]*y[i+4*Nx]+s[i+5*Nx]*y[i+5*Nx];
	
					sB[i]=0; sB[i+Nx]=0; sB[i+2*Nx]=0; sB[i+3*Nx]=0; sB[i+4*Nx]=0; sB[i+5*Nx]=0;
					for(int j=0;j<Nx;j++)
					{
						sB[i]+=s[j]*B[j*Nx*6+i]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i];
						sB[i+Nx]+=s[j]*B[j*Nx*6+i+Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+Nx];
						sB[i+2*Nx]+=s[j]*B[j*Nx*6+i+2*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+2*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+2*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+2*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+2*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+2*Nx];
						sB[i+3*Nx]+=s[j]*B[j*Nx*6+i+3*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+3*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+3*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+3*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+3*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+3*Nx];
						sB[i+4*Nx]+=s[j]*B[j*Nx*6+i+4*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+4*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+4*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+4*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+4*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+4*Nx];
						sB[i+5*Nx]+=s[j]*B[j*Nx*6+i+5*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+5*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+5*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+5*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+5*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+5*Nx];
					}
					sigma+=sB[i]*s[i]+sB[i+Nx]*s[i+Nx]+sB[i+2*Nx]*s[i+2*Nx]+sB[i+3*Nx]*s[i+3*Nx]+sB[i+4*Nx]*s[i+4*Nx]+sB[i+5*Nx]*s[i+5*Nx];
				}
				cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

				if(beta>0)//(beta>=0.2*sigma)
				{
					for(int i=0;i<Nx;i++)
					{
						Bs[i]=0; Bs[i+Nx]=0; Bs[i+2*Nx]=0; Bs[i+3*Nx]=0; Bs[i+4*Nx]=0; Bs[i+5*Nx]=0;
						for(int j=0;j<Nx;j++)
						{
							Bs[i]+=B[i*Nx*6+j]*s[j]+B[i*6*Nx+j+Nx]*s[j+Nx]+B[i*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+Nx]+=B[i+Nx*Nx*6+j]*s[j]+B[i+Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+2*Nx]+=B[i+2*Nx*Nx*6+j]*s[j]+B[i+2*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+2*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+2*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+2*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+2*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+3*Nx]+=B[i+3*Nx*Nx*6+j]*s[j]+B[i+3*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+3*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+3*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+3*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+3*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+4*Nx]+=B[i+4*Nx*Nx*6+j]*s[j]+B[i+4*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+4*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+4*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+4*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+4*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+5*Nx]+=B[i+5*Nx*Nx*6+j]*s[j]+B[i+5*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+5*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+5*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+5*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+5*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
						}
					}
					for(int i=0;i<6*Nx;i++)
					{
						for(int j=0;j<6*Nx;j++)
						{
							BssB[i*Nx*6+j]=0;
							for(int k=0;k<6*Nx;k++)	BssB[i*Nx*6+j]+=Bs[i]*s[k]*B[k*Nx*6+j];
							B[i*Nx*6+j]+=1/beta*y[i]*y[j]-1/sigma*BssB[i*Nx*6+j];
						}
					}
				}
			}
			delete[]	qx_k;
			delete[]	qy_k;
			delete[]	qz_k;
			delete[]	px_k;
			delete[]	py_k;
			delete[]	pz_k;
			delete[]	dTr_k;

			delete[]	qx_a;
			delete[]	qy_a;
			delete[]	qz_a;
			delete[]	px_a;
			delete[]	py_a;
			delete[]	pz_a;
			delete[]	w_a;
			delete[]	g_a;

			delete[]	s;
			delete[]	y;
			delete[]	sB;
			delete[]	Bs;
			delete[]	BssB;		

		}
		for(int i=0;i<Nx;i++)
		{
			E_min+=sqrt((old_qx[i]-qx[i])*(old_qx[i]-qx[i]) + (old_qy[i]-qy[i])*(old_qy[i]-qy[i]) + (old_qz[i]-qz[i])*(old_qz[i]-qz[i])
				+ (old_px[i]-px[i])*(old_px[i]-px[i]) + (old_py[i]-py[i])*(old_py[i]-py[i]) + (old_pz[i]-pz[i])*(old_pz[i]-pz[i]));

			seta[i]+=g[i];
			seta[i+Nx]+=h[i];
		}
		if(E_min<ep*1000)	r*=4;
	}

	double *fc=new double [Nx];
	int *Nv_n=new int [Nx];

	for(int i=0;i<Nx;i++)
	{
		fc[i]=-g[i];
		fc[i+Nx]=h[i];
		HYPER[i].lambda=0;
		HYPER[i].mu=0;
	}

	int nv=Nx;
	for(int i=0;i<Nx;i++)
	{
		if(fc[i+Nx]>ep)	HYPER[i].mu=0;
		else
		{
			Nv_n[nv]=i;
			nv++;
		}
	}

	double *Nrv=new double [nv];
	double *Nlv=new double [nv*nv];

	for(int i=0;i<Nx;i++)
	{
		dTr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z]);
		}
		Nrv[i]=-dTr[i];
	}

	for(int i=Nx;i<nv;i++)
	{
		dTr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dTr[i]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z]);
		}
		Nrv[i]=-dTr[i];
	}
	for(int i=0;i<nv;i++)
	{

	}


	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<nv;j++)
		{
			if(j<Nx)		Nlv[i*nv+j]=HYPER1[j*Nx+i].DgDq[A_X];
			else
			{
				int jv=Nv_n[j];
				Nlv[i*nv+jv]=-V*nx;
			}
			int jv=Nv_n[j];

			if(jv==0)	Nlv[i*nv+j]=dc0[i];
			else if(jv==1)	Nlv[i*nv+j]=dc1[i];
			else if(jv==2)	Nlv[i*nv+j]=dc2[i];
			cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
		}
		//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
		Nrv[i]=-dTr[i];
		cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
	}


	delete[]	fc;
	delete[]	Nv_n;


	int *Nv_n=new int [3];
	Nv_n[0]=0;	Nv_n[1]=0;	Nv_n[2]=0;

	//v0Å@åvéZ
	int nv=0;
	cout<<"fc="<<c0<<", "<<c1<<", "<<c2<<endl;
	if(c0>ep)	v[0]=0;
	else
	{
		Nv_n[nv]=0;
		nv++;
	}
	if(c1>ep)	v[1]=0;
	else
	{
		Nv_n[nv]=1;
		nv++;
	}
	if(c2>ep)	v[2]=0;
	else
	{
		Nv_n[nv]=2;
		nv++;
	}
	cout<<"nv="<<nv<<endl;
	if(nv>1)
	{
		double  *Nrv=new double [nv];
		double  *Nlv=new double [nv*nv];

		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nlv[i*nv+j]=dc0[i];
				else if(jv==1)	Nlv[i*nv+j]=dc1[i];
				else if(jv==2)	Nlv[i*nv+j]=dc2[i];
				cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
			}
			//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
			Nrv[i]=-dfx[i];
			cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
		}
		gauss(Nlv,Nrv,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nrv[i];
			else if(iv==1)	v[1]=Nrv[i];
			else if(iv==2)	v[2]=Nrv[i];
		}
		delete[]	Nrv;
		delete[]	Nlv;
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	}
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;

	delete[]	v;
	delete[]	Nv_n;


	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"óùò_â"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;







	delete[]	qx;
	delete[]	qy;
	delete[]	qz;
	delete[]	old_qx;
	delete[]	old_qy;
	delete[]	old_qz;
	delete[]	px;
	delete[]	py;
	delete[]	pz;
	delete[]	old_px;
	delete[]	old_py;
	delete[]	old_pz;
	delete[]	dTr;
	delete[]	w;
	delete[]	g;
	delete[]	h;
	delete[]	d;
	delete[]	B;
	delete[]	Nr;
	delete[]	seta;

}


void QP2(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	int Nx=HYPER.size();
	double *qx=new double [Nx];
	double *qy=new double [Nx];
	double *qz=new double [Nx];
	double *old_qx=new double [Nx];
	double *old_qy=new double [Nx];
	double *old_qz=new double [Nx];

	double *px=new double [Nx];
	double *py=new double [Nx];
	double *pz=new double [Nx];
	double *old_px=new double [Nx];
	double *old_py=new double [Nx];
	double *old_pz=new double [Nx];


	double Tr=0;
	double *dTr=new double [Nx*6];
	double *w=new double [Nx];
	double *g=new double [Nx];
	double *h=new double [Nx];

	double V=get_volume(&CON);
	double mi=CON.get_h_dis()*V;

	double ax=0;
	double ay=0;
	double az=0;
	double nx=0;
	double ny=0;
	double nz=1;

	double *B=new double[36*Nx*Nx];
	double *d=new double [Nx*6];	
	double *Nr=new double [Nx*6];

	for(int i=0; i<Nx; i++)
	{
		qx[i]=PART[i].r[A_X];
		qy[i]=PART[i].r[A_Y];
		qz[i]=PART[i].r[A_Z];
		px[i]=HYPER[i].p[A_X];
		py[i]=HYPER[i].p[A_Y];
		pz[i]=HYPER[i].p[A_Z];
		dTr[i]=0;
		w[i]=0;
		g[i]=0;
		h[i]=0;
	}
	for(int i=Nx; i<Nx*6; i++)	dTr[i]=0;

	int Nc=2*HYPER.size();	//çSë©èåèÇÃå¬êî

	double *seta=new double [Nc];


	for(int i=0;i<Nc;i++)			seta[i]=0;

	double r=0;
	double ep=1e-10;

	double E_min=1;
	int count_min=0;

	while(E_min>ep)
	{
		count_min++;

		for(int i=0;i<Nx;i++)
		{
			old_qx[i]=qx[i];
			old_qy[i]=qy[i];
			old_qz[i]=qz[i];
			old_px[i]=px[i];
			old_py[i]=py[i];
			old_pz[i]=pz[i];
		}
		for(int i=0;i<Nx*6;i++)
		{
			d[i]=0;
			Nr[i]=0;
			for(int j=0;j<Nx*6;j++)
			{
				if(j==i)	B[i*Nx+j]=1;
				else
				{
					B[i*Nx+j]=0;
				}
			}
		}

		double E=1;
		int count=0;

		calc_wg(CON,HYPER,HYPER1,qx,qy,qz,w,g);
		Tr=0;
		for(int i=0;i<Nx;i++)
		{
			Tr+=0.5/mi*(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-V*w[i]+0.5*r*(-g[i]+seta[i])*(-g[i]+seta[i]);			

			dTr[i]=0;
			dTr[i+Nx]=0;
			dTr[i+2*Nx]=0;
			for(int j=0;j<Nx;j++)
			{
				dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_X]*(-g[j]+seta[j]);

				dTr[i+Nx]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_Y]*(-g[j]+seta[j]);
				dTr[i+2*Nx]+=-V*(HYPER[j].stress[A_Z][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Z][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
					-r*HYPER1[j*Nx+i].DgDq[A_Z]*(-g[j]+seta[j]);			
			}

			h[i]=-V*( (qx[i]-ax)*nx + (qy[i]-ay)*ny +(qz[i]-az)*nz );	
			if(h[i]+seta[i+Nc/2]>0)
			{
				Tr+=0.5*r*(h[i]+seta[i+Nc/2])*(h[i]+seta[i+Nc/2]);
				dTr[i]+=r*(-V*nx+seta[i+Nc/2]);
				dTr[i+Nx]+=r*(-V*ny+seta[i+Nc/2]);
				dTr[i+2*Nx]+=r*(-V*nz+seta[i+Nc/2]);
			}

			dTr[i+Nx*3]=1/mi*px[i];
			dTr[i+Nx*4]=1/mi*py[i];
			dTr[i+Nx*5]=1/mi*pz[i];

			E+=sqrt(dTr[i]*dTr[i] + dTr[i+Nx]*dTr[i+Nx] + dTr[i+Nx*2]*dTr[i+Nx*2] + dTr[i+Nx*3]*dTr[i+Nx*3] + dTr[i+Nx*4]*dTr[i+Nx*4] + dTr[i+Nx*5]*dTr[i+Nx*5]);
		}
		cout<<"E0="<<E<<endl;


		if(E<ep)	break;
		else
		{
			double *qx_k=new double [Nx];
			double *qy_k=new double [Nx];
			double *qz_k=new double [Nx];
			double *px_k=new double [Nx];
			double *py_k=new double [Nx];
			double *pz_k=new double [Nx];
			double *dTr_k=new double [Nx*6];

			double *qx_a=new double [Nx];
			double *qy_a=new double [Nx];
			double *qz_a=new double [Nx];
			double *px_a=new double [Nx];
			double *py_a=new double [Nx];
			double *pz_a=new double [Nx];
			double *w_a=new double [Nx];
			double *g_a=new double [Nx];

			double *s=new double [Nx*6];
			double *y=new double [Nx*6];
			double *sB=new double [Nx*6];
			double *Bs=new double [Nx*6];
			double *BssB=new double [Nx*Nx*36];


			while(E>ep)
			{
				count++;

				for(int i=0;i<Nx;i++)
				{
					qx_k[i]=qx[i];
					qy_k[i]=qy[i];
					qz_k[i]=qz[i];
					px_k[i]=px[i];
					py_k[i]=py[i];
					pz_k[i]=pz[i];
					dTr_k[i]=dTr[i]; dTr_k[i+Nx]=dTr[i+Nx]; dTr_k[i+2*Nx]=dTr[i+2*Nx];	dTr_k[i+3*Nx]=dTr[i+3*Nx];	dTr_k[i+4*Nx]=dTr[i+4*Nx];	dTr_k[i+5*Nx]=dTr[i+5*Nx];	
					
					s[i]=0;	s[i+Nx]=0; s[i+2*Nx]=0; s[i+3*Nx]=0; s[i+4*Nx]=0; s[i+5*Nx]=0;
					y[i]=0;	y[i+Nx]=0; y[i+2*Nx]=0;	y[i+3*Nx]=0; y[i+4*Nx]=0; y[i+5*Nx]=0;
					sB[i]=0; sB[i+Nx]=0; sB[i+2*Nx]=0; sB[i+3*Nx]=0; sB[i+4*Nx]=0; sB[i+5*Nx]=0;
					Bs[i]=0; Bs[i+Nx]=0; Bs[i+2*Nx]=0; Bs[i+3*Nx]=0; Bs[i+4*Nx]=0; Bs[i+5*Nx]=0;

					for(int j=0;j<Nx*6;j++)
					{
						BssB[i*(Nx*6)+j]=0;
						BssB[(i+Nx)*(Nx*6)+j+Nx]=0;
						BssB[(i+2*Nx)*(Nx*6)+j+2*Nx]=0;
						BssB[(i+3*Nx)*(Nx*6)+j+3*Nx]=0;
						BssB[(i+4*Nx)*(Nx*6)+j+4*Nx]=0;
						BssB[(i+5*Nx)*(Nx*6)+j+5*Nx]=0;
					}
				}

				double Tr_min=Tr;
				double a_min=1e-3;
				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;

					double Tr_a=0;
					for(int i=0;i<Nx;i++)
					{
						qx_a[i]=qx[i]+d[i]*alpha;
						qy_a[i]=qy[i]+d[i+Nx]*alpha;
						qz_a[i]=qz[i]+d[i+2*Nx]*alpha;
						px_a[i]=px[i]+d[i+3*Nx]*alpha;
						py_a[i]=py[i]+d[i+4*Nx]*alpha;
						pz_a[i]=pz[i]+d[i+5*Nx]*alpha;

						Tr_a+=0.5/mi*(px_a[i]*px_a[i]+py_a[i]*py_a[i]+pz_a[i]*pz_a[i])-V*w_a[i]+0.5*r*(g_a[i]+seta[i])*(g_a[i]+seta[i]);
			
						double hi_a=-V*( (qx_a[i]-ax)*nx + (qy_a[i]-ay)*ny +(qz_a[i]-az)*nz );
						if(hi_a+seta[i+Nc/2]>0)	Tr_a+=0.5*r*(hi_a+seta[i+Nc/2])*(hi_a+seta[i+Nc/2]);
					}

					if(Tr_a<Tr_min)
					{
						Tr_min=Tr_a;
						a_min=alpha;
					}
				}
				cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
				for(int i=0;i<Nx;i++)
				{
					qx[i]+=d[i]*=a_min;
					qy[i]+=d[i+Nx]*a_min;
					qz[i]+=d[i+2*Nx]*a_min;
					px[i]+=d[i+3*Nx]*a_min;
					py[i]+=d[i+4*Nx]*a_min;
					pz[i]+=d[i+5*Nx]*a_min;
				}

				Tr=0;
				for(int i=0;i<Nx;i++)
				{
					Tr+=0.5/mi*(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-V*w[i]+0.5*r*(g[i]+seta[i])*(g[i]+seta[i]);			

					dTr[i]=0;
					dTr[i+Nx]=0;
					dTr[i+2*Nx]=0;
					for(int j=0;j<Nx;j++)
					{
						dTr[i]+=-V*(HYPER[j].stress[A_X][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_X][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_X][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_X]*(-g[j]+seta[j]);
						dTr[i+Nx]+=-V*(HYPER[j].stress[A_Y][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Y][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Y][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_Y]*(-g[j]+seta[j]);
						dTr[i+2*Nx]+=-V*(HYPER[j].stress[A_Z][A_X]*HYPER1[j*Nx+i].DgDq[A_X]+HYPER[j].stress[A_Z][A_Y]*HYPER1[j*Nx+i].DgDq[A_Y]+HYPER[j].stress[A_Z][A_Z]*HYPER1[j*Nx+i].DgDq[A_Z])
							-r*HYPER1[j*Nx+i].DgDq[A_Z]*(-g[j]+seta[j]);			
					}

					h[i]=-V*( (qx[i]-ax)*nx + (qy[i]-ay)*ny +(qz[i]-az)*nz );	
					if(h[i]+seta[i+Nc/2]>0)
					{
						Tr+=0.5*r*(h[i]+seta[i+Nc/2])*(h[i]+seta[i+Nc/2]);
						dTr[i]+=r*(-V*nx+seta[i+Nc/2]);
						dTr[i+Nx]+=r*(-V*ny+seta[i+Nc/2]);
						dTr[i+2*Nx]+=r*(-V*nz+seta[i+Nc/2]);
					}

					dTr[i+Nx*3]=1/mi*px[i];
					dTr[i+Nx*4]=1/mi*py[i];
					dTr[i+Nx*5]=1/mi*pz[i];

					E+=sqrt(dTr[i]*dTr[i] + dTr[i+Nx]*dTr[i+Nx] + dTr[i+Nx*2]*dTr[i+Nx*2] + dTr[i+Nx*3]*dTr[i+Nx*3] + dTr[i+Nx*4]*dTr[i+Nx*4] + dTr[i+Nx*5]*dTr[i+Nx*5]);
				}
				if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
				if(E<ep)	break;

				double beta=0;
				double sigma=0;

				for(int i=0;i<Nx;i++)
				{
					s[i]=qx[i]-qx_k[i];	s[i+Nx]=qy[i]-qy_k[i]; s[i+2*Nx]=qz[i]-qz_k[i]; s[i+3*Nx]=px[i]-px_k[i]; s[i+4*Nx]=py[i]-py_k[i]; s[i+5*Nx]=pz[i]-pz_k[i];
	
					y[i]=dTr[i]-dTr_k[i]; y[i+Nx]=dTr[i+Nx]-dTr_k[i+Nx]; y[i+2*Nx]=dTr[i+2*Nx]-dTr_k[i+2*Nx]; y[i+3*Nx]=dTr[i+3*Nx]-dTr_k[i+3*Nx]; y[i+4*Nx]=dTr[i+4*Nx]-dTr_k[i+4*Nx]; y[i+5*Nx]=dTr[i+5*Nx]-dTr_k[i+5*Nx];
					beta+=s[i]*y[i]+s[i+Nx]*y[i+Nx]+s[i+2*Nx]*y[i+2*Nx]+s[i+3*Nx]*y[i+3*Nx]+s[i+4*Nx]*y[i+4*Nx]+s[i+5*Nx]*y[i+5*Nx];
	
					sB[i]=0; sB[i+Nx]=0; sB[i+2*Nx]=0; sB[i+3*Nx]=0; sB[i+4*Nx]=0; sB[i+5*Nx]=0;
					for(int j=0;j<Nx;j++)
					{
						sB[i]+=s[j]*B[j*Nx*6+i]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i];
						sB[i+Nx]+=s[j]*B[j*Nx*6+i+Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+Nx];
						sB[i+2*Nx]+=s[j]*B[j*Nx*6+i+2*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+2*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+2*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+2*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+2*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+2*Nx];
						sB[i+3*Nx]+=s[j]*B[j*Nx*6+i+3*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+3*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+3*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+3*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+3*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+3*Nx];
						sB[i+4*Nx]+=s[j]*B[j*Nx*6+i+4*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+4*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+4*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+4*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+4*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+4*Nx];
						sB[i+5*Nx]+=s[j]*B[j*Nx*6+i+5*Nx]+s[(j+Nx)]*B[(j+Nx)*Nx*6+i+5*Nx]+s[(j+2*Nx)]*B[(j+2*Nx)*Nx*6+i+5*Nx]+s[(j+3*Nx)]*B[(j+3*Nx)*Nx*6+i+5*Nx]+s[(j+4*Nx)]*B[(j+4*Nx)*Nx*6+i+5*Nx]+s[(j+5*Nx)]*B[(j+5*Nx)*Nx*6+i+5*Nx];
					}
					sigma+=sB[i]*s[i]+sB[i+Nx]*s[i+Nx]+sB[i+2*Nx]*s[i+2*Nx]+sB[i+3*Nx]*s[i+3*Nx]+sB[i+4*Nx]*s[i+4*Nx]+sB[i+5*Nx]*s[i+5*Nx];
				}
				cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

				if(beta>0)//(beta>=0.2*sigma)
				{
					for(int i=0;i<Nx;i++)
					{
						Bs[i]=0; Bs[i+Nx]=0; Bs[i+2*Nx]=0; Bs[i+3*Nx]=0; Bs[i+4*Nx]=0; Bs[i+5*Nx]=0;
						for(int j=0;j<Nx;j++)
						{
							Bs[i]+=B[i*Nx*6+j]*s[j]+B[i*6*Nx+j+Nx]*s[j+Nx]+B[i*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+Nx]+=B[i+Nx*Nx*6+j]*s[j]+B[i+Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+2*Nx]+=B[i+2*Nx*Nx*6+j]*s[j]+B[i+2*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+2*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+2*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+2*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+2*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+3*Nx]+=B[i+3*Nx*Nx*6+j]*s[j]+B[i+3*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+3*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+3*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+3*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+3*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+4*Nx]+=B[i+4*Nx*Nx*6+j]*s[j]+B[i+4*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+4*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+4*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+4*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+4*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
							Bs[i+5*Nx]+=B[i+5*Nx*Nx*6+j]*s[j]+B[i+5*Nx*6*Nx+j+Nx]*s[j+Nx]+B[i+5*Nx*6*Nx+j+2*Nx]*s[j+2*Nx]+B[i+5*Nx*6*Nx+j+3*Nx]*s[j+3*Nx]+B[i+5*Nx*6*Nx+j+4*Nx]*s[j+4*Nx]+B[i+5*Nx*6*Nx+j+5*Nx]*s[j+5*Nx];
						}
					}
					for(int i=0;i<6*Nx;i++)
					{
						for(int j=0;j<6*Nx;j++)
						{
							BssB[i*Nx*6+j]=0;
							for(int k=0;k<6*Nx;k++)	BssB[i*Nx*6+j]+=Bs[i]*s[k]*B[k*Nx*6+j];
							B[i*Nx*6+j]+=1/beta*y[i]*y[j]-1/sigma*BssB[i*Nx*6+j];
						}
					}
				}
			}
			delete[]	qx_k;
			delete[]	qy_k;
			delete[]	qz_k;
			delete[]	px_k;
			delete[]	py_k;
			delete[]	pz_k;
			delete[]	dTr_k;

			delete[]	qx_a;
			delete[]	qy_a;
			delete[]	qz_a;
			delete[]	px_a;
			delete[]	py_a;
			delete[]	pz_a;
			delete[]	w_a;
			delete[]	g_a;

			delete[]	s;
			delete[]	y;
			delete[]	sB;
			delete[]	Bs;
			delete[]	BssB;		

		}
		for(int i=0;i<Nx;i++)
		{
			E_min+=sqrt((old_qx[i]-qx[i])*(old_qx[i]-qx[i]) + (old_qy[i]-qy[i])*(old_qy[i]-qy[i]) + (old_qz[i]-qz[i])*(old_qz[i]-qz[i])
				+ (old_px[i]-px[i])*(old_px[i]-px[i]) + (old_py[i]-py[i])*(old_py[i]-py[i]) + (old_pz[i]-pz[i])*(old_pz[i]-pz[i]));

			seta[i]+=g[i];
			seta[i+Nx]+=h[i];
		}
		if(E_min<ep*1000)	r*=4;
	}

	double *fc=new double [Nx];
	int *Nv_n=new int [Nx];

	for(int i=0;i<Nx;i++)
	{
		fc[i]=-g[i];
		fc[i+Nx]=h[i];
		HYPER[i].lambda=0;
		HYPER[i].mu=0;
	}

	int nv=Nx;
	for(int i=0;i<Nx;i++)
	{
		if(fc[i+Nx]>ep)	HYPER[i].mu=0;
		else
		{
			Nv_n[nv]=i;
			nv++;
		}
	}

	double *Nrv=new double [nv];
	double *Nlv=new double [nv*nv];

	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<nv;j++)
		{
			int jv=Nv_n[j];
			if(jv==0)	Nlv[i*nv+j]=dc0[i];
			else if(jv==1)	Nlv[i*nv+j]=dc1[i];
			else if(jv==2)	Nlv[i*nv+j]=dc2[i];
			cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
		}
		//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
		Nrv[i]=-dfx[i];
		cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
	}


	delete[]	fc;
	delete[]	Nv_n;


	int *Nv_n=new int [3];
	Nv_n[0]=0;	Nv_n[1]=0;	Nv_n[2]=0;

	//v0Å@åvéZ
	int nv=0;
	cout<<"fc="<<c0<<", "<<c1<<", "<<c2<<endl;
	if(c0>ep)	v[0]=0;
	else
	{
		Nv_n[nv]=0;
		nv++;
	}
	if(c1>ep)	v[1]=0;
	else
	{
		Nv_n[nv]=1;
		nv++;
	}
	if(c2>ep)	v[2]=0;
	else
	{
		Nv_n[nv]=2;
		nv++;
	}
	cout<<"nv="<<nv<<endl;
	if(nv>1)
	{
		double  *Nrv=new double [nv];
		double  *Nlv=new double [nv*nv];

		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nlv[i*nv+j]=dc0[i];
				else if(jv==1)	Nlv[i*nv+j]=dc1[i];
				else if(jv==2)	Nlv[i*nv+j]=dc2[i];
				cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
			}
			//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
			Nrv[i]=-dfx[i];
			cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
		}
		gauss(Nlv,Nrv,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nrv[i];
			else if(iv==1)	v[1]=Nrv[i];
			else if(iv==2)	v[2]=Nrv[i];
		}
		delete[]	Nrv;
		delete[]	Nlv;
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	}
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;

	delete[]	v;
	delete[]	Nv_n;


	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"óùò_â"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;







	delete[]	qx;
	delete[]	qy;
	delete[]	qz;
	delete[]	old_qx;
	delete[]	old_qy;
	delete[]	old_qz;
	delete[]	px;
	delete[]	py;
	delete[]	pz;
	delete[]	old_px;
	delete[]	old_py;
	delete[]	old_pz;
	delete[]	dTr;
	delete[]	w;
	delete[]	g;
	delete[]	h;
	delete[]	d;
	delete[]	B;
	delete[]	Nr;
	delete[]	seta;

}


void calc_wg(mpsconfig &CON, vector<mpselastic> PART, vector<hyperelastic> &HYPER, vector<hyperelastic2> &HYPER1, double **dq, double *w, double *g, double *J, double **F, double **ti_F, double **s)
{

//	cout<<"FiåvéZ";
	////FiÇÃçXêV
	int h_num=HYPER.size();
	double V=get_volume(&CON);
	double mi=V*CON.get_h_dis();
	
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];

	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		//FiÇÃåvéZ

		int Ni=HYPER[i].N;	
		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*h_num+inn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+inn].aiin[A_X],	HYPER1[i*h_num+inn].aiin[A_Y],	HYPER1[i*h_num+inn].aiin[A_Z]};
			
			fi[0][0]+=w*(PART[inn].r[A_X]+dq[inn][A_X] - PART[i].r[A_X]-dq[i][A_X])*a[A_X];
			fi[0][1]+=w*(PART[inn].r[A_X]+dq[inn][A_X] - PART[i].r[A_X]-dq[i][A_X])*a[A_Y];
			fi[0][2]+=w*(PART[inn].r[A_X]+dq[inn][A_X] - PART[i].r[A_X]-dq[i][A_X])*a[A_Z];
			
			fi[1][0]+=w*(PART[inn].r[A_Y]+dq[inn][A_Y] - PART[i].r[A_Y]-dq[i][A_Y])*a[A_X];
			fi[1][1]+=w*(PART[inn].r[A_Y]+dq[inn][A_Y] - PART[i].r[A_Y]-dq[i][A_Y])*a[A_Y];
			fi[1][2]+=w*(PART[inn].r[A_Y]+dq[inn][A_Y] - PART[i].r[A_Y]-dq[i][A_Y])*a[A_Z];

			fi[2][0]+=w*(PART[inn].r[A_Z]+dq[inn][A_Z] - PART[i].r[A_Z]-dq[i][A_Z])*a[A_X];
			fi[2][1]+=w*(PART[inn].r[A_Z]+dq[inn][A_Z] - PART[i].r[A_Z]-dq[i][A_Z])*a[A_Y];
			fi[2][2]+=w*(PART[inn].r[A_Z]+dq[inn][A_Z] - PART[i].r[A_Z]-dq[i][A_Z])*a[A_Z];
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
		
		F[i][0*DIMENSION+0]=p_Fi[0][0];	F[i][0*DIMENSION+1]=p_Fi[0][1];	F[i][0*DIMENSION+2]=p_Fi[0][2];	
		F[i][1*DIMENSION+0]=p_Fi[1][0];	F[i][1*DIMENSION+1]=p_Fi[1][1];	F[i][1*DIMENSION+2]=p_Fi[1][2];	
		F[i][2*DIMENSION+0]=p_Fi[2][0];	F[i][2*DIMENSION+1]=p_Fi[2][1];	F[i][2*DIMENSION+2]=p_Fi[2][2];	
	
		//JÇÃåvéZ
		J[i]=calc_det3(p_Fi);
	//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
		g[i]=V*(1-J[i]);
		//t_inverse_FiÇÃåvéZ
		inverse(p_Fi,DIMENSION);
		ti_F[i][0*DIMENSION+0]=p_Fi[0][0];	ti_F[i][0*DIMENSION+1]=p_Fi[1][0];	ti_F[i][0*DIMENSION+2]=p_Fi[2][0];
		ti_F[i][1*DIMENSION+0]=p_Fi[0][1];	ti_F[i][1*DIMENSION+1]=p_Fi[1][1];	ti_F[i][1*DIMENSION+2]=p_Fi[2][1];
		ti_F[i][2*DIMENSION+0]=p_Fi[0][2];	ti_F[i][2*DIMENSION+1]=p_Fi[1][2];	ti_F[i][2*DIMENSION+2]=p_Fi[2][2];
	}
	//	cout<<"----------OK"<<endl;

	for(int D=0;D>DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;




	double c10=CON.get_c10();
	double c01=CON.get_c01();

	double **dC=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	dC[D]=new double [DIMENSION];

	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	
	for(int i=0;i<h_num;i++)
	{
		double dFi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		if(J[i]<0){
			dFi[0][0]=-1/pow(-J[i],1/3)*F[i][0*DIMENSION+0];	dFi[0][1]=-1/pow(-J[i],1/3)*F[i][0*DIMENSION+1];	dFi[0][2]=-1/pow(-J[i],1/3)*F[i][0*DIMENSION+2];
			dFi[1][0]=-1/pow(-J[i],1/3)*F[i][1*DIMENSION+0];	dFi[1][1]=-1/pow(-J[i],1/3)*F[i][1*DIMENSION+1];	dFi[1][2]=-1/pow(-J[i],1/3)*F[i][1*DIMENSION+2];
			dFi[2][0]=-1/pow(-J[i],1/3)*F[i][2*DIMENSION+0];	dFi[2][1]=-1/pow(-J[i],1/3)*F[i][2*DIMENSION+1];	dFi[2][2]=-1/pow(-J[i],1/3)*F[i][2*DIMENSION+2];
		}
		else
		{
			dFi[0][0]=1/pow(J[i],1/3)*F[i][0*DIMENSION+0];	dFi[0][1]=1/pow(J[i],1/3)*F[i][0*DIMENSION+1];	dFi[0][2]=1/pow(J[i],1/3)*F[i][0*DIMENSION+2];
			dFi[1][0]=1/pow(J[i],1/3)*F[i][1*DIMENSION+0];	dFi[1][1]=1/pow(J[i],1/3)*F[i][1*DIMENSION+1];	dFi[1][2]=1/pow(J[i],1/3)*F[i][1*DIMENSION+2];
			dFi[2][0]=1/pow(J[i],1/3)*F[i][2*DIMENSION+0];	dFi[2][1]=1/pow(J[i],1/3)*F[i][2*DIMENSION+1];	dFi[2][2]=1/pow(J[i],1/3)*F[i][2*DIMENSION+2];
		}


		dC[0][0]=dFi[0][0]*dFi[0][0]+dFi[1][0]*dFi[1][0]+dFi[2][0]*dFi[2][0];
		dC[0][1]=dFi[0][0]*dFi[0][1]+dFi[1][0]*dFi[1][1]+dFi[2][0]*dFi[2][1];
		dC[0][2]=dFi[0][0]*dFi[0][2]+dFi[1][0]*dFi[1][2]+dFi[2][0]*dFi[2][2];
		dC[1][0]=dFi[0][1]*dFi[0][0]+dFi[1][1]*dFi[1][0]+dFi[2][1]*dFi[2][0];
		dC[1][1]=dFi[0][1]*dFi[0][1]+dFi[1][1]*dFi[1][1]+dFi[2][1]*dFi[2][1];
		dC[1][2]=dFi[0][1]*dFi[0][2]+dFi[1][1]*dFi[1][2]+dFi[2][1]*dFi[2][2];
		dC[2][0]=dFi[0][2]*dFi[0][0]+dFi[1][2]*dFi[1][0]+dFi[2][2]*dFi[2][0];
		dC[2][1]=dFi[0][2]*dFi[0][1]+dFi[1][2]*dFi[1][1]+dFi[2][2]*dFi[2][1];
		dC[2][2]=dFi[0][2]*dFi[0][2]+dFi[1][2]*dFi[1][2]+dFi[2][2]*dFi[2][2];


		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];
		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);

		w[i]=c10*(Ic-3)+c01*(IIc-3);

		double C[DIMENSION][DIMENSION]={{dC[0][0],dC[0][1],dC[0][2]},{dC[1][0],dC[1][1],dC[1][2]},{dC[2][0],dC[2][1],dC[2][2]}};

		inverse(dC,DIMENSION);

		if(J[i]<0){
			s[i][0*DIMENSION+0]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[0][0]) -c01*(C[0][0]-1/3*trace_dC2*dC[0][0]) );
			s[i][0*DIMENSION+1]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][1]) -c01*(C[0][1]-1/3*trace_dC2*dC[0][1]) );
			s[i][0*DIMENSION+2]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][2]) -c01*(C[0][2]-1/3*trace_dC2*dC[0][2]) );

			s[i][1*DIMENSION+0]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][0]) -c01*(C[1][0]-1/3*trace_dC2*dC[1][0]) );
			s[i][1*DIMENSION+1]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[1][1]) -c01*(C[1][1]-1/3*trace_dC2*dC[1][1]) );
			s[i][1*DIMENSION+2]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][2]) -c01*(C[1][2]-1/3*trace_dC2*dC[1][2]) );

			s[i][2*DIMENSION+0]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][0]) -c01*(C[2][0]-1/3*trace_dC2*dC[2][0]) );
			s[i][2*DIMENSION+1]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][1]) -c01*(C[2][1]-1/3*trace_dC2*dC[2][1]) );
			s[i][2*DIMENSION+2]=-1/pow(-J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[2][2]) -c01*(C[2][2]-1/3*trace_dC2*dC[2][2]) );

		}
		else
		{
			s[i][0*DIMENSION+0]=1/pow(J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[0][0]) -c01*(C[0][0]-1/3*trace_dC2*dC[0][0]) );
			s[i][0*DIMENSION+1]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][1]) -c01*(C[0][1]-1/3*trace_dC2*dC[0][1]) );
			s[i][0*DIMENSION+2]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[0][2]) -c01*(C[0][2]-1/3*trace_dC2*dC[0][2]) );

			s[i][1*DIMENSION+0]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][0]) -c01*(C[1][0]-1/3*trace_dC2*dC[1][0]) );
			s[i][1*DIMENSION+1]=1/pow(J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[1][1]) -c01*(C[1][1]-1/3*trace_dC2*dC[1][1]) );
			s[i][1*DIMENSION+2]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[1][2]) -c01*(C[1][2]-1/3*trace_dC2*dC[1][2]) );

			s[i][2*DIMENSION+0]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][0]) -c01*(C[2][0]-1/3*trace_dC2*dC[2][0]) );
			s[i][2*DIMENSION+1]=1/pow(J[i],2/3)*( (c10+c01*Ic)*( -1/3*trace_dC*dC[2][1]) -c01*(C[2][1]-1/3*trace_dC2*dC[2][1]) );
			s[i][2*DIMENSION+2]=1/pow(J[i],2/3)*( (c10+c01*Ic)*(1-1/3*trace_dC*dC[2][2]) -c01*(C[2][2]-1/3*trace_dC2*dC[2][2]) );
		}
	}
	for(int D=0;D>DIMENSION;D++)	delete[]	dC[D];
	delete[]	dC;
}

int MM_method()	
{
	////ó·ëË		ÉVÉXÉeÉÄçHäwëÊ2î≈Å@êXñkèoî≈ÅiäîÅjÅ@ââèKñ‚ëË5ÇÃ4	p.197
	int Nx=2;
	double *x=new double [Nx];
	double *old_x=new double [Nx];
	x[0]=10;
	x[1]=-10;


	double c0=0;
	double fx=0;
	double Txr=0;


	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dTxr=new double [Nx];

	double *rfx=new double [Nx*Nx];
	double *rc0=new double [Nx*Nx];
	double *rTxr=new double [Nx*Nx];

	for(int i=0;i<Nx;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dTxr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			rfx[i*Nx+j]=0;
			rc0[i*Nx+j]=0;
			rTxr[i*Nx+j]=0;
		}
	}


	double r=1;
	double ep=1e-10;
	double seta0=0, seta1=0, seta2=0, seta3=0;
	double E_min=1;
	int count_min=0;

	double *B=new double[Nx*Nx];
	double *d=new double [Nx];	
	double *Nr=new double [Nx];



	while(E_min>ep)
	{
		count_min++;
		if(count_min>5000)	break;
		for(int i=0;i<Nx;i++)
		{
			old_x[i]=x[i];
			d[i]=0;
			Nr[i]=0;
			for(int j=0;j<Nx;j++)
			{
				if(j==i)	B[i*Nx+j]=1;
				else
				{
					B[i*Nx+j]=0;
				}
			}
		}


		double E=1;
		int count=0;
		int calc_way=1;	//ç≈ìKâíTçı	Newton 0, èÄNewton 1

		if(calc_way==0)
		{
			while(E>ep)
			{
				count++;

				fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
				dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);
				rfx[0*Nx+0]=2;	rfx[0*Nx+1]=0;	rfx[1*Nx+0]=0;	rfx[1*Nx+1]=2;

				c0=x[0]*x[0]+x[1]*x[1]-1;
				dc0[0]=2*x[0];	dc0[1]=2*x[1];
				rc0[0*Nx+0]=2;	rc0[0*Nx+1]=0;	rc0[1*Nx+0]=0;	rc0[1*Nx+1]=2;

				Txr=fx;
				Txr+=0.5*r*(c0+seta0)*(c0+seta0);
	
				//cout<<"Txr="<<Txr<<endl;

				for(int i=0;i<Nx;i++)
				{
					dTxr[i]=dfx[i];
					dTxr[i]+=r*(c0+seta0)*dc0[i];
				}

				//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;

				for(int i=0;i<Nx;i++)
				{
					for(int j=0;j<Nx;j++)
					{
						rTxr[i*Nx+j]=rfx[i*Nx+j];
						rTxr[i*Nx+j]+=r*dc0[j]*dc0[i]+r*(c0+seta0)*rc0[i*Nx+j];
					}
				}
				//cout<<"rTxr="<<rTxr[0]<<", "<<rTxr[1]<<", "<<rTxr[2]<<", "<<rTxr[3]<<endl;
				//cout<<rTxr[4]<<", "<<rTxr[5]<<", "<<rTxr[6]<<", "<<rTxr[7]<<endl;
				//cout<<rTxr[8]<<", "<<rTxr[9]<<", "<<rTxr[10]<<", "<<rTxr[11]<<endl;
				//cout<<rTxr[12]<<", "<<rTxr[13]<<", "<<rTxr[14]<<", "<<rTxr[15]<<endl;

				gauss(rTxr,dTxr,Nx);
				//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;
				for(int i=0;i<Nx;i++)	x[i]-=dTxr[i];


				E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
				cout<<"E"<<count<<"="<<E<<endl;
			}
		}
		else if(calc_way==1)
		{
	
			fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
			dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

			c0=x[0]*x[0]+x[1]*x[1]-1;
			dc0[0]=2*x[0];	dc0[1]=2*x[1];

			Txr=fx+0.5*r*(c0+seta0)*(c0+seta0);			
			dTxr[0]=dfx[0]+r*(c0+seta0)*dc0[0];
			dTxr[1]=dfx[1]+r*(c0+seta0)*dc0[1];

			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
			cout<<"E0="<<E<<endl;
			if(E<ep)	return 0;
			while(E>ep)
			{
				count++;
				if(count>500)	break;
				double x_k[2]={x[0], x[1]};
				double dTxr_k[2]={dTxr[0],dTxr[1]};

				Nr[0]=dTxr[0];
				Nr[1]=dTxr[1];
				gauss(B,Nr,Nx);
				d[0]=-1*Nr[0];
				d[1]=-1*Nr[1];

//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;

				double Txr_min=Txr;
				double a_min=1e-3;

				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;
					double x0_a=x[0]+d[0]*alpha;
					double x1_a=x[1]+d[1]*alpha;

					double Txr_a=(x0_a-1)*(x0_a-1)+(x1_a-2)*(x1_a-2)+0.5*r*(x0_a*x0_a+x1_a*x1_a-1+seta0)*(x0_a*x0_a+x1_a*x1_a-1+seta0);

					if(Txr_a<Txr_min)
					{
						Txr_min=Txr_a;
						a_min=alpha;
					}
				}
//				cout<<"Txr"<<count<<"="<<Txr_min<<", alpha="<<a_min<<endl;

				double d0=d[0]*a_min;
				double d1=d[1]*a_min;
				x[0]+=d0;
				x[1]+=d1;

//				cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<endl;


				fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
				dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

				c0=x[0]*x[0]+x[1]*x[1]-1;
				dc0[0]=2*x[0];	dc0[1]=2*x[1];

				Txr=fx;
				Txr+=0.5*r*(c0+seta0)*(c0+seta0);
	
				//cout<<"Txr="<<Txr<<endl;

				for(int i=0;i<Nx;i++)
				{
					dTxr[i]=dfx[i];
					dTxr[i]+=r*(c0+seta0)*dc0[i];
				}
				
				E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
				if(count%100==0)cout<<"E"<<count<<"="<<E<<endl;
				//cout<<endl;

				if(E<ep)	break;

				double s[2]={x[0]-x_k[0],x[1]-x_k[1]};
				double y[2]={dTxr[0]-dTxr_k[0],dTxr[1]-dTxr_k[1]};

				double beta=y[0]*s[0]+y[1]*s[1];
				double sigma=(s[0]*B[0]+s[1]*B[2])*s[0]+(s[0]*B[2]+s[1]*B[3])*s[1];

//				cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

				if(beta>0)//(beta>=0.2*sigma)
				{
					double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
					double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

					B[0*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
					B[0*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
					B[1*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
					B[1*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
				}
	/*			else if(beta>0)
				{
					double seta=0.8*sigma/(sigma-beta);
					cout<<"seta="<<seta<<endl;
					double y0_s=y[0]*seta+(1-seta)*(B[0]*s[0]+B[1]*s[1]);
					double y1_s=y[1]*seta+(1-seta)*(B[2]*s[0]+B[3]*s[1]);

					beta=y0_s*s[0]+y1_s*s[1];
					double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
					double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

					B[0*Nx+0]+=1/beta*y0_s*y0_s-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
					B[0*Nx+1]+=1/beta*y0_s*y1_s-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
					B[1*Nx+0]+=1/beta*y1_s*y0_s-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
					B[1*Nx+1]+=1/beta*y1_s*y1_s-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
				}*/
			}
		}


		E_min=sqrt( (old_x[0]-x[0])*(old_x[0]-x[0]) + (old_x[1]-x[1])*(old_x[1]-x[1]));
		if(E_min<ep*1000)	r*=4;
		seta0+=c0;

		cout<<"E_min"<<count_min<<"="<<E_min<<endl<<endl;
	}

	//v0Å@åvéZ
	int nv=0;
	cout<<"fc="<<c0<<endl;
	double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	cout<<"c="<<c0<<endl;
	cout<<"v="<<v0<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"óùò_â"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/


	delete[]	x;
	delete[]	old_x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dTxr;
	delete[]	rfx;
	delete[]	rc0;
	delete[]	rTxr;
	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}

