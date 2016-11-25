#include "stdafx.h"	
void calc_nabla_laplacian(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dL,double *rL,double **dg,double **dh,double **d_dgdt,double **d_dhdt);
void calc_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void output_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,double *rT,double *dT,double *rL,double *dL,int Nx,int h_num,int count);

void calc_HYPER_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{

	////////////íËã`///////////////
	int h_num=HYPER.size();
	int Nx=h_num*4;

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
	double nG[DIMENSION]={0,0,1};
	double aG[DIMENSION]={0,0,0};

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

	ofstream f_np("p_n.csv",ios::trunc);
	ofstream f_nq("q_n.csv",ios::trunc);

	for(int i=0;i<h_num;i++)
	{
		f_np<<HYPER[i].p_n[A_X]<<","<<HYPER[i].p_n[A_Y]<<","<<HYPER[i].p_n[A_Z]<<endl;
		f_nq<<HYPER[i].q_n[A_X]<<","<<HYPER[i].q_n[A_Y]<<","<<HYPER[i].q_n[A_Z]<<endl;
		dg[i]=new double [Nx];
		dh[i]=new double [Nx];
		d_dgdt[i]=new double [Nx];
		d_dhdt[i]=new double [Nx];
	}
	f_np.close();
	f_nq.close();

	double r=0.1;
	double ep=1e-10;
	double *th_g=new double [h_num];
	double *th_h=new double [h_num];
	double *th_dgdt=new double [h_num];
	double *th_dhdt=new double [h_num];


	////////////èâä˙âªéZ///////////////
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
		}
	}
	for(int i=h_num;i<Nx;i++)
	{
		df[i]=0;
		dT[i]=0;
		dL[i]=0;
		for(int j=0;j<Nx;j++)
		{
			rf[i*Nx+j]=0;
			rT[i*Nx+j]=0;
			rL[i*Nx+j]=0;
		}
	}








	////////////QPåvéZ///////////////		
	cout<<"QP START----------------------"<<endl;		
	double E_min=1;
	int count_min=0;
	while(E_min>ep)
	{
		count_min++;

		///////ëOåvéZåãâ ï€ä«
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

			////////////////////////////////t=n+1 ïœêîåvéZ//////////////////////////////	/////t=nÅ@ïœêîåvéZåãâ ÇÕä˘ílÇ∆Ç∑ÇÈÅBâ∫ãLÇ≈ÇÃåvéZÇÕÇ»ÇµÅB
			//cout<<"calc_variables";
			calc_variables(CON,PART,HYPER,HYPER1);
			//cout<<"èIóπ"<<endl;
			//cout<<"calc_nabla_laplacian";
			calc_nabla_laplacian(CON,PART,HYPER,HYPER1,dL,rL,dg,dh,d_dgdt,d_dhdt);
			//cout<<"èIóπ"<<endl;


			////////////////////////////////ñ⁄ìIä÷êîåvéZ/////////////////////////////////
			T=0;
			L=0;
			for(int i=0;i<h_num;i++)
			{
				L+=0.5/mi*(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z])-V*HYPER[i].W;

				g[i]=V*(1-HYPER[i].J);

				dgdt[i]=0;
				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{
					int k=HYPER[i].NEI[j];
					dgdt[i]+=1/mi*(HYPER1[i*h_num+k].DgDq[A_X]*HYPER[i].p[A_X]+HYPER1[i*h_num+k].DgDq[A_Y]*HYPER[i].p[A_Y]+HYPER1[i*h_num+k].DgDq[A_Z]*HYPER[i].p[A_Z]);
				}
				dgdt[i]+=1/mi*(HYPER1[i*h_num+i].DgDq[A_X]*HYPER[i].p[A_X]+HYPER1[i*h_num+i].DgDq[A_Y]*HYPER[i].p[A_Y]+HYPER1[i*h_num+i].DgDq[A_Z]*HYPER[i].p[A_Z]);
				h[i]=-1*((PART[i].r[A_X]-aG[A_X])*nG[A_X]+(PART[i].r[A_Y]-aG[A_Y])*nG[A_Y]+(PART[i].r[A_Z]-aG[A_Z])*nG[A_Z]);
				dhdt[i]=-1/mi*(HYPER[i].p[A_X]*nG[A_X]+HYPER[i].p[A_Y]*nG[A_Y]+HYPER[i].p[A_Z]*nG[A_Z]);

				T+=0.5*r*( (g[i]+th_g[i])*(g[i]+th_g[i]) + (dgdt[i]+th_dgdt[i])*(dgdt[i]+th_dgdt[i]) );
				if(h[i]+th_h[i]>0)	T+=0.5*r*(h[i]+th_h[i])*(h[i]+th_h[i]);
				if(dhdt[i]+th_dhdt[i]>0)	T+=0.5*r*(dhdt[i]+th_dhdt[i])*(dhdt[i]+th_dhdt[i]);
			}			
			f=(L-L_n)*(L-L_n);
			T=+f;

			cout<<"L="<<L<<", Ln="<<L_n<<endl;

			for(int i=0;i<h_num;i++)
			{
				////////////å˘îzåvéZ
				df[i]=2*dL[i]*(L-L_n);	df[i+h_num]=2*dL[i+h_num]*(L-L_n);	df[i+2*h_num]=2*dL[i+2*h_num]*(L-L_n);	df[i+3*h_num]=2*dL[i+3*h_num]*(L-L_n);
				dT[i]=df[i];	dT[i+h_num]=df[i+h_num];	dT[i+2*h_num]=df[i+2*h_num];	dT[i+3*h_num]=df[i+3*h_num];

				/////////////î≠éUåvéZ
				for(int j=0;j<h_num;j++)
				{
					////////////å˘îzåvéZ
					dT[i]+=r*( dg[j][i]*(g[j]+th_g[j]) + d_dgdt[j][i]*(dgdt[j]+th_dgdt[j]) );
					dT[i+h_num]+=r*d_dgdt[j][i+h_num]*(dgdt[j]+th_dgdt[j]);
					dT[i+2*h_num]+=r*( dg[j][i+2*h_num]*(g[j]+th_g[j]) + d_dgdt[j][i+2*h_num]*(dgdt[j]+th_dgdt[j]) );
					dT[i+3*h_num]+=r*d_dgdt[j][i+3*h_num]*(dgdt[j]+th_dgdt[j]);

					if(h[j]+th_h[j]>0)
					{
						dT[i]+=r*dh[j][i]*(h[j]+th_h[j]);
						dT[i+2*h_num]+=r*dh[j][i+2*h_num]*(h[j]+th_h[j]);
					}
					if(dhdt[j]+th_dhdt[j]>0)
					{
						dT[i]+=r*d_dhdt[j][i]*(dhdt[j]+th_dhdt[j]);
						dT[i+h_num]+=r*d_dhdt[j][i+h_num]*(dhdt[j]+th_dhdt[j]);
						dT[i+2*h_num]+=r*d_dhdt[j][i+2*h_num]*(dhdt[j]+th_dhdt[j]);
						dT[i+3*h_num]+=r*d_dhdt[j][i+3*h_num]*(dhdt[j]+th_dhdt[j]);
					}

					/////////////î≠éUåvéZ
					rT[i*Nx+j]=2*(rL[i*Nx+j]*(L-L_n)+dL[i]*dL[j]);							rT[(i+h_num)*Nx+j]=2*(rL[(i+h_num)*Nx+j]*(L-L_n)+dL[i+h_num]*dL[j]);							rT[(i+2*h_num)*Nx+j]=2*(rL[(i+2*h_num)*Nx+j]*(L-L_n)+dL[i+2*h_num]*dL[j]);							rT[(i+3*h_num)*Nx+j]=2*(rL[(i+3*h_num)*Nx+j]*(L-L_n)+dL[i+3*h_num]*dL[j]);
					rT[i*Nx+j+h_num]=2*(rL[i*Nx+j+h_num]*(L-L_n)+dL[i]*dL[j+h_num]);		rT[(i+h_num)*Nx+j+h_num]=2*(rL[(i+h_num)*Nx+j+h_num]*(L-L_n)+dL[i+h_num]*dL[j+h_num]);			rT[(i+2*h_num)*Nx+j+h_num]=2*(rL[(i+2*h_num)*Nx+j+h_num]*(L-L_n)+dL[i+2*h_num]*dL[j+h_num]);		rT[(i+3*h_num)*Nx+j+h_num]=2*(rL[(i+3*h_num)*Nx+j+h_num]*(L-L_n)+dL[i+3*h_num]*dL[j+h_num]);
					rT[i*Nx+j+2*h_num]=2*(rL[i*Nx+j+2*h_num]*(L-L_n)+dL[i]*dL[j+2*h_num]);	rT[(i+h_num)*Nx+j+2*h_num]=2*(rL[(i+h_num)*Nx+j+2*h_num]*(L-L_n)+dL[i+h_num]*dL[j+2*h_num]);	rT[(i+2*h_num)*Nx+j+2*h_num]=2*(rL[(i+2*h_num)*Nx+j+2*h_num]*(L-L_n)+dL[i+2*h_num]*dL[j+2*h_num]);	rT[(i+3*h_num)*Nx+j+2*h_num]=2*(rL[(i+3*h_num)*Nx+j+2*h_num]*(L-L_n)+dL[i+3*h_num]*dL[j+2*h_num]);
					rT[i*Nx+j+3*h_num]=2*(rL[i*Nx+j+3*h_num]*(L-L_n)+dL[i]*dL[j+3*h_num]);	rT[(i+h_num)*Nx+j+3*h_num]=2*(rL[(i+h_num)*Nx+j+3*h_num]*(L-L_n)+dL[i+h_num]*dL[j+3*h_num]);	rT[(i+2*h_num)*Nx+j+3*h_num]=2*(rL[(i+2*h_num)*Nx+j+3*h_num]*(L-L_n)+dL[i+2*h_num]*dL[j+3*h_num]);	rT[(i+3*h_num)*Nx+j+3*h_num]=2*(rL[(i+3*h_num)*Nx+j+3*h_num]*(L-L_n)+dL[i+3*h_num]*dL[j+3*h_num]);

					for(int k=0;k<h_num;k++)
					{
						rT[i*Nx+j]+=r*(dg[k][i]*dg[k][j]+d_dgdt[k][i]*d_dgdt[k][j]);	rT[(i+h_num)*Nx+j]+=r*d_dgdt[k][i+h_num]*d_dgdt[k][j];					rT[(i+2*h_num)*Nx+j]+=r*(dg[k][i+2*h_num]*dg[k][j]+d_dgdt[k][i+2*h_num]*d_dgdt[k][j]);		rT[(i+3*h_num)*Nx+j]+=r*d_dgdt[k][i+3*h_num]*d_dgdt[k][j];	
						rT[i*Nx+j+h_num]+=r*d_dgdt[k][i]*d_dgdt[k][j+h_num];			rT[(i+h_num)*Nx+j+h_num]+=r*+d_dgdt[k][i+h_num]*d_dgdt[k][j+h_num];		rT[(i+2*h_num)*Nx+j+h_num]+=r*d_dgdt[k][i+2*h_num]*d_dgdt[k][j+h_num];						rT[(i+3*h_num)*Nx+j+h_num]+=r*d_dgdt[k][i+3*h_num]*d_dgdt[k][j+h_num];	
						rT[i*Nx+j+2*h_num]+=r*(dg[k][i]*dg[k][j+2*h_num]+d_dgdt[k][i]*d_dgdt[k][j+2*h_num]);	rT[(i+h_num)*Nx+j+2*h_num]+=r*d_dgdt[k][i+h_num]*d_dgdt[k][j+2*h_num];	rT[(i+2*h_num)*Nx+j+2*h_num]+=r*(dg[k][i+2*h_num]*dg[k][j+2*h_num]+d_dgdt[k][i+2*h_num]*d_dgdt[k][j+2*h_num]);		rT[(i+3*h_num)*Nx+j+2*h_num]+=r*d_dgdt[k][i+3*h_num]*d_dgdt[k][j+2*h_num];	
						rT[i*Nx+j+3*h_num]+=r*d_dgdt[k][i]*d_dgdt[k][j+3*h_num];		rT[(i+h_num)*Nx+j+3*h_num]+=r*d_dgdt[k][i+h_num]*d_dgdt[k][j+3*h_num];	rT[(i+2*h_num)*Nx+j+3*h_num]+=r*d_dgdt[k][i+2*h_num]*d_dgdt[k][j+3*h_num];					rT[(i+3*h_num)*Nx+j+3*h_num]+=r*d_dgdt[k][i+3*h_num]*d_dgdt[k][j+3*h_num];	

						if(h[k]+th_h[k]>0)
						{
							rT[i*Nx+j]+=r*dh[k][i]*dh[k][j];					rT[(i+2*h_num)*Nx+j]+=r*dh[k][i+2*h_num]*dh[k][j];
							rT[i*Nx+j+2*h_num]+=r*dh[k][i]*dh[k][j+2*h_num];	rT[(i+2*h_num)*Nx+j+2*h_num]+=r*dh[k][i+2*h_num]*dh[k][j+2*h_num];							
						}
						else if(h[k]+th_h[k]==0)
						{
							if(dh[k][i]*dh[k][j]>0)
							{
								rT[i*Nx+j]+=r*dh[k][i]*dh[k][j];					rT[(i+2*h_num)*Nx+j]+=r*dh[k][i+2*h_num]*dh[k][j];
								rT[i*Nx+j+2*h_num]+=r*dh[k][i]*dh[k][j+2*h_num];	rT[(i+2*h_num)*Nx+j+2*h_num]+=r*dh[k][i+2*h_num]*dh[k][j+2*h_num];							
							}
						}
						if(dhdt[k]+th_dhdt[k]>0)
						{
							rT[i*Nx+j]+=r*d_dhdt[k][i]*d_dhdt[k][j];					rT[(i+h_num)*Nx+j]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j];					rT[(i+2*h_num)*Nx+j]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j];					rT[(i+3*h_num)*Nx+j]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j];
							rT[i*Nx+j+h_num]+=r*d_dhdt[k][i]*d_dhdt[k][j+h_num];		rT[(i+h_num)*Nx+j+h_num]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j+h_num];		rT[(i+2*h_num)*Nx+j+h_num]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j+h_num];		rT[(i+3*h_num)*Nx+j+h_num]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j+h_num];
							rT[i*Nx+j+2*h_num]+=r*d_dhdt[k][i]*d_dhdt[k][j+2*h_num];	rT[(i+h_num)*Nx+j+2*h_num]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j+2*h_num];	rT[(i+2*h_num)*Nx+j+2*h_num]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j+2*h_num];	rT[(i+3*h_num)*Nx+j+2*h_num]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j+2*h_num];
							rT[i*Nx+j+3*h_num]+=r*d_dhdt[k][i]*d_dhdt[k][j+3*h_num];	rT[(i+h_num)*Nx+j+3*h_num]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j+3*h_num];	rT[(i+2*h_num)*Nx+j+3*h_num]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j+3*h_num];	rT[(i+3*h_num)*Nx+j+3*h_num]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j+3*h_num];
						}
						else if(dhdt[k]+th_dhdt[k]==0)
						{
							if(d_dhdt[k][i]*d_dhdt[k][j]>0)
							{
								rT[i*Nx+j]+=r*d_dhdt[k][i]*d_dhdt[k][j];					rT[(i+h_num)*Nx+j]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j];					rT[(i+2*h_num)*Nx+j]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j];					rT[(i+3*h_num)*Nx+j]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j];
								rT[i*Nx+j+h_num]+=r*d_dhdt[k][i]*d_dhdt[k][j+h_num];		rT[(i+h_num)*Nx+j+h_num]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j+h_num];		rT[(i+2*h_num)*Nx+j+h_num]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j+h_num];		rT[(i+3*h_num)*Nx+j+h_num]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j+h_num];
								rT[i*Nx+j+2*h_num]+=r*d_dhdt[k][i]*d_dhdt[k][j+2*h_num];	rT[(i+h_num)*Nx+j+2*h_num]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j+2*h_num];	rT[(i+2*h_num)*Nx+j+2*h_num]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j+2*h_num];	rT[(i+3*h_num)*Nx+j+2*h_num]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j+2*h_num];
								rT[i*Nx+j+3*h_num]+=r*d_dhdt[k][i]*d_dhdt[k][j+3*h_num];	rT[(i+h_num)*Nx+j+3*h_num]+=r*d_dhdt[k][i+h_num]*d_dhdt[k][j+3*h_num];	rT[(i+2*h_num)*Nx+j+3*h_num]+=r*d_dhdt[k][i+2*h_num]*d_dhdt[k][j+3*h_num];	rT[(i+3*h_num)*Nx+j+3*h_num]+=r*d_dhdt[k][i+3*h_num]*d_dhdt[k][j+3*h_num];
							}
						}
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

			output_data(PART,HYPER,rT,dT,rL,dL,Nx,h_num,count);

			double p_E=0;
			for(int i=0;i<Nx;i++)	p_E+=dT[i]*dT[i];
			E=sqrt(p_E);
			cout<<"	E"<<count<<"="<<E<<endl;
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
		cout<<"E_min"<<count_min<<"="<<E_min<<endl;
	}
	cout<<"OK"<<endl;


	
	///////////ÉÅÉÇÉäâï˙///////////////		
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
	}
	delete[]	dg;
	delete[]	dh;
	delete[]	d_dgdt;
	delete[]	d_dhdt;

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


	for(int l=0;l<h_num;l++)
	{
		/////////////////dLåvéZ///////////////////
		dL[l]=0;
		dL[l+h_num]=0;
		dL[l+2*h_num]=0;
		dL[l+3*h_num]=0;

		int Nl=HYPER[l].N;
		for(int i=0;i<Nl;i++)
		{
			int in=HYPER[l].NEI[i];
			dL[l]+=-0.5*Dt/mi*(HYPER[in].p[A_X]*HYPER1[l*h_num+in].DgDq_n[A_X]+HYPER[in].p[A_Y]*HYPER1[l*h_num+in].DgDq_n[A_Y]+HYPER[in].p[A_Z]*HYPER1[l*h_num+in].DgDq_n[A_Z]);			
			dL[l+h_num]+=-0.5*Dt/mi*(HYPER[in].p[A_X]*HYPER1[l*h_num+in].DgDq[A_X]+HYPER[in].p[A_Y]*HYPER1[l*h_num+in].DgDq[A_Y]+HYPER[in].p[A_Z]*HYPER1[l*h_num+in].DgDq[A_Z]);

			dL[l+2*h_num]+=0.5*Dt*Dt/mi*(
			 (HYPER[in].stress[A_X][0]*HYPER1[in*h_num+l].DgDq[0]+HYPER[in].stress[A_X][1]*HYPER1[in*h_num+l].DgDq[1]+HYPER[in].stress[A_X][2]*HYPER1[in*h_num+l].DgDq[2])*nG[A_X]
			+(HYPER[in].stress[A_Y][0]*HYPER1[in*h_num+l].DgDq[0]+HYPER[in].stress[A_Y][1]*HYPER1[in*h_num+l].DgDq[1]+HYPER[in].stress[A_Y][2]*HYPER1[in*h_num+l].DgDq[2])*nG[A_Y]
			+(HYPER[in].stress[A_Z][0]*HYPER1[in*h_num+l].DgDq[0]+HYPER[in].stress[A_Z][1]*HYPER1[in*h_num+l].DgDq[1]+HYPER[in].stress[A_Z][2]*HYPER1[in*h_num+l].DgDq[2])*nG[A_Z] );

			int Nin=HYPER[in].N;
			for(int j=0;j<Nin;j++)
			{
				int jn=HYPER[in].NEI[j];
				dL[l]+=-0.5*Dt*Dt/mi*(
				 (HYPER[jn].stress[A_X][0]*HYPER1[jn*h_num+in].DgDq[0]+HYPER[jn].stress[A_X][1]*HYPER1[jn*h_num+in].DgDq[1]+HYPER[jn].stress[A_X][2]*HYPER1[jn*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_X]
				+(HYPER[jn].stress[A_Y][0]*HYPER1[jn*h_num+in].DgDq[0]+HYPER[jn].stress[A_Y][1]*HYPER1[jn*h_num+in].DgDq[1]+HYPER[jn].stress[A_Y][2]*HYPER1[jn*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_Y]
				+(HYPER[jn].stress[A_Z][0]*HYPER1[jn*h_num+in].DgDq[0]+HYPER[jn].stress[A_Z][1]*HYPER1[jn*h_num+in].DgDq[1]+HYPER[jn].stress[A_Z][2]*HYPER1[jn*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_Z] );
			}
			dL[l]+=-0.5*Dt*Dt/mi*(
			 (HYPER[in].stress[A_X][0]*HYPER1[in*h_num+in].DgDq[0]+HYPER[in].stress[A_X][1]*HYPER1[in*h_num+in].DgDq[1]+HYPER[in].stress[A_X][2]*HYPER1[in*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_X]
			+(HYPER[in].stress[A_Y][0]*HYPER1[in*h_num+in].DgDq[0]+HYPER[in].stress[A_Y][1]*HYPER1[in*h_num+in].DgDq[1]+HYPER[in].stress[A_Y][2]*HYPER1[in*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_Y]
			+(HYPER[in].stress[A_Z][0]*HYPER1[in*h_num+in].DgDq[0]+HYPER[in].stress[A_Z][1]*HYPER1[in*h_num+in].DgDq[1]+HYPER[in].stress[A_Z][2]*HYPER1[in*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_Z] );
		}
		dL[l]+=-0.5*Dt*Dt/mi*(
		 (HYPER[in].stress[A_X][0]*HYPER1[in*h_num+in].DgDq[0]+HYPER[in].stress[A_X][1]*HYPER1[in*h_num+in].DgDq[1]+HYPER[in].stress[A_X][2]*HYPER1[in*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_X]
		+(HYPER[in].stress[A_Y][0]*HYPER1[in*h_num+in].DgDq[0]+HYPER[in].stress[A_Y][1]*HYPER1[in*h_num+in].DgDq[1]+HYPER[in].stress[A_Y][2]*HYPER1[in*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_Y]
		+(HYPER[in].stress[A_Z][0]*HYPER1[in*h_num+in].DgDq[0]+HYPER[in].stress[A_Z][1]*HYPER1[in*h_num+in].DgDq[1]+HYPER[in].stress[A_Z][2]*HYPER1[in*h_num+in].DgDq[2])*HYPER1[l*h_num+in].DgDq_n[A_Z] );

		dL[l]+=-0.5*Dt/mi*(HYPER[l].p[A_X]*HYPER1[l*h_num+l].DgDq_n[A_X]+HYPER[l].p[A_Y]*HYPER1[l*h_num+l].DgDq_n[A_Y]+HYPER[l].p[A_Z]*HYPER1[l*h_num+l].DgDq_n[A_Z]);
		dL[l+h_num]+=-0.5*Dt/mi*(HYPER[l].p[A_X]*HYPER1[l*h_num+l].DgDq[A_X]+HYPER[l].p[A_Y]*HYPER1[l*h_num+l].DgDq[A_Y]+HYPER[l].p[A_Z]*HYPER1[l*h_num+l].DgDq[A_Z]);
		dL[l+2*h_num]+=0.5*Dt*Dt/mi*(
		 (HYPER[l].stress[A_X][0]*HYPER1[l*h_num+l].DgDq[0]+HYPER[l].stress[A_X][1]*HYPER1[l*h_num+l].DgDq[1]+HYPER[l].stress[A_X][2]*HYPER1[l*h_num+l].DgDq[2])*nG[A_X]
		+(HYPER[l].stress[A_Y][0]*HYPER1[l*h_num+l].DgDq[0]+HYPER[l].stress[A_Y][1]*HYPER1[l*h_num+l].DgDq[1]+HYPER[l].stress[A_Y][2]*HYPER1[l*h_num+l].DgDq[2])*nG[A_Y]
		+(HYPER[l].stress[A_Z][0]*HYPER1[l*h_num+l].DgDq[0]+HYPER[l].stress[A_Z][1]*HYPER1[l*h_num+l].DgDq[1]+HYPER[l].stress[A_Z][2]*HYPER1[l*h_num+l].DgDq[2])*nG[A_Z] );

		dL[l+2*h_num]+=0.5*Dt/mi*(HYPER[l].p[A_X]*nG[A_X]+HYPER[l].p[A_Y]*nG[A_Y]+HYPER[l].p[A_Z]*nG[A_Z]);
		dL[l+3*h_num]=0.5*Dt/mi*(HYPER[l].p[A_X]*nG[A_X]+HYPER[l].p[A_Y]*nG[A_Y]+HYPER[l].p[A_Z]*nG[A_Z]);

	
		
		
		
		/////////////////rLåvéZ1///////////////////
		for(int k=0;k<h_num;k++)
		{
			rL[l*Nx+k]=0;			rL[(l+h_num)*Nx+k]=0;			rL[(l+2*h_num)*Nx+k]=0;			rL[(l+3*h_num)*Nx+k]=0;
			rL[l*Nx+k+h_num]=0;		rL[(l+h_num)*Nx+k+h_num]=0;		rL[(l+2*h_num)*Nx+k+h_num]=0;	rL[(l+3*h_num)*Nx+k+h_num]=0;
			rL[l*Nx+k+2*h_num]=0;	rL[(l+h_num)*Nx+k+2*h_num]=0;	rL[(l+2*h_num)*Nx+k+2*h_num]=0;	rL[(l+3*h_num)*Nx+k+2*h_num]=0;
			rL[l*Nx+k+3*h_num]=0;	rL[(l+h_num)*Nx+k+3*h_num]=0;	rL[(l+2*h_num)*Nx+k+3*h_num]=0;	rL[(l+3*h_num)*Nx+k+3*h_num]=0;
			
			dg[k][l]=0;		dg[k][l+2*h_num]=0;
			d_dgdt[k][l]=0;	d_dgdt[k][l+h_num]=0;	d_dgdt[k][l+2*h_num]=0;	d_dgdt[k][l+3*h_num]=0;
			dh[k][l]=0;		dh[k][l+2*h_num]=0;
			d_dhdt[k][l]=0;	d_dhdt[k][l+h_num]=0;	d_dhdt[k][l+2*h_num]=0;	d_dhdt[k][l+3*h_num]=0;

			int Nk=HYPER[k].N;
			for(int i=0;i<Nk;i++)
			{
				int in=HYPER[k].NEI[i];
	
				dg[k][l]+=-0.5*Dt*Dt/mi*(HYPER1[k*h_num+in].DgDq[A_X]*HYPER1[l*h_num+in].DgDq_n[A_X]+HYPER1[k*h_num+in].DgDq[A_Y]*HYPER1[l*h_num+in].DgDq_n[A_Y]+HYPER1[k*h_num+in].DgDq[A_Z]*HYPER1[l*h_num+in].DgDq_n[A_Z]);				
	
				d_dgdt[k][l+h_num]+=-0.5*Dt/mi*(HYPER1[k*h_num+in].DgDq[A_X]*HYPER1[l*h_num+in].DgDq[A_X]+HYPER1[k*h_num+in].DgDq[A_Y]*HYPER1[l*h_num+in].DgDq[A_Y]+HYPER1[k*h_num+in].DgDq[A_Z]*HYPER1[l*h_num+in].DgDq[A_Z]);				

				rL[l*Nx+k]+=0.25*Dt*Dt/mi*(HYPER1[k*h_num+in].DgDq_n[A_X]*HYPER1[l*h_num+in].DgDq_n[A_X]+HYPER1[k*h_num+in].DgDq_n[A_Y]*HYPER1[l*h_num+in].DgDq_n[A_Y]+HYPER1[k*h_num+in].DgDq_n[A_Z]*HYPER1[l*h_num+in].DgDq_n[A_Z]);

				rL[(l+h_num)*Nx+k]+=0.25*Dt*Dt/mi*(HYPER1[k*h_num+in].DgDq_n[A_X]*HYPER1[l*h_num+in].DgDq[A_X]+HYPER1[k*h_num+in].DgDq_n[A_Y]*HYPER1[l*h_num+in].DgDq[A_Y]+HYPER1[k*h_num+in].DgDq_n[A_Z]*HYPER1[l*h_num+in].DgDq[A_Z]);

				rL[l*Nx+k+h_num]+=0.25*Dt*Dt/mi*(HYPER1[l*h_num+in].DgDq_n[A_X]*HYPER1[k*h_num+in].DgDq[A_X]+HYPER1[l*h_num+in].DgDq_n[A_Y]*HYPER1[k*h_num+in].DgDq[A_Y]+HYPER1[l*h_num+in].DgDq_n[A_Z]*HYPER1[k*h_num+in].DgDq[A_Z]);
				rL[(l+h_num)*Nx+k+h_num]+=0.25*Dt*Dt/mi*(HYPER1[l*h_num+in].DgDq[A_X]*HYPER1[k*h_num+in].DgDq[A_X]+HYPER1[l*h_num+in].DgDq[A_Y]*HYPER1[k*h_num+in].DgDq[A_Y]+HYPER1[l*h_num+in].DgDq[A_Z]*HYPER1[k*h_num+in].DgDq[A_Z]);

				
				rL[(l+2*h_num)*Nx+k+2*h_num]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*
				((HYPER[in].dPIdF[A_X][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_X][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_X][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*nG[A_X]
				+(HYPER[in].dPIdF[A_Y][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*nG[A_Y]
				+(HYPER[in].dPIdF[A_Z][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*nG[A_Z])
				*(HYPER1[k*h_num+in].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+in].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+in].n0ij[A_Z]*nG[A_Z]);

				int Nin=HYPER[in].N;
				for(int j=0;j<Nin;j++)
				{
					int jn=HYPER[in].NEI[j];
					
					rL[(l+2*h_num)*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
					((HYPER[jn].dPIdF[A_X][A_X]*HYPER1[l*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_X][A_Y]*HYPER1[l*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_X][A_Z]*HYPER1[l*h_num+jn].n0ij[A_Z])*nG[A_X]
					+(HYPER[jn].dPIdF[A_Y][A_X]*HYPER1[l*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+jn].n0ij[A_Z])*nG[A_X]
					+(HYPER[jn].dPIdF[A_Z][A_X]*HYPER1[l*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+jn].n0ij[A_Z])*nG[A_X])
					*(HYPER1[in*h_num+jn].n0ij[A_X]*HYPER1[k*h_num+in].DgDq_n[A_X]+HYPER1[in*h_num+jn].n0ij[A_Y]*HYPER1[k*h_num+in].DgDq_n[A_Y]+HYPER1[in*h_num+jn].n0ij[A_Z]*HYPER1[k*h_num+in].DgDq_n[A_Z]);

					rL[l*Nx+k+2*h_num]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
					((HYPER[jn].dPIdF[A_X][A_X]*HYPER1[in*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_X][A_Y]*HYPER1[in*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_X][A_Z]*HYPER1[in*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_X]
					+(HYPER[jn].dPIdF[A_Y][A_X]*HYPER1[in*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Y][A_Y]*HYPER1[in*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Y][A_Z]*HYPER1[in*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_Y]
					+(HYPER[jn].dPIdF[A_Z][A_X]*HYPER1[in*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Z][A_Y]*HYPER1[in*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Z][A_Z]*HYPER1[in*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_Z])
					*(HYPER1[k*h_num+jn].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+jn].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+jn].n0ij[A_Z]*nG[A_Z]);

					int Nl=HYPER[l].N;
					for(int m=0;m<Nl;m++)
					{
						int mn=HYPER[l].NEI[m];
						rL[l*Nx+k]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*
						((HYPER[jn].dPIdF[A_X][A_X]*HYPER1[mn*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_X][A_Y]*HYPER1[mn*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_X][A_Z]*HYPER1[mn*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+mn].DgDq_n[A_X]
						+(HYPER[jn].dPIdF[A_Y][A_X]*HYPER1[mn*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Y][A_Y]*HYPER1[mn*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Y][A_Z]*HYPER1[mn*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+mn].DgDq_n[A_Y]
						+(HYPER[jn].dPIdF[A_Z][A_X]*HYPER1[mn*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Z][A_Y]*HYPER1[mn*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Z][A_Z]*HYPER1[mn*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+mn].DgDq_n[A_Z])
						*(HYPER1[in*h_num+jn].n0ij[A_X]*HYPER1[k*h_num+in].DgDq_n[A_X]+HYPER1[in*h_num+jn].n0ij[A_Y]*HYPER1[k*h_num+in].DgDq_n[A_Y]+HYPER1[in*h_num+jn].n0ij[A_Z]*HYPER1[k*h_num+in].DgDq_n[A_Z]);
					}
					rL[l*Nx+k]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*
					((HYPER[jn].dPIdF[A_X][A_X]*HYPER1[jn*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_X][A_Y]*HYPER1[jn*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_X][A_Z]*HYPER1[jn*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+jn].DgDq_n[A_X]
					+(HYPER[jn].dPIdF[A_Y][A_X]*HYPER1[jn*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Y][A_Y]*HYPER1[jn*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Y][A_Z]*HYPER1[jn*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+jn].DgDq_n[A_Y]
					+(HYPER[jn].dPIdF[A_Z][A_X]*HYPER1[jn*h_num+jn].n0ij[A_X]+HYPER[jn].dPIdF[A_Z][A_Y]*HYPER1[jn*h_num+jn].n0ij[A_Y]+HYPER[jn].dPIdF[A_Z][A_Z]*HYPER1[jn*h_num+jn].n0ij[A_Z])*HYPER1[l*h_num+jn].DgDq_n[A_Z])
					*(HYPER1[in*h_num+jn].n0ij[A_X]*HYPER1[k*h_num+in].DgDq_n[A_X]+HYPER1[in*h_num+jn].n0ij[A_Y]*HYPER1[k*h_num+in].DgDq_n[A_Y]+HYPER1[in*h_num+jn].n0ij[A_Z]*HYPER1[k*h_num+in].DgDq_n[A_Z]);
				}
				rL[l*Nx+k]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*
				((HYPER[in].dPIdF[A_X][A_X]*HYPER1[in*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_X][A_Y]*HYPER1[in*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_X][A_Z]*HYPER1[in*h_num+in].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_X]
				+(HYPER[in].dPIdF[A_Y][A_X]*HYPER1[in*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Y][A_Y]*HYPER1[in*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Y][A_Z]*HYPER1[in*h_num+in].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_Y]
				+(HYPER[in].dPIdF[A_Z][A_X]*HYPER1[in*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Z][A_Y]*HYPER1[in*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Z][A_Z]*HYPER1[in*h_num+in].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_Z])
				*(HYPER1[in*h_num+in].n0ij[A_X]*HYPER1[k*h_num+in].DgDq_n[A_X]+HYPER1[in*h_num+in].n0ij[A_Y]*HYPER1[k*h_num+in].DgDq_n[A_Y]+HYPER1[in*h_num+in].n0ij[A_Z]*HYPER1[k*h_num+in].DgDq_n[A_Z]);

				rL[(l+2*h_num)*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
				((HYPER[in].dPIdF[A_X][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_X][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_X][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*nG[A_X]
				+(HYPER[in].dPIdF[A_Y][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*nG[A_Y]
				+(HYPER[in].dPIdF[A_Z][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*nG[A_Z])
				*(HYPER1[in*h_num+in].n0ij[A_X]*HYPER1[k*h_num+in].DgDq_n[A_X]+HYPER1[in*h_num+in].n0ij[A_Y]*HYPER1[k*h_num+in].DgDq_n[A_Y]+HYPER1[in*h_num+in].n0ij[A_Z]*HYPER1[k*h_num+in].DgDq_n[A_Z]);

				rL[l*Nx+k+2*h_num]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
				((HYPER[in].dPIdF[A_X][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_X][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_X][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_X]
				+(HYPER[in].dPIdF[A_Y][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_Y]
				+(HYPER[in].dPIdF[A_Z][A_X]*HYPER1[l*h_num+in].n0ij[A_X]+HYPER[in].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+in].n0ij[A_Y]+HYPER[in].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+in].n0ij[A_Z])*HYPER1[l*h_num+in].DgDq_n[A_Z])
				*(HYPER1[k*h_num+in].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+in].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+in].n0ij[A_Z]*nG[A_Z]);
			}
			rL[l*Nx+k]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*
			((HYPER[k].dPIdF[A_X][A_X]*HYPER1[k*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_X][A_Y]*HYPER1[k*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_X][A_Z]*HYPER1[k*h_num+k].n0ij[A_Z])*HYPER1[l*h_num+k].DgDq_n[A_X]
			+(HYPER[k].dPIdF[A_Y][A_X]*HYPER1[k*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Y][A_Y]*HYPER1[k*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Y][A_Z]*HYPER1[k*h_num+k].n0ij[A_Z])*HYPER1[l*h_num+k].DgDq_n[A_Y]
			+(HYPER[k].dPIdF[A_Z][A_X]*HYPER1[k*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Z][A_Y]*HYPER1[k*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Z][A_Z]*HYPER1[k*h_num+k].n0ij[A_Z])*HYPER1[l*h_num+k].DgDq_n[A_Z])
			*(HYPER1[k*h_num+k].n0ij[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].n0ij[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].n0ij[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);

			rL[(l+2*h_num)*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
			((HYPER[k].dPIdF[A_X][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_X][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_X][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*nG[A_X]
			+(HYPER[k].dPIdF[A_Y][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*nG[A_Y]
			+(HYPER[k].dPIdF[A_Z][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*nG[A_Z])
			*(HYPER1[k*h_num+k].n0ij[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].n0ij[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].n0ij[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);

			rL[l*Nx+k]+=0.25*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq_n[A_X]*HYPER1[l*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].DgDq_n[A_Y]*HYPER1[l*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].DgDq_n[A_Z]*HYPER1[l*h_num+k].DgDq_n[A_Z]);
			rL[(l+h_num)*Nx+k]+=0.25*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq_n[A_X]*HYPER1[l*h_num+k].DgDq[A_X]+HYPER1[k*h_num+k].DgDq_n[A_Y]*HYPER1[l*h_num+k].DgDq[A_Y]+HYPER1[k*h_num+k].DgDq_n[A_Z]*HYPER1[l*h_num+k].DgDq[A_Z]);

			rL[(l+2*h_num)*Nx+k]+=-0.25*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq_n[A_Z]*nG[A_Z]);
			rL[(l+3*h_num)*Nx+k]=-0.25*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq_n[A_Z]*nG[A_Z]);

			rL[l*Nx+k+h_num]+=0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*HYPER1[k*h_num+k].DgDq[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*HYPER1[k*h_num+k].DgDq[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*HYPER1[k*h_num+k].DgDq[A_Z]);
			rL[(l+h_num)*Nx+k+h_num]+=0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq[A_X]*HYPER1[k*h_num+k].DgDq[A_X]+HYPER1[l*h_num+k].DgDq[A_Y]*HYPER1[k*h_num+k].DgDq[A_Y]+HYPER1[l*h_num+k].DgDq[A_Z]*HYPER1[k*h_num+k].DgDq[A_Z]);

			rL[(l+2*h_num)*Nx+k+h_num]=-0.25*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq[A_Z]*nG[A_Z]);
			rL[(l+3*h_num)*Nx+k+h_num]=rL[(l+2*h_num)*Nx+k+h_num];

			rL[l*Nx+k+2*h_num]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*
			((HYPER[k].dPIdF[A_X][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_X][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_X][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*HYPER1[l*h_num+k].DgDq_n[A_X]
			+(HYPER[k].dPIdF[A_Y][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*HYPER1[l*h_num+k].DgDq_n[A_Y]
			+(HYPER[k].dPIdF[A_Z][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*HYPER1[l*h_num+k].DgDq_n[A_Z])
			*(HYPER1[k*h_num+k].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+k].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+k].n0ij[A_Z]*nG[A_Z]);

			rL[l*Nx+k+2*h_num]+=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*nG[A_Z]);

			rL[(l+h_num)*Nx+k+2*h_num]+=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq[A_Z]*nG[A_Z]);

			rL[(l+2*h_num)*Nx+k+2*h_num]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*
			((HYPER[k].dPIdF[A_X][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_X][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_X][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*nG[A_X]
			+(HYPER[k].dPIdF[A_Y][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Y][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Y][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*nG[A_Y]
			+(HYPER[k].dPIdF[A_Z][A_X]*HYPER1[l*h_num+k].n0ij[A_X]+HYPER[k].dPIdF[A_Z][A_Y]*HYPER1[l*h_num+k].n0ij[A_Y]+HYPER[k].dPIdF[A_Z][A_Z]*HYPER1[l*h_num+k].n0ij[A_Z])*nG[A_Z])
			*(HYPER1[k*h_num+k].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+k].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+k].n0ij[A_Z]*nG[A_Z]);

			rL[l*Nx+k+3*h_num]=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*nG[A_Z]);
			rL[(l+h_num)*Nx+k+3*h_num]=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq[A_Z]*nG[A_Z]);


			/////////////////dg, dgdt, dh, dhdtåvéZ1///////////////////
			dg[k][l]+=-0.5*Dt*Dt/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[l*h_num+k].DgDq_n[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[l*h_num+k].DgDq_n[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[l*h_num+k].DgDq_n[A_Z]);				
			d_dgdt[k][l+h_num]+=-0.5*Dt/mi*(HYPER1[k*h_num+k].DgDq[A_X]*HYPER1[l*h_num+k].DgDq[A_X]+HYPER1[k*h_num+k].DgDq[A_Y]*HYPER1[l*h_num+k].DgDq[A_Y]+HYPER1[k*h_num+k].DgDq[A_Z]*HYPER1[l*h_num+k].DgDq[A_Z]);				
			

			dg[k][l+2*h_num]=0.5*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq[A_Z]*nG[A_Z]);

			d_dgdt[k][l]=dg[k][l]/Dt;		
			d_dgdt[k][l+2*h_num]=dg[k][l+2*h_num]/Dt;
			d_dgdt[k][l+3*h_num]=dg[k][l+2*h_num]/Dt;

			dh[k][l]=0.5*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*nG[A_Z]);
			d_dhdt[k][l]=dh[k][l]/Dt;
			d_dhdt[k][l+h_num]=0.5*Dt/mi*(HYPER1[l*h_num+k].DgDq[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq[A_Z]*nG[A_Z]);
		}
		/////////////////rLåvéZ2///////////////////
		rL[(l+2*h_num)*Nx+l+2*h_num]+=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
		rL[(l+3*h_num)*Nx+l+2*h_num]=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);

		rL[(l+2*h_num)*Nx+l+3*h_num]=rL[(l+3*h_num)*Nx+l+2*h_num];
		rL[(l+3*h_num)*Nx+l+3*h_num]=rL[(l+3*h_num)*Nx+l+2*h_num];

		/////////////////dg, dgdt, dh, dhdtåvéZ2///////////////////
		dh[l][l+2*h_num]=-0.5*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
		d_dhdt[l][l+2*h_num]=dh[l][l+2*h_num]/Dt;
		d_dhdt[l][l+3*h_num]=dh[l][l+2*h_num]/Dt;
	}
}



void calc_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double nG[DIMENSION]={0,0,1};


	/////////////råvéZ
	for(int i=0;i<h_num;i++)
	{
		double p_r[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			p_r[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].h_lam)*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_r[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].h_lam)*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_r[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].h_lam)*HYPER1[jn*h_num+i].DgDq_n[A_Z];
		}
		p_r[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_r[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_r[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Z];

		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*( HYPER[i].p_n[A_X]+0.5*Dt*(p_r[A_X]+nG[A_X]*HYPER[i].h_mu) );
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*( HYPER[i].p_n[A_Y]+0.5*Dt*(p_r[A_Y]+nG[A_Y]*HYPER[i].h_mu) );
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*( HYPER[i].p_n[A_Z]+0.5*Dt*(p_r[A_Z]+nG[A_Z]*HYPER[i].h_mu) );
	}

	/////////////F, J, t_inverseÇÃåvéZ
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

		int Ni=HYPER[i].N;	
		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*h_num+inn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+inn].aiin[A_X],	HYPER1[i*h_num+inn].aiin[A_Y],	HYPER1[i*h_num+inn].aiin[A_Z]};
			
			fi[0][0]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_Z];
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


	/////////////DgDq, Stress, W, S, dSdcåvéZ
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double **in_Ci=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	in_Ci[D]=new double[DIMENSION];
	double in_Ci2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double first_term[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double first_term2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double S[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dSdc[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<h_num;j++)
	{	
		////////DgDqåvéZ
		int Ni=HYPER[j].N;
		for(int l=0;l<Ni;l++)
		{			
			int k=HYPER[j].NEI[l];
			HYPER1[k*h_num+j].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[j*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[j*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[j*h_num+k].n0ij[2]);
			HYPER1[k*h_num+j].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[j*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[j*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[j*h_num+k].n0ij[2]);
			HYPER1[k*h_num+j].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[j*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[j*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[j*h_num+k].n0ij[2]);
		}
		HYPER1[j*h_num+j].DgDq[A_X]=HYPER[j].J*(HYPER[j].t_inverse_Fi[A_X][0]*HYPER1[j*h_num+j].n0ij[0]+HYPER[j].t_inverse_Fi[A_X][1]*HYPER1[j*h_num+j].n0ij[1]+HYPER[j].t_inverse_Fi[A_X][2]*HYPER1[j*h_num+j].n0ij[2]);
		HYPER1[j*h_num+j].DgDq[A_Y]=HYPER[j].J*(HYPER[j].t_inverse_Fi[A_Y][0]*HYPER1[j*h_num+j].n0ij[0]+HYPER[j].t_inverse_Fi[A_Y][1]*HYPER1[j*h_num+j].n0ij[1]+HYPER[j].t_inverse_Fi[A_Y][2]*HYPER1[j*h_num+j].n0ij[2]);
		HYPER1[j*h_num+j].DgDq[A_Z]=HYPER[j].J*(HYPER[j].t_inverse_Fi[A_Z][0]*HYPER1[j*h_num+j].n0ij[0]+HYPER[j].t_inverse_Fi[A_Z][1]*HYPER1[j*h_num+j].n0ij[1]+HYPER[j].t_inverse_Fi[A_Z][2]*HYPER1[j*h_num+j].n0ij[2]);

		////////StressåvéZ
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
	
		HYPER[j].stress[0][0]=2/J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[j].stress[0][1]=2/J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[j].stress[0][2]=2/J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[j].stress[1][0]=2/J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[j].stress[1][1]=2/J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[j].stress[1][2]=2/J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[j].stress[2][0]=2/J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[j].stress[2][1]=2/J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[j].stress[2][2]=2/J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));


		////////WåvéZ
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

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
	
		HYPER[j].W=c10*(Ic-3)+c01*(IIc-3);

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

		first_term[0][0]=1-1/3*Ic*in_Ci[0][0];	first_term[0][1]=-1/3*Ic*in_Ci[0][1];	first_term[0][2]=-1/3*Ic*in_Ci[0][2];
		first_term[1][0]=-1/3*Ic*in_Ci[1][0];	first_term[1][1]=1-1/3*Ic*in_Ci[1][1];	first_term[1][2]=-1/3*Ic*in_Ci[1][0];
		first_term[2][0]=-1/3*Ic*in_Ci[1][0];	first_term[2][1]=-1/3*Ic*in_Ci[1][0];	first_term[2][2]=1-1/3*Ic*in_Ci[2][2];

		first_term2[0][0]=first_term[0][0]*first_term[0][0]+first_term[0][1]*first_term[1][0]+first_term[0][2]*first_term[2][0];
		first_term2[0][1]=first_term[0][0]*first_term[0][1]+first_term[0][1]*first_term[1][1]+first_term[0][2]*first_term[2][1];
		first_term2[0][2]=first_term[0][0]*first_term[0][2]+first_term[0][1]*first_term[1][2]+first_term[0][2]*first_term[2][2];
		first_term2[1][0]=first_term[1][0]*first_term[0][0]+first_term[1][1]*first_term[1][0]+first_term[1][2]*first_term[2][0];
		first_term2[1][1]=first_term[1][0]*first_term[0][1]+first_term[1][1]*first_term[1][1]+first_term[1][2]*first_term[2][1];
		first_term2[1][2]=first_term[1][0]*first_term[0][2]+first_term[1][1]*first_term[1][2]+first_term[1][2]*first_term[2][2];
		first_term2[2][0]=first_term[2][0]*first_term[0][0]+first_term[2][1]*first_term[1][0]+first_term[2][2]*first_term[2][0];
		first_term2[2][1]=first_term[2][0]*first_term[0][1]+first_term[2][1]*first_term[1][1]+first_term[2][2]*first_term[2][1];
		first_term2[2][2]=first_term[2][0]*first_term[0][2]+first_term[2][1]*first_term[1][2]+first_term[2][2]*first_term[2][2];

		if(J<0)
		{
			S[0][0]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3*trace_dC2*in_Ci[0][0]) );
			S[0][1]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3*trace_dC2*in_Ci[0][1]) );
			S[0][2]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3*trace_dC2*in_Ci[0][2]) );
			S[1][0]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3*trace_dC2*in_Ci[1][0]) );
			S[1][1]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3*trace_dC2*in_Ci[1][1]) );
			S[1][2]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3*trace_dC2*in_Ci[1][2]) );
			S[2][0]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3*trace_dC2*in_Ci[2][0]) );
			S[2][1]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3*trace_dC2*in_Ci[2][1]) );		
			S[2][2]=-2*1/pow(-J,2/3)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3*trace_dC2*in_Ci[2][2]) );

			dSdc[0][0]=-2*1/pow(-J,4/3)*( c01*first_term2[0][0]-(c10+c01*Ic)*(in_Ci[0][0]-4/9*Ic*in_Ci2[0][0])+2/3*c01*(1-2/3*trace_dC2*in_Ci2[0][0]) );
			dSdc[0][1]=-2*1/pow(-J,4/3)*( c01*first_term2[0][1]-(c10+c01*Ic)*(in_Ci[0][1]-4/9*Ic*in_Ci2[0][1])-4/9*c01*trace_dC2*in_Ci2[0][1] );
			dSdc[0][2]=-2*1/pow(-J,4/3)*( c01*first_term2[0][2]-(c10+c01*Ic)*(in_Ci[0][2]-4/9*Ic*in_Ci2[0][2])-4/9*c01*trace_dC2*in_Ci2[0][2] );
			dSdc[1][0]=-2*1/pow(-J,4/3)*( c01*first_term2[1][0]-(c10+c01*Ic)*(in_Ci[1][0]-4/9*Ic*in_Ci2[1][0])-4/9*c01*trace_dC2*in_Ci2[1][0] );
			dSdc[1][1]=-2*1/pow(-J,4/3)*( c01*first_term2[1][1]-(c10+c01*Ic)*(in_Ci[1][1]-4/9*Ic*in_Ci2[1][1])+2/3*c01*(1-2/3*trace_dC2*in_Ci2[1][1]) );
			dSdc[1][2]=-2*1/pow(-J,4/3)*( c01*first_term2[1][2]-(c10+c01*Ic)*(in_Ci[1][2]-4/9*Ic*in_Ci2[1][2])-4/9*c01*trace_dC2*in_Ci2[1][2] );
			dSdc[2][0]=-2*1/pow(-J,4/3)*( c01*first_term2[2][0]-(c10+c01*Ic)*(in_Ci[2][0]-4/9*Ic*in_Ci2[2][0])-4/9*c01*trace_dC2*in_Ci2[2][0] );
			dSdc[2][1]=-2*1/pow(-J,4/3)*( c01*first_term2[2][1]-(c10+c01*Ic)*(in_Ci[2][1]-4/9*Ic*in_Ci2[2][1])-4/9*c01*trace_dC2*in_Ci2[2][1] );		
			dSdc[2][2]=-2*1/pow(-J,4/3)*( c01*first_term2[2][2]-(c10+c01*Ic)*(in_Ci[2][2]-4/9*Ic*in_Ci2[2][2])+2/3*c01*(1-2/3*trace_dC2*in_Ci2[2][2]) );
		}
		else
		{
			S[0][0]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3*trace_dC2*in_Ci[0][0]) );
			S[0][1]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3*trace_dC2*in_Ci[0][1]) );
			S[0][2]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3*trace_dC2*in_Ci[0][2]) );
			S[1][0]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3*trace_dC2*in_Ci[1][0]) );
			S[1][1]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3*trace_dC2*in_Ci[1][1]) );
			S[1][2]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3*trace_dC2*in_Ci[1][2]) );
			S[2][0]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3*trace_dC2*in_Ci[2][0]) );
			S[2][1]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3*trace_dC2*in_Ci[2][1]) );		
			S[2][2]=2*1/pow(J,2/3)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3*trace_dC2*in_Ci[2][2]) );
	
			dSdc[0][0]=2*1/pow(J,4/3)*( c01*first_term2[0][0]-(c10+c01*Ic)*(in_Ci[0][0]-4/9*Ic*in_Ci2[0][0])+2/3*c01*(12/3*trace_dC2*in_Ci2[0][0]) );
			dSdc[0][1]=2*1/pow(J,4/3)*( c01*first_term2[0][1]-(c10+c01*Ic)*(in_Ci[0][1]-4/9*Ic*in_Ci2[0][1])-4/9*c01*trace_dC2*in_Ci2[0][1] );
			dSdc[0][2]=2*1/pow(J,4/3)*( c01*first_term2[0][2]-(c10+c01*Ic)*(in_Ci[0][2]-4/9*Ic*in_Ci2[0][2])-4/9*c01*trace_dC2*in_Ci2[0][2] );
			dSdc[1][0]=2*1/pow(J,4/3)*( c01*first_term2[1][0]-(c10+c01*Ic)*(in_Ci[1][0]-4/9*Ic*in_Ci2[1][0])-4/9*c01*trace_dC2*in_Ci2[1][0] );
			dSdc[1][1]=2*1/pow(J,4/3)*( c01*first_term2[1][1]-(c10+c01*Ic)*(in_Ci[1][1]-4/9*Ic*in_Ci2[1][1])+2/3*c01*(12/3*trace_dC2*in_Ci2[1][1]) );
			dSdc[1][2]=2*1/pow(J,4/3)*( c01*first_term2[1][2]-(c10+c01*Ic)*(in_Ci[1][2]-4/9*Ic*in_Ci2[1][2])-4/9*c01*trace_dC2*in_Ci2[1][2] );
			dSdc[2][0]=2*1/pow(J,4/3)*( c01*first_term2[2][0]-(c10+c01*Ic)*(in_Ci[2][0]-4/9*Ic*in_Ci2[2][0])-4/9*c01*trace_dC2*in_Ci2[2][0] );
			dSdc[2][1]=2*1/pow(J,4/3)*( c01*first_term2[2][1]-(c10+c01*Ic)*(in_Ci[2][1]-4/9*Ic*in_Ci2[2][1])-4/9*c01*trace_dC2*in_Ci2[2][1] );		
			dSdc[2][2]=2*1/pow(J,4/3)*( c01*first_term2[2][2]-(c10+c01*Ic)*(in_Ci[2][2]-4/9*Ic*in_Ci2[2][2])+2/3*c01*(12/3*trace_dC2*in_Ci2[2][2]) );
		}
		HYPER[j].dPIdF[0][0]=S[0][0]+2*dSdc[0][0]*HYPER[j].Fi[0][0]					 +dSdc[0][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[0][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		HYPER[j].dPIdF[0][1]=S[0][1]+dSdc[0][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2*dSdc[0][1]*HYPER[j].Fi[1][1]				  +dSdc[0][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		HYPER[j].dPIdF[0][2]=S[0][2]+dSdc[0][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[0][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2*dSdc[0][2]*HYPER[j].Fi[2][2];
		HYPER[j].dPIdF[1][0]=S[1][0]+2*dSdc[1][0]*HYPER[j].Fi[0][0]					 +dSdc[1][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[1][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		HYPER[j].dPIdF[1][1]=S[1][1]+dSdc[1][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2*dSdc[1][1]*HYPER[j].Fi[1][1]				  +dSdc[1][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		HYPER[j].dPIdF[1][2]=S[1][2]+dSdc[1][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[1][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2*dSdc[1][2]*HYPER[j].Fi[2][2];
		HYPER[j].dPIdF[2][0]=S[2][0]+2*dSdc[2][0]*HYPER[j].Fi[0][0]					 +dSdc[2][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[2][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		HYPER[j].dPIdF[2][1]=S[2][1]+dSdc[2][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2*dSdc[2][1]*HYPER[j].Fi[1][1]				  +dSdc[2][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		HYPER[j].dPIdF[2][2]=S[2][2]+dSdc[2][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[2][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2*dSdc[2][2]*HYPER[j].Fi[2][2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	in_Ci[D];
	delete[]	in_Ci;


	/////////////påvéZ
	for(int i=0;i<h_num;i++)
	{
		double p_r[DIMENSION]={0,0,0};
		double p_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			p_r[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].h_lam)*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_r[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].h_lam)*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_r[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].h_lam)*HYPER1[jn*h_num+i].DgDq_n[A_Z];
	
			p_p[A_X]+=(HYPER[jn].stress[A_X][A_X]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq[A_X]+HYPER[jn].stress[A_X][A_Y]*HYPER1[jn*h_num+i].DgDq[A_Y]+HYPER[jn].stress[A_X][A_Z]*HYPER1[jn*h_num+i].DgDq[A_Z];
			p_p[A_Y]+=HYPER[jn].stress[A_Y][A_X]*HYPER1[jn*h_num+i].DgDq[A_X]+(HYPER[jn].stress[A_Y][A_Y]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq[A_Y]+HYPER[jn].stress[A_Y][A_Z]*HYPER1[jn*h_num+i].DgDq[A_Z];
			p_p[A_Z]+=HYPER[jn].stress[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq[A_X]+HYPER[jn].stress[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq[A_Y]+(HYPER[jn].stress[A_Z][A_Z]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq[A_Z];
		}
		p_r[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_r[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_r[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Z];
	
		p_p[A_X]+=(HYPER[i].stress[A_X][A_X]-HYPER[i].lam)*HYPER1[i*h_num+i].DgDq[A_X]+HYPER[i].stress[A_X][A_Y]*HYPER1[i*h_num+i].DgDq[A_Y]+HYPER[i].stress[A_X][A_Z]*HYPER1[i*h_num+i].DgDq[A_Z];
		p_p[A_Y]+=HYPER[i].stress[A_Y][A_X]*HYPER1[i*h_num+i].DgDq[A_X]+(HYPER[i].stress[A_Y][A_Y]-HYPER[i].lam)*HYPER1[i*h_num+i].DgDq[A_Y]+HYPER[i].stress[A_Y][A_Z]*HYPER1[i*h_num+i].DgDq[A_Z];
		p_p[A_Z]+=HYPER[i].stress[A_Z][A_X]*HYPER1[i*h_num+i].DgDq[A_X]+HYPER[i].stress[A_Z][A_Y]*HYPER1[i*h_num+i].DgDq[A_Y]+(HYPER[i].stress[A_Z][A_Z]-HYPER[i].lam)*HYPER1[i*h_num+i].DgDq[A_Z];

		HYPER[i].p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*( p_r[A_X]+p_p[A_X]+nG[A_X]*(HYPER[i].h_mu+HYPER[i].mu) );
		HYPER[i].p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*( p_r[A_Y]+p_p[A_Y]+nG[A_Y]*(HYPER[i].h_mu+HYPER[i].mu) );
		HYPER[i].p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*( p_r[A_Z]+p_p[A_Z]+nG[A_Z]*(HYPER[i].h_mu+HYPER[i].mu) );
	}

}


void output_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,double *rT,double *dT,double *rL,double *dL,int Nx,int h_num,int count)
{
	double Ed=1.0E10;			
	
	stringstream ss0;
	ss0<<"./rT"<<count<<".csv";
	ofstream f_rt(ss0.str(), ios::trunc);

	stringstream ss1;
	ss1<<"./dT"<<count<<".csv";
	ofstream f_dt(ss1.str(), ios::trunc);

	stringstream ss2;
	ss2<<"./rL"<<count<<".csv";
	ofstream f_rL(ss2.str(), ios::trunc);

	stringstream ss3;
	ss3<<"./dL"<<count<<".csv";
	ofstream f_dL(ss3.str(), ios::trunc);

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nx;j++)
		{
			f_rt<<rT[i*Nx+j]<<",";
			f_rL<<rL[i*Nx+j]*Ed<<",";
		}
		f_rt<<endl;
		f_rL<<endl;
		f_dt<<dT[i]<<endl;
		f_dL<<dL[i]<<endl;
	}
	f_rt.close();
	f_rL.close();
	f_dt.close();
	f_dL.close();

	stringstream ss5;
	ss5<<"./p"<<count<<".csv";
	ofstream f_p(ss5.str(), ios::trunc);
	
	stringstream ss6;
	ss6<<"./q"<<count<<".csv";
	ofstream f_q(ss6.str(), ios::trunc);

	stringstream ss7;
	ss7<<"./h_lam"<<count<<".csv";
	ofstream f_h_lam(ss7.str(), ios::trunc);
	
	stringstream ss8;
	ss8<<"./lam"<<count<<".csv";
	ofstream f_lam(ss8.str(), ios::trunc);
	
	stringstream ss9;
	ss9<<"./h_mu"<<count<<".csv";
	ofstream f_h_mu(ss9.str(), ios::trunc);
	
	stringstream ss10;
	ss10<<"./mu"<<count<<".csv";
	ofstream f_mu(ss10.str(), ios::trunc);

	for(int i=0;i<h_num;i++)
	{				
		f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		f_h_lam<<HYPER[i].h_lam<<endl;
		f_lam<<HYPER[i].lam<<endl;
		f_h_mu<<HYPER[i].h_mu<<endl;
		f_mu<<HYPER[i].mu<<endl;

	}
	f_p.close();
	f_q.close();
	f_h_lam.close();
	f_lam.close();
	f_h_mu.close();
	f_mu.close();
}









