#include "stdafx.h"	

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, double V, double mi, double Dt, double E0,vector<double > NEIw);
void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,vector<double > NEIw);
void q_variables(mpsconfig &CON, vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double **dg, double *th_g, double *h, double **dh, double *th_h, double *d, int count,int count_min,int t,double E,double En,double E0,double mi, double V,vector<double > NEIw);

void calc_n(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void p_QP(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, int t,double V, double mi, double Dt, double E0, vector<double > NEIw);
void p_variables(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num,double Dt,double mi);
void p_nab(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dE, double *G, double *H,double Dt, double V, double mi,vector<double > NEIw);
void p_lap(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *rE, double **dG, double **dH, double Dt, double V, double mi,vector<double > NEIw);
void output_data_p(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double **dG, double *G, double *th_G, double **dH, double *H, double *th_H,double *rT, double *dT, double T, double *d, int count,int count_min,int t,double E,double En, double E0,vector<double > NEIw);


void q_QP_nw1(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, double V, double mi, double Dt, double E0,vector<double > NEIw);
void q_nab_lap_nw1(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,vector<double > NEIw);
void output_data_nw1(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double **dg, double *th_g, double *h, double **dh, double *th_h, double *d, int count,int count_min,int t,double E,double En,double E0,double mi, double V,vector<double > NEIw);

void p_QP_nw1(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, int t,double V, double mi, double Dt, double E0, vector<double > NEIw);
void p_nab_nw1(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dE, double *G, double *H,double Dt, double V, double mi,vector<double > NEIw);
void p_lap_nw1(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *rE, double **dG, double **dH, double Dt, double V, double mi,vector<double > NEIw);
void output_data_p_nw1(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double **dG, double *G, double *th_G, double **dH, double *H, double *th_H,double *rT, double *dT, double T, double *d, int count,int count_min,int t,double E,double En, double E0,vector<double > NEIw);





void calc_HYPER_QP_gh(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F,vector<double > NEIw)
{
	cout<<"QP_gh start-------------------------"<<endl;


	////////////íËã`///////////////
	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
	double nG[DIMENSION]={0,0,1};

	double E0=0;	
	for(int i=0;i<h_num;i++)
	{
		E0+=0.5/mi*(HYPER[i].p0[A_X]*HYPER[i].p0[A_X]+HYPER[i].p0[A_Y]*HYPER[i].p0[A_Y]+HYPER[i].p0[A_Z]*HYPER[i].p0[A_Z])+V*HYPER[i].W0;
	}
	if(t==1)
	{
		ofstream fe0("E0.csv");
		for(int i=0;i<h_num;i++)
		{
			fe0<<0.5/mi*(HYPER[i].p0[A_X]*HYPER[i].p0[A_X]+HYPER[i].p0[A_Y]*HYPER[i].p0[A_Y]+HYPER[i].p0[A_Z]*HYPER[i].p0[A_Z])+V*HYPER[i].W0<<","<<0.5/mi*(HYPER[i].p0[A_X]*HYPER[i].p0[A_X]+HYPER[i].p0[A_Y]*HYPER[i].p0[A_Y]+HYPER[i].p0[A_Z]*HYPER[i].p0[A_Z])<<","<<V*HYPER[i].W0<<endl;
		}
		fe0.close();
	}
	cout<<"E0="<<E0<<endl;

	ofstream fq0("qn_1_QP.csv");
	ofstream fp0("pn_1_QP.csv");
	ofstream fhp0("hpn_1_QP.csv");

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
	if(Nw==1)
	{
		q_QP_nw1(CON,PART,HYPER,HYPER1,t,V,mi,Dt,E0,NEIw);
		cout<<"qçXêV"<<endl;
		calc_n(CON,PART,HYPER,HYPER1);
		p_QP_nw1(PART,HYPER,HYPER1,t,V,mi,Dt,E0,NEIw);
	}
	else
	{
		q_QP(CON,PART,HYPER,HYPER1,t,V,mi,Dt,E0,NEIw);
		cout<<"qçXêV"<<endl;
		calc_n(CON,PART,HYPER,HYPER1);
		p_QP(PART,HYPER,HYPER1,t,V,mi,Dt,E0,NEIw);
	}

	for(int i=0;i<h_num;i++)	HYPER[i].lambda=HYPER[i].lam;
	cout<<"--------------------------OK"<<endl;



}

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, double V, double mi, double Dt, double E0,vector<double > NEIw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;
	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=Nw+h_num;

	double E=1;
	double E_min=1;
	double E_sum=1;
	double ep=1.e-5;
	double ep_min=1.e-5;

	double rg=1.;
	double rh=1.;
	//if(t==6)
	//{
	//	rg=1.;
	//	rh=1000.;
	//}

	double T=0.;
	double *dT=new double [Nx];	
	double *rT=new double [Nx*Nx];

	double En=0.;
	double *dE=new double [Nx];
	double *rE=new double [Nx*Nx];
	double *p_rE=new double [Nx*Nx];


	double *g=new double [h_num];
	double **dg=new double *[h_num];
	double *h=new double [Nw];
	double **dh=new double *[Nw];

	for(int i=0;i<h_num;i++)	dg[i]=new double [Nx];
	for(int i=0;i<Nw;i++)	dh[i]=new double [Nx];

	double *th_g=new double [h_num];
	double *th_h=new double [Nw];

	double *d=new double [Nx];
	double *B=new double [Nx*Nx];


	int fw=0;
	int iw=0;
	double p_p[DIMENSION]={0,0,0};
	double hp[DIMENSION]={0,0,0};
	double W=0.;
	double Dgki_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};
	double Dgji_n[DIMENSION]={0,0,0};
	double Dgli_n[DIMENSION]={0,0,0};
	double lam=0.;
	double mu=0.;
	double hi=0.;
	double d_sum=0.;

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].lam=1.;
		HYPER[i].mu=1.;
		g[i]=0.;
		th_g[i]=0.;
		for(int j=0;j<Nx;j++)	dg[i][j]=0.;
	}
	for(int i=0;i<Nw;i++)
	{
		h[i]=0.;
		th_h[i]=0.;
		for(int j=0;j<Nx;j++)	dh[i][j]=0.;
	}
	for(int i=0;i<Nx;i++)
	{
		dT[i]=0.;
		d[i]=0.;
		dE[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			rT[i*Nx+j]=0.;
			rE[i*Nx+j]=0.;
			p_rE[i*Nx+j]=0.;
			B[i*Nx+j]=0.;
		}
	}


	////////////dh, p_rEåvéZ///////////////
	for(int i=0;i<Nw;i++)
	{
		iw=NEIw[i];
		int Ni=HYPER[iw].N;
		Dgii_n[A_Z]=HYPER1[iw*h_num+iw].DgDq_n[A_Z];
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[iw].NEI[k];
			Dgki_n[A_Z]=HYPER1[kn*h_num+iw].DgDq_n[A_Z];
			dh[i][kn]+=0.5*Dt*Dt/mi*Dgki_n[A_Z];
			p_rE[(i+h_num)*Nx+kn]-=0.25*Dt*Dt/mi*Dgki_n[A_Z];
			p_rE[kn*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgki_n[A_Z];
		}
		dh[i][iw]+=0.5*Dt*Dt/mi*Dgii_n[A_Z];
		dh[i][i+h_num]-=0.5*Dt*Dt/mi;
		p_rE[(i+h_num)*Nx+iw]-=0.25*Dt*Dt/mi*Dgii_n[A_Z];
		p_rE[iw*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgii_n[A_Z];
		p_rE[(i+h_num)*Nx+i+h_num]+=0.25*Dt*Dt/mi;
	}
	for(int i=0;i<h_num;i++)
	{
		Dgii_n[A_X]=HYPER1[i*h_num+i].DgDq_n[A_X];	Dgii_n[A_Y]=HYPER1[i*h_num+i].DgDq_n[A_Y];	Dgii_n[A_Z]=HYPER1[i*h_num+i].DgDq_n[A_Z];

		int Ni=HYPER[i].N;
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[i].NEI[k];
			Dgki_n[A_X]=HYPER1[kn*h_num+i].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[kn*h_num+i].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[kn*h_num+i].DgDq_n[A_Z];
			for(int l=0;l<Ni;l++)
			{
				int ln=HYPER[i].NEI[l];
				Dgli_n[A_X]=HYPER1[ln*h_num+i].DgDq_n[A_X];	Dgli_n[A_Y]=HYPER1[ln*h_num+i].DgDq_n[A_Y];	Dgli_n[A_Z]=HYPER1[ln*h_num+i].DgDq_n[A_Z];
				p_rE[ln*Nx+kn]+=0.25*Dt*Dt/mi*(Dgki_n[A_X]*Dgli_n[A_X]+Dgki_n[A_Y]*Dgli_n[A_Y]+Dgki_n[A_Z]*Dgli_n[A_Z]);
			}
			p_rE[kn*Nx+i]+=0.25*Dt*Dt/mi*(Dgii_n[A_X]*Dgki_n[A_X]+Dgii_n[A_Y]*Dgki_n[A_Y]+Dgii_n[A_Z]*Dgki_n[A_Z]);
			p_rE[i*Nx+kn]+=0.25*Dt*Dt/mi*(Dgki_n[A_X]*Dgii_n[A_X]+Dgki_n[A_Y]*Dgii_n[A_Y]+Dgki_n[A_Z]*Dgii_n[A_Z]);
		}
		p_rE[i*Nx+i]+=0.25*Dt*Dt/mi*(Dgii_n[A_X]*Dgii_n[A_X]+Dgii_n[A_Y]*Dgii_n[A_Y]+Dgii_n[A_Z]*Dgii_n[A_Z]);
	}






	////////////QPñ@åvéZ///////////////
	while(E_min>ep_min)
	{
		count_min++;

		for(int i=0;i<h_num;i++)	HYPER[i].old_lam=HYPER[i].lam;
		for(int i=0;i<Nw;i++)
		{
			iw=NEIw[i];
			HYPER[iw].old_mu=HYPER[iw].mu;
		}


		E=1;
		count=0;
		while(E>ep)
		{
			count++;

			q_variables(CON,PART,HYPER,HYPER1);	//Ç±Ç±Ç‹Ç≈1219
			q_nab_lap(PART,HYPER,HYPER1,dE,rE,dg,Dt,V,mi,NEIw);

			En=0.;
			for(int i=0;i<h_num;i++)
			{
				W=HYPER[i].W;
				hp[A_X]=HYPER[i].half_p[A_X];	hp[A_Y]=HYPER[i].half_p[A_Y];	hp[A_Z]=HYPER[i].half_p[A_Z];
				En+=0.5/mi*(hp[A_X]*hp[A_X]+hp[A_Y]*hp[A_Y]+hp[A_Z]*hp[A_Z])+V*W;
			}

			T=(En-E0)*(En-E0);
			//cout<<"T="<<T<<", En"<<En<<endl;
			for(int k=0;k<Nx;k++)
			{	
				dT[k]=2.*dE[k]*(En-E0);
				for(int l=0;l<Nx;l++)
				{
					rE[l*Nx+k]+=p_rE[l*Nx+k];
					rT[l*Nx+k]=2.*rE[l*Nx+k]*(En-E0)+2.*dE[k]*dE[l];
				}
			}
			//for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			for(int i=0;i<h_num;i++)
			{	
				g[i]=V*(1.-HYPER[i].J);

				if(g[i]<0)	T+=0.5*rg*(-1.*g[i]+th_g[i])*(-1.*g[i]+th_g[i]);
				else
				{
					T+=0.5*rg*(g[i]+th_g[i])*(g[i]+th_g[i]);
				}
			}
			for(int i=0;i<Nw;i++)
			{	
				int iw=NEIw[i];
				h[i]=-1.*PART[iw].r[A_Z];

				//if(h[i]<0)	T+=0.5*rh*(-1.*h[i]+th_h[i])*(-1.*h[i]+th_h[i]);
				//else
				//{
				//	T+=0.5*rh*(h[i]+th_h[i])*(h[i]+th_h[i]);
				//}
				if(h[i]+th_h[i]>0)	T+=0.5*rh*(h[i]+th_h[i])*(h[i]+th_h[i]);
			}




			for(int k=0;k<Nx;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					if(g[i]<0)	dT[k]+=-1.*rg*dg[i][k]*(-1.*g[i]+th_g[i]);
					else
					{
						dT[k]+=rg*dg[i][k]*(g[i]+th_g[i]);
					}
					//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
				}
				for(int i=0;i<Nw;i++)
				{	
					//if(h[i]<0)	dT[k]+=-1.*rh*dh[i][k]*(-1.*h[i]+th_h[i]);
					//else
					//{
					//	dT[k]+=rh*dh[i][k]*(h[i]+th_h[i]);
					//}
					if(h[i]+th_h[i]>0)	dT[k]+=rh*dh[i][k]*(h[i]+th_h[i]);
				}
				for(int l=0;l<Nx;l++)
				{
					for(int i=0;i<h_num;i++)
					{	
						rT[l*Nx+k]+=rg*dg[i][k]*dg[i][l];
					}
					for(int i=0;i<Nw;i++)
					{	
						//rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						if(h[i]+th_h[i]>0)	rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						else if(h[i]+th_h[i]>-1.e-20 && dh[i][k]*dh[i][l]>0)
						{
							rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						}
					}
				}
			}

			d_sum=0.;
			for(int i=0;i<Nx;i++)
			{
				d[i]=dT[i];
				d_sum+=fabs(dT[i]);
				for(int j=0;j<Nx;j++)
				{
					B[i*Nx+j]=rT[i*Nx+j];
				}
			}
			if(d_sum<1.e-20)	break;

			gauss(B,d,Nx);

			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=d[i];
				E_sum+=fabs(d[i]);
			}
			for(int i=0;i<Nw;i++)
			{
				iw=NEIw[i];
				HYPER[iw].mu-=d[i+h_num];
				E_sum+=fabs(d[i+h_num]);
			}

			E=E_sum;




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
				if(fw==1)
				{
					mu=HYPER[i].mu;
					HYPER[i].half_p[A_Z]+=0.5*Dt*mu;
				}
				/////////////qåvéZ
				PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
				PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
				PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
			}	
			if(count==1||count%10000==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				output_data(PART,HYPER,HYPER1,T,dT,rT,g,dg,th_g,h,dh,th_h,d,count,count_min,t,E,En,E0,mi,V,NEIw);

			}

		}

		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=fabs(HYPER[i].old_lam-HYPER[i].lam);
		for(int i=0;i<Nw;i++)
		{
			iw=NEIw[i];
			E_sum+=fabs(HYPER[iw].old_mu-HYPER[iw].mu);
		}

		E_min=E_sum;
		//if(E_min<ep_min*1000)
		//{
		//	rg*=4;
		//	rh*=4;
		//}

		for(int i=0;i<h_num;i++)
		{
			if(g[i]<0)	th_g[i]+=-1.*g[i];
			else
			{
				th_g[i]+=g[i];
			}
		}
		for(int i=0;i<Nw;i++)
		{
			if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		}


		cout<<"Emin"<<count_min<<"="<<E_min<<endl;
		//if(count_min==1||count_min%100==0)
		{
			stringstream ss0;
			ss0<<"./E_min/E"<<t<<".csv"<<endl;
			stringstream ss1;
			ss1<<"./p/hp"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss2;
			ss2<<"./q/q"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss3;
			ss3<<"./lam/lam"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss4;
			ss4<<"./T/T"<<t<<".csv"<<endl;

			if(count_min==1)
			{
				ofstream fem(ss0.str(), ios::trunc);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::trunc);
				ft<<T<<","<<En<<","<<E0<<endl;
				ft.close();
			}
			else
			{
				ofstream fem(ss0.str(), ios::app);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::app);
				ft<<T<<","<<En<<endl;
				ft.close();
			}

			ofstream fhp(ss1.str());
			ofstream fq(ss2.str());
			ofstream fl(ss3.str());
			for(int i=0;i<h_num;i++)
			{
				fhp<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
				fq<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;			
				fl<<HYPER[i].lam<<","<<HYPER[i].mu<<endl;
			}
			fhp.close();
			fq.close();
			fl.close();

		}

		//NEIw.clear();
		//for(int i=0;i<h_num;i++)
		//{
		//	hi=-1.*PART[i].r[A_Z];
		//	if(hi>0)
		//	{
		//		NEIw.push_back(i);
		//		HYPER[i].fw=1;
		//	}
		//	else
		//	{
		//		HYPER[i].fw=0;
		//	}
		//}
		//Nw=NEIw.size();
		//if(Nw>0)
		//{
		//	cout<<"ê⁄êGó±éq ";
		//	for(int i=0;i<Nw;i++)	cout<<NEIw[i]<<", ";
		//	cout<<endl;
		//}


	}

	delete[]	dE;
	delete[]	dT;
	delete[]	rE;
	delete[]	rT;
	delete[]	g;
	delete[]	h;
	for(int i=0;i<h_num;i++)	delete[]	dg[i];
	for(int i=0;i<Nw;i++)	delete[]	dh[i];
	delete[]	dg;
	delete[]	dh;
	delete[]	th_h;
	delete[]	th_g;
	delete[]	d;
	delete[]	B;
	delete[]	p_rE;

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
		if(fw==1)
		{
			mu=HYPER[i].mu;
			HYPER[i].half_p[A_Z]+=0.5*Dt*mu;
		}
		/////////////qåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}	

}

void q_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
	int h_num=HYPER.size();

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();

	/////////////pn1_2åvéZ
	double p_p[DIMENSION]={0,0,0};
	double Dgji_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};
	double lam=0.;
	double mu=0.;
	int fw=0;
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
		if(fw==1)
		{
			mu=HYPER[i].mu;
			HYPER[i].half_p[A_Z]+=0.5*Dt*mu;
		}
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
			S[0][0]=-2./pow(-J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[0][0]) + c01*(Ic-dC[0][0]-2./3.*IIc*in_Ci[0][0]) );
			S[0][1]=-2./pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][1]) + c01*(  -dC[0][1]-2./3.*IIc*in_Ci[0][1]) );
			S[0][2]=-2./pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][2]) + c01*(  -dC[0][2]-2./3.*IIc*in_Ci[0][2]) );
			S[1][0]=-2./pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][0]) + c01*(  -dC[1][0]-2./3.*IIc*in_Ci[1][0]) );
			S[1][1]=-2./pow(-J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[1][1]) + c01*(Ic-dC[1][1]-2./3.*IIc*in_Ci[1][1]) );
			S[1][2]=-2./pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][2]) + c01*(  -dC[1][2]-2./3.*IIc*in_Ci[1][2]) );
			S[2][0]=-2./pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][0]) + c01*(  -dC[2][0]-2./3.*IIc*in_Ci[2][0]) );
			S[2][1]=-2./pow(-J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][1]) + c01*(  -dC[2][1]-2./3.*IIc*in_Ci[2][1]) );	
			S[2][2]=-2./pow(-J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[2][2]) + c01*(Ic-dC[2][2]-2./3.*IIc*in_Ci[2][2]) );

			dSdc[0][0]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][0]-in_Ci[0][0])+c01*(-2.*Ic*in_Ci[0][0]+5./3.*(1.+IIc*in_Ci2[0][0])) );
			dSdc[0][1]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][1]-in_Ci[0][1])+c01*(-2.*Ic*in_Ci[0][1]+5./3.*( +IIc*in_Ci2[0][1])) );
			dSdc[0][2]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][2]-in_Ci[0][2])+c01*(-2.*Ic*in_Ci[0][2]+5./3.*( +IIc*in_Ci2[0][2])) );
			dSdc[1][0]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][0]-in_Ci[1][0])+c01*(-2.*Ic*in_Ci[1][0]+5./3.*( +IIc*in_Ci2[1][0])) );
			dSdc[1][1]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][1]-in_Ci[1][1])+c01*(-2.*Ic*in_Ci[1][1]+5./3.*(1.+IIc*in_Ci2[1][1])) );
			dSdc[1][2]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][2]-in_Ci[1][2])+c01*(-2.*Ic*in_Ci[1][2]+5./3.*( +IIc*in_Ci2[1][2])) );
			dSdc[2][0]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][0]-in_Ci[2][0])+c01*(-2.*Ic*in_Ci[2][0]+5./3.*( +IIc*in_Ci2[2][0])) );
			dSdc[2][1]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][1]-in_Ci[2][1])+c01*(-2.*Ic*in_Ci[2][1]+5./3.*( +IIc*in_Ci2[2][1])) );
			dSdc[2][2]=-4./3./pow(-J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][2]-in_Ci[2][2])+c01*(-2.*Ic*in_Ci[2][2]+5./3.*(1.+IIc*in_Ci2[2][2])) );
		}
		else
		{
			S[0][0]=2./pow(J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[0][0]) + c01*(Ic-dC[0][0]-2./3.*IIc*in_Ci[0][0]) );
			S[0][1]=2./pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][1]) + c01*(  -dC[0][1]-2./3.*IIc*in_Ci[0][1]) );
			S[0][2]=2./pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[0][2]) + c01*(  -dC[0][2]-2./3.*IIc*in_Ci[0][2]) );
			S[1][0]=2./pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][0]) + c01*(  -dC[1][0]-2./3.*IIc*in_Ci[1][0]) );
			S[1][1]=2./pow(J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[1][1]) + c01*(Ic-dC[1][1]-2./3.*IIc*in_Ci[1][1]) );
			S[1][2]=2./pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[1][2]) + c01*(  -dC[1][2]-2./3.*IIc*in_Ci[1][2]) );
			S[2][0]=2./pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][0]) + c01*(  -dC[2][0]-2./3.*IIc*in_Ci[2][0]) );
			S[2][1]=2./pow(J,2./3.)*( c10*(  -1./3.*Ic*in_Ci[2][1]) + c01*(  -dC[2][1]-2./3.*IIc*in_Ci[2][1]) );	
			S[2][2]=2./pow(J,2./3.)*( c10*(1.-1./3.*Ic*in_Ci[2][2]) + c01*(Ic-dC[2][2]-2./3.*IIc*in_Ci[2][2]) );


			dSdc[0][0]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][0]-in_Ci[0][0])+c01*(-2.*Ic*in_Ci[0][0]+5./3.*(1.+IIc*in_Ci2[0][0])) );
			dSdc[0][1]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][1]-in_Ci[0][1])+c01*(-2.*Ic*in_Ci[0][1]+5./3.*( +IIc*in_Ci2[0][1])) );
			dSdc[0][2]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[0][2]-in_Ci[0][2])+c01*(-2.*Ic*in_Ci[0][2]+5./3.*( +IIc*in_Ci2[0][2])) );
			dSdc[1][0]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][0]-in_Ci[1][0])+c01*(-2.*Ic*in_Ci[1][0]+5./3.*( +IIc*in_Ci2[1][0])) );
			dSdc[1][1]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][1]-in_Ci[1][1])+c01*(-2.*Ic*in_Ci[1][1]+5./3.*(1.+IIc*in_Ci2[1][1])) );
			dSdc[1][2]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[1][2]-in_Ci[1][2])+c01*(-2.*Ic*in_Ci[1][2]+5./3.*( +IIc*in_Ci2[1][2])) );
			dSdc[2][0]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][0]-in_Ci[2][0])+c01*(-2.*Ic*in_Ci[2][0]+5./3.*( +IIc*in_Ci2[2][0])) );
			dSdc[2][1]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][1]-in_Ci[2][1])+c01*(-2.*Ic*in_Ci[2][1]+5./3.*( +IIc*in_Ci2[2][1])) );
			dSdc[2][2]=4./3./pow(J,4./3.)*( c10*(2./3.*Ic*in_Ci2[2][2]-in_Ci[2][2])+c01*(-2.*Ic*in_Ci[2][2]+5./3.*(1.+IIc*in_Ci2[2][2])) );
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

void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,vector<double > NEIw)
{

	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;
	/////////////////èâä˙âª///////////////////
	for(int k=0;k<Nx;k++)
	{
		dE[k]=0.;
		for(int i=0;i<h_num;i++)	dg[i][k]=0.;
		for(int i=0;i<Nx;i++)	rE[i*Nx+k]=0.;
	}


	/////////////////dE, dgåvéZ///////////////////
	//////////íËã`
	double Dgkk_n[DIMENSION]={0,0,0};
	double Dgki_n[DIMENSION]={0,0,0};

	double Dgik[DIMENSION]={0,0,0};
	double Dgik_n[DIMENSION]={0,0,0};
	double Dgkk[DIMENSION]={0,0,0};
	double	p_Eik[DIMENSION]={0,0,0};
	double	p_Ekk[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		/////////////dlam E1
		Dgkk_n[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];	Dgkk_n[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];	Dgkk_n[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];
		/////////////dlam E1
		Dgkk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	Dgkk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];
		p_Ekk[A_X]=HYPER[k].stress[A_X][0]*Dgkk[0]+HYPER[k].stress[A_X][1]*Dgkk[1]+HYPER[k].stress[A_X][2]*Dgkk[2];
		p_Ekk[A_Y]=HYPER[k].stress[A_Y][0]*Dgkk[0]+HYPER[k].stress[A_Y][1]*Dgkk[1]+HYPER[k].stress[A_Y][2]*Dgkk[2];
		p_Ekk[A_Z]=HYPER[k].stress[A_Z][0]*Dgkk[0]+HYPER[k].stress[A_Z][1]*Dgkk[1]+HYPER[k].stress[A_Z][2]*Dgkk[2];
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			/////////////dlam E1
			Dgki_n[A_X]=HYPER1[k*h_num+in].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[k*h_num+in].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[k*h_num+in].DgDq_n[A_Z];
			dE[k]-=0.5*Dt/mi*(Dgki_n[A_X]*HYPER[in].half_p[A_X]+Dgki_n[A_Y]*HYPER[in].half_p[A_Y]+Dgki_n[A_Z]*HYPER[in].half_p[A_Z]);
			/////////////dlam E2
			Dgik[A_X]=HYPER1[in*h_num+k].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+k].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+k].DgDq[A_Z];
			p_Eik[A_X]=HYPER[in].stress[A_X][0]*Dgik[0]+HYPER[in].stress[A_X][1]*Dgik[1]+HYPER[in].stress[A_X][2]*Dgik[2];
			p_Eik[A_Y]=HYPER[in].stress[A_Y][0]*Dgik[0]+HYPER[in].stress[A_Y][1]*Dgik[1]+HYPER[in].stress[A_Y][2]*Dgik[2];
			p_Eik[A_Z]=HYPER[in].stress[A_Z][0]*Dgik[0]+HYPER[in].stress[A_Z][1]*Dgik[1]+HYPER[in].stress[A_Z][2]*Dgik[2];

			Dgik_n[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];	Dgik_n[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];	Dgik_n[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];

				/////////////dlam E2
				dE[jn]+=0.5*Dt*Dt/mi*(p_Eik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+p_Eik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+p_Eik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
				/////////////dlam g
				dg[in][jn]-=0.5*Dt*Dt/mi*(Dgik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+Dgik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+Dgik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);

			}
			/////////////dlam E2
			dE[k]+=0.5*Dt*Dt/mi*(p_Eik[A_X]*Dgkk_n[A_X]+p_Eik[A_Y]*Dgkk_n[A_Y]+p_Eik[A_Z]*Dgkk_n[A_Z]);
			dE[in]+=0.5*Dt*Dt/mi*(p_Ekk[A_X]*Dgik_n[A_X]+p_Ekk[A_Y]*Dgik_n[A_Y]+p_Ekk[A_Z]*Dgik_n[A_Z]);
			/////////////dlam g
			dg[in][k]-=0.5*Dt*Dt/mi*(Dgik[A_X]*Dgkk_n[A_X]+Dgik[A_Y]*Dgkk_n[A_Y]+Dgik[A_Z]*Dgkk_n[A_Z]);
			dg[k][in]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgik_n[A_X]+Dgkk[A_Y]*Dgik_n[A_Y]+Dgkk[A_Z]*Dgik_n[A_Z]);
		}
		/////////////dlam E1
		dE[k]-=0.5*Dt/mi*(Dgkk_n[A_X]*HYPER[k].half_p[A_X]+Dgkk_n[A_Y]*HYPER[k].half_p[A_Y]+Dgkk_n[A_Z]*HYPER[k].half_p[A_Z]);
		/////////////dlam E2
		dE[k]+=0.5*Dt*Dt/mi*(p_Ekk[A_X]*Dgkk_n[A_X]+p_Ekk[A_Y]*Dgkk_n[A_Y]+p_Ekk[A_Z]*Dgkk_n[A_Z]);
		/////////////dlam g
		dg[k][k]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgkk_n[A_X]+Dgkk[A_Y]*Dgkk_n[A_Y]+Dgkk[A_Z]*Dgkk_n[A_Z]);
	}


	//////////íËã`
	for(int k=0;k<Nw;k++)
	{
		int kw=NEIw[k];

		/////////////dmu E2
		Dgkk[A_X]=HYPER1[kw*h_num+kw].DgDq[A_X];	Dgkk[A_Y]=HYPER1[kw*h_num+kw].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[kw*h_num+kw].DgDq[A_Z];
		p_Ekk[A_Z]=HYPER[kw].stress[A_Z][0]*Dgkk[0]+HYPER[kw].stress[A_Z][1]*Dgkk[1]+HYPER[kw].stress[A_Z][2]*Dgkk[2];

		int Nk=HYPER[kw].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[kw].NEI[i];

			/////////////dmu E2
			Dgik[A_X]=HYPER1[in*h_num+kw].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+kw].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+kw].DgDq[A_Z];

			p_Eik[A_Z]=HYPER[in].stress[A_Z][0]*Dgik[0]+HYPER[in].stress[A_Z][1]*Dgik[1]+HYPER[in].stress[A_Z][2]*Dgik[2];
			dE[k+h_num]-=0.5*Dt*Dt/mi*p_Eik[A_Z];

			/////////////dmu g
			dg[in][k+h_num]+=0.5*Dt*Dt/mi*Dgik[A_Z];
		}
		/////////////dmu E2
		dE[k+h_num]-=0.5*Dt*Dt/mi*p_Ekk[A_Z];

		/////////////dmu E1
		dE[k+h_num]+=0.5*Dt/mi*HYPER[kw].half_p[A_Z];

		/////////////dmu g
		dg[kw][k+h_num]+=0.5*Dt*Dt/mi*Dgkk[A_Z];
	}



	/////////////////rEåvéZ///////////////////
	//////////íËã`
	/////////////rlam E1
	double Dgli_n[DIMENSION]={0,0,0};
	double Dpki[DIMENSION]={0,0,0};
	double Dpii[DIMENSION]={0,0,0};
	double *DpiDlam=new double [h_num*h_num];
	double *DpiDmu_x=new double [h_num*h_num];
	double *DpiDmu_y=new double [h_num*h_num];
	double *DpiDmu_z=new double [h_num*h_num];
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			DpiDlam[i*h_num+j]=0.;
			DpiDmu_x[i*h_num+j]=0.;
			DpiDmu_y[i*h_num+j]=0.;
			DpiDmu_z[i*h_num+j]=0.;

		}
	}
	/////////////rlam E3
	double Dgli[DIMENSION]={0,0,0};
	/////////////rlam E4
	double Dgki[DIMENSION]={0,0,0};
	double Dgii[DIMENSION]={0,0,0};
	/////////////dpidlam
	double Dgii_n[DIMENSION]={0,0,0};
	/////////////dpidmu
	double Dpik[DIMENSION]={0,0,0};
	double n0li[DIMENSION]={0,0,0};
	double n0ii[DIMENSION]={0,0,0};
	double n0ki[DIMENSION]={0,0,0};

	for(int i=0;i<h_num;i++)
	{
		/////////////dpidlama
		Dpii[A_X]=HYPER1[i*h_num+i].DpiDq[A_X];	Dpii[A_Y]=HYPER1[i*h_num+i].DpiDq[A_Y];	Dpii[A_Z]=HYPER1[i*h_num+i].DpiDq[A_Z];
		/////////////dpidmu
		n0ii[A_X]=HYPER1[i*h_num+i].n0ij[A_X];	n0li[A_Y]=HYPER1[i*h_num+i].n0ij[A_Y];	n0li[A_Z]=HYPER1[i*h_num+i].n0ij[A_Z];
		int Ni=HYPER[i].N;
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[i].NEI[k];
			/////////////dpidmu
			Dgki_n[A_X]=HYPER1[kn*h_num+i].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[kn*h_num+i].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[kn*h_num+i].DgDq_n[A_Z];
			Dpki[A_X]=HYPER1[kn*h_num+i].DpiDq[A_X];	Dpki[A_Y]=HYPER1[kn*h_num+i].DpiDq[A_Y];	Dpki[A_Z]=HYPER1[kn*h_num+i].DpiDq[A_Z];
			Dpik[A_X]=HYPER1[i*h_num+kn].DpiDq[A_X];	Dpik[A_Y]=HYPER1[i*h_num+kn].DpiDq[A_Y];	Dpik[A_Z]=HYPER1[i*h_num+kn].DpiDq[A_Z];
			n0ki[A_X]=HYPER1[kn*h_num+i].n0ij[A_X];	n0ki[A_Y]=HYPER1[kn*h_num+i].n0ij[A_Y];	n0ki[A_Z]=HYPER1[kn*h_num+i].n0ij[A_Z];
			for(int l=0;l<Ni;l++)
			{
				int ln=HYPER[i].NEI[l];
				/////////////dpidlam
				Dgli_n[A_X]=HYPER1[ln*h_num+i].DgDq_n[A_X];	Dgli_n[A_Y]=HYPER1[ln*h_num+i].DgDq_n[A_Y];	Dgli_n[A_Z]=HYPER1[ln*h_num+i].DgDq_n[A_Z];
				DpiDlam[kn*h_num+ln]+=Dpki[A_X]*Dgli_n[A_X]+Dpki[A_Y]*Dgli_n[A_Y]+Dpki[A_Z]*Dgli_n[A_Z];
				/////////////dpidmu
				n0li[A_X]=HYPER1[ln*h_num+i].n0ij[A_X];	n0li[A_Y]=HYPER1[ln*h_num+i].n0ij[A_Y];	n0li[A_Z]=HYPER1[ln*h_num+i].n0ij[A_Z];
				DpiDmu_x[kn*h_num+ln]+=Dpik[A_Z]*n0li[A_X];	DpiDmu_y[kn*h_num+ln]+=Dpik[A_Z]*n0li[A_Y];	DpiDmu_z[kn*h_num+ln]+=Dpik[A_Z]*n0li[A_Z];

			}
			/////////////dpidlama
			DpiDlam[kn*h_num+i]+=Dpki[A_X]*Dgii_n[A_X]+Dpki[A_Y]*Dgii_n[A_Y]+Dpki[A_Z]*Dgii_n[A_Z];
			DpiDlam[i*h_num+kn]+=Dpii[A_X]*Dgki_n[A_X]+Dpii[A_Y]*Dgki_n[A_Y]+Dpii[A_Z]*Dgki_n[A_Z];
			/////////////dpidmu
			DpiDmu_x[kn*h_num+i]+=Dpik[A_Z]*n0ii[A_X];	DpiDmu_y[kn*h_num+i]+=Dpik[A_Z]*n0ii[A_Y];	DpiDmu_z[kn*h_num+i]+=Dpik[A_Z]*n0ii[A_Z];
			DpiDmu_x[i*h_num+kn]+=Dpii[A_Z]*n0ki[A_X];	DpiDmu_y[i*h_num+kn]+=Dpii[A_Z]*n0ki[A_Y];	DpiDmu_z[i*h_num+kn]+=Dpii[A_Z]*n0ki[A_Z];

		}
		/////////////dpidlama
		DpiDlam[i*h_num+i]+=Dpii[A_X]*Dgii_n[A_X]+Dpii[A_Y]*Dgii_n[A_Y]+Dpii[A_Z]*Dgii_n[A_Z];
		/////////////dpidmu
		DpiDmu_x[i*h_num+i]+=Dpii[A_Z]*n0ii[A_X];	DpiDmu_y[i*h_num+i]+=Dpii[A_Z]*n0ii[A_Y];	DpiDmu_z[i*h_num+i]+=Dpii[A_Z]*n0ii[A_Z];
	}

	/////////////rlam E2
	double Dgkm_n[DIMENSION]={0,0,0};
	double n0mi[DIMENSION]={0,0,0};
	//double n0ii[DIMENSION]={0,0,0};
	//double n0li[DIMENSION]={0,0,0};
	double Dgkl_n[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		for(int i=0;i<h_num;i++)
		{
			int Ni=HYPER[i].N;
			n0ii[A_X]=HYPER1[i*h_num+i].n0ij[A_X];	n0ii[A_Y]=HYPER1[i*h_num+i].n0ij[A_Y];	n0ii[A_Z]=HYPER1[i*h_num+i].n0ij[A_Z];
			Dgki_n[A_X]=HYPER1[k*h_num+i].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[k*h_num+i].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[k*h_num+i].DgDq_n[A_Z];
			for(int l=0;l<Ni;l++)
			{
				int ln=HYPER[i].NEI[l];
				n0li[A_X]=HYPER1[ln*h_num+i].n0ij[A_X];	n0li[A_Y]=HYPER1[ln*h_num+i].n0ij[A_Y];	n0li[A_Z]=HYPER1[ln*h_num+i].n0ij[A_Z];
				Dgkl_n[A_X]=HYPER1[k*h_num+ln].DgDq_n[A_X];	Dgkl_n[A_Y]=HYPER1[k*h_num+ln].DgDq_n[A_Y];	Dgkl_n[A_Z]=HYPER1[k*h_num+ln].DgDq_n[A_Z];

				for(int m=0;m<Ni;m++)
				{
					int mn=HYPER[i].NEI[m];
					n0mi[A_X]=HYPER1[mn*h_num+i].n0ij[A_X];	n0mi[A_Y]=HYPER1[mn*h_num+i].n0ij[A_Y];	n0mi[A_Z]=HYPER1[mn*h_num+i].n0ij[A_Z];
					Dgkm_n[A_X]=HYPER1[k*h_num+mn].DgDq_n[A_X];	Dgkm_n[A_Y]=HYPER1[k*h_num+mn].DgDq_n[A_Y];	Dgkm_n[A_Z]=HYPER1[k*h_num+mn].DgDq_n[A_Z];

					rE[ln*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+ln]*(n0mi[A_X]*Dgkm_n[A_X]+n0mi[A_Y]*Dgkm_n[A_Y]+n0mi[A_Z]*Dgkm_n[A_Z]);
				}
				rE[ln*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+ln]*(n0ii[A_X]*Dgki_n[A_X]+n0ii[A_Y]*Dgki_n[A_Y]+n0ii[A_Z]*Dgki_n[A_Z]);
				rE[i*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+i]*(n0li[A_X]*Dgkl_n[A_X]+n0li[A_Y]*Dgkl_n[A_Y]+n0li[A_Z]*Dgkl_n[A_Z]);
			}
			rE[i*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+i]*(n0ii[A_X]*Dgki_n[A_X]+n0ii[A_Y]*Dgki_n[A_Y]+n0ii[A_Z]*Dgki_n[A_Z]);
		}
	}

	/////////////rmu E2
	for(int l=0;l<Nw;l++)
	{
		int lw=NEIw[l];

		for(int k=0;k<Nw;k++)
		{
			int kw=NEIw[k];
			/////////////rmu E2
			rE[(l+h_num)*Nx+(k+h_num)]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDmu_z[lw*h_num+kw];
		}
	}


	/////////////rmulam E2
	double Dplm[DIMENSION]={0,0,0};
	for(int l=0;l<Nw;l++)
	{
		int lw=NEIw[l];

		int Nl=HYPER[lw].N;
		for(int k=0;k<Nl;k++)
		{
			int kn=HYPER[lw].NEI[k];;
			for(int m=0;m<h_num;m++)
			{
				/////////////rmulam E2
				Dgkm_n[A_X]=HYPER1[kn*h_num+m].DgDq_n[A_X];	Dgkm_n[A_Y]=HYPER1[kn*h_num+m].DgDq_n[A_Y];	Dgkm_n[A_Z]=HYPER1[kn*h_num+m].DgDq_n[A_Z];
				Dplm[A_X]=DpiDmu_x[lw*h_num+m];	Dplm[A_Y]=DpiDmu_y[lw*h_num+m];	Dplm[A_Z]=DpiDmu_z[lw*h_num+m];
				rE[(l+h_num)*Nx+kn]-=0.25*Dt*Dt*Dt*Dt/mi/mi/V*(Dplm[A_X]*Dgkm_n[A_X]+Dplm[A_Y]*Dgkm_n[A_Y]+Dplm[A_Z]*Dgkm_n[A_Z]);
			}
		}
	}


	for(int k=0;k<Nw;k++)
	{
		int kw=NEIw[k];
		int Nk=HYPER[kw].N;

		for(int l=0;l<h_num;l++)
		{
			for(int i=0;i<h_num;i++)
			{
				/////////////rlammu E2
				n0ki[A_X]=HYPER1[kw*h_num+i].n0ij[A_X];	n0ki[A_Y]=HYPER1[kw*h_num+i].n0ij[A_Y];	n0ki[A_Z]=HYPER1[kw*h_num+i].n0ij[A_Z];
				rE[l*h_num+k+h_num]-=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+l]*n0ki[A_Z];
			}
		}

	}



	delete[]	DpiDlam;
	delete[]	DpiDmu_x;
	delete[]	DpiDmu_y;
	delete[]	DpiDmu_z;

}


void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double **dg, double *th_g, double *h, double **dh, double *th_h, double *d, int count,int count_min,int t,double E,double En,double E0,double mi, double V,vector<double > NEIw)
{
	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;

	stringstream ss0;
	ss0<<"./T/T_"<<t<<"_"<<count_min<<".csv";
	stringstream ss1;
	ss1<<"./T/dT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./T/rT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./lam/lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./g/g_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss7;
	ss7<<"./g/dg_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./p/hp_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./q/q_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./T/d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./E/E_"<<t<<"_"<<count_min<<".csv";
	//stringstream ss12;
	//ss12<<"./h/h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss13;
	//ss13<<"./g/J_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss14;
	//ss14<<"./g/DgDq_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss15;
	//ss15<<"./g/stress_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss16;
	ss16<<"./T/En_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss17;
	ss17<<"./h/h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss18;
	ss18<<"./h/dh_"<<t<<".csv";
	stringstream ss19;
	ss19<<"./mu/mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";


	ofstream f_dt(ss1.str());
	ofstream f_rt(ss2.str());
	ofstream f_lam(ss3.str());
	ofstream f_g(ss6.str());
	ofstream f_dg(ss7.str());
	ofstream f_p(ss8.str());
	ofstream f_q(ss9.str());
	ofstream f_d(ss10.str());
	//ofstream f_h(ss12.str());
	//ofstream f_J(ss13.str());
	//ofstream f_dgdq(ss14.str());
	//ofstream f_s(ss15.str());
	ofstream f_en(ss16.str());
	ofstream f_h(ss17.str());

	ofstream f_mu(ss19.str());

	if(count_min==1&&count==1)
	{
		ofstream f_dh(ss18.str());
		for(int i=0;i<Nw;i++)
		{
			for(int k=0;k<Nx;k++)
			{
				f_dh<<dh[i][k]<<",";
			}
			f_dh<<endl;
		}
		f_dh.close();
	}	
	if(count==1)
	{
		ofstream f_T(ss0.str(), ios::trunc);
		f_T<<count<<","<<T<<","<<En<<endl;
		f_T.close();
		ofstream f_E(ss11.str(), ios::trunc);
		f_E<<count<<","<<E<<endl;
		f_E.close();
	}	
	else
	{
		ofstream f_T(ss0.str(), ios::app);
		f_T<<count<<","<<T<<","<<En<<endl;
		f_T.close();
		ofstream f_E(ss11.str(), ios::app);
		f_E<<count<<","<<E<<endl;
		f_E.close();
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
	for(int i=0;i<h_num;i++)
	{
		f_lam<<HYPER[i].lam<<endl;
		f_mu<<HYPER[i].mu<<endl;

		for(int k=0;k<Nx;k++)
		{
			f_dg<<dg[i][k]<<",";
		}
		f_dg<<endl;

		f_g<<g[i]<<","<<th_g[i]<<endl;
		f_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
		f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		//f_J<<HYPER[i].J<<endl;

		//f_dgdq<<i<<endl;
		//int Ni=HYPER[i].N;
		//for(int j=0;j<Ni;j++)
		//{
		//	int jn=HYPER[i].NEI[j];
		//	f_dgdq<<","<<jn<<","<<HYPER1[i*h_num+jn].DgDq[A_X]<<","<<HYPER1[i*h_num+jn].DgDq[A_Y]<<","<<HYPER1[i*h_num+jn].DgDq[A_Z]<<endl;
		//}

		////f_h<<h[i]<<","<<th_h[i]<<endl;
		//f_s<<i<<","<<HYPER[i].stress[A_X][0]<<","<<HYPER[i].stress[A_X][1]<<","<<HYPER[i].stress[A_X][2]<<endl;
		//f_s<<","<<HYPER[i].stress[A_Y][0]<<","<<HYPER[i].stress[A_Y][1]<<","<<HYPER[i].stress[A_Y][2]<<endl;
		//f_s<<","<<HYPER[i].stress[A_Z][0]<<","<<HYPER[i].stress[A_Z][1]<<","<<HYPER[i].stress[A_Z][2]<<endl;
		f_en<<0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])+V*HYPER[i].W<<","<<0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])<<","<<V*HYPER[i].W<<endl;
	}
	for(int i=0;i<Nw;i++)
	{
		f_h<<h[i]<<","<<th_h[i]<<endl;
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
	//f_J.close();
	//f_dgdq.close();
	//f_s.close();
	f_en.close();
	f_h.close();
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


	double trace_dC=0.;
	double trace_dC2=0.;

	double Ic=0.;
	double IIc=0.;

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
	}
	//cout<<"----------OK"<<endl;
}


void p_QP(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, int t,double V, double mi, double Dt, double E0, vector<double > NEIw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;
	int Nw=NEIw.size();
	int h_num=HYPER.size();
	int Nx=Nw+h_num;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1.e-5;
	double ep_min=1.e-5;
	double d_sum=0.;

	double p_p[DIMENSION]={0,0,0};
	double p[DIMENSION]={0,0,0};
	double W=0.;
	double g=0.;
	double rg=1.0;
	double rh=1.0;

	double En=0.;
	double *dE=new double [Nx];	
	double *rE=new double [Nx*Nx];

	double T=0.;
	double *dT=new double [Nx];	
	double *rT=new double [Nx*Nx];

	double *d=new double [Nx];
	double *B=new double [Nx*Nx];

	double *G=new double [h_num];
	double *th_G=new double [h_num];
	double **dG=new double *[h_num];

	double *H=new double [Nw];
	double *th_H=new double [Nw];
	double **dH=new double *[Nw];

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)	dG[i]=new double [Nx];
	for(int i=0;i<Nw;i++)	dH[i]=new double [Nx];

	for(int i=0;i<h_num;i++)
	{
		G[i]=0.;
		th_G[i]=0.;
		for(int j=0;j<Nx;j++)	dG[i][j]=0.;
	}
	for(int i=0;i<Nw;i++)
	{
		H[i]=0.;
		th_H[i]=0.;
		for(int j=0;j<Nx;j++)	dH[i][j]=0.;
	}

	for(int i=0;i<Nx;i++)
	{
		dE[i]=0.;
		dT[i]=0.;
		d[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			rE[i*Nx+j]=0.;
			rT[i*Nx+j]=0.;
			B[i*Nx+j]=0.;
		}
	}

	p_lap(HYPER,HYPER1,rE,dG,dH,Dt,V,mi,NEIw);



	while(E_min>ep_min)
	{
		count_min++;

		for(int i=0;i<h_num;i++)	HYPER[i].old_lam=HYPER[i].lam;
		for(int i=0;i<Nw;i++)
		{
			int iw=NEIw[i];
			HYPER[iw].old_mu=HYPER[iw].mu;		
		}

		E=1;
		count=0;
		while(E>ep)
		{
			count++;
			p_variables(HYPER,HYPER1,h_num,Dt,mi);
			p_nab(HYPER,HYPER1,dE,G,H,Dt,V,mi,NEIw);

			En=0;
			for(int i=0;i<h_num;i++)	
			{

				p[A_X]=HYPER[i].p[A_X];	p[A_Y]=HYPER[i].p[A_Y];	p[A_Z]=HYPER[i].p[A_Z];
				W=HYPER[i].W;
				//En+=0.5/mi*(p[A_X]*p[A_X]+p[A_Y]*p[A_Y]+p[A_Z]*p[A_Z])+V*W+HYPER[i].mu*g;
				En+=0.5/mi*(p[A_X]*p[A_X]+p[A_Y]*p[A_Y]+p[A_Z]*p[A_Z])+V*W;
			}

			T=0.;
			for(int k=0;k<Nx;k++)
			{	
				dT[k]=0.;
				for(int l=0;l<Nx;l++)	rT[l*Nx+k]=0.;
			}


			T+=(En-E0)*(En-E0);
			for(int k=0;k<Nx;k++)
			{	
				dT[k]+=2.*dE[k]*(En-E0);
				for(int l=0;l<Nx;l++)	rT[l*Nx+k]+=2.*rE[l*Nx+k]*(En-E0)+2*dE[k]*dE[l];
			}
			for(int i=0;i<h_num;i++)
			{	
				if(G[i]<0)	T+=0.5*rg*(-1.*G[i]+th_G[i])*(-1.*G[i]+th_G[i]);
				else
				{
					T+=0.5*rg*(G[i]+th_G[i])*(G[i]+th_G[i]);
				}
			}
			for(int i=0;i<Nw;i++)	if(H[i]+th_H[i]>0)	T+=0.5*rh*(H[i]+th_H[i])*(H[i]+th_H[i]);

			for(int k=0;k<Nx;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					if(G[i]<0)	dT[k]+=-1.*rg*dG[i][k]*(-1.*G[i]+th_G[i]);
					else
					{
						dT[k]+=rg*dG[i][k]*(G[i]+th_G[i]);
					}
				}
				for(int i=0;i<Nw;i++)
				{	
					if(H[i]+th_H[i]>0)	dT[k]+=rh*dH[i][k]*(H[i]+th_H[i]);
				}


				for(int l=0;l<Nx;l++)
				{
					for(int i=0;i<h_num;i++)
					{	
						//rT[l*h_num+k]+=r*2.*dG[i*h_num+k]*dG[i*h_num+l]*(3*G[i]*G[i]+th_G[i]);
						rT[l*Nx+k]+=rg*dG[i][k]*dG[i][l];
						//rT[l*h_num+k]+=r*dG[i*h_num+k]*dG[i*h_num+l];
					}
					for(int i=0;i<Nw;i++)
					{	
						//rT[l*h_num+k]+=r*2.*dG[i*h_num+k]*dG[i*h_num+l]*(3*G[i]*G[i]+th_G[i]);
						if(H[i]+th_H[i]>0)	rT[l*Nx+k]+=rh*dH[i][k]*dH[i][l];
						else if (H[i]+th_H[i]>-1.e-20 && dH[i][k]*dH[i][l]>0)	rT[l*Nx+k]+=rh*dH[i][k]*dH[i][l];
						//rT[l*h_num+k]+=r*dG[i*h_num+k]*dG[i*h_num+l];
					}
				}
			}

			d_sum=0.;
			for(int i=0;i<Nx;i++)
			{
				d[i]=dT[i];
				d_sum+=fabs(dT[i]);
				for(int j=0;j<Nx;j++)	B[i*Nx+j]=rT[i*Nx+j];
			}
			if(d_sum<1.e-20)	break;
			gauss(B,d,Nx);

			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=d[i];
				E_sum+=fabs(d[i]);
			}
			for(int i=0;i<Nw;i++)
			{
				int iw=NEIw[i];
				HYPER[iw].mu-=d[i+h_num];
				E_sum+=fabs(d[i+h_num]);
			}
			E=E_sum;
			//cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;

			p_variables(HYPER,HYPER1,h_num,Dt,mi);
			if(count==1||count%10000==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				output_data_p(HYPER,HYPER1,dG,G,th_G,dH,H,th_H,rT,dT,T,d,count,count_min,t,E,En,E0,NEIw);
			}

		}
		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=fabs(HYPER[i].old_lam-HYPER[i].lam);
		for(int i=0;i<Nw;i++)
		{
			int iw=NEIw[i];
			E_sum+=fabs(HYPER[iw].old_mu-HYPER[iw].mu);
		}
		E_min=E_sum;

		//if(E_min<ep_min*1000)
		//{
		//	rg*=4;
		//	rh*=4;
		//}


		for(int i=0;i<h_num;i++)
		{
			if(G[i]<0.)	th_G[i]+=-1.*G[i];
			else
			{
				th_G[i]+=G[i];
			}
			//th_g[i]+=g[i];
			//if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		}
		for(int i=0;i<Nw;i++)	if(H[i]+th_H[i]>0.)	th_H[i]+=H[i];

		//////ÉfÅ[É^èoóÕ
		//if(count_min==1||count_min%100==0)
		{
			cout<<"Emin"<<count_min<<"="<<E_min<<endl;
			stringstream ss0;
			ss0<<"./E_min/E"<<t<<".csv"<<endl;
			stringstream ss1;
			ss1<<"./p/p"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss3;
			ss3<<"./mu/mu"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss4;
			ss4<<"./p_T/T"<<t<<".csv"<<endl;

			if(count_min==1)
			{
				ofstream fem(ss0.str(), ios::trunc);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::trunc);
				ft<<T<<","<<En<<","<<E0<<endl;
				ft.close();
			}
			else
			{
				ofstream fem(ss0.str(), ios::app);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::app);
				ft<<T<<","<<En<<endl;
				ft.close();
			}

			ofstream fp(ss1.str());
			ofstream fl(ss3.str());
			for(int i=0;i<h_num;i++)
			{
				fp<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
				fl<<HYPER[i].lam<<","<<HYPER[i].mu<<endl;
			}
			fp.close();
			fl.close();
		}


		if(count_min>c_max)	break;
	}



	for(int i=0;i<h_num;i++)	delete[]	dG[i];
	for(int i=0;i<Nw;i++)	delete[]	dH[i];

	delete[]	dE;
	delete[]	rE;
	delete[]	G;
	delete[]	dG;
	delete[]	th_G;
	delete[]	H;
	delete[]	dH;
	delete[]	th_H;
	delete[]	dT;
	delete[]	rT;
	delete[]	d;
	delete[]	B;
}


void p_nab(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dE, double *G, double *H,double Dt, double V, double mi,vector<double > NEIw)
{

	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;
	/////////////////èâä˙âª///////////////////
	for(int i=0;i<h_num;i++)	G[i]=0.;
	for(int i=0;i<Nw;i++)	H[i]=0.;
	for(int i=0;i<Nx;i++)	dE[i]=0.;


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

			dE[k]-=0.5*Dt/mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
			G[k]+=1./mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
		}
		dE[k]-=0.5*Dt/mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
		//dE[k]+=V*(1.-HYPER[k].J);
		G[k]+=1./mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
	}
	for(int k=0;k<Nw;k++)
	{
		int kw=NEIw[k];
		pk[A_Z]=HYPER[kw].p[A_Z];
		dE[k+h_num]=0.5*Dt/mi*pk[A_Z];
		H[k]=-1./mi*pk[A_Z];
	}

}


void p_lap(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *rE, double **dG, double **dH, double Dt, double V, double mi,vector<double > NEIw)
{

	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=Nw+h_num;

	/////////////////èâä˙âª///////////////////
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nx;j++)	rE[i*Nx+j]=0.;
		for(int j=0;j<Nw;j++)	dH[j][i]=0.;
		for(int j=0;j<h_num;j++)	dG[j][i]=0.;

	}


	//////////íËã`
	double Dgjj[DIMENSION]={0,0,0};
	double Dgkj[DIMENSION]={0,0,0};
	double Dglj[DIMENSION]={0,0,0};
	for(int j=0;j<h_num;j++)
	{
		int Nj=HYPER[j].N;
		Dgjj[A_X]=HYPER1[j*h_num+j].DgDq[A_X];	Dgjj[A_Y]=HYPER1[j*h_num+j].DgDq[A_Y];	Dgjj[A_Z]=HYPER1[j*h_num+j].DgDq[A_Z];

		for(int k=0;k<Nj;k++)
		{
			int kn=HYPER[j].NEI[k];

			Dgkj[A_X]=HYPER1[kn*h_num+j].DgDq[A_X];	Dgkj[A_Y]=HYPER1[kn*h_num+j].DgDq[A_Y];	Dgkj[A_Z]=HYPER1[kn*h_num+j].DgDq[A_Z];
			for(int l=0;l<Nj;l++)
			{
				int ln=HYPER[j].NEI[l];
				/////////////rlam E
				Dglj[A_X]=HYPER1[ln*h_num+j].DgDq[A_X];	Dglj[A_Y]=HYPER1[ln*h_num+j].DgDq[A_Y];	Dglj[A_Z]=HYPER1[ln*h_num+j].DgDq[A_Z];
				rE[ln*Nx+kn]+=0.25*Dt*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
				/////////////dlam G
				dG[kn][ln]-=0.5*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
			}
			/////////////rlam E
			rE[j*Nx+kn]+=0.25*Dt*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			rE[kn*Nx+j]+=0.25*Dt*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

			/////////////dlam G
			dG[kn][j]-=0.5*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			dG[j][kn]-=0.5*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

		}
		/////////////rlam E
		rE[j*Nx+j]+=0.25*Dt*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);
		/////////////dlam G
		dG[j][j]-=0.5*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);
	}

	//////////íËã`
	double Dgki[DIMENSION]={0,0,0};
	double Dgii[DIMENSION]={0,0,0};
	for(int i=0;i<Nw;i++)
	{
		int iw=NEIw[i];
		int Ni=HYPER[iw].N;
		Dgii[A_Z]=HYPER1[iw*h_num+iw].DgDq[A_Z];

		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[iw].NEI[k];
			Dgki[A_Z]=HYPER1[kn*h_num+iw].DgDq[A_Z];
			/////////////rmulam E
			rE[(i+h_num)*Nx+kn]-=0.25*Dt*Dt/mi*Dgki[A_Z];
			/////////////rlammu E
			rE[kn*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgki[A_Z];

			/////////////dmu G
			dG[kn][i+h_num]+=0.5*Dt/mi*Dgki[A_Z];
			/////////////dlam H
			dH[i][kn]+=0.5*Dt/mi*Dgki[A_Z];
		}
		/////////////rmulam E
		rE[(i+h_num)*Nx+iw]-=0.25*Dt*Dt/mi*Dgii[A_Z];
		/////////////rlammu E
		rE[iw*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgii[A_Z];
		/////////////rmu E
		rE[(i+h_num)*Nx+i+h_num]=0.25*Dt*Dt/mi;
		/////////////dmu G
		dG[iw][i+h_num]+=0.5*Dt/mi*Dgii[A_Z];
		/////////////dlam H
		dH[i][iw]+=0.5*Dt/mi*Dgii[A_Z];
		/////////////dmu H
		dH[i][i+h_num]=-0.5*Dt/mi;
	}





}


void p_variables(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num,double Dt,double mi)
{
	/////////////påvéZ
	double p_p[DIMENSION]={0,0,0};
	double DgDq_ji[DIMENSION]={0,0,0};
	double DgDq_ii[DIMENSION]={0,0,0};
	double lam=0.;
	double mu=0.;
	int fw=0.;
	for(int i=0;i<h_num;i++)
	{
		p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			lam=HYPER[jn].lam;

			DgDq_ji[A_X]=HYPER1[jn*h_num+i].DgDq[A_X];	DgDq_ji[A_Y]=HYPER1[jn*h_num+i].DgDq[A_Y];	DgDq_ji[A_Z]=HYPER1[jn*h_num+i].DgDq[A_Z];
			p_p[A_X]+=HYPER[jn].stress[A_X][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_X][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_X][A_Z]*DgDq_ji[A_Z]-lam*DgDq_ji[A_X];
			p_p[A_Y]+=HYPER[jn].stress[A_Y][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_Y][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_Y][A_Z]*DgDq_ji[A_Z]-lam*DgDq_ji[A_Y];
			p_p[A_Z]+=HYPER[jn].stress[A_Z][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_Z][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_Z][A_Z]*DgDq_ji[A_Z]-lam*DgDq_ji[A_Z];
		}
		DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
		lam=HYPER[i].lam;
		p_p[A_X]+=HYPER[i].stress[A_X][A_X]*DgDq_ii[A_X]+HYPER[i].stress[A_X][A_Y]*DgDq_ii[A_Y]+HYPER[i].stress[A_X][A_Z]*DgDq_ii[A_Z]-lam*DgDq_ii[A_X];
		p_p[A_Y]+=HYPER[i].stress[A_Y][A_X]*DgDq_ii[A_X]+HYPER[i].stress[A_Y][A_Y]*DgDq_ii[A_Y]+HYPER[i].stress[A_Y][A_Z]*DgDq_ii[A_Z]-lam*DgDq_ii[A_Y];
		p_p[A_Z]+=HYPER[i].stress[A_Z][A_X]*DgDq_ii[A_X]+HYPER[i].stress[A_Z][A_Y]*DgDq_ii[A_Y]+HYPER[i].stress[A_Z][A_Z]*DgDq_ii[A_Z]-lam*DgDq_ii[A_Z];
		HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+0.5*Dt*p_p[A_X];
		HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+0.5*Dt*p_p[A_Y];
		HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+0.5*Dt*p_p[A_Z];

		fw=HYPER[i].fw;
		if(fw==1)
		{
			mu=HYPER[i].mu;
			HYPER[i].p[A_Z]+=0.5*Dt*mu;
		}
	}
}


void output_data_p(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double **dG, double *G, double *th_G, double **dH, double *H, double *th_H,double *rT, double *dT, double T, double *d, int count,int count_min,int t,double E,double En, double E0,vector<double > NEIw)
{
	int Nw=NEIw.size();
	int h_num=HYPER.size();
	int Nx=Nw+h_num;


	stringstream ss0;
	ss0<<"./p_G/dG_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss1;
	ss1<<"./p_G/G_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./p_T/rT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./p_T/dT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss4;
	ss4<<"./p_T/T_"<<t<<"_"<<count_min<<".csv";
	stringstream ss5;
	ss5<<"./p/p_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./E/p_E_"<<t<<"_"<<count_min<<".csv";
	stringstream ss7;
	ss7<<"./p_T/d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";

	stringstream ss8;
	ss8<<"./p_H/dH_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./p_H/H_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./lam/p_lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./mu/p_mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";


	if(count==1)
	{
		ofstream f_t(ss4.str(),ios::trunc);
		f_t<<T<<","<<En<<","<<E0<<endl;
		f_t.close();
		ofstream f_e(ss6.str(),ios::trunc);
		f_e<<E<<endl;
		f_e.close();
	}
	else
	{
		ofstream f_t(ss4.str(),ios::app);
		f_t<<T<<","<<En<<endl;
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
	ofstream f_lam(ss10.str());
	ofstream f_mu(ss11.str());


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
	for(int i=0;i<h_num;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_dG<<dG[i][j]<<",";
		}
		f_dG<<endl;
		f_G<<G[i]<<","<<th_G[i]<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		f_lam<<HYPER[i].lam<<endl;
		f_mu<<HYPER[i].mu<<endl;
	}
	for(int i=0;i<Nw;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_dH<<dH[i][j]<<",";
		}
		f_dH<<endl;
		f_H<<H[i]<<","<<th_H[i]<<endl;
	}

	f_dG.close();
	f_G.close();
	f_rt.close();
	f_dt.close();
	f_p.close();
	f_d.close();
	f_dH.close();
	f_H.close();
	f_lam.close();
	f_mu.close();

}




void q_QP_nw1(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, double V, double mi, double Dt, double E0,vector<double > NEIw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;
	int h_num=HYPER.size();
	int Nx=1+h_num;

	double E=1;
	double E_min=1;
	double E_sum=1;
	double ep=1.e-5;
	double ep_min=1.e-5;

	double rg=1.;
	double rh=1.;

	double T=0.;
	double *dT=new double [Nx];	
	double *rT=new double [Nx*Nx];

	double En=0.;
	double *dE=new double [Nx];
	double *rE=new double [Nx*Nx];
	double *p_rE=new double [Nx*Nx];


	double *g=new double [h_num];
	double **dg=new double *[h_num];
	for(int i=0;i<h_num;i++)	dg[i]=new double [Nx];
	double *th_g=new double [h_num];

	double h=0.;
	double *dh=new double [Nx];
	double th_h=0.;

	double *d=new double [Nx];
	double *B=new double [Nx*Nx];


	int fw=0;
	int iw=0;
	double p_p[DIMENSION]={0,0,0};
	double hp[DIMENSION]={0,0,0};
	double W=0.;
	double Dgki_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};
	double Dgji_n[DIMENSION]={0,0,0};
	double Dgli_n[DIMENSION]={0,0,0};
	double lam=0.;
	double mu=0.;
	double hi=0.;
	double d_sum=0.;

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].lam=1.;
		HYPER[i].mu=1.;
		g[i]=0.;
		th_g[i]=0.;
		for(int j=0;j<Nx;j++)	dg[i][j]=0.;
		dh[i]=0.;
	}
	for(int i=0;i<Nx;i++)
	{
		dT[i]=0.;
		d[i]=0.;
		dE[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			rT[i*Nx+j]=0.;
			rE[i*Nx+j]=0.;
			p_rE[i*Nx+j]=0.;
			B[i*Nx+j]=0.;
		}
	}


	////////////dh, p_rEåvéZ///////////////
	for(int i=0;i<Nw;i++)
	{
		iw=NEIw[i];
		int Ni=HYPER[iw].N;
		Dgii_n[A_Z]=HYPER1[iw*h_num+iw].DgDq_n[A_Z];
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[iw].NEI[k];
			Dgki_n[A_Z]=HYPER1[kn*h_num+iw].DgDq_n[A_Z];
			dh[i][kn]+=0.5*Dt*Dt/mi*Dgki_n[A_Z];
			p_rE[(i+h_num)*Nx+kn]-=0.25*Dt*Dt/mi*Dgki_n[A_Z];
			p_rE[kn*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgki_n[A_Z];
		}
		dh[i][iw]+=0.5*Dt*Dt/mi*Dgii_n[A_Z];
		dh[i][i+h_num]-=0.5*Dt*Dt/mi;
		p_rE[(i+h_num)*Nx+iw]-=0.25*Dt*Dt/mi*Dgii_n[A_Z];
		p_rE[iw*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgii_n[A_Z];
		p_rE[(i+h_num)*Nx+i+h_num]+=0.25*Dt*Dt/mi;
	}
	for(int i=0;i<h_num;i++)
	{
		Dgii_n[A_X]=HYPER1[i*h_num+i].DgDq_n[A_X];	Dgii_n[A_Y]=HYPER1[i*h_num+i].DgDq_n[A_Y];	Dgii_n[A_Z]=HYPER1[i*h_num+i].DgDq_n[A_Z];

		int Ni=HYPER[i].N;
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[i].NEI[k];
			Dgki_n[A_X]=HYPER1[kn*h_num+i].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[kn*h_num+i].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[kn*h_num+i].DgDq_n[A_Z];
			for(int l=0;l<Ni;l++)
			{
				int ln=HYPER[i].NEI[l];
				Dgli_n[A_X]=HYPER1[ln*h_num+i].DgDq_n[A_X];	Dgli_n[A_Y]=HYPER1[ln*h_num+i].DgDq_n[A_Y];	Dgli_n[A_Z]=HYPER1[ln*h_num+i].DgDq_n[A_Z];
				p_rE[ln*Nx+kn]+=0.25*Dt*Dt/mi*(Dgki_n[A_X]*Dgli_n[A_X]+Dgki_n[A_Y]*Dgli_n[A_Y]+Dgki_n[A_Z]*Dgli_n[A_Z]);
			}
			p_rE[kn*Nx+i]+=0.25*Dt*Dt/mi*(Dgii_n[A_X]*Dgki_n[A_X]+Dgii_n[A_Y]*Dgki_n[A_Y]+Dgii_n[A_Z]*Dgki_n[A_Z]);
			p_rE[i*Nx+kn]+=0.25*Dt*Dt/mi*(Dgki_n[A_X]*Dgii_n[A_X]+Dgki_n[A_Y]*Dgii_n[A_Y]+Dgki_n[A_Z]*Dgii_n[A_Z]);
		}
		p_rE[i*Nx+i]+=0.25*Dt*Dt/mi*(Dgii_n[A_X]*Dgii_n[A_X]+Dgii_n[A_Y]*Dgii_n[A_Y]+Dgii_n[A_Z]*Dgii_n[A_Z]);
	}






	////////////QPñ@åvéZ///////////////
	while(E_min>ep_min)
	{
		count_min++;

		for(int i=0;i<h_num;i++)	HYPER[i].old_lam=HYPER[i].lam;
			iw=NEIw[0];
			HYPER[iw].old_mu=HYPER[iw].mu;


		E=1;
		count=0;
		while(E>ep)
		{
			count++;

			q_variables(CON,PART,HYPER,HYPER1);	//Ç±Ç±Ç‹Ç≈1219
			q_nab_lap(PART,HYPER,HYPER1,dE,rE,dg,Dt,V,mi,NEIw);

			En=0.;
			for(int i=0;i<h_num;i++)
			{
				W=HYPER[i].W;
				hp[A_X]=HYPER[i].half_p[A_X];	hp[A_Y]=HYPER[i].half_p[A_Y];	hp[A_Z]=HYPER[i].half_p[A_Z];
				En+=0.5/mi*(hp[A_X]*hp[A_X]+hp[A_Y]*hp[A_Y]+hp[A_Z]*hp[A_Z])+V*W;
			}

			T=(En-E0)*(En-E0);
			//cout<<"T="<<T<<", En"<<En<<endl;
			for(int k=0;k<Nx;k++)
			{	
				dT[k]=2.*dE[k]*(En-E0);
				for(int l=0;l<Nx;l++)
				{
					rE[l*Nx+k]+=p_rE[l*Nx+k];
					rT[l*Nx+k]=2.*rE[l*Nx+k]*(En-E0)+2.*dE[k]*dE[l];
				}
			}
			//for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			for(int i=0;i<h_num;i++)
			{	
				g[i]=V*(1.-HYPER[i].J);

				if(g[i]<0)	T+=0.5*rg*(-1.*g[i]+th_g[i])*(-1.*g[i]+th_g[i]);
				else
				{
					T+=0.5*rg*(g[i]+th_g[i])*(g[i]+th_g[i]);
				}
			}
			for(int i=0;i<Nw;i++)
			{	
				int iw=NEIw[i];
				h[i]=-1.*PART[iw].r[A_Z];

				//if(h[i]<0)	T+=0.5*rh*(-1.*h[i]+th_h[i])*(-1.*h[i]+th_h[i]);
				//else
				//{
				//	T+=0.5*rh*(h[i]+th_h[i])*(h[i]+th_h[i]);
				//}
				if(h[i]+th_h[i]>0)	T+=0.5*rh*(h[i]+th_h[i])*(h[i]+th_h[i]);
			}




			for(int k=0;k<Nx;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					if(g[i]<0)	dT[k]+=-1.*rg*dg[i][k]*(-1.*g[i]+th_g[i]);
					else
					{
						dT[k]+=rg*dg[i][k]*(g[i]+th_g[i]);
					}
					//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
				}
				for(int i=0;i<Nw;i++)
				{	
					//if(h[i]<0)	dT[k]+=-1.*rh*dh[i][k]*(-1.*h[i]+th_h[i]);
					//else
					//{
					//	dT[k]+=rh*dh[i][k]*(h[i]+th_h[i]);
					//}
					if(h[i]+th_h[i]>0)	dT[k]+=rh*dh[i][k]*(h[i]+th_h[i]);
				}
				for(int l=0;l<Nx;l++)
				{
					for(int i=0;i<h_num;i++)
					{	
						rT[l*Nx+k]+=rg*dg[i][k]*dg[i][l];
					}
					for(int i=0;i<Nw;i++)
					{	
						//rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						if(h[i]+th_h[i]>0)	rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						else if(h[i]+th_h[i]>-1.e-20 && dh[i][k]*dh[i][l]>0)
						{
							rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						}
					}
				}
			}

			d_sum=0.;
			for(int i=0;i<Nx;i++)
			{
				d[i]=dT[i];
				d_sum+=fabs(dT[i]);
				for(int j=0;j<Nx;j++)
				{
					B[i*Nx+j]=rT[i*Nx+j];
				}
			}
			if(d_sum<1.e-20)	break;

			gauss(B,d,Nx);

			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=d[i];
				E_sum+=fabs(d[i]);
			}
			for(int i=0;i<Nw;i++)
			{
				iw=NEIw[i];
				HYPER[iw].mu-=d[i+h_num];
				E_sum+=fabs(d[i+h_num]);
			}

			E=E_sum;




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
				if(fw==1)
				{
					mu=HYPER[i].mu;
					HYPER[i].half_p[A_Z]+=0.5*Dt*mu;
				}
				/////////////qåvéZ
				PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
				PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
				PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
			}	
			if(count==1||count%10000==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				output_data(PART,HYPER,HYPER1,T,dT,rT,g,dg,th_g,h,dh,th_h,d,count,count_min,t,E,En,E0,mi,V,NEIw);

			}

		}

		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=fabs(HYPER[i].old_lam-HYPER[i].lam);
		iw=NEIw[0];
		E_sum+=fabs(HYPER[iw].old_mu-HYPER[iw].mu);

		E_min=E_sum;
		//if(E_min<ep_min*1000)
		//{
		//	rg*=4;
		//	rh*=4;
		//}

		for(int i=0;i<h_num;i++)
		{
			if(g[i]<0)	th_g[i]+=-1.*g[i];
			else
			{
				th_g[i]+=g[i];
			}
		}
		if(h+th_h>0)	th_h+=h;

		cout<<"Emin"<<count_min<<"="<<E_min<<endl;
		//if(count_min==1||count_min%100==0)
		{
			stringstream ss0;
			ss0<<"./E_min/E"<<t<<".csv"<<endl;
			stringstream ss1;
			ss1<<"./p/hp"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss2;
			ss2<<"./q/q"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss3;
			ss3<<"./lam/lam"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss4;
			ss4<<"./T/T"<<t<<".csv"<<endl;

			if(count_min==1)
			{
				ofstream fem(ss0.str(), ios::trunc);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::trunc);
				ft<<T<<","<<En<<","<<E0<<endl;
				ft.close();
			}
			else
			{
				ofstream fem(ss0.str(), ios::app);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::app);
				ft<<T<<","<<En<<endl;
				ft.close();
			}

			ofstream fhp(ss1.str());
			ofstream fq(ss2.str());
			ofstream fl(ss3.str());
			for(int i=0;i<h_num;i++)
			{
				fhp<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
				fq<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;			
				fl<<HYPER[i].lam<<","<<HYPER[i].mu<<endl;
			}
			fhp.close();
			fq.close();
			fl.close();

		}

		//NEIw.clear();
		//for(int i=0;i<h_num;i++)
		//{
		//	hi=-1.*PART[i].r[A_Z];
		//	if(hi>0)
		//	{
		//		NEIw.push_back(i);
		//		HYPER[i].fw=1;
		//	}
		//	else
		//	{
		//		HYPER[i].fw=0;
		//	}
		//}
		//Nw=NEIw.size();
		//if(Nw>0)
		//{
		//	cout<<"ê⁄êGó±éq ";
		//	for(int i=0;i<Nw;i++)	cout<<NEIw[i]<<", ";
		//	cout<<endl;
		//}


	}

	delete[]	dE;
	delete[]	dT;
	delete[]	rE;
	delete[]	rT;
	delete[]	g;
	for(int i=0;i<h_num;i++)	delete[]	dg[i];
	delete[]	dg;
	delete[]	dh;
	delete[]	th_g;
	delete[]	d;
	delete[]	B;
	delete[]	p_rE;

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
		if(fw==1)
		{
			mu=HYPER[i].mu;
			HYPER[i].half_p[A_Z]+=0.5*Dt*mu;
		}
		/////////////qåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}	

}


void q_nab_lap_nw1(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,vector<double > NEIw)
{

	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;
	/////////////////èâä˙âª///////////////////
	for(int k=0;k<Nx;k++)
	{
		dE[k]=0.;
		for(int i=0;i<h_num;i++)	dg[i][k]=0.;
		for(int i=0;i<Nx;i++)	rE[i*Nx+k]=0.;
	}


	/////////////////dE, dgåvéZ///////////////////
	//////////íËã`
	double Dgkk_n[DIMENSION]={0,0,0};
	double Dgki_n[DIMENSION]={0,0,0};

	double Dgik[DIMENSION]={0,0,0};
	double Dgik_n[DIMENSION]={0,0,0};
	double Dgkk[DIMENSION]={0,0,0};
	double	p_Eik[DIMENSION]={0,0,0};
	double	p_Ekk[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		/////////////dlam E1
		Dgkk_n[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];	Dgkk_n[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];	Dgkk_n[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];
		/////////////dlam E1
		Dgkk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	Dgkk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];
		p_Ekk[A_X]=HYPER[k].stress[A_X][0]*Dgkk[0]+HYPER[k].stress[A_X][1]*Dgkk[1]+HYPER[k].stress[A_X][2]*Dgkk[2];
		p_Ekk[A_Y]=HYPER[k].stress[A_Y][0]*Dgkk[0]+HYPER[k].stress[A_Y][1]*Dgkk[1]+HYPER[k].stress[A_Y][2]*Dgkk[2];
		p_Ekk[A_Z]=HYPER[k].stress[A_Z][0]*Dgkk[0]+HYPER[k].stress[A_Z][1]*Dgkk[1]+HYPER[k].stress[A_Z][2]*Dgkk[2];
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			/////////////dlam E1
			Dgki_n[A_X]=HYPER1[k*h_num+in].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[k*h_num+in].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[k*h_num+in].DgDq_n[A_Z];
			dE[k]-=0.5*Dt/mi*(Dgki_n[A_X]*HYPER[in].half_p[A_X]+Dgki_n[A_Y]*HYPER[in].half_p[A_Y]+Dgki_n[A_Z]*HYPER[in].half_p[A_Z]);
			/////////////dlam E2
			Dgik[A_X]=HYPER1[in*h_num+k].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+k].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+k].DgDq[A_Z];
			p_Eik[A_X]=HYPER[in].stress[A_X][0]*Dgik[0]+HYPER[in].stress[A_X][1]*Dgik[1]+HYPER[in].stress[A_X][2]*Dgik[2];
			p_Eik[A_Y]=HYPER[in].stress[A_Y][0]*Dgik[0]+HYPER[in].stress[A_Y][1]*Dgik[1]+HYPER[in].stress[A_Y][2]*Dgik[2];
			p_Eik[A_Z]=HYPER[in].stress[A_Z][0]*Dgik[0]+HYPER[in].stress[A_Z][1]*Dgik[1]+HYPER[in].stress[A_Z][2]*Dgik[2];

			Dgik_n[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];	Dgik_n[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];	Dgik_n[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];

				/////////////dlam E2
				dE[jn]+=0.5*Dt*Dt/mi*(p_Eik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+p_Eik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+p_Eik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
				/////////////dlam g
				dg[in][jn]-=0.5*Dt*Dt/mi*(Dgik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+Dgik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+Dgik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);

			}
			/////////////dlam E2
			dE[k]+=0.5*Dt*Dt/mi*(p_Eik[A_X]*Dgkk_n[A_X]+p_Eik[A_Y]*Dgkk_n[A_Y]+p_Eik[A_Z]*Dgkk_n[A_Z]);
			dE[in]+=0.5*Dt*Dt/mi*(p_Ekk[A_X]*Dgik_n[A_X]+p_Ekk[A_Y]*Dgik_n[A_Y]+p_Ekk[A_Z]*Dgik_n[A_Z]);
			/////////////dlam g
			dg[in][k]-=0.5*Dt*Dt/mi*(Dgik[A_X]*Dgkk_n[A_X]+Dgik[A_Y]*Dgkk_n[A_Y]+Dgik[A_Z]*Dgkk_n[A_Z]);
			dg[k][in]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgik_n[A_X]+Dgkk[A_Y]*Dgik_n[A_Y]+Dgkk[A_Z]*Dgik_n[A_Z]);
		}
		/////////////dlam E1
		dE[k]-=0.5*Dt/mi*(Dgkk_n[A_X]*HYPER[k].half_p[A_X]+Dgkk_n[A_Y]*HYPER[k].half_p[A_Y]+Dgkk_n[A_Z]*HYPER[k].half_p[A_Z]);
		/////////////dlam E2
		dE[k]+=0.5*Dt*Dt/mi*(p_Ekk[A_X]*Dgkk_n[A_X]+p_Ekk[A_Y]*Dgkk_n[A_Y]+p_Ekk[A_Z]*Dgkk_n[A_Z]);
		/////////////dlam g
		dg[k][k]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgkk_n[A_X]+Dgkk[A_Y]*Dgkk_n[A_Y]+Dgkk[A_Z]*Dgkk_n[A_Z]);
	}


	//////////íËã`
	for(int k=0;k<Nw;k++)
	{
		int kw=NEIw[k];

		/////////////dmu E2
		Dgkk[A_X]=HYPER1[kw*h_num+kw].DgDq[A_X];	Dgkk[A_Y]=HYPER1[kw*h_num+kw].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[kw*h_num+kw].DgDq[A_Z];
		p_Ekk[A_Z]=HYPER[kw].stress[A_Z][0]*Dgkk[0]+HYPER[kw].stress[A_Z][1]*Dgkk[1]+HYPER[kw].stress[A_Z][2]*Dgkk[2];

		int Nk=HYPER[kw].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[kw].NEI[i];

			/////////////dmu E2
			Dgik[A_X]=HYPER1[in*h_num+kw].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+kw].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+kw].DgDq[A_Z];

			p_Eik[A_Z]=HYPER[in].stress[A_Z][0]*Dgik[0]+HYPER[in].stress[A_Z][1]*Dgik[1]+HYPER[in].stress[A_Z][2]*Dgik[2];
			dE[k+h_num]-=0.5*Dt*Dt/mi*p_Eik[A_Z];

			/////////////dmu g
			dg[in][k+h_num]+=0.5*Dt*Dt/mi*Dgik[A_Z];
		}
		/////////////dmu E2
		dE[k+h_num]-=0.5*Dt*Dt/mi*p_Ekk[A_Z];

		/////////////dmu E1
		dE[k+h_num]+=0.5*Dt/mi*HYPER[kw].half_p[A_Z];

		/////////////dmu g
		dg[kw][k+h_num]+=0.5*Dt*Dt/mi*Dgkk[A_Z];
	}



	/////////////////rEåvéZ///////////////////
	//////////íËã`
	/////////////rlam E1
	double Dgli_n[DIMENSION]={0,0,0};
	double Dpki[DIMENSION]={0,0,0};
	double Dpii[DIMENSION]={0,0,0};
	double *DpiDlam=new double [h_num*h_num];
	double *DpiDmu_x=new double [h_num*h_num];
	double *DpiDmu_y=new double [h_num*h_num];
	double *DpiDmu_z=new double [h_num*h_num];
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			DpiDlam[i*h_num+j]=0.;
			DpiDmu_x[i*h_num+j]=0.;
			DpiDmu_y[i*h_num+j]=0.;
			DpiDmu_z[i*h_num+j]=0.;

		}
	}
	/////////////rlam E3
	double Dgli[DIMENSION]={0,0,0};
	/////////////rlam E4
	double Dgki[DIMENSION]={0,0,0};
	double Dgii[DIMENSION]={0,0,0};
	/////////////dpidlam
	double Dgii_n[DIMENSION]={0,0,0};
	/////////////dpidmu
	double Dpik[DIMENSION]={0,0,0};
	double n0li[DIMENSION]={0,0,0};
	double n0ii[DIMENSION]={0,0,0};
	double n0ki[DIMENSION]={0,0,0};

	for(int i=0;i<h_num;i++)
	{
		/////////////dpidlama
		Dpii[A_X]=HYPER1[i*h_num+i].DpiDq[A_X];	Dpii[A_Y]=HYPER1[i*h_num+i].DpiDq[A_Y];	Dpii[A_Z]=HYPER1[i*h_num+i].DpiDq[A_Z];
		/////////////dpidmu
		n0ii[A_X]=HYPER1[i*h_num+i].n0ij[A_X];	n0li[A_Y]=HYPER1[i*h_num+i].n0ij[A_Y];	n0li[A_Z]=HYPER1[i*h_num+i].n0ij[A_Z];
		int Ni=HYPER[i].N;
		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[i].NEI[k];
			/////////////dpidmu
			Dgki_n[A_X]=HYPER1[kn*h_num+i].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[kn*h_num+i].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[kn*h_num+i].DgDq_n[A_Z];
			Dpki[A_X]=HYPER1[kn*h_num+i].DpiDq[A_X];	Dpki[A_Y]=HYPER1[kn*h_num+i].DpiDq[A_Y];	Dpki[A_Z]=HYPER1[kn*h_num+i].DpiDq[A_Z];
			Dpik[A_X]=HYPER1[i*h_num+kn].DpiDq[A_X];	Dpik[A_Y]=HYPER1[i*h_num+kn].DpiDq[A_Y];	Dpik[A_Z]=HYPER1[i*h_num+kn].DpiDq[A_Z];
			n0ki[A_X]=HYPER1[kn*h_num+i].n0ij[A_X];	n0ki[A_Y]=HYPER1[kn*h_num+i].n0ij[A_Y];	n0ki[A_Z]=HYPER1[kn*h_num+i].n0ij[A_Z];
			for(int l=0;l<Ni;l++)
			{
				int ln=HYPER[i].NEI[l];
				/////////////dpidlam
				Dgli_n[A_X]=HYPER1[ln*h_num+i].DgDq_n[A_X];	Dgli_n[A_Y]=HYPER1[ln*h_num+i].DgDq_n[A_Y];	Dgli_n[A_Z]=HYPER1[ln*h_num+i].DgDq_n[A_Z];
				DpiDlam[kn*h_num+ln]+=Dpki[A_X]*Dgli_n[A_X]+Dpki[A_Y]*Dgli_n[A_Y]+Dpki[A_Z]*Dgli_n[A_Z];
				/////////////dpidmu
				n0li[A_X]=HYPER1[ln*h_num+i].n0ij[A_X];	n0li[A_Y]=HYPER1[ln*h_num+i].n0ij[A_Y];	n0li[A_Z]=HYPER1[ln*h_num+i].n0ij[A_Z];
				DpiDmu_x[kn*h_num+ln]+=Dpik[A_Z]*n0li[A_X];	DpiDmu_y[kn*h_num+ln]+=Dpik[A_Z]*n0li[A_Y];	DpiDmu_z[kn*h_num+ln]+=Dpik[A_Z]*n0li[A_Z];

			}
			/////////////dpidlama
			DpiDlam[kn*h_num+i]+=Dpki[A_X]*Dgii_n[A_X]+Dpki[A_Y]*Dgii_n[A_Y]+Dpki[A_Z]*Dgii_n[A_Z];
			DpiDlam[i*h_num+kn]+=Dpii[A_X]*Dgki_n[A_X]+Dpii[A_Y]*Dgki_n[A_Y]+Dpii[A_Z]*Dgki_n[A_Z];
			/////////////dpidmu
			DpiDmu_x[kn*h_num+i]+=Dpik[A_Z]*n0ii[A_X];	DpiDmu_y[kn*h_num+i]+=Dpik[A_Z]*n0ii[A_Y];	DpiDmu_z[kn*h_num+i]+=Dpik[A_Z]*n0ii[A_Z];
			DpiDmu_x[i*h_num+kn]+=Dpii[A_Z]*n0ki[A_X];	DpiDmu_y[i*h_num+kn]+=Dpii[A_Z]*n0ki[A_Y];	DpiDmu_z[i*h_num+kn]+=Dpii[A_Z]*n0ki[A_Z];

		}
		/////////////dpidlama
		DpiDlam[i*h_num+i]+=Dpii[A_X]*Dgii_n[A_X]+Dpii[A_Y]*Dgii_n[A_Y]+Dpii[A_Z]*Dgii_n[A_Z];
		/////////////dpidmu
		DpiDmu_x[i*h_num+i]+=Dpii[A_Z]*n0ii[A_X];	DpiDmu_y[i*h_num+i]+=Dpii[A_Z]*n0ii[A_Y];	DpiDmu_z[i*h_num+i]+=Dpii[A_Z]*n0ii[A_Z];
	}

	/////////////rlam E2
	double Dgkm_n[DIMENSION]={0,0,0};
	double n0mi[DIMENSION]={0,0,0};
	//double n0ii[DIMENSION]={0,0,0};
	//double n0li[DIMENSION]={0,0,0};
	double Dgkl_n[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		for(int i=0;i<h_num;i++)
		{
			int Ni=HYPER[i].N;
			n0ii[A_X]=HYPER1[i*h_num+i].n0ij[A_X];	n0ii[A_Y]=HYPER1[i*h_num+i].n0ij[A_Y];	n0ii[A_Z]=HYPER1[i*h_num+i].n0ij[A_Z];
			Dgki_n[A_X]=HYPER1[k*h_num+i].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[k*h_num+i].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[k*h_num+i].DgDq_n[A_Z];
			for(int l=0;l<Ni;l++)
			{
				int ln=HYPER[i].NEI[l];
				n0li[A_X]=HYPER1[ln*h_num+i].n0ij[A_X];	n0li[A_Y]=HYPER1[ln*h_num+i].n0ij[A_Y];	n0li[A_Z]=HYPER1[ln*h_num+i].n0ij[A_Z];
				Dgkl_n[A_X]=HYPER1[k*h_num+ln].DgDq_n[A_X];	Dgkl_n[A_Y]=HYPER1[k*h_num+ln].DgDq_n[A_Y];	Dgkl_n[A_Z]=HYPER1[k*h_num+ln].DgDq_n[A_Z];

				for(int m=0;m<Ni;m++)
				{
					int mn=HYPER[i].NEI[m];
					n0mi[A_X]=HYPER1[mn*h_num+i].n0ij[A_X];	n0mi[A_Y]=HYPER1[mn*h_num+i].n0ij[A_Y];	n0mi[A_Z]=HYPER1[mn*h_num+i].n0ij[A_Z];
					Dgkm_n[A_X]=HYPER1[k*h_num+mn].DgDq_n[A_X];	Dgkm_n[A_Y]=HYPER1[k*h_num+mn].DgDq_n[A_Y];	Dgkm_n[A_Z]=HYPER1[k*h_num+mn].DgDq_n[A_Z];

					rE[ln*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+ln]*(n0mi[A_X]*Dgkm_n[A_X]+n0mi[A_Y]*Dgkm_n[A_Y]+n0mi[A_Z]*Dgkm_n[A_Z]);
				}
				rE[ln*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+ln]*(n0ii[A_X]*Dgki_n[A_X]+n0ii[A_Y]*Dgki_n[A_Y]+n0ii[A_Z]*Dgki_n[A_Z]);
				rE[i*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+i]*(n0li[A_X]*Dgkl_n[A_X]+n0li[A_Y]*Dgkl_n[A_Y]+n0li[A_Z]*Dgkl_n[A_Z]);
			}
			rE[i*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+i]*(n0ii[A_X]*Dgki_n[A_X]+n0ii[A_Y]*Dgki_n[A_Y]+n0ii[A_Z]*Dgki_n[A_Z]);
		}
	}

	/////////////rmu E2
	for(int l=0;l<Nw;l++)
	{
		int lw=NEIw[l];

		for(int k=0;k<Nw;k++)
		{
			int kw=NEIw[k];
			/////////////rmu E2
			rE[(l+h_num)*Nx+(k+h_num)]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDmu_z[lw*h_num+kw];
		}
	}


	/////////////rmulam E2
	double Dplm[DIMENSION]={0,0,0};
	for(int l=0;l<Nw;l++)
	{
		int lw=NEIw[l];

		int Nl=HYPER[lw].N;
		for(int k=0;k<Nl;k++)
		{
			int kn=HYPER[lw].NEI[k];;
			for(int m=0;m<h_num;m++)
			{
				/////////////rmulam E2
				Dgkm_n[A_X]=HYPER1[kn*h_num+m].DgDq_n[A_X];	Dgkm_n[A_Y]=HYPER1[kn*h_num+m].DgDq_n[A_Y];	Dgkm_n[A_Z]=HYPER1[kn*h_num+m].DgDq_n[A_Z];
				Dplm[A_X]=DpiDmu_x[lw*h_num+m];	Dplm[A_Y]=DpiDmu_y[lw*h_num+m];	Dplm[A_Z]=DpiDmu_z[lw*h_num+m];
				rE[(l+h_num)*Nx+kn]-=0.25*Dt*Dt*Dt*Dt/mi/mi/V*(Dplm[A_X]*Dgkm_n[A_X]+Dplm[A_Y]*Dgkm_n[A_Y]+Dplm[A_Z]*Dgkm_n[A_Z]);
			}
		}
	}


	for(int k=0;k<Nw;k++)
	{
		int kw=NEIw[k];
		int Nk=HYPER[kw].N;

		for(int l=0;l<h_num;l++)
		{
			for(int i=0;i<h_num;i++)
			{
				/////////////rlammu E2
				n0ki[A_X]=HYPER1[kw*h_num+i].n0ij[A_X];	n0ki[A_Y]=HYPER1[kw*h_num+i].n0ij[A_Y];	n0ki[A_Z]=HYPER1[kw*h_num+i].n0ij[A_Z];
				rE[l*h_num+k+h_num]-=0.25*Dt*Dt*Dt*Dt/mi/mi/V*DpiDlam[i*h_num+l]*n0ki[A_Z];
			}
		}

	}



	delete[]	DpiDlam;
	delete[]	DpiDmu_x;
	delete[]	DpiDmu_y;
	delete[]	DpiDmu_z;

}


void output_data_nw1(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double **dg, double *th_g, double *h, double **dh, double *th_h, double *d, int count,int count_min,int t,double E,double En,double E0,double mi, double V,vector<double > NEIw)
{
	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;

	stringstream ss0;
	ss0<<"./T/T_"<<t<<"_"<<count_min<<".csv";
	stringstream ss1;
	ss1<<"./T/dT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./T/rT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./lam/lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./g/g_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss7;
	ss7<<"./g/dg_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss8;
	ss8<<"./p/hp_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./q/q_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./T/d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./E/E_"<<t<<"_"<<count_min<<".csv";
	//stringstream ss12;
	//ss12<<"./h/h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss13;
	//ss13<<"./g/J_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss14;
	//ss14<<"./g/DgDq_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	//stringstream ss15;
	//ss15<<"./g/stress_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss16;
	ss16<<"./T/En_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss17;
	ss17<<"./h/h_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss18;
	ss18<<"./h/dh_"<<t<<".csv";
	stringstream ss19;
	ss19<<"./mu/mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";


	ofstream f_dt(ss1.str());
	ofstream f_rt(ss2.str());
	ofstream f_lam(ss3.str());
	ofstream f_g(ss6.str());
	ofstream f_dg(ss7.str());
	ofstream f_p(ss8.str());
	ofstream f_q(ss9.str());
	ofstream f_d(ss10.str());
	//ofstream f_h(ss12.str());
	//ofstream f_J(ss13.str());
	//ofstream f_dgdq(ss14.str());
	//ofstream f_s(ss15.str());
	ofstream f_en(ss16.str());
	ofstream f_h(ss17.str());

	ofstream f_mu(ss19.str());

	if(count_min==1&&count==1)
	{
		ofstream f_dh(ss18.str());
		for(int i=0;i<Nw;i++)
		{
			for(int k=0;k<Nx;k++)
			{
				f_dh<<dh[i][k]<<",";
			}
			f_dh<<endl;
		}
		f_dh.close();
	}	
	if(count==1)
	{
		ofstream f_T(ss0.str(), ios::trunc);
		f_T<<count<<","<<T<<","<<En<<endl;
		f_T.close();
		ofstream f_E(ss11.str(), ios::trunc);
		f_E<<count<<","<<E<<endl;
		f_E.close();
	}	
	else
	{
		ofstream f_T(ss0.str(), ios::app);
		f_T<<count<<","<<T<<","<<En<<endl;
		f_T.close();
		ofstream f_E(ss11.str(), ios::app);
		f_E<<count<<","<<E<<endl;
		f_E.close();
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
	for(int i=0;i<h_num;i++)
	{
		f_lam<<HYPER[i].lam<<endl;
		f_mu<<HYPER[i].mu<<endl;

		for(int k=0;k<Nx;k++)
		{
			f_dg<<dg[i][k]<<",";
		}
		f_dg<<endl;

		f_g<<g[i]<<","<<th_g[i]<<endl;
		f_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
		f_q<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		//f_J<<HYPER[i].J<<endl;

		//f_dgdq<<i<<endl;
		//int Ni=HYPER[i].N;
		//for(int j=0;j<Ni;j++)
		//{
		//	int jn=HYPER[i].NEI[j];
		//	f_dgdq<<","<<jn<<","<<HYPER1[i*h_num+jn].DgDq[A_X]<<","<<HYPER1[i*h_num+jn].DgDq[A_Y]<<","<<HYPER1[i*h_num+jn].DgDq[A_Z]<<endl;
		//}

		////f_h<<h[i]<<","<<th_h[i]<<endl;
		//f_s<<i<<","<<HYPER[i].stress[A_X][0]<<","<<HYPER[i].stress[A_X][1]<<","<<HYPER[i].stress[A_X][2]<<endl;
		//f_s<<","<<HYPER[i].stress[A_Y][0]<<","<<HYPER[i].stress[A_Y][1]<<","<<HYPER[i].stress[A_Y][2]<<endl;
		//f_s<<","<<HYPER[i].stress[A_Z][0]<<","<<HYPER[i].stress[A_Z][1]<<","<<HYPER[i].stress[A_Z][2]<<endl;
		f_en<<0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])+V*HYPER[i].W<<","<<0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])<<","<<V*HYPER[i].W<<endl;
	}
	for(int i=0;i<Nw;i++)
	{
		f_h<<h[i]<<","<<th_h[i]<<endl;
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
	//f_J.close();
	//f_dgdq.close();
	//f_s.close();
	f_en.close();
	f_h.close();
	f_mu.close();

}




void p_QP_nw1(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1, int t,double V, double mi, double Dt, double E0, vector<double > NEIw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;
	int Nw=NEIw.size();
	int h_num=HYPER.size();
	int Nx=Nw+h_num;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1.e-5;
	double ep_min=1.e-5;
	double d_sum=0.;

	double p_p[DIMENSION]={0,0,0};
	double p[DIMENSION]={0,0,0};
	double W=0.;
	double g=0.;
	double rg=1.0;
	double rh=1.0;

	double En=0.;
	double *dE=new double [Nx];	
	double *rE=new double [Nx*Nx];

	double T=0.;
	double *dT=new double [Nx];	
	double *rT=new double [Nx*Nx];

	double *d=new double [Nx];
	double *B=new double [Nx*Nx];

	double *G=new double [h_num];
	double *th_G=new double [h_num];
	double **dG=new double *[h_num];

	double *H=new double [Nw];
	double *th_H=new double [Nw];
	double **dH=new double *[Nw];

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)	dG[i]=new double [Nx];
	for(int i=0;i<Nw;i++)	dH[i]=new double [Nx];

	for(int i=0;i<h_num;i++)
	{
		G[i]=0.;
		th_G[i]=0.;
		for(int j=0;j<Nx;j++)	dG[i][j]=0.;
	}
	for(int i=0;i<Nw;i++)
	{
		H[i]=0.;
		th_H[i]=0.;
		for(int j=0;j<Nx;j++)	dH[i][j]=0.;
	}

	for(int i=0;i<Nx;i++)
	{
		dE[i]=0.;
		dT[i]=0.;
		d[i]=0.;
		for(int j=0;j<Nx;j++)
		{
			rE[i*Nx+j]=0.;
			rT[i*Nx+j]=0.;
			B[i*Nx+j]=0.;
		}
	}

	p_lap(HYPER,HYPER1,rE,dG,dH,Dt,V,mi,NEIw);



	while(E_min>ep_min)
	{
		count_min++;

		for(int i=0;i<h_num;i++)	HYPER[i].old_lam=HYPER[i].lam;
		for(int i=0;i<Nw;i++)
		{
			int iw=NEIw[i];
			HYPER[iw].old_mu=HYPER[iw].mu;		
		}

		E=1;
		count=0;
		while(E>ep)
		{
			count++;
			p_variables(HYPER,HYPER1,h_num,Dt,mi);
			p_nab(HYPER,HYPER1,dE,G,H,Dt,V,mi,NEIw);

			En=0;
			for(int i=0;i<h_num;i++)	
			{

				p[A_X]=HYPER[i].p[A_X];	p[A_Y]=HYPER[i].p[A_Y];	p[A_Z]=HYPER[i].p[A_Z];
				W=HYPER[i].W;
				//En+=0.5/mi*(p[A_X]*p[A_X]+p[A_Y]*p[A_Y]+p[A_Z]*p[A_Z])+V*W+HYPER[i].mu*g;
				En+=0.5/mi*(p[A_X]*p[A_X]+p[A_Y]*p[A_Y]+p[A_Z]*p[A_Z])+V*W;
			}

			T=0.;
			for(int k=0;k<Nx;k++)
			{	
				dT[k]=0.;
				for(int l=0;l<Nx;l++)	rT[l*Nx+k]=0.;
			}


			T+=(En-E0)*(En-E0);
			for(int k=0;k<Nx;k++)
			{	
				dT[k]+=2.*dE[k]*(En-E0);
				for(int l=0;l<Nx;l++)	rT[l*Nx+k]+=2.*rE[l*Nx+k]*(En-E0)+2*dE[k]*dE[l];
			}
			for(int i=0;i<h_num;i++)
			{	
				if(G[i]<0)	T+=0.5*rg*(-1.*G[i]+th_G[i])*(-1.*G[i]+th_G[i]);
				else
				{
					T+=0.5*rg*(G[i]+th_G[i])*(G[i]+th_G[i]);
				}
			}
			for(int i=0;i<Nw;i++)	if(H[i]+th_H[i]>0)	T+=0.5*rh*(H[i]+th_H[i])*(H[i]+th_H[i]);

			for(int k=0;k<Nx;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					if(G[i]<0)	dT[k]+=-1.*rg*dG[i][k]*(-1.*G[i]+th_G[i]);
					else
					{
						dT[k]+=rg*dG[i][k]*(G[i]+th_G[i]);
					}
				}
				for(int i=0;i<Nw;i++)
				{	
					if(H[i]+th_H[i]>0)	dT[k]+=rh*dH[i][k]*(H[i]+th_H[i]);
				}


				for(int l=0;l<Nx;l++)
				{
					for(int i=0;i<h_num;i++)
					{	
						//rT[l*h_num+k]+=r*2.*dG[i*h_num+k]*dG[i*h_num+l]*(3*G[i]*G[i]+th_G[i]);
						rT[l*Nx+k]+=rg*dG[i][k]*dG[i][l];
						//rT[l*h_num+k]+=r*dG[i*h_num+k]*dG[i*h_num+l];
					}
					for(int i=0;i<Nw;i++)
					{	
						//rT[l*h_num+k]+=r*2.*dG[i*h_num+k]*dG[i*h_num+l]*(3*G[i]*G[i]+th_G[i]);
						if(H[i]+th_H[i]>0)	rT[l*Nx+k]+=rh*dH[i][k]*dH[i][l];
						else if (H[i]+th_H[i]>-1.e-20 && dH[i][k]*dH[i][l]>0)	rT[l*Nx+k]+=rh*dH[i][k]*dH[i][l];
						//rT[l*h_num+k]+=r*dG[i*h_num+k]*dG[i*h_num+l];
					}
				}
			}

			d_sum=0.;
			for(int i=0;i<Nx;i++)
			{
				d[i]=dT[i];
				d_sum+=fabs(dT[i]);
				for(int j=0;j<Nx;j++)	B[i*Nx+j]=rT[i*Nx+j];
			}
			if(d_sum<1.e-20)	break;
			gauss(B,d,Nx);

			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=d[i];
				E_sum+=fabs(d[i]);
			}
			for(int i=0;i<Nw;i++)
			{
				int iw=NEIw[i];
				HYPER[iw].mu-=d[i+h_num];
				E_sum+=fabs(d[i+h_num]);
			}
			E=E_sum;
			//cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;

			p_variables(HYPER,HYPER1,h_num,Dt,mi);
			if(count==1||count%10000==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				output_data_p(HYPER,HYPER1,dG,G,th_G,dH,H,th_H,rT,dT,T,d,count,count_min,t,E,En,E0,NEIw);
			}

		}
		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=fabs(HYPER[i].old_lam-HYPER[i].lam);
		for(int i=0;i<Nw;i++)
		{
			int iw=NEIw[i];
			E_sum+=fabs(HYPER[iw].old_mu-HYPER[iw].mu);
		}
		E_min=E_sum;

		//if(E_min<ep_min*1000)
		//{
		//	rg*=4;
		//	rh*=4;
		//}


		for(int i=0;i<h_num;i++)
		{
			if(G[i]<0.)	th_G[i]+=-1.*G[i];
			else
			{
				th_G[i]+=G[i];
			}
			//th_g[i]+=g[i];
			//if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		}
		for(int i=0;i<Nw;i++)	if(H[i]+th_H[i]>0.)	th_H[i]+=H[i];

		//////ÉfÅ[É^èoóÕ
		//if(count_min==1||count_min%100==0)
		{
			cout<<"Emin"<<count_min<<"="<<E_min<<endl;
			stringstream ss0;
			ss0<<"./E_min/E"<<t<<".csv"<<endl;
			stringstream ss1;
			ss1<<"./p/p"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss3;
			ss3<<"./mu/mu"<<t<<"_"<<count_min<<".csv"<<endl;
			stringstream ss4;
			ss4<<"./p_T/T"<<t<<".csv"<<endl;

			if(count_min==1)
			{
				ofstream fem(ss0.str(), ios::trunc);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::trunc);
				ft<<T<<","<<En<<","<<E0<<endl;
				ft.close();
			}
			else
			{
				ofstream fem(ss0.str(), ios::app);
				fem<<E_min<<endl;
				fem.close();
				ofstream ft(ss4.str(), ios::app);
				ft<<T<<","<<En<<endl;
				ft.close();
			}

			ofstream fp(ss1.str());
			ofstream fl(ss3.str());
			for(int i=0;i<h_num;i++)
			{
				fp<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
				fl<<HYPER[i].lam<<","<<HYPER[i].mu<<endl;
			}
			fp.close();
			fl.close();
		}


		if(count_min>c_max)	break;
	}



	for(int i=0;i<h_num;i++)	delete[]	dG[i];
	for(int i=0;i<Nw;i++)	delete[]	dH[i];

	delete[]	dE;
	delete[]	rE;
	delete[]	G;
	delete[]	dG;
	delete[]	th_G;
	delete[]	H;
	delete[]	dH;
	delete[]	th_H;
	delete[]	dT;
	delete[]	rT;
	delete[]	d;
	delete[]	B;
}


void p_nab_nw1(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dE, double *G, double *H,double Dt, double V, double mi,vector<double > NEIw)
{

	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=h_num+Nw;
	/////////////////èâä˙âª///////////////////
	for(int i=0;i<h_num;i++)	G[i]=0.;
	for(int i=0;i<Nw;i++)	H[i]=0.;
	for(int i=0;i<Nx;i++)	dE[i]=0.;


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

			dE[k]-=0.5*Dt/mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
			G[k]+=1./mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
		}
		dE[k]-=0.5*Dt/mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
		//dE[k]+=V*(1.-HYPER[k].J);
		G[k]+=1./mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
	}
	for(int k=0;k<Nw;k++)
	{
		int kw=NEIw[k];
		pk[A_Z]=HYPER[kw].p[A_Z];
		dE[k+h_num]=0.5*Dt/mi*pk[A_Z];
		H[k]=-1./mi*pk[A_Z];
	}

}


void p_lap_nw1(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *rE, double **dG, double **dH, double Dt, double V, double mi,vector<double > NEIw)
{

	int h_num=HYPER.size();
	int Nw=NEIw.size();
	int Nx=Nw+h_num;

	/////////////////èâä˙âª///////////////////
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nx;j++)	rE[i*Nx+j]=0.;
		for(int j=0;j<Nw;j++)	dH[j][i]=0.;
		for(int j=0;j<h_num;j++)	dG[j][i]=0.;

	}


	//////////íËã`
	double Dgjj[DIMENSION]={0,0,0};
	double Dgkj[DIMENSION]={0,0,0};
	double Dglj[DIMENSION]={0,0,0};
	for(int j=0;j<h_num;j++)
	{
		int Nj=HYPER[j].N;
		Dgjj[A_X]=HYPER1[j*h_num+j].DgDq[A_X];	Dgjj[A_Y]=HYPER1[j*h_num+j].DgDq[A_Y];	Dgjj[A_Z]=HYPER1[j*h_num+j].DgDq[A_Z];

		for(int k=0;k<Nj;k++)
		{
			int kn=HYPER[j].NEI[k];

			Dgkj[A_X]=HYPER1[kn*h_num+j].DgDq[A_X];	Dgkj[A_Y]=HYPER1[kn*h_num+j].DgDq[A_Y];	Dgkj[A_Z]=HYPER1[kn*h_num+j].DgDq[A_Z];
			for(int l=0;l<Nj;l++)
			{
				int ln=HYPER[j].NEI[l];
				/////////////rlam E
				Dglj[A_X]=HYPER1[ln*h_num+j].DgDq[A_X];	Dglj[A_Y]=HYPER1[ln*h_num+j].DgDq[A_Y];	Dglj[A_Z]=HYPER1[ln*h_num+j].DgDq[A_Z];
				rE[ln*Nx+kn]+=0.25*Dt*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
				/////////////dlam G
				dG[kn][ln]-=0.5*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
			}
			/////////////rlam E
			rE[j*Nx+kn]+=0.25*Dt*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			rE[kn*Nx+j]+=0.25*Dt*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

			/////////////dlam G
			dG[kn][j]-=0.5*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			dG[j][kn]-=0.5*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

		}
		/////////////rlam E
		rE[j*Nx+j]+=0.25*Dt*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);
		/////////////dlam G
		dG[j][j]-=0.5*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);
	}

	//////////íËã`
	double Dgki[DIMENSION]={0,0,0};
	double Dgii[DIMENSION]={0,0,0};
	for(int i=0;i<Nw;i++)
	{
		int iw=NEIw[i];
		int Ni=HYPER[iw].N;
		Dgii[A_Z]=HYPER1[iw*h_num+iw].DgDq[A_Z];

		for(int k=0;k<Ni;k++)
		{
			int kn=HYPER[iw].NEI[k];
			Dgki[A_Z]=HYPER1[kn*h_num+iw].DgDq[A_Z];
			/////////////rmulam E
			rE[(i+h_num)*Nx+kn]-=0.25*Dt*Dt/mi*Dgki[A_Z];
			/////////////rlammu E
			rE[kn*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgki[A_Z];

			/////////////dmu G
			dG[kn][i+h_num]+=0.5*Dt/mi*Dgki[A_Z];
			/////////////dlam H
			dH[i][kn]+=0.5*Dt/mi*Dgki[A_Z];
		}
		/////////////rmulam E
		rE[(i+h_num)*Nx+iw]-=0.25*Dt*Dt/mi*Dgii[A_Z];
		/////////////rlammu E
		rE[iw*Nx+i+h_num]-=0.25*Dt*Dt/mi*Dgii[A_Z];
		/////////////rmu E
		rE[(i+h_num)*Nx+i+h_num]=0.25*Dt*Dt/mi;
		/////////////dmu G
		dG[iw][i+h_num]+=0.5*Dt/mi*Dgii[A_Z];
		/////////////dlam H
		dH[i][iw]+=0.5*Dt/mi*Dgii[A_Z];
		/////////////dmu H
		dH[i][i+h_num]=-0.5*Dt/mi;
	}





}


void output_data_p_nw1(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double **dG, double *G, double *th_G, double **dH, double *H, double *th_H,double *rT, double *dT, double T, double *d, int count,int count_min,int t,double E,double En, double E0,vector<double > NEIw)
{
	int Nw=NEIw.size();
	int h_num=HYPER.size();
	int Nx=Nw+h_num;


	stringstream ss0;
	ss0<<"./p_G/dG_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss1;
	ss1<<"./p_G/G_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss2;
	ss2<<"./p_T/rT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss3;
	ss3<<"./p_T/dT_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss4;
	ss4<<"./p_T/T_"<<t<<"_"<<count_min<<".csv";
	stringstream ss5;
	ss5<<"./p/p_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss6;
	ss6<<"./E/p_E_"<<t<<"_"<<count_min<<".csv";
	stringstream ss7;
	ss7<<"./p_T/d_"<<t<<"_"<<count_min<<"_"<<count<<".csv";

	stringstream ss8;
	ss8<<"./p_H/dH_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss9;
	ss9<<"./p_H/H_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss10;
	ss10<<"./lam/p_lam_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss11;
	ss11<<"./mu/p_mu_"<<t<<"_"<<count_min<<"_"<<count<<".csv";


	if(count==1)
	{
		ofstream f_t(ss4.str(),ios::trunc);
		f_t<<T<<","<<En<<","<<E0<<endl;
		f_t.close();
		ofstream f_e(ss6.str(),ios::trunc);
		f_e<<E<<endl;
		f_e.close();
	}
	else
	{
		ofstream f_t(ss4.str(),ios::app);
		f_t<<T<<","<<En<<endl;
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
	ofstream f_lam(ss10.str());
	ofstream f_mu(ss11.str());


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
	for(int i=0;i<h_num;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_dG<<dG[i][j]<<",";
		}
		f_dG<<endl;
		f_G<<G[i]<<","<<th_G[i]<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		f_lam<<HYPER[i].lam<<endl;
		f_mu<<HYPER[i].mu<<endl;
	}
	for(int i=0;i<Nw;i++)
	{				
		for(int j=0;j<Nx;j++)
		{
			f_dH<<dH[i][j]<<",";
		}
		f_dH<<endl;
		f_H<<H[i]<<","<<th_H[i]<<endl;
	}

	f_dG.close();
	f_G.close();
	f_rt.close();
	f_dt.close();
	f_p.close();
	f_d.close();
	f_dH.close();
	f_H.close();
	f_lam.close();
	f_mu.close();

}


