#include "stdafx.h"	

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0);
void p_QP(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0);
void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,double nG[DIMENSION]);
void q_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT,double *rL,double *dL,int h_num,int count,int count_min,int t,double E);

void calc_HYPER_QP_gh(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{
	////////////’è‹`///////////////
	int h_num=HYPER.size();
	int Nx=h_num*2;

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
	double nG[DIMENSION]={0,0,1};

	double E0=0;

	cout<<"QP start-------------------------"<<endl;
	for(int i=0;i<h_num;i++)	E0+=0.5/mi*(HYPER[i].p_n[A_X]*HYPER[i].p_n[A_X]+HYPER[i].p_n[A_Y]*HYPER[i].p_n[A_Y]+HYPER[i].p_n[A_Z]*HYPER[i].p_n[A_Z])+V*HYPER[i].W_n;
	cout<<"E0="<<E0<<endl;
	ofstream fq0("q0_QP.csv");
	ofstream fqn0("qn0_QP.csv");
	ofstream fph0("ph0_QP.csv");
	ofstream fphn0("phn0_QP.csv");
	for(int i=0;i<h_num;i++)
	{
		fq0<<HYPER[i].q_n[A_X]<<","<<HYPER[i].q_n[A_Y]<<","<<HYPER[i].q_n[A_Z]<<endl;
		fqn0<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		fph0<<HYPER[i].ph_n[A_X]<<","<<HYPER[i].ph_n[A_Y]<<","<<HYPER[i].ph_n[A_Z]<<endl;
		fphn0<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
	}
	fq0.close();
	fqn0.close();
	fph0.close();
	fphn0.close();

	////////////QPŒvŽZ///////////////		
	q_QP(CON,PART,HYPER,HYPER1,h_num,Nx,V,mi,Dt,nG,E0);

	ofstream fqn("qn_QP.csv");
	ofstream fphn("phn_QP.csv");
	for(int i=0;i<h_num;i++)
	{
		fqn<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
		fphn<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
	}
	fqn.close();
	fphn.close();

	p_QP(PART,HYPER,HYPER1,h_num,Nx,V,mi,Dt,nG,E0);

	cout<<"--------------------------OK"<<endl;
	

}

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0)
{
	////////////’è‹`///////////////
	int count=0;
	int count_min=0;
	int c_max=10000;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1e-5;
	double ep_min=1e-5;

	double r=100000;

	double En=0;
	double T=0;

	double *dE=new double [Nx];
	double *dT=new double [Nx];	

	double *rE=new double [Nx*Nx];
	double *rT=new double [Nx*Nx];

	double *g=new double [h_num];
	double *h=new double [h_num];

	double **dg=new double *[h_num];
	double **dh=new double *[h_num];

	for(int i=0;i<h_num;i++)
	{
		dg[i]=new double [Nx];
		dh[i]=new double [Nx];
	}

	double *th_g=new double [h_num];
	double *th_h=new double [h_num];
		
	////////////‰Šú‰»ŽZ///////////////
	for(int i=0;i<h_num;i++)
	{
		dE[i]=0;	dE[i+h_num]=0;
		dT[i]=0;	dT[i+h_num]=0;
		g[i]=0;
		h[i]=0;
		th_g[i]=0;
		th_h[i]=0;
		for(int j=0;j<Nx;j++)
		{
			dg[i][j]=0;
			dh[i][j]=0;
			rE[i*Nx+j]=0;	rE[(i+h_num)*Nx+j]=0;
			rT[i*Nx+j]=0;	rT[(i+h_num)*Nx+j]=0;
		}
	}

	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			dh[in][k]=0.5*Dt*Dt/mi*(HYPER1[k*h_num+i].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*nG[A_Z]);
		}
		dh[k][k+h_num]=-0.5*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
	}



	while(E_min>ep_min)
	{
		count_min++;
	
		for(int i=0;i<h_num;i++)
		{
			HYPER[i].old_lam=HYPER[i].lam;
			HYPER[i].old_mu=HYPER[i].mu;
		}
		
		E=1;
		count=0;
		while(E>ep)
		{
			count++;

			q_variables(CON,PART,HYPER,HYPER1);
			q_nab_lap(PART,HYPER,HYPER1,dE,rE,dg,Dt,V,mi,nG);

			for(int i=0;i<h_num;i++)	En=0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])+V*HYPER[i].W;

			T=(En-E0)*(En-E0);
			//cout<<"T="<<T<<", En"<<En<<endl;
			for(int k=0;k<Nx;k++)
			{	
				dT[k]=2.*dE[k]*(E-E0);
				for(int l=0;l<Nx;l++)	rT[l*Nx+k]=2.*rE[l*Nx+k]*(E-E0)+2*dE[k]*dE[l];
			}
			//for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			for(int i=0;i<h_num;i++)
			{	
				g[i]=V*(1-HYPER[i].J);
				h[i]=-1*(PART[i].r[A_X]*nG[A_X]+PART[i].r[A_Y]*nG[A_Y]+PART[i].r[A_Z]*nG[A_Z]);
				
				T+=0.5*r*(g[i]*g[i]+th_g[i])*(g[i]*g[i]+th_g[i]);
				if(h[i]+th_h[i]>0)	T+=0.5*r*(h[i]+th_h[i])*(h[i]+th_h[i]);
			}

			for(int k=0;k<Nx;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					dT[k]+=r*2.*dg[i][k]*g[i]*(g[i]*g[i]+th_g[i]);
					//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
					if(h[i]+th_h[i]>0)	dT[k]+=r*dh[i][k]*(h[i]+th_h[i]);
				}
				for(int l=0;l<Nx;l++)
				{
					for(int i=0;i<h_num;i++)
					{	
						rT[l*Nx+k]+=2.*r*g[i]*dg[i][k]*dg[i][l]*(3.*g[i]*g[i]+th_g[i]);
						if(h[i]+th_h[i]>0)	rT[l*Nx+k]+=r*dh[i][k]*dh[i][l];
						else if(h[i]+th_h[i]==0 && dh[i][k]*dh[i][l]>0)
						{
							rT[l*Nx+k]+=r*dh[i][k]*dh[i][l];
						}
					}
				}
			}
			//for(int k=0;k<Nx;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			if(count==1||count%250==0)
			{

				stringstream sg;
				sg<<"./lam/g"<<count_min<<"_"<<count<<".csv";
				ofstream fg(sg.str());
				stringstream sdg;
				sdg<<"./lam/dg"<<count_min<<"_"<<count<<".csv";
				ofstream fdg(sdg.str());

				stringstream sh;
				sh<<"./mu/h"<<count_min<<"_"<<count<<".csv";
				ofstream fh(sh.str());
				stringstream sdh;
				sdh<<"./mu/dh"<<count_min<<"_"<<count<<".csv";
				ofstream fdh(sdh.str());

				for(int i=0;i<h_num;i++)
				{
					fg<<g[i]<<","<<th_g[i]<<endl;
					fh<<h[i]<<","<<th_h[i]<<endl;
					for(int j=0;j<Nx;j++)
					{
						fdg<<dg[i][j]<<",";
						fdh<<dh[i][j]<<",";
					}
					fdg<<endl;
					fdh<<endl;
				}
				fdg.close();
				fdh.close();

				stringstream std;
				std<<"./dT/dT"<<count_min<<"_"<<count<<".csv";
				ofstream fd(std.str());
				stringstream sTr;
				sTr<<"./rT/rT"<<count_min<<"_"<<count<<".csv";
				ofstream fr(sTr.str());

				stringstream sed;
				sed<<"./dL/dE"<<count_min<<"_"<<count<<".csv";
				ofstream fde(sed.str());
				stringstream str;
				str<<"./rL/rE"<<count_min<<"_"<<count<<".csv";
				ofstream fre(str.str());

				for(int i=0;i<Nx;i++)
				{
					fd<<dT[i]<<endl;
					fde<<dE[i]<<endl;
					for(int j=0;j<Nx;j++)
					{
						fr<<rT[i*Nx+j]<<",";
						fre<<rE[i*Nx+j]<<",";
					}
					fr<<endl;
					fre<<endl;
				}
				fd.close();
				fr.close();
				fde.close();
				fre.close();

			}
			gauss(rT,dT,Nx);
			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam+=-dT[i];
				HYPER[i].mu+=-dT[i+h_num];
				E_sum+=dT[i]*dT[i]+dT[i+h_num]*dT[i+h_num];
			}
			E=sqrt(E_sum);
			if(count==1||count>600)	cout<<"E"<<count<<"="<<E<<endl;

			/////////////pn1_2ŒvŽZ
			for(int i=0;i<h_num;i++)
			{
				double p_p[DIMENSION]={0,0,0};
				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{
					int jn=HYPER[i].NEI[j];
					p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
					p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
					p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_Z];
				}
				p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
				p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
				p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Z];
				HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X]+0.5*Dt*nG[A_X]*HYPER[i].mu;
				HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y]+0.5*Dt*nG[A_Y]*HYPER[i].mu;
				HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*nG[A_Z]*HYPER[i].mu;
	
				/////////////qŒvŽZ
				PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
				PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
				PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
			}
			if(count==1||count%250==0)
			{

			stringstream stpnh;
			stpnh<<"./p/phn_QP"<<count_min<<"_"<<count<<".csv";
			ofstream fpnh(stpnh.str());
			stringstream stq;
			stq<<"./q/qn_QP"<<count_min<<"_"<<count<<".csv";
			ofstream fq(stq.str());
			for(int i=0;i<h_num;i++)
			{
				fpnh<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
				fq<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
			}
			fpnh.close();
			fq.close();
			}

		}

		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=(HYPER[i].old_lam-HYPER[i].lam)*(HYPER[i].old_lam-HYPER[i].lam)+(HYPER[i].old_mu-HYPER[i].mu)*(HYPER[i].old_mu-HYPER[i].mu);
		E_min=sqrt(E_sum);
		if(E_min<ep_min*1000)	r*=4;
		cout<<"Emin"<<count_min<<"="<<E_min<<endl;

		for(int i=0;i<h_num;i++)
		{
			th_g[i]+=g[i]*g[i];
			if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		}
	}

	delete[]	dE;
	delete[]	dT;
	delete[]	rE;
	delete[]	rT;
	delete[]	g;
	delete[]	h;
	for(int i=0;i<h_num;i++)
	{
		delete[]	dg[i];
		delete[]	dh[i];
	}
	delete[]	dg;
	delete[]	dh;


	/////////////pn1_2ŒvŽZ
	for(int i=0;i<h_num;i++)
	{
		double p_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_Z];
		}
		p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Z];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X]+0.5*Dt*nG[A_X]*HYPER[i].mu;
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y]+0.5*Dt*nG[A_Y]*HYPER[i].mu;
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*nG[A_Z]*HYPER[i].mu;
	
		/////////////qŒvŽZ
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


	/////////////pn1_2ŒvŽZ
	for(int i=0;i<h_num;i++)
	{
		double p_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int jn=HYPER[i].NEI[j];
			p_p[A_X]+=(HYPER[jn].stress_n[A_X][A_X]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_X][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_X][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_p[A_Y]+=HYPER[jn].stress_n[A_Y][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+(HYPER[jn].stress_n[A_Y][A_Y]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_Y]+HYPER[jn].stress_n[A_Y][A_Z]*HYPER1[jn*h_num+i].DgDq_n[A_Z];
			p_p[A_Z]+=HYPER[jn].stress_n[A_Z][A_X]*HYPER1[jn*h_num+i].DgDq_n[A_X]+HYPER[jn].stress_n[A_Z][A_Y]*HYPER1[jn*h_num+i].DgDq_n[A_Y]+(HYPER[jn].stress_n[A_Z][A_Z]-HYPER[jn].lam)*HYPER1[jn*h_num+i].DgDq_n[A_Z];
		}
		p_p[A_X]+=(HYPER[i].stress_n[A_X][A_X]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_X][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_X][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_p[A_Y]+=HYPER[i].stress_n[A_Y][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+(HYPER[i].stress_n[A_Y][A_Y]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Y]+HYPER[i].stress_n[A_Y][A_Z]*HYPER1[i*h_num+i].DgDq_n[A_Z];
		p_p[A_Z]+=HYPER[i].stress_n[A_Z][A_X]*HYPER1[i*h_num+i].DgDq_n[A_X]+HYPER[i].stress_n[A_Z][A_Y]*HYPER1[i*h_num+i].DgDq_n[A_Y]+(HYPER[i].stress_n[A_Z][A_Z]-HYPER[i].h_lam)*HYPER1[i*h_num+i].DgDq_n[A_Z];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_p[A_X]+0.5*Dt*nG[A_X]*HYPER[i].mu;
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_p[A_Y]+0.5*Dt*nG[A_Y]*HYPER[i].mu;
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_p[A_Z]+0.5*Dt*nG[A_Z]*HYPER[i].mu;
	
		/////////////qŒvŽZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}


	/////////////F, J, t_inverse‚ÌŒvŽZ
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


	/////////////DgDq, Stress, W, S, dSdcŒvŽZ
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double **in_Ci=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	in_Ci[D]=new double[DIMENSION];
	double in_Ci2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double first_term[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double first_term2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double sec_term[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double S[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dSdc[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	double dPIdF[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<h_num;j++)
	{	

		////////WŒvŽZ
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

	
		
		////////S, dSdcŒvŽZ
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

		first_term[0][0]=1.-1/3.*Ic*in_Ci[0][0];	first_term[0][1]=-1/3.*Ic*in_Ci[0][1];	first_term[0][2]=-1/3.*Ic*in_Ci[0][2];
		first_term[1][0]=-1/3.*Ic*in_Ci[1][0];	first_term[1][1]=1.-1/3.*Ic*in_Ci[1][1];	first_term[1][2]=-1/3.*Ic*in_Ci[1][0];
		first_term[2][0]=-1/3.*Ic*in_Ci[2][0];	first_term[2][1]=-1/3.*Ic*in_Ci[2][1];	first_term[2][2]=1.-1/3.*Ic*in_Ci[2][2];

		first_term2[0][0]=first_term[0][0]*first_term[0][0]+first_term[0][1]*first_term[1][0]+first_term[0][2]*first_term[2][0];
		first_term2[0][1]=first_term[0][0]*first_term[0][1]+first_term[0][1]*first_term[1][1]+first_term[0][2]*first_term[2][1];
		first_term2[0][2]=first_term[0][0]*first_term[0][2]+first_term[0][1]*first_term[1][2]+first_term[0][2]*first_term[2][2];
		first_term2[1][0]=first_term[1][0]*first_term[0][0]+first_term[1][1]*first_term[1][0]+first_term[1][2]*first_term[2][0];
		first_term2[1][1]=first_term[1][0]*first_term[0][1]+first_term[1][1]*first_term[1][1]+first_term[1][2]*first_term[2][1];
		first_term2[1][2]=first_term[1][0]*first_term[0][2]+first_term[1][1]*first_term[1][2]+first_term[1][2]*first_term[2][2];
		first_term2[2][0]=first_term[2][0]*first_term[0][0]+first_term[2][1]*first_term[1][0]+first_term[2][2]*first_term[2][0];
		first_term2[2][1]=first_term[2][0]*first_term[0][1]+first_term[2][1]*first_term[1][1]+first_term[2][2]*first_term[2][1];
		first_term2[2][2]=first_term[2][0]*first_term[0][2]+first_term[2][1]*first_term[1][2]+first_term[2][2]*first_term[2][2];

		sec_term[0][0]=in_Ci[0][0]*first_term[0][0]+in_Ci[0][1]*first_term[1][0]+in_Ci[0][2]*first_term[2][0];
		sec_term[0][1]=in_Ci[0][0]*first_term[0][1]+in_Ci[0][1]*first_term[1][1]+in_Ci[0][2]*first_term[2][1];
		sec_term[0][2]=in_Ci[0][0]*first_term[0][2]+in_Ci[0][1]*first_term[1][2]+in_Ci[0][2]*first_term[2][2];
		sec_term[1][0]=in_Ci[1][0]*first_term[0][0]+in_Ci[1][1]*first_term[1][0]+in_Ci[1][2]*first_term[2][0];
		sec_term[1][1]=in_Ci[1][0]*first_term[0][1]+in_Ci[1][1]*first_term[1][1]+in_Ci[1][2]*first_term[2][1];
		sec_term[1][2]=in_Ci[1][0]*first_term[0][2]+in_Ci[1][1]*first_term[1][2]+in_Ci[1][2]*first_term[2][2];
		sec_term[2][0]=in_Ci[2][0]*first_term[0][0]+in_Ci[2][1]*first_term[1][0]+in_Ci[2][2]*first_term[2][0];
		sec_term[2][1]=in_Ci[2][0]*first_term[0][1]+in_Ci[2][1]*first_term[1][1]+in_Ci[2][2]*first_term[2][1];
		sec_term[2][2]=in_Ci[2][0]*first_term[0][2]+in_Ci[2][1]*first_term[1][2]+in_Ci[2][2]*first_term[2][2];

		if(J<0)
		{
			S[0][0]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3.*trace_dC2*in_Ci[0][0]) );
			S[0][1]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3.*trace_dC2*in_Ci[0][1]) );
			S[0][2]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3.*trace_dC2*in_Ci[0][2]) );
			S[1][0]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3.*trace_dC2*in_Ci[1][0]) );
			S[1][1]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3.*trace_dC2*in_Ci[1][1]) );
			S[1][2]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3.*trace_dC2*in_Ci[1][2]) );
			S[2][0]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3.*trace_dC2*in_Ci[2][0]) );
			S[2][1]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3.*trace_dC2*in_Ci[2][1]) );		
			S[2][2]=-2.*1/pow(-J,2/3.)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3.*trace_dC2*in_Ci[2][2]) );

			dSdc[0][0]=-2.*1/pow(-J,4/3.)*( c01*first_term2[0][0]-(c10+c01*Ic)*sec_term[0][0]+1/3.*c01*(1.-5/3.*trace_dC2*in_Ci2[0][0]) );
			dSdc[0][1]=-2.*1/pow(-J,4/3.)*( c01*first_term2[0][1]-(c10+c01*Ic)*sec_term[0][1]-5/9.*c01*trace_dC2*in_Ci2[0][1] );
			dSdc[0][2]=-2.*1/pow(-J,4/3.)*( c01*first_term2[0][2]-(c10+c01*Ic)*sec_term[0][2]-5/9.*c01*trace_dC2*in_Ci2[0][2] );
			dSdc[1][0]=-2.*1/pow(-J,4/3.)*( c01*first_term2[1][0]-(c10+c01*Ic)*sec_term[1][0]-5/9.*c01*trace_dC2*in_Ci2[1][0] );
			dSdc[1][1]=-2.*1/pow(-J,4/3.)*( c01*first_term2[1][1]-(c10+c01*Ic)*sec_term[1][1]+1/3.*c01*(1.-5/3.*trace_dC2*in_Ci2[1][1]) );
			dSdc[1][2]=-2.*1/pow(-J,4/3.)*( c01*first_term2[1][2]-(c10+c01*Ic)*sec_term[1][2]-5/9.*c01*trace_dC2*in_Ci2[1][2] );
			dSdc[2][0]=-2.*1/pow(-J,4/3.)*( c01*first_term2[2][0]-(c10+c01*Ic)*sec_term[2][0]-5/9.*c01*trace_dC2*in_Ci2[2][0] );
			dSdc[2][1]=-2.*1/pow(-J,4/3.)*( c01*first_term2[2][1]-(c10+c01*Ic)*sec_term[2][1]-5/9.*c01*trace_dC2*in_Ci2[2][1] );		
			dSdc[2][2]=-2.*1/pow(-J,4/3.)*( c01*first_term2[2][2]-(c10+c01*Ic)*sec_term[2][2]+1/3.*c01*(1.-5/3.*trace_dC2*in_Ci2[2][2]) );
		}
		else
		{
			S[0][0]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][0]-c01*(dC[0][0]-1/3.*trace_dC2*in_Ci[0][0]) );
			S[0][1]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][1]-c01*(dC[0][1]-1/3.*trace_dC2*in_Ci[0][1]) );
			S[0][2]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[0][2]-c01*(dC[0][2]-1/3.*trace_dC2*in_Ci[0][2]) );
			S[1][0]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][0]-c01*(dC[1][0]-1/3.*trace_dC2*in_Ci[1][0]) );
			S[1][1]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][1]-c01*(dC[1][1]-1/3.*trace_dC2*in_Ci[1][1]) );
			S[1][2]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[1][2]-c01*(dC[1][2]-1/3.*trace_dC2*in_Ci[1][2]) );
			S[2][0]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][0]-c01*(dC[2][0]-1/3.*trace_dC2*in_Ci[2][0]) );
			S[2][1]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][1]-c01*(dC[2][1]-1/3.*trace_dC2*in_Ci[2][1]) );		
			S[2][2]=2.*1/pow(J,2/3.)*( (c10+c01*Ic)*first_term[2][2]-c01*(dC[2][2]-1/3.*trace_dC2*in_Ci[2][2]) );

			dSdc[0][0]=2.*1/pow(J,4/3.)*( c01*first_term2[0][0]-(c10+c01*Ic)*sec_term[0][0]+1/3.*c01*(1.-5/3.*trace_dC2*in_Ci2[0][0]) );
			dSdc[0][1]=2.*1/pow(J,4/3.)*( c01*first_term2[0][1]-(c10+c01*Ic)*sec_term[0][1]-5/9.*c01*trace_dC2*in_Ci2[0][1] );
			dSdc[0][2]=2.*1/pow(J,4/3.)*( c01*first_term2[0][2]-(c10+c01*Ic)*sec_term[0][2]-5/9.*c01*trace_dC2*in_Ci2[0][2] );
			dSdc[1][0]=2.*1/pow(J,4/3.)*( c01*first_term2[1][0]-(c10+c01*Ic)*sec_term[1][0]-5/9.*c01*trace_dC2*in_Ci2[1][0] );
			dSdc[1][1]=2.*1/pow(J,4/3.)*( c01*first_term2[1][1]-(c10+c01*Ic)*sec_term[1][1]+1/3.*c01*(1.-5/3.*trace_dC2*in_Ci2[1][1]) );
			dSdc[1][2]=2.*1/pow(J,4/3.)*( c01*first_term2[1][2]-(c10+c01*Ic)*sec_term[1][2]-5/9.*c01*trace_dC2*in_Ci2[1][2] );
			dSdc[2][0]=2.*1/pow(J,4/3.)*( c01*first_term2[2][0]-(c10+c01*Ic)*sec_term[2][0]-5/9.*c01*trace_dC2*in_Ci2[2][0] );
			dSdc[2][1]=2.*1/pow(J,4/3.)*( c01*first_term2[2][1]-(c10+c01*Ic)*sec_term[2][1]-5/9.*c01*trace_dC2*in_Ci2[2][1] );		
			dSdc[2][2]=2.*1/pow(J,4/3.)*( c01*first_term2[2][2]-(c10+c01*Ic)*sec_term[2][2]+1/3.*c01*(1.-5/3.*trace_dC2*in_Ci2[2][2]) );
		}
		dPIdF[0][0]=S[0][0]+2.*dSdc[0][0]*HYPER[j].Fi[0][0]					 +dSdc[0][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[0][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		dPIdF[0][1]=S[0][1]+dSdc[0][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2.*dSdc[0][1]*HYPER[j].Fi[1][1]				  +dSdc[0][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		dPIdF[0][2]=S[0][2]+dSdc[0][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[0][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2.*dSdc[0][2]*HYPER[j].Fi[2][2];
		dPIdF[1][0]=S[1][0]+2.*dSdc[1][0]*HYPER[j].Fi[0][0]					 +dSdc[1][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[1][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		dPIdF[1][1]=S[1][1]+dSdc[1][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2.*dSdc[1][1]*HYPER[j].Fi[1][1]				  +dSdc[1][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		dPIdF[1][2]=S[1][2]+dSdc[1][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[1][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2.*dSdc[1][2]*HYPER[j].Fi[2][2];
		dPIdF[2][0]=S[2][0]+2.*dSdc[2][0]*HYPER[j].Fi[0][0]					 +dSdc[2][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[2][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		dPIdF[2][1]=S[2][1]+dSdc[2][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2.*dSdc[2][1]*HYPER[j].Fi[1][1]				  +dSdc[2][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		dPIdF[2][2]=S[2][2]+dSdc[2][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[2][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2.*dSdc[2][2]*HYPER[j].Fi[2][2];

		HYPER[j].pi[0][0]=HYPER[j].Fi[0][0]*S[0][0]+HYPER[j].Fi[0][1]*S[1][0]+HYPER[j].Fi[0][2]*S[2][0];
		HYPER[j].pi[0][1]=HYPER[j].Fi[0][0]*S[0][1]+HYPER[j].Fi[0][1]*S[1][1]+HYPER[j].Fi[0][2]*S[2][1];
		HYPER[j].pi[0][2]=HYPER[j].Fi[0][0]*S[0][2]+HYPER[j].Fi[0][1]*S[1][2]+HYPER[j].Fi[0][2]*S[2][2];
		HYPER[j].pi[1][0]=HYPER[j].Fi[1][0]*S[0][0]+HYPER[j].Fi[1][1]*S[1][0]+HYPER[j].Fi[1][2]*S[2][0];
		HYPER[j].pi[1][1]=HYPER[j].Fi[1][0]*S[0][1]+HYPER[j].Fi[1][1]*S[1][1]+HYPER[j].Fi[1][2]*S[2][1];
		HYPER[j].pi[1][2]=HYPER[j].Fi[1][0]*S[0][2]+HYPER[j].Fi[1][1]*S[1][2]+HYPER[j].Fi[1][2]*S[2][2];
		HYPER[j].pi[2][0]=HYPER[j].Fi[2][0]*S[0][0]+HYPER[j].Fi[2][1]*S[1][0]+HYPER[j].Fi[2][2]*S[2][0];
		HYPER[j].pi[2][1]=HYPER[j].Fi[2][0]*S[0][1]+HYPER[j].Fi[2][1]*S[1][1]+HYPER[j].Fi[2][2]*S[2][1];
		HYPER[j].pi[2][2]=HYPER[j].Fi[2][0]*S[0][2]+HYPER[j].Fi[2][1]*S[1][2]+HYPER[j].Fi[2][2]*S[2][2];

		HYPER[j].dpidF[0][0]=S[0][0]+2.*dSdc[0][0]*HYPER[j].Fi[0][0]					 +dSdc[0][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[0][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		HYPER[j].dpidF[0][1]=S[0][1]+dSdc[0][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2.*dSdc[0][1]*HYPER[j].Fi[1][1]				  +dSdc[0][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		HYPER[j].dpidF[0][2]=S[0][2]+dSdc[0][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[0][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2.*dSdc[0][2]*HYPER[j].Fi[2][2];
		HYPER[j].dpidF[1][0]=S[1][0]+2.*dSdc[1][0]*HYPER[j].Fi[0][0]					 +dSdc[1][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[1][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		HYPER[j].dpidF[1][1]=S[1][1]+dSdc[1][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2.*dSdc[1][1]*HYPER[j].Fi[1][1]				  +dSdc[1][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		HYPER[j].dpidF[1][2]=S[1][2]+dSdc[1][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[1][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2.*dSdc[1][2]*HYPER[j].Fi[2][2];
		HYPER[j].dpidF[2][0]=S[2][0]+2.*dSdc[2][0]*HYPER[j].Fi[0][0]					 +dSdc[2][1]*(HYPER[j].Fi[1][0]+HYPER[j].Fi[0][1])+dSdc[2][2]*(HYPER[j].Fi[2][0]+HYPER[j].Fi[0][2]);
		HYPER[j].dpidF[2][1]=S[2][1]+dSdc[2][0]*(HYPER[j].Fi[0][1]+HYPER[j].Fi[1][0])+2.*dSdc[2][1]*HYPER[j].Fi[1][1]				  +dSdc[2][2]*(HYPER[j].Fi[2][1]+HYPER[j].Fi[1][2]);
		HYPER[j].dpidF[2][2]=S[2][2]+dSdc[2][0]*(HYPER[j].Fi[0][2]+HYPER[j].Fi[2][0])+dSdc[2][1]*(HYPER[j].Fi[1][2]+HYPER[j].Fi[2][1])+2.*dSdc[2][2]*HYPER[j].Fi[2][2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	in_Ci[D];
	delete[]	in_Ci;
}

void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,double nG[DIMENSION])
{

	int h_num=HYPER.size();
	int Nx=2*h_num;

	/////////////////dEŒvŽZ///////////////////
	for(int k=0;k<h_num;k++)
	{
		dE[k]=0;
		dE[k+h_num]=0;

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
			dE[k+h_num]+=-0.5*Dt*Dt/mi*((HYPER[i].pi[A_X][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].pi[A_X][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].pi[A_X][2]*HYPER1[k*h_num+i].n0ij[2])*nG[A_X]
									   +(HYPER[i].pi[A_Y][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].pi[A_Y][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].pi[A_Y][2]*HYPER1[k*h_num+i].n0ij[2])*nG[A_Y]
									   +(HYPER[i].pi[A_Z][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].pi[A_Z][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].pi[A_Z][2]*HYPER1[k*h_num+i].n0ij[2])*nG[A_Z]);
			dg[i][k+h_num]=0.5*Dt*Dt/mi*HYPER[i].J*((HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[k*h_num+i].n0ij[2])*nG[A_X]
												   +(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[k*h_num+i].n0ij[2])*nG[A_Y]
												   +(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[k*h_num+i].n0ij[2])*nG[A_Z]);
		}
		dE[k+h_num]+=0.5*Dt/mi*(HYPER[k].half_p[A_X]*nG[A_X]+HYPER[k].half_p[A_Y]*nG[A_Y]+HYPER[k].half_p[A_Z]*nG[A_Z]);
	}


	/////////////////rEŒvŽZ///////////////////
	double llE=0;
	double mlE=0;
	double lmE=0;
	double mmE=0;
	for(int l=0;l<h_num;l++)
	{
		for(int k=0;k<h_num;k++)
		{
			rE[l*Nx+k]=0;
			rE[(l+h_num)*Nx+k]=0;
			rE[l*Nx+k+h_num]=0;		
			rE[(l+h_num)*Nx+k+h_num]=0;

			mmE=0;
			for(int i=0;i<h_num;i++)
			{
				mlE=0;
				lmE=0;
				for(int m=0;m<h_num;m++)
				{
					llE=0;
					for(int s=0;s<h_num;s++)
					{
						llE+=(HYPER[i].dpidF[A_X][A_X]*HYPER1[s*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[s*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[s*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+s].DgDq_n[A_X]
						    +(HYPER[i].dpidF[A_Y][A_X]*HYPER1[s*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[s*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[s*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+s].DgDq_n[A_Y]
							+(HYPER[i].dpidF[A_Z][A_X]*HYPER1[s*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[s*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[s*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+s].DgDq_n[A_Z];
					}
					rE[l*Nx+k]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*llE*(HYPER1[m*h_num+i].n0ij[A_X]*HYPER1[k*h_num+m].DgDq_n[A_X]+HYPER1[m*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+m].DgDq_n[A_Y]+HYPER1[m*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+m].DgDq_n[A_Z]);
				
					mlE+=((HYPER[i].dpidF[A_X][A_X]*HYPER1[l*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[l*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[l*h_num+i].n0ij[A_Z])*nG[A_X]
						 +(HYPER[i].dpidF[A_Y][A_X]*HYPER1[l*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[l*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[l*h_num+i].n0ij[A_Z])*nG[A_Y]
						 +(HYPER[i].dpidF[A_Z][A_X]*HYPER1[l*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[l*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[l*h_num+i].n0ij[A_Z])*nG[A_Z])
						 *(HYPER1[m*h_num+i].n0ij[A_X]*HYPER1[k*h_num+m].DgDq_n[A_X]+HYPER1[m*h_num+i].n0ij[A_Y]*HYPER1[k*h_num+m].DgDq_n[A_Y]+HYPER1[m*h_num+i].n0ij[A_Z]*HYPER1[k*h_num+m].DgDq_n[A_Z]);

					lmE+=(HYPER[i].dpidF[A_X][A_X]*HYPER1[m*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[m*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[m*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+m].DgDq_n[A_X]
						+(HYPER[i].dpidF[A_Y][A_X]*HYPER1[m*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[m*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[m*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+m].DgDq_n[A_Y]
						+(HYPER[i].dpidF[A_Z][A_X]*HYPER1[m*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[m*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[m*h_num+i].n0ij[A_Z])*HYPER1[l*h_num+m].DgDq_n[A_Z];

				}
				rE[l*Nx+k]+=0.25*Dt*Dt/mi*(HYPER1[k*h_num+i].DgDq_n[A_X]*HYPER1[l*h_num+i].DgDq_n[A_X]+HYPER1[k*h_num+i].DgDq_n[A_Y]*HYPER1[l*h_num+i].DgDq_n[A_Y]+HYPER1[k*h_num+i].DgDq_n[A_Z]*HYPER1[l*h_num+i].DgDq_n[A_Z]);
	
				rE[(l+h_num)*Nx+k]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*mlE;
				rE[l*Nx+k+h_num]+=-0.25*Dt*Dt*Dt*Dt/mi/mi/V*lmE*(HYPER1[k*h_num+i].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+i].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].n0ij[A_Z]*nG[A_Z]);

				mmE+=((HYPER[i].dpidF[A_X][A_X]*HYPER1[l*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_X][A_Y]*HYPER1[l*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_X][A_Z]*HYPER1[l*h_num+i].n0ij[A_Z])*nG[A_X]
				  	 +(HYPER[i].dpidF[A_Y][A_X]*HYPER1[l*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Y][A_Y]*HYPER1[l*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Y][A_Z]*HYPER1[l*h_num+i].n0ij[A_Z])*nG[A_Y]
					 +(HYPER[i].dpidF[A_Z][A_X]*HYPER1[l*h_num+i].n0ij[A_X]+HYPER[i].dpidF[A_Z][A_Y]*HYPER1[l*h_num+i].n0ij[A_Y]+HYPER[i].dpidF[A_Z][A_Z]*HYPER1[l*h_num+i].n0ij[A_Z])*nG[A_Z])
					 *(HYPER1[k*h_num+i].n0ij[A_X]*nG[A_X]+HYPER1[k*h_num+i].n0ij[A_Y]*nG[A_Y]+HYPER1[k*h_num+i].n0ij[A_Z]*nG[A_Z]);

			}
			rE[(l+h_num)*Nx+k]+=-0.25*Dt*Dt/mi*(HYPER1[k*h_num+l].DgDq_n[A_X]*nG[A_X]+HYPER1[k*h_num+l].DgDq_n[A_Y]*nG[A_Y]+HYPER1[k*h_num+l].DgDq_n[A_Z]*nG[A_Z]);
			rE[l*Nx+k+h_num]+=-0.25*Dt*Dt/mi*(HYPER1[l*h_num+k].DgDq_n[A_X]*nG[A_X]+HYPER1[l*h_num+k].DgDq_n[A_Y]*nG[A_Y]+HYPER1[l*h_num+k].DgDq_n[A_Z]*nG[A_Z]);
			rE[(l+h_num)*Nx+k+h_num]+=0.25*Dt*Dt*Dt*Dt/mi/mi/V*mmE;
		}
		rE[(l+h_num)*Nx+l+h_num]+=0.25*Dt*Dt/mi*(nG[A_X]*nG[A_X]+nG[A_Y]*nG[A_Y]+nG[A_Z]*nG[A_Z]);
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











