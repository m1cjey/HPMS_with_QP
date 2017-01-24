#include "stdafx.h"	

void q_QP_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0,double **F,int Nw);
void q_variables_g(mpsconfig &CON, vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *hp_x,double *hp_y,double *hp_z,int Nw,double *g,double *dg,double *dT, double *rT,double *d_sum);
void output_data_g(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double *dg, double *th_g, double *d, int h_num,int count,int count_min,int t,double E,double En,double E0,double mi, double V,double **F);

void calc_n_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void p_QP_g(mpsconfig &CON, vector<mpselastic> &PART, vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double E0,double **F);
void p_variables_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int h_num,double Dt,double mi,double **F);
void p_nab_g(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dE, double *G, double Dt, double V, double mi);
void p_lap_g(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *rE, double *dG, double Dt, double V, double mi);
void output_data_p_g(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double *dG, double *G, double *th_G, double *rT, double *dT, double T, double *d, int Nx, int h_num,int count,int count_min,int t,double E,double En, double E0);
//void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT, double *dN,double *g, double **dg, double *h, double *th_h, int Nx, int h_num,int count,int count_min,int t,double E);


void calc_HYPER_QP_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F,int Nw)
{
	////////////íËã`///////////////
	int h_num=HYPER.size();
	int Nx=h_num*2;

	double V=get_volume(&CON);
	double mi=0.;
	double mh=CON.get_hyper_density()*V;
	double ms=CON.get_silicone_density()*V;
	double Dt=CON.get_dt();
	double nG[DIMENSION]={0,0,1.};
	double G=9.8;
	double E0_t=0.;	


	cout<<"QP_g start-------------------------"<<endl;
	//ofstream fq0("q0_QP.csv");
	//ofstream fp0("hp0_QP.csv");
	//ofstream fhp0("hp0_QP.csv");

	//ofstream fqn("qn_QP.csv");
	//ofstream fpn("pn_QP.csv");
	//ofstream fhpn("hpn_QP.csv");
	//for(int i=0;i<h_num;i++)
	//{
	//	fq0<<HYPER[i].q_n[A_X]<<","<<HYPER[i].q_n[A_Y]<<","<<HYPER[i].q_n[A_Z]<<endl;
	//	fp0<<HYPER[i].p_n[A_X]<<","<<HYPER[i].p_n[A_Y]<<","<<HYPER[i].p_n[A_Z]<<endl;
	//	fhp0<<HYPER[i].ph_n[A_X]<<","<<HYPER[i].ph_n[A_Y]<<","<<HYPER[i].ph_n[A_Z]<<endl;
	//	fqn<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;
	//	fpn<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
	//	fhpn<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
	//}
	//fq0.close();
	//fp0.close();
	//fhp0.close();
	//fqn.close();
	//fpn.close();
	//fhpn.close();

	////////////QPåvéZ///////////////		
	q_QP_g(CON,PART,HYPER,HYPER1,t,h_num,Nx,V,mi,Dt,nG,E0_t,F,Nw);

	calc_n_g(CON,PART,HYPER,HYPER1);

	p_QP_g(CON,PART,HYPER,HYPER1,t,h_num,Nx,V,mi,Dt,E0_t,F);
	for(int i=0;i<h_num;i++)	HYPER[i].lambda=HYPER[i].lam;

	for(int i=0;i<h_num;i++)
	{
		if(PART[i].toFEM==1)	mi=mh;
		else
		{
			mi=ms;
		}
		HYPER[i].E0=0.5/mi*(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z])+V*HYPER[i].W;
	}
}

void q_QP_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double nG[DIMENSION], double E0,double **F,int Nw)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1.e-5;
	double ep_min=1.e-5;
	double d_ep=1.e-20;
	double d_sum=0.;

	double p_p[DIMENSION]={0,0,0}; 
	double hp[DIMENSION]={0,0,0};
	double W=0.;
	double Efi=0.;
	double r=1.;

	double En=0.;

	double T=0.;
	double *dT=new double [h_num];	
	double *rT=new double [h_num*h_num];


	double *g=new double [h_num];
	double *dg=new double [h_num*h_num];
	double *th_g=new double [h_num];
	double *h=new double [Nw];
	double *dh=new double [h_num*h_num];
	double **rh=new double *[Nw];
	for(int i=0;i<Nw;i++)	rh[i]=new double [h_num*h_num];
	double *th_h=new double [Nw];
	double *hp_x=new double [h_num];
	double *hp_y=new double [h_num];
	double *hp_z=new double [h_num];
	double dr[DIMENSION]={0,0,0};

	double mh=V*CON.get_hyper_density();
	double ms=V*CON.get_silicone_density();
	double G=9.8;
	double p_half_p[DIMENSION]={0,0,0};
	double li=0.,lk=0;

	bool flag_FEM=CON.get_FEM_flag();
	int flag_G=CON.get_flag_G();
	double Ev=0.;
	double lam=0.;
	double Dgji_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};

	////////////èâä˙âªéZ///////////////
	for(int i=0;i<Nw;i++)
	{
		HYPER[i].lam=1.;
		HYPER[i].mu=1.;
		dT[i]=0.;
		g[i]=0.;
		h[i]=0.;
		th_g[i]=0.;
		th_h[i]=0.;
		for(int j=0;j<h_num;j++)
		{
			dg[i*h_num+j]=0.;
			rT[i*h_num+j]=0.;
			for(int k=0;k<h_num;k++)	rh[i][j*h_num+k]=0.;

		}
		hp_x[i]=0.;
		hp_y[i]=0.;
		hp_z[i]=0.;

	}
	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].lam=1.;
		HYPER[i].mu=1.;
		dT[i]=0.;
		g[i]=0.;
		th_g[i]=0.;

		for(int j=0;j<h_num;j++)
		{
			dg[i*h_num+j]=0.;
			dh[i*h_num+j]=0.;
			rT[i*h_num+j]=0.;
		}

		hp_x[i]=0.;
		hp_y[i]=0.;
		hp_z[i]=0.;

		//half_pÇÃåvéZ
		p_half_p[A_X]=0.;
		p_half_p[A_Y]=0.;
		p_half_p[A_Z]=0.;

		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{	
			int k=HYPER[i].NEI[j];
			p_half_p[A_X]+=HYPER[k].stress_n[0][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[0][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[0][2]*HYPER1[k*h_num+i].DgDq_n[2];
			p_half_p[A_Y]+=HYPER[k].stress_n[1][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[1][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[1][2]*HYPER1[k*h_num+i].DgDq_n[2];
			p_half_p[A_Z]+=HYPER[k].stress_n[2][0]*HYPER1[k*h_num+i].DgDq_n[0]+HYPER[k].stress_n[2][1]*HYPER1[k*h_num+i].DgDq_n[1]+HYPER[k].stress_n[2][2]*HYPER1[k*h_num+i].DgDq_n[2];
		}//jÇ…ä÷Ç∑ÇÈforï∂ÇÃèIÇÌÇË
		p_half_p[A_X]+=HYPER[i].stress_n[0][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[0][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[0][2]*HYPER1[i*h_num+i].DgDq_n[2];
		p_half_p[A_Y]+=HYPER[i].stress_n[1][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[1][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[1][2]*HYPER1[i*h_num+i].DgDq_n[2];
		p_half_p[A_Z]+=HYPER[i].stress_n[2][0]*HYPER1[i*h_num+i].DgDq_n[0]+HYPER[i].stress_n[2][1]*HYPER1[i*h_num+i].DgDq_n[1]+HYPER[i].stress_n[2][2]*HYPER1[i*h_num+i].DgDq_n[2];

		//èdóÕÇÃâeãø
		p_half_p[A_Z]-=G*mh;
		//îSê´çÄÇÃâeãø
		p_half_p[A_X]+=HYPER[i].vis_force[A_X];
		p_half_p[A_Y]+=HYPER[i].vis_force[A_Y];
		p_half_p[A_Z]+=HYPER[i].vis_force[A_Z];
		//é•èÍÇÃçló∂
		p_half_p[A_X]+=F[A_X][i]*V;//density;
		p_half_p[A_Y]+=F[A_Y][i]*V;//density;
		p_half_p[A_Z]+=F[A_Z][i]*V;//density;

		//âºÇÃâ^ìÆó åvéZ
		hp_x[i]=HYPER[i].p_n[A_X]+Dt*0.5*p_half_p[A_X];
		hp_y[i]=HYPER[i].p_n[A_Y]+Dt*0.5*p_half_p[A_Y];
		hp_z[i]=HYPER[i].p_n[A_Z]+Dt*0.5*p_half_p[A_Z];
	}


	while(E_min>ep_min)
	{
		count_min++;

		E=1;
		count=0;
		while(E>ep)
		{
			count++;
			////////////èâä˙âªéZ///////////////
			for(int i=0;i<Nw;i++)
			{
				for(int j=0;j<h_num;j++)
				{
					for(int k=0;k<h_num;k++)	rh[i][j*h_num+k]=0.;
				}
			}
			for(int i=0;i<h_num;i++)
			{
				for(int j=0;j<h_num;j++)
				{
					dh[j*h_num+i]=0.;
				}
			}

			q_variables_g(CON,PART,HYPER,HYPER1,hp_x,hp_y,hp_z,Nw,g,dg,dT,rT,&d_sum);

			for(int i=0;i<Nw;i++)
			{
				h[i]=(PART[i].r[A_X]-PART[i].q0[A_X])*(PART[i].r[A_X]-PART[i].q0[A_X])+
					(PART[i].r[A_Y]-PART[i].q0[A_Y])*(PART[i].r[A_Y]-PART[i].q0[A_Y])+
					(PART[i].r[A_Z]-PART[i].q0[A_Z])*(PART[i].r[A_Z]-PART[i].q0[A_Z]);
				cout<<"h"<<h[i]<<endl;
			}
			for(int i=0;i<Nw;i++)
			{
				if(h[i]>1.e-20)
				{
					int Ni=HYPER[i].N;
					for(int j=0;j<Ni;j++)
					{
						int jn=HYPER[i].NEI[j];
						//cout<<"dh="<<dh[i*h_num+jn]<<endl;
						Dgji_n[A_X]=HYPER1[jn*h_num+i].DgDq_n[A_X];
						Dgji_n[A_Y]=HYPER1[jn*h_num+i].DgDq_n[A_Y];
						Dgji_n[A_Z]=HYPER1[jn*h_num+i].DgDq_n[A_Z];
						dr[A_X]=PART[i].r[A_X]-PART[i].q0[A_X];
						dr[A_Y]=PART[i].r[A_Y]-PART[i].q0[A_Y];
						dr[A_Z]=PART[i].r[A_Z]-PART[i].q0[A_Z];

						dh[i*h_num+jn]+=-Dt*Dt/mi*(Dgji_n[A_X]*(dr[A_X])
							+Dgji_n[A_Y]*(dr[A_Y])
							+Dgji_n[A_Z]*(dr[A_Z]));
						//cout<<"dh="<<dh[i*h_num+jn]<<endl;
					}
					dh[i*h_num+i]+=-Dt*Dt/mi*(HYPER1[i*h_num+i].DgDq_n[A_X]*(PART[i].r[A_X]-PART[i].q0[A_X])
						+HYPER1[i*h_num+i].DgDq_n[A_Y]*(PART[i].r[A_Y]-PART[i].q0[A_Y])
						+HYPER1[i*h_num+i].DgDq_n[A_Z]*(PART[i].r[A_Z]-PART[i].q0[A_Z]));
					//cout<<"dh="<<dh[i*h_num+i]<<endl;
				}
			}
			for(int k=0;k<Nw;k++)
			{
				if(h[k]>1.e-20)
				{
					int Nk=HYPER[k].N;
					for(int i=0;i<Nk;i++)
					{
						int in=HYPER[k].NEI[i];
						for(int j=0;j<Nk;j++)
						{
							int jn=HYPER[k].NEI[j];
							rh[k][jn*h_num+in]+=0.5*Dt*Dt*Dt*Dt/mi/mi*(HYPER1[jn*h_num+k].DgDq_n[A_X]*HYPER1[in*h_num+k].DgDq_n[A_X]
							+HYPER1[jn*h_num+k].DgDq_n[A_Y]*HYPER1[in*h_num+k].DgDq_n[A_Y]
							+HYPER1[jn*h_num+k].DgDq_n[A_Z]*HYPER1[in*h_num+k].DgDq_n[A_Z]);
							cout<<"rh="<<rh[k][jn*h_num+in]<<endl;
						}
						rh[k][k*h_num+in]+=0.5*Dt*Dt*Dt*Dt/mi/mi*(HYPER1[k*h_num+k].DgDq_n[A_X]*HYPER1[in*h_num+k].DgDq_n[A_X]
						+HYPER1[k*h_num+k].DgDq_n[A_Y]*HYPER1[in*h_num+k].DgDq_n[A_Y]
						+HYPER1[k*h_num+k].DgDq_n[A_Z]*HYPER1[in*h_num+k].DgDq_n[A_Z]);
						rh[k][in*h_num+k]+=0.5*Dt*Dt*Dt*Dt/mi/mi*(HYPER1[in*h_num+k].DgDq_n[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]
						+HYPER1[in*h_num+k].DgDq_n[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]
						+HYPER1[in*h_num+k].DgDq_n[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);
						cout<<"rh="<<rh[k][k*h_num+in]<<endl;
						cout<<"rh="<<rh[k][in*h_num+k]<<endl;
					}
					rh[k][k*h_num+k]+=0.5*Dt*Dt*Dt*Dt/mi/mi*(HYPER1[k*h_num+k].DgDq_n[A_X]*HYPER1[k*h_num+k].DgDq_n[A_X]
					+HYPER1[k*h_num+k].DgDq_n[A_Y]*HYPER1[k*h_num+k].DgDq_n[A_Y]
					+HYPER1[k*h_num+k].DgDq_n[A_Z]*HYPER1[k*h_num+k].DgDq_n[A_Z]);
					cout<<"rh="<<rh[k][k*h_num+k]<<endl;
				}
			}

			for(int i=0;i<Nw;i++)
			{
				if(h[i]>1.e-20)
				{
					for(int j=0;j<h_num;j++)
					{
						dT[j]+=2*r*(h[i]+th_h[i])*dh[i*h_num+j];
						for(int k=0;k<h_num;k++)	rT[j*h_num+k]+=2*r*dh[i*h_num+k]*dh[i*h_num+j]+2*r*(h[i]+th_h[i])*rh[i][j*h_num+k];
					}
				}
			}


			if(d_sum<1.e-20)	break;
			gauss(rT,dT,h_num);

			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].lam-=dT[i];
				E_sum+=fabs(dT[i]);
			}
			E=E_sum;
			/////////////pn1_2åvéZ
			for(int i=Nw;i<h_num;i++)
			{
				//half_pÇÃåvéZ
				p_half_p[A_X]=0.;
				p_half_p[A_Y]=0.;
				p_half_p[A_Z]=0.;

				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{	
					int k=HYPER[i].NEI[j];
					lk=HYPER[k].lambda;
					p_half_p[A_X]+=-lk*HYPER1[k*h_num+i].DgDq_n[0];
					p_half_p[A_Y]+=-lk*HYPER1[k*h_num+i].DgDq_n[1];
					p_half_p[A_Z]+=-lk*HYPER1[k*h_num+i].DgDq_n[2];
				}//jÇ…ä÷Ç∑ÇÈforï∂ÇÃèIÇÌÇË
				li=HYPER[i].lambda;
				p_half_p[A_X]+=-li*HYPER1[i*h_num+i].DgDq_n[0];
				p_half_p[A_Y]+=-li*HYPER1[i*h_num+i].DgDq_n[1];
				p_half_p[A_Z]+=-li*HYPER1[i*h_num+i].DgDq_n[2];

				//à íuç¿ïWÇÃåvéZ
				PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt*(hp_x[i]+Dt*0.5*p_half_p[A_X])/mh;
				PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt*(hp_y[i]+Dt*0.5*p_half_p[A_Y])/mh;
				PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt*(hp_z[i]+Dt*0.5*p_half_p[A_Z])/mh;
			}


			//cout<<"	E"<<count<<"="<<E<<endl;
			//output_data(PART,HYPER,HYPER1,T,dT,rT,g,dg,th_g,d,h_num,count,count_min,t,E,En,E0,mi,V);
			//if(count==1||count%100==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				//output_data_g(PART,HYPER,HYPER1,T,dT,rT,g,dg,th_g,d,h_num,count,count_min,t,E,En,E0,mi,V,F);
			}
		}
		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=(HYPER[i].old_lam-HYPER[i].lam)*(HYPER[i].old_lam-HYPER[i].lam);
		E_min=sqrt(E_sum);
		if(E_min<ep_min*1000)	r*=4;

		for(int i=0;i<Nw;i++)
		{
			th_h[i]+=h[i];
		}

		//for(int i=0;i<h_num;i++)
		//{
		//	//th_g[i]+=g[i]*g[i];
		//	if(g[i]<0)	th_g[i]+=-1.*g[i];
		//	else
		//	{
		//		th_g[i]+=g[i];
		//	}
		//	//if(h[i]+th_h[i]>0)	th_h[i]+=h[i];
		//}

		//////ÉfÅ[É^èoóÕ
		cout<<"Emin"<<count_min<<"="<<E_min<<endl;
		//if(count_min==1||count_min%100==0)
		//{
		//	stringstream ss0;
		//	ss0<<"./E_min/E"<<t<<".csv"<<endl;
		//	stringstream ss1;
		//	ss1<<"./p/hp"<<t<<"_"<<count_min<<".csv"<<endl;
		//	stringstream ss2;
		//	ss2<<"./q/q"<<t<<"_"<<count_min<<".csv"<<endl;
		//	stringstream ss3;
		//	ss3<<"./lam/lam"<<t<<"_"<<count_min<<".csv"<<endl;
		//	stringstream ss4;
		//	ss4<<"./T/T"<<t<<".csv"<<endl;

		//	if(count_min==1)
		//	{
		//		ofstream fem(ss0.str(), ios::trunc);
		//		fem<<E_min<<endl;
		//		fem.close();
		//		ofstream ft(ss4.str(), ios::trunc);
		//		ft<<T<<","<<En<<","<<E0<<endl;
		//		ft.close();
		//	}
		//	else
		//	{
		//		ofstream fem(ss0.str(), ios::app);
		//		fem<<E_min<<endl;
		//		fem.close();
		//		ofstream ft(ss4.str(), ios::app);
		//		ft<<T<<","<<En<<endl;
		//		ft.close();
		//	}

		//	ofstream fhp(ss1.str());
		//	ofstream fq(ss2.str());
		//	ofstream fl(ss3.str());
		//	for(int i=0;i<h_num;i++)
		//	{
		//		fhp<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<endl;
		//		fq<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<endl;			
		//		fl<<HYPER[i].lam<<endl;
		//	}
		//	fhp.close();
		//	fq.close();
		//	fl.close();

		//}





		if(count_min>c_max)	break;
	}
	cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;



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
	for(int i=0;i<Nw;i++)	delete[]	rh[i];
	delete[]	rh;
	delete[]	g;
	delete[]	dg;
	delete[]	h;
	delete[]	dh;
	delete[]	dT;
	delete[]	rT;
	delete[]	th_g;
	delete[]	th_h;
	delete[]	hp_x;
	delete[]	hp_y;
	delete[]	hp_z;

}

void q_variables_g(mpsconfig &CON, vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *hp_x,double *hp_y,double *hp_z,int Nw,double *g,double *dg,double *dT, double *rT,double *d_sum)
{
	int h_num=HYPER.size();

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
	double mh=V*CON.get_hyper_density();

	/////////////pn1_2åvéZ
	double li=0.,lk=0;
	double p_half_p[DIMENSION]={0,0,0};
	for(int i=0;i<h_num;i++)
	{
		//half_pÇÃåvéZ
		p_half_p[A_X]=0.;
		p_half_p[A_Y]=0.;
		p_half_p[A_Z]=0.;

		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{	
			int k=HYPER[i].NEI[j];
			lk=HYPER[k].lambda;
			p_half_p[A_X]+=-lk*HYPER1[k*h_num+i].DgDq_n[0];
			p_half_p[A_Y]+=-lk*HYPER1[k*h_num+i].DgDq_n[1];
			p_half_p[A_Z]+=-lk*HYPER1[k*h_num+i].DgDq_n[2];
		}//jÇ…ä÷Ç∑ÇÈforï∂ÇÃèIÇÌÇË
		li=HYPER[i].lambda;
		p_half_p[A_X]+=-li*HYPER1[i*h_num+i].DgDq_n[0];
		p_half_p[A_Y]+=-li*HYPER1[i*h_num+i].DgDq_n[1];
		p_half_p[A_Z]+=-li*HYPER1[i*h_num+i].DgDq_n[2];

		//à íuç¿ïWÇÃåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt*(hp_x[i]+Dt*0.5*p_half_p[A_X])/mh;
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt*(hp_y[i]+Dt*0.5*p_half_p[A_Y])/mh;
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt*(hp_z[i]+Dt*0.5*p_half_p[A_Z])/mh;
	}



	/////////////F, J, t_inverseÇÃåvéZ
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	double J=0.;
	double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double a[DIMENSION]={0,0,0};
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].old_lam=HYPER[i].lam;

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
	}	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;


	for(int k=0;k<h_num;k++)
	{
		g[k]=V*(1-HYPER[k].J);
		dT[k]=0;

		//dE[k]=0.;
		for(int i=0;i<h_num;i++)
		{
			dg[i*h_num+k]=0.;
			//rE[i*h_num+k]=0.;
		}
	}

	/////////////////dE, dgåvéZ///////////////////
	//////////íËã`
	double Dgkk_n[DIMENSION]={0,0,0};
	double Dgik_n[DIMENSION]={0,0,0};

	double Dgik[DIMENSION]={0,0,0};
	double Dgkk[DIMENSION]={0,0,0};
	double Ji=0.,Jk=0.;

	double d_sum2=0.;

	for(int k=0;k<h_num;k++)
	{
		Dgkk_n[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];	Dgkk_n[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];	Dgkk_n[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];

		Jk=HYPER[k].J;
		Dgkk[A_X]=Jk*(HYPER[k].t_inverse_Fi[0][0]*HYPER1[k*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[0][1]*HYPER1[k*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[0][2]*HYPER1[k*h_num+k].n0ij[2]);
		Dgkk[A_Y]=Jk*(HYPER[k].t_inverse_Fi[1][0]*HYPER1[k*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[1][1]*HYPER1[k*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[1][2]*HYPER1[k*h_num+k].n0ij[2]);
		Dgkk[A_Z]=Jk*(HYPER[k].t_inverse_Fi[2][0]*HYPER1[k*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[2][1]*HYPER1[k*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[2][2]*HYPER1[k*h_num+k].n0ij[2]);

		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];

			Ji=HYPER[i].J;
			Dgik[A_X]=Ji*(HYPER[i].t_inverse_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[0][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[0][2]*HYPER1[k*h_num+i].n0ij[2]);
			Dgik[A_Y]=Ji*(HYPER[i].t_inverse_Fi[1][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[1][2]*HYPER1[k*h_num+i].n0ij[2]);
			Dgik[A_Z]=Ji*(HYPER[i].t_inverse_Fi[2][0]*HYPER1[k*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[2][1]*HYPER1[k*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);

			Dgik_n[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];	Dgik_n[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];	Dgik_n[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];

				/////////////dlam g
				dg[in*h_num+jn]-=0.5*Dt*Dt/mi*(Dgik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+Dgik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+Dgik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
				dT[jn]+=dg[in*h_num+jn]*g[in];
				d_sum2+=sqrt(dT[jn]);
			}
			/////////////dlam g
			dg[in*h_num+k]-=0.5*Dt*Dt/mi*(Dgik[A_X]*Dgkk_n[A_X]+Dgik[A_Y]*Dgkk_n[A_Y]+Dgik[A_Z]*Dgkk_n[A_Z]);
			dg[k*h_num+in]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgik_n[A_X]+Dgkk[A_Y]*Dgik_n[A_Y]+Dgkk[A_Z]*Dgik_n[A_Z]);
			dT[k]+=dg[in*h_num+k]*g[in];
			dT[in]+=dg[k*h_num+in]*g[k];
			d_sum2+=sqrt(dT[k])+sqrt(dT[in]);

		}
		/////////////dlam g
		dg[k*h_num+k]-=0.5*Dt*Dt/mi*(Dgkk[A_X]*Dgkk_n[A_X]+Dgkk[A_Y]*Dgkk_n[A_Y]+Dgkk[A_Z]*Dgkk_n[A_Z]);
		dT[k]+=dg[k*h_num+k]*g[k];
		d_sum2+=sqrt(dT[k]);
	}

	*d_sum=d_sum2;

	for(int k=0;k<h_num;k++)
	{	
		for(int l=0;l<h_num;l++)
		{	
			//dT[k]+=r*2.*dg[i*h_num+k]*g[i]*(g[i]*g[i]+th_g[i]);
			//if(g[i]<0)	dT[k]+=-1.*r*dg[i*h_num+k]*(-1.*g[i]+th_g[i]);
			//else
			//{
			//	dT[k]+=r*dg[i*h_num+k]*(g[i]+th_g[i]);
			//}
			int Nl=HYPER[l].N;
			rT[l*h_num+k]=0.;
			for(int i=0;i<Nl;i++)
			{	
				int in=HYPER[l].NEI[i];
				rT[l*h_num+k]+=dg[in*h_num+k]*dg[in*h_num+l];
				//rT[l*h_num+k]+=r*dg[i*h_num+k]*dg[i*h_num+l];
				//rT[l*h_num+k]+=r*2.*dg[i*h_num+k]*dg[i*h_num+l]*(3*g[i]*g[i]+th_g[i]);
			}
		}
	}
			for(int i=0;i<h_num;i++)
			{
				cout<<"dT"<<dT[i]<<endl;
				for(int j=0;j<h_num;j++)	cout<<"rT"<<rT[i]<<endl;
			}


}



void output_data_g(vector<mpselastic>PART,vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double T,double *dT, double *rT, double *g, double *dg, double *th_g, double *d, int h_num,int count,int count_min,int t,double E,double En,double E0,double mi, double V,double **F)
{

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
	stringstream ss13;
	ss13<<"./g/J_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss14;
	ss14<<"./g/DgDq_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss15;
	ss15<<"./g/stress_"<<t<<"_"<<count_min<<"_"<<count<<".csv";
	stringstream ss16;
	ss16<<"./T/En_"<<t<<"_"<<count_min<<"_"<<count<<".csv";


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
	ofstream f_en(ss16.str());

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

	for(int i=0;i<h_num;i++)
	{
		f_dt<<dT[i]<<endl;
		f_lam<<HYPER[i].lam<<endl;

		for(int k=0;k<h_num;k++)
		{
			f_rt<<rT[i*h_num+k]<<",";
			f_dg<<dg[i*h_num+k]<<",";
		}
		f_rt<<endl;
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
		f_en<<0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])+V*HYPER[i].W+HYPER[i].lam*(1.-HYPER[i].J)<<","<<0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])<<","<<V*HYPER[i].W<<","<<HYPER[i].lam*(1.-HYPER[i].J)<<endl;
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
	f_en.close();

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


void calc_n_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
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


void p_QP_g(mpsconfig &CON, vector<mpselastic> &PART, vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, int h_num, int Nx, double V, double mi, double Dt, double E0,double **F)
{
	////////////íËã`///////////////
	int count=0;
	int count_min=0;
	int c_max=1000;

	double E=1;
	double E_min=1;
	double E_sum=0;
	double ep=1.e-5;
	double ep_min=1.e-5;
	double d_ep=1.e-20;

	double p_p[DIMENSION]={0,0,0};
	double p[DIMENSION]={0,0,0};
	double W=0.;
	double r=1.e100;	//r=0.01ÇÕåÎç∑Ç™ëÂÇ´Ç¢
	//if(t>101) r=10.;
	double d_sum=0.;

	double En=0.;
	double *dE=new double [h_num];	
	double *rE=new double [h_num*h_num];

	double T=0.;
	double *dT=new double [h_num];	
	double *rT=new double [h_num*h_num];

	double *d=new double [h_num];
	double *B=new double [h_num*h_num];
	double g=9.8;
	double *G=new double [h_num];
	double *dG=new double [h_num*h_num];
	double *th_G=new double [h_num];
	double Efi=0.;
	double Ev=0.;
	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].lam=1.;
		HYPER[i].mu=1.;
		dE[i]=0.;
		dT[i]=0.;
		G[i]=0.;
		th_G[i]=0.;
		d[i]=0.;
		for(int j=0;j<h_num;j++)
		{
			dG[i*h_num+j]=0.;
			rT[i*h_num+j]=0.;
			rE[i*h_num+j]=0.;
			B[i*h_num+j]=0.;
		}
	}

	p_lap_g(HYPER,HYPER1,rE,dG,Dt,V,mi);



	while(E_min>ep_min)
	{
		count_min++;

		for(int i=0;i<h_num;i++)
		{
			HYPER[i].old_mu=HYPER[i].mu;
		}

		E=1;
		count=0;
		while(E>ep)
		{
			count++;
			p_variables_g(CON, PART, HYPER,HYPER1,h_num,Dt,mi,F);
			p_nab_g(HYPER,HYPER1,dE,G,Dt,V,mi);

			//En=0;
			//for(int i=0;i<h_num;i++)	
			//{

			//	p[A_X]=HYPER[i].p[A_X];	p[A_Y]=HYPER[i].p[A_Y];	p[A_Z]=HYPER[i].p[A_Z];
			//	W=HYPER[i].W;
			//	Efi=HYPER[i].Ef;
			//	Ev=HYPER[i].vis_force[A_X]*(PART[i].r[A_X]-HYPER[i].q_n[A_X])+HYPER[i].vis_force[A_Y]*(PART[i].r[A_Y]-HYPER[i].q_n[A_Y])+HYPER[i].vis_force[A_Z]*(PART[i].r[A_Z]-HYPER[i].q_n[A_Z]);				//En+=0.5/mi*(p[A_X]*p[A_X]+p[A_Y]*p[A_Y]+p[A_Z]*p[A_Z])+V*W+HYPER[i].mu*g;
			//	En+=0.5/mi*(p[A_X]*p[A_X]+p[A_Y]*p[A_Y]+p[A_Z]*p[A_Z])+V*W-mi*Ev+mi*g*PART[i].r[A_Z]+Efi;
			//}


			//T=(En-E0)*(En-E0);
			//for(int k=0;k<h_num;k++)
			//{	
			//	dT[k]=2.*dE[k]*(E-E0);
			//	for(int l=0;l<h_num;l++)	rT[l*h_num+k]=2.*rE[l*h_num+k]*(E-E0)+2*dE[k]*dE[l];
			//}
			T=0.;
			for(int k=0;k<h_num;k++)
			{	
				if(G[k]<0)	T+=0.5*r*(-1.*G[k]+th_G[k])*(-1.*G[k]+th_G[k]);
				else
				{
					T+=0.5*r*(G[k]+th_G[k])*(G[k]+th_G[k]);
				}

				dT[k]=0.;
				for(int i=0;i<h_num;i++)
				{	

					if(G[i]<0)	dT[k]+=-1.*r*dG[i*h_num+k]*(-1.*G[i]+th_G[i]);
					else
					{
						dT[k]+=r*dG[i*h_num+k]*(G[i]+th_G[i]);
					}

					rT[i*h_num+k]=0.;
					for(int s=0;s<h_num;s++)
					{	
						//rT[i*h_num+k]+=r*2.*dG[s*h_num+k]*dG[s*h_num+i]*(3*G[s]*G[s]+th_G[s]);
						rT[i*h_num+k]+=r*dG[s*h_num+k]*dG[s*h_num+i];
						//rT[i*h_num+k]+=r*dG[s*h_num+k]*dG[s*h_num+i];
					}
				}
			}

			//for(int i=0;i<h_num;i++)
			//{	
			//	G[i]/=V;
			//	//T+=0.5*r*(G[i]*G[i]+th_G[i])*(G[i]*G[i]+th_G[i]);
			//	if(G[i]<0.)T+=0.5*r*(-1.*G[i]+th_G[i])*(-1.*G[i]+th_G[i]);
			//	else
			//	{
			//		T+=0.5*r*(G[i]+th_G[i])*(G[i]+th_G[i]);
			//	}
			//}

			//for(int k=0;k<h_num;k++)
			//{	
			//	for(int i=0;i<h_num;i++)
			//	{	
			//		if(G[i]<0.)	dT[k]+=-1.*r/V*dG[i*h_num+k]*(-1.*G[i]+th_G[i]);
			//		else
			//		{
			//			dT[k]+=r/V*dG[i*h_num+k]*(G[i]+th_G[i]);
			//		}
			//		//if(G[i]<0.)	dT[k]+=-1.*r*dG[i*h_num+k]*(-1.*G[i]+th_G[i]);
			//		//else
			//		//{
			//		//	dT[k]+=r*dG[i*h_num+k]*(G[i]+th_G[i]);
			//		//}
			//		//dT[k]+=r*2.*dG[i*h_num+k]*G[i]*(G[i]*G[i]+th_G[i]);
			//	}
			//	for(int l=0;l<h_num;l++)
			//	{
			//		for(int i=0;i<h_num;i++)
			//		{	
			//			//rT[l*h_num+k]+=r*2.*dG[i*h_num+k]*dG[i*h_num+l]*(3*G[i]*G[i]+th_G[i]);
			//			rT[l*h_num+k]+=r/V/V*dG[i*h_num+k]*dG[i*h_num+l];
			//			//rT[l*h_num+k]+=r*dG[i*h_num+k]*dG[i*h_num+l];
			//		}
			//	}
			//}

			d_sum=0.;

			for(int i=0;i<h_num;i++)
			{
				d[i]=dT[i];
				d_sum+=fabs(dT[i]);
				for(int j=0;j<h_num;j++)
				{
					B[i*h_num+j]=rT[i*h_num+j];
				}
			}
			if(d_sum<1.e-20)	break;
			gauss(B,d,h_num);
			E_sum=0;
			for(int i=0;i<h_num;i++)
			{
				HYPER[i].mu-=d[i];
				E_sum+=fabs(d[i]);
			}
			E=E_sum;
			/*if(count==1||count%500==0)cout<<"E"<<count<<"="<<E<<endl;*/


			/////////////påvéZ
			p_variables_g(CON,PART,HYPER,HYPER1,h_num,Dt,mi,F);
			//if(count==1||count%100==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				output_data_p_g(HYPER,HYPER1,dG,G,th_G,rT,dT,T,d,Nx,h_num,count,count_min,t,E,En,E0);
			}

		}
		E_sum=0;
		for(int i=0;i<h_num;i++)	E_sum+=(HYPER[i].old_mu-HYPER[i].mu)*(HYPER[i].old_mu-HYPER[i].mu);
		E_min=sqrt(E_sum);
		if(E_min<ep_min*1000)	r*=4;


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

		//////ÉfÅ[É^èoóÕ
		//if(count_min==1||count_min%100==0)
		//{
		//	cout<<"Emin"<<count_min<<"="<<E_min<<endl;
		//	stringstream ss0;
		//	ss0<<"./E_min/g_E"<<t<<".csv";
		//	stringstream ss1;
		//	ss1<<"./p/g_p"<<t<<"_"<<count_min<<".csv";
		//	stringstream ss3;
		//	ss3<<"./mu/g_mu"<<t<<"_"<<count_min<<".csv";
		//	stringstream ss4;
		//	ss4<<"./p_T/g_T"<<t<<".csv";

		//	if(count_min==1)
		//	{
		//		ofstream fem(ss0.str(), ios::trunc);
		//		fem<<E_min<<endl;
		//		fem.close();
		//		ofstream ft(ss4.str(), ios::trunc);
		//		ft<<T<<","<<En<<","<<E0<<endl;
		//		ft.close();
		//	}
		//	else
		//	{
		//		ofstream fem(ss0.str(), ios::app);
		//		fem<<E_min<<endl;
		//		fem.close();
		//		ofstream ft(ss4.str(), ios::app);
		//		ft<<T<<","<<En<<endl;
		//		ft.close();
		//	}

		//	ofstream fp(ss1.str());
		//	ofstream fl(ss3.str());
		//	for(int i=0;i<h_num;i++)
		//	{
		//		fp<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		//		fl<<HYPER[i].lam<<endl;
		//	}
		//	fp.close();
		//	fl.close();
		//}


		if(count_min>c_max)	break;
	}
	cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;


	delete[]	dE;
	delete[]	rE;
	delete[]	G;
	delete[]	dG;
	delete[]	dT;
	delete[]	rT;
	delete[]	th_G;
	delete[]	d;
	delete[]	B;
}

void p_nab_g(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *dE, double *G, double Dt, double V, double mi)
{

	int h_num=HYPER.size();
	/////////////////dEåvéZ///////////////////
	for(int k=0;k<h_num;k++)
	{
		G[k]=0.;
		//dE[k]=0.;
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

			//dE[k]-=0.5*Dt/mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
			G[k]+=1./mi*(Dgki[A_X]*pi[A_X]+Dgki[A_Y]*pi[A_Y]+Dgki[A_Z]*pi[A_Z]);
		}
		//dE[k]-=0.5*Dt/mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
		//dE[k]+=V*(1.-HYPER[k].J);
		G[k]+=1./mi*(Dgkk[A_X]*pk[A_X]+Dgkk[A_Y]*pk[A_Y]+Dgkk[A_Z]*pk[A_Z]);
	}

}


void p_lap_g(vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,double *rE, double *dG, double Dt, double V, double mi)
{

	int h_num=HYPER.size();
	/////////////////dEåvéZ///////////////////
	for(int k=0;k<h_num;k++)
	{
		for(int i=0;i<h_num;i++)
		{
			dG[i*h_num+k]=0.;
			//rE[i*h_num+k]=0.;
		}
	}

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
				Dglj[A_X]=HYPER1[ln*h_num+j].DgDq[A_X];	Dglj[A_Y]=HYPER1[ln*h_num+j].DgDq[A_Y];	Dglj[A_Z]=HYPER1[ln*h_num+j].DgDq[A_Z];

				//rE[kn*h_num+ln]+=0.25*Dt*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
				dG[kn*h_num+ln]-=0.5*Dt/mi*(Dgkj[A_X]*Dglj[A_X]+Dgkj[A_Y]*Dglj[A_Y]+Dgkj[A_Z]*Dglj[A_Z]);
			}
			//rE[kn*h_num+j]+=0.25*Dt*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			//rE[j*h_num+kn]+=0.25*Dt*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

			dG[kn*h_num+j]-=0.5*Dt/mi*(Dgkj[A_X]*Dgjj[A_X]+Dgkj[A_Y]*Dgjj[A_Y]+Dgkj[A_Z]*Dgjj[A_Z]);
			dG[j*h_num+kn]-=0.5*Dt/mi*(Dgjj[A_X]*Dgkj[A_X]+Dgjj[A_Y]*Dgkj[A_Y]+Dgjj[A_Z]*Dgkj[A_Z]);

		}
		//rE[j*h_num+j]+=0.25*Dt*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);

		dG[j*h_num+j]-=0.5*Dt/mi*(Dgjj[A_X]*Dgjj[A_X]+Dgjj[A_Y]*Dgjj[A_Y]+Dgjj[A_Z]*Dgjj[A_Z]);
	}

}


void p_variables_g(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int h_num,double Dt,double mi,double **F)
{
	/////////////påvéZ
	double p_p[DIMENSION]={0,0,0};
	double DgDq_ji[DIMENSION]={0,0,0};
	double DgDq_ii[DIMENSION]={0,0,0};
	double mu=0.;
	bool flag_FEM=CON.get_FEM_flag();
	double G=9.8;

	double V=get_volume(&CON);
	double mh=V*CON.get_hyper_density();
	double ms=V*CON.get_silicone_density();

	for(int i=0;i<h_num;i++)
	{
		if(HYPER[i].fw==1)
		{
			HYPER[i].p[A_X]=0.;
			HYPER[i].p[A_Y]=0.;
			HYPER[i].p[A_Z]=0.;

		}
		else
		{
			if(PART[i].toFEM==ON)	mi=mh;
			else
			{
				mi=ms;

			}
			p_p[A_X]=0.;	p_p[A_Y]=0.;	p_p[A_Z]=0.;
			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{
				int jn=HYPER[i].NEI[j];
				mu=HYPER[jn].mu;
				DgDq_ji[A_X]=HYPER1[jn*h_num+i].DgDq[A_X];	DgDq_ji[A_Y]=HYPER1[jn*h_num+i].DgDq[A_Y];	DgDq_ji[A_Z]=HYPER1[jn*h_num+i].DgDq[A_Z];
				p_p[A_X]+=(HYPER[jn].stress[A_X][A_X]-mu)*DgDq_ji[A_X]+HYPER[jn].stress[A_X][A_Y]*DgDq_ji[A_Y]+HYPER[jn].stress[A_X][A_Z]*DgDq_ji[A_Z];
				p_p[A_Y]+=HYPER[jn].stress[A_Y][A_X]*DgDq_ji[A_X]+(HYPER[jn].stress[A_Y][A_Y]-mu)*DgDq_ji[A_Y]+HYPER[jn].stress[A_Y][A_Z]*DgDq_ji[A_Z];
				p_p[A_Z]+=HYPER[jn].stress[A_Z][A_X]*DgDq_ji[A_X]+HYPER[jn].stress[A_Z][A_Y]*DgDq_ji[A_Y]+(HYPER[jn].stress[A_Z][A_Z]-mu)*DgDq_ji[A_Z];
			}
			DgDq_ii[A_X]=HYPER1[i*h_num+i].DgDq[A_X];	DgDq_ii[A_Y]=HYPER1[i*h_num+i].DgDq[A_Y];	DgDq_ii[A_Z]=HYPER1[i*h_num+i].DgDq[A_Z];
			mu=HYPER[i].mu;
			p_p[A_X]+=(HYPER[i].stress[A_X][A_X]-mu)*DgDq_ii[A_X]+HYPER[i].stress[A_X][A_Y]	  *DgDq_ii[A_Y]+HYPER[i].stress[A_X][A_Z]			    *DgDq_ii[A_Z];
			p_p[A_Y]+=HYPER[i].stress[A_Y][A_X]		 *DgDq_ii[A_X]+(HYPER[i].stress[A_Y][A_Y]-mu)*DgDq_ii[A_Y]+HYPER[i].stress[A_Y][A_Z]			    *DgDq_ii[A_Z];
			p_p[A_Z]+=HYPER[i].stress[A_Z][A_X]		 *DgDq_ii[A_X]+HYPER[i].stress[A_Z][A_Y]	  *DgDq_ii[A_Y]+(HYPER[i].stress[A_Z][A_Z]-HYPER[i].mu)*DgDq_ii[A_Z];
			//èdóÕçÄ
			/*if(flag_G==ON)*/	p_p[A_Z]-=G*mi;
			//îSê´çÄ
			//if(flag_vis==ON)
			{
				p_p[A_X]+=mi*HYPER[i].vis_force[A_X];
				p_p[A_Y]+=mi*HYPER[i].vis_force[A_Y];
				p_p[A_Z]+=mi*HYPER[i].vis_force[A_Z];
			}
			//é•óÕçÄ
			if(flag_FEM==ON&&PART[i].toFEM==ON)
			{
				p_p[A_X]+=V*F[A_X][i];//density;
				p_p[A_Y]+=V*F[A_Y][i];//density;
				p_p[A_Z]+=V*F[A_Z][i];//density;
			}


			//â^ìÆó ÇÃçXêV
			HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*p_p[A_X];
			HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*p_p[A_Y];
			HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*p_p[A_Z];////

			//ë¨ìxÇÃçXêV
			PART[i].u[A_X]=HYPER[i].half_p[A_X]/mi;
			PART[i].u[A_Y]=HYPER[i].half_p[A_Y]/mi;
			PART[i].u[A_Z]=HYPER[i].half_p[A_Z]/mi;
			//äpâ^ìÆó ÇÃçXêV
			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];

		}
	}
}



void output_data_p_g(vector<hyperelastic>HYPER,vector<hyperelastic2>HYPER1, double *dG, double *G, double *th_G, double *rT, double *dT, double T, double *d, int Nx, int h_num,int count,int count_min,int t,double E,double En, double E0)
{
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



	for(int i=0;i<h_num;i++)
	{				
		for(int j=0;j<h_num;j++)
		{
			f_dG<<dG[i*h_num+j]<<",";
			f_rt<<rT[i*h_num+j]<<",";
		}
		f_dG<<endl;
		f_rt<<endl;
		f_G<<G[i]<<endl;
		f_dt<<dT[i]<<endl;
		f_p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
		f_d<<d[i]<<endl;
	}
	f_dG.close();
	f_G.close();
	f_rt.close();
	f_dt.close();
	f_p.close();
	f_d.close();
}









