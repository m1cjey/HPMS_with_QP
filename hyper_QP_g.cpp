#include "stdafx.h"	

void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num,double V, double mi, double Dt, double nG[DIMENSION], double E0);
void p_QP(vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int h_num,double V, double mi, double Dt, double nG[DIMENSION], double E0);
void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *rE,double **dg,double Dt, double V, double mi,double nG[DIMENSION]);
void q_variables(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);

void output_data(vector<mpselastic>PART,vector<hyperelastic>HYPER, double *rT, double *dT,double *rL,double *dL,int h_num,int count,int count_min,int t,double E);

//Hyperelstic2Ç©ÇÁÇÃà¯åpÇ¨
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

	////////////íËã`///////////////
	int h_num=HYPER.size();

	double V=get_volume(&CON);
	double mi=CON.get_hyper_density()*V;
	double Dt=CON.get_dt();
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

	////////////QPåvéZ///////////////		
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
	////////////íËã`///////////////
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

	double *dE=new double [h_num*2];
	double *dT=new double [h_num*2];	


	double *h=new double [h_num];

	double *th_h=new double [h_num];
		
	////////////èâä˙âªéZ///////////////
	for(int i=0;i<h_num;i++)
	{
		dE[i]=0;	dE[i+h_num]=0;
		dT[i]=0;	dT[i+h_num]=0;
		h[i]=0;
		th_h[i]=0;
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
			q_nab_lap(PART,HYPER,HYPER1,dE,Dt,V,mi,nG);

			En=0;
			for(int i=0;i<h_num;i++)	En+=0.5/mi*(HYPER[i].half_p[A_X]*HYPER[i].half_p[A_X]+HYPER[i].half_p[A_Y]*HYPER[i].half_p[A_Y]+HYPER[i].half_p[A_Z]*HYPER[i].half_p[A_Z])+V*HYPER[i].W;

			T=En-E0;
			//cout<<"T="<<T<<", En"<<En<<endl;
			for(int k=0;k<h_num;k++)
			{	
				dT[k]=dE[k];
			}
			//for(int k=0;k<h_num;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

			for(int i=0;i<h_num;i++)
			{	
				h[i]=-1.*(PART[i].r[A_X]*nG[A_X]+PART[i].r[A_Y]*nG[A_Y]+PART[i].r[A_Z]*nG[A_Z]);
				
				()T+=0.5*r*(h[i]*g[i]+th_g[i])*(g[i]*g[i]+th_g[i]);
			}

			for(int k=0;k<h_num;k++)
			{	
				for(int i=0;i<h_num;i++)
				{	
					dT[k]+=2.*r*dg[i][k]*g[i]*(g[i]*g[i]+th_g[i]);
					//cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<2*r*dg[i][k]*(g[i]+th_g[i])<<endl;
				}
			}
			//for(int k=0;k<h_num;k++)	cout<<"dT"<<count_min<<"_"<<count<<"_"<<k<<"="<<dT[k]<<endl;

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

			/////////////pn1_2åvéZ
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
					}//jÇ…ä÷Ç∑ÇÈforï∂ÇÃèIÇÌÇË	
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
	delete[]	th_g;

	/////////////pn1_2åvéZ
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
		}//jÇ…ä÷Ç∑ÇÈforï∂ÇÃèIÇÌÇË	
		p_half_p[A_X]+=HYPER[i].pi_n[0][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[0][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[0][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[0];
		p_half_p[A_Y]+=HYPER[i].pi_n[1][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[1][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[1][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[1];
		p_half_p[A_Z]+=HYPER[i].pi_n[2][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[2][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[2][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[2];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_half_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_half_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_half_p[A_Z];
	
		/////////////qåvéZ
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


	/////////////pn1_2åvéZ
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
		}//jÇ…ä÷Ç∑ÇÈforï∂ÇÃèIÇÌÇË	
		p_half_p[A_X]+=HYPER[i].pi_n[0][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[0][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[0][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[0];
		p_half_p[A_Y]+=HYPER[i].pi_n[1][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[1][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[1][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[1];
		p_half_p[A_Z]+=HYPER[i].pi_n[2][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].pi_n[2][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].pi_n[2][2]*HYPER1[i*h_num+i].n0ij[2]-HYPER[i].lam*HYPER1[i*h_num+i].DgDq_n[2];
		HYPER[i].half_p[A_X]=HYPER[i].p_n[A_X]+0.5*Dt*p_half_p[A_X];
		HYPER[i].half_p[A_Y]=HYPER[i].p_n[A_Y]+0.5*Dt*p_half_p[A_Y];
		HYPER[i].half_p[A_Z]=HYPER[i].p_n[A_Z]+0.5*Dt*p_half_p[A_Z];
	
		/////////////qåvéZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}


	/////////////F, J, t_inverseÇÃåvéZ
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


		for(int k=0;k<Ni;k++)
		{			
			int kn=HYPER[i].NEI[k];
			HYPER1[i*h_num+kn].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[kn*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[kn*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[kn*h_num+i].n0ij[2]);
			HYPER1[i*h_num+kn].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[kn*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[kn*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[kn*h_num+i].n0ij[2]);
			HYPER1[i*h_num+kn].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[kn*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[kn*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[kn*h_num+i].n0ij[2]);
			//cout<<"i"<<i<<"k"<<kn<<"	"<<HYPER1[kn*h_num+i].DgDq[A_X]<<","<<HYPER1[kn*h_num+i].DgDq[A_Y]<<","<<HYPER1[kn*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);

	}	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;


	/////////////W, StressåvéZ
	double trace_dC=0., trace_dC2=0.;

	double Ic=0., IIc=0.;

	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();

	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int i=0;i<h_num;i++)
	{
		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1./3.)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1./3.)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1./3.)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1./3.)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1./3.)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1./3.)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1./3.)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1./3.)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1./3.)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1./3.)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1./3.)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1./3.)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1./3.)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1./3.)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1./3.)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1./3.)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1./3.)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1./3.)*HYPER[i].Fi[2][2];
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

		trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		Ic=trace_dC;
		IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		HYPER[i].W=c10*(Ic-3)+c01*(IIc-3);

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
	
		HYPER[i].stress[0][0]=2./J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[i].stress[0][1]=2./J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[i].stress[0][2]=2./J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[i].stress[1][0]=2./J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[i].stress[1][1]=2./J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[i].stress[1][2]=2./J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[i].stress[2][0]=2./J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[i].stress[2][1]=2./J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[i].stress[2][2]=2./J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));

	}
	
		
}

void q_nab_lap(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *dE,double *dg,double Dt, double V, double mi,double nG[DIMENSION])
{
	/////////////////dEåvéZ///////////////////
	int h_num=HYPER.size();

	for(int k=0;k<h_num;k++)
	{
		dE[k]=0.;
		dE[k+h_num]=0.;
	}


	double Dgkk_n[DIMENSION]={0,0,0};
	double Dgki_n[DIMENSION]={0,0,0};


	for(int k=0;k<h_num;k++)
	{
		Dgkk_n[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];	Dgkk_n[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];	Dgkk_n[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];
		
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			Dgki_n[A_X]=HYPER1[k*h_num+in].DgDq_n[A_X];	Dgki_n[A_Y]=HYPER1[k*h_num+in].DgDq_n[A_Y];	Dgki_n[A_Z]=HYPER1[k*h_num+in].DgDq_n[A_Z];
			dE[k]-=0.5*Dt/mi*(Dgki_n[A_X]*HYPER[in].half_p[A_X]+Dgki_n[A_Y]*HYPER[in].half_p[A_Y]+Dgki_n[A_Z]*HYPER[in].half_p[A_Z]);
		}
		dE[k]-=0.5*Dt/mi*(Dgkk_n[A_X]*HYPER[k].half_p[A_X]+Dgkk_n[A_Y]*HYPER[k].half_p[A_Y]+Dgkk_n[A_Z]*HYPER[k].half_p[A_Z]);
	}

	double Dgik[DIMENSION]={0,0,0};
	double Dgik_n[DIMENSION]={0,0,0};
	double Dgkk[DIMENSION]={0,0,0};
	double Dgkk_n[DIMENSION]={0,0,0};
	double	p_Eik[DIMENSION]={0,0,0};
	double	p_Eik_n[DIMENSION]={0,0,0};
	double	p_Ekk[DIMENSION]={0,0,0};
	for(int k=0;k<h_num;k++)
	{
		Dgkk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	Dgkk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];
		p_Ekk[A_X]=HYPER[k].stress[A_X][0]*Dgkk[0]+HYPER[k].stress[A_X][1]*Dgkk[1]+HYPER[k].stress[A_X][2]*Dgkk[2];
		p_Ekk[A_Y]=HYPER[k].stress[A_Y][0]*Dgkk[0]+HYPER[k].stress[A_Y][1]*Dgkk[1]+HYPER[k].stress[A_Y][2]*Dgkk[2];
		p_Ekk[A_Z]=HYPER[k].stress[A_Z][0]*Dgkk[0]+HYPER[k].stress[A_Z][1]*Dgkk[1]+HYPER[k].stress[A_Z][2]*Dgkk[2];

		Dgkk_n[A_X]=HYPER1[k*h_num+k].DgDq_n[A_X];	Dgkk_n[A_Y]=HYPER1[k*h_num+k].DgDq_n[A_Y];	Dgkk_n[A_Z]=HYPER1[k*h_num+k].DgDq_n[A_Z];

		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			Dgik[A_X]=HYPER1[in*h_num+k].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+k].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+k].DgDq[A_Z];
			p_Eik[A_X]=HYPER[i].stress[A_X][0]*Dgik[0]+HYPER[i].stress[A_X][1]*Dgik[1]+HYPER[i].stress[A_X][2]*Dgik[2];
			p_Eik[A_Y]=HYPER[i].stress[A_Y][0]*Dgik[0]+HYPER[i].stress[A_Y][1]*Dgik[1]+HYPER[i].stress[A_Y][2]*Dgik[2];
			p_Eik[A_Z]=HYPER[i].stress[A_Z][0]*Dgik[0]+HYPER[i].stress[A_Z][1]*Dgik[1]+HYPER[i].stress[A_Z][2]*Dgik[2];

			Dgik_n[A_X]=HYPER1[in*h_num+k].DgDq_n[A_X];	Dgik_n[A_Y]=HYPER1[in*h_num+k].DgDq_n[A_Y];	Dgik_n[A_Z]=HYPER1[in*h_num+k].DgDq_n[A_Z];
			p_Eik_n[A_X]=HYPER[i].stress[A_X][0]*Dgik[0]+HYPER[i].stress[A_X][1]*Dgik[1]+HYPER[i].stress[A_X][2]*Dgik[2];
			p_Eik_n[A_Y]=HYPER[i].stress[A_Y][0]*Dgik[0]+HYPER[i].stress[A_Y][1]*Dgik[1]+HYPER[i].stress[A_Y][2]*Dgik[2];
			p_Eik_n[A_Z]=HYPER[i].stress[A_Z][0]*Dgik[0]+HYPER[i].stress[A_Z][1]*Dgik[1]+HYPER[i].stress[A_Z][2]*Dgik[2];

			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				dE[jn]+=0.5*Dt*Dt/mi*(p_Eik[A_X]*HYPER1[jn*h_num+k].DgDq_n[A_X]+p_Eik[A_Y]*HYPER1[jn*h_num+k].DgDq_n[A_Y]+p_Eik[A_Z]*HYPER1[jn*h_num+k].DgDq_n[A_Z]);
			}
			dE[k]+=0.5*Dt*Dt/mi*(p_Eik[A_X]*Dgkk_n[A_X]+p_Eik[A_Y]*Dgkk_n[A_Y]+p_Eik[A_Z]*Dgkk_n[A_Z]);
			dE[in]+=0.5*Dt*Dt/mi*(p_Ekk[A_X]*Dgik_n[A_X]+p_Ekk[A_Y]*Dgik_n[A_Y]+p_Ekk[A_Z]*Dgik_n[A_Z]);
		}
		dE[k]+=0.5*Dt*Dt/mi*(p_Ekk[A_X]*Dgkk_n[A_X]+p_Ekk[A_Y]*Dgkk_n[A_Y]+p_Ekk[A_Z]*Dgkk_n[A_Z]);
	}
	
	
	for(int k=0;k<h_num;k++)
	{
		dE[k+h_num]+=0.5*Dt/mi*(nG[A_X]*HYPER[k].half_p[A_X]+nG[A_Y]*HYPER[k].half_p[A_Y]+nG[A_Z]*HYPER[k].half_p[A_Z]);

		Dgkk[A_X]=HYPER1[k*h_num+k].DgDq[A_X];	Dgkk[A_Y]=HYPER1[k*h_num+k].DgDq[A_Y];	Dgkk[A_Z]=HYPER1[k*h_num+k].DgDq[A_Z];
		p_Ekk[A_X]=HYPER[k].stress[A_X][0]*Dgkk[0]+HYPER[k].stress[A_X][1]*Dgkk[1]+HYPER[k].stress[A_X][2]*Dgkk[2];
		p_Ekk[A_Y]=HYPER[k].stress[A_Y][0]*Dgkk[0]+HYPER[k].stress[A_Y][1]*Dgkk[1]+HYPER[k].stress[A_Y][2]*Dgkk[2];
		p_Ekk[A_Z]=HYPER[k].stress[A_Z][0]*Dgkk[0]+HYPER[k].stress[A_Z][1]*Dgkk[1]+HYPER[k].stress[A_Z][2]*Dgkk[2];
		
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			Dgik[A_X]=HYPER1[in*h_num+k].DgDq[A_X];	Dgik[A_Y]=HYPER1[in*h_num+k].DgDq[A_Y];	Dgik[A_Z]=HYPER1[in*h_num+k].DgDq[A_Z];
			p_Eik[A_X]=HYPER[i].stress[A_X][0]*Dgik[0]+HYPER[i].stress[A_X][1]*Dgik[1]+HYPER[i].stress[A_X][2]*Dgik[2];
			p_Eik[A_Y]=HYPER[i].stress[A_Y][0]*Dgik[0]+HYPER[i].stress[A_Y][1]*Dgik[1]+HYPER[i].stress[A_Y][2]*Dgik[2];
			p_Eik[A_Z]=HYPER[i].stress[A_Z][0]*Dgik[0]+HYPER[i].stress[A_Z][1]*Dgik[1]+HYPER[i].stress[A_Z][2]*Dgik[2];

			dE[k+h_num]-=0.5*Dt*Dt/mi*(p_Eik[A_X]*nG[A_X]+p_Eik[A_Y]*nG[A_Y]+p_Eik[A_Z]*nG[A_Z]);
		}
		dE[k+h_num]-=0.5*Dt*Dt/mi*(p_Ekk[A_X]*nG[A_X]+p_Ekk[A_Y]*nG[A_Y]+p_Ekk[A_Z]*nG[A_Z]);
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


////////Hyperelastic2Ç©ÇÁÇÃà¯åpÇ¨
void calc_W(mpsconfig &CON,vector<hyperelastic> &HYPER,int h_num)
{
	//	cout<<"íeê´É|ÉeÉìÉVÉÉÉãåvéZ";

	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();

	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int i=0;i<h_num;i++)
	{
		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1./3.)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1./3.)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1./3.)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1./3.)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1./3.)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1./3.)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1./3.)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1./3.)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1./3.)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1./3.)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1./3.)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1./3.)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1./3.)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1./3.)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1./3.)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1./3.)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1./3.)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1./3.)*HYPER[i].Fi[2][2];
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
	
		HYPER[i].stress[0][0]=2./J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[i].stress[0][1]=2./J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[i].stress[0][2]=2./J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[i].stress[1][0]=2./J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[i].stress[1][1]=2./J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[i].stress[1][2]=2./J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[i].stress[2][0]=2./J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[i].stress[2][1]=2./J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[i].stress[2][2]=2./J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));

	}
//	cout<<"----------OK"<<endl;
}

