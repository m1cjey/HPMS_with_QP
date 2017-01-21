#include "stdafx.h"

void calc_vis_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,int t)
{
	cout<<"”S«€ŒvŽZ";
	double v=CON.get_vis();
	double r=CON.get_h_dis();
	int h_num=HYPER.size();
	int d=3;
	double mh=CON.get_hyper_density()*get_volume(&CON);
	double ms=CON.get_silicone_density()*get_volume(&CON);
	double lambda=calclambda(CON);

	stringstream ss;
	ss<<"./Viscosity/vis"<<t<<".csv";
	ofstream fv(ss.str());

	ofstream fv_sum("Vis.csv",ios::app);
	fv_sum<<t<<",";
	double v_sum=0;
	for(int i=0;i<h_num;i++)
	{
		double pnd0=HYPER[i].pnd0;

		double p_vis[3]={0,0,0};
		double pnd=0;
		double p_lambda=0;

		for(int j=0;j<h_num;j++)
		{
			if(j!=i)
			{
				double qij[3]={PART[j].r[A_X]-PART[i].r[A_X], PART[j].r[A_Y]-PART[i].r[A_Y], PART[j].r[A_Z]-PART[i].r[A_Z]};
				double dis=sqrt((PART[j].r[A_X]-PART[i].r[A_X])*(PART[j].r[A_X]-PART[i].r[A_X])+(PART[j].r[A_Y]-PART[i].r[A_Y])*(PART[j].r[A_Y]-PART[i].r[A_Y])+(PART[j].r[A_Z]-PART[i].r[A_Z])*(PART[j].r[A_Z]-PART[i].r[A_Z]));
				if(dis<r)
				{
					double w=kernel4(r,dis);

					if(PART[i].toFEM==1)
					{
						p_vis[A_X]+=(HYPER[j].p[A_X]-HYPER[i].p[A_X])/mh*w;
						p_vis[A_Y]+=(HYPER[j].p[A_Y]-HYPER[i].p[A_Y])/mh*w;
						p_vis[A_Z]+=(HYPER[j].p[A_Z]-HYPER[i].p[A_Z])/mh*w;
					}else
					{
						p_vis[A_X]+=(HYPER[j].p[A_X]-HYPER[i].p[A_X])/ms*w;
						p_vis[A_Y]+=(HYPER[j].p[A_Y]-HYPER[i].p[A_Y])/ms*w;
						p_vis[A_Z]+=(HYPER[j].p[A_Z]-HYPER[i].p[A_Z])/ms*w;
					}

					p_lambda+=(qij[A_X]*qij[A_X]+qij[A_Y]*qij[A_Y]+qij[A_Z]*qij[A_Z])*w;
					pnd+=w;
				}
			}
		}
		double lambda=p_lambda/pnd;

		HYPER[i].vis_force[A_X]=2*v*d/lambda/pnd0*p_vis[A_X];
		HYPER[i].vis_force[A_Y]=2*v*d/lambda/pnd0*p_vis[A_Y];
		HYPER[i].vis_force[A_Z]=2*v*d/lambda/pnd0*p_vis[A_Z];

		fv<<PART[i].r[A_X]<<","<<PART[i].r[A_Y]<<","<<PART[i].r[A_Z]<<","<<HYPER[i].vis_force[A_X]<<","<<HYPER[i].vis_force[A_Y]<<","<<HYPER[i].vis_force[A_Z]<<endl;
		fv_sum<<sqrt(HYPER[i].vis_force[A_X]*HYPER[i].vis_force[A_X]+HYPER[i].vis_force[A_Y]*HYPER[i].vis_force[A_Y]+HYPER[i].vis_force[A_Z]*HYPER[i].vis_force[A_Z])<<",";
		v_sum+=sqrt(HYPER[i].vis_force[A_X]*HYPER[i].vis_force[A_X]+HYPER[i].vis_force[A_Y]*HYPER[i].vis_force[A_Y]+HYPER[i].vis_force[A_Z]*HYPER[i].vis_force[A_Z]);
	}
	fv_sum<<v_sum<<endl;
	fv_sum.close();
	fv.close();

	cout<<"---------OK"<<endl;

}

void calc_spl_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic2>&HYPER1,int h_num)
{
	double h=CON.get_h_dis();
	double r=0,s=0,theta=0;

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			theta=0;	//ŠÖ”C³15/2/8
			HYPER1[i*h_num+j].spl_f=0;
			r=sqrt(pow((PART[j].r[A_X]-PART[i].r[A_X]),2.0)+pow((PART[j].r[A_Y]-PART[i].r[A_Y]),2.0)+pow((PART[j].r[A_Z]-PART[i].r[A_Z]),2.0));
			s=r/(h/2);
			if(0<=s&&s<1)	theta=1/PI*(1-3/2*pow(s,2.0)+3/4*pow(s,3.0));
			if(1<=s&&s<2)	theta=1/PI/4*pow((2-s),3.0);
			HYPER1[i*h_num+j].spl_f=1/pow(h/2,3)*theta;	//ŠÖ”C³15/2/8
		}
	}
}


