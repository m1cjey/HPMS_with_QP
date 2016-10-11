#include "stdafx.h"
int MM_method();
void QP(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
double calc_W(double c10, double c01, double J,double Fi[3][3]);


void QP(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	int h_num=18;
	double mi=CON.get_hyper_density()*get_volume(&CON);
	double V=get_volume(&CON);
	double c10=CON.get_c10();
	double c01=CON.get_c01();

	double *W=new double [h_num];
	double *qx=new double [h_num];
	double *qy=new double [h_num];
	double *qz=new double [h_num];
	double *px=new double [h_num];
	double *py=new double [h_num];
	double *pz=new double [h_num];


	for(int i=0;i<h_num;i++)
	{
		W[i]=0;
		double J=HYPER[i].J;
		double Fi[3][3]={{HYPER[i].Fi[0][0], HYPER[i].Fi[0][1], HYPER[i].Fi[0][2]}, {HYPER[i].Fi[1][0], HYPER[i].Fi[1][1], HYPER[i].Fi[1][2]}, {HYPER[i].Fi[2][0], HYPER[i].Fi[2][1], HYPER[i].Fi[2][2]}};
		W[i]=calc_W(c10,c01,J,Fi);
	}

	double L=0;
	for(int i=0;i<h_num;i++)	L+=0.5/mi*(HYPER[i].p[A_X]*HYPER[i].p[A_X]+HYPER[i].p[A_Y]*HYPER[i].p[A_Y]+HYPER[i].p[A_Z]*HYPER[i].p[A_Z])-V*W[i];

	

	int Nx=6;
	double *x=new double [Nx];
	double *old_x=new double [Nx];
	x[0]=0;
	x[1]=0;

	double c0=0;
	double L=0;
	double Tr=0;

	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dTr=new double [Nx];

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
		
	
			fx=0.5/mi*(px*px+py*py+pz*pz)-V*W;

				(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
			dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

			c0=x[0]*x[0]+x[1]*x[1]-1;
			dc0[0]=2*x[0];	dc0[1]=2*x[1];

			Txr=fx+0.5*r*(c0+seta0)*(c0+seta0);			
			dTxr[0]=dfx[0]+r*(c0+seta0)*dc0[0];
			dTxr[1]=dfx[1]+r*(c0+seta0)*dc0[1];

			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
			cout<<"E0="<<E<<endl;
			if(E<ep)	break;
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
//				else if(beta>0)
//				{
//					double seta=0.8*sigma/(sigma-beta);
//					cout<<"seta="<<seta<<endl;
//					double y0_s=y[0]*seta+(1-seta)*(B[0]*s[0]+B[1]*s[1]);
//					double y1_s=y[1]*seta+(1-seta)*(B[2]*s[0]+B[3]*s[1]);

//					beta=y0_s*s[0]+y1_s*s[1];
//					double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
//					double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

//					B[0*Nx+0]+=1/beta*y0_s*y0_s-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
//					B[0*Nx+1]+=1/beta*y0_s*y1_s-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
//					B[1*Nx+0]+=1/beta*y1_s*y0_s-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
//					B[1*Nx+1]+=1/beta*y1_s*y1_s-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
//				}
			
		}

		E_min=sqrt( (old_x[0]-x[0])*(old_x[0]-x[0]) + (old_x[1]-x[1])*(old_x[1]-x[1]));
		if(E_min<ep*1000)	r*=4;
		seta0+=c0;

		cout<<"E_min"<<count_min<<"="<<E_min<<endl<<endl;
	}

	//v0　計算
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

	cout<<"理論解"<<endl;
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

}

int MM_method()	
{
	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
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
		int calc_way=1;	//最適解探索	Newton 0, 準Newton 1

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

	//v0　計算
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

	cout<<"理論解"<<endl;
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

double calc_W(double c10, double c01, double J,double Fi[3][3])
{
	//	cout<<"弾性ポテンシャル計算";
	double W=0;
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	
	if(J<0){
		d_Fi[0][0]=-1/pow(-J,1/3)*Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*Fi[0][2];
		d_Fi[1][0]=-1/pow(-J,1/3)*Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*Fi[1][2];
		d_Fi[2][0]=-1/pow(-J,1/3)*Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*Fi[2][2];
	}
	else
	{
		d_Fi[0][0]=1/pow(J,1/3)*Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*Fi[0][2];
		d_Fi[1][0]=1/pow(J,1/3)*Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*Fi[1][2];
		d_Fi[2][0]=1/pow(J,1/3)*Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*Fi[2][2];
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
	W=c10*(Ic-3)+c01*(IIc-3);

	return W;
}
