
void q_QP(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t, double V, double mi, double Dt, double E0,vector<double > NEIw)
{
	////////////’è‹`///////////////
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

	double rg=0.1;
	double rh=10.;


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
	double p_p[DIMENSION]={0,0,0};
	double hp[DIMENSION]={0,0,0};
	double W=0.;
	double Dgki_n[DIMENSION]={0,0,0};
	double Dgii_n[DIMENSION]={0,0,0};
	double Dgji_n[DIMENSION]={0,0,0};
	double Dgli_n[DIMENSION]={0,0,0};
	double lam=0.;
	double mu=0.;

	////////////‰Šú‰»ŽZ///////////////
	for(int i=0;i<h_num;i++)
	{
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


	////////////dh, p_rEŒvŽZ///////////////
	for(int i=0;i<Nw;i++)
	{
		int iw=NEIw[i];
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






	////////////QP–@ŒvŽZ///////////////
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

			q_variables(CON,PART,HYPER,HYPER1);	//‚±‚±‚Ü‚Å1219
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
						if(h[i]+th_h[i]>0)	rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						else if(h[i]+th_h[i]>-1.e-20 && dh[i][k]*dh[i][l]>0)
						{
							rT[l*Nx+k]+=rh*dh[i][k]*dh[i][l];
						}
					}
				}
			}

			for(int i=0;i<Nx;i++)
			{
				d[i]=dT[i];
				for(int j=0;j<Nx;j++)
				{
					B[i*Nx+j]=rT[i*Nx+j];
				}
			}


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




			/////////////pn1_2ŒvŽZ
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
				/////////////qŒvŽZ
				PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
				PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
				PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
			}	
			//if(count==1||count%100==0)
			{
				cout<<"E"<<count<<"="<<E<<", En="<<En<<endl;
				output_data(PART,HYPER,HYPER1,T,dT,rT,g,dg,th_g,h,dh,th_h,d,count,count_min,t,E,En,E0,mi,V,NEIw);

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
		if(E_min<ep_min*1000)
		{
			rg*=4;
			rh*=4;
		}

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

	/////////////pn1_2ŒvŽZ
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
		/////////////qŒvŽZ
		PART[i].r[A_X]=HYPER[i].q_n[A_X]+Dt/mi*HYPER[i].half_p[A_X];
		PART[i].r[A_Y]=HYPER[i].q_n[A_Y]+Dt/mi*HYPER[i].half_p[A_Y];
		PART[i].r[A_Z]=HYPER[i].q_n[A_Z]+Dt/mi*HYPER[i].half_p[A_Z];	
	}	

}

