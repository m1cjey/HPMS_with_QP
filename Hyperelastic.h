#ifndef HYPERELASTIC
#define HYPERELASTIC

class hyperelastic
{
public:

	unsigned ID;
	int NEI[1000];
	int N;
	double pnd0;
	int flag_wall;

	double W;
	double lambda;
	double half_p[DIMENSION];
	double stress[DIMENSION][DIMENSION];
	double differential_p[DIMENSION];
	double p[DIMENSION];
	double ang_p[DIMENSION];
	double Ai[DIMENSION][DIMENSION];
	double inverse_Ai[DIMENSION][DIMENSION];
	double t_inverse_Ai[DIMENSION][DIMENSION];
	double t_inverse_Fi[DIMENSION][DIMENSION];
	double J;
	double pnd;
	double Fi[DIMENSION][DIMENSION];
	double vis_force[DIMENSION];


	//QP—p
	double pi[DIMENSION][DIMENSION];
	double lam;
	double h_lam;
	double mu;
	double h_mu;

	double old_lam;
	double old_h_lam;
	double old_mu;
	double old_h_mu;

	double W_n;
	double stress_n[DIMENSION][DIMENSION];
	double q_n[DIMENSION];
	double p_n[DIMENSION];
	double dpidF[DIMENSION][DIMENSION];
	hyperelastic();
};

class hyperelastic2
{
public:
	unsigned ID;
	double wiin;
	double DgDq[DIMENSION];
	double aiin[DIMENSION];
	double n0ij[DIMENSION];
	double spl_f;

	double DgDq_n[DIMENSION];
	double dPIdlam;
	double dPIdmu;
	double Dg_n_Dg_n;
	hyperelastic2();
};

#endif