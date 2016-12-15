#ifndef HYPERELASTIC
#define HYPERELASTIC

class hyperelastic
{
public:

	unsigned ID;
	int NEI[60];
	int N;
	double pnd0;

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
	//double pi[DIMENSION][DIMENSION];
	//double pi_n[DIMENSION][DIMENSION];
	double lam;
	double mu;
	double E;

	double old_lam;
	double old_mu;

	double W_n;
	double stress_n[DIMENSION][DIMENSION];
	double q_n[DIMENSION];
	double ph_n[DIMENSION];
	double p_n[DIMENSION];

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
	double DpiDq[DIMENSION];
	hyperelastic2();
};

#endif