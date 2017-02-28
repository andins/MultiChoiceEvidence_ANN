
// ---------with Facilitation---------------------------------------------------------- 
#ifndef _IF_POOL_H
#define _IF_POOL_H

using namespace std;

// ------------------------------------------------------------
//    Constants
// ------------------------------------------------------------
const int maxpools = 20;
const int maxpoolneurons = 4000;//2600; yota's problem:  con piú di 5000 neurns i rates vanno a saturazione
const int maxsteps = 20000;

// Precalculation stuff
const int max_size = 15000;
const int max_broadsize = 300; // Should be a divisor of max_size

const int NMAX = 7;
const double PI = 3.1415926535;

const double a1 = -1.26551223;
const double a2 = 1.00002368;
const double a3 = .37409196;
const double a4 = .09678418;
const double a5 = -.18628806;
const double a6 = .27886087;
const double a7 = -1.13520398;
const double a8 = 1.48851587;
const double a9 = -.82215223;
const double a10 = .17087277;

// Calcium dynamics
const double VK = -80.;		// (mV) reversal Ca2+
const double gAHP = 7.5;	// (nS , 0.015 mS/cm2) 7.5
const double alphaCa = 0.1;	// (0.2 muM) increment Ca2+  0.005
const double taoCa = 50.;	// (ms) decay Ca2+ 600
const int AHP_Ca = 0;		// flag
const double taoF = 2000;	//Facilitation timeconst
const double U = 0.15;		//Facilitation (Gustavo 0.15)
// ------------------------------------------------------------

struct connweights {
	double ampa, nmda, gaba;
};

struct poolinfo {
	double VL, Vthr, Vreset;	// Resting potential, firing threshold, reset potential
	double VI, VE;		// constants?!

	double Cm;		// Membrane capacitance
	double gm;		// Membrane leak conductance
	double taurp;		// Refractory period
	double ustart;		// starting point of facilitation variable MARIO: see setup.cc, pool.cc.
	double taum;		// time constant
	double calpha, cbeta, cgamma;	// constant coefficients

// synaptic conductances
	double gAMPAext, gAMPArec, gNMDA, gGABA;
// gating variables
	double tauAMPA, tauGABA, tauNMDAdecay, tauNMDArise;
};


class pool {
	friend class brain;

      public:
	// configuration
	 pool();
	~pool();
	void initPool(int id, poolinfo * data);
	void setExternal(double nu);
	// get some parameters
	double get_gAMPArec();
	double get_gNMDA();
	void set_gAMPArec(double g);
	void set_gNMDA(double g);

	// Spiking Functions
	void resetSpikingPool();
	void calcSynapses(double time, bool snmda_f, ofstream * snmda_file);
	void calcPotentials(double time, ofstream * file, bool facilON);

	// Meanfield Functions
	double Euler_derivative(double);
	double Flow();
	void newavrV();

      private:
/* 	gsl_rng* rnd_instance; // handler for the random number generator */
	// Surrouding World
	 connweights w[maxpools];	// Incoming Connections
	int pool_id;		// Pool ID in Array of brain
	int num_neurons;	// Number of neurons
	double nuext;		// External Input

	// Pool Specific Variables
	double VL, Vthr, Vreset;	// Resting potential, firing threshold, reset
	double VI, VE;		// constants?!
	double Cm;		// Membrane capacitance
	double gm;		// Membrane leak conductance
	double taurp;		// Refractory period
	double taum;		// time constant
	double calpha, cbeta, cgamma;	// constants
	double gAMPAext, gAMPArec, gNMDA, gGABA;
	double tauAMPA, tauGABA, tauNMDAdecay, tauNMDArise;

	double ustart;		// starting point of facilitation variable MARIO: see setup.cc, pool.cc.
	// Spiking Variables
	double sAMPAext[maxpoolneurons];
	double sAMPArec[maxpoolneurons];
	double sNMDA[maxpoolneurons];
	double xNMDA[maxpoolneurons];
	double sGABA[maxpoolneurons];
	double V[maxpoolneurons];
	double dV[maxpoolneurons];
	double Ca[maxpoolneurons];
	double u[maxpoolneurons]; //Facilitation
        double sumu; 		  //Facilitation
	double lastspiketime[maxpoolneurons];
	double expdttauAMPA, expdttauNMDArise, expdttauGABA;
	
	// Spiking functions
	int extspike();
	double drand();

	// Mean Field Functions
	void initVar();
	void generateintegraldata();

	double Phi();		// old
	double Phi2();		// new
	double nerf(double z);
	double falpha();
	double fbeta();
	double tauE();
	double muE();
	double sigmaE();
	double Sx();
	double rho1();
	double rho2();
	double nx();
	double Nx();
	double nix();
	double psi(double nu);	// old
	double psi2(double nu);	// new
	double J();
	double T(int n);

	// Mean Field Variables
	double startnu;
	double startavrV;
	double avrV;
	double fak2[NMAX + 1];
	double bin[NMAX + 1][NMAX + 1];
	double alphatauNMDArisen[NMAX + 1];
	double TEext, TEAMPA, TEI;
	double crho1, crho2;
	double tauNMDA, onenutauNMDA;
	double vnx, vNx, vnix;
	double vSx, vtauE, vmuE;
	double vsigmaE, csigmaE;

	double psidata[max_size];
	double intdata[max_size];
	double intdatabroad[max_broadsize];
};

#endif
