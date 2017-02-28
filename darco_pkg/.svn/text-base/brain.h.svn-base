// Mean Field Implementation
// 
// Marco Loh, loh@in.tum.de
// ---------with Facilitation----------------------------------------------------------
#ifndef _IF_BRAIN_H
#define _IF_BRAIN_H


#include <vector>
#include <gsl/gsl_rng.h>
#include <iostream>
#include "pool.h"

class brain {

	friend class pool;

      public:
	// Constructor
	 brain();
	// Destructor
	~brain();

	/*   int flagavr; */
	/*   int trials; */
	/*   double avr[maxpools][maxsteps]; */

	//  double pp0, pp1, pp2, pp3;

	// Initialiser
	void initBrain(int num_pools, int num_neurons, int num_ext_nrns,
		       double deltat);
	void set_rseed(unsigned long int seed);
	void initPool(int pool_id, double pool_rel_size, double startnu,
		      double startavrv, double ext_rate, poolinfo * d);

	// Calculation Methods
	void Euler(double eulerdelta, int clamp_pool);
	void RK2(double eulerdelta, int clamp_pool);
	int get_fixed_rates(double eulerdelta, int clamp_pool);
	void MeanField(double eulersteps, double eulerdelta,
		       ostream * rate_output);
	int MeanFieldEffective(double, int, double, int, double);
	static double drand();

	// Pops out (and empties) the value of spikes[i]
	double pop_spikes(int pool_id);
	double pop_ufacil(int pool_id); //Facilitation
	// Outputs spiking data to file and to stdout
	void write_down(ofstream * rateoutput);
	void write_down(ofstream * rateoutput, ofstream * uoutput);
	void cout_rates();
	void print_rates(ofstream * file);
	void write_rates(ofstream * file);

	// Sets the dimensions of the avg_on_trials_matrix and others
	void init_vector_matrices(int num_bins, vector< vector<double> > &m);
	void initialize_vectors(int t);
	static int l;

	void set_window_size(int i);
	void set_window_step(int i);
	void set_time_forward(bool b);
	int get_window_size();
	int get_window_step();
	bool get_time_forward();
	
	void Spiking(int offset_time, int total_steps,
		     ofstream * rate_output, ofstream * u_output,
		     ofstream * spikes_file,
		     ofstream * snmda1, ofstream * snmda2,
                     ofstream * snmda3, bool facilON, bool self_finish);
	void resetPools();
	//  void resetavr ();

	// Adaptation Methods
	void setConnection(int pool_origin, int pool_target, double ampa,
			   double nmda, double gaba);
	void setStartPoolrate(int p, double nu);
	void setStartAvrV(int p, double v);

	double get_g_NMDA(int p);
	double get_g_AMPArec(int p);
	void set_g_NMDA(int p, double g);
	void set_g_AMPArec(int p, double g);

	// void setPoolsize(int p, double fr);
	void setExtrate(int p, double nu);
	double getExtrate(int p);
	
	double getNMDAerror();
	double getPoolrate(int p);
	void setPoolrate(int p, double nu);
	double getConnection(int pool_origin, int pool_target, char art);
	double getPoolsize(int p);

	double getFlow(int);

	void calc_mf_parameters();
	void printBrain();
	void print_pools();
	void printPool(int p);

      private:
	void setAvrV(int p, double v);

	static double poolnu[maxpools];	// pool frequencies
	static double poolf[maxpools];	// fraction of pool size of
	// pools, must add up to 1.

	static int C;		// connections to exc. cells
	static int Cext;	// external connections
	static int PoolCount;

	bool filled; // to check if the temp_matrix (a queue)
		     // has already been filled.

	static bool time_forward;
	int window_size;
	int window_step;
	
	static double dt;
	static double sumAMPArec[maxpools];
	static double sumNMDA[maxpools];
	static double sumGABA[maxpools];
	static int spikes[maxpools];
	static double sumup[maxpools];	//Facilitation

	class pool pools[maxpools];	// Pools in the brain
	double nmdaerror;	// error calculation

	// output spikes
	ostream *rateoutput;
	ostream *uoutput;

	// The matrix containing the mean spiking rate;
	// rows represent different pools and columns represent time
	vector< vector<double> > avg_on_trials;
	vector< vector<double> > avg_sumu_on_trials; //Facilitation

	// The rates recorded in intervals of window_step (2ms), which
	// will be averaged over the whole range (window_size, 50ms)
	vector< vector<double> > matrix_rates;
	vector< vector<double> > matrix_uout; //Facilitation
	
	vector<double> means;
	vector<double> means_ufacil; //Facilitation

	// Handler for the GSL random number generator
	static gsl_rng *rnd_instance;
};
#endif
