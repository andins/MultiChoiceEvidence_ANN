// Mean Field Implementation
// Marco Loh, loh@in.tum.de
// ---------with Facilitation----------------------------------------------------------
#include <iomanip>
#include <fstream>
#include <math.h>
#include <climits>
#include "pool.h"
#include "brain.h"

int brain::PoolCount = 0;
int brain::C = 0;
int brain::Cext = 0;
int brain::l = 0; // iterates over time steps for population rates
bool brain::time_forward = true;
double brain::poolnu[maxpools];
double brain::poolf[maxpools];	// coding_lvl (exc. + inh. !)
double brain::dt = 0;
double brain::sumAMPArec[maxpools];
double brain::sumNMDA[maxpools];
double brain::sumGABA[maxpools];
double brain::sumup[maxpools]; //Facilitation
int brain::spikes[maxpools];
gsl_rng *brain::rnd_instance;

// constuctor

brain::brain()
{
	const gsl_rng_type *type_rng;

	gsl_rng_env_setup();
	type_rng = gsl_rng_default;
	rnd_instance = gsl_rng_alloc(type_rng);
}

brain::~brain()
{
	gsl_rng_free(rnd_instance);
}

// void brain::resetavr(void)
// {
//      for (int i = 0; i < maxpools; i++)
//              for (int j = 0; j < maxsteps; j++)
//                      avr[i][j] = 0;
// }

void brain::initBrain(int num_pools, int number_neurons, int num_ext_nrns,
		double deltat)
{
	PoolCount = num_pools;
	C = number_neurons;
	Cext = num_ext_nrns;
	rateoutput = 0;		// null pointer
	uoutput = 0;		// null pointer

	// time step for spiking neurons!
	dt = deltat;
	int i, j;
	for (i = 0; i < num_pools; i++)	// before: maxpools
		for (j = 0; j < num_pools; j++) {
			pools[j].w[i].ampa = 0;
			pools[j].w[i].nmda = 0;
			pools[j].w[i].gaba = 0;
		}
}

void brain::set_rseed(unsigned long int seed){
	gsl_rng_set(rnd_instance, seed);
}

void brain::initPool(int p, double p_fraction, double startnu,
		double startavrv, double extrate, poolinfo * d)
{
	poolf[p] = p_fraction;
	pools[p].startnu = startnu; //for meanfield
	pools[p].startavrV = startavrv; //for meanfield
	pools[p].setExternal(extrate);
	pools[p].initPool(p, d);
}

void brain::resetPools()
{
	for (int j = 0; j < PoolCount; j++) {
		setPoolrate(j, pools[j].startnu); //for meanfield
		setAvrV(j, pools[j].startavrV); //for meanfield
		pools[j].resetSpikingPool();
		l = 0;
		filled = false;
	}
}

void brain::setConnection(int pool_origin, int pool_target, double ampa,
		double nmda, double gaba)
{
	pools[pool_target].w[pool_origin].ampa = ampa;
	pools[pool_target].w[pool_origin].nmda = nmda;
	pools[pool_target].w[pool_origin].gaba = gaba;
}

double brain::getConnection(int pool_origin, int pool_target, char art)
{
	switch (art) {
		case 'a':
			return (pools[pool_target].w[pool_origin].ampa);
		case 'n':
			return (pools[pool_target].w[pool_origin].nmda);
		case 'g':
			return (pools[pool_target].w[pool_origin].gaba);
		default:
			return (99999);	// error, you will notice ;)
	}
}

void brain::setPoolrate(int p, double nu) //for meanfield
{
	poolnu[p] = nu;
}

void brain::setAvrV(int p, double v) //for meanfield
{
	pools[p].avrV = v;
}

void brain::setExtrate(int p, double nu)
{
	pools[p].setExternal(nu);
}

double brain::getExtrate(int p)
{
	return pools[p].nuext;
}

/*void brain::setPoolsize(int p, double fr){
  poolf[p]=fr;
  }*/
void brain::setStartPoolrate(int p, double nu)
{
	pools[p].startnu = nu;
}

void brain::setStartAvrV(int p, double v) //for meanfield
{
	pools[p].startavrV = v;
}

double brain::get_g_NMDA(int p)
{
	return pools[p].get_gNMDA();
}

double brain::get_g_AMPArec(int p)
{
	return pools[p].get_gAMPArec();
}

void brain::set_g_NMDA(int p, double g)
{
	pools[p].set_gNMDA(g);
}

void brain::set_g_AMPArec(int p, double g)
{
	pools[p].set_gAMPArec(g);
}

void brain::set_window_size(int i)
{
	window_size = i;
}
void brain::set_window_step(int i)
{
	window_step = i;
}
void brain::set_time_forward(bool b) 
{
    time_forward = b;
}
int brain::get_window_size()
{
	return window_size;
}
int brain::get_window_step()
{
	return window_step;
}
bool brain::get_time_forward()
{
    return time_forward;
}
double brain::getPoolrate(int p)
{
	return poolnu[p];
}

double brain::getPoolsize(int p)
{
	return (poolf[p]);
}

double brain::getNMDAerror()
{
	return (nmdaerror);
}

//double brain::drand()  // old function by dani/marco
//{
//	double uniform_rnd =
//		static_cast<double> (gsl_rng_get(rnd_instance)) / 
//		static_cast<double> (gsl_rng_max(rnd_instance));
//	return (uniform_rnd);
//}

double brain::drand() // new function by miguel
{
	static unsigned long x=123456789, y=362436069, z=521288629;
	unsigned long t;
	x ^= x << 16;
	x ^= x >> 5;
	x ^= x << 1;
	
	t = x;
	x = y;
	y = z;
	z = t ^ x ^ y;

	return (double)z/ULONG_MAX;
}

double brain::pop_spikes(int pool_id)
{
	double f = static_cast<double> (spikes[pool_id])/static_cast<double> (pools[pool_id].num_neurons);
	spikes[pool_id] = 0;	// reset the spike counter
	return (f);
}

double brain::pop_ufacil(int pool_id)	//Facilitaion
{
	double f = static_cast<double> (sumup[pool_id])/static_cast<double> (pools[pool_id].num_neurons);
	sumup[pool_id] = 0;	// reset the spike counter
	return (f);
}

void brain::write_down(ofstream * rateoutput)
{	
	int bins_window = window_size / window_step;
	for (int j = 0; j < (l - bins_window) + 1; j++) {
		*rateoutput << fixed << setprecision(5);
		*rateoutput << setw(6) << j * window_step + window_size;
		for (int i = 0; i < PoolCount; i++) {
			*rateoutput << fixed << setprecision(5) << setw(10) << "  " << avg_on_trials[i][j];
			//*rateoutput << setw(10);
			//*rateoutput << "  " << avg_on_trials[i][j];
		}
		*rateoutput << endl;
	}
}

void brain::write_down(ofstream * rateoutput, ofstream * uoutput)
{	
	int bins_window = window_size / window_step;
	for (int j = 0; j < (l - bins_window) + 1; j++) {
		*rateoutput << fixed << setprecision(5);
		*rateoutput << setw(6) << j * window_step + window_size;
		*uoutput << fixed << setprecision(5);
		*uoutput << setw(6) << j * window_step + window_size;
		for (int i = 0; i < PoolCount; i++) {
			*rateoutput << fixed << setprecision(5) << setw(10) << "  " << avg_on_trials[i][j];
			*uoutput << fixed << setprecision(5) << setw(10) << "  " << avg_sumu_on_trials[i][j];
			//*rateoutput << setw(10);
			//*rateoutput << "  " << avg_on_trials[i][j];
		}
		*rateoutput << endl;
		*uoutput << endl;
	}
}

void brain::cout_rates()
{
	for (int j = 0; j < PoolCount; j++) {
		cout << fixed << setprecision(6);
		cout << setw(12) << "  " << getPoolrate(j);
	}
	cout << "\n";
}

void brain::write_rates(ofstream * file)
{
	for (int j = 0; j < PoolCount; j++) {
		*file << fixed << setprecision(6);
		*file << setw(12) << "  " << getPoolrate(j);
	}
	*file << "\n";
}

void brain::print_rates(ofstream * file) 
{
	for (int j = 0; j < PoolCount; j++) {
		*file << fixed << setprecision(6);
		*file << setw(12) << "  " << getPoolrate(j);
	}
	*file << "\n";
}

void brain::init_vector_matrices(int n, vector< vector<double> > &m)
{
	m.clear();
	m.resize(PoolCount);	// set the number of rows
	for (int i = 0; i < PoolCount; i++) {
		// set the number of columns and initialize
		// all elements to 0;
		m[i].resize(n); }
}

void brain::initialize_vectors(int totalspan)
{

	init_vector_matrices((totalspan - window_size) / window_step + 1, avg_on_trials);
	init_vector_matrices((totalspan - window_size) / window_step + 1, avg_sumu_on_trials); //Facilitation
	init_vector_matrices(window_size / window_step,	matrix_rates);
	init_vector_matrices(window_size / window_step,	matrix_uout); //Facilitation
	means.clear();
	means.resize(PoolCount, 0.0);
	means_ufacil.clear(); //Facilitation
	means_ufacil.resize(PoolCount, 0.0);
}

void brain::Spiking(int offset_time, int timespan, ofstream * file, ofstream * ufile,
		ofstream * spikes_file, ofstream * snmda1, ofstream * snmda2,
                ofstream * snmda3, bool facilON, bool b_self_finish=false)
{
    rateoutput = file;
    uoutput = ufile;
    // The total number of computational steps
    int total_steps = (int) (timespan / dt + 0.5);
    // The step size 
    int steps_bin = (int) (window_step / dt + 0.5);
    // The number of bins we average over
    int bins_window = window_size / window_step;
    // Output the S_{nmda} variable?
    bool snmda_f = false;
    // time in miliseconds
    double time_ms;
    // Selectivity index and relatives
    double sel_ind, rate_1, rate_2;
    int c = 0; // [sel index] counts number of steps over threshold.

    for (int j = 0; j < total_steps; j++) {
	time_ms = (j + 1) * dt + offset_time;

	double inc_status = 10.; // ANDREA: incremental time to show time flow update
	if ( fmod(time_ms, inc_status) == 0 ){
		fprintf(stderr, "Network Time: %.0f ms\r", time_ms);	// show network time flow
	}

	for (int i = 0; i < PoolCount; i++) {
	// TODO: for more than 4 pools
		if (i == 1){
		snmda_f = true;
		pools[i].calcSynapses(time_ms, snmda_f, snmda1);
	    } else if (i == 2) {
		snmda_f = true;
		pools[i].calcSynapses(time_ms, snmda_f, snmda2);
	    } else if (i == 3) {
		snmda_f = true;
		pools[i].calcSynapses(time_ms, snmda_f, snmda3);
	    } else {
		snmda_f = false;
		// Brutto: snmda1 is passed to the function
		// but not used whatsoever.
		pools[i].calcSynapses(time_ms, snmda_f, snmda1);
	    }
	}
	for (int i = 0; i < PoolCount; i++) {
	    pools[i].calcPotentials(time_ms, spikes_file, facilON);
	}
	// pop out the spike value at the end of the bin.
	if ((j + 1) % steps_bin == 0) {
	    // Push the new rates to the queue
	    for (int i = 0; i < PoolCount; i++) {
		for (int m = 0; m < bins_window - 1; m++) {
		    matrix_rates[i][m] = matrix_rates[i][m + 1];
		    matrix_uout[i][m] = matrix_uout[i][m + 1];	
		}
		matrix_rates[i][bins_window - 1] = pop_spikes(i) * (1000.0 / (double)window_step);
		matrix_uout[i][bins_window - 1] = pop_ufacil(i) * (dt / (double)window_step);
	    }
	    if( (!filled &&	(l + 1) >= bins_window) || filled) {
		filled = true;
		// Calculate the mean rate within window
		for (int i = 0; i < PoolCount; i++) {
		    means[i] = 0.0;
		    means_ufacil[i] = 0.0;	
		    for (int m = 0; m < bins_window; m++) {
			means[i] += matrix_rates[i][m] / bins_window;
			means_ufacil[i] += matrix_uout[i][m] / bins_window; //Facilitation
		    }
		}	
		for (int i = 0; i < PoolCount; i++) {
		    avg_on_trials[i][l - bins_window + 1] += means[i];
		    avg_sumu_on_trials[i][l - bins_window + 1] += means_ufacil[i];	//Facilitation
		}
		// For decision-making networks it makes no sense
		// to keep simulating after the decision has been made.
		// Here we stop the trial after making a choice.
		// TODO allow the user to set the parameters used in
		//	    the criterion.
		if(b_self_finish){
//		    TODO: max vs next? for (int i = 0; i < PoolCount-1; i++){
		    rate_1 = avg_on_trials[1][l - bins_window + 1];
		    rate_2 = avg_on_trials[2][l - bins_window + 1];
//		    }
		    sel_ind = fabs(rate_1 - rate_2) / (rate_1 + rate_2);
		    if (sel_ind > 0.5){
			c++;
			// sel_ind should be above 0.5 for at least 250ms
			if (c > 300 / window_step){
			    l++;
			    break; // Breaks the "for j<total_steps" loop
			}
		    } else {
			// Start again
			c = 0;
		    }
		}
	    }
	    l++;
	} // end of if ((j + 1) % ...
    } // end of trial
}

//==================================================================================================================
// meanfield section
//==================================================================================================================

void brain::Euler(double eulerdelta, int clamp_pool)
{
    double tempnu[maxpools];
    for (int j = 0; j < PoolCount; j++) {
		if (clamp_pool == -1 || (j != clamp_pool)) {
			tempnu[j] = getPoolrate(j) + pools[j].Euler_derivative(eulerdelta);
			setPoolrate(j, tempnu[j]);
			//cout << "Pool" << j << " Rate: " << tempnu[j] << "\n";
		}
	}
	// Updating avrV
	for (int j = 0; j < PoolCount; j++) {
		pools[j].newavrV();
		// Warning
		if (pools[j].avrV < -56.)
			nmdaerror += pow(-56. - pools[j].avrV, 2);
		if (pools[j].avrV > -48.)
			nmdaerror += pow(-48. - pools[j].avrV, 2);
	}
}

void brain::RK2(double eulerdelta, int clamp_pool)
{
	// Second order Runge-Kutta method
	double initialnu[maxpools], tempnu[maxpools];

	// Intermediate step
	// -----------------
	for (int j = 0; j < PoolCount; j++) {
		if (clamp_pool == -1 || (j != clamp_pool)) {
			initialnu[j] = getPoolrate(j);
			tempnu[j] = initialnu[j] +pools[j].Euler_derivative(eulerdelta) / 2.0;
		}

	}
	// Setting rates for the intermediate iteration
	for (int j = 0; j < PoolCount; j++) {
		if (j != clamp_pool)
			setPoolrate(j, tempnu[j]);
	}
	// Updating avrV
	for (int j = 0; j < PoolCount; j++) {
		pools[j].newavrV();
		// Warning
		if (pools[j].avrV < -56.)
			nmdaerror += pow(-56. - pools[j].avrV, 2);
		if (pools[j].avrV > -48.)
			nmdaerror += pow(-48. - pools[j].avrV, 2);
	}
	// Final step
	// ----------
	for (int j = 0; j < PoolCount; j++) {
		if (j != clamp_pool) {
			tempnu[j] = initialnu[j] +
				pools[j].Euler_derivative(eulerdelta);
		}
	}
	// Setting rates for the final iteration
	for (int j = 0; j < PoolCount; j++) {
		if (j != clamp_pool) {
			setPoolrate(j, tempnu[j]);
		}
	}
	// Updating avrV
	for (int j = 0; j < PoolCount; j++) {
		pools[j].newavrV();
		// Warning
		if (pools[j].avrV < -56.)
			nmdaerror += pow(-56. - pools[j].avrV, 2);
		if (pools[j].avrV > -48.)
			nmdaerror += pow(-48. - pools[j].avrV, 2);
	}
}

int brain::get_fixed_rates(double eulerdelta, int clamp_pool)
{
	int convergence = 0;
	int eulercount = 0;
	int n_queue = 4;

	// An array of queues, one for each pool
	double **temp_matrix = new double *[n_queue];
	for (int k = 0; k < n_queue; k++)
		temp_matrix[k] = new double[PoolCount];
	for (int k = 0; k < n_queue; k++)
		for (int j = 0; j < PoolCount; j++)
			temp_matrix[k][j] = 0.0;

	nmdaerror = 0;

	double *mean_rates = new double[PoolCount];

	while (convergence == 0) {
		eulercount++;

		// Since trajectories are not really relevant and we
		// are only interested in the fixed points, the simple
		// and imprecise Euler method will suffice.[?]
//L: 21.05.09	Laura Dempere found out, that there are numerical problems with the Euler at around 2400 Hz ext. Input (just the value we use!!)
//		--> RK2 not Euler!
		RK2(eulerdelta, clamp_pool);
		//Euler(eulerdelta, clamp_pool);

		// Every 10 (500) steps we update the array of queues,
		// and check if there is much variation among elements
		if (eulercount > 5000 and eulercount % 100 == 0) {	// 500
			// Push the new rates to the queue
			for (int k = 0; k < n_queue - 1; k++) {
				for (int j = 0; j < PoolCount; j++)
					temp_matrix[k][j] =
						temp_matrix[k + 1][j];
			}
			for (int j = 0; j < PoolCount; j++)
				temp_matrix[n_queue - 1][j] =
					getPoolrate(j);

			// Variability in each queue.
			// first, calculate the averages for every pool
			for (int j = 0; j < PoolCount; j++)
				mean_rates[j] = 0.0;

			for (int j = 0; j < PoolCount; j++) {
				double sum = 0.0;
				for (int k = 0; k < n_queue; k++)
					sum += temp_matrix[k][j];
				mean_rates[j] = sum / double (n_queue);
				if(time_forward==false)
				    cout << mean_rates[j] << "  ";
			}
			if(time_forward==false) {
			    cout << "\n";
			    cout.flush();
			}

			// Now check if ANY of the elements differs
			// from the mean in less than 1e-6
			int nout = 0;

			// if the following test is valid, the
			// fixed_point flag will be turned on
			//                         int fixed_point = 0;
			for (int j = 0; j < PoolCount; j++) {
				for (int k = 0; k < n_queue; k++) {
					if (fabs(temp_matrix[k][j] - mean_rates[j]) > 1.0e-8) {
						nout = 1;
						break;
					}
				}
				if (nout == 1)
					break;	// last exit
			}
			if (nout == 0)	// The test has succeeded
				convergence = 1;
		}
	}
	//         cout << fixed << setprecision(1);
	//         cout << setw(6) << eulercount * eulerdelta;
	//         for (int j = 1; j < PoolCount - 1; j++) { // Selective pools only
	//                 cout << fixed << setprecision(3);
	//                 cout << setw(8) << "  " << getPoolrate(j) * 1000.0;
	//              }
	//         cout << endl;

	delete[]mean_rates;

	for (int k = 0; k < n_queue; k++)
		delete[]temp_matrix[k];
	delete[]temp_matrix;

	return eulercount;
}

void brain::MeanField(double eulersteps, double eulerdelta, ostream * file)
{
	//      char d;
	int ok = 0;
	int eulercount = 0;
	nmdaerror = 0;

	rateoutput = file;
	cout_rates();
	while (ok == 0) {

		// loop termination (now like a "for" loop)
		eulercount++;
		if (eulercount == eulersteps)
			ok = 1;

		// Calculation
		RK2(eulerdelta, -1);
		//Euler(eulerdelta, -1);
		//Check for errors
		/*if (eulercount < 10 || eulercount % 100 == 0){
			cout_rates();
		}*/
			
		// Updating avrV
		for (int j = 0; j < PoolCount; j++) {
			pools[j].newavrV();
			// Warning
			if (pools[j].avrV < -56.)
				nmdaerror += pow(-56. - pools[j].avrV, 2);
			if (pools[j].avrV > -48.)
				nmdaerror += pow(-48. - pools[j].avrV, 2);
			// if (pools[j].avrV<-56 || pools[j].avrV>-48)
			//     cout <<"!" << tempnu[j] << "|"<<pools[j].avrV << "!";
		}
		// cout << endl;
	}
}

int brain::MeanFieldEffective(double eulerdelta, int pool_fix1,
		double rate_u1, int pool_fix2,
		double rate_u2)
{
	resetPools();

	double tempnu[maxpools];

	int convergence = 0;
	int eulercount = 0;
	int n_queue = 4;

	nmdaerror = 0;

	// An array of queues, one for each pool
	double **temp_matrix = new double *[n_queue];
	for (int k = 0; k < n_queue; k++)
		temp_matrix[k] = new double[PoolCount];
	for (int k = 0; k < n_queue; k++)
		for (int j = 0; j < PoolCount; j++)
			temp_matrix[k][j] = 0.0;

	// Quenched values
	setPoolrate(pool_fix1, rate_u1);
	setPoolrate(pool_fix2, rate_u2);

	double *means = new double[PoolCount];

	while (convergence == 0) {
		eulercount++;

		//                 cout << setprecision(0);
		//                 cout << setw(5) << eulercount * eulerdelta;
		for (int j = 0; j < PoolCount; j++) {
			tempnu[j] =
				getPoolrate(j) +
				pools[j].Euler_derivative(eulerdelta);
			//                         cout << fixed << setprecision(2);
			//                         cout << setw(8) << "  " << tempnu[j] * 1000;
		}
		//              cout << "[" << eulercount << "] ";

		// Setting rates for new iteration + output
		for (int j = 0; j < PoolCount; j++) {
			if (j != pool_fix1 && j != pool_fix2) {
				setPoolrate(j, tempnu[j]);
			}
			setStartPoolrate(j, tempnu[j]);
		}

		// Updating avrV
		for (int j = 0; j < PoolCount; j++) {
			pools[j].newavrV();
			// Warning
			if (pools[j].avrV < -56.)
				nmdaerror += pow(-56. - pools[j].avrV, 2);
			if (pools[j].avrV > -48.)
				nmdaerror += pow(-48. - pools[j].avrV, 2);
			// if (pools[j].avrV<-56 || pools[j].avrV>-48)
			//     cout <<"!" << tempnu[j] << "|"<<pools[j].avrV << "!";
		}
		//--------------------------------------------------
		// Convergence test
		//
		if (eulercount > 2000 and eulercount % 50 == 0) {	// 500
			// Push the new rates to the queue
			for (int k = 0; k < n_queue - 1; k++) {
				for (int j = 0; j < PoolCount; j++)
					temp_matrix[k][j] =
						temp_matrix[k + 1][j];
			}
			for (int j = 0; j < PoolCount; j++)
				temp_matrix[n_queue - 1][j] =
					getPoolrate(j);

			// Variability in each queue.
			// first, calculate the averages for every pool
			for (int j = 0; j < PoolCount; j++)
				means[j] = 0.0;

			for (int j = 0; j < PoolCount; j++) {
				double sum = 0.0;
				for (int k = 0; k < n_queue; k++)
					sum += temp_matrix[k][j];
				means[j] = sum / double (n_queue);
			}

			// Now check if ANY of the elements differs
			// from the mean in less than 1e-6
			int nout = 0;

			// if the following test is valid, the
			// fixed_point flag will be turned on
			//                         int fixed_point = 0;
			for (int j = 0; j < PoolCount; j++) {
				for (int k = 0; k < n_queue; k++) {
					if (fabs
							(temp_matrix[k][j] -
							 means[j]) > 1.0e-8) {
						nout = 1;
						break;
					}
				}
				if (nout == 1)
					break;	// last exit
			}
			if (nout == 0)	// The test has succeeded
				convergence = 1;
		}
	}
	//         cout << fixed << setprecision(1);
	//         cout << setw(6) << eulercount * eulerdelta;
	//         for (int j = 1; j < PoolCount - 1; j++) { // Selective pools only
	//                 cout << fixed << setprecision(3);
	//                 cout << setw(8) << "  " << getPoolrate(j) * 1000.0;
	//              }
	//         cout << endl;

	delete[]means;

	for (int k = 0; k < n_queue; k++)
		delete[]temp_matrix[k];
	delete[]temp_matrix;

	return eulercount;
}


double brain::getFlow(int p)
{
	return pools[p].Flow();
}

void brain::calc_mf_parameters()
{
	for (int i = 0; i < PoolCount; i++) {
		pools[i].initVar();
	}
}

//==================================================================================================================
// Print Parameters to screen
//==================================================================================================================

void brain::printBrain()
{
	cout << "====== Global Vars ======" << endl;
	cout << "C: " << C << endl;	// connections to pyr cells
	cout << "Cext: " << Cext << endl;	// external connections
	cout << "PoolCount: " << PoolCount << endl;
	cout << "Frequencies: ";
	for (int i = 0; i < PoolCount; i++)
		cout << poolnu[i] << " ";
	cout << endl;

	cout << "Poolsizes: ";
	for (int i = 0; i < PoolCount; i++)
		cout << poolf[i] << " ";
	cout << endl;

	cout << "Extrates: ";
	for (int i = 0; i < PoolCount; i++)
		cout << pools[i].nuext << " ";
	cout << endl;

	cout << "\n====== AMPA Connections =======" << endl;
	printf("  ");
	for (int i = 0; i < PoolCount; i++)
		printf("%4d", i);
	printf("\n");
	for (int i = 0; i < PoolCount; i++) {
		printf("%2d ", i);
		for (int j = 0; j < PoolCount; j++)
			printf("%4.1f", pools[j].w[i].ampa);
		printf("\n");
	}

	cout << "\n====== NMDA Connections =======" << endl;
	printf("  ");
	for (int i = 0; i < PoolCount; i++)
		printf("%4d", i);
	printf("\n");
	for (int i = 0; i < PoolCount; i++) {
		printf("%2d ", i);
		for (int j = 0; j < PoolCount; j++)
			printf("%4.1f", pools[j].w[i].nmda);
		printf("\n");
	}

	cout << "\n====== GABA Connections =======" << endl;
	printf("  ");
	for (int i = 0; i < PoolCount; i++)
		printf("%4d", i);
	printf("\n");
	for (int i = 0; i < PoolCount; i++) {
		printf("%2d ", i);
		for (int j = 0; j < PoolCount; j++)
			printf("%4.1f", pools[j].w[i].gaba);
		printf("\n");
	}
}


void brain::print_pools()
{
	cout.setf(ios::fixed, ios::showpoint);

	cout << endl << setw(15) << "VL: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			VL;
	}
	cout << endl << setw(15) << "Vthr: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			Vthr;
	}
	cout << endl << setw(15) << "Vreset: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			Vreset;
	}
	cout << endl << setw(15) << "VI: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			VI;
	}
	cout << endl << setw(15) << "VE: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			VE;
	}			// constants?!
	cout << endl << setw(15) << "avrV: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			avrV;
	};
	cout << endl << setw(15) << "Cm: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			Cm;
	};			// Membrane capacitance
	cout << endl << setw(15) << "gm: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			gm;
	};			// Membrane leak conductance
	cout << endl << setw(15) << "taurp: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			taurp;
	};			// Refractory period
	cout << endl << setw(15) << "taum: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			taum;
	};			// time constant
	cout << endl << setw(15) << "calpha: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			calpha;
	}
	cout << endl << setw(15) << "cbeta: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			cbeta;
	}
	cout << endl << setw(15) << "cgamma: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			cgamma;
	}

	// synaptic conductances
	cout << endl << setw(15) << "gAMPAexp: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			gAMPAext;
	}
	cout << endl << setw(15) << "gAMPArec: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			gAMPArec;
	}
	cout << endl << setw(15) << "gNMDA: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			gNMDA;
	}
	cout << endl << setw(15) << "gGABA: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			gGABA;
	}

	// gating variables
	cout << endl << setw(15) << "tauAMPA: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			tauAMPA;
	}
	cout << endl << setw(15) << "tauGABA: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			tauGABA;
	}
	cout << endl << setw(15) << "tauNMDAdecay: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			tauNMDAdecay;
	}
	cout << endl << setw(15) << "tauNMDArise: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			tauNMDArise;
	}
	cout << endl << setw(15) << "nuext: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].nuext * 1000;	// overall external input (Hz)
	}

	cout << "\n\n ------- Mean Field Functions --------- ";
	cout << endl << setw(15) << "Phi: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			Phi();
	}
	cout << endl << setw(15) << "falpha: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			falpha();
	}
	cout << endl << setw(15) << "fbeta: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			fbeta();
	}
	cout << endl << setw(15) << "tauE: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			tauE();
	}
	cout << endl << setw(15) << "muE: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			muE();
	}
	cout << endl << setw(15) << "sigmaE: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			sigmaE();
	}
	cout << endl << setw(15) << "Sx: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			Sx();
	}
	cout << endl << setw(15) << "rho1: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			rho1();
	}
	cout << endl << setw(15) << "rho2: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			rho2();
	}
	cout << endl << setw(15) << "nx: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			nx();
	}
	cout << endl << setw(15) << "Nx: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			Nx();
	}
	cout << endl << setw(15) << "J: ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			J();
	}
	// cout << "fak(5): " <<  pools[p].fak(5) << endl;
	cout << endl << setw(15) << "psi(nuE): ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			psi(poolnu[i]);
	}
	cout << endl << setw(15) << "T(2): ";
	for (int i = 0; i < PoolCount; i++) {
		cout << fixed << setprecision(5) << setw(10) << pools[i].
			T(2);
	}
	cout << endl;
}

// bunch of output for testing purposes
void brain::printPool(int p)
{
	cout << "\n====== Local Vars Pool " << p << " =======" << endl;
	cout << "VL: " << pools[p].VL << endl;
	cout << "Vthr: " << pools[p].Vthr << endl;
	cout << "Vreset: " << pools[p].Vreset << endl;
	cout << "VI: " << pools[p].VI << endl;
	cout << "VE: " << pools[p].VE << endl;	// constants?!
	cout << "avrV: " << pools[p].avrV << endl;
	cout << "Cm: " << pools[p].Cm << endl;	// Membrane capacitance
	cout << "gm: " << pools[p].gm << endl;	// Membrane leak conductance
	cout << "taurp: " << pools[p].taurp << endl;	// Refractory period
	cout << "taum: " << pools[p].taum << endl;	// time constant
	cout << "calpha: " << pools[p].calpha << endl;
	cout << "cbeta: " << pools[p].cbeta << endl;
	cout << "cgamma: " << pools[p].cgamma << endl;

	// synaptic conductances
	cout << "gAMPAexp: " << pools[p].gAMPAext << endl;
	cout << "gAMPArec: " << pools[p].gAMPArec << endl;
	cout << "gNMDA: " << pools[p].gNMDA << endl;
	cout << "gGABA: " << pools[p].gGABA << endl;

	// gating variables
	cout << "tauAMPA: " << pools[p].tauAMPA << endl;
	cout << "tauGABA: " << pools[p].tauGABA << endl;
	cout << "tauNMDAdecay: " << pools[p].tauNMDAdecay << endl;
	cout << "tauNMDArise: " << pools[p].tauNMDArise << endl;
	cout << "Nr: " << pools[p].pool_id << endl;
	cout << "nuext: " << pools[p].nuext << endl;	// overall external input (Hz)

	cout << " ------- Functions --------- " << endl;
	cout << "Phi: " << pools[p].Phi() << endl;
	cout << "falpha: " << pools[p].falpha() << endl;
	cout << "fbeta: " << pools[p].fbeta() << endl;
	cout << "tauE: " << pools[p].tauE() << endl;
	cout << "muE: " << pools[p].muE() << endl;
	cout << "sigmaE: " << pools[p].sigmaE() << endl;
	cout << "Sx: " << pools[p].Sx() << endl;
	cout << "rho1: " << pools[p].rho1() << endl;
	cout << "rho2: " << pools[p].rho2() << endl;
	cout << "nx: " << pools[p].nx() << endl;
	cout << "Nx: " << pools[p].Nx() << endl;
	cout << "J: " << pools[p].J() << endl;
	// cout << "fak(5): " <<  pools[p].fak(5) << endl;
	cout << "psi(nuE): " << pools[p].psi(poolnu[p]) << endl;
	cout << "T(2): " << pools[p].T(2) << endl;
}
