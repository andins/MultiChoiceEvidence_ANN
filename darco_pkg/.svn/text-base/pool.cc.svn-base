// ---------with Facilitation---------------------------------------------------------- 
#include <math.h>
#include <iostream>
#include <fstream>

#include "brain.h"

pool::pool()
{
}

pool::~pool()
{
}

void pool::setExternal(double nu)
{
	// Hz
	nuext = nu;
}

void pool::initPool(int i, poolinfo * pi)
{
	pool_id = i;
//	cout << "pool: " << i << endl;
       	num_neurons =
	    static_cast<int> (brain::C * brain::poolf[pool_id] + 0.5);
//	cout << "pool[" << i << "]= " << num_neurons << endl;
	// to make the program more readable
	VL = pi->VL;
	Vthr = pi->Vthr;
	Vreset = pi->Vreset;
	VE = pi->VE;
	VI = pi->VI;
	Cm = pi->Cm;
	gm = pi->gm;
	taurp = pi->taurp;
	ustart = pi->ustart; // MARIO: see setup.cc, pool.h
	taum = pi->taum;
	calpha = pi->calpha;
	cbeta = pi->cbeta;
	cgamma = pi->cgamma;

	gAMPAext = pi->gAMPAext;
	gAMPArec = pi->gAMPArec;
	gNMDA = pi->gNMDA;
	gGABA = pi->gGABA;
	tauAMPA = pi->tauAMPA;
	tauGABA = pi->tauGABA;
	tauNMDAdecay = pi->tauNMDAdecay;
	tauNMDArise = pi->tauNMDArise;

	// precalculations
	// bin = n over k = n!/((n-k)!*k!) needed for Tn (appendix2 Brunel and Wang)
	bin[0][0] = bin[1][0] = bin[1][1] = 1.0;
	fak2[0] = fak2[1] = 1.0;	// fak2 = faculty = n!
	alphatauNMDArisen[0] = 1.0;	// (-1.0* calpha * tauNMDArise)^n
	alphatauNMDArisen[1] = -1.0 * calpha * tauNMDArise;

	for (int n = 2; n <= NMAX; n++) {
		fak2[n] = fak2[n - 1] * (double) (n);
		bin[n][0] = 1.0;
		for (int k = 1; k < n; k++)
			bin[n][k] = bin[n - 1][k - 1] + bin[n - 1][k];
			bin[n][n] = 1.0;
			alphatauNMDArisen[n] = -1.0 * calpha * tauNMDArise *
		    	alphatauNMDArisen[n - 1];
	}

	TEext = (gAMPAext * brain::Cext * tauAMPA) / gm;
	TEAMPA = (gAMPArec * brain::C * tauAMPA) / gm;
	TEI = (gGABA * brain::C * tauGABA) / gm;
	tauNMDA = calpha * tauNMDArise * tauNMDAdecay;

	crho1 = (gNMDA * brain::C) / gm;
	crho2 = cbeta * crho1;
	csigmaE = (pow(gAMPAext, 2) * brain::Cext * pow(tauAMPA, 2)) / pow(gm * taum, 2);

	for (int i = 0; i < max_size; i++) {
		// Precalculated function for values up to 300Hz
		psidata[i] = psi(i * 300e-3 / max_size);
	}
	generateintegraldata();

	// Spiking
	srand(0);

	expdttauAMPA = exp(-brain::dt / tauAMPA);
	expdttauNMDArise = exp(-brain::dt / tauNMDArise);
	expdttauGABA = exp(-brain::dt / tauGABA);

	resetSpikingPool();
}

double pool::get_gAMPArec()
{
	return gAMPArec;
}

double pool::get_gNMDA()
{
	return gNMDA;
}

void pool::set_gAMPArec(double g)
{
	gAMPArec = g;
}

void pool::set_gNMDA(double g)
{
	gNMDA = g;
}

void pool::resetSpikingPool()
{
	for (int i = 0; i < maxpoolneurons; i++){ 
		sAMPAext[i] = 0;
		sAMPArec[i] = 0;
		sNMDA[i] = 0;
		xNMDA[i] = 0;
		sGABA[i] = 0;
		Ca[i] = 0;
		u[i]= ustart; // MARIO: see the code in setup.cc
		V[i] = Vreset + brain::drand() * (Vthr - Vreset);
		dV[i] = 0;
		lastspiketime[i] = taurp + 1;
	}
}
// double pool::brain::drand()
// {
// //   double uniform_rnd = gsl_rng_get(rnd_instance) /
// //           gsl_rng_max(rnd_instance);
// // //        double uniform_rnd = brain::C / (brain::Cext + 1);
// //   return(uniform_rnd);

//      //pseudo brain::drand
//      int i = rand();
//      int j = RAND_MAX;
//      double result = double (i) / double (j);
//      //cout << i<<" "<< j <<" " <<result << endl;
//      return (result);
// }

int pool::extspike()
{
	// nuext in Hz!
	if (brain::drand() < (1e-3 * nuext * brain::Cext * brain::dt))
		return 1;
	else
		return 0;
}

// Improvement potentials
// diverse precalculations, eg. exp(...)!!
// only ampa,nmda OR gaba

void pool::calcSynapses(double time_ms, bool snmda_f, ofstream * snmda_file)
{
	brain::sumAMPArec[pool_id] = 0;
	brain::sumNMDA[pool_id] = 0;
	brain::sumGABA[pool_id] = 0;
	// Take a sample of NMDA activating variables
	// every ´sample_ms' miliseconds.
	int sample_ms = 100;
	
	if (snmda_f == true && snmda_file->is_open() &&
		int(time_ms) % sample_ms == 0){
	    * snmda_file << time_ms << " ";
	}
	for (int i = 0; i < num_neurons; i++) {	// Loop over all neurons in the pool
		// External Spikes can be seen in "advance"
		sAMPAext[i] = sAMPAext[i] * expdttauAMPA;
		if (extspike() == 1) {
			sAMPAext[i] += exp(-brain::dt * brain::drand() / tauAMPA); //why is here brain::drand()? So that spike has fired sometime 
		}							//between the last timestep and now? 
		
		double dsold = (-sNMDA[i] / tauNMDAdecay)+ calpha * xNMDA[i] * (1 - sNMDA[i]);
		double sNMDAeuler = sNMDA[i] + brain::dt * dsold;

		sAMPArec[i] = sAMPArec[i] * expdttauAMPA;
		xNMDA[i] = xNMDA[i] * expdttauNMDArise;
		sGABA[i] = sGABA[i] * expdttauGABA;
		sNMDA[i] = sNMDA[i] + brain::dt * 0.5 * (dsold - (sNMDAeuler / tauNMDAdecay) +  calpha * xNMDA[i] * (1 - sNMDAeuler));
		if (snmda_f == true && snmda_file->is_open() &&
			int(time_ms) % sample_ms == 0){
		    *snmda_file << sNMDA[i] << " "; 
		}
		brain::sumAMPArec[pool_id] += sAMPArec[i]*u[i]; //Facilitation #####################################################################
		brain::sumNMDA[pool_id] += sNMDA[i]*u[i];	//Facilitation
		brain::sumGABA[pool_id] += sGABA[i];
		//if (i == 275 && pool_id == 1 && time_ms < 10.0){
		//	cout << "spikesAMPArec: " << spikesAMPArec[0][i] << "	spikexNMDA: " << spikexNMDA[0][i] << "	spikesGABA: " << spikesGABA[0][i] << "\n";}
	}
	*snmda_file << endl;
}


void pool::calcPotentials(double time, ofstream * file, bool facilON)
{
	//Calculate Incoming Synaptic Potentials
	double inputAMPArec = 0;
	double inputNMDA = 0;
	double inputGABA = 0;
	for (int i = 0; i < brain::PoolCount; i++) {
		inputAMPArec += w[i].ampa * brain::sumAMPArec[i];
		inputNMDA += w[i].nmda * brain::sumNMDA[i];
		inputGABA += w[i].gaba * brain::sumGABA[i];
	}
	sumu=0;		//Facilitation
	for (int i = 0; i < num_neurons; i++) {	// Loop over all neurons
		Ca[i] = Ca[i] - brain::dt * Ca[i] / taoCa;
		if (facilON == 1){
		u[i]=u[i]+brain::dt*(U-u[i])/taoF; // ####################################################
                sumu=sumu+u[i];
		}
		if (lastspiketime[i] > taurp) {
			//Heun Algorithm for new potential
			//Ca[i] = Ca[i] - brain::dt * Ca[i] / taoCa; Comment for Facilitation
			double Veuler = V[i] + brain::dt * dV[i];
			V[i] = 	V[i] + brain::dt * 0.5 * 		// 0.5 from Heun algorithm
				(dV[i] - (gm / Cm) * (Veuler - VL) 
				- (gAMPAext * (Veuler - VE) * sAMPAext[i] 
				+  gAMPArec * (Veuler - VE) * inputAMPArec 
				+  gNMDA * (Veuler - VE) * inputNMDA / (1. + cgamma * exp(-cbeta*Veuler))
				+  gGABA * (Veuler - VI) * inputGABA 
				+  AHP_Ca * gAHP * Ca[i] * (Veuler - VK)) / Cm);
			//V[i] = Veuler;
			//Derivative at this point - used in the next step...
			dV[i] = - (gm / Cm) * (V[i] - VL)
			    	- (gAMPAext * (V[i] - VE) * sAMPAext[i]
			       	+  gAMPArec * (V[i] - VE) * inputAMPArec
			       	+  gNMDA * (V[i] - VE) * inputNMDA / (1. + cgamma * exp(-cbeta*V[i]))		
				+  gGABA * (V[i] - VI) * inputGABA 
				+  AHP_Ca * gAHP * Ca[i] * (V[i] - VK)) / Cm;
			// Spike Generation
			if (V[i] > Vthr) {
// General raster
//                                 if ( i < 10 && file->is_open() ){
//                                         *file << time << "\t"
//                                                       << pool_id * 10 + i << "\n";
//                                 }
				// Typical electrophysiological data: measures of one single
				// neuron across several trials.
				// Select neuron #2 in the first and second selective pool.
//<				if ((pool_id == 1 || pool_id == 2)
//<				    && i == 2 && file->is_open()) {
//<					*file << time << "\t" << pool_id <<
//<					    "\n";
//<				}
				if ( file->is_open() && 
                                      (   (pool_id == 0) && (i < 20) 
                                       || (pool_id > 0) && (pool_id < brain::PoolCount-1) && (i < 10)
                                       || (pool_id == brain::PoolCount-1) && (i < 60) ) ) {
					*file << time << "\t" << pool_id 
                                            << "\t" << i << "\n";
				}
				brain::spikes[pool_id]++;
				lastspiketime[i] = 0;
				sAMPArec[i] += 1.;
				xNMDA[i] += 1.;
				sGABA[i] += 1.;
				V[i] = Vreset;
				Ca[i] += alphaCa;
				if (facilON == 1){
					u[i] += U*(1-u[i]); //Facilitation ·###########################################
					//u[i] = 1;
				}
			}
		}
		lastspiketime[i] += brain::dt;
	}
	brain::sumup[pool_id]=brain::sumup[pool_id]+sumu; //Facilitation
}

//==================================================================================================================
// meanfield section
//==================================================================================================================

void pool::initVar()
{
	//init variables - order is important!
	vnx = nx();
	vNx = Nx();
	vnix = nix();
	vSx = Sx();
	vtauE = tauE();
	vmuE = muE();
	vsigmaE = sigmaE();
}

double pool::Euler_derivative(double eulerdelta)
{
	double phi2;

	// don't forget to init variables before calling this function
	initVar();

	phi2 = Phi2(); //Phi2: return (1e3 / (taurp + vtauE * sum));
	//cout << "Phi2: " << phi2 << "	";
	//cout << "poolnu: " << brain::poolnu[pool_id] << "	";
	//cout << "vtauE: " << vtauE << "		";
	if(phi2 == -1)
	    brain::time_forward=true;
	if(brain::time_forward)
	    // Returns in Hz!
	    return (eulerdelta * (-1 * brain::poolnu[pool_id] + phi2)
		    / vtauE);
	else
	    // Time is reversed (useful to find unstable points)
	    return (eulerdelta * (brain::poolnu[pool_id] - phi2)
		    / vtauE);

}

void pool::newavrV()
{
	avrV =
	    vmuE - (Vthr -
		    Vreset) * 1.0e-3 * brain::poolnu[pool_id] * vtauE;
// <-relaxation, nearly factor 2 in speed!
}

//======================= 

double pool::Flow()
{

	// don't forget to init variables before calling this function
	initVar();

	return (Phi2());
}

double pool::Phi2()
{
	// Should return rate values in Hz!

	double b = fbeta();	//start point
	double a = falpha();	//end point

	const double min_val = -20.0;
	const double max_val = +10.0;
	const double step = (max_val - min_val) / max_size;

	if ((b < min_val) || (a > max_val)) {
		cout << a << "\t" << b << "\t" << 
		    "Phi2: Integral Range exceeds precalculated values." <<
		    endl;
		//return Phi(); // Use the slow function instead of exiting
		return -1;
	}

	int xi = static_cast<int> ((b - min_val) / step);	// position in array
	int yi = static_cast<int> ((a - min_val) / step);	// position in array

	double xf = ((b - min_val) - xi * step) / step;	// fraction between 0,1
	double yf = ((a - min_val) - yi * step) / step;	// fraction between 0,1

	double sum = 0.0;
	int start, end;

	if (xf > 0.5) {
		sum += intdata[xi + 1] * (1.5 - xf);
		start = xi + 2;
	} else {
		sum += intdata[xi] * (0.5 - xf);
		start = xi + 1;
	}

	if (yf > 0.5) {
		sum += intdata[yi + 1] * (yf - 0.5);
		end = yi;
	} else {
		sum += intdata[yi] * (yf + 0.5);
		end = yi - 1;
	}

	int broad_step = max_size / max_broadsize;

	int inter1 = start + (broad_step - (start % broad_step));
	int inter2 = end - (end % broad_step);

	int broad1 = static_cast<int> (start / broad_step) + 1;
	int broad2 = static_cast<int> (end / broad_step);

	for (int i = start; i < inter1; i++)
		sum += intdata[i];

	for (int i = broad1; i < broad2; i++)
		sum += intdatabroad[i];

	for (int i = inter2; i <= end; i++)
		sum += intdata[i];


	sum *= sqrt(PI) * 0.002;
	return (1e3 / (taurp + vtauE * sum));
}


//old (slow) phi function
double pool::Phi()
{
	double a = falpha();
	double b = fbeta();

	double z;
	int N = 1000;
	double sum = 0;

	for (int i = 0; i <= N; i++) {
		z = b + (a - b) * i / N;
		if (i == 0 | i == N)
			sum += 0.5 * nerf(z);
		else
			sum += nerf(z);
	}

	sum *= sqrt(PI) * fabs(a - b) / N; // last term is because of numerical integration

	return (1e3 / (taurp + vtauE * sum));
}


//=======================
double pool::falpha()
{
	return (Vthr - vmuE) / vsigmaE * (1.0 + 0.5 * tauAMPA / vtauE) +
	    1.03 * sqrt(tauAMPA / vtauE) - 0.5 * tauAMPA / vtauE;
}

double pool::fbeta()
{
	return (Vreset - vmuE) / vsigmaE;
}

double pool::nerf(double z)
{				/* function exp(z^2)(erf(z)+1) */
	double t, ef, at;
	double w;
	w = fabs(z);
	t = 1. / (1. + 0.5 * w);
	at = a1 + t * (a2 +
		       t * (a3 +
			    t * (a4 +
				 t * (a5 +
				      t * (a6 +
					   t * (a7 +
						t * (a8 +
						     t * (a9 +
							  t * a10))))))));
	ef = t * exp(at);
	if (z > 0.)
		ef = 2. * exp(w * w) - ef;
	return (ef);
}

void pool::generateintegraldata() // necessary for intdata precalculation needed in Phi2()
{
	int c = 0;
	int d = 0;
	for (int i = 0; i < max_size; i++) {
		intdata[i] = nerf(-20.0 + i * 30.0 / max_size);
		if (c < max_size / max_broadsize)
			intdatabroad[d] += intdata[i];
		else {
			c = 0;
			d++;
			intdatabroad[d] = intdata[i];
		}
		c++;
	}
}

//======================= necessary for eq(8) = diffeq of V
double pool::muE()
{
	double res;
	if (VE == 0)
		res = (VL + rho2() * vNx * avrV + TEI * vnix * VI) / vSx;
	else
		res = (VL + (TEext * 1e-3 * nuext + TEAMPA * vnx + rho1() * vNx) * VE 
			  + rho2() * vNx * avrV + TEI * vnix * VI) / vSx;
	return res;
}

double pool::tauE()
{
	return Cm / (gm * vSx);
}

double pool::sigmaE()
{
	//sqrt((pow(gAMPAext,2)*
	// pow(avrV-VE,2)*Cext*nuext*pow(tauAMPA,2)*tauE())/pow(gm*taum,2));
	return sqrt(pow(avrV - VE, 2) * vtauE * csigmaE * nuext * 1e-3);
}

double pool::Sx()
{
	return (1.0 + TEext * nuext * 1e-3 + TEAMPA * vnx +
		(rho1() + rho2()) * vNx + TEI * vnix);
}

//======================= necessary for some of the above equations
double pool::J()
{
	return 1 + cgamma * exp(-cbeta * avrV);
}

double pool::rho1()
{
	//(gNMDA * CE)/(gm * J());
	return (crho1 / J());
}

double pool::rho2()
{
	//cbeta * (gNMDA * CE * (avrV-VE)* (J()-1)) / (gm * pow(J(),2));
	return crho2 * (avrV - VE) * (J() - 1) / pow(J(), 2);
}

//======================= (initial values for nuE=nx, psi(nuE)=Nx and nuI=nix)

double pool::nx()
{
	double res = 0;
	for (int i = 0; i < brain::PoolCount; i++) {
		res +=
		    brain::poolf[i] * w[i].ampa * 1.0e-3 *
		    brain::poolnu[i];
// incoming spikes from all pools
	}
	return res;
}

double pool::Nx()
{
	double res = 0;
	for (int i = 0; i < brain::PoolCount; i++) {
		res +=
		    brain::poolf[i] * w[i].nmda * psi2(1.0e-3 *
						       brain::poolnu[i]);
	}
	return res;
}

double pool::nix()
{
	double res = 0;
	for (int i = 0; i < brain::PoolCount; i++) {
		if (w[i].gaba != 0)	// taking very few inhib. pools into account
			res += brain::poolf[i] * w[i].gaba *
			    1.0e-3 * brain::poolnu[i];
// incoming spikes from all pools
	}
	return res;
}

//======================= necessary calculations for nonlinear NMDA Nx()
double pool::psi2(double nu)
{
	// Precalculation for nu < 300Hz
	double steps = max_size / 300e-3; // 50000
	// nu here is in kHz!
	int a = static_cast<int> (nu * steps);
	if (nu < 300e-3) // Use precalculated values
		return (psidata[a] + (psidata[a + 1] - psidata[a])
			* (nu - (a / steps)) * steps);
	else // Calculate it directly
		return psi(nu);
}

double pool::psi(double nu)
{
	// nu in kHz!
	double sum1 = 0;
	onenutauNMDA = 1 + nu * tauNMDA;
	for (int n = 1; n < 7; n++) {
		//sum1 += pow(-1* calpha * tauNMDArise,n) * T(n,nu) / fak2[n+1];
		sum1 += alphatauNMDArisen[n] * T(n) / fak2[n + 1];
	}
	return (nu * tauNMDA) / onenutauNMDA * (1 + (1 / onenutauNMDA) * sum1);
}

double pool::T(int n)
{
	// depends on psi for onenutauNMDA and therefore on it's nu!!
	double sum = 0;
	for (int k = 0; k <= n; k++) {
		if (k % 2 == 0)
			sum += bin[n][k] * (tauNMDArise * onenutauNMDA) /
			       (tauNMDArise * onenutauNMDA + k * tauNMDAdecay);
		else
			sum -= bin[n][k] * (tauNMDArise * onenutauNMDA) /
			       (tauNMDArise * onenutauNMDA + k * tauNMDAdecay);
	}
	return sum;
}
