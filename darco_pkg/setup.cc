#include "brain.h"
#include "setup.h"
#include <iostream> // MARIO
#include <fstream> // MARIO
#include <string>// MARIO
#include <iomanip>// MARIO

brain pfc;
struct poolinfo pool_inh, pool_exc;

void set_pool_data()
{
	//
	// Brunel & Wang.
	// Conductances fixed for a network with 1000 neurons.
	// These are automatically rescaled according to the
	// number of neurons.
	//

	// -------------------------
	//    Inhibitory pool
	// -------------------------

	// MARIO loading param. from an external file:
	string line;
	ifstream neuronsfile ("neurons.ini");
        if (neuronsfile.is_open())	{
		// ***  Neurons ***
		getline(neuronsfile,line);
		neuronsfile >> pool_inh.VL;	// resting potential (mV)
		neuronsfile >> pool_inh.Vthr;	// firing threshold (mV)
		neuronsfile >> pool_inh.Vreset;	// reset potential (mV)
		neuronsfile >> pool_inh.VI;
		neuronsfile >> pool_inh.VE;
	 	neuronsfile>>pool_inh.Cm;	// membrane capacitance (pF)
		neuronsfile>>pool_inh.gm;	// leak conductance (nS)
		neuronsfile>>pool_inh.taurp;	// refractory period (ms)
		neuronsfile>>pool_inh.taum;	// membrane time constant (ms)
		neuronsfile>>pool_inh.calpha;
		neuronsfile>>pool_inh.cbeta;
		neuronsfile>>pool_inh.cgamma;

		// *** Synaptic conductances ***
		neuronsfile>>pool_inh.gAMPAext;
		neuronsfile>>pool_inh.gAMPArec;
		neuronsfile>>pool_inh.gNMDA;
		neuronsfile>>pool_inh.gGABA; 

		// --- gating variables ---
		neuronsfile>>pool_inh.tauAMPA;
		neuronsfile>>pool_inh.tauGABA;
		neuronsfile>>pool_inh.tauNMDAdecay;
		neuronsfile>>pool_inh.tauNMDArise;


		// --- facilitation variables ---
		neuronsfile>>pool_inh.ustart;

		// Excitatory populations:
		// ***  Neurons ***
		neuronsfile >> pool_exc.VL;	// resting potential (mV)
		neuronsfile >> pool_exc.Vthr;	// firing threshold (mV)
		neuronsfile >> pool_exc.Vreset;	// reset potential (mV)
		neuronsfile >> pool_exc.VI;
		neuronsfile >> pool_exc.VE;
		neuronsfile>>pool_exc.Cm;	// membrane capacitance (pF)
		neuronsfile>>pool_exc.gm;	// leak conductance (nS)
		neuronsfile>>pool_exc.taurp;	// refractory period (ms)
		neuronsfile>>pool_exc.taum;	// membrane time constant (ms)
		neuronsfile>>pool_exc.calpha;
		neuronsfile>>pool_exc.cbeta;
		neuronsfile>>pool_exc.cgamma;

		// *** Synaptic conductances ***
		neuronsfile>>pool_exc.gAMPAext;
		neuronsfile>>pool_exc.gAMPArec;
		neuronsfile>>pool_exc.gNMDA;
		neuronsfile>>pool_exc.gGABA; 

		// --- gating variables ---
		neuronsfile>>pool_exc.tauAMPA;
		neuronsfile>>pool_exc.tauGABA;
		neuronsfile>>pool_exc.tauNMDAdecay;
		neuronsfile>>pool_exc.tauNMDArise;

		
		// --- facilitation variables ---
		neuronsfile>>pool_exc.ustart;

		neuronsfile.close();

	}
	// MARIO if the external file doesn't exist use the default values:
	else {
		cout << "Unable to open file neurons.ini\n";

		// ***  Neurons ***
		pool_inh.VL = -70.0;	// resting potential (mV)
		pool_inh.Vthr = -50.0;	// firing threshold (mV)
		pool_inh.Vreset = -55.0;	// reset potential (mV)
		pool_inh.VI = -70.0;
		pool_inh.VE = 0.0;
		
		pool_inh.Cm = 200.0;	// membrane capacitance (pF)
		pool_inh.gm = 20.0;	// leak conductance (nS)
		pool_inh.taurp = 1.0;	// refractory period (ms)
		pool_inh.taum = 10.0;	// membrane time constant (ms)
		pool_inh.calpha = 0.5;
		pool_inh.cbeta = 0.062;
		pool_inh.cgamma = 0.2801120448;

		// *** Synaptic conductances ***
		pool_inh.gAMPAext = 1.62;
		pool_inh.gAMPArec = 0.081;
		pool_inh.gNMDA = 0.258;
		pool_inh.gGABA = 2.0; 

		// --- gating variables ---
		pool_inh.tauAMPA = 2.0;
		pool_inh.tauGABA = 5.0;
		pool_inh.tauNMDAdecay = 100.0;
		pool_inh.tauNMDArise = 2.0;

		// -------------------------
		//    Excitatory pools
		// -------------------------

		// *** Neurons ***
		pool_exc.VL = -70.0;
		pool_exc.Vthr = -50.0;
		pool_exc.Vreset = -55.0;
		pool_exc.VI = -70.0;
		pool_exc.VE = 0.0;

		pool_exc.Cm = 500.0; // pF
		pool_exc.gm = 25.0; // nS
		pool_exc.taurp = 2.0; // ms
		pool_exc.taum = 20.0; // ms
		pool_exc.calpha = 0.5; // ms^{-1}
		pool_exc.cbeta = 0.062;
		pool_exc.cgamma = 0.2801120448;

		// *** Synaptic conductances ***
		pool_exc.gAMPAext = 2.08;
		pool_exc.gAMPArec = 0.104 ;
		pool_exc.gNMDA = 0.327 ;
		pool_exc.gGABA = 2.0*1.3; 

		// --- gating variables ---
		pool_exc.tauAMPA = 2.0;
		pool_exc.tauGABA = 5.0;
		pool_exc.tauNMDAdecay = 100.0;
		pool_exc.tauNMDArise = 2.0;

		// MARIO:
		// --- facilitation variables ---
		pool_exc.ustart = 0.3;
		pool_inh.ustart = 0.6;
	}//MARIO
}

void setup(int num_neurons, int num_pools, int num_ext_neurns,
	   double time_step, double inh_percent,
	   double coding_lvl, double g_ratio, double ratio)
{
	set_pool_data();
	// Scaling ratio
	if (ratio!=true){
		pool_exc.gAMPArec = pool_exc.gAMPArec *ratio;
		pool_exc.gNMDA = pool_exc.gNMDA*ratio;
		pool_exc.gGABA = pool_exc.gGABA*ratio; 
		pool_inh.gAMPArec = pool_inh.gAMPArec*ratio;
		pool_inh.gNMDA = pool_inh.gNMDA*ratio;
		pool_inh.gGABA = pool_inh.gGABA*ratio; 
	}
	pfc.initBrain(num_pools, num_neurons, num_ext_neurns, time_step);

//	cout << "num pools: " << num_pools << "; num neur: " << num_neurons << endl; // DEBUG MARIO

        // If we want to regulate the NMDA contribution (using -N option, in a
	// range from 0 to 1). 0 means no modification between the parameters
	// given by B&W. 1 means entire AMPA contribution, by keeping the
	// average external excitatory inputs equal to the average recurrent
	// excitatory inputs.
	if (g_ratio != 0.0) {
		cerr << "NMDA/AMPA ratio has been modified...\n";
		// 10, the ratio of the NMDA over the AMPA component
		// in terms of the charge entry near the threshold, -55mv
		// (See Brunel & Wang, p.66)
		pool_inh.gNMDA *= (1 - g_ratio);
		pool_exc.gNMDA *= (1 - g_ratio);
		pool_inh.gAMPArec *= (1 + 10 * g_ratio);
		pool_exc.gAMPArec *= (1 + 10 * g_ratio);
	}
	// MARIO loading param. from an external file:
	string line;
	ifstream pfcfile ("pfc.ini");
	if (coding_lvl==2){ // or a condition to discriminate between the condition with one layer
		// and the condition with, for example, multipe inhibitory population.
	        if (pfcfile.is_open())	{
			getline(pfcfile,line);
			char linePfc[200];
			while(pfcfile.getline(linePfc,200)) {
				int Cinput, pid, num_neurons_p, exc_flag;
				Cinput=sscanf(linePfc, "%d %d %d",&pid, &num_neurons_p,&exc_flag);
				float perc_neur=(float) num_neurons_p/num_neurons;
//				cout << "P[" << pid << "], perc neur: " << perc_neur << endl; // DEBUG MARIO
				if (exc_flag==1)
					pfc.initPool(pid,perc_neur, 3.0, -52.5, 3.0, &pool_exc);
				else
					pfc.initPool(pid,perc_neur, 9.0, -52.5, 3.0, &pool_inh);
			}
			pfcfile.close();
		}
		else{   // if pfc.ini doesn't exist something went wrog
	//		disp(ERROR)
			cout << "No pfc.ini file. Check the coherence of your inin files.\n";
			exit(1);
		}
	
	}
	else
	{
		// initialize the inhibitory pool
		//cout << "perc neur" << inh_percent << endl; // DEBUG MARIO
		pfc.initPool(0, inh_percent, 9.0, -52.5, 3.0, &pool_inh);
		// renormalize coding_lvl to the whole network (inh. + exc.)
		double global_coding_lvl = coding_lvl * (1.0 - inh_percent);
		//cout << "perc neur selec" << global_coding_lvl << endl; // DEBUG MARIO
		// initialize the excitatory selective pools
		for (int i = 1; i < num_pools - 1; i++) {
			pfc.initPool(i, global_coding_lvl, 3.0, -52.5, 3.0,&pool_exc);
		}
		// initialitze the excitatory non-selective pool
		float cod2 = (1 - inh_percent) * (1 - coding_lvl * (num_pools - 2));
		//cout << "perc neur non-selec" << cod2 << endl; // DEBUG MARIO
		pfc.initPool(num_pools - 1, cod2, 3.0, -52.5, 3.0, &pool_exc);
	}
}

void set_all_connections(int num_pools, double coding_lvl, double w_plus, double w_I, double w_T, string connectivityFileName)
{	//w_T has no use here for 2 choices
	// ----------------------------------------
	// setting the connections
	// ----------------------------------------
	// MARIO loading param. from an external file:
	string line;
	ifstream connfile (connectivityFileName.c_str());
        if (connfile.is_open())	{
		// ***  connections ***
		// set all the connectivity to zero!
		for (int i = 0; i < num_pools; i++) {
	 		for (int j = 0; j < num_pools; j++) {
				pfc.setConnection(j, i, 0, 0, 0);
			}
		}

		getline(connfile,line);
		char lineCon[200];
		while(connfile.getline(lineCon,200)) {
//			cout << "[ " << lineCon << " ]" << endl;
			int pPre, pPost, Cinput;
			float wCon,wAmpaCon,wGabaCon;
			Cinput=sscanf(lineCon, "%d %d %f %f %f",&pPost, &pPre, &wCon,&wAmpaCon,&wGabaCon);
			if (Cinput!=5) {
				Cinput=sscanf(lineCon, "%d %d %f %f",&pPost, &pPre, &wCon,&wAmpaCon);
				if (Cinput!=4) {
					Cinput=sscanf(lineCon, "%d %d %f",&pPost, &pPre, &wCon);
					if (Cinput!=3){
						cout << "error" << endl;
					}
					else{
					       	if (wCon>0){
							pfc.setConnection(pPre, pPost, wCon, wCon,0); // 3 infos
						}
						else{
							pfc.setConnection(pPre, pPost, 0, 0,-wCon); // 3 infos 
						}
						//cout << "pre:" << pPre << " post:" << pPost << ' ' << wCon << endl;
					}
				}
				else{
					pfc.setConnection(pPre, pPost, wCon, wAmpaCon,0); // 4 infos
				}
			}
			else
				pfc.setConnection(pPre, pPost, wCon, wAmpaCon,wGabaCon); // 5 infos
		}
		connfile.close();
	}
	// MARIO if the external file doesn't exist use the default values:
	else {
		cout << "Unable to open file connectivity.ini\n";
		double w_minus = (1. - coding_lvl * (2.*w_T + w_plus)) / (1. - coding_lvl);

		// inhibitory projections: baseline
		for (int i = 0; i < num_pools; i++) {
			pfc.setConnection(0, i, 0, 0, w_I);
		}

		// excitatory projections
		for (int i = 1; i < num_pools; i++) {
			// excitatory to inhibitory: baseline
			pfc.setConnection(i, 0, 1.0, 1.0, 0);

			// excitatory to non-selective pool: baseline
			// TODO: Shouldn't it be wminus?
			pfc.setConnection(i, num_pools - 1, 1.0, 1.0, 0);

			// excitatory to selective pools: depressed
			for (int j = 1; j < num_pools - 1; j++)
				pfc.setConnection(i, j, w_minus, w_minus, 0);
				//pfc.setConnection(i, j, 0.8, 0.8, 0);
		}

		// recurrent potentiation for selective pools
		// (overwriting previous assignments)
		for (int i = 1; i < num_pools - 1; i++) {
			pfc.setConnection(i, i, w_plus, w_plus, 0);
		}

	} // MARIO
	
}

