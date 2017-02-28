/*=============================================================================
 *
 *  spiking: a leaky Integrate & Fire network simulator, using Brunel & Wang
 *	     model.
 *_____________________________________________________________________________
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>		// for atoi() and atof()
#include <string>
#include <unistd.h>		// for getopt()
#include <stdio.h>		// for exp
#include <math.h>		// for exp

#include "brain.h"
#include "setup.h"

void help_message_spkg(char *program)
{
	cerr << "usage: " << program
	    << "\n\t [-h (for help)]\n"
	    << "\t [-c connectivity_ratio]\n"  //MARIO set a ratio to proportional change all the w/j
	    << "\t [-F facilitation on]\n"  //ANDREA 
	    << "\t [-n number_neurons]\n"
	    << "\t [-j potentiated_synaptic_wPlus]\n"
	    << "\t [-W inhibitory connections]\n"
	    << "\t [-T added tuned connections]\n"
	    << "\t [-o output_file]\n"
	    << "\t [-t number_of_trials]\n"
	    << "\t [-r random_seed]\n"
	    << "\t [-w window_time_size_in_ms]\n"
	    << "\t [-s output_spikes_file]\n"
	    << "\t [-e number_external_neurons]\n"
	    << "\t [-a number_of_averages]\n"
	    << "\t [-C connectivityFileName]\n"
	    << "\t [-p protocolFileName]\n"
	    << "\t [-S output_filename_for_Snmda]\n"
	    << "\t [-D should the trial finish when a decision is made?]\n"
	    << "\t [-N g_AMPA_factor_and_suppress_NMDA]\n" 
	    << "\t [-L trial total simulation time]" << endl;
}

/*-----------------------------------------------------------------------------
 *  read_options  --  reads the command line options, and implements them.
 *
 *  note: when the help option -h is invoked, the return value is set to false
 *        to prevent further execution of the main program; likewise, if an
 *        unknown option is used, the return value is set to false.
 *--------------------------------------------------------------------------- */
bool read_options_spkg(int argc, char *argv[], int &num_neurons,
		       double &w_I, double &w_T, unsigned long int
		       &seed, int &trials, int &num_ext_neurons, string &
		       output_file, int &window_size, double
		       &wPlus, string & output_spikes, double & g_ratio, string & snmda_filename, int & trial_life,
		       bool & b_self_finish, bool & facilON, bool & verboseON,string & protocolFileName,
		       double & ratio, string & connectivityFileName) {
	int c;
	// Colons in the third getopt argument allow args for the options
	while ((c = getopt(argc, argv, "hDFvn:l:W:T:r:t:e:o:w:j:s:z:S:N:L:p:c:C:?")) != -1)
		switch (c) {
		case 'h':
			help_message_spkg(argv[0]);
			return false;	// execution should stop after help
		case 'D':
			b_self_finish = true;
			break;
		case 'F':
			facilON = true;
			break;
		case 'n':
			num_neurons = atoi(optarg);
			break;
		case 'p':
			protocolFileName = optarg;
			break;
		case 'C':
			connectivityFileName = optarg;
			break;
		case 'W':
			w_I = atof(optarg);
			break;
		case 'T':
			w_T = atof(optarg);
			break;
		case 'r':
			seed = atol(optarg);
			break;
		case 't':
			trials = atoi(optarg);
			break;
		case 'e':
			num_ext_neurons = atoi(optarg);
			break;
		case 'o':
			output_file = optarg;
			break;
		case 'w':
			window_size = atoi(optarg);
			break;
		case 'j':
			wPlus = atof(optarg);
			break;
		case 's':
			output_spikes = optarg;
			break;
		case 'S':
			snmda_filename = optarg;
			break;
		case 'N':
			g_ratio = atof(optarg);
			break;
		case 'L':
			trial_life = atoi(optarg);
			break;
		case 'v':
			verboseON = true;
			break;
		case 'c':
			ratio = atof(optarg);
			break;
		case '?':
			help_message_spkg(argv[0]);
			return false;	// execution should stop after error
		}

	return true;		// ready to continue program execution
}

inline std::string stringify_extension(int i)
{
	// Not very neat nor simple:
	ostringstream temp;
	temp.width(4);
	temp.fill('0');
	temp << i;
	return temp.str();
}

int main(int argc, char *argv[])
{
	// ---------------------------
	//   Initialize pfc
	// ---------------------------
	double time_step = 0.02;	// in ms
	double inh_percent;	// inhibitory population (in %).

	int num_pools=0; 		//decision making 2 alternatives
	// Initialize arguments to appropriate default values
	// These values may be changed as command-line options
	int num_neurons=0;	// option -n 1000
	int num_ext_neurons=800;	// option -e 800
      	double coding_lvl;

	bool facilON = false;	//ANDREA facilitation switch: if set to 0 remember to set ustart to 1 in neurons.ini
	bool verboseON = false; // MARIO output off at the command line if is not -v 
	string protocolFileName = "protocol.ini"; // MARIO protocolFileName can change at the command line
	string connectivityFileName = "connectivity.ini"; // MARIO connectivityFileName can change at the command line

	double ratio = 1;	// MARIO the default value for the ratio 
	string output_file = "trial";	// option -o data.dat
	string output_uout = "ufacil";	// option -o uout.dat
	string output_spikes = "";	// option -s spikes.dat
	string snmda_filename = "";	// option -S snmda.dat
	
	int window_size = 50;	// option -w 50
	int window_step = 5;
	double w_T = 0.0;		// additional tuned connections	
	double g_ratio = 0.0; // option -N 4.0			// (if > 0.5 NMDA is switched off)
	double tau2 = 15.; // 15.0; //for second exponential 
	double tot_ext_rate = 2400.0; 

	unsigned long int seed = 0;
	bool b_self_finish = false;

	char * descr = getenv("GSL_RNG_SEED");
	if(descr){
	    seed = strtoul(descr,0,0);
	}
	
	
	// ---------------------------
	//   Variables
	// -----------------------
	int trial_life = 1000; //  trial total simulation time in milliseconds option -L 5000 
	double w_I= 0.97; // W_i inhibicio
	double wPlus = 2.17;	// option -j 1.7
	int trials =1;		// option -t 1
	
	// MARIO loading param. from an external file:
	string line;
	ifstream pfcfile ("pfc.ini");
        if (pfcfile.is_open())	{
		getline(pfcfile,line);
		char linePfc[200];
		int pid, Cinput, num_neurons_p, exc_flag, num_inh_neurons=0,TEMPnum_neurons,TEMPnum_pools;
		while(pfcfile.getline(linePfc,200)) {
/*			pch = strtok (linePfc," ,-");// DEBUG MARIO
		      	while (pch != NULL)
	    		{
	//			printf ("%s\n",pch);
				pch = strtok (NULL, " ,-");
				cout << "pch" << pch << endl;
			}*/

			Cinput=sscanf(linePfc, "%d %lf %d %lf",&TEMPnum_neurons, &inh_percent, &TEMPnum_pools, &coding_lvl);
//			cout << "Cinput: " << Cinput << endl; MARIO DEBUG
			if (Cinput!=4)
		       	{
				Cinput=sscanf(linePfc, "%d %d %d",&pid, &num_neurons_p, &exc_flag);
				if (Cinput!=3) {
					cout << "Bad definition of pool.ini structrure.\n";
					exit(1);
				}
				else{
					num_pools += 1;
					num_neurons += num_neurons_p;
					if (exc_flag==false)
						num_inh_neurons+=num_neurons_p;	
				}
			}
			else {
				num_pools = TEMPnum_pools;
				num_neurons = TEMPnum_neurons;
				break;
			}
		}
		if (num_inh_neurons>0){
			coding_lvl = 2;
			inh_percent = (float) num_inh_neurons /num_neurons;
//			cout << "#pools= " << num_pools << "; inh_percent " << inh_percent << endl; // DEBUG MARIO
		}
	}
	// MARIO if the external file doesn't exist use the default values:
	else {
		cout << "Unable to open file pfc.ini\n";
        	inh_percent = 0.20;	// inhibitory population (in %).
		num_pools = 4; 		//decision making 2 alternatives
		// Initialize arguments to appropriate default values
		// These values may be changed as command-line options
		num_neurons = 1000;	// option -n 1000
		// num_pools_with_input = 2; //option -p4	
		coding_lvl = 0.1;
//		bool LayerStruct = false;
		
	} // MARIO

	cout << "coding_lvl: " << coding_lvl  << endl;
	
   	ofstream fs("parametros.txt");    
       // comments in the output file parametros
		fs << endl
		   <<  trials <<  endl
		    << wPlus <<  endl
		    << w_I <<  endl
		    << num_neurons << endl
		    <<  coding_lvl <<  endl;
   		fs.close();
   
        tot_ext_rate = tot_ext_rate / (double) num_ext_neurons;
	if (!read_options_spkg
	    (argc, argv, num_neurons, 
	     w_I, w_T, seed, trials,
	     num_ext_neurons, output_file, window_size, 
	     wPlus, output_spikes, g_ratio,
	     snmda_filename, trial_life, b_self_finish, facilON, verboseON, protocolFileName,ratio, connectivityFileName))
		return 1;	// halt criterion detected by read_options()

	// TODO: Check the input
	if (g_ratio < 0 || g_ratio > 1){
		cerr << "The argument for the -N option must be in the range [0,1].\n";
		exit(1);
	}	
//	cout << "\nNumber of pools with input: " << num_pools_with_input << endl;		
	setup(num_neurons, num_pools, num_ext_neurons, time_step,inh_percent, coding_lvl, g_ratio, ratio);	
	
	if (facilON==false){
		cout << "no facilitation" << endl;
		pool_exc.ustart	= 1;
		pool_inh.ustart	= 1;
	}
	
	// --------------------------------------------------
	// Set the connections
	// --------------------------------------------------
	set_all_connections(num_pools, coding_lvl, wPlus, w_I, w_T,connectivityFileName);	
	
	
	
	// ------------------------------------------------------------
	// Start the simulation
	// ------------------------------------------------------------
	pfc.set_rseed(seed);
	cout << "\nRandom seed: " << seed << endl;
	
	pfc.set_window_size(window_size);
	pfc.set_window_step(window_step);
	
	ofstream data, uout, spdata, snmda1, snmda2, snmda3;

	// *********************************
	// *** Spiking Neuron simulation ***
	// *********************************
	for (int k = 0; k < trials; k++) {
		//to change b during trials
		cout << "k_value= " << k << endl; 

		string output_file_ext = output_file + "_" +    stringify_extension(k + 1);
		data.open(output_file_ext.c_str(), ios::out);
		
		string output_uout_ext = output_uout + "_" +   stringify_extension(k + 1);
		if (facilON == true){
		uout.open(output_uout_ext.c_str(), ios::out);
		}

		string output_spikes_ext = output_spikes + "_" +  stringify_extension(k + 1);
		if (output_spikes != "") {
			spdata.open(output_spikes_ext.c_str(), ios::out);
		}
		if (snmda_filename != "") {
		    // For NMDA component stuff
		    string snmda_filename_ext = snmda_filename + "_" +stringify_extension(k + 1);
		    string snmda_filename_ext_1 = snmda_filename_ext +"_1";
		    string snmda_filename_ext_2 = snmda_filename_ext +"_2";
		    string snmda_filename_ext_3 = snmda_filename_ext +"_3";
		    snmda1.open(snmda_filename_ext_1.c_str(), ios::out);
		    snmda2.open(snmda_filename_ext_2.c_str(), ios::out);
		    snmda3.open(snmda_filename_ext_3.c_str(), ios::out);
		}		    
		cout << "Simulating trial " << k + 1 << "...";
		cout.flush();

		// comments in the output file
		data << "%\n% num_neurons: " << num_neurons
		    << ", coding level: " << coding_lvl
  		    << ", extrate: " << tot_ext_rate
		    << ", num pools: " << num_pools
		    << "% w_{+}: " << wPlus	
		    << ", w inhibitory: " << w_I
		    << ", w_T tuned: " << w_T 
		    << ", AN: " << g_ratio	
		    << ", tau2: " << tau2 << ".\n"
		    << "%" << endl;
		
		uout << "%\n% num_neurons: " << num_neurons
		    << ", coding level: " << coding_lvl
  		    << ", extrate: " << tot_ext_rate
		    << ", num pools: " << num_pools
		    << "% w_{+}: " << wPlus	
		    << ", w inhibitory: " << w_I
		    << ", w_T tuned: " << w_T 
		    << ", AN: " << g_ratio	
		    << ", tau2: " << tau2 << ".\n"
		    << "%" << endl;
		//resetting variables
		
		pfc.resetPools();
		// Reset and resize all the matrices used to keep
		// the firing rates and to store their average values.
		pfc.initialize_vectors(trial_life);
		//----------------------------------------------------------
		//Time course of input 
		//----------------------------------------------------------
		// before stimulus (0.5s), from 0 to 0.5s
		for(int pp=0; pp < num_pools; pp++){
             		pfc.setExtrate(pp, tot_ext_rate);
        	}		
	
		// ANDREA loading param. from an external file:
		string line;
		ifstream protocolfile (protocolFileName.c_str()); 
	        if (protocolfile.is_open())	{
			// ***  protocol.ini***
			getline(protocolfile,line);
			char lineProt[200];
			int poolStim, StimDur=0, Cinput;
			trial_life=0;
			while(protocolfile.getline(lineProt,200)) {
				float NewNu;
				StimDur =0;
				Cinput=sscanf(lineProt, "%d %d %f",&StimDur, &poolStim, &NewNu);
				trial_life+=StimDur;
			}
			cout << "\ntrial_life= " << trial_life  << endl;
			pfc.initialize_vectors(trial_life);
			protocolfile.close();
			ifstream protocolfile (protocolFileName.c_str());
			getline(protocolfile,line);
			int begTime=0;
			while(protocolfile.getline(lineProt,200)) {
				float NewNu;
				StimDur =0;
				Cinput=sscanf(lineProt, "%d %d %f",&StimDur, &poolStim, &NewNu);
				pfc.setExtrate(poolStim,NewNu);
//				cout << "Cinput " << Cinput << endl; // DEBUG MARIO
				if (StimDur!=0){
					pfc.Spiking(begTime,StimDur, &data, &uout, &spdata, &snmda1,&snmda2, &snmda3, facilON, false);
					begTime=begTime+StimDur;
				}
			}
			protocolfile.close();
			cout << "\n" << endl;
		}
		// ANDREA if the external file doesn't exist use the default values:
		else {
			cout << "Unable to open file protocol.ini\n";
			pfc.Spiking(0 ,trial_life , &data, &uout, &spdata, &snmda1,&snmda2, &snmda3, facilON, false);
		}

	        if (facilON == false){		//ANDREA if facilitation is enabled use proper function to write down u[i] time course	
			pfc.write_down(&data);
		}
		else {
			pfc.write_down(&data, &uout);
		}

		// Output
		if (verboseON == true){
                         pfc.printBrain();
                         pfc.print_pools();
		}
		cout << " Data from trial " << k + 1 << " saved in '"
		     << output_file_ext << "'." << endl;

		if ( facilON == true ){
		cout << " Facilitation data from trial " << k + 1 << " saved in '"
			<< output_uout_ext << "'." << endl;
		}

		data.close();
		if (uout.is_open()){
		uout.close();
		}
		if (spdata.is_open()) {
			cout << "Spikes from one neuron in each of the " << 
			    num_pools - 2 << " selective pools, trial " << k +
			    1 << ", saved in file '" <<
			    output_spikes_ext << "'." << endl;
			spdata.close();
		}
		if (snmda_filename != "") {
//			cout << "NMDA component for selective populations "
//			    << "1 and 2 saved into files '"
//			    << snmda_filename_ext_1 << "' and '"
//			    << snmda_filename_ext_2 << "'." << endl;
			snmda1.close();
			snmda2.close();
			snmda3.close();
		}
		cout.flush();
	} // Trial loop
	return 0;
}
