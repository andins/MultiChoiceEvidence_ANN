#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>		// for atoi() and atof()
#include <string>
#include <unistd.h>		// for getopt()

#include "brain.h"
#include "setup.h"

void help_message_mfield(char *program)
{
	cerr << "usage: " << program
	    << "\n\t [-h (for help)] \n"
	    << "\t [-l lambda_average_in_Hz]\n"
	    << "\t [-b lambda_bias_in_Hz]\n"
	    << "\t [-u lambda_uncertain_in_Hz]\n"
	    << "\t [-o output_file]\n"
	    << "\t [-w window_time_size_in_ms]\n"
	    << "\t [-j potentiated_synaptic_efficacy]\n"
	    << "\t [-T enhanced_synapses(neighboring_pools)]\n"
	    << "\t [-s number_of_euler_steps]\n"
	    << "\t [-e number_external_neurons]\n"
	    << "\t [-i init_rate_of_1st_pool]\n"
	    << "\t [-k init_rate_of_2nd_pool]\n"
	    << "\t [-g init_rate_of_3rd_pool]\n"
	    << "\t [-z init_rate_of_0th_pool]\n"
	    << "\t [-x init_rate_of_4th_pool]\n"
	    << "\t [-N g_AMPA_factor_and_suppress_NMDA]" << endl;
}


bool read_options_mfield(int argc, char *argv[], double &l_average,
			 double &l_bias, double &l_uncertain, int &num_ext_neurns,
			 string & output_file, int output_flag,
			 int &window_size, double &wPlus, double &w_T, double &eps,
			 int &total_steps, double &initial_rate1,
			 double &initial_rate2, double &initial_rate3, double &initial_rate0, double &initial_rate4, int &self_termination,
			 double &g_ratio)
{
	int c;
	// Colons in the third getopt argument allow args for the options
	while ((c = getopt(argc, argv, "hl:b:u:e:o:w:j:T:E:s:i:k:g:z:x:N:?")) != -1)
		switch (c) {
		case 'h':
			help_message_mfield(argv[0]);
			return false;	// execution should stop after help
		case 'l':
			l_average = atof(optarg);
			break;
		case 'b':
			l_bias = atof(optarg);
			break;
		case 'u':
			l_uncertain = atof(optarg);
			break;
		case 'e':
			num_ext_neurns = atoi(optarg);
			break;
		case 'o':
			output_file = optarg;
			output_flag = 1;
			break;
		case 'w':
			window_size = atoi(optarg);
			break;
		case 'j':
			wPlus = atof(optarg);
			break;
		case 'T':
			w_T = atof(optarg);
			break;
		case 'E':
			eps = atof(optarg);
			break;
		case 's':
			total_steps = atoi(optarg);
			self_termination = 0;	// explicit nmber of steps given
			break;
		case 'i':
			initial_rate1 = atof(optarg);
			break;
		case 'k':
			initial_rate2 = atof(optarg);
			break;
		case 'g':
			initial_rate3 = atof(optarg);
			break;
		case 'z':
			initial_rate0 = atof(optarg);
			break;
		case 'x':
			initial_rate4 = atof(optarg);
			break;
		case 'N':
			g_ratio = atof(optarg);
			break;
		case '?':
			help_message_mfield(argv[0]);
			return false;	// execution should stop after error
		}

	return true;		// ready to continue program execution
}

int main(int argc, char *argv[])
{
	// ---------------------------
	//   Initialize pfc
	// ---------------------------
	double time_step = 0.02;	// in ms
	double inh_percent = 0.20;	// inhibitory population (over 1).
	double coding_lvl = 0.20;	// coding level of every sel. pool

	int num_pools=0; 		//decision making 2 alternatives

	int num_neurons = 1000;

	// Initialize arguments to appropriate default values
	// These values may be changed as command-line options
	int num_ext_neurns = 800;	// option -e 800
	double l_average = 30.0;	// option -l 0.0
	double l_bias = 0.0;	// option -b 3.0
	double l_uncertain = 0.0;	// option -u 0.0

	bool verboseON = false; // ANDREA output off at the command line if is not -v 
        string protocolFileName = "protocol.ini"; // ANDREA protocolFileName can change at the command line
        string connectivityFileName = "connectivity.ini"; // ANDREA connectivityFileName can change at the command line

        double ratio = 1;       // ANDREA the default value for the ratio

	int output_flag = 0;	// option -o typed or not
	string output_file = "data.dat";	// option -o data.dat

	int window_size = 50;	// option -w 50
	int total_steps = 4000;	// option -s 4000
	int self_termination = 1;	// this flag is on unless -s option
	double wPlus = 1.5;	// option -j 1.6
	double w_T = 0.0;	// option -T 0.01
	double w_I = 1.0;
	double eps = 0.0;		// enhanced competition option -E
	double initial_rate0 = 19.0;	// option -z 3.0
	double initial_rate1 = 3.0;	// option -i 3.0
	double initial_rate2 = 3.0;	// option -k 3.0
	double initial_rate3 = 3.0;	// option -g 3.0
	double initial_rate4 = 3.0;	// option -x 3.0
	double g_ratio = 0.0;

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

//	cout << "coding_lvl: " << coding_lvl  << endl;
	

	if (!read_options_mfield(argc, argv, l_average, l_bias, l_uncertain,
				 num_ext_neurns, output_file,
				 output_flag, window_size, wPlus, w_T, eps,
				 total_steps, initial_rate1,
				 initial_rate2, initial_rate3, initial_rate0, initial_rate4, self_termination,
				 g_ratio))
		return 1;	// halt criterion detected by read_options()


	// TODO: Check the input
	if (g_ratio < 0 || g_ratio > 1){
		cerr << "The argument for the -N option must be in the range [0,1].\n";
		exit(1);
	}

	setup(num_neurons, num_pools, num_ext_neurns, time_step,
	      inh_percent, coding_lvl, g_ratio, ratio);


	// --------------------------------------------------
	// Set the connections
	// --------------------------------------------------
	set_all_connections(num_pools, coding_lvl, wPlus, w_I, w_T,connectivityFileName);	

	
	// ------------------------------------------------------------
	// Start the simulation
	// ------------------------------------------------------------
	ofstream data;
	data.open(output_file.c_str(), ios::out);

	if (output_flag == 1) {
		data.open(output_file.c_str(), ios::out);
	}

	pfc.setExtrate(1, (2400.0 + l_average + l_bias) /
		       (double) num_ext_neurns);
	pfc.setExtrate(2, (2400.0 + l_average - l_bias) /
		       (double) num_ext_neurns);
	pfc.setExtrate(3, (2400.0 + l_uncertain) /
		       (double) num_ext_neurns);
	pfc.setStartPoolrate(0, initial_rate0);
	pfc.setStartPoolrate(1, initial_rate1);
	pfc.setStartPoolrate(2, initial_rate2);
	pfc.setStartPoolrate(3, initial_rate3);
	pfc.setStartPoolrate(4, initial_rate4);
	pfc.resetPools();
	pfc.calc_mf_parameters();
	//                 pfc.printBrain();
	//                 pfc.print_pools();

 
//	if (output_flag == 1) {						******** OLD FUNCTIONS (THE PROGRAM DOESN'T CHANGE THE output_flag WHEN
//		pfc.MeanField(total_steps, 0.1, &data);						-o ARGUMENT IS TYPED)[???]	***************	
//		cout << "\n Mean Field results saved into the file '"	***********************************************************************
//		    << output_file << "'.\n";
//		data.close();
//	} else {
		if (self_termination == 1) {
			int counts;
			counts = pfc.get_fixed_rates(0.1, -1);	// Evolves to stationarity
			cout << fixed << setprecision(1);
			cout << setw(4) << counts * 0.1;
			pfc.cout_rates();
			data << fixed << setprecision(1);
			data << setw(4) << counts * 0.1;
			pfc.write_rates(&data);
			data.close();
		} else {	// use the given number of steps
			pfc.MeanField(total_steps, 0.1, &cout);
			cout << fixed << setprecision(1);
			cout << setw(4) << total_steps * 0.1;
			pfc.cout_rates();
			data << fixed << setprecision(1);
			data << setw(4) << total_steps * 0.1;
			pfc.write_rates(&data);
			data.close();
		}
//	}
	return 0;
}
