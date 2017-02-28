code to simulate a spiking network (e.g. Wang 2002)

the program is spread in many files

The main file to simulate the spiking network is darco.cc.
This file reads from some .ini files that need to be in the same directory (explained in 
INI_FILES.txt).
The program accept command line options.
The program will generate a text file with firing-rates. Column 1 is time in ms and
subsequent columns are the firing-rate time series for each pool.
Also a file called parametros.txt is generated but this can be just ignored.

Another main file is mfield_darco.cc, which calculate the mean-field equivalent of the 
spiking network. Its output is given through standard output and consist of convergence time
and the firing-rate of each pool once in the attractor.
Initial conditions (and other parameters) can be modified with command line options.
