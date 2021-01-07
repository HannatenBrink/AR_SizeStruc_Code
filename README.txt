The command 'make' will compile the program "RunIBM".

To run, the program needs one or two files, a .cvf file with all the parameters of the model, and an (optional) .isf file with the initial state of the population.
If there is no isf-file, the population will be initialized with 20 individuals, 10 newborn and 5 adults of each sex. The resource densities are in this case set equal to the Rmax values. In case of clonal reproduction, all individuals are females. 
In the absence of an ISF-file, the initial phenotypes are determined by the cvf file. Note that for the mating type, only values of -1, 0, and 1 are possible.
In the presence of an ISF-file, make sure that the ISF-file has the same number of alleles as indicated in the cvf-file, otherwise, initialisation will not work properly.
 
The command './RunIBM NAME' will run the IBM, where NAME is the name of the cvf (& isf) file. 
For example, './RunIBM Default' will run the ecological dynamics for 10000 days, for the parameters specified in Default.cvf. Since the mutation probabilities are zero in this file, there is no evolution. An additional (optional) argument determines the random seed of the run, in case this optional argument is not used, the current time is used to get a random seed (saved in the log file).  

The first row of the ISF file has 9 values, the starttime, the size of the system (in m^3), and the densities of the resources, separated by a tab.
The rest of the file contains the data for the initial population. Each line represents one individual with its phenotype & genotype. 

The IBM creates a minimum of 3 output files, a logfile (NAME_log.txt) with parameter values and settings, one file with the population & resource dynamics over time (NAME_time.txt), and one file with the state of the population (NAME.esf) at the end of the run. The parameters "write full pop" and "write matefile" determine how often the program will write to the files NAME_alltraits.txt, and NAME_mate.txt respectively. The NAME_alltraits.txt contains the phenotype & genotype of all individuals in the population at a certain time. The NAME_mate.txt file contains the mating choice of females at a certain time (only in case of sexual reproduction). 

The file with the state of the population (.esf) has the same structure as the .isf file, with the first line of the file containing the environmental values, and the rest of the file containing the information of the population (phenotype & genotype). These .esf files can be used for a new run.

Stepsize is super important! Too big and you'll get strange dynamics. It depends on model parameters which stepsize works. A stepsize of 0.01 seems to always work well but for some parameters 0.1 or 1 works. Since bigger stepsize speeds up calculations considerably, it is worthwhile to first do some testruns to determine the stepsize. 

(Not yet in this repository): Deterministic dynamics in absence of evolution can be found using the ebttool. The files determining the ecological dynamics called XXXX.
  



