The command 'make' will compile the program "RunIBM".

To run, the program needs one or two files, a .csf file with all the parameters of the model, and an (optional) .isf file with the initial state of the population.
If there is no isf-file, the population will be initialized with 20 individuals, 10 newborn and 10 adults individuals of each sex. The resource densities are in this case set equal to the Rmax values.
In the absence of an ISF-file, the initial phenotypes are determined by the cvf file. Note that for the mating type, only values of -1, 0, and 1 are possible.
In the presence of an ISF-file, make sure that the ISF-file has the same number of alleles as is indicated in the csf-file. 
 
The command './RunIBM NAME' will run the IBM, where NAME is the name of the cvs (& isf) file. 
For example, './RunIBM Default' will run the ecological dynamics for 10000 days, for the parameters specified in Default.csv. Since the mutation probabilities are zero in this file, there is no evolution. An additional (optional) argument determines the random seed of the run, in case this optional argument is not used, the current time is used to get a random seed.  

The first row of the ISF file has 8 values, the starttime, the size of the system (in m^3), and the densities of the resources, separated by a tab.
The rest of the file contains the data for the initial population. Each line represents one individual with its phenotype & genotype. 

The IBM creates a minimum of 2 output files, one file with the population & trait dynamics over time (NAME_time.txt), and one file with the state of the population (NAME.ESF) at the end of the run. The parameters "Write full pop", "write trait pop", and "write matefile", determine how often the program will write to the files NAME_alltraits.txt, NAME_traits.txt, and NAME_mate.txt respectively. The NAME_alltraits.txt contains the phenotype & genotype of all individuals in the population at a certain time, the NAME_traits.txt contains the phenotype of all individuals in the population at a certain time. The NAME_mate.txt file contains the mating choice of females at a certain time. 


The file with the state of the population (.ESF) has the same structure as the .isf file, with the first line of the file containing the environmental values, and the rest of the file containing the information of the population. These .esf files can be used for a new run.
  



