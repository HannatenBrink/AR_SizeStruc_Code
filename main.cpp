///////INCLUDE FILES//////////////
#include "main.h"

auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//without user input, the time is used as the random seed//

//#define TIMECHECK //In case I want to get an idea about how long calculations take

/*---------------------START OF MAIN-----------------------------------------*/
int main(int argc, char* argv[]) {


/*------------------------------Check input-----------------------------------------*/
  if(argc < 2 ){
    std::cerr << "Usage: " << argv[0] << " Name of the csv/isf file" << ", random seed (optional)" << std::endl;
    exit(1);
  }

/*-------------------Use second argument as filename for output------------------------------*/
  const string Filename = argv[1];
  cout << "Start of run: " << Filename << std::endl;

/*--------------------------------Create random number generators-----------------------------------------------------------------------*/
  //auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  if(argc == 3) {
    seed = atof(argv[2]);
  } //Set seed via command line (optional)
  mt19937 mt_rand(seed);


/*---------------------------Create output file to log info----------------------------------------*/

  //Logfile, containing error messages & random seed//
  ss<<Filename<<"_log"<<".txt";
  LogFile.open(ss.str());
  ss.str("");
  LogFile << "The random seed for this run is: " << seed << endl;
  LogFile << "Get a random number between 0 - 1 to check seed: " << unif(mt_rand) << endl;


  /*---------------------------Read Cvf file----------------------------------------*/
  LogFile << "Reading the CVF file" << endl;
  ss<<Filename<<".cvf";
  ifstream input(ss.str());

  if(!input){
    std::cerr<<"The parameterfile " <<ss.str()<< " does not exist.\n";
    LogFile << "No calculations done, cvf file with parameter values is missing" << endl;
    exit(1);
  }

  string line;

  const string search = "[ 0]"; // test variable to search in file
  unsigned int curLine = 0; //Var to keep track of line numbers
  unsigned int curLine2 = 0; //Var2 to keep track of line numbers
  int Start_read = 0;

  LogFile << endl << "________________________________________" << '\n' <<"Parameters of this run" << endl;

  std::vector<std::string> Parameter_names {"Rho" , "RJmax" , "RAmax" , "ThetaB" ,
  "Theta1" , "Theta2" , "Theta3" , "Theta4" ,"Theta5" , "Theta6" , "TauB" ,
  "Tau1", "Tau2" , "Tau3" , "Tau4" , "Tau5" , "Tau6" ,
  "Amax" , "qpow" , "h" , "npow" , "alphapar" , "kmet" , "pmain" , "u" , "eta" , "epsilon" ,
  "mu_b" , "Xi" , "M_B" , "M_Mat" , "M_Shift" , "s_ass" , "s_diss", "N_eco" , "N_neutral" , "N_mating" ,
  "mut_std" , "Mut_rate_eco" , "Mut_rate_mate" , "Mut_rate_neut", "delta_t" , "volume" , "clonal" , "ini_eco" , "ini_mate" , "ini_neutral", "discrete_shift"};

  int f = 0;
  while(getline(input, line))
    {
    curLine++;
    vector<string> lines;
    if (Start_read == 0){
      split(line, '\t', lines);
      string b = lines.back();
      double AStand = strtof((b).c_str(),0);
      Setting.push_back(AStand); //Read in time settings data in vector 'Setting'
      if (line.find(search) != string::npos) {
        Start_read = 1; //Start with reading parameter values
      }}
      if (Start_read == 1){
        curLine2++;
        vector<string> row_values;
        split(line, '\t', row_values);
        string a =  row_values.back();
        double aNumero=strtof((a).c_str(),0);
        Parameter.push_back(aNumero); //Add parameter value to the vector 'Parameter'
        LogFile << Parameter_names[f] << ": " << '\t' << aNumero << '\t';
        f += 1;
      }
    }

  LogFile << endl << endl;
  LogFile << "Simulation will run for " << Max_time << " days" << endl;
  LogFile << "Stepsize in the run is " << delta_t << endl;
  LogFile << "Writing to timefile every " << Output_time << " days" << endl;
  LogFile << "Writing full pop file every " << Pop_Output << " days" << endl;
  LogFile << "Writing mating file every " << MateFile << " days" << endl;
  ss.str("");
  input.close();
  if (clonal) {
    LogFile << "Reproduction is clonal" << endl;
    if(N_mating>0 || N_neutral > 0) {
      LogFile << "Error: Mating genotypes should be set to zero in cvf file. Exit run..." << endl;
      cout << "Error: Mating genotypes should be set to zero in cvf file. Exit run..." << endl;
      exit(1);
    }} else {
    LogFile << "Reproduction is sexual" << endl;
    if(N_mating==0) {
      LogFile << "Error: Sexual reproduction is only possible in case N_mating > 0. Exit run..." << endl;
      cout << "Error: Sexual reproduction is only possible in case N_mating > 0. Exit run..." << endl;
      exit(1);
    }
  }
  LogFile << endl << "________________________________________" << endl << endl;

  /*---------------------------Create output files----------------------------------------*/

  //File with mate choice data//
  ss<<Filename<<"_mate"<<".txt";
  Matefile.open(ss.str());
  ss.str("");

  //File with resource & consumer dynamics over time//
  ss<<Filename<<"_time"<<".txt";
  Timefile.open(ss.str());
  ss.str("");


  //File containing all individuals and their genotypes&phenotypes in the population//
  ss<<Filename<<"_alltraits"<<".txt";
  FullTraitfile.open(ss.str());
  ss.str("");

  //Population at the end of run//
  ss<<Filename<<".esf";
  Endfile.open(ss.str());
  ss.str("");


 /*--------------------------------Define variables-------------------------------------------*/
   MutNorm.param(std::normal_distribution<double>(0, mut_std).param());
   unif.param(std::uniform_real_distribution<double>(0.0, 1.0).param());
   Surv_age.param(exponential_distribution<double>((mu_b)).param());

   double Tot = 0;
   int mate;
   double Starttime_init = 0;
   vector<double> Feeding;    //Feeding by consumers
   Time = 0; //The time
   double b1 = (ThetaB + Theta1) / 2;
   double b2 = (Theta1 + Theta2) / 2;
   double b3 = (Theta2 + Theta3) / 2;
   double b4 = (Theta3 + Theta4) / 2;
   double b5 = (Theta4 + Theta5) / 2;
   double b6 = (Theta5 + Theta6) / 2;
   SpeciesDiv.push_back(b6);
   SpeciesDiv.push_back(b5);
   SpeciesDiv.push_back(b4);
   SpeciesDiv.push_back(b3);
   SpeciesDiv.push_back(b2);
   SpeciesDiv.push_back(b1);
   double newtime = Max_time / delta_t;


/*------------------------Initialize the environment with ISF file-----------------------------*/
  LogFile << "Initializing the population..." << endl;
  #ifdef TIMECHECK
  auto start = std::chrono::high_resolution_clock::now();
  #endif
  ss<<Filename<<".isf";
  ifstream InitPop(ss.str());
  ss.str("");

  if(!InitPop){
    LogFile<< "The ISF-file " <<ss.str() << " does not exist. Initialization via the program. Start with " << N_ini * 4 <<  " individuals\n";
    Init_Env(mt_rand);
  } else {
    LogFile << "Initialization occurs via the isf file\n";
    Time = Starttime_init;
    double RB_init, R1_init, R2_init, R3_init, R4_init, R5_init, R6_init, Volume_init;
    double sex_init, age_init, size_init, mating_trait_init, eco_trait_init, neutral_trait_init, MaxAge_init;
    double dummy,  idnr,reprobuf;
    int id, Mature, Matings,  nrStarve;
    int i = 0, z = 0;
    int MinItems = 10 + 2 * N_neutral + 2 * N_mating + 2 * N_eco;
    if(N_mating > 0) {
      MinItems += 1;
    }
    if(N_neutral > 0) {
      MinItems +=1;
    }
    if(N_eco > 0) {
      MinItems +=1;
    }
    string data;
    while (getline(InitPop, line)){
      stringstream linestream(line);
      if (i == 1) {
        linestream >> Starttime_init >> Volume_init >> RB_init >> R1_init >> R2_init >> R3_init >> R4_init >> R5_init >> R6_init;
        Time = Starttime_init/delta_t;
        Resource RB(RB_init, RJmax, volume, ThetaB, TauB, "RB");
        Resource R1(R1_init, RAmax, volume, Theta1, Tau1, "R1");
        Resource R2(R2_init, RAmax, volume, Theta2, Tau2, "R2");
        Resource R3(R3_init, RAmax, volume, Theta3, Tau3, "R3");
        Resource R4(R4_init, RAmax, volume, Theta4, Tau4, "R4");
        Resource R5(R5_init, RAmax, volume, Theta5, Tau5, "R5");
        Resource R6(R6_init, RAmax, volume, Theta6, Tau6, "R6");
        AllFood.push_back(RB);
        RmaxChange.push_back(0.0);
        AllFood.push_back(R1);
        RmaxChange.push_back(0.0);
        AllFood.push_back(R2);
        RmaxChange.push_back(0.0);
        AllFood.push_back(R3);
        RmaxChange.push_back(0.0);
        AllFood.push_back(R4);
        RmaxChange.push_back(0.0);
        AllFood.push_back(R5);
        RmaxChange.push_back(0.0);
        AllFood.push_back(R6);
        RmaxChange.push_back(0.0);
      }
      else if (i > 3) {
        if (i == 5){
          auto countline = line;
          stringstream countstream(countline);
          while(getline(countstream, data, '\t')){
            z+=1;
          }
        if (z != MinItems) {
          cout << "The .isf file should have " << MinItems << " variables per individual." << endl;
          cout << "However, there are " << z << " items. The isf file is therefore incorrect." << endl;
          LogFile << "The .isf file should have " << MinItems << " variables per individual." << endl;
          LogFile << "However, there are " << z << " items. The isf file is therefore incorrect." << endl;
          exit(1);}
        }
        vector <double> mating_f;
        vector <double> eco_f;
        vector <double> neutral_f;
        vector <double> mating_m;
        vector <double> eco_m;
        vector <double> neutral_m;
        if(N_neutral && N_mating){
          linestream >> idnr >> sex_init >> age_init >> size_init >>  id >> mating_trait_init >>
              neutral_trait_init >> eco_trait_init >> Mature >> Matings >> reprobuf >> MaxAge_init >> nrStarve;
          for(int j = 0; j < N_mating; ++j) {
                  linestream >> dummy;
                  mating_f.push_back(dummy);
              }
          for(int j = 0; j < N_mating; ++j) {
                  linestream >> dummy;
                  mating_m.push_back(dummy);
              }
          for(int j = 0; j < N_neutral; ++j) {
                  linestream >> dummy;
                  neutral_f.push_back(dummy);
              }
          for(int j = 0; j < N_neutral; ++j) {
                  linestream >> dummy;
                  neutral_m.push_back(dummy);
              }
          for(int j = 0; j < N_eco; ++j) {
                  linestream >> dummy;
                  eco_f.push_back(dummy);
              }
          for(int j = 0; j < N_eco; ++j) {
                  linestream >> dummy;
                  eco_m.push_back(dummy);
              }
            unique_ptr<Individual> IndivPtr(new Individual(idnr, sex_init, age_init, size_init, Mature, Matings, reprobuf, MaxAge_init,  nrStarve, mating_f, mating_m, neutral_f, neutral_m,
              eco_m, eco_f, AllFood));
              if(Mature) {
                if(sex_init==1){
                Femalevec.push_back(move(IndivPtr));
              } else {
                Malevec.push_back(move(IndivPtr));
              }
              } else {
                Juvvec.push_back(move(IndivPtr));
              }
         }
        else if(N_mating){
            linestream >>  idnr >> sex_init >> age_init >> size_init >> id >> mating_trait_init >>
                eco_trait_init >> Mature >> Matings >> reprobuf >> MaxAge_init >> nrStarve;
          for(int j = 0; j < N_mating; ++j) {
                  linestream >> dummy;
                  mating_f.push_back(dummy);
              }
          for(int j = 0; j < N_mating; ++j) {
                  linestream >> dummy;
                  mating_m.push_back(dummy);
              }
          for(int j = 0; j < N_eco; ++j) {

                  linestream >> dummy;
                  eco_f.push_back(dummy);
              }
          for(int j = 0; j < N_eco; ++j) {
                  linestream >> dummy;
                  eco_m.push_back(dummy);
              }
          unique_ptr<Individual> IndivPtr(new Individual(idnr, sex_init, age_init, size_init, Mature, Matings, reprobuf, MaxAge_init,   nrStarve, mating_f, mating_m,
                eco_m, eco_f,  AllFood));
          if(Mature) {
            if(sex_init==1){
            Femalevec.push_back(move(IndivPtr));
          } else {
            Malevec.push_back(move(IndivPtr));
          }
          } else {
            Juvvec.push_back(move(IndivPtr));
          }
         }
        else { //id
           linestream >> idnr >> sex_init >> age_init >> size_init >>  id >>
               eco_trait_init >> Mature >> Matings >> reprobuf >> MaxAge_init >> nrStarve;
           for(int j = 0; j < N_eco; ++j) {
                   linestream >> dummy;
                   eco_f.push_back(dummy);
               }
           for(int j = 0; j < N_eco; ++j) {
                   linestream >> dummy;
                   eco_m.push_back(dummy);
               }
               //cout << "Individual info " << i << '\t' << idnr << '\t' << age_init << '\t' << size_init << '\t' << Mature << '\t' << Matings << '\t' << reprobuf << '\t' << MaxAge_init << '\t' << nrStarve << '\t' << eco_f[0] << '\t' << eco_m[0] << endl;
               unique_ptr<Individual> IndivPtr(new Individual(idnr, sex_init, age_init, size_init, Mature, Matings, reprobuf, MaxAge_init, nrStarve, eco_f, eco_m,   AllFood));
               if(Mature) {
                 Femalevec.push_back(move(IndivPtr));
               } else {
                 Juvvec.push_back(move(IndivPtr));
               }
          }
      }
      i += 1;
    }
  }

  #ifdef TIMECHECK
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  LogFile << "Initialization of the population takes " << duration.count() << " microseconds" << endl;
  #endif

  LogFile << "The initial population size is " << Femalevec.size()+Malevec.size()+Juvvec.size() << endl;
  LogFile << "...End of initializing the population" << endl;
  LogFile << endl << "________________________________________" << endl << endl;


/*-----------------------------------Give names to the files-----------------------------------*/

  LogFile << "Add names of variables to the files, and initialize variables..." << endl;

  Timefile << "Time" << '\t' << "Starve" << '\t' << "OverCons" << '\t' <<
   "PopSize" << '\t' << "PopMass"   << '\t' << "Females"  << '\t' << "Males" << '\t' <<
   "PopSize1" << '\t' << "PopMass1" << '\t' << "Females1" << '\t' << "Males1" << '\t' <<
   "PopSize2" << '\t' << "PopMass2" << '\t' << "Females2" << '\t' << "Males2" << '\t' <<
   "PopSize3" << '\t' << "PopMass3" << '\t' << "Females3" << '\t' << "Males3" << '\t' <<
   "PopSize4" << '\t' << "PopMass4" << '\t' << "Females4" << '\t' << "Males4" << '\t' <<
   "PopSize5" << '\t' << "PopMass5" << '\t' << "Females5" << '\t' << "Males5" << '\t' <<
   "PopSize6" << '\t' << "PopMass6" << '\t' << "Females6" << '\t' << "Males6" << '\t' <<
   "PopSize7" << '\t' << "PopMass7" << '\t' << "Females7" << '\t' << "Males7" << '\t';
  for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
    print_resourceName(Timefile, *it_r);
    Timefile << "\t";}
  #ifdef TIMECHECK
  Timefile << "FeedingDur" << '\t' <<  "Remove_Juv" << '\t' <<  "Remove_Death" << '\t' << "ReproductionDur";
  #endif
  Timefile << endl;

  FullTraitfile << "Time" << "\t";
  print_individualnames(FullTraitfile);
  FullTraitfile << "\n" << endl;

  Matefile <<"Time" <<"\t" << "f_age" << "\t" << "f_size" << "\t" << "f_eco" << "\t" << "f_neu" << "\t" <<
    "f_mating" <<"\t" << "f_reprobuf" << '\t' << "f_ID" << '\t' << "f_IDNR" << '\t' <<
    "m_age" << "\t" << "m_size" << "\t" <<
    "m_eco" << "\t" << "m_neu" << "\t" << "m_mating" << '\t' << "m_reprobuf" << '\t'<<  "m_id" << '\t' << "m_IDNR" << endl;

/*------------------------------Set variables to raise exceptions------------------------------------*/
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sig_to_exception;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;
  sigaction(SIGINT, &sigIntHandler, NULL);
LogFile << "...Finished with adding names of variables to the files, and initializing variables" << endl;
LogFile << "___________________________________________" << endl;



/*------------------------Start of the loop over time-----------------------------------------------*/
 LogFile << "...Start with the simulation" << endl;

 try
 {
   while (Time <= newtime)  {
     #ifdef TIMECHECK
     LogFile << "\n \n \nTime is " << Time * delta_t << endl;
     #endif


     //Write FULL POPstate output to file//
     if ((Pop_Output > 0) && ((round(fmod(T_POP_FULL, (Pop_Output / delta_t))) == 0)|| (round(fmod(T_POP_FULL, (Pop_Output / delta_t)) - (Pop_Output / delta_t))  == 0))){
       LogFile << "Writing full output (phenotype & genotype) to file at time " << Time * delta_t << endl << endl;
       #ifdef TIMECHECK
       start = std::chrono::high_resolution_clock::now();
       #endif
       T_POP_FULL = 0;
       for (auto&& it_f : Femalevec){
         FullTraitfile << Time * delta_t << "\t";
         print_individual(FullTraitfile, *it_f);
         FullTraitfile << endl;
       }
       for (auto&& it_f : Malevec){
         FullTraitfile << Time * delta_t << "\t";
         print_individual(FullTraitfile, *it_f);
         FullTraitfile << endl;
       }
       for (auto&& it_f : Juvvec){
         FullTraitfile << Time * delta_t << "\t";
         print_individual(FullTraitfile, *it_f);
         FullTraitfile << endl;
       }
       #ifdef TIMECHECK
       stop = std::chrono::high_resolution_clock::now();
       duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
       LogFile << "Writing output to the full output takes " << duration.count() << " microseconds" << endl;
       #endif
     }

     //Write output to time file//
     if ((Output_time > 0) && ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) || (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
       #ifdef TIMECHECK
       start = std::chrono::high_resolution_clock::now();
       #endif
       T_Time = 0;
         Pop_Density = Malevec.size() + Femalevec.size() + Juvvec.size();
         Pop_Mass = 0;
         Females = 0;
         Males = 0;

         Pop_Density1 = 0;
         Pop_Mass1 = 0;
         Females1 = 0;
         Males1 = 0;

         Pop_Density2 = 0;
         Pop_Mass2 = 0;
         Females2 = 0;
         Males2 = 0;

         Pop_Density3 = 0;
         Pop_Mass3 = 0;
         Females3 = 0;
         Males3 = 0;

         Pop_Density4 = 0;
         Pop_Mass4 = 0;
         Females4 = 0;
         Males4 = 0;

         Pop_Density5 = 0;
         Pop_Mass5 = 0;
         Females5 = 0;
         Males5 = 0;

         Pop_Density6 = 0;
         Pop_Mass6 = 0;
         Females6 = 0;
         Males6 = 0;

         Pop_Density7 = 0;
         Pop_Mass7 = 0;
         Females7 = 0;
         Males7 = 0;

       for (auto&& it_f : Femalevec){
         Pop_Mass += it_f->size;
         Females += 1;
         switch(it_f->SpeciesID) {
           case 1:
           Pop_Mass1 += it_f->size;
           Pop_Density1 += 1;
           Females1 += 1;
           break;
           case 2:
           Pop_Mass2 += it_f->size;
           Pop_Density2 += 1;
           Females2 += 1;
           break;
           case 3:
           Pop_Mass3 += it_f->size;
           Pop_Density3 += 1;
           Females3 += 1;
           break;
           case 4:
           Pop_Mass4 += it_f->size;
           Pop_Density4 += 1;
           Females4 += 1;
           break;
           case 5:
           Pop_Mass5 += it_f->size;
           Pop_Density5 += 1;
           Females5 += 1;
           break;
           case 6:
           Pop_Mass6 += it_f->size;
           Pop_Density6 += 1;
           Females6 += 1;
           break;
           case 7:
           Pop_Mass7 += it_f->size;
           Pop_Density7 += 1;
           Females7 += 1;
           break;

         }
       }
       for (auto&& it_f : Malevec){
         Pop_Mass += it_f->size;
         Males += 1;
         switch(it_f->SpeciesID) {
           case 1:
           Pop_Mass1 += it_f->size;
           Pop_Density1 += 1;
           Males1 += 1;
           break;
           case 2:
           Pop_Mass2 += it_f->size;
           Pop_Density2 += 1;
           Males2 += 1;
           break;
           case 3:
           Pop_Mass3 += it_f->size;
           Pop_Density3 += 1;
           Males3 += 1;
           break;
           case 4:
           Pop_Mass4 += it_f->size;
           Pop_Density4 += 1;
           Males4 += 1;
           break;
           case 5:
           Pop_Mass5 += it_f->size;
           Pop_Density5 += 1;
           Males5 += 1;
           break;
           case 6:
           Pop_Mass6 += it_f->size;
           Pop_Density6 += 1;
           Males6 += 1;
           break;
           case 7:
           Pop_Mass7 += it_f->size;
           Pop_Density7 += 1;
           Males7 += 1;
           break;

         }
       }
       for (auto&& it_f : Juvvec){
         Pop_Mass += it_f->size;
         switch(it_f->SpeciesID) {
           case 1:
           Pop_Mass1 += it_f->size;
           Pop_Density1 += 1;
           break;
           case 2:
           Pop_Mass2 += it_f->size;
           Pop_Density2 += 1;
           break;
           case 3:
           Pop_Mass3 += it_f->size;
           Pop_Density3 += 1;
           break;
           case 4:
           Pop_Mass4 += it_f->size;
           Pop_Density4 += 1;
           break;
           case 5:
           Pop_Mass5 += it_f->size;
           Pop_Density5 += 1;
           break;
           case 6:
           Pop_Mass6 += it_f->size;
           Pop_Density6 += 1;
           break;
           case 7:
           Pop_Mass7 += it_f->size;
           Pop_Density7 += 1;
           break;

         }
       }


       Timefile << Time * delta_t << '\t' << TotStarv << '\t' << OverCons << '\t'//remove the starve
       << Pop_Density << '\t'  << Pop_Mass  << '\t'  << Females <<'\t' << Males <<'\t'
       << Pop_Density1 << '\t' << Pop_Mass1 << '\t'  << Females1 <<'\t' << Males1 <<'\t'
       << Pop_Density2 << '\t' << Pop_Mass2 << '\t'  << Females2 <<'\t' << Males2 <<'\t'
       << Pop_Density3 << '\t' << Pop_Mass3 << '\t'  << Females3 <<'\t' << Males3 <<'\t'
       << Pop_Density4 << '\t' << Pop_Mass4 << '\t'  << Females4 <<'\t' << Males4 <<'\t'
       << Pop_Density5 << '\t' << Pop_Mass5 << '\t'  << Females5 <<'\t' << Males5 <<'\t'
       << Pop_Density6 << '\t' << Pop_Mass6 << '\t'  << Females6 <<'\t' << Males6 <<'\t'
       << Pop_Density7 << '\t' << Pop_Mass7 << '\t'  << Females7 <<'\t' << Males7 <<'\t' ;

       for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
           print_resourceDensity(Timefile, *it_r);
           Timefile << "\t";
       }
       //Timefile << endl;
       #ifdef TIMECHECK
       stop = std::chrono::high_resolution_clock::now();
       duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
       LogFile << "Writing output to the timefile takes " << duration.count() << " microseconds" << endl;
       #endif
     }

 /*--------------------------Feeding, Growing----------------------------------------------*/
   //Set feeding vector to zero//
   #ifdef TIMECHECK
   start = std::chrono::high_resolution_clock::now();
   #endif
   for(it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
     Feeding.push_back(0);
   }
   #ifdef TIMECHECK
   stop = std::chrono::high_resolution_clock::now();
   duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
   LogFile << "Emptying the feeding vector " << duration.count() << " microseconds" << endl;
   #endif

   #ifdef TIMECHECK
   auto start3 = std::chrono::high_resolution_clock::now();
   auto stop3 = std::chrono::high_resolution_clock::now();
   auto dur3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
   auto tot3 = dur3.count();
   tot3 = 0;
   start = std::chrono::high_resolution_clock::now();
   #endif
   TotStarv = 0;
   OverCons = 0;
   //Feed and grow and die//
   for (auto&& it_f : Femalevec){
     it_f->R_Intake(AllFood, Feeding, mt_rand);
   }
   for (auto&& it_f : Malevec){
     it_f->R_Intake(AllFood, Feeding, mt_rand);
   }
   for (auto&& it_f : Juvvec){
     it_f->R_Intake(AllFood, Feeding, mt_rand);
     if(it_f->Mature){
       if(!clonal){
         if(unif(mt_rand) > 0.5) {it_f->sex = 1;} else {it_f->sex = 0;}
         if(it_f->sex==1){
           Femalevec.push_back(std::move(it_f));
         } else {
           Malevec.push_back(std::move(it_f));
         }
       } else {
         it_f->sex = 1;
         Femalevec.push_back(std::move(it_f));
       }
     }
   }

   #ifdef TIMECHECK
   stop = std::chrono::high_resolution_clock::now();
   duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
   LogFile << "Total feeding takes " << duration.count() << " microseconds" << endl;
   if ((Output_time > 0) &&
      ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
      (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
   Timefile << duration.count() << '\t';}
   #endif
   #ifdef TIMECHECK
   start = std::chrono::high_resolution_clock::now();
   #endif
   Juvvec.erase(std::remove(begin(Juvvec), end(Juvvec), nullptr),
            end(Juvvec));
   #ifdef TIMECHECK
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            LogFile << "Removing empty juv takes " << duration.count() << " microseconds" << endl;
            if ((Output_time > 0) &&
               ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
               (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
            Timefile << duration.count() << '\t';}
   #endif

 /*------------------------Growth of Resources-------------------------------------------*/
   #ifdef TIMECHECK
   start = std::chrono::high_resolution_clock::now();
   #endif

   for(i = 0, it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r, ++i){it_r->Growth(Feeding[i], RmaxChange[i]);}
   Feeding.clear(); //empty the feeding vector

   #ifdef TIMECHECK
   stop = std::chrono::high_resolution_clock::now();
   duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
   LogFile << "Change in resource densities takes " << duration.count() << " microseconds" << endl;
   #endif

 /*-------------------------------Dying---------------------------------------------------------*/
   #ifdef TIMECHECK
   start = std::chrono::high_resolution_clock::now();
   #endif

   #ifdef TIMECHECK
   stop = std::chrono::high_resolution_clock::now();
   duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
   LogFile << "Dying takes " << duration.count() << " microseconds" << endl;
   #endif

   #ifdef TIMECHECK
   start = std::chrono::high_resolution_clock::now();
   #endif

   Femalevec.erase(std::remove_if(Femalevec.begin(), Femalevec.end(), [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Femalevec.end());
   Malevec.erase(std::remove_if(Malevec.begin(), Malevec.end(), [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Malevec.end());
   Juvvec.erase(std::remove_if(Juvvec.begin(), Juvvec.end(), [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Juvvec.end());


   #ifdef TIMECHECK
   stop = std::chrono::high_resolution_clock::now();
   duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
   LogFile << "Removing death animals takes " << duration.count() << " microseconds" << endl;
   LogFile << "Size of the new population is " << Femalevec.size() + Malevec.size() << endl;
   if ((Output_time > 0) &&
      ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
      (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
   Timefile << duration.count() << '\t';}
   #endif

   //Exit when population is extinct
   if(Femalevec.size()<=0 && Malevec.size()<=0 && Juvvec.size()<=0){
     std::cerr<<"The consumer population is extinct at time " << Time * delta_t << "\n";
     exit(1);
   }

 /*------------------------Mating----CHANGE ---------------------------*/
   if(!clonal){//sexual reproduction
     #ifdef TIMECHECK
     start = std::chrono::high_resolution_clock::now();
     auto start2 = std::chrono::high_resolution_clock::now();
     auto stop2 = std::chrono::high_resolution_clock::now();
     auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
     auto Tot2 = dur2.count();
     Tot2 = 0;
     #endif
     //maybe better to shuffle it depending on reproductive buffer? In this way maybe less mate limitation?
     sort(Femalevec.begin(), Femalevec.end(), cmp_by_repro); //sort
     //sort(Malevec.begin(), Malevec.end(), cmp_by_repro); //sorting of males not necessary
   for (auto&& it_f : Femalevec) {
     if (it_f->Fecund) {
     vector<double> cumsum;
     Tot = 0;
     #ifdef TIMECHECK
     start2 = std::chrono::high_resolution_clock::now();
     #endif
     for (auto&& it_m : Malevec){
       it_m->MateProb(*it_f, cumsum, Tot); //calculate probabilities
     }
     #ifdef TIMECHECK
     stop2 = std::chrono::high_resolution_clock::now();
     dur2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
     Tot2 += dur2.count();
     #endif
     RandomVal = unif(mt_rand) * Tot;
     mate = weighted_random_known_sums_floats(cumsum, RandomVal, cumsum.size()); //roulette wheel selection with binary search
     auto it_m = next(Malevec.begin(), mate);
     if (((*it_m)->Fecund) & (Tot > 0)) {
       it_f->Matings += 1;
       (*it_m)->Matings += 1;

     //Write mate choice to matefile
     if ((MateFile > 0) && ((round(fmod(T_Mate, (MateFile / delta_t))) == 0)|| (round(fmod(T_Mate, (MateFile / delta_t)) - (MateFile / delta_t))  == 0))) {
         Matefile << Time * delta_t << "\t" <<
         it_f->age << "\t"
         << it_f->size << "\t" <<
         it_f->ecological_trait << "\t"
         << it_f->neutral_trait << "\t" <<
         it_f->mating_trait << "\t" <<
         it_f->repro_buffer << "\t" <<
         it_f->SpeciesID << "\t" <<
         it_f->IDNR << "\t" <<
         (*it_m)->age << "\t" <<
         (*it_m)->size << "\t" <<
         (*it_m)->ecological_trait << "\t" <<
         (*it_m)->neutral_trait << "\t" <<
         (*it_m)->mating_trait << "\t" <<
         (*it_m)->repro_buffer << "\t" <<
         (*it_m)->SpeciesID << "\t" <<
         (*it_m)->IDNR << "\t" <<
         endl;}
         it_f->SexualRepro(**it_m, mt_rand); //reproduce
     } else if((MateFile > 0) && ((round(fmod(T_Mate, (MateFile / delta_t))) == 0)|| (round(fmod(T_Mate, (MateFile / delta_t)) - (MateFile / delta_t))  == 0))) {
         Matefile << Time * delta_t << "\t" <<
         it_f->age << "\t" <<
         it_f->size << "\t" <<
         it_f->ecological_trait << "\t" <<
         it_f->neutral_trait << "\t" <<
         it_f->mating_trait << "\t" <<
         it_f->repro_buffer << "\t" <<
         it_f->SpeciesID << "\t" <<
         it_f->IDNR << "\t" <<
         0 << "\t" << 0 << "\t" <<
         0 << "\t" <<
         0 << "\t" <<
         0 << "\t" <<
         0 << "\t" <<
         0 << "\t" <<
         0 << endl;}
       }}

   if ((MateFile > 0) && ((round(fmod(T_Mate, (MateFile / delta_t))) == 0)||(round(fmod(T_Mate, (MateFile / delta_t)) - (MateFile / delta_t))  == 0))) {
       T_Mate = 0;
       LogFile << "Writing mating combinations to matefile at time " << Time * delta_t << endl << endl;
     }

   #ifdef TIMECHECK
   stop = std::chrono::high_resolution_clock::now();
   duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
   LogFile << "Sexual reproduction takes " << duration.count() << " microseconds" << endl;
   LogFile << "Of which " << Tot2 << " microseconds is calculating mating probabilities" << endl;
   if ((Output_time > 0) &&
      ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
      (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
   Timefile << duration.count();}
   #endif
   }
   else {//clonal reproduction
     #ifdef TIMECHECK
     start = std::chrono::high_resolution_clock::now();
     #endif
     for (auto&& it_f : Femalevec){
       if (it_f->Fecund) {
         it_f->ClonalRepro(mt_rand);
         it_f->Matings += 1;}
     }
     for (auto&& it_f : Malevec){ //will not be used but just for clarity.
       if (it_f->Fecund) {
         it_f->ClonalRepro(mt_rand);
         it_f->Matings += 1;}
     }

     #ifdef TIMECHECK
     stop = std::chrono::high_resolution_clock::now();
     duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
     LogFile << "Clonal reproduction takes " << duration.count() << " microseconds" << endl;
     if ((Output_time > 0) &&
        ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
        (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
     Timefile << duration.count();}
     #endif
   }

   if ((Output_time > 0) &&
      ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
      (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
       Timefile << endl;}

 /*------------------------Increase time-------------------------------------------------*/
   Time += 1;
   T_Time +=1;
   T_POP_FULL += 1;
   T_POP += 1;
   T_Mate += 1;

 }

 LogFile << "... Finished with the simulation" << endl << endl;

 /*------------------------End of the timeloop-------------------------------------------*/

 /*------------------------Create last output--------------------------------------------*/


 Endfile << "Time" << "\t" << "Volume" << "\t";
 for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
   print_resourceName(Endfile, *it_r);
   Endfile << "\t";
 }
 Endfile << endl;

 Endfile << Time * delta_t << "\t" << volume << "\t";
 for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
   print_resourceDensity(Endfile, *it_r);
   Endfile << "\t";
 }
 Endfile << endl << endl;

 print_individualnames(Endfile);
 Endfile << endl;

 for (auto&& it_f : Femalevec){
   print_individual(Endfile, *it_f);
   Endfile << endl;
 }
 for (auto&& it_f : Malevec){
   print_individual(Endfile, *it_f);
   Endfile << endl;
 }
 for (auto&& it_f : Juvvec){
   print_individual(Endfile, *it_f);
   Endfile << endl;
 }


 LogFile << "\nSimulation ended at time " << Time * delta_t - delta_t << endl;
 cout << "Simulation ended at time " << Time * delta_t - delta_t << endl;

 //close files
   Endfile.close();
   FullTraitfile.close();
   Timefile.close();
   LogFile.close();
   Matefile.close();
 }

 catch(InterruptException& e) //Try to write to output in case of interruption
     {
       LogFile  << "Caught signal " << e.S << endl;
       LogFile << "Writing current output to the .esf file" << endl;

       //Resource Names to file///
       Endfile << "Time" << "\t" << "Volume" << "\t";
       for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
         print_resourceName(Endfile, *it_r);
         Endfile << "\t";
       }
       Endfile << endl;
       //Resource data to file//
       Endfile << Time*delta_t << "\t" << volume << "\t";
       for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
         print_resourceDensity(Endfile, *it_r);
         Endfile << "\t";
       }
       Endfile << endl;

       //individual names to file///
       print_individualnames(Endfile);
       Endfile << "\n" << endl;
       //individual
       for (auto&& it_f : Femalevec){
         print_individual(Endfile, *it_f);
         Endfile << endl;
       }
       for (auto&& it_f : Malevec){
         print_individual(Endfile, *it_f);
         Endfile << endl;
       }
       for (auto&& it_f : Juvvec){
         print_individual(Endfile, *it_f);
         Endfile << endl;
       }

       LogFile  << "Simulation ended at time " << Time * delta_t << endl;
       cout << "Simulation ended at time " << Time * delta_t << endl;

       Endfile.close();
       FullTraitfile.close();
       Timefile.close();
       LogFile.close();
       Matefile.close();

       return 1;
     }
   return 0;
 }
