///////INCLUDE FILES//////////////
#include "main.h"


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
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  if(argc == 3) {
    seed = atof(argv[2]);
  } //Set seed via command line (optional)

/*---------------------------Create output files----------------------------------------*/
  //File with mate choice data//
  ofstream Matefile;
  ss<<Filename<<"_mate"<<".txt";
  Matefile.open(ss.str());
  ss.str("");

  //File with resource & consumer dynamics over time//
  ofstream Timefile;
  ss<<Filename<<"_time"<<".txt";
  Timefile.open(ss.str());
  ss.str("");

  //File containing all individuals and their phenotypes in the population//
  ofstream Traitfile;
  ss<<Filename<<"_traits"<<".txt";
  Traitfile.open(ss.str());
  ss.str("");

  //File containing all individuals and their genotypes&phenotypes in the population//
  ofstream FullTraitfile;
  ss<<Filename<<"_alltraits"<<".txt";
  FullTraitfile.open(ss.str());
  ss.str("");

  //Population at the end of run//
  ofstream Endfile;
  ss<<Filename<<".esf";
  Endfile.open(ss.str());
  ss.str("");

  //Logfile, containing error messages & random seed//
  ofstream LogFile;
  ss<<Filename<<"_log"<<".txt";
  LogFile.open(ss.str());
  ss.str("");

  LogFile << "The random seed for this run is: " << seed << endl;
  mt19937 mt_rand(seed);

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
  "mu_b" , "Xi" , "M_B" , "M_Mat" , "M_Shift" , "s_ass" , "N_eco" , "N_neutral" , "N_mating" ,
  "mut_std" , "Mut_rate" , "Mut_rate_dialic" , "delta_t" , "volume" , "clonal" , "ini_eco" , "ini_mate" , "ini_neutral"};

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
  LogFile << "Simulation will run for " << Max_time << " years" << endl;
  LogFile << "Stepsize in the run is " << delta_t << endl;
  LogFile << "Writing to timefile every " << Output_time << " days" << endl;
  LogFile << "Writing full pop file every " << Pop_Output << " days" << endl;
  LogFile << "Writing mating file every " << MateFile << " days" << endl;
  ss.str("");
  input.close();
  if (clonal) {
    LogFile << "Reproduction is clonal" << endl;} else {
    LogFile << "Reproduction is sexual" << endl;
    if(N_mating==0) {
      LogFile << "Error: Sexual reproduction is only possible in case N_mating > 0. Exit run..." << endl;
      cout << "Error: Sexual reproduction is only possible in case N_mating > 0. Exit run..." << endl;
      exit(1);
    }
  }
  LogFile << endl << "________________________________________" << endl << endl;

 /*--------------------------------Define variables-------------------------------------------*/
   MutNorm.param(std::normal_distribution<double>(0,mut_std).param());
   unif.param(std::uniform_real_distribution<double>(0.0, 1.0).param());
   Surv_age.param(exponential_distribution<double>((mu_b)).param());

   double Tot = 0;
   int mate;
   double Starttime_init = 0;
   vector<double> Feeding;    //Feeding by consumers
   double Time = 0; //The time
   double b1 = (Theta1 + Theta2) / 2;
   double b2 = (Theta2 + Theta3) / 2;
   double b3 = (Theta3 + Theta4) / 2;
   double b4 = (Theta4 + Theta5) / 2;
   double b5 = (Theta5 + Theta6) / 2;
   SpeciesDiv.push_back(b5);
   SpeciesDiv.push_back(b4);
   SpeciesDiv.push_back(b3);
   SpeciesDiv.push_back(b2);
   SpeciesDiv.push_back(b1);
   double newtime = Max_time / delta_t + delta_t;


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
    Init_Env();
  } else {
    LogFile << "Initialization occurs via the isf file\n";

    double RB_init, R1_init, R2_init, R3_init, R4_init, R5_init, R6_init, Volume_init;
    double age_init, size_init, mating_trait_init, eco_trait_init, neutral_trait_init, MaxAge_init;
    double dummy;
    int sex_init, id, Mature, Matings;
    int i = 0;
    while (getline(InitPop, line)){
      stringstream linestream(line);
      if (i == 1) {
        linestream >> Starttime_init >> Volume_init >> RB_init >> R1_init >> R2_init >> R3_init >> R4_init >> R5_init >> R6_init;
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
      if (i > 3) {
        vector <int> mating_f;
        vector <double> eco_f;
        vector <double> neutral_f;
        vector <int> mating_m;
        vector <double> eco_m;
        vector <double> neutral_m;
        if(N_neutral && N_mating){
          linestream >> age_init >> size_init >> sex_init >> id >> mating_trait_init >>
              neutral_trait_init >> eco_trait_init >> Mature >> Matings >> MaxAge_init;
          //cout << "Age_init is " << age_init << endl;
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
          if(sex_init == 1) {
            unique_ptr<Individual> IndivPtr(new Individual(age_init, size_init, sex_init, MaxAge_init, Mature, mating_f, mating_m, neutral_f, neutral_m,
              eco_m, eco_f, AllFood));
              Femalesvec.push_back(move(IndivPtr));
          } else {
            unique_ptr<Individual> IndivPtr(new Individual(age_init, size_init, sex_init, MaxAge_init, Mature, mating_f, mating_m, neutral_f, neutral_m,
              eco_m, eco_f, AllFood));
            Malesvec.push_back(move(IndivPtr));
          }
         }
        else if(N_mating){
          linestream >> age_init >> size_init >> sex_init >> id >> mating_trait_init >>
              eco_trait_init >> Mature >> Matings >> MaxAge_init;
          //cout << "Age_init is " << age_init << endl;
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
          unique_ptr<Individual> IndivPtr(new Individual(age_init, size_init, sex_init, MaxAge_init, Mature, mating_f, mating_m,
                eco_m, eco_f,  AllFood));
          if(sex_init == 1) {
              Femalesvec.push_back(move(IndivPtr));
          } else {
              Malesvec.push_back(move(IndivPtr));
          }
         }
        else {
           linestream >> age_init >> size_init >> sex_init >> id >>
               eco_trait_init >> Mature >> Matings >> MaxAge_init;
           for(int j = 0; j < N_eco; ++j) {
                   linestream >> dummy;
                   eco_f.push_back(dummy);
               }
           for(int j = 0; j < N_eco; ++j) {
                   linestream >> dummy;
                   eco_m.push_back(dummy);
               }
               unique_ptr<Individual> IndivPtr(new Individual(age_init, size_init, sex_init, MaxAge_init, Mature, eco_f, eco_m,   AllFood));
           if(sex_init == 1) {
               Femalesvec.push_back(move(IndivPtr));
           } else {
              Malesvec.push_back(move(IndivPtr));
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

  LogFile << "The initial population size is " << Femalesvec.size() + Malesvec.size() << endl;
  LogFile << "...End of initializing the population" << endl;
  LogFile << endl << "________________________________________" << endl << endl;


/*-----------------------------------Give names to the files-----------------------------------*/

  LogFile << "Add names of variables to the files, and initialize variables..." << endl;

  Timefile << "Time" << '\t' << "Females" << '\t' << "FemaleMass" << '\t' <<
   "Males" << '\t' << "MaleMass" << '\t' << "Adults" << '\t';
  for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
    print_resourceName(Timefile, *it_r);
    Timefile << "\t";}
  #ifdef TIMECHECK
  Timefile << "FeedingDur" << '\t' <<  "Remove_Death" << '\t' << "ReproductionDur";
  #endif
  Timefile << endl;

  FullTraitfile << "Time" << "\t";
  print_individualnames(FullTraitfile);
  FullTraitfile << "\n" << endl;

  Traitfile << "Time" << "\t";
  print_traitsindividualnames(Traitfile);
  Traitfile << endl;

  Matefile <<"Time" <<"\t" << "f_age" << "\t" << "f_size" << "\t" << "f_eco" << "\t" << "f_neu" << "\t" << "f_mating" <<"\t" << "m_age" << "\t" << "m_size" << "\t" <<
    "m_eco" << "\t" << "m_neu" << "\t" << "m_mating" << endl;

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
    Time = Starttime_init;

    while (Time <= newtime)  {
      #ifdef TIMECHECK
      LogFile << "\n \n \nTime is " << Time * delta_t << endl;
      #endif

      //Write phenotypic traits only to traitfile//
      if ((Pop_TraitOutput> 0) && ((round(fmod(T_POP, (Pop_TraitOutput / delta_t))) == 0) || ((round(fmod(T_POP,   (Pop_TraitOutput / delta_t)) - (Pop_TraitOutput / delta_t))  == 0)))) {
      LogFile << "Writing phenotypes to traitfile at time " << Time * delta_t << endl << endl;
      #ifdef TIMECHECK
      start = std::chrono::high_resolution_clock::now();
      #endif
      T_POP = 0;
      for (auto&& it_f : Femalesvec){
        Traitfile << Time * delta_t << "\t";
        print_traitsindividual(Traitfile, *it_f);
        Traitfile << endl;
      }
      for (auto&& it_m : Malesvec){
        Traitfile << Time * delta_t << "\t";
        print_traitsindividual(Traitfile, *it_m);
        Traitfile << endl;
      }
      #ifdef TIMECHECK
      stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      LogFile << "Writing output to the traitfile took " << duration.count() << " microseconds" << endl;
      #endif
    }

      //Write FULL POPstate output to file//
      if ((Pop_Output > 0) && ((round(fmod(T_POP_FULL, (Pop_Output / delta_t))) == 0)|| (round(fmod(T_POP_FULL, (Pop_Output / delta_t)) - (Pop_Output / delta_t))  == 0))){
        LogFile << "Writing full output (phenotype & genotype) to file at time " << Time * delta_t << endl << endl;
        #ifdef TIMECHECK
        start = std::chrono::high_resolution_clock::now();
        #endif
        T_POP_FULL = 0;    
        for (auto&& it_f : Femalesvec){
          FullTraitfile << Time * delta_t << "\t";
          print_individual(FullTraitfile, *it_f);
          FullTraitfile << endl;
        }
        for (auto&& it_m : Malesvec){
          FullTraitfile << Time * delta_t << "\t";
          print_individual(FullTraitfile, *it_m);
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
        double Female_Density = Femalesvec.size();
        double Male_Density = Malesvec.size();
        double Female_Mass = 0;
        double Male_Mass = 0;
        double Adults = 0;
        for (auto&& it_f : Femalesvec){
          Female_Mass += it_f->size;
          if (it_f->Mature){
            Adults += 1;
          }
        }
        for (auto&& it_m : Malesvec){
          Male_Mass += it_m -> size;
          if (it_m->Mature){
            Adults += 1;
          }
        }

        Timefile << Time * delta_t << '\t' << Female_Density << '\t'
        << Female_Mass << '\t'  << Male_Density << '\t' << Male_Mass << '\t' << Adults <<'\t' ;
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
    //Feed and grow and die//
    for (auto&& it_f : Femalesvec){
      it_f->R_Intake(AllFood, Feeding);
    }
    for (auto&& it_m : Malesvec){
      it_m->R_Intake(AllFood, Feeding);
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




/*------------------------Growth of Resources-------------------------------------------*/
    #ifdef TIMECHECK
    start = std::chrono::high_resolution_clock::now();
    #endif

    for(i = 0, it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r, ++i){
        it_r->Growth(Feeding[i], RmaxChange[i]);
    }
    Feeding.clear(); //empty the feeding vector

    #ifdef TIMECHECK
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    LogFile << "Change in resource densities takes " << duration.count() << " microseconds" << endl;
    #endif

/*------------------------Dying---------------------------------------------------------*/

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

    Femalesvec.erase(std::remove_if(Femalesvec.begin(), Femalesvec.end(),
    [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Femalesvec.end());

    Malesvec.erase(std::remove_if(Malesvec.begin(), Malesvec.end(),
    [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Malesvec.end());


    #ifdef TIMECHECK
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    LogFile << "Removing death animals takes " << duration.count() << " microseconds" << endl;
    if ((Output_time > 0) &&
       ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
       (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
    Timefile << duration.count() << '\t';}
    #endif

/*------------------------Mating-------------------------------*/
    if(!clonal){
      #ifdef TIMECHECK
      start = std::chrono::high_resolution_clock::now();
      auto start2 = std::chrono::high_resolution_clock::now();
      auto stop2 = std::chrono::high_resolution_clock::now();
      auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
      auto Tot2 = dur2.count();
      Tot2 = 0;
      #endif
    for (auto&& it_f : Femalesvec) {
      if (it_f->Fecund) {
      vector<double> cumsum;
      Tot = 0;

      #ifdef TIMECHECK
      start2 = std::chrono::high_resolution_clock::now();
      #endif
      for (auto&& it_m : Malesvec){
        it_m->MateProb(*it_f, cumsum, Tot); //calculate probabilities
      }
      #ifdef TIMECHECK
      stop2 = std::chrono::high_resolution_clock::now();
      dur2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
      Tot2 += dur2.count();
      #endif
      RandomVal = unif(mt_rand) * Tot;
      mate = weighted_random_known_sums_floats(cumsum ,RandomVal, cumsum.size()); //roulette wheel selection with binary search
      auto it_m = next(Malesvec.begin(), mate);
      it_f->Matings += 1;
      (*it_m)->Matings += 1;
      //Write mate choice to matefile
      if ((MateFile > 0) && ((round(fmod(T_Mate, (MateFile / delta_t))) == 0)|| (round(fmod(T_Mate, (MateFile / delta_t)) - (MateFile / delta_t))  == 0))) {
          Matefile << Time * delta_t << "\t" << it_f->age << "\t" << it_f->size << "\t" <<
          it_f->ecological_trait << "\t"
          << it_f->neutral_trait << "\t" <<
          it_f->mating_trait << "\t" <<
          (*it_m)->age << "\t" << (*it_m)->size << "\t" <<
          (*it_m)->ecological_trait << "\t" <<
          (*it_m)->neutral_trait << "\t" <<
          (*it_m)->mating_trait << endl;
      }
      it_f->SexualRepro(**it_m);
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
    else {
      #ifdef TIMECHECK
      start = std::chrono::high_resolution_clock::now();
      #endif
      for (auto&& it_f : Femalesvec){
        if (it_f->Fecund) {
          it_f->ClonalRepro();
          it_f->Matings += 1;}
      }
      for (auto&& it_m : Malesvec){
        if (it_m->Fecund) {
          it_m->ClonalRepro();
          it_m->Matings +=1;}
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


/*------------------------Add juveniles to the pop------------------------------------*/
    #ifdef TIMECHECK
    start = std::chrono::high_resolution_clock::now();
    #endif

    Femalesvec.insert(Femalesvec.end(),
    std::make_move_iterator(JuvFemalesvec.begin()),
    std::make_move_iterator(JuvFemalesvec.end()));

     Malesvec.insert(Malesvec.end(),
     std::make_move_iterator(JuvMalesvec.begin()),
     std::make_move_iterator(JuvMalesvec.end()));

    #ifdef TIMECHECK
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    LogFile << "Adding newborns to the population takes " << duration.count() << " microseconds" << endl;
    #endif


    //Clear the newborn vector
    JuvFemalesvec.clear();
    JuvMalesvec.clear();


/*------------------------Create output-------------------------------------------------*/
    //Exit when population is extinct
    if(Femalesvec.size()<=0 && Malesvec.size()<=0){
      std::cerr<<"The consumer population is extinct at time " << Time * delta_t << "\n";
      exit(1);
    }

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

  for (auto&& it_f : Femalesvec){
    print_individual(Endfile, *it_f);
    Endfile << endl;
  }

  for (auto&& it_m : Malesvec){
    print_individual(Endfile, *it_m);
    Endfile << endl;
  }

  LogFile << "\nSimulation ended at time " << Time * delta_t - delta_t << endl;
  cout << "Simulation ended at time " << Time * delta_t - delta_t << endl;

    Endfile.close();
    Traitfile.close();
    FullTraitfile.close();
    Timefile.close();
    LogFile.close();
    Matefile.close();
}

catch(InterruptException& e) //Try to write to output in case of interruption
      {
        LogFile  << "Caught signal " << e.S << endl;
        LogFile << "Writing current output to the .esf file" << endl;
        Endfile << "Time" << "\t" << "Volume" << "\t";

        for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
          print_resourceName(Endfile, *it_r);
          Endfile << "\t";
        }
        Endfile << endl;

        Endfile << Time*delta_t << "\t" << volume << "\t";

        for (it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r){
          print_resourceDensity(Endfile, *it_r);
          Endfile << "\t";
        }
        Endfile << endl;

        print_individualnames(Endfile);
        Endfile << "\n" << endl;
        for (auto&& it_f : Femalesvec){
          print_individual(Endfile, *it_f);
          Endfile << endl;
        }
        for (auto&& it_m : Malesvec){
          print_individual(Endfile, *it_m);
          Endfile << endl;
        }

        LogFile  << "Simulation ended at time " << Time * delta_t << endl;
        cout << "Simulation ended at time " << Time * delta_t << endl;

        Endfile.close();
        Traitfile.close();
        FullTraitfile.close();
        Timefile.close();
        LogFile.close();
        Matefile.close();

        return 1;
      }

    return 0;
  }
