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
  "mu_b" , "Xi" , "M_B" , "M_Mat" , "M_Shift" , "s_ass" , "s_diss", "N_eco" , "N_neutral" , "N_mating" ,
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
    }
     } else {
    LogFile << "Reproduction is sexual" << endl;
    if(N_mating==0) {
      LogFile << "Error: Sexual reproduction is only possible in case N_mating > 0. Exit run..." << endl;
      cout << "Error: Sexual reproduction is only possible in case N_mating > 0. Exit run..." << endl;
      exit(1);
    }
  }
  LogFile << endl << "________________________________________" << endl << endl;

 /*--------------------------------Define variables-------------------------------------------*/
   MutNorm.param(std::normal_distribution<double>(0, mut_std).param());
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
    Init_Env();
  } else {
    LogFile << "Initialization occurs via the isf file\n";
    //The new runs will have as first element the ID of the individual, and in addition reprobuffer
    //Has to change in the code
    Time = Starttime_init;
    double RB_init, R1_init, R2_init, R3_init, R4_init, R5_init, R6_init, Volume_init;
    double age_init, size_init, mating_trait_init, eco_trait_init, neutral_trait_init, MaxAge_init;
    double dummy, idnr, reprobuf, nrStarve;
    int id, Mature, Matings;
    int i = 0;
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
      if (i > 3) {
        //vector <int> mating_f;
        vector <double> mating_f;
        vector <double> eco_f;
        vector <double> neutral_f;
        //vector <int> mating_m;
        vector <double> mating_m;
        vector <double> eco_m;
        vector <double> neutral_m;
        if(N_neutral && N_mating){
          linestream >> idnr >> age_init >> size_init >>  id >> mating_trait_init >>
              neutral_trait_init >> eco_trait_init >> Mature >> Matings >> reprobuf >> MaxAge_init;
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
            unique_ptr<Individual> IndivPtr(new Individual(idnr, age_init, size_init, Matings, reprobuf, MaxAge_init, Mature,  mating_f, mating_m, neutral_f, neutral_m,
              eco_m, eco_f, AllFood));
              if(Mature) {
                Advec.push_back(move(IndivPtr));
              } else {
                Juvvec.push_back(move(IndivPtr));
              }
         }
        else if(N_mating){
          //linestream >>  age_init >> size_init >> sex_init >> id >> mating_trait_init >>
            //  eco_trait_init >> Mature >> Matings >> MaxAge_init;
            linestream >>  idnr >> age_init >> size_init >> id >> mating_trait_init >>
                eco_trait_init >> Mature >> Matings >> reprobuf >> MaxAge_init;
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
          unique_ptr<Individual> IndivPtr(new Individual(idnr, age_init, size_init, Matings, reprobuf, MaxAge_init, Mature, mating_f, mating_m,
                eco_m, eco_f,  AllFood)); //already changed.
          if(Mature) {
            Advec.push_back(move(IndivPtr));
          } else {
            Juvvec.push_back(move(IndivPtr));
          }
         }
        else {
           linestream >> idnr >> age_init >> size_init >>  id >>
               eco_trait_init >> Mature >> Matings >> reprobuf >> MaxAge_init >> nrStarve;
           for(int j = 0; j < N_eco; ++j) {
                   linestream >> dummy;
                   eco_f.push_back(dummy);
               }
           for(int j = 0; j < N_eco; ++j) {
                   linestream >> dummy;
                   eco_m.push_back(dummy);
               }
               unique_ptr<Individual> IndivPtr(new Individual(idnr, age_init, size_init, Matings, reprobuf, MaxAge_init, Mature, eco_f, eco_m,   AllFood));
               if(Mature) {
                 Advec.push_back(move(IndivPtr));
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

  LogFile << "The initial population size is " << Advec.size()+Juvvec.size()
    << endl;
  LogFile << "...End of initializing the population" << endl;
  LogFile << endl << "________________________________________" << endl << endl;


/*-----------------------------------Give names to the files-----------------------------------*/

  LogFile << "Add names of variables to the files, and initialize variables..." << endl;

  Timefile << "Time" << '\t' << "Starve" << '\t' << "OverCons" << '\t' <<
   "PopSize" << '\t' << "PopMass" << '\t' << "Adults" << '\t' <<
   "PopSize1" << '\t' << "PopMass1" << '\t' << "Adults1" << '\t' <<
   "PopSize2" << '\t' << "PopMass2" << '\t' << "Adults2" << '\t' <<
   "PopSize3" << '\t' << "PopMass3" << '\t' << "Adults3" << '\t' <<
   "PopSize4" << '\t' << "PopMass4" << '\t' << "Adults4" << '\t' <<
   "PopSize5" << '\t' << "PopMass5" << '\t' << "Adults5" << '\t' <<
   "PopSize6" << '\t' << "PopMass6" << '\t' << "Adults6" << '\t' <<
   "PopSize7" << '\t' << "PopMass7" << '\t' << "Adults7" << '\t';
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
    //Time = Starttime_init;

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
      for (auto&& it_f : Juvvec){
        Traitfile << Time * delta_t << "\t";
        print_traitsindividual(Traitfile, *it_f);
        Traitfile << endl;
      }
      for (auto&& it_f : Advec){
        Traitfile << Time * delta_t << "\t";
        print_traitsindividual(Traitfile, *it_f);
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
        for (auto&& it_f : Advec){
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
          Pop_Density = Advec.size() + Juvvec.size();
          Pop_Mass = 0;
          Adults = 0;

          Pop_Density1 = 0;
          Pop_Mass1 = 0;
          Adults1 = 0;

          Pop_Density2 = 0;
          Pop_Mass2 = 0;
          Adults2 = 0;

          Pop_Density3 = 0;
          Pop_Mass3 = 0;
          Adults3 = 0;

          Pop_Density4 = 0;
          Pop_Mass4 = 0;
          Adults4 = 0;

          Pop_Density5 = 0;
          Pop_Mass5 = 0;
          Adults5 = 0;

          Pop_Density6 = 0;
          Pop_Mass6 = 0;
          Adults6 = 0;

          Pop_Density7 = 0;
          Pop_Mass7 = 0;
          Adults7 = 0;

        for (auto&& it_f : Advec){
          Pop_Mass += it_f->size;
          Adults += 1;
          switch(it_f->SpeciesID) {
            case 1:
            Pop_Mass1 += it_f->size;
            Pop_Density1 += 1;
            Adults1 += 1;
            break;
            case 2:
            Pop_Mass2 += it_f->size;
            Pop_Density2 += 1;
            Adults2 += 1;
            break;
            case 3:
            Pop_Mass3 += it_f->size;
            Pop_Density3 += 1;
            Adults3 += 1;
            break;
            case 4:
            Pop_Mass4 += it_f->size;
            Pop_Density4 += 1;
            Adults4 += 1;
            break;
            case 5:
            Pop_Mass5 += it_f->size;
            Pop_Density5 += 1;
            Adults5 += 1;
            break;
            case 6:
            Pop_Mass6 += it_f->size;
            Pop_Density6 += 1;
            Adults6 += 1;
            break;
            case 7:
            Pop_Mass7 += it_f->size;
            Pop_Density7 += 1;
            Adults7 += 1;
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
        << Pop_Density << '\t'  << Pop_Mass  << '\t'  << Adults <<'\t'
        << Pop_Density1 << '\t' << Pop_Mass1 << '\t'  << Adults1 <<'\t'
        << Pop_Density2 << '\t' << Pop_Mass2 << '\t'  << Adults2 <<'\t'
        << Pop_Density3 << '\t' << Pop_Mass3 << '\t'  << Adults3 <<'\t'
        << Pop_Density4 << '\t' << Pop_Mass4 << '\t'  << Adults4 <<'\t'
        << Pop_Density5 << '\t' << Pop_Mass5 << '\t'  << Adults5 <<'\t'
        << Pop_Density6 << '\t' << Pop_Mass6 << '\t'  << Adults6 <<'\t'
        << Pop_Density7 << '\t' << Pop_Mass7 << '\t'  << Adults7 <<'\t' ;

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
    for (auto&& it_f : Advec){
      it_f->R_Intake(AllFood, Feeding);
    }

    //cout << "Before growth Size of the juvvec is " << Juvvec.size() << " Adult size is " << Advec.size() << endl;
    for (auto&& it_f : Juvvec){
      it_f->R_Intake(AllFood, Feeding);
      if(it_f->Mature){
        Advec.push_back(std::move(it_f)); //move mature individuals to adults
      }
    }
    //cout << "After growth, Size of the juvvec is " << Juvvec.size() << " Adult size is " << Advec.size() << endl;

    //cout << "After removing, Size of the juvvec is " << Juvvec.size() << " Adult size is " << Advec.size() << endl;



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

    for(i = 0, it_r = AllFood.begin(); it_r != AllFood.end(); ++it_r, ++i){
        //print_resourceName(cout, *it_r);
        //cout << endl << "Before eating the total density is: ";
        //print_resourceDensity(cout, *it_r);
        it_r->Growth(Feeding[i], RmaxChange[i]);
        //cout << endl;
        //cout << "Resource eaten is " << Feeding[i]/volume << endl << "After feeding the total density is: ";
        //print_resourceDensity(cout, *it_r);
        //cout << endl;
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

    Advec.erase(std::remove_if(Advec.begin(), Advec.end(), [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Advec.end());
    Juvvec.erase(std::remove_if(Juvvec.begin(), Juvvec.end(), [](std::unique_ptr<Individual> &tst) { return tst->Is_dead;}), Juvvec.end());



    #ifdef TIMECHECK
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    LogFile << "Removing death animals takes " << duration.count() << " microseconds" << endl;
    LogFile << "Size of the new population is " << Advec.size() << endl;
    if ((Output_time > 0) &&
       ((round(fmod(T_Time, (Output_time / delta_t))) <= 0) ||
       (round(fmod(T_Time, (Output_time / delta_t)) - (Output_time / delta_t))  == 0))) {
    Timefile << duration.count() << '\t';}
    #endif

    //Exit when population is extinct
    if(Advec.size()<=0 && Juvvec.size()<=0){
      std::cerr<<"The consumer population is extinct at time " << Time * delta_t << "\n";
      exit(1);
    }

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
      //shuffle(Advec.begin(), Advec.end(), mt_rand); //necessary because females that choose first have higher prob to find a fecund mate
      //maybe better to shuffle it depending on reproductive buffer? In this way maybe less mate limitation?
      //downside is that a female that can't find a mate because she has a rare trait, might ultimately have a rare mating opportunity.
      sort(Advec.begin(), Advec.end(), cmp_by_repro); //sort
    for (auto&& it_f : Advec) {
      if (it_f->Fecund) {
      vector<double> cumsum;
      Tot = 0;
      //cout << "Female eco trait: " << it_f->ecological_trait << " mating trait: " << it_f->mating_trait << endl;
      #ifdef TIMECHECK
      start2 = std::chrono::high_resolution_clock::now();
      #endif
      //int test = 0;
      for (auto&& it_m : Advec){
        if(it_m == it_f){
         it_m->matingProb = 0;
         cumsum.push_back(Tot);
         //cout << "Male " << test << " fecund: " << it_m->Fecund << " ecotrait female: " << it_f->ecological_trait << " ecotrait: " << it_m->ecological_trait << " its mating prob: " << it_m->matingProb << endl;
        }
        else {
        it_m->MateProb(*it_f, cumsum, Tot); //calculate probabilities
        //cout << "Male " << test << " fecund: " << it_m->Fecund << " ecotrait female: " << it_f->ecological_trait << " ecotrait: " << it_m->ecological_trait << " its mating prob: " << it_m->matingProb << endl;
        }
        //test += 1;
      }
      #ifdef TIMECHECK
      stop2 = std::chrono::high_resolution_clock::now();
      dur2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
      Tot2 += dur2.count();
      #endif
      RandomVal = unif(mt_rand) * Tot;
      mate = weighted_random_known_sums_floats(cumsum, RandomVal, cumsum.size()); //roulette wheel selection with binary search
      auto it_m = next(Advec.begin(), mate);
      //cout << "Ecotrait female: " << it_f->ecological_trait << " Chosen Male ecotrait: " << (*it_m)->ecological_trait << " its mating prob: " << (*it_m)->matingProb << endl;
      if(((*it_m)->Fecund) & (Tot > 0)) {
      it_f->Matings += 1;
      (*it_m)->Matings += 1;
      it_f->SexualRepro(**it_m); //reproduce
      //Write mate choice to matefile
      if ((MateFile > 0) && ((round(fmod(T_Mate, (MateFile / delta_t))) == 0)|| (round(fmod(T_Mate, (MateFile / delta_t)) - (MateFile / delta_t))  == 0))) {
          Matefile << Time * delta_t << "\t" << it_f->age << "\t" << it_f->size << "\t" <<
          it_f->ecological_trait << "\t"
          << it_f->neutral_trait << "\t" <<
          it_f->mating_trait << "\t" <<
          it_f->repro_buffer << "\t" <<
          it_f->SpeciesID << "\t" <<
          it_f->IDNR << "\t" <<
          (*it_m)->age << "\t" << (*it_m)->size << "\t" <<
          (*it_m)->ecological_trait << "\t" <<
          (*it_m)->neutral_trait << "\t" <<
          (*it_m)->mating_trait << "\t" <<
          (*it_m)->repro_buffer << "\t" <<
          (*it_m)->SpeciesID << "\t" <<
          (*it_m)->IDNR << "\t" <<
          endl;
      }
    } else if ((MateFile > 0) && ((round(fmod(T_Mate, (MateFile / delta_t))) == 0)|| (round(fmod(T_Mate, (MateFile / delta_t)) - (MateFile / delta_t))  == 0))) {

          Matefile << Time * delta_t << "\t" << it_f->age << "\t" << it_f->size << "\t" <<
          it_f->ecological_trait << "\t"
          << it_f->neutral_trait << "\t" <<
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
          0 << endl;
      }
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
      for (auto&& it_f : Advec){
        if (it_f->Fecund) {
          it_f->ClonalRepro();
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


/*------------------------Add juveniles to the pop------------------------------------*/
    /*#ifdef TIMECHECK
    start = std::chrono::high_resolution_clock::now();
    #endif*/

  /*  Advec.insert(Advec.end(),
    std::make_move_iterator(Juvvec.begin()),
    std::make_move_iterator(Juvvec.end()));*/


    /*#ifdef TIMECHECK
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    LogFile << "Adding newborns to the population takes " << duration.count() << " microseconds" << endl;
    #endif*/


    //Clear the newborn vector
    //Juvvec.clear();



/*------------------------Create output-------------------------------------------------*/

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

  for (auto&& it_f : Advec){
    print_individual(Endfile, *it_f);
    Endfile << endl;
  }
  for (auto&& it_f : Juvvec){
    print_individual(Endfile, *it_f);
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
        for (auto&& it_f : Advec){
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
        Traitfile.close();
        FullTraitfile.close();
        Timefile.close();
        LogFile.close();
        Matefile.close();

        return 1;
      }
    return 0;
  }
