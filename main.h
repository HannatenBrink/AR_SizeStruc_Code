#include "Resource.h"
#include "Pars.h"
#include "Individual.h"

#include <set>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <stdio.h>
#include <cstdio>
#include <random>
#include <algorithm>
#include <exception>
#include <signal.h>
#include <functional>


using namespace std;
stringstream ss;

//Max exponent value
double MAX_EXP = 50;
//When is something zero?
double eps = pow(10,-9);
double N_ini = 5; //initial number of individuals per age group (total is 4*N_ini)
int Nr_Res = 6;

double Time;
double phi; //Ontogenetic shift
double psi;
double Total;
double IntakeTot;
double Value;
int Offspring;
double multi;
double dif;
double RandomVal;
int TotStarv = 0;
double OverCons = 0;

double Pop_Density = Femalevec.size() + Malevec.size();
double Pop_Mass = 0;
double Females = 0;
double Males = 0;

double Pop_Density1 = 0;
double Pop_Mass1 = 0;
double Females1 = 0;
double Males1 = 0;

double Pop_Density2 = 0;
double Pop_Mass2 = 0;
double Females2 = 0;
double Males2 = 0;

double Pop_Density3 = 0;
double Pop_Mass3 = 0;
double Females3 = 0;
double Males3 = 0;

double Pop_Density4 = 0;
double Pop_Mass4 = 0;
double Females4 = 0;
double Males4 = 0;

double Pop_Density5 = 0;
double Pop_Mass5 = 0;
double Females5 = 0;
double Males5 = 0;

double Pop_Density6 = 0;
double Pop_Mass6 = 0;
double Females6 = 0;
double Males6 = 0;

double Pop_Density7 = 0;
double Pop_Mass7 = 0;
double Females7 = 0;
double Males7 = 0;



uniform_real_distribution<double> unif;
uniform_int_distribution<> IDNR_Rand(1, 100000);
normal_distribution<double> MutNorm;
exponential_distribution<double> Surv_age;
vector<int>::iterator int_it; //Iterator over integers
vector<double>::iterator int_doub; //Iterator over doubles
vector<Resource>::iterator it_r;  //Iterator over resources
int i, j;

double T_POP = 0; //For output
double T_Time = 0; //For output
double T_POP_FULL = 0; //for output
double T_Mate = 0; //for output of mate file

ofstream LogFile;
ofstream Matefile;
ofstream Timefile;
ofstream Endfile;
ofstream FullTraitfile;

//Vectors for time-settings and parameters
std::vector<double> Setting(0);
std::vector<double> Parameter(0);
std::vector<double> RmaxChange; //Change in Rmaxvalues//
//vectors of food and individuals
vector<unique_ptr<Individual>> Femalevec;
vector<unique_ptr<Individual>> Malevec;
vector<unique_ptr<Individual>> Juvvec;
vector<Resource> AllFood;



//Function to split a line. Needed to read the cvf file
void split(const std::string &s, char delim, std::vector<std::string> &elems){
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);}
  }


//Function to get a random number from a vector of cumulative probabilities//
inline int weighted_random_known_sums_floats(std::vector<double> &cumulative_rates, double &r, int num_elements) {
if(num_elements == 1) return 0;

if(num_elements == 2) {
  if(cumulative_rates[0] < r) return 1;
  else return 0;
}

int lowGuess = 0;
int highGuess = num_elements - 1;
int guess;

while(highGuess >= lowGuess) {
  guess = (lowGuess + highGuess) / 2;
  if(cumulative_rates[guess] < r) lowGuess = guess + 1;
  else if(guess > 0 && cumulative_rates[guess - 1] > r) highGuess = guess - 1;
  else return guess;
}
return guess;
}

//Species identities///
std::vector<double> SpeciesDiv;

//Printing functions//
inline std::ostream& print_individual(std::ostream& os, const Individual& s){
  if(N_neutral && N_mating){
  os << s.IDNR << '\t' << s.sex << "\t" << s.age << "\t" << s.size << "\t"  << s.SpeciesID << "\t"
  << s.mating_trait << "\t" << s.neutral_trait << "\t" << s.ecological_trait << "\t" <<
      s.Mature << "\t" << s.Matings << "\t" << s.repro_buffer << "\t" << s.MaxAge << "\t" << s.Starve << "\t";
  }
  else if(N_mating){
     os << s.IDNR << '\t' << s.sex << '\t' << s.age << "\t" << s.size << "\t"  << s.SpeciesID << "\t"
      <<  s.mating_trait << "\t" << s.ecological_trait << "\t" << s.Mature << "\t" << s.Matings << "\t"
      << s.repro_buffer << "\t" << s.MaxAge << "\t" << s.Starve << "\t";
} else {
    os << s.IDNR << '\t' << s.sex << '\t' << s.age << "\t" << s.size << "\t"  << s.SpeciesID << "\t"
       << s.ecological_trait << "\t" << s.Mature << "\t"
       << s.Matings << "\t" << s.repro_buffer << "\t" << s.MaxAge << "\t" << s.Starve << "\t";
  }
  for(i = 0; i < N_mating; ++i){
    os << s.mating_trait_alleles_f[i] << "\t";
  }
  for(i = 0; i < N_mating; ++i){
    os << s.mating_trait_alleles_m[i] << "\t";
  }
  for(i = 0; i < N_neutral; ++i){
    os << s.neutral_trait_alleles_f[i] << "\t";
  }
  for(i = 0; i < N_neutral; ++i){
    os << s.neutral_trait_alleles_m[i] << "\t";
  }
  for(i = 0; i < N_eco; ++i){
    os << s.ecological_trait_alleles_f[i] << "\t";
  }
  for(i = 0; i < N_eco; ++i){
    os << s.ecological_trait_alleles_m[i] << "\t";
  }
  return os;
}
inline std::ostream& print_individualnames(std::ostream& os) {
  if(N_neutral && N_mating){
    os << "ID_NR" << '\t' << "sex" << '\t' << "age" << "\t" << "size" << "\t"  << "ID" << "\t"
      << "m_trait" << "\t" << "n_trait" << "\t" << "eco_trait" << "\t"
      << "Mature" << "\t" <<"Matings" << "\t" << "ReproBuffer" << "\t" << "MaxAge" << "\t" << "Starve" << '\t';
  } else if (N_mating){
      os << "ID_NR" << '\t' << "sex" << '\t' << "age" << "\t" << "size" << "\t"  << "ID" << "\t"
         << "m_trait" << "\t"  << "eco_trait" << "\t" <<  "Mature" << "\t" << "Matings" << "\t"
         << "ReproBuffer" << "\t" << "MaxAge" << "\t" << "Starve" << '\t';
  }
  else {
    os << "ID_NR" << '\t' << "sex" << '\t' << "age" << "\t" << "size" << "\t"  << "ID" << "\t"
      <<  "eco_trait" << "\t" <<  "Mature" << "\t" <<"Matings" << "\t" << "ReproBuffer"
      << "\t" <<" MaxAge" << "\t" << "Starve" << '\t';
}
  for(i = 0; i < N_mating; ++i){
    os << "m_A_f" << i << "\t";
  }
  for(i = 0; i < N_mating; ++i){
    os << "m_A_m" << i << "\t";
  }
  for(i = 0; i < N_neutral; ++i){
    os << "n_A_f" << i << "\t";
  }
  for(i = 0; i < N_neutral; ++i){
    os << "n_A_m" << i << "\t";
  }
  for(i = 0; i < N_eco; ++i){
    os << "eco_A_f" << i << "\t";
  }
  for(i = 0; i < N_eco; ++i){
    os << "eco_A_m" << i << "\t";
  }
  return os;
}
inline std::ostream& print_traitsindividual(std::ostream& os, const Individual& s){
  if(N_neutral && N_mating){
    os << s.IDNR << '\t' << s.sex << '\t' << s.age << "\t" << s.size << "\t"  << s.SpeciesID << "\t"
      << s.mating_trait << "\t" << s.neutral_trait << "\t" << s.ecological_trait << "\t"
      << s.Mature << "\t" << s.Matings << "\t" << s.repro_buffer << '\t' << s.MaxAge << '\t' << s.Starve << '\t';
  } else if (N_mating) {
      os << s.IDNR << '\t' << s.sex << '\t' << s.age << "\t" << s.size << "\t"  << s.SpeciesID << "\t"
        << s.mating_trait << "\t" << s.ecological_trait << "\t" << s.Mature << "\t" <<
        s.Matings << "\t" << s.repro_buffer << '\t' << s.MaxAge << '\t' << s.Starve << '\t';
  }
  else {
    os << s.IDNR << '\t' << s.sex << '\t' << s.age << "\t" << s.size << "\t"  << s.SpeciesID << "\t"
    << s.ecological_trait << "\t" << s.Mature << "\t"
    << s.Matings << "\t" << s.repro_buffer << '\t' << s.MaxAge << '\t' << s.Starve; //REMOVE Starve
}
  return os;
}
inline std::ostream& print_traitsindividualnames(std::ostream& os){
  if(N_neutral && N_mating){
  os << "ID_NR" << '\t' << "sex" << '\t' << "age" << "\t" << "size" << "\t"  << "SpeciesID" << "\t"
   << "m_trait" << "\t" << "n_trait" << "\t"
  << "eco_trait "<< "\t" << "Mature" << "\t" << "Matings" << "\t" << "ReproBuffer" << "\t" << "MaxAge" << "\t" << "Starve" << std::endl;
  } else if (N_mating){
      os << "ID_NR" << '\t' << "sex" << '\t' << "age" << "\t" << "size" << "\t"  << "SpeciesID" << "\t"
      << "m_trait" << "\t"
        << "eco_trait " << "\t" <<   "Mature" << "\t" << "Matings" << "\t" << "ReproBuffer" << '\t' << "MaxAge" << "\t" << "Starve" << std::endl;
  }
  else {
  os << "ID_NR" << '\t' << "sex" << '\t' << "age" << "\t" << "size" << "\t"  << "SpeciesID" << "\t"
    << "eco_trait " << "\t" <<  "Mature" << "\t" << "Matings" << "\t" << "ReproBuffer" << '\t' << "MaxAge" << '\t' << "Starve" << std::endl;
}
  return os;
}
inline std::ostream& print_resource(std::ostream& os, const Resource& s){
    os << "The density of resource " << s.Name <<"equals "<< s.Density << " gram per m3. The volume of the system is " << s.Volume
    << " m3. Therefore, the total amount of this resource is " <<s.Density * s.Volume * 0.001 << " gram." << std::endl;
    return os;
}
inline std::ostream& print_resourceDensity(std::ostream& os, const Resource& s){
    os << s.Density;
    return os;
}
inline std::ostream& print_resourceDensityTotal(std::ostream& os, const Resource& s){
    os << s.Density * s.Volume;
    return os;
}
inline std::ostream& print_resourceName(std::ostream& os, const Resource& s){
    os << s.Name;
    return os;
}
inline std::ostream& print_resourceInfo(std::ostream& os, const Resource& s){
    os << s.Name << "\t" << s.Density << "\t" << s.OptTrait << "\t" << s.Rmax << "\t" << std::endl;
    return os;
}
inline std::ostream& print_resourceInfo_name(std::ostream& os){
    os << "Name" << "\t" << "Density" << "\t" << "OptTrait"<< "\t" << "Rmax" << "\t" << std::endl;
    return os;
}

inline bool IsMarkedToDelete(const Individual & o);

inline bool IsMarkedToDelete(const Individual & o)
{
  return o.Is_dead > 0;
}

inline bool IsMarkedToMove(const Individual & o);

inline bool IsMarkedToMove(const Individual & o)
{
  return o.Mature > 0;
}

class InterruptException : public std::exception
{
public:
  InterruptException(int s) : S(s) {}
  int S;
};

void sig_to_exception(int s)
{
  throw InterruptException(s);
}

/*--------------Function to close program when receiving a signal------------------*/

void signalHandler( int signum ){
   cout << "Interrupt signal (" << signum << ") received.\n";
   cout << "Stop IBMrun and write output to the .esf file" << endl;

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
   if(Malevec.size() > 0){for (auto&& it_f = Malevec.begin(); it_f != Malevec.end(); ++it_f){
     print_individual(Endfile, **it_f);
     Endfile << endl;
   }}
   if(Femalevec.size() > 0){for (auto&& it_f = Femalevec.begin(); it_f != Femalevec.end(); ++it_f){
     print_individual(Endfile, **it_f);
     Endfile << endl;
   }}

   if(Juvvec.size() > 0){
   for (auto&& it_f = Juvvec.begin(); it_f != Juvvec.end(); ++it_f){
     print_individual(Endfile, **it_f);
     Endfile << endl;
   }}



   LogFile << "\nSimulation ended at time " << Time * delta_t - delta_t << endl;
   cout << "Simulation ended at time " << Time * delta_t - delta_t << endl;

 //close files
     Endfile.close();
     FullTraitfile.close();
     Timefile.close();
     LogFile.close();
     Matefile.close();
   exit(signum);
}




/*-----Init the env and pop in absence of isf file----------*/

void Init_Env(std::mt19937& mt_rand_ref) {
/*------------------------Initialize the environment-------------------------------*/
    Resource RB(RJmax, volume, ThetaB, TauB, "RB");
    Resource R1(RAmax, volume, Theta1, Tau1, "R1");
    Resource R2(RAmax, volume, Theta2, Tau2, "R2");
    Resource R3(RAmax, volume, Theta3, Tau3, "R3");
    Resource R4(RAmax, volume, Theta4, Tau4, "R4");
    Resource R5(RAmax, volume, Theta5, Tau5, "R5");
    Resource R6(RAmax, volume, Theta6, Tau6, "R6");
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
    std::cout << "Create resources" << std::endl;


/*------------------------Initialize the pop-------------*/
vector<double> mate_traits_f_ini;
vector<double> mate_traits_m_ini;
vector<double> neutral_traits_f_ini;
vector<double> neutral_traits_m_ini;
vector<double> ecological_traits_f_ini;
vector<double> ecological_traits_m_ini;


for(i = 0; i < N_mating; ++i){
  mate_traits_f_ini.push_back(ini_mate);
  mate_traits_m_ini.push_back(ini_mate);
}

for(i = 0; i < N_neutral; ++i){
  neutral_traits_f_ini.push_back(ini_neutral/(2 * N_neutral));
  neutral_traits_m_ini.push_back(ini_neutral/(2 * N_neutral));
}

for(i = 0; i < N_eco; ++i){
  ecological_traits_f_ini.push_back(ini_eco/(2 * N_eco));
  ecological_traits_m_ini.push_back(ini_eco/(2 * N_eco));
}

std::cout << "Add " << N_ini * 4 << " individuals to the init pop" << std::endl;
for(int k = 0; k < N_ini; ++k){

  //juveniles//
  unique_ptr<Individual> IndivPtr1(new Individual(mt_rand_ref, (k+1),mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
  ecological_traits_f_ini, ecological_traits_m_ini,  AllFood));
  unique_ptr<Individual> IndivPtr2(new Individual(mt_rand_ref, (k+1)+4,mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
  ecological_traits_f_ini, ecological_traits_m_ini,  AllFood));

  //adults (IDnr, sex, age, size, mature, matings, reprodbuf, maxage, starve_nr)
  unique_ptr<Individual> IndivPtr3(new Individual((k+1)+8, 1, 0, 100,  1, 0, 0, 100, 0, mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
  ecological_traits_f_ini, ecological_traits_m_ini,  AllFood));
  unique_ptr<Individual> IndivPtr4(new Individual((k+1)+12, 0, 0, 100,  1, 0, 0, 100, 0, mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
  ecological_traits_f_ini, ecological_traits_m_ini,  AllFood));

  Juvvec.push_back(move(IndivPtr1));
  Juvvec.push_back(move(IndivPtr2));
  Femalevec.push_back(move(IndivPtr3));
  Malevec.push_back(move(IndivPtr4));
}
}

bool cmp_by_repro(const std::unique_ptr<Individual> &a, const std::unique_ptr<Individual> &b)
{
    return a->repro_buffer > b->repro_buffer;
}
