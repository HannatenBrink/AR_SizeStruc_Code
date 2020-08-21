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


using namespace std;
stringstream ss;

//Max exponent value
double MAX_EXP = 50;
//When is something zero?
double eps = 10e-9;
double N_ini = 5; //initial number of individuals per sex & age group (total is 4*N_ini)

mt19937 mt_rand;
uniform_real_distribution<double> unif;
normal_distribution<double> MutNorm;
vector<int>::iterator int_it; //Iterator over integers
vector<double>::iterator int_doub; //Iterator over doubles
vector<Resource>::iterator it_r;  //Iterator over resources
int i, j;

double T_POP = 0; //For output
double T_Time = 0; //For output
double T_POP_FULL = 0; //for output
double T_Mate = 0; //for output

//Vectors for time-settings and parameters
std::vector<double> Setting(0);
std::vector<double> Parameter(0);
std::vector<double> RmaxChange; //Change in Rmaxvalues//
//vectors of food and individuals
vector<Individual> Females;
vector<Individual> Males;
vector<Individual> Females2;
vector<Individual> Males2;
vector<Individual> JuvFemales;
vector<Individual> JuvMales;
vector<Resource> JFood;
vector<Resource> AFood;

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
//decide the species identity//
int Nr_Res = 6;

std::vector<double> SpeciesDiv;

//Printing functions//
inline std::ostream& print_individual(std::ostream& os, const Individual& s){
  if(N_neutral && N_mating){
  os << s.age << "\t" << s.size << "\t" << s.sex << "\t" << s.SpeciesID << "\t"
  << s.mating_trait << "\t" << s.neutral_trait << "\t" << s.ecological_trait << "\t" <<
      s.Mature << "\t" <<
       s.Matings << "\t";
  }
  else if(N_mating){
     os << s.age << "\t" << s.size << "\t" << s.sex << "\t" << s.SpeciesID << "\t"
      <<  s.mating_trait << "\t" << s.ecological_trait << "\t" <<
      s.Mature << "\t" <<
         s.Matings << "\t";
} else {
  os << s.age << "\t" << s.size << "\t" << s.sex << "\t" << s.SpeciesID << "\t"
  <<   s.ecological_trait << "\t" <<
  s.Mature << "\t" <<
     s.Matings << "\t";
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
  os << "age" << "\t" << "size" << "\t" << "sex" << "\t" << "ID" << "\t"
  << "m_trait" << "\t" << "n_trait" << "\t" << "eco_trait" << "\t" <<
     "Mature" << "\t" <<
       "Matings" << "\t";
  } else if (N_mating){
      os << "age" << "\t" << "size" << "\t" << "sex" << "\t" << "ID" << "\t"
         << "m_trait" << "\t"  << "eco_trait" << "\t" <<  "Mature" << "\t" << "Matings" << "\t";
  }

  else {
  os << "age" << "\t" << "size" << "\t" << "sex" << "\t" << "ID" << "\t"
    <<  "eco_trait" << "\t" <<  "Mature" << "\t" <<"Matings" << "\t";
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
  os << s.age << "\t" << s.size << "\t" << s.sex << "\t" << s.SpeciesID << "\t"
  << s.mating_trait << "\t" << s.neutral_trait << "\t" << s.ecological_trait << "\t" <<
  s.Mature << "\t" <<
      s.Matings;
  } else if (N_mating) {
      os << s.age << "\t" << s.size << "\t" << s.sex << "\t" << s.SpeciesID << "\t"
      << s.mating_trait << "\t" << s.ecological_trait << "\t" << s.Mature << "\t" <<
        s.Matings;
  }
  else {
  os << s.age << "\t" << s.size << "\t" << s.sex << "\t" << s.SpeciesID << "\t"
  << s.ecological_trait << "\t" << s.Mature << "\t" <<
    s.Matings;
}
  return os;
}
inline std::ostream& print_traitsindividualnames(std::ostream& os){
  if(N_neutral && N_mating){
  os << "age" << "\t" << "size" << "\t" << "sex" << "\t" << "SpeciesID" << "\t"
   << "m_trait" << "\t" << "n_trait" << "\t"
  << "eco_trait "<< "\t" << "Mature" << "\t" << "Matings" << std::endl;
  } else if (N_mating){
      os << "age" << "\t" << "size" << "\t" << "sex" << "\t" << "SpeciesID" << "\t"
      << "m_trait" << "\t"
        << "eco_trait " << "\t" <<   "Mature" << "\t" << "Matings" << std::endl;
  }
  else {
  os << "age" << "\t" << "size" << "\t" << "sex" << "\t" << "SpeciesID" << "\t"
    << "eco_trait " << "\t" <<  "Mature" << "\t" << "Matings" << std::endl;
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

inline bool IsMarkedToStay(const Individual & o);

inline bool IsMarkedToStay(const Individual & o)
{
  return o.Is_dead == 0;
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

/*-----Init the env and pop in absence of isf file----------*/

void Init_Env() {
/*------------------------Initialize the environment-------------------------------*/
    Resource RB(RJmax, volume, ThetaB, TauB, "RB");
    Resource R1(RAmax, volume, Theta1, Tau1, "R1");
    Resource R2(RAmax, volume, Theta2, Tau2, "R2");
    Resource R3(RAmax, volume, Theta3, Tau3, "R3");
    Resource R4(RAmax, volume, Theta4, Tau4, "R4");
    Resource R5(RAmax, volume, Theta5, Tau5, "R5");
    Resource R6(RAmax, volume, Theta6, Tau6, "R6");
    JFood.push_back(RB);
    RmaxChange.push_back(0.0);
    AFood.push_back(R1);
    RmaxChange.push_back(0.0);
    AFood.push_back(R2);
    RmaxChange.push_back(0.0);
    AFood.push_back(R3);
    RmaxChange.push_back(0.0);
    AFood.push_back(R4);
    RmaxChange.push_back(0.0);
    AFood.push_back(R5);
    RmaxChange.push_back(0.0);
    AFood.push_back(R6);
    RmaxChange.push_back(0.0);
    std::cout << "Create resources" << std::endl;


/*------------------------Initialize the pop-------------*/
vector<int> mate_traits_f_ini;
vector<int> mate_traits_m_ini;
vector<double> neutral_traits_f_ini;
vector<double> neutral_traits_m_ini;
vector<double> ecological_traits_f_ini;
vector<double> ecological_traits_m_ini;

int p = 1;
if(N_mating && (ini_mate != 0 && ini_mate != 1 && ini_mate != -1)){
        std::cerr << "Initial average mating trait should be -1, 0, or 1" << std::endl;
        exit(1);
    }

if(ini_mate == 0) {
  for(i = 0; i < N_mating; ++i){
  p *= -1;
  mate_traits_f_ini.push_back(p);
  mate_traits_m_ini.push_back(-p);
}
} else {
  for(i = 0; i < N_mating; ++i){
  mate_traits_f_ini.push_back(ini_mate);
  mate_traits_m_ini.push_back(ini_mate);
}
}

for(i = 0; i < N_neutral; ++i){
  neutral_traits_f_ini.push_back(ini_neutral/(2 * N_neutral));
  neutral_traits_m_ini.push_back(ini_neutral/(2 * N_neutral));
}

for(i = 0; i < N_eco; ++i){
  ecological_traits_f_ini.push_back(ini_eco/(2 * N_eco));
  ecological_traits_m_ini.push_back(ini_eco/(2 * N_eco));
}

Individual MaleJuv(mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
ecological_traits_f_ini, ecological_traits_m_ini, 0, JFood, AFood);
Individual FemaleJuv(mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
ecological_traits_f_ini, ecological_traits_m_ini, 1, JFood, AFood);

Individual MaleAd(0, 100, 0, mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
ecological_traits_f_ini, ecological_traits_m_ini, JFood, AFood);
Individual FemaleAd(0, 100, 1, mate_traits_f_ini, mate_traits_m_ini, neutral_traits_f_ini, neutral_traits_m_ini,
ecological_traits_f_ini, ecological_traits_m_ini, JFood, AFood);

std::cout << "Add " << N_ini * 4 << " individuals to the init pop" << std::endl;
for(i = 0; i < N_ini; ++i){
  Females.push_back(FemaleJuv);
  Females.push_back(FemaleAd);
  Males.push_back(MaleJuv);
  Males.push_back(MaleAd);
}
}
