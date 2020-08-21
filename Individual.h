#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <numeric>
#include <random>
#include <math.h>
#include "Pars.h"
#include "Resource.h"

extern int i;
extern int j;
extern std::vector<int>::iterator int_it;
extern std::vector<double>::iterator int_doub;
extern std::mt19937 mt_rand;
extern std::uniform_real_distribution<double> unif;
extern std::normal_distribution<double> MutNorm;
extern std::vector<Resource>::iterator it_r;
extern std::ostream& print_individualnames(std::ostream&);
extern std::ostream& print_traitsindividualnames(std::ostream&);
extern std::vector<double> SpeciesDiv;
extern int Nr_Res;
extern std::string arrname[];

class Individual {

/*------------------Friends---------------------------------------*/
  friend std::ostream& print_individual(std::ostream& os, const Individual& s);
  friend std::ostream& print_traitsindividual(std::ostream& os, const Individual& s);

public:

  /*------------------Constructor---------------------------------------*/

  //Old individuals no mut//
  Individual(double age, double size, const int sex,
    std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
    std::vector<double> neutral_traits_f, std::vector<double> neutral_traits_m,
    std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
    std::vector<Resource>& JFood, std::vector<Resource>& AFood)
    : size(size), age(age), sex(sex),
    mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
    neutral_trait_alleles_f(neutral_traits_f), neutral_trait_alleles_m(neutral_traits_m),
    ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
     {
      repro_buffer = 0;
      Starve = false;
      Fecund = false;
      Is_dead = false;
      Mature = false;
      ecological_trait = 0;
      if (N_neutral) {
        neutral_trait = std::accumulate(neutral_trait_alleles_f.begin(), neutral_trait_alleles_f.end(), 0.0);
        neutral_trait += std::accumulate(neutral_trait_alleles_m.begin(), neutral_trait_alleles_m.end(), 0.0);}
      mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0)/(2*N_mating);
      mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0)/(2*N_mating);
      ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0);
      ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0);
      //std::cout << "Traits (neutral mate eco) " << neutral_trait << " " << mating_trait << " " << ecological_trait << std::endl;


      //Determine species identity//
      for (i = 0; i < Nr_Res - 1; ++i){
        SpeciesID = 1;
        if (ecological_trait >= SpeciesDiv[i]){
          SpeciesID = Nr_Res - i;
          break;
        }
      }
      //These are correct//
      for(i = 0, it_r = JFood.begin(); it_r < JFood.end(); ++it_r, ++i) {
          double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait, 2)/(2*it_r->Tau*it_r->Tau));
          attack_constants.push_back(Const);
      }
      for(it_r = AFood.begin(); it_r < AFood.end(); ++it_r, ++i) {
        double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
        attack_constants.push_back(Const);
       }
    }

    Individual(double age, double size, const int sex,
      std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
      std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
      std::vector<Resource>& JFood, std::vector<Resource>& AFood)
      : size(size), age(age), sex(sex),
      mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
      ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
       {
        repro_buffer = 0;
        Starve = false;
        Fecund = false;
        Is_dead = false;
        Mature = false;
        ecological_trait = 0;
        mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0)/(2*N_mating);
        mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0)/(2*N_mating);
        ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0);
        ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0);


        //Determine species identity//
        for (i = 0; i < Nr_Res - 1; ++i){
          SpeciesID = 1;
          if (ecological_trait >= SpeciesDiv[i]){
            SpeciesID = Nr_Res - i;
            break;
          }
        }
        //These are correct//
        for(i = 0, it_r = JFood.begin(); it_r < JFood.end(); ++it_r, ++i) {
            double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait, 2)/(2*it_r->Tau*it_r->Tau));
            attack_constants.push_back(Const);
        }
        for(it_r = AFood.begin(); it_r < AFood.end(); ++it_r, ++i) {
          double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
          attack_constants.push_back(Const);
         }
      }

    Individual(double age, double size, const int sex,
      std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
      std::vector<Resource>& JFood, std::vector<Resource>& AFood)
      : size(size), age(age), sex(sex),
      ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
       {
        repro_buffer = 0;
        Starve = false;
        Fecund = false;
        Is_dead = false;
        Mature = false;
        ecological_trait = 0;
       ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0);
        ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0);

        //Determine species identity//
        for (i = 0; i < Nr_Res - 1; ++i){
          SpeciesID = 1;
          if (ecological_trait >= SpeciesDiv[i]){
            SpeciesID = Nr_Res - i;
            break;
          }
        }
        //These are correct//
        for(i = 0, it_r = JFood.begin(); it_r < JFood.end(); ++it_r, ++i) {
            double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait, 2)/(2*it_r->Tau*it_r->Tau));
            attack_constants.push_back(Const);
        }
        for(it_r = AFood.begin(); it_r < AFood.end(); ++it_r, ++i) {
          double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
          attack_constants.push_back(Const);
         }
      }


  //Newborn no mut//
  Individual(std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
    std::vector<double> neutral_traits_f, std::vector<double> neutral_traits_m,
    std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
    const int sex, std::vector<Resource>& JFood, std::vector<Resource>& AFood)
    : sex(sex), mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
    neutral_trait_alleles_f(neutral_traits_f), neutral_trait_alleles_m(neutral_traits_m),
    ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
     {
      age = 0;
      size = size_birth;
      repro_buffer = 0;
      Starve = false;
      Fecund = false;
      Is_dead = false;
      Mature = false;
      if (N_neutral) {
        neutral_trait = std::accumulate(neutral_trait_alleles_f.begin(), neutral_trait_alleles_f.end(), 0.0);
        neutral_trait += std::accumulate(neutral_trait_alleles_m.begin(), neutral_trait_alleles_m.end(), 0.0);}
      mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0)/(2*N_mating);
      mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0)/(2*N_mating);
      ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0);
      ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0);

      //Determine species identity//
      for (i = 0; i < Nr_Res - 1; ++i){
        SpeciesID = 1;
        if (ecological_trait >= SpeciesDiv[i]){
          SpeciesID = Nr_Res - i;
          break;
        }
      }

      for(i = 0, it_r = JFood.begin(); it_r < JFood.end(); ++it_r, ++i) {
          double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
          attack_constants.push_back(Const);
      }
      for(it_r = AFood.begin(); it_r < AFood.end(); ++it_r, ++i) {
        double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
        attack_constants.push_back(Const);
       }
    }

  //Newborn with mut//
  Individual(std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
    std::vector<double> neutral_traits_f, std::vector<double> neutral_traits_m,
    std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
    const int sex, std::vector<Resource>& JFood, std::vector<Resource>& AFood, const int mut)
    : sex(sex),mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
    neutral_trait_alleles_f(neutral_traits_f), neutral_trait_alleles_m(neutral_traits_m),
    ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
    {
      age = 0;
      size = size_birth;
      repro_buffer = 0;
      Starve = false;
      Fecund = false;
      Is_dead = false;
      Mature = false;
      this->Mate_mut();
      this->Eco_mut();
      this->Neutral_mut();
      if (N_neutral) {
        neutral_trait = std::accumulate(neutral_trait_alleles_f.begin(), neutral_trait_alleles_f.end(), 0.0);
        neutral_trait += std::accumulate(neutral_trait_alleles_m.begin(), neutral_trait_alleles_m.end(), 0.0);}
      mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0)/(2*N_mating);
      mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0)/(2*N_mating);
      ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0);
      ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0);

      //Determine species identity//
      for (i = 0; i < Nr_Res - 1; ++i){
        SpeciesID = 1;
        if (ecological_trait >= SpeciesDiv[i]){
          SpeciesID = Nr_Res - i;
          break;
        }
      }

      for(i = 0, it_r = JFood.begin(); it_r < JFood.end(); ++it_r, ++i) {
          double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
          attack_constants.push_back(Const);
      }
      for(it_r = AFood.begin(); it_r < AFood.end(); ++it_r, ++i) {
        double Const = Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau));
        attack_constants.push_back(Const);
       }
    }


  /*------------------Deconstructor---------------------------------------*/
  ~Individual() {};

  /*------------------Member functions declarations---------------------------------------*/
  void R_Intake(std::vector<Resource>&, std::vector<Resource>&); //Food intake
  void Grow();
  void Die();

  /*------------------Mutations---------------------------------------*/
  void Mate_mut();

  void Neutral_mut();

  void Eco_mut();

  /*------------------Mating probability---------------------------------------*/
  Individual& MateProb(const Individual&, std::vector<double> &, double &, double &);
  void MateProb2(const Individual&, std::vector<double> &, double &, double &);


  /*------------------Mating function---------------------------------------*/
  void Mating(Individual&);
  void ClonalMating();

  /*-------------------Sorting based on age---------------------------------*/
  bool operator <(Individual const & IndividualObj)const;


  /*------------------Data members public---------------------------------------*/
  double mating_trait;
  double repro_buffer;
  std::vector<double> Intake;
  bool Is_dead;
  double matingProb;
  bool Fecund;
  double size;
    int SpeciesID;
    bool Starve;
    double NetProd;
    double age;
    double ecological_trait;
    double neutral_trait;
    double Matings = 0;
    bool Mature;
  private:
  /*------------------Data members private---------------------------------------*/

  int sex;
  std::vector<int> mating_trait_alleles_f;
  std::vector<int> mating_trait_alleles_m;
  std::vector<double> neutral_trait_alleles_f;
  std::vector<double> neutral_trait_alleles_m;
  std::vector<double> ecological_trait_alleles_f;
  std::vector<double> ecological_trait_alleles_m;
  std::vector<double> attack_constants;

  };


  extern std::vector<Individual> Females;
  extern std::vector<Individual> Males;
  extern std::vector<Individual> JuvFemales;
  extern std::vector<Individual> JuvMales;
  extern std::vector<Resource> JFood;
  extern std::vector<Resource> AFood;



#endif
