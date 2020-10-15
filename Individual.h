#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <numeric>
#include <random>
#include <functional>
#include <math.h>
#include "Pars.h"
#include "Resource.h"
#include <memory>
#include <algorithm>
#include "cmath"

extern int i;
extern int j;
extern std::vector<int>::iterator int_it; //iterator over integers
extern std::vector<double>::iterator int_doub; //iterator over doubles
extern std::mt19937 mt_rand; //random number gen
extern std::uniform_real_distribution<double> unif; //uniform dist
extern std::uniform_int_distribution<> IDNR_Rand;
extern std::normal_distribution<double> MutNorm; //normal dist
extern std::exponential_distribution<double> Surv_age; //exponential dist
extern std::vector<Resource>::iterator it_r; //iterator over resources
extern std::ostream& print_individualnames(std::ostream&); //print function
extern std::ostream& print_traitsindividualnames(std::ostream&); //print function
extern std::vector<double> SpeciesDiv;
extern int Nr_Res; //number of resources
extern double multi; //temp variable to hold multiplication factor
extern double phi; //Ontogenetic diet shift
extern double psi; //Allocation to maturation
extern double Total; //temp var to calculate food intake
extern double IntakeTot; //Variable that holds total food intake of indiviudal
extern double Value; //temp variable to calculate food intake
extern int Offspring; //Number of offspring
extern double dif; //Difference in trait between two parents
extern double eps;
extern int TotStarv;


class Individual {

/*------------------Friends---------------------------------------*/
  friend std::ostream& print_individual(std::ostream& os, const Individual& s);
  friend std::ostream& print_individualnames(std::ostream& os);
  friend std::ostream& print_traitsindividual(std::ostream& os, const Individual& s);

public:

  /*------------------Constructor---------------------------------------*/

    //Old individuals no mut, all 3 genotypes (assortative mating based on neutral trait)//
    Individual(double idnr, double age, double size,  int Mating, double reprobuf, double mxage, bool Mat,
    std::vector<double> mate_traits_f, std::vector<double> mate_traits_m,
    std::vector<double> neutral_traits_f, std::vector<double> neutral_traits_m,
    std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
    std::vector<Resource>& AllFood)
      : IDNR(idnr), size(size), age(age),  Mature(Mat),  Matings(Mating), repro_buffer(reprobuf),
      mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
      neutral_trait_alleles_f(neutral_traits_f), neutral_trait_alleles_m(neutral_traits_m),
      ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
       {
        repro_buffer = 0;
        //Starve = false;
        Starve = 0;
        Fecund = false;
        Is_dead = false;
        Mature = false;
        ecological_trait = 0;
        if (N_neutral) {
          neutral_trait = std::accumulate(neutral_trait_alleles_f.begin(), neutral_trait_alleles_f.end(), 0.0f);
          neutral_trait += std::accumulate(neutral_trait_alleles_m.begin(), neutral_trait_alleles_m.end(), 0.0f);}

        mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0f)/(2*N_mating);
        mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0f)/(2*N_mating);
        ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0f);
        ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0f);

        //strength of assortative mating
        //AssM = -0.5 * pow((pow(mating_trait, 2)/s_ass), 2);
        AssM = s_ass/pow(mating_trait, 2);

        //determine maximum age of an individual
        MaxAge = mxage;

        //Determine species identity//
        for (i = 0; i < Nr_Res - 1; ++i){
          SpeciesID = 1;
          if (ecological_trait >= SpeciesDiv[i]){
            SpeciesID = Nr_Res - i;
            break;
          }
        }

        //Feeding efficiencies
        for(it_r = AllFood.begin(); it_r < AllFood.end(); ++it_r) {
          if(it_r->Tau>=2000){
            attack_constants.push_back(Amax);
          } else {
          attack_constants.push_back(Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau)));
         }
      }
    }

    //Old individuals no mut, 2 genotypes (assortative mating based on eco trait)//
    Individual(double idnr, double age, double size,  int Mating, double reprobuf, double mxage, bool Mat,
      //std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
      std::vector<double> mate_traits_f, std::vector<double> mate_traits_m,
      std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
      std::vector<Resource>& AllFood)
      : IDNR(idnr), size(size), age(age),  Mature(Mat), Matings(Mating),repro_buffer(reprobuf),
      mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
      ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
       {
        //Starve = false;
          Starve = 0;
        Fecund = false;
        Is_dead = false;
        Mature = false;
        ecological_trait = 0;
        mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0f)/(2*N_mating);
        mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0f)/(2*N_mating);
        ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0f);
        ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0f);

        //strength of assortative mating
        if(mating_trait > 0){
          AssM = s_ass/pow(mating_trait,2);
        } else {
          AssM = s_diss/pow(mating_trait,2);
        }
        //AssM = -0.5 * pow((pow(mating_trait, 2)/s_ass), 2);
        //AssM = s_ass/pow(mating_trait,2);

        //determine maximum age of an individual
        MaxAge = mxage;

        //Determine species identity//
        for (i = 0; i < Nr_Res - 1; ++i){
          SpeciesID = 1;
          if (ecological_trait >= SpeciesDiv[i]){
            SpeciesID = Nr_Res - i;
            break;
          }
        }

        //Feeding efficiencies
         for(it_r = AllFood.begin(); it_r < AllFood.end(); ++it_r) {
           if(it_r->Tau>=2000){
             attack_constants.push_back(Amax);
           } else {
           attack_constants.push_back(Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau)));
          }}
      }

    //Old individuals no mut, 1 genotype (clonal repro)//
    Individual(double idnr, double age, double size,  int Mating, double reprobuf, double mxage, bool Mat,
      std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
      std::vector<Resource>& AllFood)
      : IDNR(idnr), size(size), age(age),  Mature(Mat), Matings(Mating),repro_buffer(reprobuf),
      ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
       {
        repro_buffer = 0;
        //Starve = false;
          Starve = 0;
        Fecund = false;
        Is_dead = false;
        Mature = false;
        ecological_trait = 0;
        ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0f);
        ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0f);

        //strenght of assortative mating
        //AssM = -0.5 * pow((pow(mating_trait, 2)/s_ass), 2);
        AssM = s_ass/pow(mating_trait,2);
        //determine maximum age of an individual
        MaxAge = mxage;

        //Determine species identity//
        for (i = 0; i < Nr_Res - 1; ++i){
          SpeciesID = 1;
          if (ecological_trait >= SpeciesDiv[i]){
            SpeciesID = Nr_Res - i;
            break;
          }
        }
        //These are correct//

        for(it_r = AllFood.begin(); it_r < AllFood.end(); ++it_r) {
          if(it_r->Tau>=2000){
            attack_constants.push_back(Amax);
          } else {
          attack_constants.push_back(Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau)));
         }
        }
      }


    //Newborn no mut//
    Individual(
      std::vector<double> mate_traits_f, std::vector<double> mate_traits_m,
      //std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
      std::vector<double> neutral_traits_f, std::vector<double> neutral_traits_m,
      std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
       std::vector<Resource>& AllFood)
      :  mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
      neutral_trait_alleles_f(neutral_traits_f), neutral_trait_alleles_m(neutral_traits_m),
      ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
       {
        age = 0;
        Matings = 0;
        size = size_birth;
        repro_buffer = 0;
        //Starve = false;
          Starve = 0;
        Fecund = false;
        Is_dead = false;
        Mature = false;
        IDNR = IDNR_Rand(mt_rand);
        if (N_neutral) {
          neutral_trait = std::accumulate(neutral_trait_alleles_f.begin(), neutral_trait_alleles_f.end(), 0.0f);
          neutral_trait += std::accumulate(neutral_trait_alleles_m.begin(), neutral_trait_alleles_m.end(), 0.0f);}
        mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0f)/(2*N_mating);
        mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0f)/(2*N_mating);
        ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0f);
        ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0f);

        //strength of assortative mating
        //AssM = -0.5 * pow((pow(mating_trait, 2)/s_ass), 2);
        AssM = s_ass/pow(mating_trait,2);

        //determine maximum age of an individual
        MaxAge = Surv_age(mt_rand);

        //Determine species identity//
        for (i = 0; i < Nr_Res - 1; ++i){
          SpeciesID = 1;
          if (ecological_trait >= SpeciesDiv[i]){
            SpeciesID = Nr_Res - i;
            break;
          }
        }

        //feeding efficiency
        for(it_r = AllFood.begin(); it_r < AllFood.end(); ++it_r) {
          if(it_r->Tau>=2000){
            attack_constants.push_back(Amax);
          } else {
          attack_constants.push_back(Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau)));
         }
        }

      }

  //Newborn with mut (it will mutate its traits first)//
  Individual(
    std::vector<double> mate_traits_f, std::vector<double> mate_traits_m,
    //std::vector<int> mate_traits_f, std::vector<int> mate_traits_m,
    std::vector<double> neutral_traits_f, std::vector<double> neutral_traits_m,
    std::vector<double> ecological_traits_f, std::vector<double> ecological_traits_m,
     std::vector<Resource>& AllFood, const int mut)
    : mating_trait_alleles_f(mate_traits_f), mating_trait_alleles_m(mate_traits_m),
    neutral_trait_alleles_f(neutral_traits_f), neutral_trait_alleles_m(neutral_traits_m),
    ecological_trait_alleles_f(ecological_traits_f), ecological_trait_alleles_m(ecological_traits_m)
    {
      age = 0;
      Matings = 0;
      size = size_birth;
      repro_buffer = 0;
      //Starve = false;
        Starve = 0; 
      Fecund = false;
      Is_dead = false;
      Mature = false;
      IDNR = IDNR_Rand(mt_rand);
      if(mut == 1){
      this->Mate_mut();
      this->Eco_mut();
      this->Neutral_mut();
      }
      if (N_neutral) {
        neutral_trait = std::accumulate(neutral_trait_alleles_f.begin(), neutral_trait_alleles_f.end(), 0.0f);
        neutral_trait += std::accumulate(neutral_trait_alleles_m.begin(), neutral_trait_alleles_m.end(), 0.0f);}
      mating_trait = std::accumulate(mating_trait_alleles_f.begin(), mating_trait_alleles_f.end(), 0.0f)/(2*N_mating);
      mating_trait += std::accumulate(mating_trait_alleles_m.begin(), mating_trait_alleles_m.end(), 0.0f)/(2*N_mating);
      ecological_trait = std::accumulate(ecological_trait_alleles_f.begin(), ecological_trait_alleles_f.end(), 0.0f);
      ecological_trait += std::accumulate(ecological_trait_alleles_m.begin(), ecological_trait_alleles_m.end(), 0.0f);

      //strenght of assortative mating
      //AssM = -0.5 * pow((pow(mating_trait, 2)/s_ass), 2);
      AssM = s_ass/pow(mating_trait,2);

      //determine maximum age of an individual
      MaxAge = Surv_age(mt_rand);

      //Determine species identity//
      for (i = 0; i < Nr_Res - 1; ++i){
        SpeciesID = 1;
        if (ecological_trait >= SpeciesDiv[i]){
          SpeciesID = Nr_Res - i;
          break;
        }
      }

      //feeding efficiency
      for(it_r = AllFood.begin(); it_r < AllFood.end(); ++it_r) {
        if(it_r->Tau>=2000){
          attack_constants.push_back(Amax);
        } else {
        attack_constants.push_back(Amax * exp(-pow(ecological_trait-it_r->OptTrait,2)/(2*it_r->Tau*it_r->Tau)));
       }
      }

    }


  /*------------------Deconstructor---------------------------------------*/
  ~Individual() {};

  /*------------------Member functions declarations---------------------------------------*/
  void R_Intake(std::vector<Resource>&, std::vector<double>&); //Food intake


  /*------------------Mutations---------------------------------------*/
//  void Mate_mut_diallic();

  void Mate_mut();

  void Neutral_mut();

  void Eco_mut();

  /*------------------Mating probability---------------------------------------*/
  Individual& MateProb(const Individual&, std::vector<double> &, double &);

  /*------------------Mating function---------------------------------------*/
  void SexualRepro(Individual&);
  void ClonalRepro();

  /*-------------------Sorting based on age ---------------------------------*/
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
  double IDNR;
  double Starve;
  double NetProd;
  double age;
  double ecological_trait;
  double neutral_trait;
  double Matings;
  bool Mature;
  double MaxAge;
  double AssM;

  private:
  /*------------------Data members private---------------------------------------*/

  std::vector<double> mating_trait_alleles_f;
  std::vector<double> mating_trait_alleles_m;
  std::vector<double> neutral_trait_alleles_f;
  std::vector<double> neutral_trait_alleles_m;
  std::vector<double> ecological_trait_alleles_f;
  std::vector<double> ecological_trait_alleles_m;
  std::vector<double> attack_constants;

  };


  extern std::vector<Resource> AllFood;
  extern std::vector<std::unique_ptr<Individual>> Advec;
  extern std::vector<std::unique_ptr<Individual>> Juvvec;


#endif
