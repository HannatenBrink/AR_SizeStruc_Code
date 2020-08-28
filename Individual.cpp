#include "Individual.h"


/*---------------------------------------Resource intake and growth---------------------------------------*/
void Individual::R_Intake(std::vector<Resource>& AllResource, std::vector<double>& FeedVec) {

  phi = 1 - 1/(1 + exp(-(this->size - m_shift)));
  Intake.clear();
  Starve = false;
  Fecund = false;

  for(i = 0, it_r = AllResource.begin(); it_r < AllResource.end(); ++it_r, ++i){
    if(i > 0){
    Value = (1 - phi) * this->attack_constants[i] * pow(this->size, qpow) * it_r->Density;
    Total += Value;
    Intake.push_back(Value);
  } else {
    Value = phi * this->attack_constants[i] * pow(this->size, qpow) * it_r->Density;
    Total = Value;
    Intake.push_back(Value);
  }}

  multi = hpar* pow(this->size, npow) / (Total + hpar* pow(this->size, npow));
  IntakeTot = multi * Total;

  std::transform(Intake.begin(), Intake.end(), Intake.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, multi));

  std::transform (FeedVec.begin(), FeedVec.end(), Intake.begin(), FeedVec.begin(), std::plus<double>());

  this->NetProd = alphapar * IntakeTot - kmet * pow(this->size, pmain);

  //Growth and death//
  age += delta_t; //increase in age
  if(age >= MaxAge){ //Old and therefore dead
    this->Is_dead = 1;
  }
  else if(NetProd > 0) {  //growth
    psi = pow((1 + pow(this->size/M_Mat,-u)),-1) * pow((eta * this->size/M_Mat), (1-npow));
    this->size += (1 - psi) * this->NetProd * delta_t;
    this->repro_buffer += psi * this->NetProd * epsilonpar * pow(size_birth, -1) * 0.5 * delta_t;
    if(this->repro_buffer >= 1){
      Fecund = true;
      Mature = true;
    }
  } else {  //starvation mortality
    this->Starve = true;
    if (- this->NetProd / (Xi * this->size) * delta_t > unif(mt_rand)){
        this->Is_dead = 1;
  }
}

}


/*--------------------Calculate mating probability---------------------------------------*/
//Following Bolnick & Doebeli//
//where parameter s is the importance of assortative mating//
//AM depends on either a neutral trait or the ecological trait//
Individual& Individual::MateProb(const Individual& female, std::vector<double> &vec, double &total){
      if (this->Fecund){
      if (female.mating_trait == 0){
        this->matingProb = 1;
      } else if (female.mating_trait > 0) {
        if (N_neutral == 0) {
          dif = pow((female.ecological_trait - this->ecological_trait),2);
        } else {
          dif = pow((female.neutral_trait - this->neutral_trait),2);
        }
        this->matingProb = exp(female.AssM * dif);
      } else {
        if (N_neutral == 0) {
          dif = pow(2 - std::abs(female.ecological_trait - this->ecological_trait),2);
        } else {
          dif = pow(2 - std::abs(female.neutral_trait - this->neutral_trait),2);
        }
        this->matingProb = exp(female.AssM * dif);
      }}
      else {
        this->matingProb = 0;
      }
      total += this->matingProb;
      vec.push_back(total);
      return *this;
    }


/*------------------Mating function---------------------------------------*/


void Individual::SexualRepro(Individual& male) {
  //decide upon number of offspring and clear the repro buffer//
  //print_individualnames(std::cout);
  //std::cout << std::endl;

  Offspring = std::trunc(this->repro_buffer/1) + std::trunc(male.repro_buffer/1);
  this->repro_buffer = std::fmod(this->repro_buffer, 1);
  male.repro_buffer = std::fmod(male.repro_buffer, 1);
  this->Fecund = false;
  male.Fecund = false;
  /*std::cout << "Mother and father traits: " << std::endl;
  print_individual(std::cout, *this);
  std::cout << std::endl;
  print_individual(std::cout, male);
  std::cout << std::endl << std::endl;*/
  //for loop over number of offspring//
  for(int k = 0; k < Offspring; ++k){
    std::vector<int> mating_alleles_mother;
    std::vector<int> mating_alleles_father;
    std::vector<double> neutral_alleles_mother;
    std::vector<double> neutral_alleles_father;
    std::vector<double> eco_alleles_mother;
    std::vector<double> eco_alleles_father;

    //get the mating alleles of the mother//
    for(j = 0; j < N_mating; ++j) {
      if(unif(mt_rand) > 0.5) {
        mating_alleles_mother.push_back(this->mating_trait_alleles_f[j]);} else {
        mating_alleles_mother.push_back(this->mating_trait_alleles_m[j]);
       }
    }
    //get the mating alleles of the father//
    for(j = 0; j < N_mating; ++j) {
      if(unif(mt_rand) > 0.5) {
        mating_alleles_father.push_back(male.mating_trait_alleles_f[j]);}
        else {
        mating_alleles_father.push_back(male.mating_trait_alleles_m[j]);
      }
    }
    //get the neutral alleles of the mother//
    for(j = 0; j < N_neutral; ++j) {
      if(unif(mt_rand) > 0.5) {neutral_alleles_mother.push_back(this->neutral_trait_alleles_f[j]);} else {
        neutral_alleles_mother.push_back(this->neutral_trait_alleles_m[j]);
      }
    }
    //get the neutral alleles of the father//
    for(j = 0; j < N_neutral; ++j) {
      if(unif(mt_rand) > 0.5) {neutral_alleles_father.push_back(male.neutral_trait_alleles_f[j]);} else {
        neutral_alleles_father.push_back(male.neutral_trait_alleles_m[j]);
      }
    }
    //get the eco alleles of the mother//
    for(j = 0; j < N_eco; ++j) {
      if(unif(mt_rand) > 0.5) {eco_alleles_mother.push_back(this->ecological_trait_alleles_f[j]);} else {
        eco_alleles_mother.push_back(this->ecological_trait_alleles_m[j]);
      }
    }
    //get the eco alleles of the father//
    for(j = 0; j < N_eco; ++j) {
      if(unif(mt_rand) > 0.5) {eco_alleles_father.push_back(male.ecological_trait_alleles_f[j]);} else {
        eco_alleles_father.push_back(male.ecological_trait_alleles_m[j]);
      }
    }


    //determine the sex of the newborn//
    if(unif(mt_rand) > 0.5) {sex_off = 1;} else {sex_off = 0;}

    //Create a newborn, mutate the alleles, and add to the juvenile population//
    std::unique_ptr<Individual> IndivPtr(new Individual(mating_alleles_mother, mating_alleles_father,
      neutral_alleles_mother, neutral_alleles_father,
      eco_alleles_mother, eco_alleles_father,
      sex_off,  AllFood, 1));
    //print_individual(std::cout, *IndivPtr);
    //std::cout << std::endl;
    if(sex_off == 1){
      JuvFemalesvec.push_back(move(IndivPtr));
    } else {
      JuvMalesvec.push_back(move(IndivPtr));
    }
  }
}

void Individual::ClonalRepro() {
  //decide upon number of offspring and clear the repro buffer//
  Offspring = std::trunc(this->repro_buffer/1);
  this->repro_buffer = std::fmod(this->repro_buffer, 1);
  this->Fecund = false;
  //for loop over number of offspring//
  for(int k = 0; k < Offspring; ++k){
    //determine the sex of the newborn//
    if(unif(mt_rand) > 0.5) {sex_off = 1;} else {sex_off = 0;}
    //create newborn, mutate alleles, and add to the population
      std::unique_ptr<Individual> IndivPtr(new Individual(this->mating_trait_alleles_f, this->mating_trait_alleles_m,
        this->neutral_trait_alleles_f, this->neutral_trait_alleles_m,
        this->ecological_trait_alleles_f , this->ecological_trait_alleles_m,
        sex_off,  AllFood, 1));
    if(sex_off == 1){
      JuvFemalesvec.push_back(move(IndivPtr));
    } else {
      JuvMalesvec.push_back(move(IndivPtr));
    }
  }
}

/*---------------------------------------Mutate---------------------------------------*/
inline void Individual::Mate_mut_diallic(){
  for(int_it = this->mating_trait_alleles_f.begin(); int_it != this->mating_trait_alleles_f.end(); ++int_it) {
    if(unif(mt_rand) < mut_rate_di) {
      *int_it *= -1;
    }
  }
  for(int_it = this->mating_trait_alleles_m.begin(); int_it != this->mating_trait_alleles_m.end(); ++int_it) {
    if(unif(mt_rand) < mut_rate_di) {
      *int_it *= -1;
    }
  }
}

inline void Individual::Mate_mut(){
  for(int_it = this->mating_trait_alleles_f.begin(); int_it != this->mating_trait_alleles_f.end(); ++int_it) {
    if(unif(mt_rand) < mut_rate) {
      *int_doub += MutNorm(mt_rand);
      if(*int_doub > 1){
        *int_doub = 1;
      } else if(*int_doub < -1){
        *int_doub = -1;
      }
    }
  }
  for(int_it = this->mating_trait_alleles_m.begin(); int_it != this->mating_trait_alleles_m.end(); ++int_it) {
    if(unif(mt_rand) < mut_rate) {
    *int_doub += MutNorm(mt_rand);
    if(*int_doub > 1){
      *int_doub = 1;
    } else if(*int_doub < -1){
      *int_doub = -1;
    }
    }
  }
}

inline void Individual::Neutral_mut(){
  for(int_doub = this->neutral_trait_alleles_f.begin(); int_doub != this->neutral_trait_alleles_f.end(); ++int_doub) {
    if(unif(mt_rand) < mut_rate) {
      *int_doub += MutNorm(mt_rand);
    }
  }
  for(int_doub = this->neutral_trait_alleles_m.begin(); int_doub != this->neutral_trait_alleles_m.end(); ++int_doub) {
    if(unif(mt_rand) < mut_rate) {
      *int_doub += MutNorm(mt_rand);
    }
  }
}

inline void Individual::Eco_mut() {
  for(int_doub = ecological_trait_alleles_f.begin(); int_doub != ecological_trait_alleles_f.end(); ++int_doub) {
    if(unif(mt_rand) < mut_rate) {
      *int_doub += MutNorm(mt_rand);
    }}
    for(int_doub = this->ecological_trait_alleles_m.begin(); int_doub != this->ecological_trait_alleles_m.end(); ++int_doub) {
      if(unif(mt_rand) < mut_rate) {
        *int_doub += MutNorm(mt_rand);
      }}
    }

/*-----Sorting operator----*/
inline bool Individual::operator <(Individual const& IndividualObj)const
	{
		return age < IndividualObj.age;
	}
