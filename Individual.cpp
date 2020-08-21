#include "Individual.h"


/*---------------------------------------Resource intake---------------------------------------*/
void Individual::R_Intake(std::vector<Resource>& JResource, std::vector<Resource>& AResource) {
  double phi = 1 - 1/(1 + exp(-(this->size - m_shift))); //Ontogenetic shift
  double Total = 0;
  double IntakeTot = 0;
  double Value;
  Intake.clear();
  this->Starve = false;
  //Size-specific attack rates are correct//
  for(i = 0, it_r = JResource.begin(); it_r < JResource.end(); ++it_r, ++i){
    //Intake[i] = 0;
    Value = phi * this->attack_constants[i] * pow(this->size, qpow) * it_r->Density;
    Total += Value;
    Intake.push_back(Value);
  }
  for(it_r = AResource.begin(); it_r < AResource.end(); ++it_r, ++i){
    //Intake[i] = 0;
    Value = (1 - phi) * this->attack_constants[i] * pow(this->size, qpow) * it_r->Density;
    Total += Value;
    Intake.push_back(Value);
  }

  IntakeTot = h * pow(this->size, npow) * Total / (Total + h * pow(this->size, npow)); //Total Intake is correct
  for(i = 0, it_r = JResource.begin(); it_r < JResource.end(); ++it_r, ++i){
    Intake[i] *= h * pow(this->size, npow) / (Total + h * pow(this->size, npow)); //correct
  }
  for(it_r = AResource.begin(); it_r < AResource.end(); ++it_r, ++i){
    Intake[i] *= h * pow(this->size, npow) / (Total + h * pow(this->size, npow)); //correct
  }

  this->NetProd = alphapar * IntakeTot - kmet * pow(this->size, pmain); //correct
  if(NetProd < 0) {
    this->Starve = true;
  }
}

/*---------------------------------------Growth---------------------------------------*/
void Individual::Grow(){
   age += delta_t;
   if(!Starve) {
     double psi = pow((1 + pow(this->size/M_Mat,-u)),-1) * pow((eta * this->size/M_Mat), (1-npow));
     this->size += (1 - psi) * this->NetProd * delta_t;
       this->repro_buffer += psi * this->NetProd * epsilonpar * pow(size_birth, -1) * 0.5 * delta_t;
     if(this->repro_buffer >= 1){
       Fecund = true;
       Mature = true;
     } else {Fecund = false;}
   }
}

/*---------------------------------------Death---------------------------------------*/
void Individual:: Die() {
  double prob = unif(mt_rand);
  double mu_starv = 0;
  if (this->Starve) {
    mu_starv = - this->NetProd / (Xi * this->size);
  }
  double mu_total = (mu_b + mu_starv);
  if (mu_total * delta_t > prob){
      this->Is_dead = 1;
  }
}

/*--------------------Calculate mating probability---------------------------------------*/
//Following Bolnick & Doebeli//
//where parameter s is the importance of assortative mating//
//input fac is the scaling of the normal distribution//
//AM depends on either a neutral trait or the ecological trait//
Individual& Individual::MateProb(const Individual& female, std::vector<double> &vec, double &total, double &fac){
      double dif;
      if (this->Fecund){
      if (female.mating_trait == 0){
        this->matingProb = 1;
      } else if (female.mating_trait > 0) {
        if (N_neutral == 0) {
          dif = pow((female.ecological_trait - this->ecological_trait),2);
        } else {
          dif = pow((female.neutral_trait - this->neutral_trait),2);
        }
        this->matingProb = exp(fac * dif);
      } else {
        if (N_neutral == 0) {
          dif = pow(2 - std::abs(female.ecological_trait - this->ecological_trait),2);
        } else {
          dif = pow(2 - std::abs(female.neutral_trait - this->neutral_trait),2);
        }
        this->matingProb = exp(fac * dif);
      }}
      else {
        this->matingProb = 0;
      }
      total += this->matingProb;
      vec.push_back(total);
      return *this;
    }


/*------------------Mating function---------------------------------------*/

void Individual::Mating(Individual& male) {
  int sex;
  //decide upon number of offspring and clear the repro buffer//
  int Offspring = std::trunc(this->repro_buffer/1) + std::trunc(male.repro_buffer/1);
  //std::cout << "Total Offspring is " << Offspring << std::endl;
  this->repro_buffer = std::fmod(this->repro_buffer, 1);
  //std::cout << "Remaining female offspring is " << this->repro_buffer << std::endl;
  male.repro_buffer = std::fmod(male.repro_buffer, 1);
  this->Fecund = false;
  male.Fecund = false;
  //std::cout << "Remaining male offspring is " << male.repro_buffer << std::endl;
  //for loop over number of offspring//
  for(int k = 0; k < Offspring; ++k){
    //std::cout << "Offspring number " << i << std::endl;
    std::vector<int> mating_alleles_mother;
    std::vector<int> mating_alleles_father;
    std::vector<double> neutral_alleles_mother;
    std::vector<double> neutral_alleles_father;
    std::vector<double> eco_alleles_mother;
    std::vector<double> eco_alleles_father;

    //get the mating alleles of the mother//
    //std::cout << "Get mating alleles of the mom" << std::endl;
    for(j = 0; j < N_mating; ++j) {
      if(unif(mt_rand) > 0.5) {mating_alleles_mother.push_back(this->mating_trait_alleles_f[j]);} else {
        mating_alleles_mother.push_back(this->mating_trait_alleles_m[j]);
      }
      int_it = mating_alleles_mother.end() - 1;
      //std::cout << "The new mating allele of the mom is is: " << *int_it << std::endl;
    }
    //get the mating alleles of the father//
    //std::cout << "Get mating alleles of the dad" << std::endl;
    for(j = 0; j < N_mating; ++j) {
      if(unif(mt_rand) > 0.5) {mating_alleles_father.push_back(male.mating_trait_alleles_f[j]);} else {
        mating_alleles_father.push_back(male.mating_trait_alleles_m[j]);
      }
      int_it = mating_alleles_father.end() - 1;
      //std::cout << "The new mating allele of the dad is is: " << *int_it << std::endl;
    }
    //get the neutral alleles of the mother//
    //  std::cout << "Get neutral alleles of the mom" << std::endl;
    for(j = 0; j < N_neutral; ++j) {
      if(unif(mt_rand) > 0.5) {neutral_alleles_mother.push_back(this->neutral_trait_alleles_f[j]);} else {
        neutral_alleles_mother.push_back(this->neutral_trait_alleles_m[j]);
      }
      int_doub = neutral_alleles_mother.end() - 1;
    //  std::cout << "The new neutral allele of the mom is is: " << *int_it << std::endl;
    }
    //std::cout << "Get neutral alleles of the dad" << std::endl;
    for(j = 0; j < N_neutral; ++j) {
      if(unif(mt_rand) > 0.5) {neutral_alleles_father.push_back(male.neutral_trait_alleles_f[j]);} else {
        neutral_alleles_father.push_back(male.neutral_trait_alleles_m[j]);
      }
      int_doub = neutral_alleles_father.end() - 1;
    //  std::cout << "The new neutral allele of the father is is: " << *int_it << std::endl;
    }
    //get the eco alleles of the mother//
    //std::cout << "Get eco alleles of the mom" << std::endl;
    for(j = 0; j < N_eco; ++j) {
      if(unif(mt_rand) > 0.5) {eco_alleles_mother.push_back(this->ecological_trait_alleles_f[j]);} else {
        eco_alleles_mother.push_back(this->ecological_trait_alleles_m[j]);
      }
      int_doub = eco_alleles_mother.end() - 1;
    //  std::cout << "The new eco allele of the mom is is: " << *int_doub << std::endl;
    }
    //get the eco alleles of the father//
    //std::cout << "Get eco alleles of the dad" << std::endl;
    for(j = 0; j < N_eco; ++j) {
      if(unif(mt_rand) > 0.5) {eco_alleles_father.push_back(male.ecological_trait_alleles_f[j]);} else {
        eco_alleles_father.push_back(male.ecological_trait_alleles_m[j]);
      }
      int_doub = eco_alleles_father.end() - 1;
    //  std::cout << "The new eco allele of the dad is is: " << *int_doub << std::endl;
    }


    //determine the sex of the newborn//
    if(unif(mt_rand) > 0.5) {sex = 1;} else {sex = 0;}
    Individual NewBorn(mating_alleles_mother, mating_alleles_father,
      neutral_alleles_mother, neutral_alleles_father,
      eco_alleles_mother, eco_alleles_father,
      sex, JFood, AFood, 1);
    if(sex){
      JuvFemales.push_back(NewBorn);
      //std::cout << "New female" << std::endl;
      //print_traitsindividualnames(std::cout);
      //print_traitsindividual(std::cout, NewBorn);
    } else {
      JuvMales.push_back(NewBorn);
      //std::cout << "New male" << std::endl;
      //print_traitsindividualnames(std::cout);
      //print_traitsindividual(std::cout, NewBorn);
    }
  }
}

void Individual::ClonalMating() {
  int sex;
  //decide upon number of offspring and clear the repro buffer//
  int Offspring = std::trunc(this->repro_buffer/1);
  this->repro_buffer = std::fmod(this->repro_buffer, 1);
  this->Fecund = false;
  //for loop over number of offspring//
  for(int k = 0; k < Offspring; ++k){
    //determine the sex of the newborn//
    if(unif(mt_rand) > 0.5) {sex = 1;} else {sex = 0;}
    Individual NewBorn(this->mating_trait_alleles_f, this->mating_trait_alleles_m,
      this->neutral_trait_alleles_f, this->neutral_trait_alleles_m,
      this->ecological_trait_alleles_f , this->ecological_trait_alleles_m,
      sex, JFood, AFood, 1);
    if(sex){
      JuvFemales.push_back(NewBorn);
    } else {
      JuvMales.push_back(NewBorn);
    }
  }
}

/*---------------------------------------Mutate---------------------------------------*/
inline void Individual::Mate_mut(){
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
      //std::cout << "Old val " << *int_doub << std::endl;
      *int_doub += MutNorm(mt_rand);
      //std::cout << "New val " << *int_doub << std::endl;
    }}
    for(int_doub = this->ecological_trait_alleles_m.begin(); int_doub != this->ecological_trait_alleles_m.end(); ++int_doub) {
      if(unif(mt_rand) < mut_rate) {
        //std::cout << "Old val " << *int_doub << std::endl;
        *int_doub += MutNorm(mt_rand);
        //std::cout << "New val " << *int_doub << std::endl;
      }}
    }

/*-----Sorting operator----*/
inline bool Individual::operator <(Individual const& IndividualObj)const
	{
		return age < IndividualObj.age;
	}
