#include "Individual.h"


/*---------------------------------------Resource intake and growth---------------------------------------*/
void Individual::R_Intake(std::vector<Resource>& AllResource, std::vector<double>& FeedVec, std::mt19937& mt_rand_ref) {
  if(discrete_shift==1){
    if (this->size <= m_shift)
    phi = 1;
    else {
      phi = 0;
    }
  }
  else {
  phi = 1 - 1/(1 + exp(-(this->size - m_shift)));
    }
  Intake.clear();
  Fecund = false;

  /*it_r = AllResource.begin();
  Value = phi * this->attack_constants[0] * pow(this->size, qpow) * it_r->Density;
  Total = Value;
  Intake.push_back(Value);

  for(i=1, it_r = AllResource.begin()+1; it_r != AllResource.end(); ++it_r, ++i){
    Value = (1 - phi) * this->attack_constants[i] * pow(this->size, qpow) * it_r->Density;
    Total += Value;
    Intake.push_back(Value);
  }*/

  for(i = 0, it_r = AllResource.begin(); it_r != AllResource.end(); ++it_r, ++i){
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
    //std::cout << "With an age of " << age <<  " and a maxage of " << MaxAge << " the individual is too old " << std::endl;
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
    this->Starve += 1;
    TotStarv += 1;
    if (- this->NetProd / (Xi * this->size) * delta_t > unif(mt_rand_ref)){
        this->Is_dead = 1;
        //std::cout << "with a net loss of " << this->NetProd << " probability to die is " << - this->NetProd / (Xi * this->size) * delta_t << std::endl;
  }
}


}


/*--------------------Calculate mating probability---------------------------------------*/
//Following Doebeli & Dieckmann//
//where parameter s_ass/s_diss is the importance of assortative mating//
//AM depends on either a neutral trait or the ecological trait//
//I think it is fine like this, but it might be prone to weird matings in case of biased sex ratio/
//In that case, make the minimum mating prob bigger, and track how often a female can't find a mate//
//if it gets a lot, let it mate with a male that is less good//
//not sure if this will happen because of the sorting//
Individual& Individual::MateProb(const Individual& female, std::vector<double> &vec, double &total){
      if(this->Fecund)  {
      if (female.mating_trait != 0){
        if (N_neutral == 0) {
        dif = 1/(female.AssM*sqrt(2*M_PI)) * exp(-0.5*pow(((this->ecological_trait - female.ecological_trait)/female.AssM),2));} else {
        dif = 1/(female.AssM*sqrt(2*M_PI)) * exp(-0.5*pow(((this->neutral_trait - female.neutral_trait)/female.AssM),2));
      }
      if (female.mating_trait > 0){
        this->matingProb = dif;
      } else if (female.mating_trait < 0){
        this->matingProb = 1 - dif;
      }
    } else {
      this->matingProb = 1;
    }
      if(this->matingProb < pow(10,-10)){
      this->matingProb = 0;
    }
  } else {
    this->matingProb = 0;
  }
      //std::cout << "mating probability for a " << this->Fecund<< " male with eco trait " << this->ecological_trait << " and female with eco trait " << female.ecological_trait << " and mating trait " << female.mating_trait << " equals " << this->matingProb << '\n' << std::endl;
      total += this->matingProb;
      vec.push_back(total);
      return *this;
    }


/*------------------Mating function---------------------------------------*/


void Individual::SexualRepro(Individual& male, std::mt19937& mt_rand_ref) {
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
    //std::vector<int> mating_alleles_mother;
    //std::vector<int> mating_alleles_father;
    std::vector<double> mating_alleles_mother;
    std::vector<double> mating_alleles_father;
    std::vector<double> neutral_alleles_mother;
    std::vector<double> neutral_alleles_father;
    std::vector<double> eco_alleles_mother;
    std::vector<double> eco_alleles_father;

    //get the mating alleles of the mother//
    for(j = 0; j < N_mating; ++j) {
      if(unif(mt_rand_ref) > 0.5) {
        mating_alleles_mother.push_back(this->mating_trait_alleles_f[j]);} else {
        mating_alleles_mother.push_back(this->mating_trait_alleles_m[j]);
       }
    }
    //get the mating alleles of the father//
    for(j = 0; j < N_mating; ++j) {
      if(unif(mt_rand_ref) > 0.5) {
        mating_alleles_father.push_back(male.mating_trait_alleles_f[j]);}
        else {
        mating_alleles_father.push_back(male.mating_trait_alleles_m[j]);
      }
    }
    //get the neutral alleles of the mother//
    for(j = 0; j < N_neutral; ++j) {
      if(unif(mt_rand_ref) > 0.5) {neutral_alleles_mother.push_back(this->neutral_trait_alleles_f[j]);} else {
        neutral_alleles_mother.push_back(this->neutral_trait_alleles_m[j]);
      }
    }
    //get the neutral alleles of the father//
    for(j = 0; j < N_neutral; ++j) {
      if(unif(mt_rand_ref) > 0.5) {neutral_alleles_father.push_back(male.neutral_trait_alleles_f[j]);} else {
        neutral_alleles_father.push_back(male.neutral_trait_alleles_m[j]);
      }
    }
    //get the eco alleles of the mother//
    for(j = 0; j < N_eco; ++j) {
      if(unif(mt_rand_ref) > 0.5) {eco_alleles_mother.push_back(this->ecological_trait_alleles_f[j]);} else {
        eco_alleles_mother.push_back(this->ecological_trait_alleles_m[j]);
      }
    }
    //get the eco alleles of the father//
    for(j = 0; j < N_eco; ++j) {
      if(unif(mt_rand_ref) > 0.5) {eco_alleles_father.push_back(male.ecological_trait_alleles_f[j]);} else {
        eco_alleles_father.push_back(male.ecological_trait_alleles_m[j]);
      }
    }


    //Create a newborn, mutate the alleles, and add to the juvenile population//
    std::unique_ptr<Individual> IndivPtr(new Individual(mt_rand_ref, mating_alleles_mother, mating_alleles_father,
      neutral_alleles_mother, neutral_alleles_father,
      eco_alleles_mother, eco_alleles_father,
      AllFood, 1));
    //print_individual(std::cout, *IndivPtr);
    //std::cout << std::endl;
      Juvvec.push_back(move(IndivPtr));

  }
}

void Individual::ClonalRepro(std::mt19937& mt_rand_ref) {
  //decide upon number of offspring and clear the repro buffer//
  Offspring = std::trunc(this->repro_buffer/1);
  this->repro_buffer = std::fmod(this->repro_buffer, 1);
  this->Fecund = false;
  //for loop over number of offspring//
  for(int k = 0; k < Offspring; ++k){

    //create newborn, mutate alleles, and add to the population
      std::unique_ptr<Individual> IndivPtr(new Individual(mt_rand_ref, this->mating_trait_alleles_f, this->mating_trait_alleles_m,
        this->neutral_trait_alleles_f, this->neutral_trait_alleles_m,
        this->ecological_trait_alleles_f , this->ecological_trait_alleles_m,
    AllFood, 1));
      Juvvec.push_back(move(IndivPtr));

  }
}

/*---------------------------------------Mutate---------------------------------------*/
/*inline void Individual::Mate_mut_diallic(){
  for(int_it = this->mating_trait_alleles_f.begin(); int_it != this->mating_trait_alleles_f.end(); ++int_it) {
    if(unif(mt_rand_ref) < mut_rate_di) {
      *int_it *= -1;
    }
  }
  for(int_it = this->mating_trait_alleles_m.begin(); int_it != this->mating_trait_alleles_m.end(); ++int_it) {
    if(unif(mt_rand_ref) < mut_rate_di) {
      *int_it *= -1;
    }
  }
}*/

inline void Individual::Mate_mut(std::mt19937& mt_rand_ref){
  for(int_doub = this->mating_trait_alleles_f.begin(); int_doub != this->mating_trait_alleles_f.end(); ++int_doub) {
    if(unif(mt_rand_ref) < mut_rate_mate) {
      *int_doub += MutNorm(mt_rand_ref);

      if(*int_doub > 1){
        *int_doub = 1;
      } else if(*int_doub < -1){
        *int_doub = -1;
      }
    }
  }
  for(int_doub = this->mating_trait_alleles_m.begin(); int_doub != this->mating_trait_alleles_m.end(); ++int_doub) {
    if(unif(mt_rand_ref) < mut_rate_mate) {
    *int_doub += MutNorm(mt_rand_ref);
    if(*int_doub > 1){
      *int_doub = 1;
    } else if(*int_doub < -1){
      *int_doub = -1;
    }
    }
  }
}

inline void Individual::Neutral_mut(std::mt19937& mt_rand_ref){

  for(int_doub = this->neutral_trait_alleles_f.begin(); int_doub != this->neutral_trait_alleles_f.end(); ++int_doub) {
    if(unif(mt_rand_ref) < mut_rate_neut) {
      *int_doub += MutNorm(mt_rand_ref);
    }
  }
  for(int_doub = this->neutral_trait_alleles_m.begin(); int_doub != this->neutral_trait_alleles_m.end(); ++int_doub) {
    if(unif(mt_rand_ref) < mut_rate_neut) {
      *int_doub += MutNorm(mt_rand_ref);
    }
  }
}

inline void Individual::Eco_mut(std::mt19937& mt_rand_ref) {

  for(int_doub = ecological_trait_alleles_f.begin(); int_doub != ecological_trait_alleles_f.end(); ++int_doub) {
    if(unif(mt_rand_ref) < mut_rate_eco) {
      *int_doub += MutNorm(mt_rand_ref);
      //*int_doub = abs(*int_doub); //This is not necessary I think?
    }}
    for(int_doub = this->ecological_trait_alleles_m.begin(); int_doub != this->ecological_trait_alleles_m.end(); ++int_doub) {
      if(unif(mt_rand_ref) < mut_rate_eco) {
        *int_doub += MutNorm(mt_rand_ref);
        //*int_doub = abs(*int_doub);
      }}
    }

/*-----Sorting operators---*/
inline bool Individual::operator <(Individual const& IndividualObj)const
	{
		return age < IndividualObj.age;
	}
