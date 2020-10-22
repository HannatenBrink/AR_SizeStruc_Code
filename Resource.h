#ifndef RESOURCE_H
#define RESOURCE_H
#include <iostream>
#include "Pars.h"
#include <memory>

extern double OverCons;

class Resource {
  friend std::ostream& print_resource(std::ostream&, const Resource&);
  friend std::ostream& print_resourceDensity(std::ostream&, const Resource&);
  friend std::ostream& print_resourceDensityTotal(std::ostream&, const Resource&);
  friend std::ostream& print_resourceInfo(std::ostream&, const Resource&);

  public:
  //Constructors//

  Resource(double dens, const double MaxDens, const double Vol, const double Trait, const double tau, const std::string Rname)
     : Density(dens), Rmax(MaxDens), Volume(Vol), OptTrait(Trait), Tau(tau), Name(Rname) {}

  Resource(const double MaxDens, const double Vol, const double Trait, const double tau, const std::string Rname)
      : Rmax(MaxDens), Volume(Vol), OptTrait(Trait), Tau(tau), Name(Rname) {Density = MaxDens;}

  //Deconstructor//
  ~Resource() {};

  //data members//
  double Density;
  double Delay = 0;
  const double Rmax;
  const double Volume;
  const double OptTrait;
  const double Tau;
  const std::string Name;


  //member function declaration// //correct
  inline Resource& Growth(double &popIntake, double &Rmaxch) {
    /*Change in the resource due to growth and due to consumption
    It is possible to change the rmax value over time
    useful in case you want to do a bifurcation over the rmax value of a resource
    note that this is implemented as an INCREASE in the rmax value*/
    double Delta_R = rho * ((this->Rmax + Rmaxch) - this->Density)- popIntake/this->Volume;
    this->Density += Delta_R * delta_t;
    this->Density = std::max(0.0, this->Density);
    return *this;
  }

  //member function declaration// new one with delay in resource
  /*inline Resource& Growth(double &popIntake, double &Rmaxch) {
    double Delta_R = rho * ((this->Rmax + Rmaxch) - this->Density) - popIntake/this->Volume;
    this->Density += Delta_R * delta_t + this->Delay;
    if(this->Density < 0){
      this->Delay = this->Density;
      OverCons  += 1;
    }
    this->Density = std::max(0.0, this->Density);
    return *this;
  }*/


};


#endif
