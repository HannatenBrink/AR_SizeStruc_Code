#ifndef PARS_H
#define PARS_H
#include <vector>
#include <memory>

extern std::vector<double> Setting;
extern std::vector<double> Parameter;

/*------------All Parameter definitions-----------------------*/
#define rho          Parameter[ 0]  //Resource overturn
#define RJmax        Parameter[ 1]  //Max resource of juveniles
#define RAmax        Parameter[ 2]  //max resource of Adults

#define ThetaB       Parameter[ 3]  //Optimal traits
#define Theta1       Parameter[ 4]
#define Theta2       Parameter[ 5]
#define Theta3       Parameter[ 6]
#define Theta4       Parameter[ 7]
#define Theta5       Parameter[ 8]
#define Theta6       Parameter[ 9]

#define TauB         Parameter[10]  //CurveWidths
#define Tau1         Parameter[11]
#define Tau2         Parameter[12]
#define Tau3         Parameter[13]
#define Tau4         Parameter[14]
#define Tau5         Parameter[15]
#define Tau6         Parameter[16]

#define Amax         Parameter[17]
#define qpow         Parameter[18]
#define hpar         Parameter[19]
#define npow         Parameter[20]
#define alphapar     Parameter[21]
#define kmet         Parameter[22]
#define pmain        Parameter[23]
#define u            Parameter[24]
#define eta          Parameter[25]

#define epsilonpar   Parameter[26]
#define mu_b         Parameter[27]
#define Xi           Parameter[28]
#define size_birth   Parameter[29]
#define M_Mat        Parameter[30]
#define m_shift      Parameter[31]


#define s_ass        Parameter[32]
#define N_eco        Parameter[33]
#define N_neutral    Parameter[34] //If this one is zero, ecological character determines assortative mating. Otherwise the neutral trait
#define N_mating     Parameter[35]
#define mut_std      Parameter[36]
#define mut_rate     Parameter[37]
#define mut_rate_di  Parameter[38]

#define delta_t      Parameter[39]
#define volume       Parameter[40]

#define clonal       Parameter[41]

#define ini_eco      Parameter[42]
#define ini_mate     Parameter[43]
#define ini_neutral  Parameter[44]



/*------------All settings definitions-----------------------*/
#define Max_time        Setting[0]      //Maximum timesteps
#define Output_time     Setting[1]      //How often write to timefile
#define Pop_Output      Setting[2]      //How often write a full pop file
#define Pop_TraitOutput Setting[3]      //How often write traits only
#define MateFile        Setting[4]      //How often write to matefile

#endif