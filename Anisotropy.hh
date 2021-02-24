/*
* $Id: Anisotropy.hh 20.02.2021 A Fijalkowska $
* 
* \file Anisotropy.hh
* \brief Class calculates angular distribution of gamma rays emitted from 
*  oriented state. 
* Input data:
* I0 - spin of mother (polarised) state
* am - population of nuclear m-states, avaliable m-states are -I..I, size of am vector is 2*I0+1
* Ii - spin of primary state
* If - spin of final state
* L1 - multipolarity of the transition
* L2 - second multipolarity of the transition
* Alghoritm calculate transitions that mix NO more than two multipolarities
* mixRatio - multipole mixing ratio
*/
#ifndef Anisotropy_h
#define Anisotropy_h 1

#include "TF1.h"
#include <vector>

class Anisotropy
{
  public:
     Anisotropy(double I0, std::vector<double> mVal, std::vector<double> am, double Ii, double If, double L1, double L2, double mixRatio);
     TF1* GetThetaDistr() {return thetaDistr;}
    ~Anisotropy();


  private: 
     double Wigner3j(std::vector<double> J123, std::vector<double> M123);
     double Wigner6j(std::vector<double> J123, std::vector<double> J456);
     double BParam(double I0, std::vector<double> m, std::vector<double> am, double lambda);
     double FParam(double Ii, double If, double L1, double L2, double lambda);
     double UCoef(double I0, double Ii, double L, double lambda);//if
     double AParam(double Ii, double If, double L1, double L2, double mixratio, double lambda);
     TF1* MakeThetaDistr(std::vector<double> angularCoeff );

     double TriangleCoefficient(double a, double b, double c);
     double Factorial(double x);
     bool IsVectorHalfInteger(std::vector<double> data);
     bool IsVectorInteger(std::vector<double> data);
     TF1* thetaDistr;
};

#endif // Anisotropy_h
