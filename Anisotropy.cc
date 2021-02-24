#include <fstream>
#include "TMath.h"
#include "Anisotropy.hh"

Anisotropy::Anisotropy(double I0, std::vector<double> mVal, std::vector<double> am, double Ii, double If, double L1, double L2, double mixRatio)
{
  R__LOAD_LIBRARY(libMathMore);  
  double lambdaMax=2.*L1;
  std::vector<double> angularCoeff(lambdaMax+1);
   for(int i=0; i<=lambdaMax; ++i)
   {
      double A = AParam(Ii, If, L1, L2, mixRatio, i);
      double B = BParam(I0, mVal, am, i);
      double U = UCoef(I0, Ii, 1, i);//L=1 GT transition
      angularCoeff.at(i)=A*B*U;   
   }
   thetaDistr = MakeThetaDistr(angularCoeff);
}

Anisotropy::~Anisotropy()
{
   delete thetaDistr;
}



double Anisotropy::BParam(double I0, std::vector<double> m, std::vector<double> am, double lambda)
{
   if(m.size() != am.size())
   {
      std::cout << "double BParam error: invalid input parameters" << std::endl;
      return 0;
   }
   //TODO B=0 for lambdra>2m
   
   double sum=0;
   for(int i=0; i!= m.size(); ++i)
   { 
      std::vector<double> WignerUp = {I0, I0, lambda};
      std::vector<double> WignerDown = {-m.at(i), m.at(i), 0};
      sum+=TMath::Power(-1, I0+m.at(i))*Wigner3j(WignerUp, WignerDown)*am.at(i);
   }
   return sum*TMath::Power((2*lambda+1)*(2*I0+1), 0.5);
}

double Anisotropy::FParam(double Ii, double If, double L1, double L2, double lambda)
{
   std::vector<double> Wigner3jUp = {L1, L2, lambda};
   std::vector<double> Wigner3jDown = {1, -1, 0};
   std::vector<double> Wigner6jUp = {L1, L2, lambda};
   std::vector<double> Wigner6jDown = {Ii, Ii, If};
   double fParam = TMath::Power(-1, If+Ii+1) 
                   * TMath::Power( (2*lambda+1)*(2*L1+1)*(2*L2+1)*(2*Ii+1), 0.5 )
                   * Wigner3j(Wigner3jUp,Wigner3jDown)
                   * Wigner6j(Wigner6jUp, Wigner6jDown);
   return fParam;
}

double Anisotropy::UCoef(double I0, double Ii, double L, double lambda)
{
   std::vector<double> Wigner6jUp = {I0, I0, lambda};
   std::vector<double> Wigner6jDown = {Ii, Ii, L};
   return TMath::Power(-1, I0+Ii+L+lambda)
          *TMath::Power((2*I0+1)*(2*Ii+1), 0.5)
          *Wigner6j(Wigner6jUp, Wigner6jDown);
}

double Anisotropy::AParam(double Ii, double If, double L1, double L2, double mixratio, double lambda)
{
   return(FParam( Ii, If, L1, L2, lambda)+2.*mixratio*FParam(Ii, If, L2, L2, lambda) +2.*mixratio*FParam(Ii, If, L1, L2, lambda) )/(1+mixratio*mixratio);
}

TF1* Anisotropy::MakeThetaDistr(std::vector<double> angularCoeff )
{   
   TF1* thetaDistr = new TF1("thetaDistr","1.+[0]+ [1]*cos(x) + [2]*0.5*(3*cos(x)**2-1)  + [3]*0.5*(5*cos(x)**3-3*cos(x)) + [4]*1./8.*(35.*cos(x)**4-30.*cos(x)**2+3.) + [5]*1./8.*(63.*cos(x)**5-70.*cos(x)**3+15.*cos(x) + [6]*1./16.*(231*cos(x)**6-315.*cos(x)**4+105*cos(x)**2-5.))",0, TMath::Pi());
   //thetaDistr->SetParameter(0,angularCoeff.at(0));
   thetaDistr->SetParameter(0,0);
   int i=1;
   for(i; i!=angularCoeff.size(); ++i)
   {
      thetaDistr->SetParameter(i,angularCoeff.at(i));
   }
   
   for(i; i!=7; ++i)
   {
      thetaDistr->SetParameter(i,0);
   }
   
   return thetaDistr;

}


double Anisotropy::Wigner6j(std::vector<double> J123, std::vector<double> J456)
{
   double j1 = J123.at(0); 
   double j2 = J123.at(1); 
   double j3 = J123.at(2);
   double j4 = J456.at(0); 
   double j5 = J456.at(1); 
   double j6 = J456.at(2);

   double tri1 = TriangleCoefficient(j1,j2,j3);
   double tri2 = TriangleCoefficient(j1,j5,j6);
   double tri3 = TriangleCoefficient(j4,j2,j6);
   double tri4 = TriangleCoefficient(j4,j5,j3);

   if (tri1==0||tri2==0||tri3==0||tri4==0)
     return 0;

   // Finding the range of summation in the Racah formula.
   std::vector<double> a(4);
   a.at(0) = j1 + j2 + j3;
   a.at(1) = j1 + j5 + j6;
   a.at(2) = j4 + j2 + j6;
   a.at(3) = j4 + j5 + j3;

   double rangei =*max_element(a.begin(), a.end());    
   std::vector<double> k(3);
   k.at(0) = j1 + j2 + j4 + j5;
   k.at(1) = j2 + j3 + j5 + j6;
   k.at(2) = j3 + j1 + j6 + j4;
   double rangef = *min_element(k.begin(), k.end());
   double Wigner = 0;
   
   for(double t=rangei; t<=rangef; ++t)
   {
      double denominator = Factorial(t-a.at(0)) * Factorial(t-a.at(1))
                          *Factorial(t-a.at(2)) * Factorial(t-a.at(3))
                          *Factorial(k.at(0)-t) * Factorial(k.at(1)-t) * Factorial(k.at(2)-t);
      Wigner += TMath::Power(-1, t)*Factorial(t+1)/denominator;     
   }   
   Wigner = TMath::Power(tri1*tri2*tri3*tri4,0.5)*Wigner;
   return Wigner;
}

double Anisotropy::Wigner3j(std::vector<double> J123, std::vector<double> M123)
{
   double j1 = J123.at(0); 
   double j2 = J123.at(1); 
   double j3 = J123.at(2);
   double m1 = M123.at(0); 
   double m2 = M123.at(1); 
   double m3 = M123.at(2);

// Input error checking
   if ( j1 < 0 || j2 < 0 || j3 < 0 )
   {
       std::cout << "The j must be non-negative" << std::endl;
       return 0;
   } 
   if (!IsVectorInteger(J123) && !IsVectorHalfInteger(J123))
   {
       std::cout << "The J vector has to be integer or half-integer" << std::endl;
       return 0;
   }
   if (!IsVectorInteger(M123) && !IsVectorHalfInteger(M123))
   {
       std::cout << "The J vector has to be integer or half-integer" << std::endl;
       return 0;
   }
   if ((IsVectorInteger(J123) && !IsVectorInteger(M123)) || (IsVectorHalfInteger(J123) && !IsVectorHalfInteger(M123)))
   {
       std::cout << "The J vector doesn't match to M vector" << std::endl;
       return 0;
   }


 // Selection rules
   if ( (j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ))//j3 out of interval
   {
      std::cout << "j3 out of interval " << std::endl;
      return 0;
   }
      
   if( abs(m1 + m2 + m3) > 0.0001 )//non-conserving angular momentum
   {
      std::cout << "non-conserving angular momentum " << std::endl;
      return 0;
   }
      
   if( abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3 )//m is larger than j
   {
      std::cout << "m is larger than j " << std::endl;
      return 0;
   }

   // Evaluation
   double firstPower = j1-j2-m3;
   double firstTerm = TMath::Power(-1, firstPower);
   double secondTerm = TMath::Power(TriangleCoefficient(j1, j2, j3), 0.5);
   
   double thirdTerm=TMath::Power(Factorial(j1+m1)*Factorial(j1-m1)
                         *Factorial(j2+m2)*Factorial(j2-m2)
                         *Factorial(j3+m3)*Factorial(j3-m3), 0.5);

   double t1 = j2 - m1 - j3;
   double t2 = j1 + m2 - j3;
   double t3 = j1 + j2 - j3;
   double t4 = j1 - m1;
   double t5 = j2 + m2;

   double tMin1 = TMath::Max(t1, t2);
   double tmin = TMath::Max(0., tMin1);
   double tmax = TMath::Min(t3, TMath::Min(t4, t5));
   double t;
   double i=tmin;
   double lastTerm = 0;
   while(i <= tmax)
   {
       t=i;
       double x=Factorial(t)*Factorial(t-t1)
                *Factorial(t-t2)*Factorial(t3-t)
                *Factorial(t4-t)*Factorial(t5-t);
       lastTerm += pow(-1, t)/x;         
       i++;
   }
   return firstTerm*secondTerm*thirdTerm*lastTerm;
}
     
double Anisotropy::TriangleCoefficient(double a, double b, double c)
{
   double triangleCoefficient = Factorial(a+b-c)*Factorial(a-b+c)*Factorial(-a+b+c)/Factorial(a+b+c+1);
   return triangleCoefficient;
}
double Anisotropy::Factorial(double x)
{
   if(x==0 || x==1)
      return 1;
   double fact=1;
   double iter = x;
   while(iter >1)
   {
      fact*=iter;
      iter--;
   }
   return fact;
}

bool Anisotropy::IsVectorHalfInteger(std::vector<double> data)
{
   for(int i = 0; i!=data.size(); ++i)
   {
      if(std::floor(2.*data.at(i)) != 2.*data.at(i))
         return false; 
   }
   return true;
}
     
bool Anisotropy::IsVectorInteger(std::vector<double> data)
{
   for(int i = 0; i!=data.size(); ++i)
   {
      if(std::floor(data.at(i)) != data.at(i))
         return false; 
   }
   return true;
}
