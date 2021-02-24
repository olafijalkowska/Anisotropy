//script calculate random direction of beta particle 
/*Input parameters:
beta energy in MeV unit
polarisation
spin of initial level
delta spin (0 or 1)
decay type:
decayType=0 - FF decay
decayType=1 - GT decay
example:
root -l "RandDistrWithPolar.cxx(2, 0.4, 1, 1)"
*/
enum TransitionType { FF,GT };

double FindAsymmertyFactor(TransitionType transitionType, int deltaI=0, int initialSpin = 0);
TF1* FindAsymDistrFunc(double asymFactor, double velocity, double polarisation);
double FindVelocity(double energyinMeV);
std::vector<double> FindRandomDir(TF1 *asymDistrFunc, TRandom* randomGen);


void RandDistrWithPolar(double energyinMeV, double polarisation, int initialSpin, int deltaSpin, int decayType=1)
{
   double asymFactor;
   if(decayType == 0)
       asymFactor=FindAsymmertyFactor(FF, deltaSpin, initialSpin);
   else if(decayType == 1)
       asymFactor=FindAsymmertyFactor(GT, deltaSpin, initialSpin);
   else
   {
       std::cout << " " << std::endl;
       
   }    
   double velocity = FindVelocity(energyinMeV);
   TF1* asymFunct = FindAsymDistrFunc(asymFactor, velocity, polarisation);
   int nrOfPoints = 10000;
   double x[nrOfPoints];
   double y[nrOfPoints];
   double z[nrOfPoints];
   TRandom* randomGen = new TRandom();
   TH1F* phiHisto = new TH1F("phiHisto", "phiHisto", 1000, 0, 2.*TMath::Pi());
   for(int i=0; i!=nrOfPoints; ++i)
   { 
      std::vector<double> distr = FindRandomDir(asymFunct, randomGen);
      x[i]=distr.at(0);
      y[i]=distr.at(1);
      z[i]=distr.at(2);
      double phi=distr.at(3);                   
      phiHisto->Fill(phi);
   }
   gStyle->SetOptStat(0);
   //phiHisto->Draw();
   TGraph2D *directionGraph = new TGraph2D(nrOfPoints, x, y, z);
   directionGraph->SetMarkerStyle(6);
   directionGraph->Draw("p0");
}




double FindAsymmertyFactor(TransitionType transitionType, int deltaI=0, int initialSpin = 0)
{
   if(deltaI != 0 && deltaI != -1 && deltaI != 1)
   {
       std::cout << "FindAsymmertyFactor designed only for allowed beta transitions. "
                 << "Possible Delta I = -1, 0, 1" << std::endl;
   }
   double asymFact;
   switch(transitionType)
   {
       case FF: 
           asymFact = 0.;   
           break;
       case GT: 
           if(deltaI == -1)
               asymFact = -1.;
           if(deltaI == 1)
           {
              asymFact = (double)initialSpin/(initialSpin+1.);
           }
           if(deltaI == 0)
           {
              asymFact = -1./(initialSpin+1.);
           }
           break;
   }
   return asymFact;
}


TF1* FindAsymDistrFunc(double asymFactor, double velocity, double polarisation)
{
   TF1 *asymDistrFunc = new TF1("asymDistrFunc","1+[0]*cos(x)",0,2.*TMath::Pi());
   asymDistrFunc->SetParameter(0,asymFactor*velocity*polarisation);
   return asymDistrFunc;
}

double FindVelocity(double energyinMeV)
{
//Return electron velocity in c unit (beta).
   double electronMass = 0.511;//in MeV unit
   double beta = pow((1-(electronMass*electronMass/(energyinMeV*energyinMeV))), 0.5);
   return beta;
}

std::vector<double> FindRandomDir(TF1 *asymDistrFunc, TRandom* randomGen)
{
   
   double cosTheta = randomGen->Uniform (-1., 1);
   double sinTheta = sqrt( 1.0 - cosTheta * cosTheta );
   double phi = asymDistrFunc->GetRandom();
   //double phi = randomGen->Uniform (0, 2.*TMath::Pi());
   double randomXaim = cos(phi) * sinTheta;
   double randomYaim = sin(phi) * sinTheta;
   double randomZaim = cosTheta;
   std::vector<double> direction{randomXaim, randomYaim, randomZaim, phi};	
   return direction;
}

//zrób historham rozkładu, wylosuj z histograma phi i narysuj x,y,z wylosowanej liczby
