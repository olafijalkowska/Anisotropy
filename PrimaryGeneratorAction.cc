#include <fstream>
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Anisotropy.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	SetUp();
	SetUpGammaAssymetry();
}

void PrimaryGeneratorAction::SetUp( void )
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun( n_particle );

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* gammaPD = particleTable->FindParticle("gamma");
  particleGun->SetParticleDefinition(gammaPD);
  particleGun->SetParticleEnergy( 500.0 * keV );
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

}

void PrimaryGeneratorAction::SetUpBetaAssymetry( void )
{
    betaEnergy = 2*MeV;//optionally, one may assume v/c=1 and set beta energy later
    double polarisation = 0.4; 
    int initialSpin = 1; 
    int deltaSpin = 1; 
    TransitionType transitionType = TransitionType::GT;
    
    double velocity = FindVelocity(betaEnergy);
    double assymetryFactor = FindAsymmertyFactor(transitionType, deltaSpin, initialSpin);
    betaAsymDistr = FindAsymDistrFunc(assymetryFactor,velocity, polarisation);
  
}

void PrimaryGeneratorAction::SetUpGammaAssymetry()
{     
   double I0=3;
   std::vector<double> mVal= {-3, -2, -1, 0, 1, 2, 3};
   std::vector<double> am = {0.002, 0.002, 0.004, 0.08, 0.1667, 0.3333, 0.5};  
   double Ii=2;
   double If=0;
   double L1=2;
   double L2=2; 
   double mixRatio = 0;

   Anisotropy* gammaAnisotropy = new Anisotropy(I0, mVal, am, Ii, If, L1, L2, mixRatio);
   gammaAsymDistr = gammaAnisotropy->  GetThetaDistr(); 
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent )
{
	//Position
	G4ThreeVector startPos( 0.0*cm, 0.0*cm, 0.0*cm );
	
	//Direction
	G4ThreeVector direction( 1.0, 0.0, 0.0 );
	GenerateAsymDirection(&direction, gammaAsymDistr);
	
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* gammaPD = particleTable->FindParticle("gamma");
	
	particleGun->SetParticleDefinition(gammaPD);
	//particleGun->SetParticleEnergy(betaEnergy);	
	particleGun->SetParticlePosition( startPos );
	particleGun->SetParticleMomentumDirection( direction );
	particleGun->GeneratePrimaryVertex(anEvent);
		
}



void  PrimaryGeneratorAction::GenerateAsymDirection(G4ThreeVector* direction, TF1* asymDistr)
{
   double cosTheta = ( G4UniformRand() - 0.5 ) * 2.0;
   double sinTheta = sqrt( 1.0 - cosTheta * cosTheta );
   double phi = asymDistr->GetRandom();
   //double phi = randomGen->Uniform (0, 2.*TMath::Pi());
   double randomXaim = cos(phi) * sinTheta;
   double randomYaim = sin(phi) * sinTheta;
   double randomZaim = cosTheta;
   *direction = G4ThreeVector (randomXaim, randomYaim, randomZaim);	   
}


TF1* PrimaryGeneratorAction::FindAsymDistrFunc(double asymFactor, double velocity, double polarisation)
{
   TF1 *asymDistrFunc = new TF1("asymDistrFunc","1+[0]*cos(x)",0,2.*TMath::Pi());
   asymDistrFunc->SetParameter(0,asymFactor*velocity*polarisation);
   return asymDistrFunc;
}

double PrimaryGeneratorAction::FindAsymmertyFactor(PrimaryGeneratorAction::TransitionType transitionType, int deltaI, int initialSpin)
{
   if(deltaI != 0 && deltaI != -1 && deltaI != 1)
   {
       std::cout << "FindAsymmertyFactor designed only for allowed beta transitions. "
                 << "Possible Delta I = -1, 0, 1" << std::endl;
   }
   double asymFact;
   switch(transitionType)
   {
       case TransitionType::FF: 
           asymFact = 0.;   
           break;
       case TransitionType::GT: 
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

double PrimaryGeneratorAction::FindVelocity(double energyinMeV)
{
//Return electron velocity in c unit (beta).
   double electronMass = 0.511;//in MeV unit
   double beta = pow((1-(electronMass*electronMass/(energyinMeV*energyinMeV))), 0.5);
   return beta;
}
    
    
