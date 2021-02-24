//w kolejnym kroku wrzuc zliczenia z histogramu
enum TransitionType { FF,GT };


double FindAsymmertyFactor(TransitionType transitionType, int deltaI=0, int primarySpin = 0)
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
              asymFact = (double)primarySpin/(primarySpin+1.);
           }
           if(deltaI == 0)
           {
              asymFact = -1./(primarySpin+1.);
           }
           break;
   }
   return asymFact;
}

double FindAsymDistr(double theta, double asymFactor, double velocity, double polarisation)
{
   return 1.+asymFactor*velocity*polarisation*cos(theta);
}

void PlotPolar()
{
   //TCanvas * CPol = new TCanvas("CPol","TGraphPolar Example",500,500);
   double asymFactor=FindAsymmertyFactor(GT, 1, 1);
   std::cout << "asym fact " << asymFactor << std::endl;
   double velocity = 1.;//in c units
   double polarisation = 0.4;//example
   int nrOfPoints=360;
   Double_t theta[nrOfPoints];
   Double_t radiusNoPolar[nrOfPoints];
   Double_t radiusPolar[nrOfPoints];
   //Double_t radiusPolar[nrOfPoints];
 
   for (int i=0; i<nrOfPoints; i++) {
      theta[i]   = (double)(i+1)*(2.*TMath::Pi()/nrOfPoints);
      radiusNoPolar[i] = 0.25;
      radiusPolar[i]  = FindAsymDistr(theta[i], asymFactor, velocity, polarisation);
   }
 
   TGraphPolar * noPolarGraph = new TGraphPolar(nrOfPoints, theta, radiusNoPolar);
   noPolarGraph->SetTitle("Beta distribution tests");
   noPolarGraph->SetLineColor(kRed);
   noPolarGraph->SetLineWidth(3);
   //noPolarGraph->SetMaxRadial (1.5);	
   TGraphPolar * polarGraph = new TGraphPolar(nrOfPoints, theta, radiusPolar);
   polarGraph->SetLineColor(kBlue);
   polarGraph->SetLineWidth(3);
   
   noPolarGraph->Draw();
   polarGraph->Draw("same");
   // Update, otherwise GetPolargram returns 0
   //CPol->Update();
   //noPolarGraph->GetPolargram()->SetToRadian();
 
   //return CPol;
}

