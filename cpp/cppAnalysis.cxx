#include "lcio.h"

#include "IO/LCReader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include <string>
#include <iostream>

using namespace std;
using namespace lcio;

int main(int argc, char** argv ){

  // read file names from command line (only argument) 
  if( argc < 2) {
    cout << " usage:  readmcparticles input-file" << endl << endl;
    exit(1);
  }
 

  int nEvents = 0;
  TH1D hist("massHist", "Z mass", 10, 85, 95);
  
  for (int nfiles=1; nfiles< argc; nfiles++)
  {
      string FILEN = argv[nfiles];
      cout << "Opening File " << nfiles << " "<<FILEN << endl;
      
      LCReader* lcReader = LCFactory::getInstance()->createLCReader();
      lcReader->open( FILEN );
      
      LCEvent* evt;
      
      
      //----------- the event loop -----------
      while ( (evt = lcReader->readNextEvent()) != 0 ) {
          TLorentzVector mu1;
          TLorentzVector mu2;
          LCCollection* col = evt->getCollection("PandoraPFOs");
          for (size_t i=0,N=col->getNumberOfElements(); i < N; ++i) {
              ReconstructedParticle * p = (ReconstructedParticle*) col->getElementAt(i);
              if (abs(p->getType()) != 13) continue;
              const double e = p->getEnergy();
              const double* m = p->getMomentum();
              if (e > mu1.E()) {
                  mu1.SetPxPyPzE(m[0], m[1], m[2], e);
              } else if (e > mu2.E()) {
                  mu2.SetPxPyPzE(m[0], m[1], m[2], e);
              }
          }
          hist.Fill((mu1+mu2).M());
          nEvents++;
      }
      lcReader->close();
      delete lcReader;
      
  }
  // -------- end of event loop -----------
  
  hist.Fit("gaus");
  TF1 *fitfunction = hist.GetFunction("gaus");
  // Extract the Z Mass
  double m= 0.0;
  double m_error= 0.0;
  if (fitfunction!=NULL)
  {
     m=fitfunction->GetParameter(1);
     m_error=fitfunction->GetParError(1);
  }
  cout << " read " << nEvents << "  <m> = " << m << "+/-" << m_error<< endl;
  
  return 0;
}

