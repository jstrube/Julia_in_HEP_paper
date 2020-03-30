#include "lcio.h"

#include "IO/LCReader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TFile.h"

#include "libFoxWo.hh"

#include <string>
#include <iostream>
#include <vector>

using namespace std;
using namespace lcio;


int main(int argc, char** argv ){

  // read file names from command line (only argument) 
  if( argc < 2) {
    cout << " usage:  cppFoxWo input-files" << endl << endl;
    exit(1);
  }
 
  vector <ReconstructedParticle *> particleList;
  
  
  
  int nEvents = 0;
  double h10,h20,h30,h40,thrust;
  
  
  TFile* rootfile =new TFile("foxwo.root","RECREATE");
  if ( rootfile->IsOpen() )
  {
      printf("ROOT output file opened successfully\n");
  }
  else
  {
      exit(-1);
  }
  
  TH1D hist1("H10", "H10", 50, 0.0, 1.0);
  TH1D hist2("H20", "H20", 50, 0.0, 1.0);
  TH1D hist3("H30", "H30", 50, 0.0, 1.0);
  TH1D hist4("H40", "H40", 50, 0.0, 1.0); 
  TH1D hist5("Thrust", "Thrust", 50, 0.5, 1.0);
  
  for (int nfiles=1; nfiles< argc; nfiles++)
  {
      string FILEN = argv[nfiles];
      cout << "Opening File " << nfiles << " "<<FILEN << endl;
      
      LCReader* lcReader = LCFactory::getInstance()->createLCReader();
      lcReader->open( FILEN );
      
      LCEvent* evt;
      
      
      //----------- the event loop -----------
      while ( (evt = lcReader->readNextEvent()) != 0 ) 
      {
          LCCollection* col = evt->getCollection("PandoraPFOs");
     
          
          for (size_t i=0,N=col->getNumberOfElements(); i < N; ++i) 
          {
              ReconstructedParticle * p = (ReconstructedParticle*) col->getElementAt(i);
              particleList.push_back(p);              
          }
          foxwo(particleList,h10,h20,h30,h40);
          Thrust(particleList,thrust);
          //cout <<"FOXWO RES" << h10<<" "<<h20<<endl;
          //printf("%.8f %.8f %.8f %.8f %.8f\n", h10, h20, h30, h40, thrust);
          hist1.Fill(h10);
          hist2.Fill(h20);
          hist3.Fill(h30);
          hist4.Fill(h40);
          hist5.Fill(thrust);
          particleList.clear();
          nEvents++;
      }
      lcReader->close();
      delete lcReader;
      
  }
  // -------- end of event loop -----------
  rootfile->Write();
  rootfile->Close();
  delete rootfile;

  return 0;
}

