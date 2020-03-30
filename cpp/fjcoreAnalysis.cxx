#include "lcio.h"

#include "IO/LCReader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "fjcore.hh"
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
    TH1D hist("nJetsHist", "nJetsHist", 20, 0, 20);
    // create a jet definition: 
    // a jet algorithm with a given radius parameter
    //----------------------------------------------------------
    double R = 0.7;
	double sum = 0;
    fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, R);
    for (int nfiles=1; nfiles< argc; nfiles++) {
        string FILEN = argv[nfiles];
        cout << "Opening File " << nfiles << " "<<FILEN << endl;
        
        LCReader* lcReader = LCFactory::getInstance()->createLCReader();
        lcReader->open( FILEN );
        
        LCEvent* evt;
        //----------- the event loop -----------
        while ( (evt = lcReader->readNextEvent()) != 0 ) {
            vector<fjcore::PseudoJet> input_particles;
            LCCollection* col = evt->getCollection("PandoraPFOs");
            for (size_t i=0,N=col->getNumberOfElements(); i < N; ++i) {
                ReconstructedParticle * p = (ReconstructedParticle*) col->getElementAt(i);
                const double e = p->getEnergy();
                const double* m = p->getMomentum();
                input_particles.push_back(fjcore::PseudoJet(m[0], m[1], m[2], e)); 
            }
            // run the jet clustering with the above jet definition
            //----------------------------------------------------------
            fjcore::ClusterSequence clust_seq(input_particles, jet_def);
            // get the resulting jets ordered in pt
            //----------------------------------------------------------
            double ptmin = 1.0;
            vector<fjcore::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
            for (auto j : inclusive_jets) {
                //printf("%.8f %.8f %.8f %.8f\n", j.px(), j.py(), j.pz(), j.e());
                sum += j.px() + j.py() + j.pz() + j.e();
            }
            //printf("%zu %zu\n", input_particles.size(), inclusive_jets.size());
            nEvents++;
        }
        lcReader->close();
        delete lcReader;    
    }
	printf("%.8f\n", sum);
    return 0;
}

