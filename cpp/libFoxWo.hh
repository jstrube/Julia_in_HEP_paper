#pragma once

#include "lcio.h"

#include "IO/LCReader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <vector>

int foxwo(std::vector <EVENT::ReconstructedParticle *> & particles,double & h10,double & h20, double & h30, double & h40);

//a simple implementation avoiding vectors should make it simpler to port to Python and Julia
int CrossProduct(double **v,int i, int j, double * result);

int Thrust(std::vector <EVENT::ReconstructedParticle *> & particles,double & thrust);

