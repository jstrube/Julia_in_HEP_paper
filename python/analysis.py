# -*- coding: utf-8 -*-
from pyLCIO import IOIMPL
import sys
import ROOT
import time
ROOT.gROOT.SetBatch()

hist = ROOT.TH1D("invMass", "Z mass", 10, 85, 95)
iEvent = 0
start = time.time()
for file in range(1,len(sys.argv)):
    print("Opening ",sys.argv[file])
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(sys.argv[file])
    
    for event in reader:
        firstmu = ROOT.TLorentzVector()
        secondmu = ROOT.TLorentzVector()
        pfos = event.getCollection("PandoraPFOs")
        for particle in pfos:
            if abs(particle.getType()) != 13:
                continue
            energy = particle.getEnergy()
            p = particle.getMomentum()
            if energy > firstmu.E():
                firstmu.SetPxPyPzE(p[0], p[1], p[2], energy)
            elif energy > secondmu.E():
                secondmu.SetPxPyPzE(p[0], p[1], p[2], energy)
            if firstmu.E() == 0 or secondmu.E() == 0:
                continue
        hist.Fill((firstmu + secondmu).M())
        iEvent += 1
                    
print("execution of the event loop took", time.time()-start, "seconds")
print "looped over ", iEvent, " events"
hist.Fit("gaus")
fitfunction =  hist.GetFunction("gaus")
m=fitfunction.GetParameter(1);
m_err=fitfunction.GetParError(1);
print m, " ", m_err

