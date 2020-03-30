# -*- coding: utf-8 -*-
from pyLCIO import IOIMPL
import sys
import ROOT
import time
from FoxWo import foxWo, Thrust
ROOT.gROOT.SetBatch()

hist = ROOT.TH1D("invMass", "Z mass", 10, 85, 95)
iEvent = 0
start = time.time()
for file in range(1,len(sys.argv)):
    print("Opening ",sys.argv[file])
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(sys.argv[file])
    
    for event in reader:
        if iEvent >= 2339:
            break
        if iEvent % 100 == 0:
            print("Time for", iEvent, "events:", time.time()-start, "seconds")
        firstmu = ROOT.TLorentzVector()
        secondmu = ROOT.TLorentzVector()
        pfos = event.getCollection("PandoraPFOs")
        foxWo(pfos)
        Thrust(pfos)	
        iEvent += 1
                    
print("execution of the event loop took", time.time()-start, "seconds")
print "looped over ", iEvent, " events"

