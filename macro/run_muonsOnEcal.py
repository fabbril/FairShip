#!/usr/bin/env python
import ROOT,os,sys,getopt,time
import shipunit as u
import shipRoot_conf
from ShipGeoConfig import ConfigRegistry
#import ShipGeoConfig
decayVolumeLength    = 5000 # cm 
meanMomentum = 10.0 # GeV
meanTheta    = 0.01 # in radiants
mcEngine     = "TGeant4"
simEngine    = "Pythia8"  # "Genie" # Ntuple
nEvents      = 10
radius       = 200   # cm
inclusive    = False  # True = all processes if False only ccbar -> HNL
deepCopy     = False  # False = copy only stable particles to stack, except for HNL events
eventDisplay = False
inputFile    = None
theSeed      = int(10000 * time.time() % 10000000)

print "FairShip setup for",simEngine,"to produce",nEvents,"events"
if simEngine == "Ntuple" and not inputFile :
  print 'input file required if simEngine = Ntuple'
if simEngine == "Pythia6" and not inputFile :
  print 'pythia6 requires inputfile'
ROOT.gRandom.SetSeed(theSeed)  # this should be propagated via ROOT to Pythia8 and Geant4VMC
shipRoot_conf.configure()
#ship_geo = ShipGeoConfig.Config().loadpy("$FAIRSHIP/geometry/geometry_config.py")
ship_geo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py",muShieldDesign=2,targetOpt=5)

# Output file name
outFile ="geant_test_newmudec_g4_10.root"

# rm older files !!! 
os.system("rm params.root")
# Parameter file name
parFile="params.root"

# In general, the following parts need not be touched
# ========================================================================

# -----Timer--------------------------------------------------------
timer = ROOT.TStopwatch()
timer.Start()
# ------------------------------------------------------------------------

# -----Create simulation run----------------------------------------
run = ROOT.FairRunSim()
run.SetName(mcEngine)  # Transport engine
run.SetOutputFile(outFile)  # Output file
run.SetPythiaDecayer(ROOT.kTRUE)
# only with TGeant4?
#run.SetUserConfig("g4Config.C") # user configuration file default g4Config.C 
rtdb = run.GetRuntimeDb() 
# -----Create geometry----------------------------------------------
import shipDet_conf
shipDet_conf.configure(run,ship_geo)
# -----Create PrimaryGenerator--------------------------------------
primGen=ROOT.FairPrimaryGenerator()
myGen=ROOT.MuonBeamGen(13,1)
myGen.SetPMean(meanMomentum)
myGen.SetMeanTheta(meanTheta)
myGen.SetRZ(radius,-decayVolumeLength/2)
myGen.SetEleFrac(0.0)
primGen.AddGenerator(myGen)
run.SetGenerator(primGen)

# ------------------------------------------------------------------------

#---Store the visualization info of the tracks, this make the output file very large!!
#--- Use it only to display but not for production!
if eventDisplay: run.SetStoreTraj(ROOT.kTRUE)
else:            run.SetStoreTraj(ROOT.kFALSE)
# -----Initialize simulation run------------------------------------
run.Init()
fStack = ROOT.gMC.GetStack()
#fStack.SetEnergyCut(.3*u.MeV)
# ------------------------------------------------------------------------
if simEngine != "Genie" and simEngine != "Ntuple":  
# -----Runtime database---------------------------------------------
 kParameterMerged = ROOT.kTRUE
 parOut = ROOT.FairParRootFileIo(kParameterMerged)
 parOut.open(parFile)
 rtdb.setOutput(parOut)
 rtdb.saveOutput()   # for the moment, it blocks when using Genie, no idea why
 rtdb.printParamContexts()
# -----Start run----------------------------------------------------
run.Run(nEvents)
# ------------------------------------------------------------------------
run.CreateGeometryFile("geofile_full.root")  
# -----Finish-------------------------------------------------------
timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
print ' ' 
print "Macro finished succesfully." 
print "Output file is ",  outFile 
print "Parameter file is ",parFile
print "Real time ",rtime, " s, CPU time ",ctime,"s"

# ------------------------------------------------------------------------
