# example for accessing smeared hits and fitted tracks
import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry

PDG = ROOT.TDatabasePDG.Instance()
inputFile = 'geant_mudec_g4_newgeo_10k_rec.root'
inputFile2 = 'ana_mudec_g4_newgeo_10k.root'
dy = 10.

try:
        opts, args = getopt.getopt(sys.argv[1:], "f:A:Y:i", [])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter file name'
        sys.exit()
for o, a in opts:
        print o, a
        if o in ("-f"):
            inputFile = a
        if o in ("-Y"): 
            dy = float(a)

if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')
  try:
    dy = float( tmp[1]+'.'+tmp[2] )
  except:
    dy = None
else:
 inputFile = 'ship.'+str(dy)+'Pythia8-TGeant4_rec.root'
 inputFile = 'geant_mudec_g4_newgeo_10k_rec.root'
print dy
print inputFile
  
f     = ROOT.TFile(inputFile)
f2    = ROOT.TFile(inputFile2)
sTree = f.cbmsim
sTree2 = f2.decay

# init geometry and mag. field
ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py", Yheight = dy )
# -----Create geometry----------------------------------------------
import shipDet_conf
run = ROOT.FairRunSim()
modules = shipDet_conf.configure(run,ShipGeo)

tgeom = ROOT.TGeoManager("Geometry", "Geane geometry")
geofile = inputFile.replace('ship.','geofile_full.').replace('_rec.','.')
gMan  = tgeom.Import(geofile)
geoMat =  ROOT.genfit.TGeoMaterialInterface()
ROOT.genfit.MaterialEffects.getInstance().init(geoMat)

bfield = ROOT.genfit.BellField(ShipGeo.Bfield.max, ShipGeo.Bfield.z,2)
fM = ROOT.genfit.FieldManager.getInstance()
fM.init(bfield)

ev     = ut.PyListOfLeaves()
leaves = sTree.GetListOfLeaves()
names  = ut.setAttributes(ev, leaves)

h = {}
ut.bookHist(h,'delPOverP','delP / P',100,0.,50.,100,-0.2,0.2)
ut.bookHist(h,'delPOverP2','delP / P chi2<25',100,0.,50.,100,-0.2,0.2)
ut.bookHist(h,'chi2','chi2 after trackfit',100,0.,1000.)
#ut.bookHist(h,'IP','Impact Parameter',100,0.,10.)
#ut.bookHist(h,'Doca','Doca between two tracks',100,0.,10.)
#ut.bookHist(h,'IP0','Impact Parameter to target',100,0.,100.)
#ut.bookHist(h,'IP0/mass','Impact Parameter to target vs mass',100,0.,2.,100,0.,100.)
#ut.bookHist(h,'HNL','reconstructed Mass',100,0.,2.)
ut.bookHist(h,'meas','number of measurements',25,-0.5,24.5)
ut.bookHist(h,'delEOverP','delE / P',80,0.,0.36)
ut.bookHist(h,'P_ele','momentum of selected electrons',50,0.,50.)
ut.bookHist(h,'posZ_ele','Z position of selected electrons',50,1500.,1510.)
ut.bookHist(h,'posXY_ele','X Y position of selected electrons',50,-300.,300.,50,-300.,300.)
ut.bookHist(h,'P_all','momentum of all electrons',50,0.,50.)
ut.bookHist(h,'posZ_all','Z position of all electrons ',50,1500.,1510.)
ut.bookHist(h,'posXY_all','X Y position of all electrons',50,-300.,300.,50,-300.,300.)
ut.bookHist(h,'pdgCode','PDG codes',40,-20.,20)

def fitSingleGauss(x,ba=None,be=None):
    name    = 'myGauss_'+x 
    myGauss = h[x].GetListOfFunctions().FindObject(name)
    if not myGauss:
       if not ba : ba = h[x].GetBinCenter(1) 
       if not be : be = h[x].GetBinCenter(h[x].GetNbinsX()) 
       bw    = h[x].GetBinWidth(1) 
       mean  = h[x].GetMean()
       sigma = h[x].GetRMS()
       norm  = h[x].GetEntries()*0.3
       myGauss = ROOT.TF1(name,'[0]*'+str(bw)+'/([2]*sqrt(2*pi))*exp(-0.5*((x-[1])/[2])**2)+[3]',4)
       myGauss.SetParameter(0,norm)
       myGauss.SetParameter(1,mean)
       myGauss.SetParameter(2,sigma)
       myGauss.SetParameter(3,1.)
       myGauss.SetParName(0,'Signal')
       myGauss.SetParName(1,'Mean')
       myGauss.SetParName(2,'Sigma')
       myGauss.SetParName(3,'bckgr')
    h[x].Fit(myGauss,'','',ba,be) 


def makePlots():
   ut.bookCanvas(h,key='strawanalysis',title='Distance to wire and mean nr of hits',nx=1200,ny=600,cx=2,cy=1)
   cv = h['strawanalysis'].cd(1)
   h['disty'].Draw()
   h['distu'].Draw('same')
   h['distv'].Draw('same')
   cv = h['strawanalysis'].cd(2)
   h['meanhits'].Draw()
   ut.bookCanvas(h,key='fitresults',title='Fit Results',nx=1600,ny=1200,cx=2,cy=2)
   cv = h['fitresults'].cd(1)
   h['delPOverP'].Draw('box')
   cv = h['fitresults'].cd(2)
   cv.SetLogy(1)
   h['chi2'].Draw()
   cv = h['fitresults'].cd(3)
   h['delPOverP_proj'] = h['delPOverP'].ProjectionY()
   ROOT.gStyle.SetOptFit(11111)
   h['delPOverP_proj'].Draw()
   h['delPOverP_proj'].Fit('gaus')
   cv = h['fitresults'].cd(4)
   h['delPOverP2_proj'] = h['delPOverP2'].ProjectionY()
   h['delPOverP2_proj'].Draw()
   h['delPOverP2_proj'].Fit('gaus')
   h['fitresults'].Print('fitresults.gif')
   ut.bookCanvas(h,key='fitresults2',title='Fit Results',nx=1600,ny=1200,cx=2,cy=2)
   print 'finished with first canvas'
   cv = h['fitresults2'].cd(1)
   h['P_ele'].Draw()
   cv = h['fitresults2'].cd(2)
   h['posXY_ele'].Draw()
   #h['IP0'].Draw()
   cv = h['fitresults2'].cd(3)
   h['delEOverP'].Draw()
   h['delEOverP'].Fit('gaus','', '', 0.15, 0.198)
   # h['HNL'].Draw()
   #h['HNL'].Fit('gaus')
   cv = h['fitresults2'].cd(4)
   h['posZ_ele'].Draw()
   #h['IP0/mass'].Draw('box')
   h['fitresults2'].Print('fitresults2.gif')
   ut.bookCanvas(h,key='fitresults3',title='Fit Results',nx=1600,ny=1200,cx=2,cy=2)
   print 'finished with second canvas'
   cv = h['fitresults3'].cd(1)
   h['P_all'].Draw()
   cv = h['fitresults3'].cd(2)
   h['posXY_all'].Draw()
   cv = h['fitresults3'].cd(3)
   h['posZ_all'].Draw()
   cv = h['fitresults3'].cd(4)
   h['pdgCode'].Draw()
   h['fitresults3'].Print('fitresults3.gif')
   print 'finished making plots'


def myVertex(t1,t2,PosDir):
 # closest distance between two tracks
   V=0
   for i in range(3):   V += PosDir[t1][1](i)*PosDir[t2][1](i)
   S1=0
   for i in range(3):   S1 += (PosDir[t1][0](i)-PosDir[t2][0](i))*PosDir[t1][1](i)
   S2=0
   for i in range(3):   S2 += (PosDir[t1][0](i)-PosDir[t2][0](i))*PosDir[t2][1](i)
   l = (S2-S1*V)/(1-V*V)
   x2 = PosDir[t2][0](0)+l*PosDir[t2][1](0)
   y2 = PosDir[t2][0](1)+l*PosDir[t2][1](1)
   z2 = PosDir[t2][0](2)+l*PosDir[t2][1](2)
   x1 = PosDir[t1][0](0)+l*PosDir[t1][1](0)
   y1 = PosDir[t1][0](1)+l*PosDir[t1][1](1)
   z1 = PosDir[t1][0](2)+l*PosDir[t1][1](2)
   dist = ROOT.TMath.Sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
   return (x1+x2)/2.,(y1+y2)/2.,(z1+z2)/2.,dist


# start event loop
def myEventLoop():
 nEvents = sTree.GetEntries()
 nEvents2 = sTree2.GetEntries()
 count = 0
 ntrack = 0
 notfitt = 0
 for n in range(nEvents):
  count+=1
  rc = sTree.GetEntry(n)
  rc2 = sTree2.GetEntry(n)
  eloss = sTree2.ElEloss
  eid = sTree2.eid
# make some straw hit analysis
  hitlist = {}
  for ahit in sTree.strawtubesPoint:
     detID = ahit.GetDetectorID()
     top = ROOT.TVector3()
     bot = ROOT.TVector3()
     modules["Strawtubes"].StrawEndPoints(detID,bot,top)
     dw  = ahit.dist2Wire()
     if abs(top.y())==abs(bot.y()): h['disty'].Fill(dw)
     if abs(top.y())>abs(bot.y()): h['distu'].Fill(dw)
     if abs(top.y())<abs(bot.y()): h['distv'].Fill(dw)
#
     trID = ahit.GetTrackID()
     if not trID < 0 :
      if hitlist.has_key(trID):  hitlist[trID]+=1
      else:  hitlist[trID]=1
  for tr in hitlist:  h['meanhits'].Fill(hitlist[tr])
  key = 0
  fittedTracks = {}
  for atrack in ev.FitTracks.GetObject():
   fitStatus   = atrack.getFitStatus()
   nmeas = atrack.getNumPoints()
   h['meas'].Fill(nmeas)
   if not fitStatus.isFitConverged() :
     notfitt += 1
    #print "notfitt =", notfitt
    #continue
   fittedTracks[key] = atrack
# needs different study why fit has not converged, continue with fitted tracks
   chi2        = fitStatus.getChi2()/nmeas
   fittedState = atrack.getFittedState()
   h['chi2'].Fill(chi2)
   h['measVSchi2'].Fill(atrack.getNumPoints(),chi2)
   P = fittedState.getMomMag()
   mcPartKey = sTree.fitTrack2MC[key]
   mcPart    = sTree.MCTrack[mcPartKey]
   if not mcPart : continue
   Ptruth    = mcPart.GetP()
   delPOverP = (Ptruth - P)/Ptruth
   h['delPOverP'].Fill(Ptruth,delPOverP)
   if chi2>25: continue
   h['delPOverP2'].Fill(Ptruth,delPOverP)
# try measure impact parameter
   trackDir = fittedState.getDir()
   trackPos = fittedState.getPos()
   vx = ROOT.TVector3()
   mcPart.GetStartVertex(vx)
   t = 0
   for i in range(3):   t += trackDir(i)*(vx(i)-trackPos(i)) 
   dist = 0
   for i in range(3):   dist += (vx(i)-trackPos(i)-t*trackDir(i))**2
   dist = ROOT.TMath.Sqrt(dist)
#h['IP'].Fill(dist)
   key+= 1
   h['pdgCode'].Fill(mcPart.GetPdgCode())
   h['posZ_all'].Fill(trackPos(2))
   h['posXY_all'].Fill(trackPos(0), trackPos(1))
   h['P_all'].Fill(P)
                       
   #if mcPartKey == eid :
   h['delEOverP'].Fill(eloss/P)
   h['P_ele'].Fill(P)
   h['posZ_ele'].Fill(trackPos(2))
   h['posXY_ele'].Fill(trackPos(0), trackPos(1))

# ---
# go for 2-track combinations
  if len(fittedTracks) == 2:
     LV  = {}
     PosDir = {} 
     for tr in fittedTracks:
      xx  = fittedTracks[tr].getFittedState()
      PosDir[tr] = [xx.getPos(),xx.getDir()]
     keys = fittedTracks.keys()
     t1,t2 = keys[0],keys[1] 
     xv,yv,zv,doca = myVertex(t1,t2,PosDir)
     h['Doca'].Fill(dist)  
     HNLPos = ROOT.TVector3(xv,yv,zv)
     for tr in fittedTracks:
      xx  = fittedTracks[tr].getFittedState()
      # make a new rep
      rep = ROOT.genfit.RKTrackRep(xx.getPDG())
      pos = xx.getPos()
      mom = xx.getMom()
      state = ROOT.genfit.StateOnPlane(rep)
      rep.setPosMom(state, pos, mom)
      origPlane = state.getPlane()
      origState = ROOT.genfit.StateOnPlane(state)
      point = ROOT.TVector3(1.,5.,3.)
      try:
       rep.extrapolateToPoint(state, HNLPos, False)
      except:
        print 'extrap did not worked'
      LV[tr] = ROOT.TLorentzVector()
      mass = PDG.GetParticle(xx.getPDG()).Mass()
      mom = rep.getMom(state)  
      E = ROOT.TMath.Sqrt( mass*mass + mom.Mag2() )
      LV[tr].SetPxPyPzE(mom.x(),mom.y(),mom.z(),E)
     HNL = LV[t1]+LV[t2]
     tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
     t = 0
     for i in range(3):   t += HNL(i)/HNL.P()*(tr(i)-HNLPos(i)) 
     dist = 0
     for i in range(3):   dist += (tr(i)-HNLPos(i)-t*HNL(i)/HNL.P())**2
     dist = ROOT.TMath.Sqrt(dist)
     #h['IP0'].Fill(dist)
     #h['IP0/mass'].Fill(HNL.M(),dist)
#h['HNL'].Fill(HNL.M())


def access2SmearedHits():
 key = 0
 for ahit in ev.SmearedHits.GetObject():
   print ahit[0],ahit[1],ahit[2],ahit[3],ahit[4],ahit[5],ahit[6]
   # follow link to true MCHit
   mchit   = TrackingHits[key]
   mctrack =  MCTracks[mchit.GetTrackID()]
   print mchit.GetZ(),mctrack.GetP(),mctrack.GetPdgCode()
   key+=1

myEventLoop()
makePlots()


