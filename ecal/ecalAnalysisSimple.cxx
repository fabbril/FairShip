#include "ecalAnalysisSimple.h"

#include "FairRootManager.h"

#include "ecalStructure.h"
#include "ecalCell.h"
#include "ecalPoint.h"
#include "ShipMCTrack.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <iostream>
#include <fstream>
#include <list>

using namespace std;

/** Loop procedure **/
void ecalAnalysisSimple::Exec(Option_t* option)
{
    Int_t n=fTracks->GetEntries();
    Int_t npl=fPLTracks->GetEntries();
    Int_t i;
    int trkID = -2;
    ecalPoint* t;
    ecalPoint* tl;
    ecalCell* cell;
    ecalCell* mcell;
    TVector3 m;
    list<ecalCell*> cells;
    list<ecalCell*>::const_iterator p;
    //  TVector3 m1;
    fEv++;
    InitTree();
    
    //
    // Get MCTrack information
    //    
    Int_t nMC=fMCTracks->GetEntries();
    ShipMCTrack* mutrk = (ShipMCTrack*) fMCTracks->At(3);
    fspdg = mutrk->GetPdgCode();
    fsx = mutrk->GetStartX();
    fsy = mutrk->GetStartY();
    fsz = mutrk->GetStartZ();
    fspx = mutrk->GetPx();
    fspy = mutrk->GetPy();
    fspz = mutrk->GetPz();
    fdpdg=0; fdx=0.; fdy=0.; fdz=-3000.; fdpx=0.; fdpy=0.; fdpz=0.; fde=0., fdtime=-10., fdecay=false, isDecay =false;
    eDecTrackId=0;
    ShipMCTrack* etrk = (ShipMCTrack*) fMCTracks->At(0);
    fdpdg = etrk->GetPdgCode();
    fdx = etrk->GetStartX();
    fdy = etrk->GetStartY();
    fdz = etrk->GetStartZ();
    fdpx = etrk->GetPx();
    fdpy = etrk->GetPy();
    fdpz = etrk->GetPz();
    fde = etrk->GetEnergy();
    fdtime = etrk->GetStartT();
    
    fElEloss=0.0;
    cout << "EcalPointLite: npl = " << npl << endl;
    Int_t n_matching_point = 0;
    for (int ipl = 0; ipl < npl; ipl++) {
        tl=(ecalPoint*)fPLTracks->At(ipl);
        cout <<"tl->GetTrackID() = " << tl->GetTrackID() << endl;
        if(tl->GetTrackID() == 0){
            // cout << "Match!" << endl;
            n_matching_point++;
            fElEloss += tl->GetEnergyLoss();
        }
    }
    cout << "n matching point:" << n_matching_point << endl;
    fTree->Fill();
    
    cout << "EcalPoint: n = " << n << endl;
    for(i=0;i<n;i++) {
        isDecay = false;
        t=(ecalPoint*)fTracks->At(i);
        //cout<< "Ecal: trackID =" <<t->GetTrackID()<< " PgdCode = " << t->GetPdgCode() << endl;
        if(t->GetTrackID() == 0) {
            isDecay = true;
            cout<< "Ecal Point! i =" << i << endl;
        }
        fX=t->GetX();
        fY=t->GetY();
        t->Momentum(m);
        fP=m.Mag();
        fPX=m.Px();
        fPY=m.Py();
        fPZ=m.Pz();
        
        //    m1=m.Unit();
        cell=fStr->GetCell(fX, fY);
        if (!cell) {
            cout<< "EcalPoint: W no cell was found"<< endl;
            continue;
        }
        mcell=cell;
        cell->GetNeighborsList(cells);
        for(p=cells.begin();p!=cells.end();++p)
            if ((*p)->GetEnergy()>mcell->GetEnergy())
                mcell=(*p);
        
        mcell->GetNeighborsList(cells);
        for(p=cells.begin();p!=cells.end();++p)
            if ((*p)->GetEnergy()>mcell->GetEnergy())
                break;
        
        if (p!=cells.end()) {
            cout<< "EcalPoint: W a neighbor cell with higher energy was found"<< endl;
            continue;
        }
        
        fCX=mcell->GetCenterX();
        fCY=mcell->GetCenterY();
        fCE=mcell->GetEnergy();
        fCellNum=mcell->GetCellNumber();
        fADC=mcell->ADC();
        fOE=fCE;
        for(p=cells.begin();p!=cells.end();++p)
            fOE+=(*p)->GetEnergy();
        
        eTree->Fill();
    }
    
    
}

void ecalAnalysisSimple::InitTree()
{
    if (eTree) return;
    eTree=new TTree("calib", "calib");
    eTree->Branch("px", &fPX, "px/D");
    eTree->Branch("py", &fPY, "py/D");
    eTree->Branch("pz", &fPZ, "pz/D");
    eTree->Branch("p" , &fP , "p/D");
    eTree->Branch("x" , &fX , "x/D");
    eTree->Branch("y" , &fY , "y/D");
    eTree->Branch("cx", &fCX, "cx/D");
    eTree->Branch("cy", &fCY, "cy/D");
    eTree->Branch("ce", &fCE, "ce/D");
    eTree->Branch("oe", &fOE, "oe/D");
    eTree->Branch("ev", &fEv, "ev/I");
    eTree->Branch("cn", &fCellNum, "cn/I");
    eTree->Branch("adc", &fADC, "adc/I");
    eTree->Branch("EcalD", &isDecay, "isDecay/B");
    
    if (fTree) return;
    fTree=new TTree("decay", "decay");
    fTree->Branch("Spdg", &fspdg, "spdg/I");
    fTree->Branch("SX", &fsx, "sx/D");
    fTree->Branch("SY", &fsy, "sy/D");
    fTree->Branch("SZ", &fsz, "sz/D");
    fTree->Branch("SPX", &fspx, "spx/D");
    fTree->Branch("SPY", &fspy, "spy/D");
    fTree->Branch("SPZ", &fspz, "spz/D");
    fTree->Branch("Dpdg", &fdpdg, "dpdg/I");
    fTree->Branch("DX", &fdx, "dx/D");
    fTree->Branch("DY", &fdy, "dy/D");
    fTree->Branch("DZ", &fdz, "dz/D");
    fTree->Branch("DPX", &fdpx, "dpx/D");
    fTree->Branch("DPY", &fdpy, "dpy/D");
    fTree->Branch("DPZ", &fdpz, "dpz/D");
    fTree->Branch("DE", &fde, "de/D");
    fTree->Branch("DTIME", &fde, "dtime/D");
    fTree->Branch("Ddecay", &fdecay, "ddecay/B");
    fTree->Branch("eid", &eDecTrackId, "eid/I");
    fTree->Branch("mid", &fmotherid, "mid/I");
    fTree->Branch("ElEloss" , &fElEloss , "ElEloss/D");
    
}

ecalAnalysisSimple::ecalAnalysisSimple(const char* name, const Int_t iVerbose)
: FairTask(name, iVerbose),
fTree(NULL),
eTree(NULL),
fX(0.),
fY(0.),
fCX(0.),
fCY(0.),
fP(0.),
fCE(0.),
fOE(0.),
fPX(0.),
fPY(0.),
fPZ(0.),
fEv(0),
fCellNum(0),
fADC(0),
fspdg(0),
fsx(0.),
fsy(0.),
fsz(0.),
fspx(0.),
fspy(0.),
fspz(0.),
fdpdg(0),
fdx(0.),
fdy(0.),
fdz(0.),
fdpx(0.),
fdpy(0.),
fdpz(0.),
fde(0.),
fdtime(0.),
fElEloss(0.),
fdecay(false),
isDecay(false),
fmotherid(3),
fStr(NULL),
fTracks(NULL),
fPLTracks(NULL),
fMCTracks(NULL)
{
}

ecalAnalysisSimple::ecalAnalysisSimple()
: FairTask(),
fTree(NULL),
eTree(NULL),
fX(0.),
fY(0.),
fCX(0.),
fCY(0.),
fP(0.),
fCE(0.),
fOE(0.),
fPX(0.),
fPY(0.),
fPZ(0.),
fEv(0),
fCellNum(0),
fADC(0),
fspdg(0),
fsx(0.),
fsy(0.),
fsz(0.),
fspx(0.),
fspy(0.),
fspz(0.),
fdpdg(0),
fdx(0.),
fdy(0.),
fdz(0.),
fdpx(0.),
fdpy(0.),
fdpz(0.),
fde(0.),
fdtime(0.),
fElEloss(0.),
fdecay(false),
isDecay(false),
fmotherid(3),
fStr(NULL),
fTracks(NULL),
fPLTracks(NULL),
fMCTracks(NULL)
{
}

/** Initing routine **/
InitStatus ecalAnalysisSimple::Init()
{
    FairRootManager* fManager=FairRootManager::Instance();
    fStr=(ecalStructure*)fManager->GetObject("EcalStructure");
    if (!fStr)
    {
        Fatal("Init()", "Can't find calorimeter structure. ");
        return kFATAL;
    }
    
    // EcalPoint fTracks
    fTracks=(TClonesArray*)fManager->GetObject("EcalPoint");
    if (!fTracks)
    {
        Fatal("Init()", "Can't find array of reconstructed tracks. ");
        return kFATAL;
    }
    // EcalPointLite  fPLTracks
    
    fPLTracks=(TClonesArray*)fManager->GetObject("EcalPointLite");
    if (!fPLTracks)
    {
        Fatal("Init()", "Can't find array of reconstructed tracks in EcalPointLite. ");
        return kFATAL;
    }
    // MC fMCTracks
    fMCTracks=(TClonesArray*)fManager->GetObject("MCTrack");
    if (!fMCTracks)
    {
        Fatal("Init()", "Can't find array of MC tracks. ");
        return kFATAL;
    }
    return kSUCCESS;
}

/** Finishing routine **/
void ecalAnalysisSimple::Finish()
{
    if (fTree)
        fTree->Write();
    if (eTree)
        eTree->Write();
}

ClassImp(ecalAnalysisSimple)
