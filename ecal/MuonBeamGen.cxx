/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
// -------------------------------------------------------------------------
// -----            MuonBeamGen source file                        -----
// -----          Created 09/09/04  by Yu.Kharlov
// -------------------------------------------------------------------------

/* $Id: MuonBeamGen.cxx,v 1.4 2006/07/18 09:28:06 prokudin Exp $ */

/* History of cvs commits:
 *
 * $Log: MuonBeamGen.cxx,v $
 * Revision 1.4  2006/07/18 09:28:06  prokudin
 * Should be * instead /
 *
 * Revision 1.3  2006/07/14 11:23:57  kharlov
 * Add protection for simultaneously set ranges; split vertex and kinematics ranges
 *
 * Revision 1.2  2006/03/29 16:25:50  kharlov
 * New functionality added
 *
 */

#include "MuonBeamGen.h"
#include "FairPrimaryGenerator.h"

#include "TRandom.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TGenPhaseSpace.h"


// ------------------------------------------------------------------------
MuonBeamGen::MuonBeamGen()
:  FairGenerator(),
fPDGType(13),
fPMean(0),fThetaMean(0),fR(0),fZ(0),fEleFrac(0),
fDebug(0)
{
  // Default constructor
}

// ------------------------------------------------------------------------
MuonBeamGen::MuonBeamGen(Int_t pdgid, Int_t mult,Int_t debug) :
FairGenerator(),
fPDGType(13),
fPMean(0),fThetaMean(0),fR(0),fZ(0),fEleFrac(0)
{
    fDebug=debug;
 
    // Constructor. Set default kinematics limits
    //  SetPhiRange  ();
}
// -------------------------------------------------------------------------
Bool_t  MuonBeamGen::Init()
{
    // Initialize generator
    
    
    if (fPMean<0.01)
        Fatal("Init()","MuonBeamGen: Mean momentum must be positive and meaningfull!");
    if (fR<0.00000001) {
        Fatal("Init()","MuonBeamGen: Radius must be positive!");
    }
    if (fThetaMean<0) {
        Fatal("Init()","MuonBeamGen: ThetaMean must be positive or null");
    }
    if (fEleFrac<0|| fEleFrac>1) {
        Fatal("Init()","MuonBeamGen: Electron fraction should be between 0 and 1");
    }
    
     // Check for particle type
    TDatabasePDG* pdgBase = TDatabasePDG::Instance();
    TParticlePDG* particle = pdgBase->GetParticle(fPDGType);
    
    if (!particle) {
        Fatal("MuonBeamGen","PDG code %d not defined.",fPDGType);
    } else {
        fPDGMass = particle->Mass();
        fPDGLifetime = particle->Lifetime();
        
        std::cout << " MuonBeamGen::Init() ------------  " << std::endl;
        printf("particle->Mass() = %f \n", fPDGMass);
        printf("particle->Lifetime() = %f \n", fPDGLifetime);
        fPDGLifetime=2197.03;
        return kTRUE;
    }
    
    return kTRUE;
}

// ------------------------------------------------------------------------
Bool_t MuonBeamGen::ReadEvent(FairPrimaryGenerator* primGen)
{
    // Generate one event: produce primary particles emitted from one vertex.
    // Primary particles are distributed uniformly along
    // those kinematics variables which were limitted by setters.
    // if SetCosTheta() function is used, the distribution will be uniform in
    // cos(theta)
    
    
    Double32_t mom=0, theta=0, phi=0, pt=0, px, py, pz;
    Double32_t u,v,ruvq, ruv;
    do {
        u = gRandom->Uniform(-1., 1.);
        v = gRandom->Uniform(-1., 1.);
        ruvq = u*u+v*v;
    } while (ruvq>1.);
    ruv = sqrt(ruvq); // now u/ruv and v/ruv are cos and sin of phi
    
    mom = gRandom->Exp(fPMean-0.01)+0.01;
    //  mom = gRandom->Uniform(fPMean-0.01,fPMean+0.01);
    theta = gRandom->Exp(fThetaMean);
    pz = mom*cos(theta);
    pt = mom*sin(theta);
    px = pt*u/ruv;
    py = pt*v/ruv;

    static double const	c_0 = 299792458;
    Double_t gamma = sqrt(1.+pow(mom/(fPDGMass*c_0),2));
    
    // choose particle type
    int charge, ptype;
    charge = -1; // assume negative particles
    if( gRandom->Uniform(-1., 1.)>0 ) charge = 1;
    ptype = 13; // assume negative muon
    if( fEleFrac>0 && gRandom->Uniform(0, 1.)<fEleFrac ) ptype = 11;
    
    if( charge>0 ) ptype = -ptype; // ptype>0 means charge<0!!
    
    fPDGType = ptype;
    
    // find a point on a disk of radius fR
    Double32_t r, rq, x, y, z, b;
    
    rq=2*fR*fR;
    do {
        x = gRandom->Uniform(-fR, fR);
        y = gRandom->Uniform(-fR, fR);
        rq = x*x+y*y;
    } while (rq>fR*fR);
    r = sqrt(rq);
    // find a decay Z-position
//    do {
//        z = gRandom->Uniform(fZ,-fZ);
//        b = gRandom->Uniform(0.,1.);
//    } while (b > exp((z-fZ)*fPDGMass/(pz*gamma*fPDGLifetime)));
    
    do {
        b = gRandom->Uniform(0.,1.);
        // std::cout << " mean Lifetime = " << pz*gamma*fPDGLifetime/fPDGMass<< std::endl;
        // std::cout << "fZ = " << fZ << "  z = " << fZ-log(b)*(pz*gamma*fPDGLifetime)/fPDGMass<< std::endl;
        z = fZ-log(b)*(pz*gamma*fPDGLifetime)/fPDGMass;
    } while (z>-fZ);
    //std::cout << " z = " << z << "  [fZ,-fZ] = [" << fZ << ","<< -fZ << "]"<<std::endl;
    // z = fZ;   // now x, y, z contains the primary vertex.
             
    // Generate particle
    if (fDebug)
        printf("BoxGen: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n",
               ptype, px, py, pz, x, y, z);
  
    TGenPhaseSpace event;
    TLorentzVector W(px,py,pz,sqrt(px*px+py*py+pz*pz+fPDGMass*fPDGMass));
    Int_t    elPDG = 11, nu_elPDG = -12, nu_muPDG=14;
    if (ptype<0) {
        elPDG = -elPDG;
        nu_elPDG = -nu_elPDG;
        nu_muPDG = -nu_muPDG;
        }
    
    if(fDebug) std::cout<< " elPDG=" <<elPDG << " \t nu_elPDG = " << nu_elPDG << "\t nu_muPDG = " << nu_muPDG << std::endl;

    TDatabasePDG* pdgBase = TDatabasePDG::Instance();
    TParticlePDG* electron = pdgBase->GetParticle(elPDG);
    TParticlePDG* nu_el = pdgBase->GetParticle(nu_elPDG);
    TParticlePDG* nu_mu = pdgBase->GetParticle(nu_muPDG);

    Double_t masses[3] ={electron->Mass(), nu_el->Mass(), nu_mu->Mass()};
    event.SetDecay(W,3,masses);
    event.Generate();
 
    if (fDebug){
        printf("muon mass: %f \n",fPDGMass);
        for (int i=0; i<3; i++)  printf("masses: %f \n", masses[i]);
        std::cout << "Muon TLorentzVector W = ";
        W.Print();
        std::cout << "mass =" << W.M() <<" P= " << W.Px() << W. Py()<< W.Pz() << std::endl;
        std::cout << "Electron TLorentzVector =";   event.GetDecay(0)->Print();
        std::cout << event.GetDecay(0)->M() << std::endl;
        printf("SetDecay: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n",
               elPDG,   event.GetDecay(0)->Px(),event.GetDecay(0)->Py(),event.GetDecay(0)->Pz(),x,y,z);
        printf("SetDecay: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n",nu_elPDG,event.GetDecay(1)->Px(),event.GetDecay(1)->Py(),event.GetDecay(1)->Pz(),x,y,z);
        printf("SetDecay: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n",nu_muPDG,event.GetDecay(2)->Px(),event.GetDecay(2)->Py(),event.GetDecay(2)->Pz(),x,y,z);
        printf("SetDecay: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n",ptype,  -px, -py, -pz, x, y, z);
    }
    
    primGen->AddTrack(elPDG,   event.GetDecay(0)->Px(),event.GetDecay(0)->Py(),event.GetDecay(0)->Pz(),x,y,z);
    primGen->AddTrack(nu_elPDG,event.GetDecay(1)->Px(),event.GetDecay(1)->Py(),event.GetDecay(1)->Pz(),x,y,z);
    primGen->AddTrack(nu_muPDG,event.GetDecay(2)->Px(),event.GetDecay(2)->Py(),event.GetDecay(2)->Pz(),x,y,z);
    primGen->AddTrack(ptype,  -px, -py, -pz, x, y, z);
    return kTRUE;
             }
// ------------------------------------------------------------------------


ClassImp(MuonBeamGen)
