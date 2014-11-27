/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *  
 *                  copied verbatim in the file "LICENSE"                       *
 ******************************************************************************* */

#ifndef MUONBEAMGEN_H
#define MUONBEAMGEN_H

#include "FairGenerator.h"              // for FairGenerator

#include "Rtypes.h"                     // for Double32_t, Bool_t, kTRUE, etc

class FairPrimaryGenerator;

class MuonBeamGen : public FairGenerator
{
  public:

    /** Default constructor. **/
    MuonBeamGen();

    /** Constructor with PDG-ID, multiplicity
     **@param pdgid Particle type (PDG encoding)
     **@param mult  Multiplicity (default is 1)
     **/
    MuonBeamGen(Int_t pdgid, Int_t mult=1, Int_t debug=0);

    /** Destructor **/
    virtual ~MuonBeamGen() {};

    /** Modifiers **/
    //    void SetPDGType      (Int_t pdg)  {fPDGType = pdg;  };

    //    void SetMultiplicity (Int_t mult) {fMult    = mult; };
 
    void SetPMean(Double32_t pmean=0) {fPMean=pmean;}

    void SetMeanTheta   (Double32_t thetamean=0) {fThetaMean=thetamean; };

    void SetRZ   (Double32_t r=0, Double32_t z=0) {
      fR=r;
      fZ=z;
    }

    void SetEleFrac   (Double32_t eleFrac=0) {fEleFrac=eleFrac; };

    void SetDebug(Bool_t debug=0) {fDebug = debug;}

       
    /** Initializer **/
    Bool_t Init();

    /** Creates an event with given type and multiplicity.
     **@param primGen  pointer to the FairPrimaryGenerator
     **/
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

  private:
    Int_t      fPDGType;             // Particle type (PDG encoding)
    // Int_t      fMult;                // Multiplicity
    Double32_t fPDGMass;             // Particle mass [GeV]
    Double32_t fPDGLifetime;         // Particle lifetime [s]
    
    Double32_t fPMean;               // Mean momentum range in lab system
    Double32_t fThetaMean;           // Mean Polar angle range in lab system
                                     //   radiant
    Double32_t fR, fZ;               // Point vertex coordinates [cm]
                                     // disk at Z or radius R
    Double32_t fEleFrac;             // Electron/positron fractions
    Bool_t     fDebug;               // Debug switch

    ClassDef(MuonBeamGen,3);

    
};


#endif
