#ifndef ECALANALYSISCALIB_H
#define ECALANALYSISCALIB_H

#include "FairTask.h"

#include "TString.h"

#include <list>

class TTree;
class ecalStructure;
class TClonesArray;

class ecalAnalysisSimple : public FairTask
{
public:
  ecalAnalysisSimple(const char* name, const Int_t iVerbose);

  /** Default constructor **/
  ecalAnalysisSimple();

  /** Initing routine **/
  virtual InitStatus Init();

  /** Loop procedure **/
  virtual void Exec(Option_t* option);

  /** Finishing routine **/
  virtual void Finish();

  /** Destructor **/
  virtual ~ecalAnalysisSimple() {};

private:
  TTree* fTree;
  TTree* eTree;
  Double_t fX;
  Double_t fY;
  Double_t fCX;
  Double_t fCY;
  Double_t fP;
  Double_t fCE;
  Double_t fOE;
  Double_t fPX;
  Double_t fPY;
  Double_t fPZ;
  Int_t fEv;
  Int_t fCellNum;
  Int_t fADC;
  Int_t fspdg;
  Double_t fsx,fsy,fsz,fspx,fspy,fspz;
  Int_t fdpdg;
  Double_t fdx,fdy,fdz,fdpx,fdpy,fdpz,fde,fdtime;
  Bool_t    fdecay;
  Bool_t    isDecay;
  Int_t fmotherid;
  Int_t eDecTrackId;
  Double_t fElEloss;
    
  void InitTree();
  ecalStructure* fStr;
  TClonesArray* fTracks;
  TClonesArray* fPLTracks;
  TClonesArray* fMCTracks;

  ecalAnalysisSimple(const ecalAnalysisSimple&);
  ecalAnalysisSimple operator=(const ecalAnalysisSimple&);

  ClassDef(ecalAnalysisSimple,1)
};

#endif 

