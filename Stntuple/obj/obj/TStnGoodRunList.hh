#ifndef STNTUPLE_TStnGoodRunList
#define STNTUPLE_TStnGoodRunList
//-----------------------------------------------------------------------------
//  definition of a good run list for STNTUPLE
//  Author:    P.Murat (CDF/FNAL)
//  Date:      Jun 14 2002
//-----------------------------------------------------------------------------
#include "TNamed.h"
#include "TArrayI.h"

class  TFile;
class  TTree;
class  TStnRunSummary;

class TStnGoodRunList : public TNamed {
//------------------------------------------------------------------------------
// data members / supposed to have name "GoodRunList"
// natural to assume that there is one and only one per job...
//------------------------------------------------------------------------------
protected:
  Float_t           fMinGoodRunLum;	// min lumi for a good run
  TArrayI           fListOfRuns;	// list of good run ranges
  Int_t             fNRanges;		// 
  Int_t             fUseUncheckedRuns;  //
  TFile*            fFile;		// file with good run ntuple
  TTree*            fTree;		// good run ntuple
  TStnRunSummary*   fRunSummary;	// !
  Double_t          fCurrentEntry;	// !
  Double_t          fNEntries;		// !
  Int_t             fMinRunNumber;	// !
  Int_t             fMaxRunNumber;	// !
  Int_t             fElectronFlag;      // ! various requirements
  Int_t             fMuonFlag;          // !
  Int_t             fSiliconFlag;       // !

  Int_t            (*fGoodRunRoutine)(int RunNumber, int RunSection, int Mask); // !
//------------------------------------------------------------------------------
//  function members
//------------------------------------------------------------------------------
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TStnGoodRunList(const char* Name     = "GoodRunList",
		  const char* Filename = NULL);

  virtual ~TStnGoodRunList();
//-----------------------------------------------------------------------------
// accessors 
//-----------------------------------------------------------------------------
  Float_t  MinGoodRunLum   () const { return fMinGoodRunLum; }
  TArrayI* ListOfRuns      ()       { return &fListOfRuns; }
  Int_t    NRanges         () const { return fNRanges    ; }
  Int_t    UseUncheckedRuns() const { return fUseUncheckedRuns; }

  Int_t    MinRunNumber(Int_t Range) { return fListOfRuns.At(2*Range  ); }
  Int_t    MaxRunNumber(Int_t Range) { return fListOfRuns.At(2*Range+1); }

  TStnRunSummary* GetRunSummary(Int_t RunNumber);
  TStnRunSummary* GetEntry     (Int_t RunNumber);

  Double_t   NEntries() const { return fNEntries; }
//-----------------------------------------------------------------------------
// setters/modifiers
//-----------------------------------------------------------------------------
  void     SetMinGoodRunLum   (Float_t  Lum) { fMinGoodRunLum    = Lum; }
  void     SetUseUncheckedRuns(Int_t    Use) { fUseUncheckedRuns = Use; }
  void     SetListOfRuns      (Int_t*  List);
  Int_t    Init               (const char* Filename = 0);
//-----------------------------------------------------------------------------
//  good run list routines - those are static
//-----------------------------------------------------------------------------
				// call with RunNumber<0 sets print level

  void  SetGoodRunRoutine(Int_t (*F)(Int_t,Int_t,Int_t)) { 
    fGoodRunRoutine = F;
    if (fGoodRunRoutine) fGoodRunRoutine(-1,-1,0);
  }
				// by default good run checking is not set

  Int_t         GoodRun    (Int_t RunNumber, Int_t RunSection = -1, Int_t Mask = 0) { 
    return (! fGoodRunRoutine) ? 1 : fGoodRunRoutine(RunNumber,RunSection,Mask);
  }

  Int_t  ElectronFlag() { return fElectronFlag; }
  Int_t  MuonFlag    () { return fMuonFlag;     }
  Int_t  SiliconFlag () { return fSiliconFlag;  }
//-----------------------------------------------------------------------------
// various concrete implementations of GRL
//-----------------------------------------------------------------------------
protected:
  static int GoodRunListEtf       (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_WZ_PRD   (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_DQM_V6   (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_DQM_V7   (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_DQM_V13  (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_DQM_V27  (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_DQM_V32  (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_DQM_V34  (int RunNumber, int RunSection, int Mask);
  static int GoodRunList_MC_1001  (int RunNumber, int RunSection, int Mask);
  static int DefaultGoodRunRoutine(int RunNumber, int RunSection, int Mask);
public:
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void   Clear(Option_t* Opt = "");
  virtual void   Print(Option_t* Opt = "") const;

  ClassDef(TStnGoodRunList,1)

};

#endif
