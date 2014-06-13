#ifndef STNTUPLE_TStnHeaderBlock
#define STNTUPLE_TStnHeaderBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE event header
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Oct 31 2000
// 
//-----------------------------------------------------------------------------
#include "TStnDataBlock.hh"

class TStnEvent;

class TStnHeaderBlock : public TStnDataBlock {

  friend Int_t StntupleInitMu2eHeaderBlock     (TStnDataBlock*, AbsEvent* , int);
  friend Int_t StntupleInitMu2eHeaderBlockLinks(TStnDataBlock*, AbsEvent* , int);

public:
  Int_t             fVersion;
  Int_t             fEventNumber;
  Int_t             fRunNumber;
  Int_t             fSectionNumber;	// section number within the run
  Int_t             fMcFlag;		// MC flag, 0 for real data
  Int_t             fGoodRun;		// run flag
  Int_t             fBrCode;		// 
  Int_t             fGoodTrig;
  Int_t             fTrigWord;		// 
  Int_t             fNTracks;           //
  Int_t             fCpu;               // packed word with processing time
  Float_t           fInstLum;		// instantaneous luminosity
  TString           fStnVersion;        // like dev_243_16

//------------------------------------------------------------------------------
//  now the data which are not a part of STNTUPLE
//------------------------------------------------------------------------------
					// number/run number for the last 
					// printed event

  Int_t             fLastNumber;	// !
  Int_t             fLastRunNumber;	// ! 
//------------------------------------------------------------------------------
//  function members
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnHeaderBlock();
  virtual ~TStnHeaderBlock();
					// ****** accessors 

  Int_t  EventNumber  () const { return fEventNumber;   }
  Int_t  RunNumber    () const { return fRunNumber;     }
  Int_t  SectionNumber() const { return fSectionNumber; }
  Int_t  McFlag       () const { return fMcFlag;        }
  Int_t  NTracks      () const { return fNTracks;       }

  Float_t InstLum     () const { return fInstLum;       }

  Float_t CpuTime     () const { return ((fCpu>>8)/10.0);   } // in s
  Float_t CpuSpeed    () const { return ((fCpu&0xFF)/5.0); } // in GHz
  const TString& StnVersion () const { return fStnVersion;    }

					// ****** setters/modifiers

					// ****** overloaded functions of 
					// TObject
  void   Clear(Option_t* opt = "");
  void   Print(Option_t* opt = "") const;

					// ****** schema evolution
  void   ReadV50(TBuffer& R__b);


  ClassDef(TStnHeaderBlock,1)	       // Mu2e STNTUPLE event header
};

#endif
