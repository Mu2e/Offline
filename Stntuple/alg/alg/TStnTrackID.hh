#ifndef murat_TStnTrackID
#define murat_TStnTrackID
///////////////////////////////////////////////////////////////////////////////
// class representing Track ID cuts
// Author:    P. Murat, borrowed from CDF code
// Date:      
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TString.h"

class  TBuffer;

class TStnTrack;
class TH1F;

class TStnTrackID: public TNamed {
public:
//-----------------------------------------------------------------------------
//  try 
//-----------------------------------------------------------------------------
  enum { 
    kNActiveBit        = 0x1 <<  0,     // 0x00000001
    kFitConsBit        = 0x1 <<  1,     // 0x00000002
    kT0Bit             = 0x1 <<  2,     // 0x00000004
    kT0ErrBit          = 0x1 <<  3,     // 0x00000008
    kFitMomErrBit      = 0x1 <<  4,     // 0x00000010
    kTanDipBit         = 0x1 <<  5,     // 0x00000020
    kD1Bit             = 0x1 <<  6,     // 0x00000040
    kD2Bit             = 0x1 <<  7      // 0x00000080
  };

  enum { 
    kNFreeInts         =  4,
    kNFreeFloats       = 10
  };

  struct Hist_t {
					// tight object cuts
    TH1F*    fNActive     [5];
    TH1F*    fFitCons     [5];
    TH1F*    fT0          [5];
    TH1F*    fT0Err       [5];
    TH1F*    fFitMomErr   [5];
    TH1F*    fTanDip      [5];
    TH1F*    fD1          [5];
    TH1F*    fD2          [5];
					// summary histogram
    TH1F*    fFailedBits;
    TH1F*    fPassed;
  };

protected:
  Int_t      fUseMask;
  Int_t      fMinNActive;
  Int_t      fMaxNActive;
  Int_t      fInteger[kNFreeInts];		// for future

  Float_t    fMinFitCons;		// 
  Float_t    fMinT0;
  Float_t    fMaxT0Err;

  Float_t    fMaxFitMomErr;
  Float_t    fMinTanDip;
  Float_t    fMaxTanDip;

  Float_t    fMinD1;			// R1 - minimal radius, asymmetric cut
  Float_t    fMaxD1;

  Float_t    fMinD2;			// R2 - maximal track radius, asymmetric cut
  Float_t    fMaxD2;

  Float_t    fFloat[kNFreeFloats];	// spare words, added in V5

  void*      fEOR;		// ! end of record
//-----------------------------------------------------------------------------
//  methods
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor

  TStnTrackID(const char* Name = "Anonymous");
  virtual ~TStnTrackID();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Float_t MinFitCons   () const { return fMinFitCons;   }
  Float_t MinT0        () const { return fMinT0;        }
  Int_t   MinNActive   () const { return fMinNActive;   }
  Int_t   MaxNActive   () const { return fMaxNActive;   }
  Float_t MaxT0Err     () const { return fMaxT0Err;     }
  Float_t MaxFitMomErr () const { return fMaxFitMomErr; }
  Float_t MinTanDip    () const { return fMinTanDip;    }
  Float_t MaxTanDip    () const { return fMaxTanDip;    }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    SetMinFitCons  (Float_t FitCons) { fMinFitCons    = FitCons; }
  void    SetMinT0       (Float_t T0     ) { fMinT0         = T0;      }
  void    SetMinNActive  (Int_t   N      ) { fMinNActive    = N;       }
  void    SetMaxNActive  (Int_t   N      ) { fMaxNActive    = N;       }
  void    SetMaxT0Err    (Float_t T0Err  ) { fMaxT0Err      = T0Err;   }
  void    SetMaxFitMomErr(Float_t MomErr ) { fMaxFitMomErr  = MomErr;  }

  void    SetMinTanDip   (Float_t TanDip ) { fMinTanDip     = TanDip;  }
  void    SetMaxTanDip   (Float_t TanDip ) { fMaxTanDip     = TanDip;  }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  virtual Int_t  IDWord     (TStnTrack* Track);
  virtual Int_t  TightIDWord(TStnTrack* Track);
  virtual Int_t  LooseIDWord(TStnTrack* Track);

                                        // Mode=1: tight ID, =2: loose ID

  void FillHistograms(Hist_t& Hist, TStnTrack* Track, Int_t Mode=1);
//-----------------------------------------------------------------------------
//  overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Print  (Option_t*  Option = "") const ;
//-----------------------------------------------------------------------------
//  schema evolution
//  2009-12-24: current version - V5
//-----------------------------------------------------------------------------
//   void  ReadV1(TBuffer& R__b);
//   void  ReadV2(TBuffer& R__b);
//   void  ReadV3(TBuffer& R__b);
//   void  ReadV4(TBuffer& R__b);

  ClassDef(TStnTrackID,0)

};

#endif
