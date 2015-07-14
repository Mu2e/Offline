//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnPid.hh,v 1.1 2014/06/13 06:14:48 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:14:48 $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_obj_TStnPid_hh
#define Stntuple_obj_TStnPid_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

class TStnTrack;

namespace mu2e {
  class AvikPIDProduct;
}

class TStnPid : public TObject {

  enum {
    kNFreeInts     = 10,		// V1
    kNFreeFloats   = 10			// V1
  };

public:
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
  int                       fTrackNumber;     // track number
  int                       fNMatched;
  int                       fNMatchedAll;
  int                       fNUsedOsEle;      // for fSumAvikOsEle
  int                       fNUsedOsMuo;
  int                       fNUsedSsEle;      // for fDrdsSsEle
  int                       fNUsedSsMuo;
  int                       fInt[kNFreeInts];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
  float                     fLogDedxProbEle;
  float                     fLogDedxProbMuo;

  float                     fDrdsVadimEle;
  float                     fDrdsVadimEleErr;
  float                     fDrdsVadimMuo;
  float                     fDrdsVadimMuoErr;
  
  float                     fSumAvikEle;
  float                     fSumAvikMuo;

  float                     fSq2AvikEle;
  float                     fSq2AvikMuo;

  float                     fDrdsOsEle;
  float                     fDrdsOsEleErr;
  float                     fDrdsOsMuo;
  float                     fDrdsOsMuoErr;
  
  float                     fSumAvikOsEle;
  float                     fSumAvikOsMuo;

  float                     fDrdsSsEle;
  float                     fDrdsSsEleErr;
  float                     fDrdsSsMuo;
  float                     fDrdsSsMuoErr;
  float                     fFloat[kNFreeFloats];
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
  const mu2e::AvikPIDProduct*  fAvikPid;  //!
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnPid(int i = -1);
  ~TStnPid();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     TrackNumber() { return fTrackNumber; }

  float   DrdsVadimEle   () { return fDrdsVadimEle;    }
  float   DrdsVadimEleErr() { return fDrdsVadimEleErr; }
  float   DrdsVadimMuo   () { return fDrdsVadimMuo;    }
  float   DrdsVadimEleMuo() { return fDrdsVadimMuoErr; }

					// the rest accessors - to come
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* Opt = "") ;
  virtual void Print(Option_t* Opt = "") const ;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
//  void ReadV1(TBuffer& R__b);

  ClassDef(TStnPid,1)
};

#endif
