//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TDisk.hh,v 1.1 2014/06/13 05:30:42 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 05:30:42 $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_obj_TDisk_hh
#define Stntuple_obj_TDisk_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>
#include <math.h>

#include "TMath.h"
#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TVector2.h"

#include "Stntuple/base/THexIndex.hh"

class TStnCrystal;

namespace mu2e {
  class Disk;
}

class TDisk : public TObject {
public:
  int                fSectionID;
					
  TObjArray*         fListOfCrystals;	// list of crystals
  int                fNCrystals;

  int                fNRings;
  //  int                fNInsideTot;
  int                fChannelOffset;             // offset of the first channel

  int*               fFirst   ;	                // index of the first crystal in the ring
  int*               fNCrystalsPerRing ;          // number of crystals in the i-th ring
  int*               fNInside ;
  
  double             fRMin;
  double             fRMax;
  double             fZ0;
  double             fHexSize;
  double             fDeadSpace;	// wrapper+shell
  double             fMinFraction;	// min fraction of area to add a crystal
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TDisk();
  TDisk(int      SectionID  , 
	double   RMin       , 
	double   RMax       , 
	double   Z0         ,
	double   HexSize    , 
	double   DeadSpace  , 
	double   MinFraction);

  ~TDisk();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int            SectionID() { return fSectionID; }
  TObjArray*     ListOfCrystals() { return fListOfCrystals; }
  int            NCrystals()      { return fListOfCrystals->GetEntriesFast(); }
  int            First    (int I) { return fFirst[I]; }
  TStnCrystal*   Crystal  (int I) { return (TStnCrystal*) fListOfCrystals->UncheckedAt(I); }


  void           GetHexIndex(int I, THexIndex* HexIndex);
				// for crystal number I return its ring number
  int            GetRing (int I);
  int            GetRing (THexIndex* Index);

  double         GetRMin() { return fRMin; }
  double         GetRMax() { return fRMax; }

  int            GetNRings   () { return fNRings;    }
  int            GetNCrystals() { return fNCrystals; }

  void           GetPosition(int I          , TVector2* Pos);
  void           GetPosition(THexIndex* Index, TVector2* Pos);

  int            GetNCrystalsPerRing(int Ring) { return fNCrystalsPerRing[Ring]; }
  int            GetNInside         (int Ring) { return fNInside         [Ring]; }

  //  int       GetNInsideTot() { return fNInsideTot;    }

  double         GetActiveArea() { return GetCrystalArea()*fNCrystals; }
  double         GetTotalArea () { return (fRMax*fRMax-fRMin*fRMin)*TMath::Pi() ; }

  double         GetRadius(int I);
  double         GetRadius(THexIndex* Index);

  double         GetCrystalArea() { return  fHexSize*fHexSize*sqrt(3.)/2.; }

  // calculates inside fraction  
  int            IsInside(THexIndex* Index, double *Fraction);

  int            ChannelOffset() { return fChannelOffset; }

//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void           SetChannelOffset(int Offset) { fChannelOffset = Offset; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void   Paint(Option_t* Opt = "");
  virtual void   Clear(Option_t* Opt = "") ;
  virtual void   Print(Option_t* Opt = "") const ;

  ClassDef(TDisk,0)
    
};

#endif
