//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnCluster.hh,v 1.1 2014/06/13 06:14:48 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:14:48 $
//
// Contact person Pavel Murat
//
#ifndef murat_inc_TStnCluster_hh
#define murat_inc_TStnCluster_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

//namespace murat {

class TStnTrack;

namespace mu2e {
  class CaloCluster;
}

class TStnCluster : public TObject {

  enum {
    kNFreeIntsV1   = 10,		// V1
    kNFreeFloatsV1 = 10,		// V1

    kNFreeInts     = 10,		// V2
    kNFreeFloats   =  3			// V2
  };

public:
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
  int                       fNumber;          // index in the list of reconstructed clusters
  int                       fDiskID;	      // 
  int                       fNCrystals;       //
  int                       fNCr1     ;       // N crystals above 1 MeV
  int                       fTrackNumber;     // closest track in TStnTrackBlock
  int                       fIx1;	      // [row, column] or [x1,x2] for a disk
  int                       fIx2;
  int                       fInt[kNFreeInts];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
  float                     fX;
  float                     fY;
  float                     fZ;
  float                     fYMean;
  float                     fZMean;
  float                     fSigY;      // cluster width in Y 
  float                     fSigZ;	// cluster width in Z/X (fSigZ is reused)
  float                     fSigR;      // don't know what it is
  float                     fEnergy   ;
  float                     fTime     ; 
  float                     fFrE1     ; // e1/etotal
  float                     fFrE2     ; // (e1+e2)/etotal
  float                     fSigE1    ;
  float                     fSigE2    ;
//-----------------------------------------------------------------------------
// 7 words added in version 2
//-----------------------------------------------------------------------------
  float                     fSigXX;    // sums over crystals
  float                     fSigXY;
  float                     fSigYY;
  float                     fNx   ;    // cluster direction
  float                     fNy   ;
  float                     fFloat[kNFreeFloats];
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
  const mu2e::CaloCluster*  fCaloCluster;  //!
  TStnTrack*                fClosestTrack; //!
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnCluster(int i = -1);
  ~TStnCluster();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     DiskID     () { return fDiskID  ; }
  int     Ix1        () { return fIx1; }
  int     Ix2        () { return fIx2; }
  int     NCrystals  () { return fNCrystals; }
  int     Number     () { return fNumber; }
  int     TrackNumber() { return fTrackNumber; }

  float   Energy     () { return fEnergy; }
  float   Time       () { return fTime  ; }

  float   Nx         () { return fNx;    }
  float   Ny         () { return fNy;    }

  float   SigX       () { return fSigZ;  }
  float   SigY       () { return fSigY;  }
  float   SigZ       () { return fSigZ;  }

  float   SigXX      () { return fSigXX; }
  float   SigXY      () { return fSigXY; }
  float   SigYY      () { return fSigYY; }

  TStnTrack* ClosestTrack() { return fClosestTrack; }

  void    SetNumber(int I) { fNumber = I; } // 

//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* Opt = "") ;
  virtual void Print(Option_t* Opt = "") const ;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer& R__b);

  ClassDef(TStnCluster,2)
};

#endif
