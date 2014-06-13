//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnCluster.hh,v 1.1 2014/06/13 05:30:42 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 05:30:42 $
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
    kNFreeInts   = 10,
    kNFreeFloats = 10
  };

public:
  int                       fNumber;          // track index in the list of reconstructed clusters
  int                       fDiskID;	      // 
  int                       fNCrystals;       //
  int                       fNCr1     ;       // above 1 MeV
  int                       fTrackNumber;     // closest track in TStnTrackBlock
  int                       fIx1;	      // [row, column] or [x1,x2] for a disk
  int                       fIx2;
  int                       fInt[kNFreeInts];
					      // float part
  float                     fX;
  float                     fY;
  float                     fZ;
  float                     fYMean;
  float                     fZMean;
  float                     fSigY;
  float                     fSigZ;
  float                     fSigR;
  float                     fEnergy;
  float                     fTime     ; 
  float                     fFrE1     ; // e1/etotal
  float                     fFrE2     ; // (e1+e2)/etotal
  float                     fSigE1    ;
  float                     fSigE2    ;
  float                     fFloat[kNFreeFloats];
//-----------------------------------------------------------------------------
// transient variables
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

  TStnTrack* ClosestTrack() { return fClosestTrack; }

  void    SetNumber(int I) { fNumber = I; } // 

  void    Print(Option_t* Opt) const ;
					// I/O: Mu2e version 1
  ClassDef(TStnCluster,1)

};

#endif
