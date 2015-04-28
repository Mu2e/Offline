//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnTrack.hh,v 1.3 2014/09/13 01:16:28 murat Exp $
// $Author: murat $
// $Date: 2014/09/13 01:16:28 $
//
// Contact person Pavel Murat
//
#ifndef murat_inc_TStnTrack_hh
#define murat_inc_TStnTrack_hh

// Mu2e includes.

// #include "CLHEP/Geometry/HepPoint.h"
// #include "CLHEP/Vector/ThreeVector.h"
// #include "CLHEP/Matrix/SymMatrix.h"
// #include "CLHEP/Matrix/Vector.h"
// 
// #include "TrkBase/TrkHelixUtils.hh"
// #include "TrkBase/HelixParams.hh"
// #include "KalmanTrack/KalHit.hh"
// #include "KalmanTests/inc/KalRepCollection.hh"

// storable objects (data products)
// #include "RecoDataProducts/inc/StrawHitCollection.hh"
// #include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
// #include "RecoDataProducts/inc/CaloHitCollection.hh"
// #include "RecoDataProducts/inc/CaloClusterCollection.hh"
// #include "MCDataProducts/inc/GenParticleCollection.hh"


// #include "TrackCaloMatching/inc/TrackClusterLink.hh"
// #include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
// #include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Stntuple/base/TBitset.hh"
// #include "TStnCluster.hh"

					// 'KalRep' is a BaBar class
class KalRep ;
class TStnCluster;

namespace mu2e {
  class CaloCluster;
  class TrkToCaloExtrapol;
}

// namespace murat {

class TStnTrack : public TObject {

  enum {
    kNFreeIntsV1   =  9,
    kNFreeFloatsV1 = 10,

    kNFreeIntsV2   =  9,
    kNFreeFloatsV2 =  6,

    kNFreeIntsV3   =  9,
    kNFreeFloatsV3 =  6,

    kNFreeIntsV4   =  7,		// hopefully, correct
    kNFreeFloatsV4 =  6,		// hopefully, correct

    kNFreeIntsV5   =  7,
    kNFreeFloatsV5 =  6,

    kNFreeIntsV6   =  7,
    kNFreeFloatsV6 =  6,

    kNFreeIntsV7   =  7,
    kNFreeFloatsV7 =  5,

    kNFreeInts     =  7,
    kNFreeFloats   =  4
  };

  enum { kMaxNLayers = 88 }; // 22x2*2

					// added 3 chi*2's
public:
  enum { kNVanesV4 = 4,
	 kNVanesV5 = 4,
	 kNDisks   = 2			// there are only two disks, no vane support
  };

  struct InterData_t {
    int          fID;			// = -1 if no intersection
    int          fClusterIndex;         // cluster index in the list of clusters
    float        fTime;			// extrapolated track time, not corrected by _dtOffset
    float        fEnergy;		// closest cluster energy
    float        fXTrk;
    float        fYTrk;
    float        fZTrk;
    float        fNxTrk;		// track direction cosines in the intersection point
    float        fNyTrk;
    float        fNzTrk;
    float        fXCl;			// cluster coordinates
    float        fYCl;
    float        fZCl;
    float        fDx;			// TRK-CL
    float        fDy;			// TRK-CL
    float        fDz;
    float        fDt;			// TRK-CL , _corrected_ by _dTOffset (!)
    float        fDu;			// added in V6
    float        fDv;			// added in V6
    float        fChi2Match;		// track-cluster match chi&^2 (coord)
    float        fChi2Time;		// track-cluster match chi&^2 (time)
    float        fPath;			// track path in the disk
    float        fIntDepth;             // assumed interaction depth, added in V6;
    const mu2e::CaloCluster*       fCluster;
    const mu2e::TrkToCaloExtrapol* fExtrk;
  };
    
  TLorentzVector            fMomentum;         // this assumes DELE fit hypothesis
  
  TBitset                   fHitMask;	       // bit #i: 1 if there is a hit 
  TBitset                   fExpectedHitMask;   // bit #i: 1 if expect to have a hit at this Z

  int                       fNumber;       // track index in the list of reconstructed tracks
  int                       fNHyp;          // number of hyp's with successfull fits
  int                       fBestHyp[2];    // hypothesis with the best chi2/ndof
  int                       fIDWord;	    // now - for selection "C"
  
  int                       fNActive;	    // NWrong << 16 + NActive
  int                       fVaneID;	    // 
  int                       fDiskID;	    // 
  int                       fPdgCode;       // PDF code of the particle produced most hits
  int                       fNGoodMcHits;   // Nhits produced by the associated MC particle
  int                       fPartID;        // MC particle ID (number in the list)
  int                       fNMcStrawHits;  // Nhits by associated particle in the straw tracker
  int                       fAlgorithmID;   // bit-packed : (alg_mask << 16 ) | best
  int                       fInt[kNFreeInts];     // provision for future expension
  
  float                     fChi2;
  float                     fChi2C;          // *NOT USED* 
  float                     fFitCons;
  float                     fT0;
  float                     fT0Err;
  float                     fFitMomErr;
  float                     fTanDip;	 // at Z=Z0
  float                     fP;		 // total momentum in the first point
  float                     fCharge;
  float                     fPt;	 // transverse momentum in the first point
  float                     fD0;	 // at Z=Z0
  float                     fZ0;

  float                     fPStOut;    // MC momentum in the VD ST_Out 
  float                     fPFront;    // MC momentum in the VD front of the tracker

  float                     fClusterE;	// energy of the associated cluster
  float                     fDt;
  float                     fEp;
  float                     fDx;	// about 0 for vanes
  float                     fDy;	// 
  float                     fDz;	// about 0 for disks
  
  float                     fEleLogLHCal;         // log likelihood of the electron hypothesis
  float                     fMuoLogLHCal;	  // log lilelihood of the muon     hypothesis

  float                     fRSlope;	          // timing residual slope dres(T)/dZ
  float                     fRSlopeErr;

  float                     fLogLHRXs;            // XSlope-only-based likelihood

  float                     fEleLogLHDeDx;        // dE/dX LH calculated by Vadim based 
  float                     fMuoLogLHDeDx;        // 
  float                     fX1;	          // momentum defined at Z1
  float                     fY1;
  float                     fZ1;
  float                     fP0;                  // momentum defined at Z0
  float                     fP2;                  // momentum defined at Z0
  float                     fC0;	          // curvature at Z0
  float                     fPhi0;	          // phi0 at Z0 **specify in V8 , no I/O changes***
  float                     fFloat[kNFreeFloats]; // provision for future I/O expansion

  InterData_t               fDisk[kNDisks];       // track intersections with disks
//-----------------------------------------------------------------------------
//  transient data members, all the persistent ones should go above
//-----------------------------------------------------------------------------
  const mu2e::TrkToCaloExtrapol* fExtrk;              //! 
  const mu2e::CaloCluster*       fClosestCaloCluster; //!

  InterData_t*                   fVMinS;	      //! intersection with min S
  InterData_t*                   fVMaxEp;	      //! intersection with max E/P
  KalRep*                        fKalRep[4];          //! different fits, sequence: e-, e+, mu-, mu+

  int                            fITmp[ 5];           //!
  Float_t                        fTmp [10];           //! for temporary analysis needs
  int                            fNHPerStation[50];   //! currently - 44 x 3
  TStnCluster*                   fCluster;            //! 
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TStnTrack(int i = -1);
  ~TStnTrack();

  TLorentzVector* Momentum() { return &fMomentum; }

  TBitset*        HitMask        () { return &fHitMask; }
  TBitset*        ExpectedHitMask() { return &fExpectedHitMask; }
  
  int    Number   () const { return fNumber; }
  int    NActive  () const { return (fNActive      ) & 0xffff; }
  int    NWrong   () const { return (fNActive >> 16) & 0xffff; }
  int    NClusters();
  int    NMcStrawHits() const { return fNMcStrawHits; }
  int    NGoodMcHits () const { return fNGoodMcHits; }

  int    AlgorithmID() const { return fAlgorithmID; }
  int    BestAlg    () const { return fAlgorithmID & 0xffff; }
  int    AlgMask    () const { return (fAlgorithmID >> 16) & 0xffff; }

  float  T0       () const { return fT0;    }
  float  T0Err    () const { return fT0Err; }
  float  FitCons  () const { return fFitCons; }
  float  FitMomErr() const { return fFitMomErr; }
  float  TanDip   () const { return fTanDip; }

  float  D0       () const { return fD0; }
  float  Z0       () const { return fZ0; }
  float  Dt       () const { return fDt; }
  float  Ep       () const { return fEp; }
  float  Dy       () const { return fDy; }
  float  Dz       () const { return fDz; }
  float  Chi2Dof  () const { return fChi2/(NActive()-5+1.e-12) ; }
  float  ClusterE () const { return fClusterE;    }

  float  EleLogLHCal() const { return fEleLogLHCal; }
  float  MuoLogLHCal() const { return fMuoLogLHCal; }
  float  LogLHRCal  () const { return fEleLogLHCal-fMuoLogLHCal; }

  float  RSlope   () const { return fRSlope;    }
  float  RSlopeErr() const { return fRSlopeErr; }
  float  XSlope   () const { return fRSlope/fRSlopeErr; }

					// track-only likelihood based on timing residuals

  float  EleLogLHDeDx() const { return fEleLogLHDeDx; }
  float  MuoLogLHDeDx() const { return fMuoLogLHDeDx; }
  float  LogLHRDeDx  () const { return fEleLogLHDeDx-fMuoLogLHDeDx; }
  float  LogLHRXs    () const { return fLogLHRXs; }

  float  LogLHRTrk   () const { return LogLHRXs()+LogLHRDeDx(); }

  TStnCluster*   Cluster() const { return fCluster; }

  float  P () const { return fP; }
  float  Pt() const { return fPt; }
  Int_t  GetMomentum  (TLorentzVector* Momentum) ;
  
  KalRep*   GetKalRep() { return fKalRep[0]; }
  
  float  Phi0 () const { return fFloat[0];}
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void   SetNumber      (int I ) { fNumber      = I ; }
  void   SetAlgorithmID (int ID) { fAlgorithmID = ID; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear (Option_t* Opt = "") ;
  virtual void Print (Option_t* Opt = "") const ;

//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV4(TBuffer& R__b);
  void ReadV5(TBuffer& R__b);
  void ReadV6(TBuffer& R__b);

  ClassDef(TStnTrack,8)

};

#endif
