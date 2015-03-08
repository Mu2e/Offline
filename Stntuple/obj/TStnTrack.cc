///////////////////////////////////////////////////////////////////////////////
//
// V6 : add fIntDepth to TStnTrack::InterData_t
// V7 : introduce  curvature fC0, use one of the floats for that
///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TMatrix.h"
#include "obj/TStnTrack.hh"

namespace {
  //  const double BF = 1.4116 ;  // CDF case
  const double BF = 1.0 ;         // MU2E case
}

//namespace murat {
ClassImp(TStnTrack)
//-----------------------------------------------------------------------------
// Read an object of class TStnTrack (version 3).
// v2 didn't store the intersection data
//-----------------------------------------------------------------------------
void TStnTrack::ReadV4(TBuffer &R__b) {

  struct InterDataV4_t {
    int          fID;			// = -1 if no intersection
    float        fTime;			// track time
    float        fEnergy;		// closest cluster energy
    float        fXTrk;
    float        fYTrk;
    float        fZTrk;
    float        fXCl;
    float        fYCl;
    float        fZCl;
    float        fDx;			// TRK-CL
    float        fDy;			// TRK-CL
    float        fDz;
    float        fDt;			// TRK-CL
    float        fNxTrk;		// track direction cosines in the intersection point
    float        fNyTrk;
    float        fNzTrk;
    const mu2e::CaloCluster*       fCluster;
    const mu2e::TrkToCaloExtrapol* fExtrk;
  };

  struct TStnTrackDataV4_t {
    TLorentzVector            fMomentum;         // this assumes DELE fit hypothesis
    
    TBitset                   fHitMask;	       // bit #i: 1 if there is a hit 
    TBitset                   fExpectedHitMask;   // bit #i: 1 if expect to have a hit at this Z

    int                       fNumber;       // track index in the list of reconstructed tracks
    int                       fNHyp;          // number of hyp's with successfull fits
    int                       fBestHyp[2];    // hypothesis with the best chi2/ndof
    int                       fIDWord;	    // now - for selection "C"
    
    int                       fNActive;	    // total number of hits
    int                       fVaneID;	    // 
    int                       fDiskID;	    // 
    int                       fPdgCode;       // PDF code of the particle produced most hits
    int                       fNGoodMcHits;   // Nhits produced by the associated MC particle
    int                       fPartID;        // MC particle ID (number in the list)
    int                       fNMcStrawHits;  // Nhits by associated particle in the straw tracker
    int                       fAlgorithmID;   // bit-packed : (alg_mask << 16 ) | best
    int                       fInt[kNFreeIntsV4];     // provision for future expension
  
    float                     fChi2;
    float                     fChi2C;          // calculated...
    float                     fFitCons;
    float                     fT0;
    float                     fT0Err;
    float                     fFitMomErr;
    float                     fTanDip;
    float                     fP;		 // total momentum
    float                     fCharge;
    float                     fPt;	 // transverse momentum
    float                     fD0;
    float                     fZ0;
    
    float                     fPStOut;    // MC momentum in the VD ST_Out 
    float                     fPFront;    // MC momentum in the VD front of the tracker
    
    float                     fClusterE;	// energy of the associated cluster
    float                     fDt;
    float                     fEp;
    float                     fDx;	// about 0 for vanes
    float                     fDy;	// 
    float                     fDz;	// about 0 for disks
    
    float                     fEleLogLHCal;  // log likelihood of the electron hypothesis
    float                     fMuoLogLHCal;	// log lilelihood of the muon     hypothesis
    
    float                     fRSlope;	    // timing residual slope dres(T)/dZ
    float                     fRSlopeErr;
    
    float                     fLogLHRXs;      // XSlope-only-based likelihood
    
    float                     fEleLogLHDeDx;  // dE/dX LH calculated by Vadim based 
    float                     fMuoLogLHDeDx;  // 
    float                     fX1;	    // momentum defined at Z1
    float                     fY1;
    float                     fZ1;
    float                     fP0;            // momentum defined at Z0
    float                     fP2;            // momentum defined at Z0
    float                     fFloat[kNFreeFloatsV4]; // provision for future I/O expansion
    InterDataV4_t             fVane [kNVanesV4];      // intersection data are saved 
  };
  
  InterDataV4_t     vane;

  TStnTrackDataV4_t data;
  
  int            nwi, nwf, nwf_vint, imins, imaxep;
  
  nwi      = ((int*  ) &data.fChi2   ) - &data.fNumber;
  nwf      = ((float*) &data.fVane   ) - &data.fChi2  ;
  nwf_vint = ((float*) &vane.fCluster) - &vane.fTime  ;
    
  fMomentum.Streamer(R__b);
  fHitMask.Streamer(R__b);
  fExpectedHitMask.Streamer(R__b);

  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fChi2  ,nwf);

  fNumber           = data.fNumber      ;
  fNHyp             = data.fNHyp        ;        // number of hyp's with successfull fits
  fBestHyp[0]       = data.fBestHyp[0]  ;  // hypothesis with the best chi2/ndof
  fBestHyp[1]       = data.fBestHyp[1]  ;  // hypothesis with the best chi2/ndof
  fIDWord           = data.fIDWord      ;	    // now - for selection "C"
  fNActive          = data.fNActive     ;	    // total number of hits
  fVaneID           = data.fVaneID      ;	    // 
  fDiskID           = data.fDiskID      ;	    // 
  fPdgCode          = data.fPdgCode     ;     // PDF code of the particle produced most hits
  fNGoodMcHits      = data.fNGoodMcHits ; // 
  fPartID           = data.fPartID      ;          // MC particle ID (number in the list)
  fNMcStrawHits     = data.fNMcStrawHits;
  fAlgorithmID      = data.fAlgorithmID ;
					      // floats
  fChi2             = data.fChi2        ;
  fChi2C            = data.fChi2C       ;      // calculated...
  fFitCons          = data.fFitCons     ;
  fT0               = data.fT0          ;
  fT0Err            = data.fT0Err       ;
  fFitMomErr        = data.fFitMomErr   ;
  fTanDip           = data.fTanDip      ;
  fP                = data.fP           ;		 // total momentum
  fCharge           = data.fCharge      ;
  fPt               = data.fPt          ;	 // transverse momentum
  fD0               = data.fD0          ;
  fZ0               = data.fZ0          ;
  fPStOut           = data.fPStOut      ;    // MC momentum in the VD ST_Out 
  fPFront           = data.fPFront      ;    // MC momentum in the VD front of the tracker
  fClusterE         = data.fClusterE    ;	// energy of the associated cluster
  fDt               = data.fDt          ;
  fEp               = data.fEp          ;
  fDx               = data.fDx          ;	// about 0 for vanes
  fDy               = data.fDy          ;	// 
  fDz               = data.fDz          ;	// about 0 for disks
  fEleLogLHCal      = data.fEleLogLHCal ;  // log likelihood of the electron hypothesis
  fMuoLogLHCal      = data.fMuoLogLHCal ;	// log lilelihood of the muon     hypothesis
  fRSlope           = data.fRSlope      ;	    // timing residual slope dres(T)/dZ
  fRSlopeErr        = data.fRSlopeErr   ;
  fLogLHRXs         = data.fLogLHRXs    ;      // XSlope-only-based likelihood
  fEleLogLHDeDx     = data.fEleLogLHDeDx;  // dE/dX LH calculated by Vadim based 
  fMuoLogLHDeDx     = data.fMuoLogLHDeDx;  // 
  fX1               = data.fX1          ;	    // momentum defined at Z1
  fY1               = data.fY1          ;
  fZ1               = data.fZ1          ;
  fP0               = data.fP0          ;            // momentum defined at Z0
  fP2               = data.fP2          ;            // momentum defined at Z0
//-----------------------------------------------------------------------------
// read intersection info 
//-----------------------------------------------------------------------------
  R__b >> imins;
  R__b >> imaxep;

  for (int i=0; i<kNVanesV4; i++) {
    R__b >> vane.fID;
    R__b.ReadFastArray(&vane.fTime,nwf_vint);

    if (i < kNDisks) { 
					 // there are only 2 disks anyway...for a long time
      fDisk[i].fID       = vane.fID    ; // = -1 if no intersection
      fDisk[i].fClusterIndex = -1.;
      fDisk[i].fTime     = vane.fTime  ; // track time
      fDisk[i].fEnergy   = vane.fEnergy; // closest cluster energy
      fDisk[i].fXTrk     = vane.fXTrk  ;
      fDisk[i].fYTrk     = vane.fYTrk  ;
      fDisk[i].fZTrk     = vane.fZTrk  ;
      fDisk[i].fXCl      = vane.fXCl   ;
      fDisk[i].fYCl      = vane.fYCl   ;
      fDisk[i].fZCl      = vane.fZCl   ;
      fDisk[i].fDx       = vane.fDx    ; // TRK-CL
      fDisk[i].fDy       = vane.fDy    ; // TRK-CL
      fDisk[i].fDz       = vane.fDz    ; // 
      fDisk[i].fDt       = vane.fDt    ; // TRK-CL
      fDisk[i].fDu       = -1.e6       ; // added in V6
      fDisk[i].fDv       = -1.e6       ; // added in V6
      fDisk[i].fNxTrk    = vane.fNxTrk ; // track direction cosines in the intersection point
      fDisk[i].fNyTrk    = vane.fNyTrk ;
      fDisk[i].fNzTrk    = vane.fNzTrk ;
      fDisk[i].fChi2Match= -1.e6       ; // added in V5, undefined in V4
      fDisk[i].fChi2Time = -1.e6       ; // added in V6, undefined in V4
      fDisk[i].fPath     = -1.e6       ; // added in V5
      fDisk[i].fIntDepth = -1.e6       ; // added in V6
      fDisk[i].fCluster  = NULL;
      fDisk[i].fExtrk    = NULL;
    }
  }

  if (imins >= 0) fVMinS = &fDisk[imins];
  else            fVMinS = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fDisk[imaxep];
  else             fVMaxEp = NULL;
//-----------------------------------------------------------------------------
// introdiced in V7, initialize to produce nonsence results if used
//-----------------------------------------------------------------------------
  fC0 = -1.e12;
}


//-----------------------------------------------------------------------------
// Read an object of class TStnTrack (version 3).
// v2 didn't store the intersection data
//-----------------------------------------------------------------------------
void TStnTrack::ReadV5(TBuffer &R__b) {

  struct InterDataV5_t {
    int          fID;			// = -1 if no intersection
    float        fTime;			// track time
    float        fEnergy;		// closest cluster energy
    float        fXTrk;
    float        fYTrk;
    float        fZTrk;
    float        fXCl;
    float        fYCl;
    float        fZCl;
    float        fDx;			// TRK-CL
    float        fDy;			// TRK-CL
    float        fDz;
    float        fDt;			// TRK-CL
    float        fNxTrk;		// track direction cosines in the intersection point
    float        fNyTrk;
    float        fNzTrk;
    float        fChi2Match;		// track-cluster match chi&^2 
    float        fPath;
    const mu2e::CaloCluster*       fCluster;
    const mu2e::TrkToCaloExtrapol* fExtrk;
  };
    
  struct TStnTrackDataV5_t {
    TLorentzVector            fMomentum;         // this assumes DELE fit hypothesis

    TBitset                   fHitMask;	       // bit #i: 1 if there is a hit 
    TBitset                   fExpectedHitMask;   // bit #i: 1 if expect to have a hit at this Z

    int                       fNumber;       // track index in the list of reconstructed tracks
    int                       fNHyp;          // number of hyp's with successfull fits
    int                       fBestHyp[2];    // hypothesis with the best chi2/ndof
    int                       fIDWord;	    // now - for selection "C"
    
    int                       fNActive;	    // total number of hits
    int                       fVaneID;	    // 
    int                       fDiskID;	    // 
    int                       fPdgCode;       // PDF code of the particle produced most hits
    int                       fNGoodMcHits;   // Nhits produced by the associated MC particle
    int                       fPartID;        // MC particle ID (number in the list)
    int                       fNMcStrawHits;  // Nhits by associated particle in the straw tracker
    int                       fAlgorithmID;   // bit-packed : (alg_mask << 16 ) | best
    int                       fInt[kNFreeIntsV5];     // provision for future expension
  
    float                     fChi2;
    float                     fChi2C;          // calculated...
    float                     fFitCons;
    float                     fT0;
    float                     fT0Err;
    float                     fFitMomErr;
    float                     fTanDip;
    float                     fP;		 // total momentum
    float                     fCharge;
    float                     fPt;	 // transverse momentum
    float                     fD0;
    float                     fZ0;
    
    float                     fPStOut;    // MC momentum in the VD ST_Out 
    float                     fPFront;    // MC momentum in the VD front of the tracker
    
    float                     fClusterE;	// energy of the associated cluster
    float                     fDt;
    float                     fEp;
    float                     fDx;	// about 0 for vanes
    float                     fDy;	// 
    float                     fDz;	// about 0 for disks
    
    float                     fEleLogLHCal;  // log likelihood of the electron hypothesis
    float                     fMuoLogLHCal;	// log lilelihood of the muon     hypothesis
    
    float                     fRSlope;	    // timing residual slope dres(T)/dZ
    float                     fRSlopeErr;
    
    float                     fLogLHRXs;      // XSlope-only-based likelihood
    
    float                     fEleLogLHDeDx;  // dE/dX LH calculated by Vadim based 
    float                     fMuoLogLHDeDx;  // 
    float                     fX1;	    // momentum defined at Z1
    float                     fY1;
    float                     fZ1;
    float                     fP0;            // momentum defined at Z0
    float                     fP2;            // momentum defined at Z0
    float                     fFloat[kNFreeFloatsV5]; // provision for future I/O expansion
    
    InterData_t               fVane[kNVanesV5];  // intersection data are saved 
  };

  InterDataV5_t     vane;

  TStnTrackDataV5_t data;
  
  int               nwi, nwf, nwf_vint, imins, imaxep;
  
  nwi      = ((int*  ) &data.fChi2   ) - &data.fNumber;
  nwf      = ((float*) &data.fVane   ) - &data.fChi2  ;
  nwf_vint = ((float*) &vane.fCluster) - &vane.fTime  ;
    
  fMomentum.Streamer(R__b);
  fHitMask.Streamer(R__b);
  fExpectedHitMask.Streamer(R__b);

  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fChi2  ,nwf);


  fNumber       = data.fNumber      ; 
  fNHyp         = data.fNHyp        ;          
  fBestHyp[0]   = data.fBestHyp[0]  ;
  fBestHyp[1]   = data.fBestHyp[1]  ;
  fIDWord       = data.fIDWord      ;	
		  	       
  fNActive      = data.fNActive     ;	
  fVaneID       = data.fVaneID      ;	
  fDiskID       = data.fDiskID      ;	
  fPdgCode      = data.fPdgCode     ;       
  fNGoodMcHits  = data.fNGoodMcHits ;   
  fPartID       = data.fPartID      ;        
  fNMcStrawHits = data.fNMcStrawHits;  
  fAlgorithmID  = data.fAlgorithmID ;   
					// floats
  fChi2         = data.fChi2        ;
  fChi2C        = data.fChi2C       ;         
  fFitCons      = data.fFitCons     ;
  fT0           = data.fT0          ;
  fT0Err        = data.fT0Err       ;
  fFitMomErr    = data.fFitMomErr   ;
  fTanDip       = data.fTanDip      ;
  fP            = data.fP           ;		
  fCharge       = data.fCharge      ;
  fPt           = data.fPt          ;	 
  fD0           = data.fD0          ;
  fZ0           = data.fZ0          ;
		  	       
  fPStOut       = data.fPStOut      ;    // MC momentum in the VD ST_Out 
  fPFront       = data.fPFront      ;    // MC momentum in the VD front of the tracker
  		  	       
  fClusterE     = data.fClusterE    ;	// energy of the associated cluster
  fDt           = data.fDt          ;
  fEp           = data.fEp          ;
  fDx           = data.fDx          ;	// about 0 for vanes
  fDy           = data.fDy          ;	// 
  fDz           = data.fDz          ;	// about 0 for disks
		  	       
  fEleLogLHCal  = data.fEleLogLHCal ;  // log likelihood of the electron hypothesis
  fMuoLogLHCal  = data.fMuoLogLHCal ;	// log lilelihood of the muon     hypothesis
  		  	       
  fRSlope       = data.fRSlope      ;	    // timing residual slope dres(T)/dZ
  fRSlopeErr    = data.fRSlopeErr   ;
  		  	       
  fLogLHRXs     = data.fLogLHRXs    ;      // XSlope-only-based likelihood
		  	       
  fEleLogLHDeDx = data.fEleLogLHDeDx;  // dE/dX LH calculated by Vadim based 
  fMuoLogLHDeDx = data.fMuoLogLHDeDx;  // 
  fX1           = data.fX1          ;	    // momentum defined at Z1
  fY1           = data.fY1          ;
  fZ1           = data.fZ1          ;
  fP0           = data.fP0          ;            // momentum defined at Z0
  fP2           = data.fP2          ;            // momentum defined at Z0

//-----------------------------------------------------------------------------
// read intersection info 
//-----------------------------------------------------------------------------
  R__b >> imins;
  R__b >> imaxep;

  for (int i=0; i<kNVanesV5; i++) {
    R__b >> vane.fID;
    R__b.ReadFastArray(&vane.fTime,nwf_vint);

					 // for a long time there are only 2 disks...
    if (i < kNDisks) { 
      fDisk[i].fID       = vane.fID    ; // = -1 if no intersection
      fDisk[i].fClusterIndex = -1.;
      fDisk[i].fTime     = vane.fTime  ; // track time
      fDisk[i].fEnergy   = vane.fEnergy; // closest cluster energy
      fDisk[i].fXTrk     = vane.fXTrk  ;
      fDisk[i].fYTrk     = vane.fYTrk  ;
      fDisk[i].fZTrk     = vane.fZTrk  ;
      fDisk[i].fXCl      = vane.fXCl   ;
      fDisk[i].fYCl      = vane.fYCl   ;
      fDisk[i].fZCl      = vane.fZCl   ;
      fDisk[i].fDx       = vane.fDx    ; // TRK-CL
      fDisk[i].fDy       = vane.fDy    ; // TRK-CL
      fDisk[i].fDz       = vane.fDz    ;
      fDisk[i].fDt       = vane.fDt    ; // TRK-CL
      fDisk[i].fDu       = -1.e6       ; // added in V6
      fDisk[i].fDv       = -1.e6       ; // added in V6
      fDisk[i].fNxTrk    = vane.fNxTrk ; // track direction cosines in the intersection point
      fDisk[i].fNyTrk    = vane.fNyTrk ;
      fDisk[i].fNzTrk    = vane.fNzTrk ;
      fDisk[i].fChi2Match= -1.e6       ; // added in V5, undefined in V4
      fDisk[i].fChi2Time = -1.e6       ; // added in V6, undefined in V4
      fDisk[i].fPath     = -1.e6       ; // added in V5
      fDisk[i].fIntDepth = -1.e6       ; // added in V6
      fDisk[i].fCluster  = NULL        ;
      fDisk[i].fExtrk    = NULL        ;
    }
  }

  if (imins >= 0) fVMinS   = &fDisk[imins];
  else            fVMinS   = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fDisk[imaxep];
  else             fVMaxEp = NULL;
//-----------------------------------------------------------------------------
// introdiced in V7, initialize to produce nonsence results if used
//-----------------------------------------------------------------------------
  fC0 = -1.e12;
}


//-----------------------------------------------------------------------------
// Read an object of class TStnTrack (version 6).
//-----------------------------------------------------------------------------
void TStnTrack::ReadV6(TBuffer &R__b) {

  struct InterDataV6_t {
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
    
  struct TStnTrackDataV6_t {
    TLorentzVector            fMomentum;         // this assumes DELE fit hypothesis

    TBitset                   fHitMask;	       // bit #i: 1 if there is a hit 
    TBitset                   fExpectedHitMask;   // bit #i: 1 if expect to have a hit at this Z

    int                       fNumber;       // track index in the list of reconstructed tracks
    int                       fNHyp;          // number of hyp's with successfull fits
    int                       fBestHyp[2];    // hypothesis with the best chi2/ndof
    int                       fIDWord;	    // now - for selection "C"
    
    int                       fNActive;	    // total number of hits
    int                       fVaneID;	    // 
    int                       fDiskID;	    // 
    int                       fPdgCode;       // PDF code of the particle produced most hits
    int                       fNGoodMcHits;   // Nhits produced by the associated MC particle
    int                       fPartID;        // MC particle ID (number in the list)
    int                       fNMcStrawHits;  // Nhits by associated particle in the straw tracker
    int                       fAlgorithmID;   // bit-packed : (alg_mask << 16 ) | best
    int                       fInt[kNFreeIntsV6];     // provision for future expension
  
    float                     fChi2;
    float                     fChi2C;          // calculated...
    float                     fFitCons;
    float                     fT0;
    float                     fT0Err;
    float                     fFitMomErr;
    float                     fTanDip;
    float                     fP;		 // total momentum
    float                     fCharge;
    float                     fPt;	 // transverse momentum
    float                     fD0;
    float                     fZ0;
    
    float                     fPStOut;    // MC momentum in the VD ST_Out 
    float                     fPFront;    // MC momentum in the VD front of the tracker
    
    float                     fClusterE;	// energy of the associated cluster
    float                     fDt;
    float                     fEp;
    float                     fDx;	// about 0 for vanes
    float                     fDy;	// 
    float                     fDz;	// about 0 for disks
    
    float                     fEleLogLHCal;  // log likelihood of the electron hypothesis
    float                     fMuoLogLHCal;	// log lilelihood of the muon     hypothesis
    
    float                     fRSlope;	    // timing residual slope dres(T)/dZ
    float                     fRSlopeErr;
    
    float                     fLogLHRXs;      // XSlope-only-based likelihood
    
    float                     fEleLogLHDeDx;  // dE/dX LH calculated by Vadim based 
    float                     fMuoLogLHDeDx;  // 
    float                     fX1;	    // momentum defined at Z1
    float                     fY1;
    float                     fZ1;
    float                     fP0;            // momentum defined at Z0
    float                     fP2;            // momentum defined at Z0
    float                     fFloat[kNFreeFloatsV6]; // provision for future I/O expansion
    
    InterData_t               fVane[kNDisks];  // intersection data are saved 
  };

  InterDataV6_t     disk;

  TStnTrackDataV6_t data;
  
  int               nwi, nwf, nwf_vint, imins, imaxep;
  
  nwi      = ((int*  ) &data.fChi2   ) - &data.fNumber;
  nwf      = ((float*) &data.fVane   ) - &data.fChi2  ;
  nwf_vint = ((float*) &disk.fCluster) - &disk.fTime  ;
    
  fMomentum.Streamer(R__b);
  fHitMask.Streamer(R__b);
  fExpectedHitMask.Streamer(R__b);

  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fChi2  ,nwf);


  fNumber       = data.fNumber      ; 
  fNHyp         = data.fNHyp        ;          
  fBestHyp[0]   = data.fBestHyp[0]  ;
  fBestHyp[1]   = data.fBestHyp[1]  ;
  fIDWord       = data.fIDWord      ;	
		  	       
  fNActive      = data.fNActive     ;	
  fVaneID       = data.fVaneID      ;	
  fDiskID       = data.fDiskID      ;	
  fPdgCode      = data.fPdgCode     ;       
  fNGoodMcHits  = data.fNGoodMcHits ;   
  fPartID       = data.fPartID      ;        
  fNMcStrawHits = data.fNMcStrawHits;  
  fAlgorithmID  = data.fAlgorithmID ;   
					// floats
  fChi2         = data.fChi2        ;
  fChi2C        = data.fChi2C       ;         
  fFitCons      = data.fFitCons     ;
  fT0           = data.fT0          ;
  fT0Err        = data.fT0Err       ;
  fFitMomErr    = data.fFitMomErr   ;
  fTanDip       = data.fTanDip      ;
  fP            = data.fP           ;		
  fCharge       = data.fCharge      ;
  fPt           = data.fPt          ;	 
  fD0           = data.fD0          ;
  fZ0           = data.fZ0          ;
		  	       
  fPStOut       = data.fPStOut      ;    // MC momentum in the VD ST_Out 
  fPFront       = data.fPFront      ;    // MC momentum in the VD front of the tracker
  		  	       
  fClusterE     = data.fClusterE    ;	// energy of the associated cluster
  fDt           = data.fDt          ;
  fEp           = data.fEp          ;
  fDx           = data.fDx          ;	// about 0 for vanes
  fDy           = data.fDy          ;	// 
  fDz           = data.fDz          ;	// about 0 for disks
		  	       
  fEleLogLHCal  = data.fEleLogLHCal ;  // log likelihood of the electron hypothesis
  fMuoLogLHCal  = data.fMuoLogLHCal ;	// log lilelihood of the muon     hypothesis
  		  	       
  fRSlope       = data.fRSlope      ;	    // timing residual slope dres(T)/dZ
  fRSlopeErr    = data.fRSlopeErr   ;
  		  	       
  fLogLHRXs     = data.fLogLHRXs    ;      // XSlope-only-based likelihood
		  	       
  fEleLogLHDeDx = data.fEleLogLHDeDx;  // dE/dX LH calculated by Vadim based 
  fMuoLogLHDeDx = data.fMuoLogLHDeDx;  // 
  fX1           = data.fX1          ;	    // momentum defined at Z1
  fY1           = data.fY1          ;
  fZ1           = data.fZ1          ;
  fP0           = data.fP0          ;            // momentum defined at Z0
  fP2           = data.fP2          ;            // momentum defined at Z0
//-----------------------------------------------------------------------------
// read intersection info 
//-----------------------------------------------------------------------------
  R__b >> imins;
  R__b >> imaxep;

  for (int i=0; i<kNDisks; i++) {
    R__b >> disk.fID;
    R__b >> disk.fClusterIndex;
    R__b.ReadFastArray(&disk.fTime,nwf_vint);

					 // for a long time there are only 2 disks...
    if (i < kNDisks) { 
      fDisk[i].fID       = disk.fID    ; // = -1 if no intersection
      fDisk[i].fClusterIndex = disk.fClusterIndex;
      fDisk[i].fTime     = disk.fTime  ; // track time
      fDisk[i].fEnergy   = disk.fEnergy; // closest cluster energy
      fDisk[i].fXTrk     = disk.fXTrk  ;
      fDisk[i].fYTrk     = disk.fYTrk  ;
      fDisk[i].fZTrk     = disk.fZTrk  ;
      fDisk[i].fXCl      = disk.fXCl   ;
      fDisk[i].fYCl      = disk.fYCl   ;
      fDisk[i].fZCl      = disk.fZCl   ;
      fDisk[i].fDx       = disk.fDx    ; // TRK-CL
      fDisk[i].fDy       = disk.fDy    ; // TRK-CL
      fDisk[i].fDz       = disk.fDz    ;
      fDisk[i].fDt       = disk.fDt    ; // TRK-CL
      fDisk[i].fDu       = disk.fDu    ; // TRK-CL, added in V6
      fDisk[i].fDv       = disk.fDv    ; // TRK-CL, added in V6
      fDisk[i].fNxTrk    = disk.fNxTrk ; // track direction cosines in the intersection point
      fDisk[i].fNyTrk    = disk.fNyTrk ;
      fDisk[i].fNzTrk    = disk.fNzTrk ;
      fDisk[i].fChi2Match= disk.fChi2Match; 
      fDisk[i].fChi2Time = disk.fChi2Time; 
      fDisk[i].fPath     = disk.fPath  ;
      fDisk[i].fIntDepth = disk.fIntDepth;
      fDisk[i].fCluster  = NULL        ;
      fDisk[i].fExtrk    = NULL        ;
    }
  }

  if (imins >= 0) fVMinS   = &fDisk[imins];
  else            fVMinS   = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fDisk[imaxep];
  else             fVMaxEp = NULL;
//-----------------------------------------------------------------------------
// undefined in V6, set to produce nonsence results if used
//-----------------------------------------------------------------------------
  fC0 = -1.e12;
}



//-----------------------------------------------------------------------------
void TStnTrack::Streamer(TBuffer& R__b) {

  int nwi, nwf, nwf_vint, imins, imaxep;

  nwi      = ((int*  ) &fChi2            ) - &fNumber;
  nwf      = ((float*) &fDisk            ) - &fChi2;
  nwf_vint = ((float*) &fDisk[0].fCluster) - &fDisk[0].fTime;

  if (R__b.IsReading()) {
//-----------------------------------------------------------------------------
// read TStnTrack, Mu2e:V1
//-----------------------------------------------------------------------------
    Version_t R__v = R__b.ReadVersion(); 

    if      (R__v < 4) {
      printf(" >>> ERROR: TStnTrack::Streamer can't read old data version = %i. BAIL OUT\n",R__v);
    }
    else if (R__v == 4) ReadV4(R__b);
    else if (R__v == 5) ReadV5(R__b);
    else if (R__v == 6) ReadV6(R__b);
    else {
//-----------------------------------------------------------------------------
// current version: v8 
// I/O for V7 is exactly the same
//-----------------------------------------------------------------------------
      fMomentum.Streamer(R__b);
      fHitMask.Streamer(R__b);
      fExpectedHitMask.Streamer(R__b);

      R__b.ReadFastArray(&fNumber,nwi);
      R__b.ReadFastArray(&fChi2,nwf);
					// read intersection info
      R__b >> imins;
      R__b >> imaxep;

      for (int i=0; i<kNDisks; i++) {
	R__b >> fDisk[i].fID;
	R__b >> fDisk[i].fClusterIndex;
	R__b.ReadFastArray(&fDisk[i].fTime,nwf_vint);
      }
    
      if (imins >= 0) fVMinS = &fDisk[imins];
      else            fVMinS = NULL;

      if (imaxep >= 0) fVMaxEp = &fDisk[imaxep];
      else             fVMaxEp = NULL;
    }
  }
  else {
//-----------------------------------------------------------------------------
// write track data out
//-----------------------------------------------------------------------------
    R__b.WriteVersion(TStnTrack::IsA());
    fMomentum.Streamer(R__b);
    fHitMask.Streamer(R__b);
    fExpectedHitMask.Streamer(R__b);

    R__b.WriteFastArray(&fNumber,nwi);
    R__b.WriteFastArray(&fChi2,nwf);
					// write out intersection info
    if (fVMinS != NULL) imins  = fVMinS-fDisk;
    else                imins  = -1;

    if (fVMaxEp != NULL) imaxep = fVMaxEp-fDisk;
    else                 imaxep = -1;

    R__b << imins;
    R__b << imaxep;
    
    for (int i=0; i<kNDisks; i++) {
      R__b << fDisk[i].fID;
      R__b << fDisk[i].fClusterIndex;
      R__b.WriteFastArray(&fDisk[i].fTime,nwf_vint);
    }
  }
}

//_____________________________________________________________________________
TStnTrack::TStnTrack(Int_t Number) : TObject (),
				     fHitMask(kMaxNLayers),
				     fExpectedHitMask(kMaxNLayers)
{
    // 'Number' can be -1 ...

  
  Clear();

  fNumber = Number;
}


//_____________________________________________________________________________
TStnTrack::~TStnTrack() {
}


//_____________________________________________________________________________
int TStnTrack::NClusters() {
  int ncl(0);
  for (int iv=0; iv<kNDisks; iv++) {
    if (fDisk[iv].fCluster != NULL) {
      ncl++;
    }
  }
  return ncl;
}

//_____________________________________________________________________________
void TStnTrack::Clear(Option_t* Opt) {
  fNumber    = -1;

  fHitMask.Clear();
  fExpectedHitMask.Clear();

  fKalRep[0] = 0;
  fKalRep[1] = 0;
  fKalRep[2] = 0;
  fKalRep[3] = 0;
  
  fExtrk              = NULL;
  fClosestCaloCluster = NULL;
  fCluster            = NULL;
    
  fNHyp           = -1;
  fBestHyp[0]     = -1;
  fBestHyp[1]     = -1;

  for (int i=0; i<kNDisks; i++) {
    fDisk[i].fID           = -1;
    fDisk[i].fClusterIndex = -1;
    fDisk[i].fEnergy       = -1.;

    fDisk[i].fXTrk         = -1.e6;
    fDisk[i].fYTrk         = -1.e6;
    fDisk[i].fZTrk         = -1.e6;
    fDisk[i].fTime         = -1.e6;
    fDisk[i].fNxTrk        = -1.e6;
    fDisk[i].fNyTrk        = -1.e6;
    fDisk[i].fNzTrk        = -1.e6;
			   
    fDisk[i].fXCl          = -1.e6;
    fDisk[i].fYCl          = -1.e6;
    fDisk[i].fZCl          = -1.e6;
			   
    fDisk[i].fDx           = -1.e6;
    fDisk[i].fDy           = -1.e6;
    fDisk[i].fDz           = -1.e6;
    fDisk[i].fDt           = -1.e6;
			   
    fDisk[i].fDu           = -1.e6;
    fDisk[i].fDv           = -1.e6;
			   
    fDisk[i].fChi2Match    = 1.e6;
    fDisk[i].fChi2Time     = 1.e6;
    fDisk[i].fPath         = 1.e6;
    fDisk[i].fIntDepth     = 1.e6;
    fDisk[i].fCluster      = NULL;
    fDisk[i].fExtrk        = NULL;
  }
  
  fVMinS        = NULL;
  fVMaxEp       = NULL;
					// by definition, LogLH < 0
  fEleLogLHCal  =  1.;
  fMuoLogLHCal  =  1.;
  fEleLogLHDeDx =  1.;
  fMuoLogLHDeDx =  1.;
  fRSlope       = -1.e6;
  fRSlopeErr    = -1.;
					// for the same reason, E/P > 0
  fEp           = -1.;
}

//_____________________________________________________________________________
void TStnTrack::Print(Option_t* opt) const {
  printf("TStnTrack::Print <WARNING> Not implemented yet\n");
}

// } // end namespace
