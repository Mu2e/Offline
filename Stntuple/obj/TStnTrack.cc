//
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
// Read an object of class TStnTrack (version 1).
// v1 didn't store the track direction
//-----------------------------------------------------------------------------
void TStnTrack::ReadV1(TBuffer &R__b) {

  struct InterDataV1_t {
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
    const mu2e::CaloCluster*       fCluster;
    const mu2e::TrkToCaloExtrapol* fExtrk;
  };

  struct TStnTrackDataV1_t {
    int                       fNumber;       // track index in the list of reconstructed tracks
    int                       fNHPerStation[40];
    int                       fNHyp;        // number of hyp's with successfull fits
    int                       fBestHyp[2];  // hypothesis with the best chi2/ndof
    int                       fIDWord;	    // now - for selection "C"
  
    int                       fNActive;	    // total number of hits
    int                       fVaneID;	    // 
    int                       fDiskID;	    // 
    int                       fPdgCode;     // PDF code of the particle produced most hits
    int                       fNGoodMcHits; // 
    int                       fPartID;          // MC particle ID (number in the list)
    int                       fInt[kNFreeIntsV1];     // provision for future expension
    
    float                     fChi2;
    float                     fChi2C;      // calculated...
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
    float                     fFloat[kNFreeFloatsV1]; // provision for future I/O expansion
    InterDataV1_t             fVane[kNVanes];  // intersection data are saved 
  };

  TStnTrackDataV1_t data;
  InterDataV1_t     vane;
  int               nwi, nwf, nwf_vint, imins, imaxep;
  
  nwi      = ((int*  ) &data.fChi2   ) - &data.fNumber;
  nwf      = ((float*) &data.fVane   ) - &data.fChi2  ;
  nwf_vint = ((float*) &vane.fCluster) - &vane.fTime  ;
    
  fMomentum.Streamer(R__b);

  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fChi2  ,nwf);

  fNumber           = data.fNumber          ;
  memcpy(fNHPerStation,data.fNHPerStation,40*sizeof(int));
  fNHyp             = data.fNHyp            ;        // number of hyp's with successfull fits
  fBestHyp[2]       = data.fBestHyp[2]      ;  // hypothesis with the best chi2/ndof
  fIDWord           = data.fIDWord          ;	    // now - for selection "C"
  fNActive          = data.fNActive         ;	    // total number of hits
  fVaneID           = data.fVaneID          ;	    // 
  fDiskID           = data.fDiskID          ;	    // 
  fPdgCode          = data.fPdgCode         ;     // PDF code of the particle produced most hits
  fNGoodMcHits      = data.fNGoodMcHits     ; // 
  fPartID           = data.fPartID          ;          // MC particle ID (number in the list)
//

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
//-----------------------------------------------------------------------------
// read intersection info
//-----------------------------------------------------------------------------
  R__b >> imins;
  R__b >> imaxep;

  for (int i=0; i<kNVanes; i++) {
    R__b >> fVane[i].fID;
    R__b.ReadFastArray(&fVane[i].fTime,nwf_vint);
//-----------------------------------------------------------------------------
// set undefined fields to zero
//-----------------------------------------------------------------------------
    fVane[i].fNxTrk     = -1.e6;
    fVane[i].fNyTrk     = -1.e6;
    fVane[i].fNzTrk     = -1.e6;
    fVane[i].fChi2Match = -1.e6;
    fVane[i].fPath      = -1.e6;
    fVane[i].fCluster  = NULL;
    fVane[i].fExtrk    = NULL;
  }

  if (imins >= 0) fVMinS = &fVane[imins];
  else            fVMinS = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fVane[imaxep];
  else             fVMaxEp = NULL;
//-----------------------------------------------------------------------------
// undefined in V1:
//-----------------------------------------------------------------------------
  fP0 = fP;
  fX1 = 1e12;
  fY1 = 1e12;
  fZ1 = fZ0;
//-----------------------------------------------------------------------------
// undefined in V2:
// set bitset bits to 1 to preserve efficiency of all cuts
//-----------------------------------------------------------------------------
  for (int i=0; i<TStnTrack::kMaxNLayers; i++) {
    fHitMask.SetBit(i);
    fExpectedHitMask.SetBit(i);
  }
}

//-----------------------------------------------------------------------------
// Read an object of class TStnTrack (version 2).
// v1 didn't store the track direction
//-----------------------------------------------------------------------------
void TStnTrack::ReadV2(TBuffer &R__b) {

  struct InterDataV2_t {
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

  struct TStnTrackDataV2_t {
    TLorentzVector            fMomentum;         // this assumes DELE fit hypothesis

    int                       fNumber;       // track index in the list of reconstructed tracks
    int                       fNHPerStation[40];
    int                       fNHyp;        // number of hyp's with successfull fits
    int                       fBestHyp[2];  // hypothesis with the best chi2/ndof
    int                       fIDWord;	    // now - for selection "C"
    
    int                       fNActive;	    // total number of hits
    int                       fVaneID;	    // 
    int                       fDiskID;	    // 
    int                       fPdgCode;     // PDF code of the particle produced most hits
    int                       fNGoodMcHits; // 
    int                       fPartID;          // MC particle ID (number in the list)
    int                       fInt[kNFreeIntsV2];     // provision for future expension
    
    float                     fChi2;
    float                     fChi2C;      // calculated...
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
    float                     fFloat[kNFreeFloatsV2]; // provision for future I/O expansion
    InterDataV2_t             fVane[kNVanes];  // intersection data are saved 
  };

  InterDataV2_t     vane;

  TStnTrackDataV2_t data;
  
  int            nwi, nwf, nwf_vint, imins, imaxep;
  
  nwi      = ((int*  ) &data.fChi2   ) - &data.fNumber;
  nwf      = ((float*) &data.fVane   ) - &data.fChi2  ;
  nwf_vint = ((float*) &vane.fCluster) - &vane.fTime  ;
    
  fMomentum.Streamer(R__b);

  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fChi2  ,nwf);

  fNumber           = data.fNumber          ;
  memcpy(fNHPerStation,data.fNHPerStation,40*sizeof(int));
  fNHyp             = data.fNHyp            ; // number of hyp's with successfull fits
  fBestHyp[0]       = data.fBestHyp[0]      ; // hypothesis with the best chi2/ndof
  fBestHyp[1]       = data.fBestHyp[1]      ; // hypothesis with the best chi2/ndof
  fIDWord           = data.fIDWord          ; // now - for selection "C"
  fNActive          = data.fNActive         ; // total number of hits
  fVaneID           = data.fVaneID          ; // 
  fDiskID           = data.fDiskID          ; // 
  fPdgCode          = data.fPdgCode         ; // PDF code of the particle produced most hits
  fNGoodMcHits      = data.fNGoodMcHits     ; // 
  fPartID           = data.fPartID          ; // MC particle ID (number in the list)
//

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

  for (int i=0; i<kNVanes; i++) {
    R__b >> vane.fID;
    R__b.ReadFastArray(&vane.fTime,nwf_vint);

    fVane[i].fID       = vane.fID    ;			// = -1 if no intersection
    fVane[i].fTime     = vane.fTime  ;			// track time
    fVane[i].fEnergy   = vane.fEnergy;		// closest cluster energy
    fVane[i].fXTrk     = vane.fXTrk  ;
    fVane[i].fYTrk     = vane.fYTrk  ;
    fVane[i].fZTrk     = vane.fZTrk  ;
    fVane[i].fXCl      = vane.fXCl   ;
    fVane[i].fYCl      = vane.fYCl   ;
    fVane[i].fZCl      = vane.fZCl   ;
    fVane[i].fDx       = vane.fDx    ;			// TRK-CL
    fVane[i].fDy       = vane.fDy    ;			// TRK-CL
    fVane[i].fDz       = vane.fDz    ;
    fVane[i].fDt       = vane.fDt    ;			// TRK-CL
    fVane[i].fNxTrk    = vane.fNxTrk ;		// track direction cosines in the intersection point
    fVane[i].fNyTrk    = vane.fNyTrk ;
    fVane[i].fNzTrk    = vane.fNzTrk ;
    fVane[i].fChi2Match= -1.e6;		// added in V5
    fVane[i].fPath     = -1.e6       ;	// added in V5
    fVane[i].fCluster  = NULL;
    fVane[i].fExtrk    = NULL;
  }

  if (imins >= 0) fVMinS = &fVane[imins];
  else            fVMinS = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fVane[imaxep];
  else             fVMaxEp = NULL;
//-----------------------------------------------------------------------------
// undefined in V2:
// set bitset bits to 1 to preserve efficiency of all cuts
//-----------------------------------------------------------------------------
  for (int i=0; i<TStnTrack::kMaxNLayers; i++) {
    fHitMask.SetBit(i);
    fExpectedHitMask.SetBit(i);
  }
}


//-----------------------------------------------------------------------------
// Read an object of class TStnTrack (version 3).
// v2 didn't store the intersection data
//-----------------------------------------------------------------------------
void TStnTrack::ReadV3(TBuffer &R__b) {

  struct InterDataV3_t {
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

  struct TStnTrackDataV3_t {
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
    int                       fNGoodMcHits;   // 
    int                       fPartID;        // MC particle ID (number in the list)
    int                       fInt[kNFreeIntsV3];     // provision for future expension
    
    float                     fChi2;
    float                     fChi2C;      // calculated...
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
    float                     fFloat[kNFreeFloatsV3]; // provision for future I/O expansion
    
    InterDataV3_t               fVane[kNVanes];  // intersection data are saved 
  };
  
  InterDataV3_t     vane;

  TStnTrackDataV3_t data;
  
  int            nwi, nwf, nwf_vint, imins, imaxep;
  
  nwi      = ((int*  ) &data.fChi2   ) - &data.fNumber;
  nwf      = ((float*) &data.fVane   ) - &data.fChi2  ;
  nwf_vint = ((float*) &vane.fCluster) - &vane.fTime  ;
    
  fMomentum.Streamer(R__b);
  fHitMask.Streamer(R__b);
  fExpectedHitMask.Streamer(R__b);

  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fChi2  ,nwf);

  fNumber           = data.fNumber          ;
  fNHyp             = data.fNHyp            ;        // number of hyp's with successfull fits
  fBestHyp[0]       = data.fBestHyp[0]      ;  // hypothesis with the best chi2/ndof
  fBestHyp[1]       = data.fBestHyp[1]      ;  // hypothesis with the best chi2/ndof
  fIDWord           = data.fIDWord          ;	    // now - for selection "C"
  fNActive          = data.fNActive         ;	    // total number of hits
  fVaneID           = data.fVaneID          ;	    // 
  fDiskID           = data.fDiskID          ;	    // 
  fPdgCode          = data.fPdgCode         ;     // PDF code of the particle produced most hits
  fNGoodMcHits      = data.fNGoodMcHits     ; // 
  fPartID           = data.fPartID          ;          // MC particle ID (number in the list)
//

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

  for (int i=0; i<kNVanes; i++) {
    R__b >> vane.fID;
    R__b.ReadFastArray(&vane.fTime,nwf_vint);

    fVane[i].fID       = vane.fID    ;			// = -1 if no intersection
    fVane[i].fTime     = vane.fTime  ;			// track time
    fVane[i].fEnergy   = vane.fEnergy;		// closest cluster energy
    fVane[i].fXTrk     = vane.fXTrk  ;
    fVane[i].fYTrk     = vane.fYTrk  ;
    fVane[i].fZTrk     = vane.fZTrk  ;
    fVane[i].fXCl      = vane.fXCl   ;
    fVane[i].fYCl      = vane.fYCl   ;
    fVane[i].fZCl      = vane.fZCl   ;
    fVane[i].fDx       = vane.fDx    ;			// TRK-CL
    fVane[i].fDy       = vane.fDy    ;			// TRK-CL
    fVane[i].fDz       = vane.fDz    ;
    fVane[i].fDt       = vane.fDt    ;			// TRK-CL
    fVane[i].fNxTrk    = vane.fNxTrk ;		// track direction cosines in the intersection point
    fVane[i].fNyTrk    = vane.fNyTrk ;
    fVane[i].fNzTrk    = vane.fNzTrk ;
    fVane[i].fChi2Match= -1.e6;		// added in V5
    fVane[i].fPath     = -1.e6;		// added in V5
    fVane[i].fCluster  = NULL;
    fVane[i].fExtrk    = NULL;
  }

  if (imins >= 0) fVMinS = &fVane[imins];
  else            fVMinS = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fVane[imaxep];
  else             fVMaxEp = NULL;
//-----------------------------------------------------------------------------
// undefined in V3: fAlgorithmID - mask << 16 | best ; (set best to 0)
// it was on eof unused integers
//-----------------------------------------------------------------------------
  fAlgorithmID = (1 << 16) | 0x0000;
}


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
    int                       fInt[kNFreeInts];     // provision for future expension
  
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
    float                     fFloat[kNFreeFloats]; // provision for future I/O expansion
    InterDataV4_t             fVane[kNVanes];  // intersection data are saved 
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

  fNumber           = data.fNumber          ;
  fNHyp             = data.fNHyp            ;        // number of hyp's with successfull fits
  fBestHyp[0]       = data.fBestHyp[0]      ;  // hypothesis with the best chi2/ndof
  fBestHyp[1]       = data.fBestHyp[1]      ;  // hypothesis with the best chi2/ndof
  fIDWord           = data.fIDWord          ;	    // now - for selection "C"
  fNActive          = data.fNActive         ;	    // total number of hits
  fVaneID           = data.fVaneID          ;	    // 
  fDiskID           = data.fDiskID          ;	    // 
  fPdgCode          = data.fPdgCode         ;     // PDF code of the particle produced most hits
  fNGoodMcHits      = data.fNGoodMcHits     ; // 
  fPartID           = data.fPartID          ;          // MC particle ID (number in the list)
  fNMcStrawHits     = data.fNMcStrawHits    ;
  fAlgorithmID      = data.fAlgorithmID     ;
//

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

  for (int i=0; i<kNVanes; i++) {
    R__b >> vane.fID;
    R__b.ReadFastArray(&vane.fTime,nwf_vint);

    fVane[i].fID       = vane.fID    ;			// = -1 if no intersection
    fVane[i].fTime     = vane.fTime  ;			// track time
    fVane[i].fEnergy   = vane.fEnergy;		// closest cluster energy
    fVane[i].fXTrk     = vane.fXTrk  ;
    fVane[i].fYTrk     = vane.fYTrk  ;
    fVane[i].fZTrk     = vane.fZTrk  ;
    fVane[i].fXCl      = vane.fXCl   ;
    fVane[i].fYCl      = vane.fYCl   ;
    fVane[i].fZCl      = vane.fZCl   ;
    fVane[i].fDx       = vane.fDx    ;			// TRK-CL
    fVane[i].fDy       = vane.fDy    ;			// TRK-CL
    fVane[i].fDz       = vane.fDz    ;
    fVane[i].fDt       = vane.fDt    ;			// TRK-CL
    fVane[i].fNxTrk    = vane.fNxTrk ;		// track direction cosines in the intersection point
    fVane[i].fNyTrk    = vane.fNyTrk ;
    fVane[i].fNzTrk    = vane.fNzTrk ;
    fVane[i].fChi2Match= -1.e6       ;  // added in V5, undefined in V4
    fVane[i].fPath     = -1.e6       ;	// added in V5
    fVane[i].fCluster  = NULL;
    fVane[i].fExtrk    = NULL;
  }

  if (imins >= 0) fVMinS = &fVane[imins];
  else            fVMinS = NULL;
  
  if (imaxep >= 0) fVMaxEp = &fVane[imaxep];
  else             fVMaxEp = NULL;
}


//-----------------------------------------------------------------------------
void TStnTrack::Streamer(TBuffer& R__b) {

  int nwi, nwf, nwf_vint, imins, imaxep;

  nwi      = ((int*  ) &fChi2            ) - &fNumber;
  nwf      = ((float*) &fVane            ) - &fChi2;
  nwf_vint = ((float*) &fVane[0].fCluster) - &fVane[0].fTime;

  if (R__b.IsReading()) {
//-----------------------------------------------------------------------------
// read TStnTrack, Mu2e:V1
//-----------------------------------------------------------------------------
    Version_t R__v = R__b.ReadVersion(); 

    if      (R__v == 1) ReadV1(R__b);
    else if (R__v == 2) ReadV2(R__b);
    else if (R__v == 3) ReadV3(R__b);
    else if (R__v == 4) ReadV4(R__b);
    else {
//-----------------------------------------------------------------------------
// current version: v3
//-----------------------------------------------------------------------------
      fMomentum.Streamer(R__b);
      fHitMask.Streamer(R__b);
      fExpectedHitMask.Streamer(R__b);

      R__b.ReadFastArray(&fNumber,nwi);
      R__b.ReadFastArray(&fChi2,nwf);
					// read intersection info
      R__b >> imins;
      R__b >> imaxep;

      for (int i=0; i<kNVanes; i++) {
	R__b >> fVane[i].fID;
	R__b.ReadFastArray(&fVane[i].fTime,nwf_vint);
      }
    
      if (imins >= 0) fVMinS = &fVane[imins];
      else            fVMinS = NULL;

      if (imaxep >= 0) fVMaxEp = &fVane[imaxep];
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
    if (fVMinS != NULL) imins  = fVMinS-fVane;
    else                imins  = -1;

    if (fVMaxEp != NULL) imaxep = fVMaxEp-fVane;
    else                 imaxep = -1;

    R__b << imins;
    R__b << imaxep;
    
    for (int i=0; i<kNVanes; i++) {
      R__b << fVane[i].fID;
      R__b.WriteFastArray(&fVane[i].fTime,nwf_vint);
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
  for (int iv=0; iv<kNVanes; iv++) {
    if (fVane[iv].fCluster != NULL) {
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
  
  fExtrk          = NULL;
  fClosestCluster = NULL;
    
  fNHyp           = -1;
  fBestHyp[0]     = -1;
  fBestHyp[1]     = -1;

  for (int i=0; i<kNVanes; i++) {
    fVane[i].fID      = -1;
    fVane[i].fEnergy  = -1.;
    fVane[i].fXTrk    = -1.e6;
    fVane[i].fYTrk    = -1.e6;
    fVane[i].fZTrk    = -1.e6;
    fVane[i].fTime    = -1.e6;
    fVane[i].fXCl     = -1.e6;
    fVane[i].fYCl     = -1.e6;
    fVane[i].fZCl     = -1.e6;
    fVane[i].fDx      = -1.e6;
    fVane[i].fDy      = -1.e6;
    fVane[i].fDz      = -1.e6;
    fVane[i].fDt      = -1.e6;
    fVane[i].fNxTrk   = -1.e6;
    fVane[i].fNyTrk   = -1.e6;
    fVane[i].fNzTrk   = -1.e6;
    fVane[i].fCluster = NULL;
    fVane[i].fExtrk   = NULL;
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




