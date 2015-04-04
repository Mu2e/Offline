///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <cmath>

#include "TH1F.h"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/obj/TStnTrack.hh"

ClassImp(TStnTrackID)

//______________________________________________________________________________
void TStnTrackID::Streamer(TBuffer &R__b) {
// //  // Stream an object of class TStnTrackID.
// //
  UInt_t R__s, R__c;
// //
// //  int nwi = ((Int_t*  )&fMinFitCons)-&fUseMask;
// //  int nwf = ((Float_t*)&fEOR       )-&fMinFitCons  ;
// //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
// //    TNamed::Streamer(R__b);
// //    R__b.ReadFastArray(&fUseMask   ,nwi);
// //    R__b.ReadFastArray(&fMinFitCons,nwf);
// //    R__b.CheckByteCount(R__s, R__c, TStnTrackID::IsA());
  } 
  else {
//-----------------------------------------------------------------------------
//  write section
//-----------------------------------------------------------------------------
    R__c = R__b.WriteVersion(TStnTrackID::IsA(), kTRUE);
// //    TNamed::Streamer(R__b);
// //    R__b.WriteFastArray(&fUseMask   ,nwi);
// //    R__b.WriteFastArray(&fMinFitCons,nwf);
// //    R__b.SetByteCount(R__c, kTRUE);
  }
}

//-----------------------------------------------------------------------------
// default: "Set C" cuts , as defined in 
// http://mu2e-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=2447;filename=Resolutions.pdf;version=2
//
// 2014-04-02: update according to 
//    http://mu2e-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=3996;filename=cutset_c.pdf;version=2
//-----------------------------------------------------------------------------
TStnTrackID::TStnTrackID(const char* Name): TNamed(Name,Name) {
  
  fUseMask         = 0xffffffff;
  fMinFitCons      = 2.e-3;
  fMinT0           = 700.;		// ns
  fMinNActive      = 25;
  fMaxT0Err        = 0.9;  		// in ns
  fMaxFitMomErr    = 0.25;  		// in MeV
  fMinTanDip       = tan(M_PI/6.);	// 0.5773
  fMaxTanDip       = 1.0;  
  fMinD1           = -80.;		// in mm
  fMaxD1           = 105.;

  fMinD2           = -1.e6;
  fMaxD2           =  1.e6;
				// initialize spare words to zero

  for (int i=0; i< 5; i++) fInteger[i] = 0;
  for (int i=0; i<10; i++) fFloat  [i] = 0;
}


//_____________________________________________________________________________
TStnTrackID::~TStnTrackID() {
}

//-----------------------------------------------------------------------------
//  need to implement the D0 cuts
//-----------------------------------------------------------------------------
int TStnTrackID::IDWord(TStnTrack* Track) {
  Int_t      id_word = 0;
  int        nactive;
  double     fcons, d0, t0, t0_err, fitmom_err, tan_dip, two_over_omega, rmax;
  
  //  KalRep* trk = Track->fKalRep[0];

  fcons      = Track->fFitCons;
  t0         = Track->fT0;
  t0_err     = Track->fT0Err;
  nactive    = Track->NActive();
  tan_dip    = Track->fTanDip;
  fitmom_err = Track->fFitMomErr;
  d0         = Track->fD0;
					// so far - kludge
  two_over_omega = 0. ;			// 2/c0 - signed diameter
  rmax           = d0+two_over_omega;

  if (fcons            < fMinFitCons        ) id_word |= kFitConsBit ;
  if (t0               < fMinT0             ) id_word |= kT0Bit ;
  if (t0_err           > fMaxT0Err          ) id_word |= kT0ErrBit ;
  if (nactive          < fMinNActive        ) id_word |= kNActiveBit ;
  if (fitmom_err       > fMaxFitMomErr      ) id_word |= kFitMomErrBit ;
  if (tan_dip          < fMinTanDip         ) id_word |= kTanDipBit ;
  if (tan_dip          > fMaxTanDip         ) id_word |= kTanDipBit ;

  if (d0               < fMinD1             ) id_word |= kD1Bit ;
  if (d0               > fMaxD1             ) id_word |= kD1Bit ;

  if (rmax             < fMinD2             ) id_word |= kD2Bit ;
  if (rmax             > fMaxD2             ) id_word |= kD2Bit ;

  return  (id_word & fUseMask);
}

//_____________________________________________________________________________
int TStnTrackID::TightIDWord(TStnTrack* Track) {
  int  id_word;
  id_word = IDWord(Track);
  return id_word;
}


//-----------------------------------------------------------------------------
// so far, no use
//-----------------------------------------------------------------------------
int TStnTrackID::LooseIDWord(TStnTrack* Track) {
  return IDWord(Track);
}


//_____________________________________________________________________________
void TStnTrackID::FillHistograms(Hist_t& Hist, TStnTrack* Track, Int_t Mode) {
//-----------------------------------------------------------------------------
// distributions for ID variables (with n-1 ID cuts applied). 
// Plot 3 sets of histograms:
// --------------------------
// set 1: for all the muons
// set 2: for the parameter given all the rest cuts were successfull
// set 3&4: for the cuts in sequence they are applied
// Mode = 1: tight ID (default) 
//      = 2: loose ID
//-----------------------------------------------------------------------------
  int     id_word(0xffffffff);
 
  if      (Mode == 1) id_word = IDWord     (Track);
  else if (Mode == 2) id_word = LooseIDWord(Track);
//-----------------------------------------------------------------------------
//  1. number of points
//-----------------------------------------------------------------------------
  Hist.fNActive[0]->Fill(Track->NActive());
  if ((id_word & ~kNActiveBit) == 0) Hist.fNActive[1]->Fill(Track->NActive());
  if (id_word == 0) Hist.fNActive[4]->Fill(Track->NActive());

//    iso1 = Muo->Iso()/pt;
//  
//    Hist.fIso1[0]->Fill(iso1);
//    if ((id_word & ~kIso1Bit) == 0) Hist.fIso1[1]->Fill(iso1);
//    if (id_word == 0) Hist.fIso1[4]->Fill(iso1);
//  //-----------------------------------------------------------------------------
//  //  2. calorimetry
//  //-----------------------------------------------------------------------------
//    Hist.fHadEnergy[0]->Fill(Muo->HadEnergy());
//    if ((id_word & ~kEHadBit) == 0) Hist.fHadEnergy[1]->Fill(Muo->HadEnergy());
//    if (id_word == 0) Hist.fHadEnergy[4]->Fill(Muo->HadEnergy());
//  
//    Hist.fEmEnergy[0]->Fill(Muo->EmEnergy());
//    if ((id_word & ~kEEmBit) == 0) Hist.fEmEnergy[1]->Fill(Muo->EmEnergy());
//    if (id_word == 0) Hist.fEmEnergy[4]->Fill(Muo->EmEnergy());
//  //-----------------------------------------------------------------------------
//  //  7. CMU residuals
//  //-----------------------------------------------------------------------------
//    if (Muo->HasCmuStub()) {
//      Hist.fCmuDelX[0]->Fill(Muo->CmuDelX());
//      if ((id_word & ~kCmuDelXBit) == 0) Hist.fCmuDelX[1]->Fill(Muo->CmuDelX());
//      if (id_word == 0) Hist.fCmuDelX[4]->Fill(Muo->CmuDelX());
//  
//      Hist.fCmuDelZ[0]->Fill(Muo->CmuDelZ());
//      if ((id_word & ~kCmuDelZBit) == 0) Hist.fCmuDelZ[1]->Fill(Muo->CmuDelZ());
//      if (id_word == 0) Hist.fCmuDelZ[4]->Fill(Muo->CmuDelZ());
//    }
//  //-----------------------------------------------------------------------------
//  //  8. CMP residuals
//  //-----------------------------------------------------------------------------
//    if (Muo->HasCmpStub()) {
//      Hist.fCmpDelX[0]->Fill(Muo->CmpDelX());
//      if ((id_word & ~kCmpDelXBit) == 0) Hist.fCmpDelX[1]->Fill(Muo->CmpDelX());
//      if (id_word == 0) Hist.fCmpDelX[4]->Fill(Muo->CmpDelX());
//    }
//  //-----------------------------------------------------------------------------
//  //  9. CMX residuals
//  //-----------------------------------------------------------------------------
//    if (Muo->HasCmxStub()) {
//      Hist.fCmxDelX[0]->Fill(Muo->CmxDelX());
//      if ((id_word & ~kCmxDelXBit) == 0) Hist.fCmxDelX[1]->Fill(Muo->CmxDelX());
//      if (id_word == 0) Hist.fCmxDelX[4]->Fill(Muo->CmxDelX());
//    }
//  //-----------------------------------------------------------------------------
//  // 10. fill track histograms - cut on corrected D0
//  //-----------------------------------------------------------------------------
//    Hist.fTIso[0]->Fill(Muo->TIso());
//    if ((id_word & ~kTIsoBit) == 0) Hist.fTIso[1]->Fill(Muo->TIso());
//    if (id_word == 0) Hist.fTIso[4]->Fill(Muo->TIso());
//  
//    double d0_corr;
//    int    n_ax_hits, n_st_hits;
//  
//    int it = Muo->TrackNumber(); 
//    if (it < 0) {
//  //-----------------------------------------------------------------------------
//  // muon track is not defined
//  //-----------------------------------------------------------------------------
//      d0_corr   = 1.e6;
//      n_ax_hits = -1;
//      n_st_hits = -1;
//    }
//    else {
//  //-----------------------------------------------------------------------------
//  // muon has a track
//  //-----------------------------------------------------------------------------
//      trk       = Muo->Track();
//      d0_corr   = TStntuple::CorrectedD0(trk);
//      n_ax_hits = trk->NCotHitsAx();
//      n_st_hits = trk->NCotHitsSt();
//    }
//  
//    Hist.fTrackD0[0]->Fill(d0_corr);
//    if ((id_word & ~kTrackD0Bit) == 0) Hist.fTrackD0[1]->Fill(d0_corr);
//    if (id_word == 0) Hist.fTrackD0[4]->Fill(d0_corr);
//  
//    Hist.fTrackZ0[0]->Fill(Muo->Z0());
//    if ((id_word & ~kTrackZ0Bit) == 0) Hist.fTrackZ0[1]->Fill(Muo->Z0());
//    if (id_word == 0) Hist.fTrackZ0[4]->Fill(Muo->Z0());
//  
//    Hist.fCotNAxHits[0]->Fill(n_ax_hits);
//    if ((id_word & ~kCotNAxBit) == 0) Hist.fCotNAxHits[1]->Fill(n_ax_hits);
//    if (id_word == 0) Hist.fCotNAxHits[4]->Fill(n_ax_hits);
//  
//    Hist.fCotNStHits[0]->Fill(n_st_hits);
//    if ((id_word & ~kCotNStBit) == 0) Hist.fCotNStHits[1]->Fill(n_st_hits);
//    if (id_word == 0) Hist.fCotNStHits[4]->Fill(n_st_hits);
//  //-----------------------------------------------------------------------------
//  //  single histogram showing how often every particular cut failed
//  //-----------------------------------------------------------------------------
//    for (int bit=0; bit<32; bit++) {
//      if (((id_word >> bit) & 0x1) == 1) {
//        Hist.fFailedBits->Fill(bit);
//      }
//    }
//    Hist.fPassed->Fill(id_word == 0);
//  //-----------------------------------------------------------------------------
//  //  ***** now histogram effect of cuts in the order they are applied
//  //-----------------------------------------------------------------------------
//    Hist.fDetector[2]->Fill(Muo->Detector());
//    Hist.fDetEta  [2]->Fill(dteta);
//    if ((id_word & kDetectorBit) != 0)                        goto END;
//    Hist.fDetector[3]->Fill(Muo->Detector());
//    Hist.fDetEta  [3]->Fill(dteta);
//  
//    Hist.fTrackPt[2]->Fill(pt);
//    if ((id_word & kTrackPtBit) != 0)                         goto END;
//    Hist.fTrackPt[3]->Fill(pt);
//  
//    Hist.fEIso[2]->Fill(Muo->Iso());
//    if ((id_word & kIsoBit) != 0)                             goto END;
//    Hist.fEIso[3]->Fill(Muo->Iso());
//  
//    Hist.fIso1[2]->Fill(Muo->TIso());
//    if ((id_word & kIso1Bit) != 0)                            goto END;
//    Hist.fIso1[3]->Fill(Muo->TIso());
//  
//    Hist.fHadEnergy[2]->Fill(Muo->HadEnergy());
//    if ((id_word & kEHadBit) != 0)                            goto END;
//    Hist.fHadEnergy[3]->Fill(Muo->HadEnergy());
//  
//    Hist.fEmEnergy[2]->Fill(Muo->EmEnergy());
//    if ((id_word & kEEmBit) != 0)                             goto END;
//    Hist.fEmEnergy[3]->Fill(Muo->EmEnergy());
//  
//    if (Muo->HasCmuStub()) {
//      Hist.fCmuDelX[2]->Fill(Muo->CmuDelX());
//      if ((id_word & kCmuDelXBit) != 0)                       goto END;
//      Hist.fCmuDelX[3]->Fill(Muo->CmuDelX());
//  
//      Hist.fCmuDelZ[2]->Fill(Muo->CmuDelZ());
//      if ((id_word & kCmuDelZBit) != 0)                       goto END;
//      Hist.fCmuDelZ[3]->Fill(Muo->CmuDelZ());
//    }
//  
//    Hist.fTrackD0[2]->Fill(d0_corr);
//    if ((id_word & kTrackD0Bit) != 0)                         goto END;
//    Hist.fTrackD0[3]->Fill(d0_corr);
//  
//    Hist.fTrackZ0[2]->Fill(Muo->Z0());
//    if ((id_word & kTrackZ0Bit) != 0)                         goto END;
//    Hist.fTrackZ0[3]->Fill(Muo->Z0());
//  
//    Hist.fCotNAxHits[2]->Fill(n_ax_hits);
//    if ((id_word & kCotNAxBit) != 0)                          goto END;
//    Hist.fCotNAxHits[3]->Fill(n_ax_hits);
//  
//    Hist.fCotNStHits[2]->Fill(n_st_hits);
//    if ((id_word & kCotNStBit) != 0)                          goto END;
//    Hist.fCotNStHits[3]->Fill(n_st_hits);
//  
//    Hist.fTrackBcPt[2]->Fill(bcpt);
//    if ((id_word & kTrackBcPtBit) != 0)                       goto END;
//    Hist.fTrackBcPt[3]->Fill(bcpt);
//  
//   END:;
}


//_____________________________________________________________________________
void TStnTrackID::Print(const char* Opt) const {
  printf("-----------------------------------------------------\n");
  printf("      track ID cuts                                  \n");
  printf("-----------------------------------------------------\n");
  printf(" bit  0: fMinFitCons     = %12.4f\n",fMinFitCons   );
  printf(" bit  1: fMinT0          = %12.4f\n",fMinT0        );
  printf(" bit  2: fMinNActive     = %12.i\n" ,fMinNActive   );
  printf(" bit  3: fMaxT0Err       = %12.4f\n",fMaxT0Err     );
  printf(" bit  4: fMaxFitMomErr   = %12.4f\n",fMaxFitMomErr );
  printf(" bit  5: fTanDip         = %12.4f < tan(dip)   < %12.4f\n",fMinTanDip,fMaxTanDip);
  printf(" bit  6: fD1             = %12.4f < D0         < %12.4f\n",fMinD1    ,fMaxD1    );
  printf(" bit  7: fD2             = %12.4f < D0+2/omega < %12.4f\n",fMinD2    ,fMaxD2    );
}


