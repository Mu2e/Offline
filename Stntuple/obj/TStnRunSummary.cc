///////////////////////////////////////////////////////////////////////////////
//  STNTUPLE run summary table to be stored in the STNTUPLE DB subdirectory
//
//  Author:    P.Murat (CDF/FNAL)
//  Date:      Nov 03 2001
//  2009-05-21: add arrays[5]  of run section ranges
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/obj/TStnRunSummary.hh"
#include <values.h>


ClassImp(TStnRunSummary)

//_____________________________________________________________________________
void TStnRunSummary::ReadV1(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// added in v2
//-----------------------------------------------------------------------------
  fTriggerTableTag = -1;
  fL1Accepts = -1;
  fL1Accepts = -1;
  fL1Accepts = -1;
//-----------------------------------------------------------------------------
// added in v3
//-----------------------------------------------------------------------------
  fBeamStatus    = 0;
  fGoodRunStatus = 0;
  fOfflineStatus = 0;
//-----------------------------------------------------------------------------
// added in v4, all runs are good by default
//-----------------------------------------------------------------------------
  fStatusWord    = 1;
  fOnlineLumiRS  = 0;
  fOfflineLumiRS = 0;
//-----------------------------------------------------------------------------
// added in v5
//-----------------------------------------------------------------------------
  fB0InitLumi    = 0.;
  fB0TermLumi    = 0.;
//-----------------------------------------------------------------------------
// read added in V6 low and high runsections
// assumption that there is only one good runsection range per run had to be
// revisited in V7
// including V8 changes: good runs: 1 RS range per run.... can be more...
//-----------------------------------------------------------------------------
  for (int i=1; i<kMaxRsRanges; i++) {
    fRsRange[i][0]        = 0;
    fRsRange[i][1]        = INT_MAX;
    fGoodOnlineLumiRS [i] = 0;
    fGoodOfflineLumiRS[i] = 0;
  }

  fNRsRanges            = 1;
  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  for (int i=0; i<fNRsRanges; i++) {
    fStatusWordRs[i] = fStatusWord;
    fOfflStatusRs[i] = fOfflineStatus;
  }
}

//_____________________________________________________________________________
void TStnRunSummary::ReadV2(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b >> fTriggerTableTag;
  R__b >> fL1Accepts;
  R__b >> fL2Accepts;
  R__b >> fL3Accepts;
  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// added in v3
//-----------------------------------------------------------------------------
  fBeamStatus    = 0;
  fGoodRunStatus = 0;
  fOfflineStatus = 0;
//-----------------------------------------------------------------------------
// added in v4
//-----------------------------------------------------------------------------
  fStatusWord    = 1;
  fOnlineLumiRS  = 0;
  fOfflineLumiRS = 0;
//-----------------------------------------------------------------------------
// added in v5
//-----------------------------------------------------------------------------
  fB0InitLumi    = 0.;
  fB0TermLumi    = 0.;
//-----------------------------------------------------------------------------
// read added in V6 low and high runsections
// assumption that there is only one good runsection range per run had to be
// revisited in V7
// including V8 changes: good runs: 1 RS range per run.... can be more...
//-----------------------------------------------------------------------------
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]        = 0;
    fRsRange[i][1]        = INT_MAX;
    fGoodOnlineLumiRS [i] = 0;
    fGoodOfflineLumiRS[i] = 0;
  }

  fNRsRanges            = 1;
  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  for (int i=0; i<fNRsRanges; i++) {
    fStatusWordRs[i] = fStatusWord;
    fOfflStatusRs[i] = fOfflineStatus;
  }
}

//_____________________________________________________________________________
void TStnRunSummary::ReadV3(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b >> fTriggerTableTag;
  R__b >> fL1Accepts;
  R__b >> fL2Accepts;
  R__b >> fL3Accepts;
  R__b >> fBeamStatus;
  R__b >> fGoodRunStatus;
  R__b >> fOfflineStatus;
  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// added in v4
//-----------------------------------------------------------------------------
  fStatusWord    = 1;
  fOnlineLumiRS  = 0;
  fOfflineLumiRS = 0;
//-----------------------------------------------------------------------------
// added in v5
//-----------------------------------------------------------------------------
  fB0InitLumi    = 0.;
  fB0TermLumi    = 0.;
//-----------------------------------------------------------------------------
// read added in V6 low and high runsections
// assumption that there is only one good runsection range per run had to be
// revisited in V7
// including V8 changes: good runs: 1 RS range per run.... can be more...
//-----------------------------------------------------------------------------
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]        = 0;
    fRsRange[i][1]        = INT_MAX;
    fGoodOnlineLumiRS [i] = 0;
    fGoodOfflineLumiRS[i] = 0;
  }

  fNRsRanges            = 1;
  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  for (int i=0; i<fNRsRanges; i++) {
    fStatusWordRs[i] = fStatusWord;
    fOfflStatusRs[i] = fOfflineStatus;
  }
}

//_____________________________________________________________________________
void TStnRunSummary::ReadV4(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fOnlineLumiRS;
  R__b >> fOfflineLumiRS;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b >> fTriggerTableTag;
  R__b >> fL1Accepts;
  R__b >> fL2Accepts;
  R__b >> fL3Accepts;
  R__b >> fBeamStatus;
  R__b >> fGoodRunStatus;
  R__b >> fOfflineStatus;
  R__b >> fStatusWord;
  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// added in v4
//-----------------------------------------------------------------------------
  fStatusWord    = 1;
  fOnlineLumiRS  = 0;
  fOfflineLumiRS = 0;
//-----------------------------------------------------------------------------
// added in v5
//-----------------------------------------------------------------------------
  fB0InitLumi    = 0.;
  fB0TermLumi    = 0.;
//-----------------------------------------------------------------------------
// read added in V6 low and high runsections
// assumption that there is only one good runsection range per run had to be
// revisited in V7
// including V8 changes: good runs: 1 RS range per run.... can be more...
//-----------------------------------------------------------------------------
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]        = 0;
    fRsRange[i][1]        = INT_MAX;
    fGoodOnlineLumiRS [i] = 0;
    fGoodOfflineLumiRS[i] = 0;
  }

  fNRsRanges            = 1;
  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  for (int i=0; i<fNRsRanges; i++) {
    fStatusWordRs[i] = fStatusWord;
    fOfflStatusRs[i] = fOfflineStatus;
  }
}

//_____________________________________________________________________________
void TStnRunSummary::ReadV5(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fOnlineLumiRS;
  R__b >> fOfflineLumiRS;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b >> fTriggerTableTag;
  R__b >> fL1Accepts;
  R__b >> fL2Accepts;
  R__b >> fL3Accepts;
  R__b >> fBeamStatus;
  R__b >> fGoodRunStatus;
  R__b >> fOfflineStatus;
  R__b >> fStatusWord;
  R__b >> fB0InitLumi;
  R__b >> fB0TermLumi;
  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// read added in V6 low and high runsections
// assumption that there is only one good runsection range per run had to be
// revisited in V7
// including V8 changes: good runs: 1 RS range per run.... can be more...
//-----------------------------------------------------------------------------
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]        = 0;
    fRsRange[i][1]        = INT_MAX;
    fGoodOnlineLumiRS [i] = 0;
    fGoodOfflineLumiRS[i] = 0;
  }

  fNRsRanges            = 1;
  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  for (int i=0; i<fNRsRanges; i++) {
    fStatusWordRs[i]   = fStatusWord;
    fOfflStatusRs[i]   = fOfflineStatus;
  }
}

//_____________________________________________________________________________
void TStnRunSummary::ReadV6(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
//-----------------------------------------------------------------------------
// V6 version of TStnRunSummary
//-----------------------------------------------------------------------------
//  class TStnRunSummary: public TObject {
//  public:
//  protected:
//    TString fObjName;                  // name of this for list searching
//  
//  protected:
//    Int_t     fRunNumber;
//    Int_t     fMyronMode;
//    Int_t     fL1Early;
//    Float_t   fLumiTev;
//    Float_t   fLumiTape;
//    Float_t   fOnlineLumiRS;
//    Float_t   fOfflineLumiRS;
//    Float_t   fLivetime;
//    TString   fDate;
//    TString   fTime;
//    TString   fTriggerTableName;
//    Int_t     fTriggerTableTag;
//    Int_t     fL1Accepts;
//    Int_t     fL2Accepts;
//    Int_t     fL3Accepts;
//    Int_t     fBeamStatus;
//    Int_t     fGoodRunStatus;
//    Int_t     fOfflineStatus;
//    Int_t     fStatusWord;
//    Float_t   fB0InitLumi;	// initial inst lumi (in 1e30), added in v5
//    Float_t   fB0TermLumi;        // added in v5
//    Int_t     fLowRunSection;     // added in V6
//    Int_t     fHighRunSection;    // added in V6
//-----------------------------------------------------------------------------
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fOnlineLumiRS;
  R__b >> fOfflineLumiRS;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b >> fTriggerTableTag;
  R__b >> fL1Accepts;
  R__b >> fL2Accepts;
  R__b >> fL3Accepts;
  R__b >> fBeamStatus;
  R__b >> fGoodRunStatus;
  R__b >> fOfflineStatus;
  R__b >> fStatusWord;
  R__b >> fB0InitLumi;
  R__b >> fB0TermLumi;
//-----------------------------------------------------------------------------
// read added in V6 low and high runsections
// assumption that there is only one good runsection range per run had to be
// revisited in V7
//-----------------------------------------------------------------------------
  R__b >> fRsRange[0][0];
  R__b >> fRsRange[0][1];
  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// added in v7
// including V8 changes: good runs: 1 RS range per run.... can be more...
//-----------------------------------------------------------------------------
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]        = 0;
    fRsRange[i][1]        = INT_MAX;
    fGoodOnlineLumiRS [i] = 0;
    fGoodOfflineLumiRS[i] = 0;
  }

  fNRsRanges            = 1;
  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  for (int i=0; i<kMaxRsRanges; i++) {
    fStatusWordRs[i] = fStatusWord;
    fOfflStatusRs[i] = fOfflineStatus;
  }
}

//_____________________________________________________________________________
void TStnRunSummary::ReadV7(TBuffer& R__b, UInt_t R__s, UInt_t R__c) {
//-----------------------------------------------------------------------------
// V6 version of TStnRunSummary
//-----------------------------------------------------------------------------
//  class TStnRunSummary: public TObject {
//  public:
//  protected:
//    TString fObjName;                  // name of this for list searching
//  
//  protected:
//    Int_t     fRunNumber;
//    Int_t     fMyronMode;
//    Int_t     fL1Early;
//    Float_t   fLumiTev;
//    Float_t   fLumiTape;
//    Float_t   fOnlineLumiRS;
//    Float_t   fOfflineLumiRS;
//    Float_t   fLivetime;
//    TString   fDate;
//    TString   fTime;
//    TString   fTriggerTableName;
//    Int_t     fTriggerTableTag;
//    Int_t     fL1Accepts;
//    Int_t     fL2Accepts;
//    Int_t     fL3Accepts;
//    Int_t     fBeamStatus;
//    Int_t     fGoodRunStatus;
//    Int_t     fOfflineStatus;
//    Int_t     fStatusWord;
//    Float_t   fB0InitLumi;	// initial inst lumi (in 1e30), added in v5
//    Float_t   fB0TermLumi;        // added in v5
//    Int_t     fLowRunSection;     // added in V6
//    Int_t     fHighRunSection;    // added in V6
//-----------------------------------------------------------------------------
  TObject::Streamer(R__b);
  R__b >> fRunNumber;
  R__b >> fMyronMode;
  R__b >> fL1Early;
  R__b >> fLumiTev;
  R__b >> fLumiTape;
  R__b >> fOnlineLumiRS;
  R__b >> fOfflineLumiRS;
  R__b >> fLivetime;
  fDate.Streamer(R__b);
  fTime.Streamer(R__b);
  fTriggerTableName.Streamer(R__b);
  R__b >> fTriggerTableTag;
  R__b >> fL1Accepts;
  R__b >> fL2Accepts;
  R__b >> fL3Accepts;
  R__b >> fBeamStatus;
  R__b >> fGoodRunStatus;
  R__b >> fOfflineStatus;
  R__b >> fStatusWord;
  R__b >> fB0InitLumi;
  R__b >> fB0TermLumi;
  R__b >> fNRsRanges;
  for (int i=0; i<fNRsRanges; i++) {
    R__b >> fRsRange[i][0];
    R__b >> fRsRange[i][1];
  }

  fGoodOnlineLumiRS [0] = fOnlineLumiRS ;
  fGoodOfflineLumiRS[0] = fOfflineLumiRS;

  R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
//-----------------------------------------------------------------------------
// v7: status words are the same for all the RS ranges
// including V8 changes: good runs: 1 RS range per run.... can be more...
// kludge: if fNsRunRanges == 0, fStatusWordRs[0] still needs to be defined
//-----------------------------------------------------------------------------
  fStatusWordRs[0] = fStatusWord;
  fOfflStatusRs[0] = fOfflineStatus;

  for (int i=0; i<fNRsRanges; i++) {
    fStatusWordRs[i] = fStatusWord;
    fOfflStatusRs[i] = fOfflineStatus;
  }

}


//______________________________________________________________________________
void TStnRunSummary::Streamer(TBuffer &R__b) {
   // Stream an object of class TStnRunSummary.

  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
//-----------------------------------------------------------------------------
// read branch, check version
//-----------------------------------------------------------------------------
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
    if      (R__v == 1) ReadV1(R__b,R__s,R__c);
    else if (R__v == 2) ReadV2(R__b,R__s,R__c);
    else if (R__v == 3) ReadV3(R__b,R__s,R__c);
    else if (R__v == 4) ReadV4(R__b,R__s,R__c);
    else if (R__v == 5) ReadV5(R__b,R__s,R__c);
    else if (R__v == 6) ReadV6(R__b,R__s,R__c);
    else if (R__v == 7) ReadV7(R__b,R__s,R__c);
    else {
//-----------------------------------------------------------------------------
// V6: 2006-09-19
// V7: 2009-05-21 
// V8: 2010-03-21 (current version): different RS ranges can have different
//                                   good run bit settings
//-----------------------------------------------------------------------------
      TObject::Streamer(R__b);
      R__b >> fRunNumber;
      R__b >> fMyronMode;
      R__b >> fL1Early;
      R__b >> fLumiTev;
      R__b >> fLumiTape;
      R__b >> fOnlineLumiRS;
      R__b >> fOfflineLumiRS;
      R__b >> fLivetime;
      fDate.Streamer(R__b);
      fTime.Streamer(R__b);
      fTriggerTableName.Streamer(R__b);
      R__b >> fTriggerTableTag;
      R__b >> fL1Accepts;
      R__b >> fL2Accepts;
      R__b >> fL3Accepts;
      R__b >> fBeamStatus;
      R__b >> fGoodRunStatus;
      R__b >> fOfflineStatus;
      R__b >> fStatusWord;
      R__b >> fB0InitLumi;
      R__b >> fB0TermLumi;
      R__b >> fNRsRanges;

      R__b.ReadFastArray(&fRsRange[0][0],2*fNRsRanges);

				        // added in V8

      R__b.ReadFastArray(fGoodOnlineLumiRS ,fNRsRanges);
      R__b.ReadFastArray(fGoodOfflineLumiRS,fNRsRanges);
      R__b.ReadFastArray(fStatusWordRs     ,fNRsRanges);
      R__b.ReadFastArray(fOfflStatusRs     ,fNRsRanges);

      R__b.CheckByteCount(R__s, R__c, TStnRunSummary::IsA());
    }
  } 
  else {
//-----------------------------------------------------------------------------
// write branch, current version: V8
//-----------------------------------------------------------------------------
    R__c = R__b.WriteVersion(TStnRunSummary::IsA(), kTRUE);
    TObject::Streamer(R__b);
    R__b << fRunNumber;
    R__b << fMyronMode;
    R__b << fL1Early;
    R__b << fLumiTev;
    R__b << fLumiTape;
    R__b << fOnlineLumiRS;
    R__b << fOfflineLumiRS;
    R__b << fLivetime;
    fDate.Streamer(R__b);
    fTime.Streamer(R__b);
    fTriggerTableName.Streamer(R__b);
    R__b << fTriggerTableTag;
    R__b << fL1Accepts;
    R__b << fL2Accepts;
    R__b << fL3Accepts;
    R__b << fBeamStatus;
    R__b << fGoodRunStatus;
    R__b << fOfflineStatus;
    R__b << fStatusWord;
    R__b << fB0InitLumi;
    R__b << fB0TermLumi;
//-----------------------------------------------------------------------------
// write out good run section ranges
//-----------------------------------------------------------------------------
    R__b << fNRsRanges;
    R__b.WriteFastArray(&fRsRange[0][0],2*fNRsRanges);

				        // added in V8

    R__b.WriteFastArray(fGoodOnlineLumiRS ,fNRsRanges);
    R__b.WriteFastArray(fGoodOfflineLumiRS,fNRsRanges);
    R__b.WriteFastArray(fStatusWordRs     ,fNRsRanges);
    R__b.WriteFastArray(fOfflStatusRs     ,fNRsRanges);

    R__b.SetByteCount(R__c, kTRUE);
  }
}

//_____________________________________________________________________________
TStnRunSummary::TStnRunSummary():
  fObjName("TStnRunSummary"),
  fDate(""),
  fTime(""),
  fTriggerTableName("")
{
  fRunNumber        = -1;
  fMyronMode        = 0;
  fL1Early          = 0;
  fLumiTev          = -1;
  fLumiTape         = -1;
  fOnlineLumiRS     = -1;
  fOfflineLumiRS    = -1;
  fLivetime         = -1;
  fTriggerTableTag  = -1;
  fL1Accepts        = 0;
  fL2Accepts        = 0;
  fL3Accepts        = 0;
  fBeamStatus       = 0;
  fGoodRunStatus    = 0;
  fOfflineStatus    = 0;
  fStatusWord       = 0;
  fB0InitLumi       = -1;
  fB0TermLumi       = -1;
//-----------------------------------------------------------------------------
// added in v7, replacing V6
//-----------------------------------------------------------------------------
  fNRsRanges        = -1;
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]  = -1;
    fRsRange[i][1]  = -1;
    fGoodOnlineLumiRS[i]  = -1;
    fGoodOfflineLumiRS[i] = -1;
  }
}

//_____________________________________________________________________________
TStnRunSummary::TStnRunSummary(const TStnRunSummary& Rs)
                                        :fObjName("RunSummary") {
  fRunNumber        = Rs.fRunNumber;
  fMyronMode        = Rs.fMyronMode;
  fL1Early          = Rs.fL1Early;
  fLumiTev          = Rs.fLumiTev;
  fLumiTape         = Rs.fLumiTape;
  fOnlineLumiRS     = Rs.fOnlineLumiRS;
  fOfflineLumiRS    = Rs.fOfflineLumiRS;
  fLivetime         = Rs.fLivetime;
  fDate             = Rs.fDate;
  fTime             = Rs.fTime;
  fTriggerTableName = Rs.fTriggerTableName;
  fTriggerTableTag  = Rs.fTriggerTableTag;
  fL1Accepts        = Rs.fL1Accepts;
  fL2Accepts        = Rs.fL2Accepts;
  fL3Accepts        = Rs.fL3Accepts;
  fBeamStatus       = Rs.fBeamStatus;
  fGoodRunStatus    = Rs.fGoodRunStatus;
  fOfflineStatus    = Rs.fOfflineStatus;
  fStatusWord       = Rs.fStatusWord;
  fB0InitLumi       = -1;
  fB0TermLumi       = -1;
//-----------------------------------------------------------------------------
// added in v7, replacing V6
//-----------------------------------------------------------------------------
  fNRsRanges        = -1;
  for (int i=0; i<kMaxRsRanges; i++) {
    fRsRange[i][0]  = -1;
    fRsRange[i][1]  = -1;
    fGoodOnlineLumiRS [i]  = -1;
    fGoodOfflineLumiRS[i] = -1;
  }
}

//_____________________________________________________________________________
TStnRunSummary::~TStnRunSummary() {
}

//_____________________________________________________________________________
void TStnRunSummary::Clear(Option_t* Opt) {
}

//_____________________________________________________________________________
void TStnRunSummary::Print(Option_t* Opt) const {

  int         n;
  char        buffer[200];
  char*       tail;

  if (strstr(Opt,"banner") || (Opt[0] == 0)) {
    printf("------------------------------------------------------------------------");
    printf("----------------------------------");
    printf("----------------------------------------------------------------------\n");
    printf("   Run [rs_lo:rs_hi]:G.M.L   Trigger Table     ");
    printf("    B0InitL    L(TeV)   L(tape)  Lonl(RS)  Lofl(RS)   Lonl(RS)  Lofl(RS)");
    printf(" Livetime     Date      Time        L1A       L2A     L3A  \n");
    printf("                     R.                                                    ");
    printf("                            (good)    (good)\n");
    printf("------------------------------------------------------------------------");
    printf("----------------------------------");
    printf("----------------------------------------------------------------------\n");
  }

  tail = (char*) fTriggerTableName.Data();

  int x;

  if (strstr(Opt,"data") || (Opt[0] == 0)) {
    n = strlen(tail);
    if (n > 22) n = 22;
    strncpy(buffer,tail,n);
    buffer[n] = 0;
    tail = tail+n;

    
    if (fRsRange[0][1] == INT_MAX) x = 99999;
    else                           x = fRsRange[0][1];

    printf("%7i[%5i:%5i]:%1i.%1i.%1i  %-22s %7.1f %9.3f %9.3f %9.3f %9.3f ",
	   fRunNumber,fRsRange[0][0],x,
	   fStatusWordRs[0]&0x1,
	   fMyronMode,
	   fL1Early,
	   buffer,
	   fB0InitLumi,
	   fLumiTev,
	   fLumiTape,
	   fOnlineLumiRS,
	   fOfflineLumiRS);

    printf(" %9.3f %9.3f ",
	   fGoodOnlineLumiRS [0],
	   fGoodOfflineLumiRS[0]);

    printf("%7.1f   %10s %8s %10i %8i %7i\n",
	   fLivetime,
	   fDate.Data(),
	   fTime.Data(),
	   fL1Accepts,
	   fL2Accepts,
	   fL3Accepts);

//-----------------------------------------------------------------------------
// print end of the trigger table name
//-----------------------------------------------------------------------------
    while (*tail) {
      n = strlen(tail);
      if (n > 22) n = 22;
      strncpy(buffer,tail,n);
      buffer[n] = 0;
      tail = tail+n;
      printf("                %-22s\n",buffer);
    }
//-----------------------------------------------------------------------------
// print end of the run section ranges, teh first one has already been printed
//-----------------------------------------------------------------------------
    for (int i=1; i<fNRsRanges; i++) {

      if (fRsRange[i][1] == INT_MAX) x = 99999;
      else                           x = fRsRange[0][1];

      printf("%7s[%5i:%5i] %78s"," ",fRsRange[i][0],x," ");
      printf(" %9.3f %9.3f\n",
	     fGoodOnlineLumiRS [i],
	     fGoodOfflineLumiRS[i]);
    }
  }
}

