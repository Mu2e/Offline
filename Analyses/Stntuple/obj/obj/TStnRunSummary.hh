#ifndef STNTUPLE_TStnRunSummary
#define STNTUPLE_TStnRunSummary
//-----------------------------------------------------------------------------
// description of the run summary information to be stored in DB section
// of the STNTUPLE
// 
// Author:    P. Murat (CDF/FNAL)
// Date:      Nov 3 2001
//-----------------------------------------------------------------------------

#include "TObject.h"
#include "TString.h"

class TStnRunSummary: public TObject {
public:
  enum {
    kGoodRunBit       = 0x1 <<  0,
    kCheckedRunBit    = 0x1 <<  2,
    kBeamStatusBit    = 0x1 <<  3,
    kOfflineStatusBit = 0x1 <<  4,
    kClcStatusBit     = 0x1 <<  5,
    kCalStatusBit     = 0x1 <<  6,
    kSmxStatusBit     = 0x1 <<  7,
    kCotStatusBit     = 0x1 <<  8,
    kSvxStatusBit     = 0x1 <<  9,
    kIslStatusBit     = 0x1 << 10,
    kL00StatusBit     = 0x1 << 11,
    kCmuStatusBit     = 0x1 << 12,
    kCmpStatusBit     = 0x1 << 13,
    kCmxStatusBit     = 0x1 << 14,
    kImuStatusBit     = 0x1 << 15,
    kTofStatusBit     = 0x1 << 16,
    kL1tStatusBit     = 0x1 << 17,
    kL2tStatusBit     = 0x1 << 18,
    kL3tStatusBit     = 0x1 << 19,
    kSvtStatusBit     = 0x1 << 20,
    kMplStatusBit     = 0x1 << 21,
    kBscStatusBit     = 0x1 << 22,
    kCesStatusBit     = 0x1 << 23,
    kPesStatusBit     = 0x1 << 24
  };

  enum { kAllStatusBits = 0x01ffffff };

  enum { kMaxRsRanges   = 5          };
protected:
  TString fObjName;                  // name of this for list searching

protected:
  Int_t     fRunNumber;
  Int_t     fMyronMode;
  Int_t     fL1Early;
  Float_t   fLumiTev;
  Float_t   fLumiTape;
  Float_t   fOnlineLumiRS;
  Float_t   fOfflineLumiRS;
  Float_t   fLivetime;
  TString   fDate;
  TString   fTime;
  TString   fTriggerTableName;
  Int_t     fTriggerTableTag;
  Int_t     fL1Accepts;
  Int_t     fL2Accepts;
  Int_t     fL3Accepts;
  Int_t     fBeamStatus;
  Int_t     fGoodRunStatus;
  Int_t     fOfflineStatus;                      // 2010-03-25: not used since V8
  Int_t     fStatusWord;	                 // 2010-03-25: not used since V8
  Float_t   fB0InitLumi;	                 // initial inst lumi (in 1e30), added in v5
  Float_t   fB0TermLumi;                         // added in v5
  Int_t     fNRsRanges;                          // V7
  Int_t     fRsRange          [kMaxRsRanges][2]; // V7: up to 5 rs ranges [low,high]
  Float_t   fGoodOnlineLumiRS [kMaxRsRanges];    // V7
  Float_t   fGoodOfflineLumiRS[kMaxRsRanges];    // V7
  Int_t     fStatusWordRs     [kMaxRsRanges];    // V8: each RS range can have its status
  Int_t     fOfflStatusRs     [kMaxRsRanges];    // V8: each RS range can have its status
//-----------------------------------------------------------------------------
//  methods
//-----------------------------------------------------------------------------
public:
  TStnRunSummary();
  TStnRunSummary(const TStnRunSummary& Rs);
  ~TStnRunSummary();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t    RunNumber        () { return fRunNumber         ; }
  Int_t    MyronMode        () { return fMyronMode         ; }
  Int_t    L1Early          () { return fL1Early           ; }
  Float_t  LumiTev          () { return fLumiTev           ; }
  Float_t  LumiTape         () { return fLumiTape          ; }
  Int_t    L1Accepts        () { return fL1Accepts         ; }
  Int_t    L2Accepts        () { return fL2Accepts         ; }
  Int_t    L3Accepts        () { return fL3Accepts         ; }
  Float_t  Livetime         () { return fLivetime          ; }
  Int_t    BeamStatusWord   () { return fBeamStatus        ; }
  Int_t    GoodRunStatus    () { return fGoodRunStatus     ; }
  Float_t  OnlineLumiRS     () { return fOnlineLumiRS      ; }
  Float_t  OfflineLumiRS    () { return fOfflineLumiRS     ; }
  Float_t  B0InitLumi       () { return fB0InitLumi        ; }
  Float_t  B0TermLumi       () { return fB0TermLumi        ; }
  TString& TriggerTableName () { return fTriggerTableName  ; }
  Int_t    TriggerTableTag  () { return fTriggerTableTag   ; }
  TString& Date             () { return fDate              ; }
  TString& Time             () { return fTime              ; }

  Int_t    NRsRanges        () { return fNRsRanges         ; }
  Int_t    LowRS     (Int_t I) { return fRsRange[I][0]     ; }
  Int_t    HighRS    (Int_t I) { return fRsRange[I][1]     ; }
  Float_t  GoodOnlineLumiRS (int I) { return fGoodOnlineLumiRS [I]; }
  Float_t  GoodOfflineLumiRS(int I) { return fGoodOfflineLumiRS[I]; }
//-----------------------------------------------------------------------------
// various general, detector, trigger and offline status bits (0=BAD, 1=OK)
// list is not complete...
//-----------------------------------------------------------------------------
  Int_t   StatusWord       (Int_t I=0) { return   fStatusWordRs[I]           ; }
  Int_t   OfflineStatusWord(Int_t I=0) { return   fOfflStatusRs[I]           ; }
  Int_t   IsGood           (Int_t I=0) { return ((fStatusWordRs[I]      ) & 0x1); }
  Int_t   IsChecked        (Int_t I=0) { return ((fStatusWordRs[I] >>  2) & 0x1); }
  Int_t   BeamStatusBit    (Int_t I=0) { return ((fStatusWordRs[I] >>  3) & 0x1); }
  Int_t   OfflineStatusBit (Int_t I=0) { return ((fStatusWordRs[I] >>  4) & 0x1); }
  Int_t   ClcStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >>  5) & 0x1); }
  Int_t   CalStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >>  6) & 0x1); }
  Int_t   SmxStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >>  7) & 0x1); }
  Int_t   CotStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >>  8) & 0x1); }
  Int_t   SvxStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >>  9) & 0x1); }
  Int_t   IslStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 10) & 0x1); }
  Int_t   L00StatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 11) & 0x1); }
  Int_t   CmuStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 12) & 0x1); }
  Int_t   CmpStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 13) & 0x1); }
  Int_t   CmxStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 14) & 0x1); }
  Int_t   ImuStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 15) & 0x1); }
  Int_t   TofStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 16) & 0x1); }
  Int_t   L1tStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 17) & 0x1); }
  Int_t   L2tStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 18) & 0x1); }
  Int_t   L3tStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 19) & 0x1); }
  Int_t   SvtStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 20) & 0x1); }
  Int_t   MplStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 21) & 0x1); }
  Int_t   BscStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 22) & 0x1); }
  Int_t   CesStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 23) & 0x1); }
  Int_t   PesStatusBit     (Int_t I=0) { return ((fStatusWordRs[I] >> 24) & 0x1); }
//-----------------------------------------------------------------------------
// offline status bits
//-----------------------------------------------------------------------------
  Int_t   ClcOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >>  5) & 0x1); }
  Int_t   CalOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >>  6) & 0x1); }
  Int_t   SmxOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >>  7) & 0x1); }
  Int_t   CotOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >>  8) & 0x1); }
  Int_t   SvxOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >>  9) & 0x1); }
  Int_t   IslOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 10) & 0x1); }
  Int_t   L00OfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 11) & 0x1); }
  Int_t   CmuOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 12) & 0x1); }
  Int_t   CmpOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 13) & 0x1); }
  Int_t   CmxOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 14) & 0x1); }
  Int_t   ImuOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 15) & 0x1); }
  Int_t   TofOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 16) & 0x1); }
  Int_t   L1tOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 17) & 0x1); }
  Int_t   L2tOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 18) & 0x1); }
  Int_t   L3tOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 19) & 0x1); }
  Int_t   SvtOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 20) & 0x1); }
  Int_t   MplOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 21) & 0x1); }
  Int_t   BscOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 22) & 0x1); }
  Int_t   CesOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 23) & 0x1); }
  Int_t   PesOfflineBit    (Int_t I=0) { return ((fOfflStatusRs[I] >> 24) & 0x1); }
//-----------------------------------------------------------------------------
// average instantaneous luminosity: integrated luminosity / livetime
// number of buckets in the Tevatron ring: 159, number of filled ones: 36
// so in reality average bunch spacing 132ns*159/36 = 583ns is 47% larger 
// than "conventional" number of 396ns
// fLivetime is already in seconds, converted from bunchcrossings
//-----------------------------------------------------------------------------
  Float_t AverageInstLumi() {
    return 1.0e33*fOfflineLumiRS/fLivetime;
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    SetRunNumber        (Int_t       Run ) { fRunNumber         = Run;  }
  void    SetMyronMode        (Int_t       Mode) { fMyronMode         = Mode; }
  void    SetL1Early          (Int_t       L1  ) { fL1Early           = L1;   }
  void    SetLumiTev          (Float_t     Lumi) { fLumiTev           = Lumi; }
  void    SetLumiTape         (Float_t     Lumi) { fLumiTape          = Lumi; }
  void    SetB0InitLumi       (Float_t     Lumi) { fB0InitLumi        = Lumi; }
  void    SetB0TermLumi       (Float_t     Lumi) { fB0TermLumi        = Lumi; }
  void    SetLivetime         (Float_t     Time) { fLivetime          = Time; }
  void    SetDate             (const char* Date) { fDate              = Date; }
  void    SetTime             (const char* Time) { fTime              = Time; }
  void    SetL1Accepts        (Int_t       N   ) { fL1Accepts         = N   ; }
  void    SetL2Accepts        (Int_t       N   ) { fL2Accepts         = N   ; }
  void    SetL3Accepts        (Int_t       N   ) { fL3Accepts         = N   ; }
  void    SetBeamStatus       (Int_t       Word) { fBeamStatus        = Word; }
  void    SetNRsRanges        (Int_t       N   ) { fNRsRanges         = N   ; }
  void    SetGoodRunStatus    (Int_t       Word) { fGoodRunStatus     = Word; }
  void    SetOnlineLumiRS     (Float_t     Lumi) { fOnlineLumiRS      = Lumi; }
  void    SetOfflineLumiRS    (Float_t     Lumi) { fOfflineLumiRS     = Lumi; }
  void    SetTriggerTableName (const char* Name) { fTriggerTableName  = Name; }
  void    SetTriggerTableTag  (Int_t        Tag) { fTriggerTableTag   = Tag;  }

  void    SetGoodOnlineLumiRS (Float_t Lumi, int I=0) { fGoodOnlineLumiRS [I] = Lumi; }
  void    SetGoodOfflineLumiRS(Float_t Lumi, int I=0) { fGoodOfflineLumiRS[I] = Lumi; }
  void    SetStatusWord       (int     Word, int I=0) { fStatusWordRs[I] = Word; }
  void    SetOfflineStatusWord(int     Word, int I=0) { fOfflStatusRs[I] = Word; }

  void    SetRsRange          (Int_t I, Int_t LowRs, Int_t HighRs) { 
    fRsRange[I][0] = LowRs; fRsRange[I][1] = HighRs; 
  }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  const char* GetName() const { return fObjName.Data(); }
  void        Clear(Option_t* option = "");
  void        Print(Option_t* option = "") const;   // *MENU* 
//-----------------------------------------------------------------------------
// schema evolution
// v6: adds low and high runsections
// v7: replaces low and high runsections with arrays [5]
//-----------------------------------------------------------------------------
  void      ReadV1(TBuffer& R__b, UInt_t R__s, UInt_t R__c);
  void      ReadV2(TBuffer& R__b, UInt_t R__s, UInt_t R__c);
  void      ReadV3(TBuffer& R__b, UInt_t R__s, UInt_t R__c);
  void      ReadV4(TBuffer& R__b, UInt_t R__s, UInt_t R__c);
  void      ReadV5(TBuffer& R__b, UInt_t R__s, UInt_t R__c);
  void      ReadV6(TBuffer& R__b, UInt_t R__s, UInt_t R__c);
  void      ReadV7(TBuffer& R__b, UInt_t R__s, UInt_t R__c);

  ClassDef(TStnRunSummary,8)
};

#endif

