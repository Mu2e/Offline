#ifndef ILCDCHRECOPARAM_H
#define ILCDCHRECOPARAM_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with DCH reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

class IlcDCHRecoParam : public TObject
{
 public: 
  IlcDCHRecoParam();
  virtual ~IlcDCHRecoParam();
  Double_t GetCtgRange() const     { return fCtgRange;}
  Double_t GetMaxSnpTracker() const{ return fMaxSnpTracker;}
  Double_t GetMaxSnpTrack() const  { return fMaxSnpTrack;}
  //
  Int_t    GetFirstBin() const     { return fFirstBin;}
  Int_t    GetLastBin() const      { return fLastBin;}
  void     SetTimeBinRange(Int_t first, Int_t last){ fFirstBin = first; fLastBin = last;}
  Bool_t   GetCalcPedestal() const { return fBCalcPedestal;}
  Bool_t   GetDoUnfold() const     { return fBDoUnfold;}
  Float_t  GetDumpAmplitudeMin() const  { return fDumpAmplitudeMin;}
  Float_t  GetMaxNoise() const     { return fMaxNoise;}  
  //
  Bool_t   GetDoKinks() const      { return fBKinkFinder;}
  Float_t  GetMaxC()    const      { return fMaxC;}
  Bool_t   GetSpecialSeeding() const { return fBSpecialSeeding;}
  Bool_t   GetBYMirror() const { return fBYMirror;}
  static   IlcDCHRecoParam *GetLowFluxParam();        // make reco parameters for low  flux env.
  static   IlcDCHRecoParam *GetHighFluxParam();       // make reco parameters for high flux env. 
  static   IlcDCHRecoParam *GetLaserTestParam(Bool_t bPedestal);  // special setting for laser 
  static   IlcDCHRecoParam *GetCosmicTestParam(Bool_t bPedestal); // special setting for cosmic  
  //
 protected:
  Double_t fCtgRange;        // +-fCtgRange is the ctg(Theta) window used for clusterization and tracking (MI) 
  Double_t fMaxSnpTracker;   // max sin of local angle  - for DCH tracker
  Double_t fMaxSnpTrack;     // max sin of local angle  - for track 
  Bool_t   fBYMirror;        // mirror of the y - pad coordinate 
  //
  //   clusterer parameters
  //
  Int_t    fFirstBin;        // first time bin used by cluster finder
  Int_t    fLastBin;         // last time bin  used by cluster finder 
  Bool_t   fBCalcPedestal;   // calculate Pedestal
  Bool_t   fBDoUnfold;       // do unfolding of clusters
  Float_t  fDumpAmplitudeMin; // minimal amplitude of signal to be dumped 
  Float_t  fMaxNoise;        // maximal noise sigma on pad to be used in cluster finder
  //
  //
  Float_t  fMaxC;            // maximal curvature for tracking
  Bool_t   fBSpecialSeeding; // special seeding with big inclination angles allowed (for Cosmic and laser)
  Bool_t   fBKinkFinder;     // do kink finder reconstruction
  ClassDef(IlcDCHRecoParam, 1)
};


#endif
