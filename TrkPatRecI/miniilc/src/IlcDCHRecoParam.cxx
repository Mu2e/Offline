/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with DCH reconstruction parameters                                  //
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "IlcDCHRecoParam.h"

//ClassImp(IlcDCHRecoParam)




//_____________________________________________________________________________
IlcDCHRecoParam::IlcDCHRecoParam():
  fCtgRange(100.05),       
  fMaxSnpTracker(0.999),
  fMaxSnpTrack(0.999),
  fBYMirror(kTRUE),
  fFirstBin(0),
  fLastBin(-1),
  fBCalcPedestal(kFALSE),
  fBDoUnfold(kTRUE),
  fDumpAmplitudeMin(100),
  fMaxNoise(3.),
  fMaxC(0.3),
  fBSpecialSeeding(kFALSE),
  fBKinkFinder(kTRUE)
{
  //
  // constructor
  //
}

//_____________________________________________________________________________
IlcDCHRecoParam::~IlcDCHRecoParam() 
{
  //
  // destructor
  //  
}




IlcDCHRecoParam *IlcDCHRecoParam::GetLowFluxParam(){
  //
  // make default reconstruction  parameters for low  flux env.
  //
  IlcDCHRecoParam *param = new IlcDCHRecoParam;
  param->fCtgRange = 10;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  return param;
}

IlcDCHRecoParam *IlcDCHRecoParam::GetHighFluxParam(){
  //
  // make reco parameters for high flux env.
  //
  IlcDCHRecoParam *param = new IlcDCHRecoParam;
  param->fCtgRange = 10.05;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  return param;
}

IlcDCHRecoParam *IlcDCHRecoParam::GetLaserTestParam(Bool_t bPedestal){
  //
  // special setting for laser
  //
  IlcDCHRecoParam *param = new IlcDCHRecoParam;
  param->fCtgRange = 10.05;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  param->fBCalcPedestal = bPedestal;
  param->fBDoUnfold     = kFALSE;
  param->fDumpAmplitudeMin = 150;
  param->fBKinkFinder   = kFALSE;
  param->fMaxSnpTracker = 0.98;
  param->fMaxC          = 0.02;
  param->fBSpecialSeeding = kTRUE;
  param->fBYMirror      = kFALSE;
  return param;
}

IlcDCHRecoParam *IlcDCHRecoParam::GetCosmicTestParam(Bool_t bPedestal){
  //
  // special setting for cosmic 
  // 
  IlcDCHRecoParam *param = new IlcDCHRecoParam;
  param->fCtgRange = 10.05;    // full DCH
  param->fFirstBin = 60;
  param->fLastBin  = 1000;
  param->fBCalcPedestal = bPedestal;
  param->fBDoUnfold     = kFALSE;
  param->fBSpecialSeeding = kTRUE;
  param->fMaxC          = 0.07;
  param->fBKinkFinder   = kFALSE;
  param->fBYMirror      = kFALSE;
  return param;
}



