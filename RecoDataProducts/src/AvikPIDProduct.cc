//
//  $Id: 
//  $Author: 
//  $Date: 
//
//  cloned from RecoDataProducts/src/PIDProduct.cc by Vadim Rusu
//

#include "RecoDataProducts/inc/AvikPIDProduct.hh"

using namespace std;

namespace mu2e {

  // Constructors
  AvikPIDProduct::AvikPIDProduct() {}

  AvikPIDProduct::AvikPIDProduct(const AvikPIDProduct & p) {
    _eleTrkID         = p._eleTrkID;
    _muoTrkID         = p._muoTrkID;
    _logDedxProbEle   = p._logDedxProbEle;
    _logDedxProbMuo   = p._logDedxProbMuo;
    _drdsVadimEle     = p._drdsVadimEle;
    _drdsVadimEleErr  = p._drdsVadimEleErr;
    _drdsVadimMuo     = p._drdsVadimMuo;
    _drdsVadimMuoErr  = p._drdsVadimMuoErr;
    _nMatched         = p._nMatched;
    _nMatchedAll      = p._nMatchedAll;
    _sumAvikEle       = p._sumAvikEle   ; 
    _sumAvikMuo       = p._sumAvikMuo   ; 
    _sq2AvikEle       = p._sq2AvikEle;
    _sq2AvikMuo       = p._sq2AvikMuo;
    _drdsOsEle        = p._drdsOsEle    ;
    _drdsOsEleErr     = p._drdsOsEleErr ;
    _drdsOsMuo        = p._drdsOsMuo    ;
    _drdsOsMuoErr     = p._drdsOsMuoErr ;
    _nUsedSsEleH      = p._nUsedSsEleH;
    _nUsedSsMuoH      = p._nUsedSsMuoH;
    _drdsSsEle        = p._drdsSsEle    ;
    _drdsSsEleErr     = p._drdsSsEleErr ;
    _drdsSsMuo        = p._drdsSsMuo    ;
    _drdsSsMuoErr     = p._drdsSsMuoErr ;
    _nUsedOsEleH      = p._nUsedOsEleH;
    _nUsedOsMuoH      = p._nUsedOsMuoH;
    _sumAvikOsEle     = p._sumAvikOsEle ;
    _sumAvikOsMuo     = p._sumAvikOsMuo ;
    _nUsedOsEleD      = p._nUsedOsEleD;
    _nUsedOsMuoD      = p._nUsedOsMuoD;
  }

  // operator overloading
  AvikPIDProduct& AvikPIDProduct::operator= (const AvikPIDProduct & p) {
    _eleTrkID         = p._eleTrkID;
    _muoTrkID         = p._muoTrkID;
    _logDedxProbEle   = p._logDedxProbEle;
    _logDedxProbMuo   = p._logDedxProbMuo;
    _drdsVadimEle     = p._drdsVadimEle;
    _drdsVadimEleErr  = p._drdsVadimEleErr;
    _drdsVadimMuo     = p._drdsVadimMuo;
    _drdsVadimMuoErr  = p._drdsVadimMuoErr;
    _nMatched         = p._nMatched;
    _nMatchedAll      = p._nMatchedAll;
    _sumAvikEle       = p._sumAvikEle   ; 
    _sumAvikMuo       = p._sumAvikMuo   ; 
    _sq2AvikEle       = p._sq2AvikEle;
    _sq2AvikMuo       = p._sq2AvikMuo;
    _drdsOsEle        = p._drdsOsEle    ;
    _drdsOsEleErr     = p._drdsOsEleErr ;
    _drdsOsMuo        = p._drdsOsMuo    ;
    _drdsOsMuoErr     = p._drdsOsMuoErr ;
    _nUsedSsEleH      = p._nUsedSsEleH;
    _nUsedSsMuoH      = p._nUsedSsMuoH;
    _drdsSsEle        = p._drdsSsEle    ;
    _drdsSsEleErr     = p._drdsSsEleErr ;
    _drdsSsMuo        = p._drdsSsMuo    ;
    _drdsSsMuoErr     = p._drdsSsMuoErr ;
    _nUsedOsEleH      = p._nUsedOsEleH;
    _nUsedOsMuoH      = p._nUsedOsMuoH;
    _sumAvikOsEle     = p._sumAvikOsEle ;
    _sumAvikOsMuo     = p._sumAvikOsMuo ;
    _nUsedOsEleD      = p._nUsedOsEleD;
    _nUsedOsMuoD      = p._nUsedOsMuoD;

    return (*this);
  }


  void AvikPIDProduct::clear () {
    _eleTrkID         = -1;
    _muoTrkID         = -1;
    _logDedxProbEle   = -1.;
    _logDedxProbMuo   = -1.;
    _drdsVadimEle     = -1.;
    _drdsVadimEleErr  = -1.;
    _drdsVadimMuo     = -1.;
    _drdsVadimMuoErr  = -1.;
    _nMatched         = -1 ;
    _nMatchedAll      = -1 ;
    _sumAvikEle       = -1.; 
    _sumAvikMuo       = -1.; 
    _sq2AvikEle       = -1.;
    _sq2AvikMuo       = -1.;
    _drdsOsEle        = -1.;
    _drdsOsEleErr     = -1.;
    _drdsOsMuo        = -1.;
    _drdsOsMuoErr     = -1.;
    _nUsedSsEleH       = -1;
    _nUsedSsMuoH       = -1;
    _drdsSsEle        = -1.;
    _drdsSsEleErr     = -1.;
    _drdsSsMuo        = -1.;
    _drdsSsMuoErr     = -1.;
    _nUsedOsEleH       = -1;
    _nUsedOsMuoH       = -1;
    _sumAvikOsEle     = 1.e6;
    _sumAvikOsMuo     = 1.e6;
    _nUsedOsEleD       = -1;
    _nUsedOsMuoD       = -1;
  }

  void  AvikPIDProduct::init(int     EleTrkID        , int   MuoTrkID        ,
			     float   LogDedxProbEle  , float LogDedxProbMuo  , 
			     float   DrdsVadimEle    , float DrdsVadimEleErr ,
			     float   DrdsVadimMuo    , float DrdsVadimMuoErr ,
			     int     NMatched        , int   NMatchedAll     ,
			     float   SumAvikEle      , float SumAvikMuo      ,
			     float   Sq2AvikEle      , float Sq2AvikMuo      ,
			     float   DrdsOsEle       , float DrdsOsEleErr    ,
			     float   DrdsOsMuo       , float DrdsOsMuoErr    ,
			     int     NUsedSsEleH     , int   NUsedSsMuoH     ,    
			     float   DrdsSsEle       , float DrdsSsEleErr    ,
			     float   DrdsSsMuo       , float DrdsSsMuoErr    ,
			     int     NUsedOsEleH     , int   NUsedOsMuoH     ,    
			     float   SumAvikOsEle    , float SumAvikOsMuo    ,
			     int     NUsedOsEleD     , int   NUsedOsMuoD     )
  {
    _eleTrkID         = EleTrkID       ;
    _muoTrkID         = MuoTrkID       ;
    _logDedxProbEle   = LogDedxProbEle ;
    _logDedxProbMuo   = LogDedxProbMuo ;
    _drdsVadimEle     = DrdsVadimEle   ;
    _drdsVadimEleErr  = DrdsVadimEleErr;
    _drdsVadimMuo     = DrdsVadimMuo   ;
    _drdsVadimMuoErr  = DrdsVadimMuoErr;
    _nMatched         = NMatched;
    _nMatchedAll      = NMatchedAll    ;
    _sumAvikEle       = SumAvikEle     ; 
    _sumAvikMuo       = SumAvikMuo     ; 
    _sq2AvikEle       = Sq2AvikEle;
    _sq2AvikMuo       = Sq2AvikMuo;
    _drdsOsEle        = DrdsOsEle      ;
    _drdsOsEleErr     = DrdsOsEleErr   ;
    _drdsOsMuo        = DrdsOsMuo      ;
    _drdsOsMuoErr     = DrdsOsMuoErr   ;
    _nUsedSsEleH      = NUsedSsEleH    ;
    _nUsedSsMuoH      = NUsedSsMuoH    ;
    _drdsSsEle        = DrdsSsEle      ;
    _drdsSsEleErr     = DrdsSsEleErr   ;
    _drdsSsMuo        = DrdsSsMuo      ;
    _drdsSsMuoErr     = DrdsSsMuoErr   ;
    _nUsedOsEleH      = NUsedOsEleH;
    _nUsedOsMuoH      = NUsedOsMuoH;
    _sumAvikOsEle     = SumAvikOsEle ;
    _sumAvikOsMuo     = SumAvikOsMuo ;
    _nUsedOsEleD      = NUsedOsEleD;
    _nUsedOsMuoD      = NUsedOsMuoD;
  }



} // end namespace mu2e


