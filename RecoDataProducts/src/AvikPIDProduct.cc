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
    _trkID            = p._trkID;
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
    _nUsedSsEle       = p._nUsedSsEle;
    _nUsedSsMuo       = p._nUsedSsMuo;
    _drdsSsEle        = p._drdsSsEle    ;
    _drdsSsEleErr     = p._drdsSsEleErr ;
    _drdsSsMuo        = p._drdsSsMuo    ;
    _drdsSsMuoErr     = p._drdsSsMuoErr ;
    _nUsedOsEle       = p._nUsedOsEle;
    _nUsedOsMuo       = p._nUsedOsMuo;
    _sumAvikOsEle     = p._sumAvikOsEle ;
    _sumAvikOsMuo     = p._sumAvikOsMuo ;
  }

  // operator overloading
  AvikPIDProduct& AvikPIDProduct::operator= (const AvikPIDProduct & p) {
    _trkID            = p._trkID;
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
    _nUsedSsEle       = p._nUsedSsEle;
    _nUsedSsMuo       = p._nUsedSsMuo;
    _drdsSsEle        = p._drdsSsEle    ;
    _drdsSsEleErr     = p._drdsSsEleErr ;
    _drdsSsMuo        = p._drdsSsMuo    ;
    _drdsSsMuoErr     = p._drdsSsMuoErr ;
    _nUsedOsEle       = p._nUsedOsEle;
    _nUsedOsMuo       = p._nUsedOsMuo;
    _sumAvikOsEle     = p._sumAvikOsEle ;
    _sumAvikOsMuo     = p._sumAvikOsMuo ;

    return (*this);
  }


  void AvikPIDProduct::clear () {
    _trkID            = -1;
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
    _nUsedSsEle       = -1;
    _nUsedSsMuo       = -1;
    _drdsSsEle        = -1.;
    _drdsSsEleErr     = -1.;
    _drdsSsMuo        = -1.;
    _drdsSsMuoErr     = -1.;
    _nUsedOsEle       = -1;
    _nUsedOsMuo       = -1;
    _sumAvikOsEle     = -1.;
    _sumAvikOsMuo     = -1.;
  }

  void  AvikPIDProduct::init(int     TrkID,
			     float   LogDedxProbEle  , float LogDedxProbMuo  , 
			     float   DrdsVadimEle    , float DrdsVadimEleErr ,
			     float   DrdsVadimMuo    , float DrdsVadimMuoErr ,
			     int     NMatched        , int   NMatchedAll     ,
			     float   SumAvikEle      , float SumAvikMuo      ,
			     float   Sq2AvikEle      , float Sq2AvikMuo      ,
			     float   DrdsOsEle       , float DrdsOsEleErr    ,
			     float   DrdsOsMuo       , float DrdsOsMuoErr    ,
			     int     NUsedSsEle      , int   NUsedSsMuo      ,    
			     float   DrdsSsEle       , float DrdsSsEleErr    ,
			     float   DrdsSsMuo       , float DrdsSsMuoErr    ,
			     int     NUsedOsEle      , int   NUsedOsMuo      ,    
			     float   SumAvikOsEle    , float SumAvikOsMuo
			     )
  {
    _trkID            = TrkID          ;
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
    _nUsedSsEle       = NUsedSsEle     ;
    _nUsedSsMuo       = NUsedSsMuo     ;
    _drdsSsEle        = DrdsSsEle      ;
    _drdsSsEleErr     = DrdsSsEleErr   ;
    _drdsSsMuo        = DrdsSsMuo      ;
    _drdsSsMuoErr     = DrdsSsMuoErr   ;
    _nUsedOsEle       = NUsedOsEle;
    _nUsedOsMuo       = NUsedOsMuo;
    _sumAvikOsEle     = SumAvikOsEle ;
    _sumAvikOsMuo     = SumAvikOsMuo ;
  }



} // end namespace mu2e


