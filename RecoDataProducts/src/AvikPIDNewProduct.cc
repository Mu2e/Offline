//
//  $Id: 
//  $Author: 
//  $Date: 
//
//  cloned from RecoDataProducts/src/PIDProduct.cc by Vadim Rusu
//

#include "RecoDataProducts/inc/AvikPIDNewProduct.hh"

using namespace std;

namespace mu2e {

  // Constructors
  AvikPIDNewProduct::AvikPIDNewProduct() {}

  AvikPIDNewProduct::AvikPIDNewProduct(const AvikPIDNewProduct & p) {
    _trkID          = p._trkID;
    _nMatched       = p._nMatched;
    _nMatchedAll    = p._nMatchedAll;
    _nUsedOsH       = p._nUsedOsH;
    _nUsedSsH       = p._nUsedSsH;
    _nUsedOsD       = p._nUsedOsD;

    _logDedxProbEle = p._logDedxProbEle;
    _logDedxProbMuo = p._logDedxProbMuo;
    _drdsVadim      = p._drdsVadim;
    _drdsVadimErr   = p._drdsVadimErr;
    _drdsOs         = p._drdsOs    ;
    _drdsOsErr      = p._drdsOsErr ;
    _drdsSs         = p._drdsSs    ;
    _drdsSsErr      = p._drdsSsErr ;

    _sumAvik        = p._sumAvik   ; 
    _sq2Avik        = p._sq2Avik;
    _sumAvikOs      = p._sumAvikOs ;
  }

  // operator overloading
  AvikPIDNewProduct& AvikPIDNewProduct::operator= (const AvikPIDNewProduct & p) {
    _trkID          = p._trkID;
    _nMatched       = p._nMatched;
    _nMatchedAll    = p._nMatchedAll;
    _nUsedOsH       = p._nUsedOsH;
    _nUsedSsH       = p._nUsedSsH;
    _nUsedOsD       = p._nUsedOsD;

    _logDedxProbEle = p._logDedxProbEle;
    _logDedxProbMuo = p._logDedxProbMuo;

    _drdsVadim      = p._drdsVadim;
    _drdsVadimErr   = p._drdsVadimErr;
    _drdsOs         = p._drdsOs    ;
    _drdsOsErr      = p._drdsOsErr ;
    _drdsSs         = p._drdsSs    ;
    _drdsSsErr      = p._drdsSsErr ;

    _sumAvik        = p._sumAvik   ; 
    _sq2Avik        = p._sq2Avik;
    _sumAvikOs      = p._sumAvikOs ;

    return (*this);
  }


  void AvikPIDNewProduct::clear () {
    _trkID          = -1;
    _nMatched       = -1 ;
    _nMatchedAll    = -1 ;
    _nUsedOsH       = -1;
    _nUsedSsH       = -1;
    _nUsedOsD       = -1;

    _logDedxProbEle = -1.;
    _logDedxProbMuo = -1.;

    _drdsVadim      = -1.;
    _drdsVadimErr   = -1.;
    _drdsOs         = -1.;
    _drdsOsErr      = -1.;
    _drdsSs         = -1.;
    _drdsSsErr      = -1.;

    _sumAvik        = -1.; 
    _sq2Avik        = -1.;
    _sumAvikOs      = 1.e6;
  }

  void  AvikPIDNewProduct::init(int     TrkID           , 
				int     NMatched        , int   NMatchedAll  ,
				int     NUsedOsH     , int   NUsedSsH     ,    
				int     NUsedOsD     ,

				float   LogDedxProbEle, float   LogDedxProbMuo, 

				float   DrdsVadim    , float DrdsVadimErr ,
				float   DrdsOs       , float DrdsOsErr    ,
				float   DrdsSs       , float DrdsSsErr    ,

				float   SumAvik      , 
				float   Sq2Avik      , 
				float   SumAvikOs    )
  {
    _trkID         = TrkID       ;
    _nMatched      = NMatched    ;
    _nMatchedAll   = NMatchedAll ;
    _nUsedOsH      = NUsedOsH    ;
    _nUsedSsH      = NUsedSsH    ;
    _nUsedOsD      = NUsedOsD    ;

    _logDedxProbEle = LogDedxProbEle;
    _logDedxProbMuo = LogDedxProbMuo;

    _drdsVadim     = DrdsVadim   ;
    _drdsVadimErr  = DrdsVadimErr;

    _drdsOs        = DrdsOs      ;
    _drdsOsErr     = DrdsOsErr   ;
    _drdsSs        = DrdsSs      ;
    _drdsSsErr     = DrdsSsErr   ;
    _sumAvik       = SumAvik     ; 
    _sq2Avik       = Sq2Avik     ;
    _sumAvikOs     = SumAvikOs   ;
  }



} // end namespace mu2e


