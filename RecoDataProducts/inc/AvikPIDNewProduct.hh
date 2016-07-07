///////////////////////////////////////////////////////////////////////////////
//  $Id: 
//  $Author: 
//  $Date: 
///////////////////////////////////////////////////////////////////////////////
#ifndef RecoDataProducts_AvikPIDNewProduct_HH
#define RecoDataProducts_AvikPIDNewProduct_HH

#include <utility>

namespace mu2e {


  class AvikPIDNewProduct {

  private:
    int    _trkID;			// track ID
    int    _nMatched;
    int    _nMatchedAll;
    int    _nUsedOsH;
    int    _nUsedSsH;
    int    _nUsedOsD;
					// dE/dX - electron and muon hypotheses
    float  _logDedxProbEle;
    float  _logDedxProbMuo;
					// Vadim's dr/ds - all hits
    float  _drdsVadim;
    float  _drdsVadimErr;
    float  _drdsOs;			// dr/ds, OS doublets only
    float  _drdsOsErr;
    float  _drdsSs;			// dr/ds, SS doublets only
    float  _drdsSsErr;
					// Avik's part
    float  _sumAvik;			// sum of Avik's terms
    float  _sq2Avik;	                // sum of local slope residuals, OS float ts 
    float  _sumAvikOs;

  public:

    AvikPIDNewProduct(); 
    AvikPIDNewProduct (const AvikPIDNewProduct & p) ; 
    ~AvikPIDNewProduct() {} 
    AvikPIDNewProduct & operator = (const AvikPIDNewProduct & p) ;

    void   clear() ;

    void   init(int     TrkID         , int   NMatched      , int   NMatchedAll,
		int     NUsedSsH      , int   NUsedOsH      , int   NUsedOsD   ,
		float   LogDedxProbEle, float LogDedxProbMuo, 
		float   DrdsVadim     , float DrdsVadimErr  ,
		float   DrdsOs        , float DrdsOsErr     ,
		float   DrdsSs        , float DrdsSsErr     ,
		float   SumAvik       , float Sq2Avik       , float SumAvikOs);

    int    trkID       () const { return _trkID; }
    int    nMatched    () const { return _nMatched; }
    int    nMatchedAll () const { return _nMatchedAll; }
    int    nUsedSsH    () const { return _nUsedSsH; }
    int    nUsedOsH    () const { return _nUsedOsH; }
    int    nUsedOsD    () const { return _nUsedOsD; }

    float  logDedxProbEle () const { return _logDedxProbEle ; }
    float  logDedxProbMuo () const { return _logDedxProbMuo ; }

    float  drdsVadim   () const { return _drdsVadim   ; }
    float  drdsVadimErr() const { return _drdsVadimErr; }

    float  drdsOs      () const { return _drdsOs    ; }
    float  drdsOsErr   () const { return _drdsOsErr ; }
    float  drdsSs      () const { return _drdsSs    ; }
    float  drdsSsErr   () const { return _drdsSsErr ; }

    float  sumAvik     () const { return _sumAvik;   }
    float  sq2Avik     () const { return _sq2Avik;   }
    float  sumAvikOs   () const { return _sumAvikOs; }
  };



} // end namespace mu2e


#endif
