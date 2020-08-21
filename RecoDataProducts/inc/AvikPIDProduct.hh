//
//
//  Original author Avik Laha

#ifndef RecoDataProducts_AvikPIDProduct_HH
#define RecoDataProducts_AvikPIDProduct_HH

#include <utility>

namespace mu2e {


  class AvikPIDProduct {

  private:
    int    _eleTrkID;			// electron track ID
    int    _muoTrkID;			// muon     track ID
					// dE/dX
    float  _logDedxProbEle;
    float  _logDedxProbMuo;
					// Vadim's dr/ds
    float  _drdsVadimEle;
    float  _drdsVadimEleErr;

    float  _drdsVadimMuo;
    float  _drdsVadimMuoErr;
					// Avik's part
    float  _sumAvikEle;			// sum of Avik's terms
    float  _sumAvikMuo;			// 

    int    _nMatched;
    int    _nMatchedAll;

    float  _sq2AvikEle;		// sum of local slope residuals, OS float ts 
    float  _sq2AvikMuo;		// 
					// global dR/dS , SS and OS float ts only; (R = residual)
    float  _drdsOsEle;
    float  _drdsOsEleErr;

    float  _drdsOsMuo;
    float  _drdsOsMuoErr;

    int    _nUsedOsEleH;
    int    _nUsedOsMuoH;

    float  _sumAvikOsEle;
    float  _sumAvikOsMuo;

    int    _nUsedOsEleD;
    int    _nUsedOsMuoD;

    float  _drdsSsEle;
    float  _drdsSsEleErr;

    float  _drdsSsMuo;
    float  _drdsSsMuoErr;

    int    _nUsedSsEleH;
    int    _nUsedSsMuoH;

  public:

    AvikPIDProduct(); 
    AvikPIDProduct (const AvikPIDProduct & p) ; 
    ~AvikPIDProduct() {} 
    AvikPIDProduct & operator = (const AvikPIDProduct & p) ;

    void   clear() ;

    void   init(int     EleTrkID        , int   MuoTrkID        ,
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
		int     NUsedOsEleD     , int   NUsedOsMuoD	);
  int                       fNUsedOsEleH;      // for fDrdsSsEle
  int                       fNUsedOsMuoH;

    int    eleTrkID       () const { return _eleTrkID; }
    int    muoTrkID       () const { return _muoTrkID; }

    float  logDedxProbEle () const { return _logDedxProbEle ; }
    float  logDedxProbMuo () const { return _logDedxProbMuo ; }

    float  drdsVadimEle   () const { return _drdsVadimEle; }
    float  drdsVadimEleErr() const { return _drdsVadimEleErr; }
    float  drdsVadimMuo   () const { return _drdsVadimMuo; }
    float  drdsVadimMuoErr() const { return _drdsVadimMuoErr; }

    int    nMatched       () const { return _nMatched; }
    int    nMatchedAll    () const { return _nMatchedAll; }
    float  sumAvikEle     () const { return _sumAvikEle; }
    float  sumAvikMuo     () const { return _sumAvikMuo; }
    float  sq2AvikEle     () const { return _sq2AvikEle; }
    float  sq2AvikMuo     () const { return _sq2AvikMuo; }

    float  drdsOsEle      () const { return _drdsOsEle ; }
    float  drdsOsEleErr   () const { return _drdsOsEleErr ; }
    float  drdsOsMuo      () const { return _drdsOsMuo ; }
    float  drdsOsMuoErr   () const { return _drdsOsMuoErr ; }

    int    nUsedSsEleH    () const { return _nUsedSsEleH; }
    int    nUsedSsMuoH    () const { return _nUsedSsMuoH; }

    float  drdsSsEle      () const { return _drdsSsEle ; }
    float  drdsSsEleErr   () const { return _drdsSsEleErr ; }
    float  drdsSsMuo      () const { return _drdsSsMuo ; }
    float  drdsSsMuoErr   () const { return _drdsSsMuoErr ; }

    int    nUsedOsEleH    () const { return _nUsedOsEleH; }
    int    nUsedOsMuoH    () const { return _nUsedOsMuoH; }

    float  sumAvikOsEle   () const { return _sumAvikOsEle; }
    float  sumAvikOsMuo   () const { return _sumAvikOsMuo; }

    int    nUsedOsEleD    () const { return _nUsedOsEleD; }
    int    nUsedOsMuoD    () const { return _nUsedOsMuoD; }
  };



} // end namespace mu2e


#endif
