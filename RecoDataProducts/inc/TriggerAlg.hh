#ifndef RecoDataProducts_TriggerAlg_hh
#define RecoDataProducts_TriggerAlg_hh

#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>


namespace mu2e {

  struct TriggerAlgDetail 
  {
    typedef unsigned int mask_type;

    enum bit_type {unbiased            , 
		   minBiasStrawDigiCount, largeStrawDigiCount, 
		   minBiasCaloDigiCount , largeCaloDigiCount ,
		   caloMVACE           , caloMVAPhoton     , 
		   caloMVAMixedCE      ,
		   caloLHCE            , caloLHPhoton       , 
		   caloCalibCosmic     , caloCalibLaser     ,
		   tprTimeClusterDeM   , tprTimeClusterDeP  ,  
		   tprHelixDeM         , tprHelixDeP        ,  
		   tprSeedDeM          , tprSeedDeP         ,  
		   cprTimeClusterDeM   , cprTimeClusterDeP  ,  
		   cprHelixDeM         , cprHelixDeP        ,  
		   cprSeedDeM          , cprSeedDeP         ,  
		   tprHelixIPADeM      , tprHelixIPADeP};

    static const std::string& typeName();
    static const std::map<std::string,mask_type>& bitNames();
    static mask_type bit_to_mask(bit_type b) {return 1<<b;}
  };

  typedef BitMap<TriggerAlgDetail> TriggerAlg;
}


#endif

