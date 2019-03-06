#include "RecoDataProducts/inc/TriggerAlg.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& TriggerAlgDetail::typeName()
  {
     static std::string type("TriggerAlg");
     return type;
  }

  std::map<std::string,TriggerAlgDetail::mask_type> const& TriggerAlgDetail::bitNames()
  {
     static std::map<std::string,mask_type> bitnames;
     if (bitnames.size()==0)
     {
         bitnames[std::string("Unbiased")]             = bit_to_mask(unbiased);
         bitnames[std::string("MinBiasSrawDigiCount")] = bit_to_mask(minBiasSrawDigiCount);
         bitnames[std::string("LargeStrawDigiCount")]  = bit_to_mask(largeStrawDigiCount);
         bitnames[std::string("MinBiasCaloDigiCount")] = bit_to_mask(minBiasCaloDigiCount);
         bitnames[std::string("LargeCaloDigiCount")]   = bit_to_mask(largeCaloDigiCount);
         bitnames[std::string("CaloMVACE")]            = bit_to_mask(caloMVACE);
         bitnames[std::string("CaloMVAMixedCE")]       = bit_to_mask(caloMVAMixedCE);
         bitnames[std::string("CaloMVAPhoton")]        = bit_to_mask(caloMVAPhoton);
         bitnames[std::string("CaloLHCE")]             = bit_to_mask(caloLHCE);
         bitnames[std::string("CaloLHPhoton")]         = bit_to_mask(caloLHPhoton);
         bitnames[std::string("CaloCalibCosmic")]      = bit_to_mask(caloCalibCosmic);
         bitnames[std::string("CaloCalibLaser")]       = bit_to_mask(caloCalibLaser);
         bitnames[std::string("TprTimeClusterDeM")]    = bit_to_mask(tprTimeClusterDeM);
         bitnames[std::string("TprTimeClusterDeP")]    = bit_to_mask(tprTimeClusterDeP);
         bitnames[std::string("TprHelixDeM")]          = bit_to_mask(tprHelixDeM);
         bitnames[std::string("TprHelixDeP")]          = bit_to_mask(tprHelixDeP);
         bitnames[std::string("TprSeedDeM")]           = bit_to_mask(tprSeedDeM);
         bitnames[std::string("TprSeedDeP")]           = bit_to_mask(tprSeedDeP);
         bitnames[std::string("CprTimeClusterDeM")]    = bit_to_mask(cprTimeClusterDeM);
         bitnames[std::string("CprTimeClusterDeP")]    = bit_to_mask(cprTimeClusterDeP);
         bitnames[std::string("CprHelixDeM")]          = bit_to_mask(cprHelixDeM);
         bitnames[std::string("CprHelixDeP")]          = bit_to_mask(cprHelixDeP);
         bitnames[std::string("CprSeedDeM")]           = bit_to_mask(cprSeedDeM);
         bitnames[std::string("CprSeedDeP")]           = bit_to_mask(cprSeedDeP);
         bitnames[std::string("TprHelixIPADeM")]       = bit_to_mask(tprHelixIPADeM);
         bitnames[std::string("TprHelixIPADeP")]       = bit_to_mask(tprHelixIPADeP);
     }
     return bitnames;
  }

}
