#include "RecoDataProducts/inc/TriggerFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& TriggerFlagDetail::typeName()
  {
     static std::string type("TriggerFlag");
     return type;
  }

  std::map<std::string,TriggerFlagDetail::mask_type> const& TriggerFlagDetail::bitNames()
  {
     static std::map<std::string,mask_type> bitnames;
     if (bitnames.size()==0)
     {
         bitnames[std::string("PrescaleRandom")]  = bit_to_mask(prescaleRandom);
         bitnames[std::string("HelixTrkOk")]      = bit_to_mask(helixTrkOk);
         bitnames[std::string("SeedTrkOk")]      = bit_to_mask(seedTrkOk);
         bitnames[std::string("KalTrkOk")]      = bit_to_mask(kalTrkOk);
         bitnames[std::string("CaloClusterOK")]   = bit_to_mask(caloClusterOk);
         bitnames[std::string("AnotherTrigger")]  = bit_to_mask(AnotherTrigger);
     }
     return bitnames;
  }

}
