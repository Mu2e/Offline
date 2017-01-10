#ifndef RecoDataProducts_TriggerFlag_hh
#define RecoDataProducts_TriggerFlag_hh

#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>


namespace mu2e {

   struct TriggerFlagDetail 
   {
       typedef unsigned mask_type;

       enum bit_type {prescaleRandom=0, helixTrkOk=1, caloClusterOk=2, AnotherTrigger=3};

       static const std::string& typeName();
       static const std::map<std::string,mask_type>& bitNames();
       static mask_type bit_to_mask(bit_type b) {return 1<<b;}
   };

   typedef BitMap<TriggerFlagDetail> TriggerFlag;

}


#endif

