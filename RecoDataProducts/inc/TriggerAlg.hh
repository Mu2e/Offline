#ifndef RecoDataProducts_TriggerAlg_hh
#define RecoDataProducts_TriggerAlg_hh

#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>


namespace mu2e {

  struct TriggerAlgDetail 
  {
    typedef unsigned int mask_type;
    
    enum bit_type { trig0,  trig1,   trig2,  trig3,
		    trig4,  trig5,   trig6,  trig7,
		    trig8,  trig9,  trig10, trig11,
		    trig12, trig13, trig14, trig15,
		    trig16, trig17, trig18, trig19,
		    trig20, trig21, trig22, trig23,
		    trig24, trig25, trig26, trig27,
		    trig28, trig29, trig30, trig31 };

    static const std::string& typeName();
    static const std::map<std::string,mask_type>& bitNames();
    static mask_type bit_to_mask(bit_type b) {return 1<<b;}
  };

  typedef BitMap<TriggerAlgDetail> TriggerAlg;
}


#endif

