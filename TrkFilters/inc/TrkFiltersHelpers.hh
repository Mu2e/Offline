// Helper templates for extracting from fhicl::ParameterSet objects of
// types that are not directly supported by fhicl.  After including this
// file you will be able to do, for example
//
//    pset.get<PhiPrescalingParams>("prescaler")
//
// in you code.  The fhicl file syntax for specifying vectors is
//
//    prescaler : { amplitude : 0 frequency : 0 phase : 0 }
//
// Gianantonio Pezzullo, 2019

#ifndef TrkFilters_inc_TrkFiltersHelpers_hh
#define TrkFilters_inc_TrkFiltersHelpers_hh

// #include <set>
// #include <string>
// #include <vector>
#include "fhiclcpp/ParameterSet.h"
//#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

namespace art { class InputTag; }

namespace mu2e { class PhiPrescalingParams;}

namespace fhicl {
  
  // Read the parameters value need to define the phiPrescaling Function (sinusoidal function)
  // key : { amplitude : <amplitude>  frequencey : <frequency>  phase : <phase> }
  template<> bool ParameterSet::get_if_present<mu2e::PhiPrescalingParams>(std::string const & key, mu2e::PhiPrescalingParams&value) const;

}

#endif/*TrkFilters_inc_TrkFiltersHelpers_hh*/
