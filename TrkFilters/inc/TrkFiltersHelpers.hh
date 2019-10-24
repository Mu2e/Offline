// Helper templates for extracting from fhicl::ParameterSet objects of
// types that are not directly supported by fhicl.  After including this
// file you will be able to do, for example
//
//    pset.get<CLHEP::Hep3Vector>("position")
//
// in you code.  The fhicl file syntax for specifying vectors is
//
//    position : [ 1.1, 2.2, 3.3 ]
//
// Andrei Gaponenko, 2012

#ifndef TrkFilters_inc_TrkFiltersHelpers_hh
#define TrkFilters_inc_TrkFiltersHelpers_hh

#include <set>
#include <string>
#include <vector>
#include "fhiclcpp/ParameterSet.h"

namespace art { class InputTag; }

class PhiPrescalingParams;

namespace fhicl {
  
  // Read the parameters value need to define the phiPrescaling Function (sinusoidal function)
  // key : { amplitude : <amplitude>  frequencey : <frequency>  phase : <phase> }
  template<> bool ParameterSet::get_if_present<PhiPrescalingParams>(std::string const & key, PhiPrescalingParams&value) const;

}

#endif/*TrkFilters_inc_TrkFiltersHelpers_hh*/
