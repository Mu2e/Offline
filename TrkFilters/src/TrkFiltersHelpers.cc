// Gianantonio Pezzullo, 2019

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "TrkFilters/inc/TrkFiltersHelpers.hh"

//#include <iostream>
//#include <sstream>
//#include <stdexcept>

//#include <vector>

#include "TrkFilters/inc/PhiPrescalingParams.hh"

#include "canvas/Utilities/InputTag.h"

template<>
bool fhicl::ParameterSet::get_if_present<mu2e::PhiPrescalingParams>(std::string const & key, mu2e::PhiPrescalingParams& value) const {
  fhicl::ParameterSet pset;
  const bool present = get_if_present<fhicl::ParameterSet>(key,pset);
  float amplitude = pset.get<float>("amplitude");
  float frequency = pset.get<float>("frequency");
  float phase     = pset.get<float>("phase"    );

  value = mu2e::PhiPrescalingParams(amplitude, frequency, phase);
  return present;
}
