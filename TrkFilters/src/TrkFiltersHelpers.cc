// Andrei Gaponenko, 2012

#include "TrkFilters/inc/TrkFiltersHelpers.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <vector>

#include "TrkFilters/inc/PhiPrescalingParams.hh"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

#include "canvas/Utilities/InputTag.h"

template<>
bool fhicl::ParameterSet::get_if_present<PhiPrescalingParams>(std::string const & key, PhiPrescalingParams& value) const {
  fhicl::ParameterSet pset;
  const bool present = get_if_present<fhicl::ParameterSet>(key,pset);
  float amplitude = pset.get<float>("amplitude");
  float frequency = pset.get<float>("frequency");
  float phase     = pset.get<float>("phase"    );

  value = PhiPrescalingParams(amplitude, frequency, phase);
  return present;
}
