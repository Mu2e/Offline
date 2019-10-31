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

#ifndef GeneralUtilities_inc_ParameterSetHelpers_hh
#define GeneralUtilities_inc_ParameterSetHelpers_hh

#include <set>
#include <string>
#include <vector>
#include "fhiclcpp/ParameterSet.h"

namespace art { class InputTag; }
namespace CLHEP { class Hep3Vector; }
namespace mu2e { class PhiPrescalingParams; }

class Binning;

namespace fhicl {

  // Converting list of strings into vector of art::InputTag's
  template<> std::vector<art::InputTag> 
  ParameterSet::get<std::vector<art::InputTag>>(std::string const & key ) const;

  template<> std::vector<art::InputTag> 
  ParameterSet::get<std::vector<art::InputTag>>(std::string const & key, std::vector<art::InputTag> const& default_value ) const;

  // Converting list of ints into set of ints
  template<> std::set<int>
  ParameterSet::get<std::set<int>>(std::string const & key ) const;

  template<> std::set<int>
  ParameterSet::get<std::set<int>>(std::string const & key, std::set<int> const& default_value ) const;

  // Converting list of doubles to CLHEP::Hep3Vector
  template<> bool ParameterSet::get_if_present<CLHEP::Hep3Vector>(std::string const & key, CLHEP::Hep3Vector& value) const;

  // Read a binning for a histogram, represented by:
  // key : { n : <number of points>  low : <low edge of low bin>  high : <high edge of high bin> }
  template<> bool ParameterSet::get_if_present<Binning>(std::string const & key, Binning& value) const;

  //the function below should became a template function
  mu2e::PhiPrescalingParams getPhiPrescalerParams(fhicl::ParameterSet const&pset, std::string const& key);
}

#endif/*GeneralUtilities_inc_ParameterSetHelpers_hh*/
