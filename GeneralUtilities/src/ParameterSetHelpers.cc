// Andrei Gaponenko, 2012

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"

#include "GeneralUtilities/inc/Binning.hh"
#include "GeneralUtilities/inc/PhiPrescalingParams.hh"

#include "canvas/Utilities/InputTag.h"

//-----------------------------------------------------------------
// std::vector<std::string> ====> std::vector<art::InputTag>
//-----------------------------------------------------------------
template<>
std::vector<art::InputTag> fhicl::ParameterSet::get<std::vector<art::InputTag>>( std::string const & key ) const {

  std::vector<std::string> stringset;
  const bool present = get_if_present<std::vector<std::string>>( key, stringset );
 
  if ( !present ) throw fhicl::exception(cant_find, key);

  std::vector<art::InputTag> value;
  value.reserve( stringset.size() );

  for ( const auto& is : stringset ) value.push_back( is ); // implicit conversion to art::InputTag
  return value;
}

//-----------------------------------------------------------------
// std::vector<std::string> ====> std::vector<art::InputTag>
//-----------------------------------------------------------------
template<>
std::vector<art::InputTag> fhicl::ParameterSet::get<std::vector<art::InputTag>>( std::string const & key, std::vector<art::InputTag> const & default_value ) const {

  std::vector<std::string> stringset;
  const bool present = get_if_present<std::vector<std::string>>( key, stringset );
  
  if ( !present) return default_value;

  std::vector<art::InputTag> value;
  value.reserve( stringset.size() );

  for ( const auto& is : stringset ) value.push_back( is ); // implicit conversion to art::InputTag
  return value;
}

//-----------------------------------------------------------------
// std::vector<std::string> ====> std::set<int>
//-----------------------------------------------------------------
template<>
std::set<int> fhicl::ParameterSet::get<std::set<int>>( std::string const & key ) const {

  std::vector<int> intset;
  const bool present = get_if_present<std::vector<int>>( key, intset );
 
  if ( !present ) throw fhicl::exception(cant_find, key);

  std::set<int> value;
  for ( const auto& is : intset ) {
    if ( !value.insert( is ).second ) { 
      std::ostringstream os;
      os << "Value << " << is << " >> is duplicated in input list for FHiCL parameter: << " << key << " >> !\n";
      throw std::runtime_error(os.str());
    }
  }
  return value;
}

//-----------------------------------------------------------------
// std::vector<std::string> ====> std::set<int>
//-----------------------------------------------------------------
template<>
std::set<int> fhicl::ParameterSet::get<std::set<int>>( std::string const & key, std::set<int> const & default_value ) const {

  std::vector<int> intset;
  const bool present = get_if_present<std::vector<int>>( key, intset );
  
  if ( !present) return default_value;

  std::set<int> value;
  for ( const auto& is : intset ) {
    if ( !value.insert( is ).second ) { 
      std::ostringstream os;
      os << "Value << " << is << " >> is duplicated in input list for FHiCL parameter: << " << key << " >> !\n";
      throw std::runtime_error(os.str());
    }
  }
  return value;
}

//-----------------------------------------------------------------
// std::vector<double> ====> CLHEP::Hep3Vector
//-----------------------------------------------------------------
template<>
bool fhicl::ParameterSet::get_if_present<CLHEP::Hep3Vector>(std::string const & key, CLHEP::Hep3Vector& value) const {
  std::vector<double> val;
  const bool present = get_if_present<std::vector<double> >(key,val);
  if(present) {
    if(val.size() != 3) {
      std::ostringstream os;
      os<<"Error converting std::vector<double> to CLHEP::Hep3Vector for key \""<<key<<"\": wrong input size = "<<val.size();
      throw std::runtime_error(os.str());
    }
    value = CLHEP::Hep3Vector(val[0], val[1], val[2]);
  }
  return present;
}

template<>
bool fhicl::ParameterSet::get_if_present<Binning>(std::string const & key, Binning& value) const {
  fhicl::ParameterSet pset;
  const bool present = get_if_present<fhicl::ParameterSet>(key,pset);
  double low   = pset.get<double>("low");
  double high  = pset.get<double>("high");
  int    nbins = pset.get<int>("nbins");
  value = Binning(nbins,low,high);
  return present;
}


mu2e::PhiPrescalingParams fhicl::getPhiPrescalerParams(fhicl::ParameterSet const&pset, std::string const& key){
  fhicl::ParameterSet scalerPset = pset.get<fhicl::ParameterSet>(key);

  float amplitude = scalerPset.get<float>("amplitude");
  float frequency = scalerPset.get<float>("frequency");
  float phase     = scalerPset.get<float>("phase"    );
  
  return mu2e::PhiPrescalingParams(amplitude, frequency, phase);
}
