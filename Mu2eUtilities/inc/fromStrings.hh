#ifndef Mu2eUtilities_fromStrings_hh
#define Mu2eUtilities_fromStrings_hh
//
// A free function template that translates from std::vector<string>
// to a std::vector of an enum type specified as the template type argument.
//
// Additional functions to support the ParticleDataList concepts of display name
// code name and alias.
//
// Original Author Rob Kutschke
//

#include <string>
#include <vector>

namespace mu2e {

  // For any type that is constructable from an std::string,
  // take an std::vector of strings and return a std::vector of that type.
  template<typename T>
  std::vector<T> fromStrings( std::vector<std::string> const& vs ){
    std::vector<T> r;
    r.reserve(vs.size());
    for ( auto const& s : vs ){
      r.emplace_back(s);
    }
    return r;
  }

  // The functions below assume familiarly with GlobalConstantsService/inc/ParticleDataList.hh
  //
  // The above function template will allow you to specify particles using its "code name".
  // The functions belows will allow you to specify PDGCodes using the name or the alias.

  // Optimization to avoid repeated creatation of Handles to ParticleDataList
  mu2e::PDGCode::type PDGCodefromString( std::string const& s, mu2e::ParticleDataList const& pdt ){
    mu2e::ParticleData const& p = pdt.particle(s);
    return (static_cast<mu2e::PDGCode::type>(p.id()));
  }

  mu2e::PDGCode::type PDGCodefromString( std::string const& s ){
    mu2e::GlobalConstantsHandle<mu2e::ParticleDataList> pdt;
    return PDGCodefromString( s, *pdt);
  }

  // Specialziation of the templated function fromStrings for PDGCode::type.
  // Cannot do this PDGCode without maintenance on that class.
  std::vector<mu2e::PDGCode::type> PDGCodesfromStrings( std::vector<std::string> const& v){
    mu2e::GlobalConstantsHandle<mu2e::ParticleDataList> pdt;
    std::vector<mu2e::PDGCode::type> r;
    for ( auto const& s : v){
      r.emplace_back( PDGCodefromString( s, *pdt) );
    }
    return r;
  }

} // end namespace mu2e

#endif
