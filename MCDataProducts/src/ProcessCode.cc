//
// An enum-matched-to-names class used to indicate why a SimParticle was
// created and why it stopped. The class contains enum entries for all
// physics processes known in G4; it also contains an enum entry to indicate
// that the particle is a primary particle and other enum entries to
// indicate that a particle was killed in one of the user actions written by G4.
//
//
// Original author Rob Kutschke
//

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Offline/MCDataProducts/inc/ProcessCode.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {

  const char* ProcessCode::_name[] = { PROCESSCODE_NAMES };

  BOOST_STATIC_ASSERT(sizeof(ProcessCode::_name)/sizeof(char*) == ProcessCode::lastEnum);

  void ProcessCode::printAll( std::ostream& ost){
    ost << "List of known process codes: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      ost << setw(3) << i << " " << _name[i] << std::endl;
    }
  }

  ProcessCode ProcessCode::findByName ( std::string const& name){

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<size(); ++i ){
      if ( _name[i] == name ){
        return ProcessCode(enum_type(i));
      }
    }
    return ProcessCode(unknown);
  }

  // Return a vector of the codes that are mu2e specific.
  std::vector<ProcessCode> ProcessCode::mu2eCodes(){
    std::vector<ProcessCode> codes;
    for ( size_t i=0; i<size(); ++i ){
      string name = _name[i];
      if ( name.substr(0,4) == "mu2e" ){
        codes.push_back(ProcessCode(enum_type(i)));
      }
    }
    return codes;
  }

  bool ProcessCode::isFromProtonBeam(const ProcessCode& pcode) {
    static std::vector<ProcessCode::enum_type> pbeamcodes = {
      ProcessCode::mu2eMuonCaptureAtRest,
      ProcessCode::mu2eMuonCaptureAtRest,
      ProcessCode::mu2eMuonDecayAtRest,
      ProcessCode::mu2eCeMinusEndpoint,
      ProcessCode::mu2eCeMinusLeadingLog,
      ProcessCode::mu2eCePlusEndpoint,
      ProcessCode::mu2eDIOLeadingLog,
      ProcessCode::mu2eInternalRMC,
      ProcessCode::mu2eExternalRMC,
      ProcessCode::mu2eFlateMinus,
      ProcessCode::mu2eFlatePlus,
      ProcessCode::mu2eFlatPhoton,
      ProcessCode::mu2eCePlusLeadingLog,
      ProcessCode::mu2ePionCaptureAtRest,
      ProcessCode::mu2eExternalRPC,
      ProcessCode::mu2eInternalRPC,
      ProcessCode::mu2eCaloCalib,
      ProcessCode::mu2ePienu,
      ProcessCode::mu2eGammaConversion,
      ProcessCode::mu2eAntiproton,
      ProcessCode::mu2eFlatMuMinus,
      ProcessCode::mu2eFlatMuPlus };
    for(auto pbeamcode : pbeamcodes){
      if(pcode._id == pbeamcode) return true;
    }
    return false;
  }



}
