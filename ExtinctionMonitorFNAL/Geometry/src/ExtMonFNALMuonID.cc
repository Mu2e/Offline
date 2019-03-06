//Source file for new Muon ID Detector

//Jackson Waters, 2018

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMuonID.hh"


#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib_except/exception.h"

//================================================================
namespace mu2e {

mu2e::ExtMonFNALMuonID::ExtMonFNALMuonID()
  : m_motherTransverseHalfSize()
  , nominalMomentum_()
{}

  /*std::ostream& operator<<(std::ostream& os, const ExtMonFNALMuonID& muid) {
  return os<<"ExtMonFNALMuonID(halfsize="<<muid.outerHalfSize()
	    <<" )";
	    }*/
//================================================================
} //namespace mu2e
