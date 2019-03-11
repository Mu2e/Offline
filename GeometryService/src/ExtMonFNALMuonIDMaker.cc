// Jackson Waters, 2018
//Maker source file for Muon ID detector geometry

// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMuonID.hh"
#include "GeometryService/inc/ExtMonFNALMuonIDMaker.hh"

// C++ includes
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

//Framework includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"



//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl; 
#define AGDEBUG(stuff)

//using namespace std;

namespace mu2e {

   std::unique_ptr<ExtMonFNALMuonID> ExtMonFNALMuonIDMaker::make(const SimpleConfig& c) {
     std::unique_ptr<ExtMonFNALMuonID> muid ( new ExtMonFNALMuonID() );

     return muid;
   }

  ExtMonFNALMuonID
  ExtMonFNALMuonIDMaker::read(const SimpleConfig& c, const std::string& prefix,
                             const CLHEP::HepRotation& muonIDInRotationInMu2e, // of the input arm of ref trajectory
			      const CLHEP::Hep3Vector& refTrajmuonIDEntranceInMu2e,
                             double nominalMomentum)
  {
    ExtMonFNALMuonID muid;
    muid.muonIDRotationInMu2e_ = muonIDInRotationInMu2e;
    muid.nominalMomentum_ = nominalMomentum;

    c.getVectorDouble(prefix + ".motherTransverseHalfSize", muid.m_motherTransverseHalfSize, 2);
 
    // distance between the points
   //const double MuonIDEntranceToBendPointDistance = muid.outerHalfSize_[2]/cos(trackBendHalfAngle);

    muid.refPointInMu2e_ = refTrajmuonIDEntranceInMu2e;
    muid.m_motherStartZ = c.getDouble(prefix+".motherStartZ");
    muid.m_motherEndZ = c.getDouble(prefix+".motherEndZ");
   // muid.geometricCenterInMu2e_ = muid.refPointInMu2e_ + muid.muonIDRotationInMu2e_ * CLHEP::Hep3Vector(0, 0, -2000);
      
   if(c.getInt("extMonFNAL.verbosityLevel") > 0) {
     std::cout<<"ExtMonFNALMuonID "<<prefix<<": refTrajMuonIDEntranceInMu2e = "<< refTrajmuonIDEntranceInMu2e<<std::endl;
     std::cout<<"ExtMonFNALMuonID "<<prefix<<": refPointInMu2e = "<<muid.refPointInMu2e_<<std::endl;
     std::cout<<"ExtMonFNALMuonID "<<prefix<<": Momentum = "<<muid.nominalMomentum_<<std::endl;
     std::cout<<"ExtMonFNALMuonID "<<prefix<<": motherTransverseHalfSize[0] = "<<muid.m_motherTransverseHalfSize[0]<<std::endl;  
     std::cout<<"ExtMonFNALMuonID "<<prefix<<": motherTransverseHalfSize[1] = "<<muid.m_motherTransverseHalfSize[1]<<std::endl; 
  }

    return muid;

    
  } 

} // namespace mu2e

