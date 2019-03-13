//Building the Muon ID detector-using code from ExtMonFNALMagnet.hh and DetectorSolenoidShieldingMaker.hh
//As of now it is a rectangular box with half sizes, rotation, and center defined

//Jackson Waters, 2018

#ifndef EXTMONFNALMUONID_HH
#define EXTMONFNALMUONID_HH

#include "Mu2eInterfaces/inc/Detector.hh"

#include <ostream>
#include <memory>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
  
namespace mu2e {

  namespace ExtMonFNAL { class ExtMonMaker; }
  namespace ExtMonFNAL { class ExtMon; }

  class ExtMonFNALMuonIDMaker;

  class ExtMonFNALMuonID : virtual public Detector {
  public:
    //Transverse size and Z-position similar to plane stack mother volume
    const std::vector<double>& motherTransverseHalfSize() const { return m_motherTransverseHalfSize; }
    const double motherStartZ() const { return m_motherStartZ; }
    const double motherEndZ() const { return m_motherEndZ; }
    //Positioning of Muon ID Detector
    CLHEP::Hep3Vector refPointInMu2e() const { return refPointInMu2e_; }
    CLHEP::HepRotation const& muonIDRotationInMu2e() const { return muonIDRotationInMu2e_; }
           
    double nominalMomentum() const { return nominalMomentum_; }

 ;

  private:
    ExtMonFNALMuonID();
    // An initialized instance of this class should be obtained via ExtMonFNALMuonIDMaker
    // Included these first two classes to make it similar to ExtMonFNALPlaneStack files
    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNAL::ExtMonMaker;
    friend class ExtMonFNALMuonIDMaker;


    // Mother Volume based off of ExtMonFNALPlaneStack
    std::vector<double> m_motherTransverseHalfSize;
    double m_motherStartZ;
    double m_motherEndZ;
    
    //Momentum, Reference point, and Rotation

    double nominalMomentum_;

    CLHEP::Hep3Vector  refPointInMu2e_;
    CLHEP::HepRotation muonIDRotationInMu2e_;
  

  };

  // std::ostream& operator<<(std::ostream& os, const ExtMonFNALMuonID& muid);

  }// namespace mu2e*/

#endif/*EXTMONFNALMUONID_HH*/

