#ifndef MCDataProducts_ScorerConfigSummary_hh
#define MCDataProducts_ScorerConfigSummary_hh

//
// Information about scorers managed by Geant4
//
// Original author Bertrand Echenard
//

#include "Offline/DataProducts/inc/GenVector.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include <string>


namespace mu2e {

  struct ScorerConfigSummary {

    // A default c'tor is required for ROOT.
    ScorerConfigSummary():
      name_(),
      nbinsX_(0u),
      nbinsY_(0u),
      nbinsZ_(0u),
      halfSize_(),
      center_()
    {}

    ScorerConfigSummary(const std::string& name, unsigned nbinsX, unsigned nbinsY, unsigned nbinsZ,
                        const CLHEP::Hep3Vector& halfSize, const CLHEP::Hep3Vector& center):
      name_(name),
      nbinsX_(nbinsX),
      nbinsY_(nbinsY),
      nbinsZ_(nbinsZ),
      halfSize_(halfSize),
      center_(center)
    {}

    const     std::string& name()     const {return name_;}
    unsigned  nbinsX()                const {return nbinsX_;}
    unsigned  nbinsY()                const {return nbinsY_;}
    unsigned  nbinsZ()                const {return nbinsZ_;}
    const     XYZVectorF&  halfSize() const {return halfSize_;}
    const     XYZVectorF&  center()   const {return center_;}


  private:
    std::string  name_;
    unsigned     nbinsX_;
    unsigned     nbinsY_;
    unsigned     nbinsZ_;
    XYZVectorF   halfSize_;
    XYZVectorF   center_;
  };

  typedef std::vector<mu2e::ScorerConfigSummary> ScorerConfigSummaryCollection;
}

#endif
