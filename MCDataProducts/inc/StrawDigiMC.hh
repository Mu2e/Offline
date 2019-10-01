#ifndef MCDataProducts_StrawDigiMC_hh
#define MCDataProducts_StrawDigiMC_hh
//
//  Summary of MC information used to create a StrawDigi
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "DataProducts/inc/StrawId.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
//art includes
#include "canvas/Persistency/Common/Ptr.h"
// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  class StrawDigiMC{

  public:

    StrawDigiMC();
    // construct from hitlets
    StrawDigiMC(StrawId sid, double wetime[2], 
    CLHEP::HepLorentzVector cpos[2], 
    art::Ptr<StepPointMC> stepMC[2], std::vector<art::Ptr<StepPointMC> > const& stepMCs);

    // copy construcors
    StrawDigiMC(const StrawDigiMC& rhs); // default: don't copy art::Ptrs
    StrawDigiMC(const StrawDigiMC& rhs, art::Ptr<StepPointMC> stepMC[2], std::vector<art::Ptr<StepPointMC> > const& stepMCs);

    // Accessors
    StrawId strawId() const { return _strawid; }
    double wireEndTime(StrawEnd strawend) const { return _wetime[strawend]; }

    CLHEP::HepLorentzVector const& clusterPosition(StrawEnd strawend) const { return _cpos[strawend]; }
    double driftDistance(StrawEnd strawend) const;
    double distanceToMid(StrawEnd strawend) const;
    art::Ptr<StepPointMC> const&  stepPointMC(StrawEnd strawend) const { return _stepMC[strawend]; }
    std::vector<art::Ptr<StepPointMC> > const& stepPointMCs() const { return _stepMCs; }
    StrawEnd earlyEnd() const { return (_wetime[StrawEnd::cal] < _wetime[StrawEnd::hv]) ?  StrawEnd::cal : StrawEnd::hv; }
    art::Ptr<StepPointMC> const&  earlyStepPointMC() const { return stepPointMC(earlyEnd()); }

    double energySum() const;// sum of all MC true energy deposited
    // energy sum of particle that triggered the discriminator
    double triggerEnergySum(StrawEnd strawend) const;
    // check if this digi was created by cross-talk
    bool isCrossTalk(StrawEnd strawend) const;
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:
    StrawId  _strawid;      // Straw sid
  // the following should be an std::array<,2>, but that's not supported by CINT: FIXME!!
    CLHEP::HepLorentzVector _cpos[2]; // Positions of the clusters responsible
    // for the TDC firings on each end (can be the same)
    double _wetime[2]; // times at the wire ends of the signals which fired the TDC.  This needs double precision as neutrons can take seconds (!) to generate signals
    // in ns from event window marker

    art::Ptr<StepPointMC> _stepMC[2];	  // Ptr into StepPointMC collection of step which
    // triggered the TDC
    std::vector<art::Ptr<StepPointMC> > _stepMCs; // All StepPointMCs which contributed to the waveform
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawDigiMC const& hit){
    hit.print(ost,false);
    return ost;
  }
   typedef std::vector<mu2e::StrawDigiMC> StrawDigiMCCollection;

} // namespace mu2e

#endif /* MCDataProducts_StrawDigiMC_hh */
