#ifndef MCDataProducts_StrawDigiMC_hh
#define MCDataProducts_StrawDigiMC_hh
//
//  Summary of MC information used to create a StrawDigi
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
//art includes
#include "canvas/Persistency/Common/Ptr.h"
// C++ includes
#include <iostream>
#include <vector>
#include <array>

// Mu2e includes

namespace mu2e {

  class StrawDigiMC{

  public:

    StrawDigiMC();
    // construct from hitlets
    StrawDigiMC(StrawId sid, double wetime[2], 
	CLHEP::HepLorentzVector cpos[2], 
	art::Ptr<StrawGasStep> sgs[2], art::Ptr<StepPointMC> stepMC[2]);
    // temporary legacy construtctor fot testing FIXME!
    StrawDigiMC(StrawId sid, double wetime[2], 
	CLHEP::HepLorentzVector cpos[2], art::Ptr<StepPointMC> stepMC[2], std::vector<art::Ptr<StepPointMC> > const& stepmcs);

    // copy construcors
    StrawDigiMC(const StrawDigiMC& rhs); // default: don't copy art::Ptrs
    StrawDigiMC(const StrawDigiMC& rhs, art::Ptr<StepPointMC> stepMC[2] );
    // Accessors
    StrawId strawId() const { return _strawid; }
    double wireEndTime(StrawEnd strawend) const { return _wetime[strawend]; }

    CLHEP::HepLorentzVector const& clusterPosition(StrawEnd strawend) const { return _cpos[strawend]; }
    double driftDistance(StrawEnd strawend) const;
    double distanceToMid(StrawEnd strawend) const;
    art::Ptr<StepPointMC> const&  stepPointMC(StrawEnd strawend) const { return _stepMC[strawend]; }
    std::array<art::Ptr<StepPointMC>,StrawEnd::nends> const&  stepPointMCs() const { return _stepMC; }
    art::Ptr<StrawGasStep> const&  strawGasStep(StrawEnd strawend) const { return _sgs[strawend]; }
    StrawEnd earlyEnd() const { return (_wetime[StrawEnd::cal] < _wetime[StrawEnd::hv]) ?  StrawEnd::cal : StrawEnd::hv; }
    art::Ptr<StepPointMC> const&  earlyStepPointMC() const { return stepPointMC(earlyEnd()); }
    art::Ptr<StrawGasStep> const&  earlyStrawGasStep() const { return strawGasStep(earlyEnd()); }

    double energySum() const;// sum of all MC true energy deposited
    // energy sum of particle that triggered the discriminator
    double triggerEnergySum(StrawEnd strawend) const;
    // check if this digi was created by cross-talk
    bool isCrossTalk(StrawEnd strawend) const;
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:
    StrawId  _strawid;      // Straw sid
    std::array<CLHEP::HepLorentzVector,StrawEnd::nends> _cpos; // Positions of the clusters responsible
    // for the TDC firings on each end (can be the same)
    std::array<double,StrawEnd::nends> _wetime; // times at the wire ends of the signals which fired the TDC.  This needs double precision as neutrons can take seconds (!) to generate signals
    // in ns from event window marker

    std::array<art::Ptr<StrawGasStep>,StrawEnd::nends> _sgs;	  // Ptr into StrawGasStep collection of step which
    std::array<art::Ptr<StepPointMC>,StrawEnd::nends> _stepMC;	  // Ptr into StepPointMC collection of step which
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawDigiMC const& hit){
    hit.print(ost,false);
    return ost;
  }
   typedef std::vector<mu2e::StrawDigiMC> StrawDigiMCCollection;

} // namespace mu2e

#endif /* MCDataProducts_StrawDigiMC_hh */
