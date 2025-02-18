#ifndef MCDataProducts_StrawDigiMC_hh
#define MCDataProducts_StrawDigiMC_hh
//
//  Summary of MC information used to create a StrawDigi.  Everything is referenced by the threshold digitization end
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/MCDataProducts/inc/DigiProvenance.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
//art includes
#include "canvas/Persistency/Common/Ptr.h"
// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <string>

namespace mu2e {
  class StrawDigiMC{

  public:

    typedef art::Ptr<StrawGasStep> SGSP;
    typedef std::array<SGSP,StrawEnd::nends> SGSPA;
    typedef std::array<XYZVectorF,StrawEnd::nends> PA;
    typedef std::array<float,StrawEnd::nends> FA;

    StrawDigiMC();
    // construct from hitlets
    StrawDigiMC(StrawId sid, PA cpos, FA ctime, FA wetime, SGSPA sgs, DigiProvenance::enum_type=DigiProvenance::Simulation);
    // construct as noise
    StrawDigiMC(StrawId sid, bool isnoise);

    // use compuater copy construcors
    StrawDigiMC(const StrawDigiMC& rhs, SGSPA sgsp ); // update the Ptrs
    StrawDigiMC(const StrawDigiMC& rhs, DigiProvenance::enum_type provenance ); // update validity

    // Accessors
    StrawId strawId() const { return _strawid; }

    DigiProvenance provenance() const { return _provenance; }
    bool containsSimulation() const;
    bool isNoise() const { return _isnoise; }

    SGSP const&  strawGasStep(StrawEnd strawend) const { return _sgspa[strawend]; }
    SGSPA const&  strawGasSteps() const { return _sgspa; }
    SGSPA&        strawGasSteps()       { return _sgspa; }
    StrawEnd earlyEnd() const { return (_wtime[StrawEnd::cal] < _wtime[StrawEnd::hv]) ?  StrawEnd::cal : StrawEnd::hv; }
    SGSP const&  earlyStrawGasStep() const { return strawGasStep(earlyEnd()); }
    double energySum() const;// sum of all MC true energy deposited

    XYZVectorF const& clusterPos(StrawEnd strawend) const { return _cpos[strawend]; }
    // note that all the time functions below have time offsets already applied!
    float clusterTime(StrawEnd strawend) const { return _ctime[strawend]; } // time the ion cluster was created
    float wireEndTime(StrawEnd strawend) const { return _wtime[strawend]; } // time the signal reaches the end of the wire
    // energy sum of particle that triggered the discriminator
    double triggerEnergySum(StrawEnd strawend) const;
    // check if this digi was created by cross-talk
    bool isCrossTalk(StrawEnd strawend) const;
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

// legacy functions
    CLHEP::HepLorentzVector clusterPosition(StrawEnd strawend) const { return CLHEP::HepLorentzVector(GenVector::Hep3Vec(_cpos[strawend]),_ctime[strawend]); }

  private:
    StrawId  _strawid;      // Straw sid
    PA _cpos; // Positions of the trigger clusters
    FA _ctime; // times of the trigger clusters
    FA _wtime; // times at the wire ends of the signals which fired the TDC.
    SGSPA _sgspa; // StrawGasStep that triggered each end
    DigiProvenance _provenance; // level of association with any MC truth info
    bool _isnoise; // true if this digi was produced as noise; false otherwise
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawDigiMC const& hit){
    hit.print(ost,false);
    return ost;
  }
   typedef std::vector<mu2e::StrawDigiMC> StrawDigiMCCollection;

} // namespace mu2e

#endif /* MCDataProducts_StrawDigiMC_hh */
