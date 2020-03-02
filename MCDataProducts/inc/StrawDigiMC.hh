#ifndef MCDataProducts_StrawDigiMC_hh
#define MCDataProducts_StrawDigiMC_hh
//
//  Summary of MC information used to create a StrawDigi.  Everything is referenced by the thresold digitization end
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/XYZVec.hh"
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
    typedef art::Ptr<StrawGasStep> SGSP;
    typedef std::array<SGSP,StrawEnd::nends> SGSPA;
    typedef std::array<XYZVec,StrawEnd::nends> PA;
    typedef std::array<float,StrawEnd::nends> FA;

    StrawDigiMC();
    // construct from hitlets
    StrawDigiMC(StrawId sid, PA cpos, FA ctime, FA wetime, SGSPA sgs);

    // use compuater copy construcors
    StrawDigiMC(const StrawDigiMC& rhs, SGSPA sgsp ); // update the Ptrs
    // Accessors
    StrawId strawId() const { return _strawid; }

    SGSP const&  strawGasStep(StrawEnd strawend) const { return _sgspa[strawend]; }
    SGSPA const&  strawGasSteps() const { return _sgspa; }
    StrawEnd earlyEnd() const { return (_wtime[StrawEnd::cal] < _wtime[StrawEnd::hv]) ?  StrawEnd::cal : StrawEnd::hv; }
    SGSP const&  earlyStrawGasStep() const { return strawGasStep(earlyEnd()); }
    double energySum() const;// sum of all MC true energy deposited

    XYZVec const& clusterPos(StrawEnd strawend) const { return _cpos[strawend]; }
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
    CLHEP::HepLorentzVector clusterPosition(StrawEnd strawend) const { return CLHEP::HepLorentzVector(Geom::Hep3Vec(_cpos[strawend]),_ctime[strawend]); }

  private:
    StrawId  _strawid;      // Straw sid
    PA _cpos; // Positions of the trigger clusters
    FA _ctime; // times of the trigger clusters
    FA _wtime; // times at the wire ends of the signals which fired the TDC.
    SGSPA _sgspa; // StrawGasStep that triggered each end
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawDigiMC const& hit){
    hit.print(ost,false);
    return ost;
  }
   typedef std::vector<mu2e::StrawDigiMC> StrawDigiMCCollection;

} // namespace mu2e

#endif /* MCDataProducts_StrawDigiMC_hh */
