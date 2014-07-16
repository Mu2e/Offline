#ifndef MCDataProducts_StrawDigiMC_hh
#define MCDataProducts_StrawDigiMC_hh
//
//  Summary of MC information used to create a StrawDigi
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "DataProducts/inc/StrawIndex.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
//art includes
#include "art/Persistency/Common/Ptr.h"
// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  class StrawDigiMC{

  public:

    StrawDigiMC();
    // construct from hitlets
    StrawDigiMC(StrawIndex index, double wetime[2], CLHEP::HepLorentzVector cpos[2], 
    art::Ptr<StepPointMC> stepMC[2], std::vector<art::Ptr<StepPointMC> > const& stepMCs);

    // Accessors
    StrawIndex strawIndex() const { return _strawIndex; }
    double wireEndTime(StrawDigi::TDCChannel itdc) const { 
      return _wetime[static_cast<size_t>(itdc)]; }

    CLHEP::HepLorentzVector const& clusterPosition(StrawDigi::TDCChannel itdc) const { 
      return _cpos[static_cast<size_t>(itdc)]; }

    double driftDistance(StrawDigi::TDCChannel itdc) const;
    double distanceToMid(StrawDigi::TDCChannel itdc) const;
    art::Ptr<StepPointMC> const&  stepPointMC(StrawDigi::TDCChannel itdc) const {
      return _stepMC[static_cast<size_t>(itdc)]; }

    bool hasTDC(StrawDigi::TDCChannel itdc) const { return !stepPointMC(itdc).isNull();} 

    std::vector<art::Ptr<StepPointMC> > const& stepPointMCs() const { return _stepMCs; }

    double energySum() const;// sum of all MC true energy deposited
    // energy sum of particle that triggered the discriminator
    double triggerEnergySum(StrawDigi::TDCChannel itdc=StrawDigi::zero) const;
    // check if this digi was created by cross-talk
    bool isCrossTalk(StrawDigi::TDCChannel itdc) const;
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:
    StrawIndex  _strawIndex;      // Straw index
  // the following should be an std::array<,2>, but that's not supported by CINT: FIXME!!
    CLHEP::HepLorentzVector _cpos[2]; // Positions of the clusters responsible
    // for the TDC firings on each end (can be the same)
    double _wetime[2]; // times at the wire ends of the signals which fired the TDC
    art::Ptr<StepPointMC> _stepMC[2];	  // Ptr into StepPointMC collection of step which
    // triggered the TDC
    std::vector<art::Ptr<StepPointMC> > _stepMCs; // All StepPointMCs which contributed to the waveform
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawDigiMC const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* MCDataProducts_StrawDigiMC_hh */
