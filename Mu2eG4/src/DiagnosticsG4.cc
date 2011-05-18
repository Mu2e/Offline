//
// A place to make diagnostic histograms, tables etc for G4.
// This is called by G4_plugin at appropriate times.
//
// $Id: DiagnosticsG4.cc,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Despite having member function names reminiscent of those in
//    module classes, this class is not a module.
//

// Framework includes
#include "art/Framework/Core/Run.h"

#include "Mu2eG4/inc/DiagnosticsG4.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"

#include "TH1F.h"

namespace mu2e {

  // Cannot do much at c'tor time.
  DiagnosticsG4::DiagnosticsG4(){}

  void DiagnosticsG4::beginJob(){}
  void DiagnosticsG4::endJob(){}

  void DiagnosticsG4::beginRun( art::Run const &run, PhysicalVolumeHelper const& volInfo ){}
  void DiagnosticsG4::endRun(art::Run const& run){}

  void DiagnosticsG4::analyze( StatusG4                  const& status,
                               SimParticleCollection     const& sims,
                               StepPointMCCollection     const& trkSteps,
                               StepPointMCCollection     const& calSteps,
                               StepPointMCCollection     const& calROSteps,
                               StepPointMCCollection     const& crvSteps,
                               StepPointMCCollection     const& foilSteps,
                               StepPointMCCollection     const& vdSteps,
                               PointTrajectoryCollection const& trajectories,
                               PhysicalVolumeHelper      const& volInfo ){}

}  // end namespace mu2e
