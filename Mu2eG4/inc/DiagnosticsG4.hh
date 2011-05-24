#ifndef Mu2eG4_DiagnosticsG4_hh
#define Mu2eG4_DiagnosticsG4_hh
//
// A place to make diagnostic histograms, tables etc for G4.
// This is called by G4_plugin at appropriate times.
//
// $Id: DiagnosticsG4.hh,v 1.5 2011/05/24 17:19:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:19:03 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Despite having member function names reminiscent of those in
//    module classes, this class is not a module.
//

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"

// Forward declarations.
class TH1F;

namespace art{
  class Run;
}

namespace mu2e {

  class DiagnosticsG4{

  public:

    DiagnosticsG4();
    // Accept compiler generated d'tor.  Class is not copyable; see private section.

    void beginJob();
    void endJob();

    void beginRun( art::Run const &run, PhysicalVolumeHelper const& volInfo );
    void endRun(art::Run const& run);

    void analyze( StatusG4                  const& status,
                  SimParticleCollection     const& sims,
                  StepPointMCCollection     const& trkSteps,
                  StepPointMCCollection     const& calSteps,
                  StepPointMCCollection     const& calROSteps,
                  StepPointMCCollection     const& crvSteps,
                  StepPointMCCollection     const& foilSteps,
                  StepPointMCCollection     const& vdSteps,
                  PointTrajectoryCollection const& trajectories,
                  PhysicalVolumeHelper      const& volInfo);


  private:

    // Not copyable or assignable.  These will not be implemented.
    DiagnosticsG4 ( DiagnosticsG4 const& rhs );
    DiagnosticsG4& operator=(DiagnosticsG4 const& rhs);

    // ROOT owns the pointees; do not delete the pointee.
    TH1F* hG4time_;
    TH1F* htrackTime_;
    TH1F* hSimSize_;

  }; // end class Diagnostics G4

} // end namespace mu2e

#endif /* Mu2eG4_DiagnosticsG4_hh */
