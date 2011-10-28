//
// A place to make diagnostic histograms, tables etc for G4.
// This is called by G4_plugin at appropriate times.
//
// $Id: DiagnosticsG4.cc,v 1.3 2011/10/28 18:47:06 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:06 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Despite having member function names reminiscent of those in
//    module classes, this class is not a module.
//

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes
#include "Analyses/inc/DiagnosticsG4.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

// ROOT includes
#include "TH1F.h"

// C++ includes
#include <iomanip>

using namespace std;

namespace mu2e {

  // Cannot do much at c'tor time.
  DiagnosticsG4::DiagnosticsG4():
    hStatusValue_(0),
    hNG4Tracks_(0),
    hNG4Tracks1Sup_(0),
    hNKilledStep_(0),
    hRealTime_(0),
    hRealTimeWide_(0),
    hCPUTime_(0),
    hCPUTimeWide_(0),
    hNTrkSteps_(0),
    hNCalSteps_(0),
    hNCalROSteps_(0),
    hNCRVSteps_(0),
    hNFoilSteps_(0),
    hNVDetSteps_(0),
    hNTrajectories_(0),
    hNPhysVolumes_(0){
  }

  // Book histograms in the root directory for the current module.
  void DiagnosticsG4::book(){
    art::ServiceHandle<art::TFileService> tfs;
    book(*tfs);
  }

  // Book histograms in the subdirectory, given by the relativePath; that path is
  // relative to the root TFileDirectory for the current module.
  void DiagnosticsG4::book( std::string const& relativePath ){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( relativePath.c_str() );
    book (tfdir);
  }

  // Book histograms.
  void DiagnosticsG4::book(art::TFileDirectory& tfdir ){

    // Histogram binning definitions are sensible when the event generator creates a single particle.
    hStatusValue_   = tfdir.make<TH1F>( "hStatusValue",   "Non-Zero Values of the G4 Status",          20,  0.,    20.   );
    hNG4Tracks_     = tfdir.make<TH1F>( "hNG4Tracks",     "Number of Tracks Created by G4",           200,  0.,  2000.   );
    hNG4Tracks1Sup_ = tfdir.make<TH1F>( "hNG4Tracks1Sup", "Number of Tracks Created by G4, 1 Suppressed",
                                                                                                     200,  1.,   401.   );
    hNKilledStep_   = tfdir.make<TH1F>( "hNKilledStep",   "Number Killed by Step Limit",              100,  0.,   100.   );
    hRealTime_      = tfdir.make<TH1F>( "hRealTime",      "Wall Clock time/event (10 ms ticks);(ms)",  50,  0.,   500.   );
    hRealTimeWide_  = tfdir.make<TH1F>( "hRealTimeWide",  "Wall Clock time/event (10 ms ticks);(ms)", 100,  0., 10000.   );
    hCPUTime_       = tfdir.make<TH1F>( "hCPUTime",       "CPU  time/event(10 ms ticks);(ms)",         50,  0.,   500.   );
    hCPUTimeWide_   = tfdir.make<TH1F>( "hCPUTimeWide",   "CPU  time/event(10 ms ticks);(ms)",        100,  0., 10000.   );
    hNTrkSteps_     = tfdir.make<TH1F>( "hNTrkSteps",     "N Tracker StepPointMCs",                   200,  1.,  1201    );
    hNCalSteps_     = tfdir.make<TH1F>( "hNCalSteps",     "N Calorimeter Crystal StepPointMCs",       300,  1.,   301    );
    hNCalROSteps_   = tfdir.make<TH1F>( "hNCalROSteps",   "N Calorimeter Readout StepPointMCs",        50,  1.,    51    );
    hNCRVSteps_     = tfdir.make<TH1F>( "hNCRVSteps",     "N CRV StepPointMCs",                       200,  1.,   201    );
    hNFoilSteps_    = tfdir.make<TH1F>( "hNFoilSteps",    "N Stopping Target Foil StepPointMCs",      200,  1.,   201    );
    hNVDetSteps_    = tfdir.make<TH1F>( "hNVDetSteps",    "N Virtual Detector StepPointMCs",          200,  1.,   201    );
    hNTrajectories_ = tfdir.make<TH1F>( "hNTrajectories", "Number of saved trajectories",              50,  1.,    51    );
    hNPhysVolumes_  = tfdir.make<TH1F>( "hNPhysVolumes",  "Number of Physical Volumes",               200,  1., 20001    );
  }

  void DiagnosticsG4::fillStatus( StatusG4 const& status){
    if ( status.status() != 0 ){
      hStatusValue_->Fill(status.status());
    }
    hNG4Tracks_->Fill(status.nG4Tracks());
    if ( status.nG4Tracks() > 1 ) hNG4Tracks1Sup_->Fill(status.nG4Tracks());
    if ( status.nKilledStepLimit() > 0 ) { hNKilledStep_->Fill(status.nKilledStepLimit()); }
  }

  void DiagnosticsG4::fill( StatusG4                     const& status,
                            SimParticleCollection        const& sims,
                            StepPointMCCollection        const& trkSteps,
                            StepPointMCCollection        const& calSteps,
                            StepPointMCCollection        const& calROSteps,
                            StepPointMCCollection        const& crvSteps,
                            StepPointMCCollection        const& foilSteps,
                            StepPointMCCollection        const& vdSteps,
                            PointTrajectoryCollection    const& trajectories,
                            PhysicalVolumeInfoCollection const& volInfo ){

    fillStatus(status);

    // Convert times to ms;
    double cpu = status.cpuTime()*1000.;
    double real = status.realTime()*1000.;

    hRealTime_    ->Fill(real);
    hRealTimeWide_->Fill(real);
    hCPUTime_     ->Fill(cpu);
    hCPUTimeWide_ ->Fill(cpu);

    hNTrkSteps_     ->Fill(trkSteps.size());
    hNCalSteps_     ->Fill(calSteps.size());
    hNCalROSteps_   ->Fill(calROSteps.size());
    hNCRVSteps_     ->Fill(crvSteps.size());
    hNFoilSteps_    ->Fill(foilSteps.size());
    hNVDetSteps_    ->Fill(vdSteps.size());
    hNTrajectories_ ->Fill(trajectories.size());
    hNPhysVolumes_  ->Fill(volInfo.size());

  }

}  // end namespace mu2e
