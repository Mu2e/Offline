#ifndef Analyses_DiagnosticsG4_hh
#define Analyses_DiagnosticsG4_hh
//
// A place to make diagnostic histograms, tables etc for G4.
// This is called by G4_plugin at appropriate times.
//
// $Id: DiagnosticsG4.hh,v 1.1 2013/05/30 18:40:35 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/30 18:40:35 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Despite having member function names reminiscent of those in
//    module classes, this class is not a module.
//

// Mu2e includes
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// Art includes
#include "art_root_io/TFileDirectory.h"

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

    // Book histograms at the root TFileDirectory for the current module.
    void book( );

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    void book( std::string const& relativePath );

    // Book histograms in the specified TFileDirectory.
    void book( art::TFileDirectory& tfdir );

    // Fill just the status part.
    void fillStatus( StatusG4 const& status);

    // Fill all information.
    void fill( StatusG4                     const& status,
               SimParticleCollection        const& sims,
               StepPointMCCollection        const& trkSteps,
               StepPointMCCollection        const& calSteps,
               StepPointMCCollection        const& calROSteps,
               StepPointMCCollection        const& crvSteps,
               StepPointMCCollection        const& foilSteps,
               StepPointMCCollection        const& vdSteps,
               StepPointMCCollection        const& extMonUCITofSteps,
               PointTrajectoryCollection    const& trajectories,
               PhysicalVolumeInfoCollection const& volInfo);

    // nullptr args are allowed
    void fill( StatusG4                     const* status,
               SimParticleCollection        const* sims,
               StepPointMCCollection        const* trkSteps,
               StepPointMCCollection        const* calSteps,
               StepPointMCCollection        const* calROSteps,
               StepPointMCCollection        const* crvSteps,
               StepPointMCCollection        const* foilSteps,
               StepPointMCCollection        const* vdSteps,
               StepPointMCCollection        const* extMonUCITofSteps,
               PointTrajectoryCollection    const* trajectories,
               PhysicalVolumeInfoCollection const* volInfo);

    // Fill information for PA only
    void fillPA ( StepPointMCCollection        const& paSteps);

    // nullptr arg is allowed
    void fillPA ( StepPointMCCollection        const* paSteps);

  private:

    // Not copyable or assignable.  These will not be implemented.
    DiagnosticsG4 ( DiagnosticsG4 const& rhs );
    DiagnosticsG4& operator=(DiagnosticsG4 const& rhs);

    // ROOT owns the pointees; do not delete the pointee.
    TH1F* hStatusValue_;
    TH1F* hNG4Tracks_;
    TH1F* hNG4TracksLog_;
    TH1F* hNG4Tracks1Sup_;
    TH1F* hNKilledStep_;
    TH1F* hRealTime_;
    TH1F* hRealTimeWide_;
    TH1F* hCPUTime_;
    TH1F* hCPUTimeWide_;
    TH1F* hNTrkSteps_;
    TH1F* hNCalSteps_;
    TH1F* hNCalROSteps_;
    TH1F* hNCRVSteps_;
    TH1F* hNFoilSteps_;
    TH1F* hNVDetSteps_;
    TH1F* hNExtMonUCITofSteps_;
    TH1F* hNTrajectories_;
    TH1F* hNPhysVolumes_;
    TH1F* hNPASteps_;


  }; // end class Diagnostics G4

} // end namespace mu2e

#endif /* Analyses_DiagnosticsG4_hh */
