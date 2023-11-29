// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke

#include "Offline/ExtinctionMonitorFNAL/Analyses/inc/EMFRawHitHistograms.hh"

#include "Offline/RecoDataProducts/inc/ExtMonFNALRawHit.hh"

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

#include <sstream>
#include <cmath>

namespace mu2e {

  EMFRawHitHistograms::EMFRawHitHistograms(const fhicl::ParameterSet& pset)
    : hitClock_()
    , hitToT_()
    , numClockTicksPerDebuncherPeriod_(GlobalConstantsHandle<PhysicsParams>()->getNominalDAQTicks())
  {}

  // Book histograms in the subdirectory, given by the relativePath; that path is
  // relative to the root TFileDirectory for the current module.
  void EMFRawHitHistograms::book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath)
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir(relativePath.c_str());
    book (extmon, tfdir);
  }

  // Book the histograms.
  void EMFRawHitHistograms::book(const ExtMonFNAL::ExtMon& extmon, const art::TFileDirectory& tfdir) {

    hitClock_ = tfdir.make<TH1D>("hitClock", "Hit clock, all hits",
                                 numClockTicksPerDebuncherPeriod_,
                                 -0.5, numClockTicksPerDebuncherPeriod_ - 0.5);

    hitToT_ =   tfdir.make<TH1D>("hitToT", "Hit time over threshold, all hits", 16, -0.5, 15.5);


    const unsigned nmodules = extmon.nmodules();
    ExtMonFNALModuleIdConverter con(extmon);
    for(unsigned mid = 0; mid < nmodules; ++mid) {
      for(unsigned ix = 0; ix < extmon.module().nxChips(); ++ix) {
        for(unsigned iy = 0; iy < extmon.module().nyChips(); ++iy) {
          ExtMonFNALChipId cid(con.moduleId(ExtMonFNALModuleDenseId(mid)), ix, iy);

          std::ostringstream osname;
          osname<<"chipOccupancy_"<< cid;
          std::ostringstream ostitle;
          ostitle<<"Occupancy for "<<cid;
          chipOccupancy_[cid] = tfdir.make<TH2D>(osname.str().c_str(),
                                                 ostitle.str().c_str(),
                                                 extmon.chip().nColumns(), -0.5, extmon.chip().nColumns() - 0.5,
                                                 extmon.chip().nRows(), -0.5, extmon.chip().nRows() - 0.5
                                                 );

          chipOccupancy_[cid]->SetOption("colz");
        }
      }
    }

  } // end EMFRawHitHistograms::book()

  void EMFRawHitHistograms::fill(const ExtMonFNALRawHitCollection& hits){

    for(ExtMonFNALRawHitCollection::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
      hitClock_->Fill(hit->clock());
      hitToT_->Fill(hit->tot());
      const ExtMonFNALPixelId& pix = hit->pixelId();
      chipOccupancy_[pix.chip()]->Fill(pix.col(), pix.row());

    }

  } // end EMFRawHitHistograms::fill()

} // end namespace mu2e
