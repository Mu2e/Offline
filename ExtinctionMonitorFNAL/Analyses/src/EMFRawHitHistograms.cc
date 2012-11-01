// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke

#include "ExtinctionMonitorFNAL/Analyses/inc/EMFRawHitHistograms.hh"

#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "ConditionsService/inc/ExtMonFNALConditions.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

#include <sstream>
#include <cmath>

namespace mu2e {

  EMFRawHitHistograms::EMFRawHitHistograms(const fhicl::ParameterSet& pset)
    : foldTimeToMicrobunch_(pset.get<bool>("foldTimeToMicrobunch"))
    , numMicrobunchTicks_()
    , hitClock_()
    , hitToT_()
  {}

  // Book histograms in the subdirectory, given by the relativePath; that path is
  // relative to the root TFileDirectory for the current module.
  void EMFRawHitHistograms::book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath)
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
    book (extmon, tfdir);
  }

  // Book the histograms.
  void EMFRawHitHistograms::book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir) {

    if(foldTimeToMicrobunch_) {
      ConditionsHandle<ExtMonFNALConditions> condEMF("ignored");
      ConditionsHandle<AcceleratorParams> condAcc("ignored");
      numMicrobunchTicks_ = rint(condAcc->deBuncherPeriod/condEMF->clockTick());
    }

    hitClock_ = tfdir.make<TH1D>("hitClock", "Hit clock, all hits", 100, -20.5, 79.5);
    hitToT_ =   tfdir.make<TH1D>("hitToT", "Hit time over threshold, all hits", 16, -0.5, 15.5);

    // assumes one sensor per plane
    const unsigned nsensors = extmon.up().nplanes() + extmon.dn().nplanes();
    for(unsigned sid = 0; sid < nsensors; ++sid) {
      for(unsigned ix = 0; ix < extmon.sensor().nxChips(); ++ix) {
        for(unsigned iy = 0; iy < extmon.sensor().nyChips(); ++iy) {
          ExtMonFNALChipId cid(ExtMonFNALSensorId(sid), ix, iy);

          std::ostringstream osname;
          osname<<"chipOccupancy_"<<sid<<"_"<<iy<<"_"<<ix;
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

      int clk = foldTimeToMicrobunch_ ? (hit->clock() % numMicrobunchTicks_): hit->clock();
      hitClock_->Fill(clk);
      hitToT_->Fill(hit->tot());
      const ExtMonFNALPixelId& pix = hit->pixelId();
      chipOccupancy_[pix.chip()]->Fill(pix.col(), pix.row());
    }

  } // end EMFRawHitHistograms::fill()

} // end namespace mu2e
