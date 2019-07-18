// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke

#include "ExtinctionMonitorFNAL/Analyses/inc/EMFRecoClusterHistograms.hh"

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

#include <sstream>

namespace mu2e {

  EMFRecoClusterHistograms::EMFRecoClusterHistograms()
    : extmon_(0)
    , perPlaneClusterMultiplicity_()
    , perClusterHitMultiplicity_()
    , clusterClock_()
    , clusterWidth_()
    , clusterPosition_()
  {}

  // Book histograms in the subdirectory, given by the relativePath; that path is
  // relative to the root TFileDirectory for the current module.
  void EMFRecoClusterHistograms::book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath)
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
    book (extmon, tfdir);
  }

  // Book the histograms.
  void EMFRecoClusterHistograms::book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir) {

    extmon_ = &extmon;

    const unsigned nplanes = extmon.up().nplanes() + extmon.dn().nplanes();

    perPlaneClusterMultiplicity_ = tfdir.make<TH2D>("clustersPerPlane", "Plane number vs per-plane cluster multiplicity",
                                                    100, -0.5, 99.5, nplanes, -0.5, nplanes-0.5);

    perPlaneClusterMultiplicity_->SetOption("colz");
    perPlaneClusterMultiplicity_->GetXaxis()->SetTitle("multiplicity");
    perPlaneClusterMultiplicity_->GetYaxis()->SetTitle("plane number");

    perClusterHitMultiplicity_   = tfdir.make<TH2D>("hitsPerCluster", "Plane number vs per-cluster hit multiplicity",
                                                     10, -0.5, 9.5, nplanes, -0.5, nplanes-0.5);

    perClusterHitMultiplicity_->SetOption("colz");
    perClusterHitMultiplicity_->GetXaxis()->SetTitle("multiplicity");
    perClusterHitMultiplicity_->GetYaxis()->SetTitle("plane number");

    clusterClock_   = tfdir.make<TH2D>("clusterTime", "Plane number vs cluster clock",
                                       110, -10.5, 99.5, nplanes, -0.5, nplanes-0.5);

    clusterClock_->SetOption("colz");
    clusterClock_->GetXaxis()->SetTitle("clock");
    clusterClock_->GetYaxis()->SetTitle("plane number");

    clusterWidth_   = tfdir.make<TH2D>("clusterWidth", "Y vs X cluster width", 5, 0.5, 5.5, 5, 0.5, 5.5);
    clusterWidth_->SetOption("colz");
    clusterWidth_->GetXaxis()->SetTitle("xWidth, pixels");
    clusterWidth_->GetYaxis()->SetTitle("yWidth, pixels");

    // Bin cluster position histograms according to pixel size
    const unsigned nx = extmon.module().nxChips() * extmon.chip().nColumns();
    const double   lx = nx * extmon.chip().xPitch();
 
    clusterPosition_.resize(nplanes);
    for(unsigned plane=0; plane < nplanes; ++plane) {
      const unsigned ny = extmon.plane(plane).nModules()*extmon.module().nyChips() * extmon.chip().nRows();
      const double   ly = ny * extmon.chip().yPitch();
      std::ostringstream osname;
      osname<<"clusterPosition"<<plane;
      std::ostringstream ostitle;
      ostitle<<"Cluster position in plane "<<plane;
      clusterPosition_[plane] = tfdir.make<TH2D>(osname.str().c_str(),
                                                 ostitle.str().c_str(),
                                                 nx, -lx/2, +lx/2,
                                                 ny, -ly/2, +ly/2
                                                 );

      clusterPosition_[plane]->SetOption("colz");
    }

  } // end EMFRecoClusterHistograms::book()

  void EMFRecoClusterHistograms::fill(const ExtMonFNALRecoClusterCollection& ccoll){

    for(unsigned plane = 0; plane < ccoll.nplanes(); ++plane) {

      ExtMonFNALRecoClusterCollection::PlaneClusters pc = ccoll.clusters(plane);

      perPlaneClusterMultiplicity_->Fill(pc.size(), plane);

      for(unsigned i=0; i<pc.size(); ++i) {

        perClusterHitMultiplicity_->Fill(pc[i].raw()->hits().size(), plane);
        clusterClock_->Fill(pc[i].clock(), plane);
        clusterWidth_->Fill(pc[i].yWidth(), pc[i].xWidth());
        CLHEP::Hep3Vector stackPos(pc[i].position());
        clusterPosition_[plane]->Fill(stackPos.x(), stackPos.y());

      }
    }

  } // end EMFRecoClusterHistograms::fill()

} // end namespace mu2e
