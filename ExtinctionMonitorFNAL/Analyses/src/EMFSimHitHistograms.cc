// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke

#include "ExtinctionMonitorFNAL/Analyses/inc/EMFSimHitHistograms.hh"

#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

#include <sstream>

namespace mu2e {

  // Book histograms in the subdirectory, given by the relativePath; that path is
  // relative to the root TFileDirectory for the current module.
  void EMFSimHitHistograms::book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath)
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
    book (extmon, tfdir);
  }

  // Book the histograms.
  void EMFSimHitHistograms::book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir) {


    // assumes one sensor per plane
    const unsigned nsensors = extmon.up().nplanes() + extmon.dn().nplanes();
    for(unsigned sensor = 0; sensor < nsensors; ++sensor) {
      ExtMonFNALSensorId sid(sensor);
      std::ostringstream osname;
      osname<<"hitOnSensor_"<<sensor;
      std::ostringstream ostitle;
      ostitle<<"Local hit position for sensor "<<sid;

      hitPosition_[sid] = tfdir.make<TH2D>(osname.str().c_str(),
                                           ostitle.str().c_str(),
                                           400, -extmon.sensor().halfSize()[0], +extmon.sensor().halfSize()[0],
                                           400, -extmon.sensor().halfSize()[1], +extmon.sensor().halfSize()[1]
                                           );

      hitPosition_[sid]->SetOption("colz");
    }

  } // end EMFSimHitHistograms::book()

  void EMFSimHitHistograms::fill(const ExtMonFNALSimHitCollection& coll) {
    for(ExtMonFNALSimHitCollection::const_iterator i = coll.begin(); i != coll.end(); ++i) {
      hitPosition_[i->sensorId()]->Fill(i->localStartPosition().x(), i->localStartPosition().y());
    }
  } // end EMFSimHitHistograms::fill()

} // end namespace mu2e
