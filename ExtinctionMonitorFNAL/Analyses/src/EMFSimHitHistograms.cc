// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke

#include "ExtinctionMonitorFNAL/Analyses/inc/EMFSimHitHistograms.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"

#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "cetlib_except/exception.h"

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
    hitTimes_ = tfdir.make<TH2D>("hitTimes", "Module plane vs hit time",
                                 400, -0.5, 399.5, /* ns resolution */
                                 3, -0.5, 3-0.5
                                 );
    hitTimes_->SetOption("colz");

    moduleHits_ = tfdir.make<TH1D>("mhits","Hits per module", 150, 0., 20);
    
    energyDeposit_ = tfdir.make<TH1D>("eion", "Ionizing energy deposit", 150, 0., 0.150);
    energyDeposit_->GetXaxis()->SetTitle("energy [MeV]");
 
    ExtMonFNALModuleIdConverter con(extmon);

    for(unsigned iplane=0; iplane<extmon.nplanes(); ++iplane) {
         
      ExtMonFNALPlane plane = extmon.plane(iplane);
           
      for(unsigned imod = 0; imod < plane.nModules(); ++imod) {
        ExtMonFNALModuleId mid = con.getModuleId(iplane,imod);
     
        std::ostringstream osname;
        osname<<"hitOnModule_"<<mid;
        std::ostringstream ostitle;
        ostitle<<"Local hit position for module "<<mid;
        
        hitPosition_[mid] = tfdir.make<TH2D>(osname.str().c_str(),
                                             ostitle.str().c_str(),
                                             400, -extmon.module().sensorHalfSize()[0], +extmon.module().sensorHalfSize()[0],
                                             400, -extmon.module().sensorHalfSize()[1], +extmon.module().sensorHalfSize()[1]
                                             );
        
        hitPosition_[mid]->SetOption("colz");
        
      } // for (unsigned imod..)
    } // for (unsigned iplane..)

  } // end EMFSimHitHistograms::book()
  
  void EMFSimHitHistograms::fill(const ExtMonFNAL::ExtMon& extmon, const ExtMonFNALSimHitCollection& coll) {
    ExtMonFNALModuleIdConverter con(extmon);
    for(ExtMonFNALSimHitCollection::const_iterator i = coll.begin(); i != coll.end(); ++i) {
      moduleHits_->Fill(con.denseModuleNumber(i->moduleId()).number());
      hitTimes_->Fill(i->startTime(), i->moduleId().plane());
      energyDeposit_->Fill(i->ionizingEnergyDeposit());
      if(hitPosition_[i->moduleId()] != NULL)
        hitPosition_[i->moduleId()]->Fill(i->localStartPosition().x(), i->localStartPosition().y());
      else throw cet::exception("RANGE")<<"module " << i->moduleId() << " not within range of booked histograms\n";
      
    }
  } // end EMFSimHitHistograms::fill()

} // end namespace mu2e
