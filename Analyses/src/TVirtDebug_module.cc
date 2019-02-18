//
// Diagnostics on virtual detectors related to the Tracker.
//
//  $Id: TVirtDebug_module.cc,v 1.1 2013/12/20 20:05:12 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2013/12/20 20:05:12 $
//
// Original author Ivan Logashenko
//

#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include <iostream>
#include <string>
#include <set>

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  class TVirtDebug : public art::EDAnalyzer {
  public:

    typedef vector<int> Vint;

    explicit TVirtDebug(fhicl::ParameterSet const& pset);

    void beginRun( const art::Run& r);
    void analyze(const art::Event& e);

  private:

    // Name of the VD and TVD StepPoint collections
    std::string  _vdStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

  };

  TVirtDebug::TVirtDebug(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _vdStepPoints(StepInstanceName(StepInstanceName::virtualdetector).name()),
    _nAnalyzed(0),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run"))
  {}


  void TVirtDebug::beginRun(const art::Run& run) {
    GeomHandle<VirtualDetector> vdg;

    static std::set<int> vds = {
      VirtualDetectorId::TT_Mid,         VirtualDetectorId::TT_MidInner,  // 11, 12
      VirtualDetectorId::TT_FrontHollow, VirtualDetectorId::TT_FrontPA,   // 13, 14
      VirtualDetectorId::TT_Back,                                         // 15 Back (downstream) face of the tracker.
      VirtualDetectorId::TT_OutSurf,     VirtualDetectorId::TT_InSurf,    // 22-23:   external and internal surface of the Tracker envelope
    };

    for ( int id: vds ){
      VirtualDetectorId vid(id);

      CLHEP::Hep3Vector gpos = vdg->getGlobal(id);
      CLHEP::Hep3Vector lpos = vdg->getLocal(id);
      cout << "VDet: "
           << vid  << " "
           << gpos << " "
           << lpos
           << endl;
    }

  }

  void TVirtDebug::analyze(const art::Event& event) {

    ++_nAnalyzed;

    GeomHandle<VirtualDetector> vdg;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_vdStepPoints,hits);

    static std::set<int> vd_save = {
      VirtualDetectorId::TT_Mid,         VirtualDetectorId::TT_MidInner,  // 11, 12
      VirtualDetectorId::TT_FrontHollow, VirtualDetectorId::TT_FrontPA,   // 13, 14
      VirtualDetectorId::TT_Back,                                         // 15 Back (downstream) face of the tracker.
      VirtualDetectorId::TT_OutSurf,     VirtualDetectorId::TT_InSurf,    // 22-23:   external and internal surface of the Tracker envelope
    };

    for ( size_t i=0; i<hits->size(); ++i ){

      const StepPointMC& hit = (*hits)[i];

      int id = hit.volumeId();

      VirtualDetectorId vid(id);

      // If virtual detector id is not in the list - skip it
      if( vd_save.size()>0 && vd_save.find(id) == vd_save.end() ) continue;

      const CLHEP::Hep3Vector& gpos = hit.position();
      const CLHEP::Hep3Vector& lpos = gpos-vdg->getGlobal(id);

      cout << "VD Step: "
           << vid << " "
           << gpos << " "
           << lpos << " | "
           << lpos.perp()
           << endl;
    } // end loop over hits.

  } // end analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TVirtDebug);
