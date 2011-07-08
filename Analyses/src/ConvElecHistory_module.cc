//
// A module to follow the conversion electron in the events
//
// $Id: ConvElecHistory_module.cc,v 1.1 2011/07/08 22:14:02 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/07/08 22:14:02 $
//
// Original author Gianni Onorato
//

#include "CLHEP/Units/PhysicalConstants.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "Mu2eG4/inc/ConvElecUtilities.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>

using namespace std;

namespace mu2e {


  class ConvElecHistory : public art::EDAnalyzer {
  public:
    explicit ConvElecHistory(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _tNtup(0),
      _nBadG4Status(0),
      _nOverflow(0),
      _nKilled(0),
      _totalcputime(0),
      _totalrealtime(0)
    {
    }
    virtual ~ConvElecHistory() {
    }
    virtual void beginJob();
    virtual void endJob();

    void analyze(art::Event const& e );

  private:

    void doTracker(art::Event const& evt);


    // Diagnostic level
    int _diagLevel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    double _minimumEnergy; //minimum energy deposition of hits

    TNtuple* _tNtup;

    int _nBadG4Status, _nOverflow, _nKilled;
    float _totalcputime, _totalrealtime;

  };


  void ConvElecHistory::beginJob( ) {
  }

  void ConvElecHistory::analyze(art::Event const& evt ) {

    static int ncalls(0);
    ++ncalls;

    art::Handle<StatusG4> g4StatusHandle;
    evt.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
        << "Aborting ConvElecHistory::analyze due to G4 status\n"
        << g4Status;
      return;
    }

    if (g4Status.overflowSimParticles()) {
      ++_nOverflow;
      mf::LogError("G4")
        << "Aborting ConvElecHistory::analyze due to overflow of particles\n"
        << g4Status;
      return;
    }

    if (g4Status.nKilledStepLimit() > 0) {
      ++_nKilled;
      mf::LogError("G4")
        << "Aborting ConvElecHistory::analyze due to nkilledStepLimit reached\n"
        << g4Status;
      return;
    }

    _totalcputime += g4Status.cpuTime();
    _totalrealtime += g4Status.realTime();

    if (ncalls == 1) {

      art::ServiceHandle<art::TFileService> tfs;

      // evt:run: (2)
      //gentime:genx:geny:genz:gene:genp:gencosth:genphi:genfoil (9)
      //nhit:totedep:firsthitx:hirsthity:firsthitz:firsthitp:firsthitcosth:firsthitphi:firsthittime (9)
      //xdead:ydead:zdead:timedead:voldead:ndau (6)

      _tNtup        = tfs->make<TNtuple>( "ConvElec", "ConvElec Info", "evt:run:gentime:genx:geny:genz:gene:genp:gencosth:genphi:genfoil:nhit:totedep:hit1x:hit1y:hit1z:hit1p:hit1costh:hit1phi:hit1time:xdead:ydead:zdead:timedead:voldead:ndau");
   }

    doTracker(evt);
    
  } // end of analyze

  void ConvElecHistory::endJob() {
    cout << "ConvElecHistory::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
	 << "\nConvElecHistory::endJob Number of overflow events "
         << "due to too many particles in G4: "
         << _nOverflow
	 << "\nConvElecHistory::endJob Number of events with killed particles "
         << "due to too many steps in G4: "
         << _nKilled
	 << "\nConvElecHistory::endJob total CpuTime "
         << _totalcputime
	 << "\nConvElecHistory::endJob total RealTime "
         << _totalrealtime
         << endl;
  }


  void ConvElecHistory::doTracker(art::Event const& evt) {

    ConvElecUtilities CEUt(evt,_generatorModuleLabel,
                           _g4ModuleLabel, _trackerStepPoints);

    if (!CEUt.hasStepPointMC()) return;

    float tntpArray[27];
    int idx(0);
    tntpArray[idx++] = evt.id().event();
    tntpArray[idx++] = evt.run();
    
    const GenParticle& genCE = CEUt.genConvElec();

    tntpArray[idx++] = genCE.time();
    tntpArray[idx++] = genCE.position().x();
    tntpArray[idx++] = genCE.position().y();
    tntpArray[idx++] = genCE.position().z();
    tntpArray[idx++] = genCE.momentum().e();
    tntpArray[idx++] = genCE.momentum().vect().mag();
    tntpArray[idx++] = genCE.momentum().cosTheta();
    tntpArray[idx++] = genCE.momentum().phi();

    GeomHandle<Target> target;
    float nfoil = 0;
    for (int i=0; i<target->nFoils(); ++i) {
      TargetFoil const& foil = target->foil(i);
      if (genCE.position().z() >= foil.center().z()-foil.halfThickness() &&
          genCE.position().z() <= foil.center().z()+foil.halfThickness() )
        nfoil = i;
    }
    tntpArray[idx++] = nfoil;
    tntpArray[idx++] = CEUt.hasStepPointMC();
    tntpArray[idx++] = CEUt.totEDep();

    StepPointMC const & hit1 = CEUt.firstHit();

    tntpArray[idx++] = hit1.position().x();
    tntpArray[idx++] = hit1.position().y();
    tntpArray[idx++] = hit1.position().z();
    tntpArray[idx++] = hit1.momentum().mag();
    tntpArray[idx++] = hit1.momentum().cosTheta();
    tntpArray[idx++] = hit1.momentum().phi();
    tntpArray[idx++] = hit1.time();

    const SimParticle& simCE = CEUt.simConvElec(); 

    tntpArray[idx++] = simCE.endPosition().x();
    tntpArray[idx++] = simCE.endPosition().y();
    tntpArray[idx++] = simCE.endPosition().z();
    tntpArray[idx++] = simCE.endGlobalTime();
    tntpArray[idx++] = simCE.endVolumeIndex();
    tntpArray[idx++] = simCE.daughters().size();
    tntpArray[idx++] = simCE.daughters().size();

    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);

    PhysicalVolumeInfo const& volInfo = volumes->at(simCE.endVolumeIndex());

    cout << "Event " << evt.id().event() << " : \nConversion Electron "
         << "dead in " << simCE.endPosition() << " in the volume "
         << volInfo.name() << '\t' << simCE.endVolumeIndex() << " because of "
         << simCE.stoppingCode() << endl;


    _tNtup->Fill(tntpArray);

  } // end of doTracker
  
}

using mu2e::ConvElecHistory;
DEFINE_ART_MODULE(ConvElecHistory);

