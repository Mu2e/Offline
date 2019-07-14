//
// Module to read SimParticles and create an ntuple
//
// Original author KLGenser based on previous readback modules
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "art/Framework/Principal/Provenance.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

namespace mu2e {

  class SimParticleAnalyzer : public art::EDAnalyzer {
  public:

    typedef SimParticleCollection::key_type key_type;

    explicit SimParticleAnalyzer(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _verbosityLevel(pset.get<int>("verbosityLevel",0)),
      _ntpssp(0),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run"))
    {
    }

    virtual ~SimParticleAnalyzer() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;
    int _verbosityLevel;

    TNtuple* _ntpssp;

    // mu2e::SimParticlemv_g4run__stoppedMuonsSingleStage

    // DataType_ModuleLabel_InstanceName_ProcessName


    // Module label of the g4 module that produced the particles
    std::string _g4ModuleLabel;

  };

  void SimParticleAnalyzer::beginJob(){

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntpssp = tfs->make<TNtuple>( "ntpssp", "SimParticle ntuple",
                                  "run:evt:trk:gid:pdg:pk:ppdg:"
                                  "sgt:spt:svid:sg4stat:scc:rcc:sx:sy:sz:spx:spy:spz:"
                                  "egt:ept:evid:eg4stat:esc:ex:ey:ez:epx:epy:epz:"
                                  "eKE:nstep"
                                  );
  }

  void SimParticleAnalyzer::beginRun(art::Run const& run){

  }

  void SimParticleAnalyzer::analyze(const art::Event& event) {

    ++_nAnalyzed;
    
    // ntuple buffer.
    float nt[_ntpssp->GetNvar()];

    art::Handle<SimParticleCollection> simPCH;
    event.getByLabel(_g4ModuleLabel, simPCH);

    const SimParticleCollection& simPC = *simPCH;

    if (_verbosityLevel >0) {
      cout << "SimParticleCollection has " << simPC.size() << " particles" << endl;
    }

    for (const auto& simPMVO : simPC) {

      const mu2e::SimParticle& simP = simPMVO.second;

      long unsigned int pkey = 0;
      int ppdg = 0;
      long unsigned int gid(0);

      art::Ptr<SimParticle> const& pptr = simP.parent();
      art::Ptr<GenParticle> const& gptr = simP.genParticle();

      if(pptr) {
        pkey = pptr->id().asUint();
        ppdg = pptr->pdgId();
      }

      if(gptr) {
        gid  = gptr->generatorId().id();
        ppdg = gptr->pdgId();
      }
      // Fill the ntuple.

      nt[ 0] = event.id().run();
      nt[ 1] = event.id().event();
      nt[ 2] = simP.id().asUint();
      nt[ 3] = gid;
      nt[ 4] = simP.pdgId();
      nt[ 5] = pkey;
      nt[ 6] = ppdg;
      nt[ 7] = simP.startGlobalTime();
      nt[ 8] = simP.startProperTime();
      nt[ 9] = simP.startVolumeIndex();
      nt[10] = simP.startG4Status();
      nt[11] = simP.creationCode().id();
      nt[12] = simP.originParticle().creationCode().id();
      nt[13] = simP.startPosition().x();
      nt[14] = simP.startPosition().y();
      nt[15] = simP.startPosition().z();
      nt[16] = simP.startMomentum().x();
      nt[17] = simP.startMomentum().y();
      nt[18] = simP.startMomentum().z();
      nt[19] = simP.endGlobalTime();
      nt[20] = simP.endProperTime();
      nt[21] = simP.endVolumeIndex();
      nt[22] = simP.endG4Status();
      nt[23] = simP.stoppingCode().id();
      nt[24] = simP.endPosition().x();
      nt[25] = simP.endPosition().y();
      nt[26] = simP.endPosition().z();
      nt[27] = simP.endMomentum().x();
      nt[28] = simP.endMomentum().y();
      nt[29] = simP.endMomentum().z();
      nt[30] = simP.endKineticEnergy();
      nt[31] = simP.nSteps();

      _ntpssp->Fill(nt);

    } // end loop over simparticles

  }

}  // end namespace mu2e

using mu2e::SimParticleAnalyzer;
DEFINE_ART_MODULE(SimParticleAnalyzer);
