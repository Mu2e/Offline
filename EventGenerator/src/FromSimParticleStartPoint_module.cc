// Read a SimParticle collection and create a GenParticleCollection from the start point of the former.
// This module is originally designed to restart anti-proton simulation on production target
// from its production vertex.
//
// Author Zhengyun You
//
//

#include <iostream>
#include <string>
#include <memory>
#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TNtuple.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/GenSimParticleLink.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "EventGenerator/inc/GeneratorBase.hh"
#include "SeedService/inc/SeedService.hh"

using namespace std;

namespace mu2e {

  typedef vector<int> Vpdg;
  typedef vector<string> Vstring;

  class FromSimParticleStartPoint : public art::EDProducer {
  public:
    explicit FromSimParticleStartPoint(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& event) override;
    void setNtuplaInfo(art::Event& event, GenParticle& gen,TNtuple*& ntup);

    string _inModuleLabel;
    Vpdg _inPdgId;
    Vstring _inVolumes;
    ProcessCode _inProcessCodeToLook;
    GenId _outGenIdToCreate;
    bool _doHistograms;
    int _diagLevel;
    bool _firstevent;
    vector<unsigned> _selVolumes;
    TNtuple* _ntup;
    void readVolumesToSelect(art::Event& event);

  };

  FromSimParticleStartPoint::FromSimParticleStartPoint(fhicl::ParameterSet const& pset):
    EDProducer{pset},
    _inModuleLabel(pset.get<std::string>("inputG4ModuleLabel","g4run")),
    _inPdgId(pset.get<Vpdg>("inputPdgIds", {})),
    _inVolumes(pset.get<Vstring>("inputVolumes", {})),
    _inProcessCodeToLook{ProcessCode::findByName(pset.get<string>("inputProcessCodeDrop"))},
    _outGenIdToCreate{GenId::findByName(pset.get<string>("outputGenId"))},
    _doHistograms(pset.get<bool>("doHistograms",true)),
    _diagLevel(pset.get<int>("diagLevel", -1)),
    _firstevent(true)
  {
    produces<GenParticleCollection>();
    produces<GenSimParticleLink>();

    if ( _doHistograms ){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "FromSPStartPoint" );

      _ntup = tfs->make<TNtuple>("fromspsp","Generator information",
                                 "evt:E:x:y:z:p:costh:phi:t:pdgId:GenId");
    }
  }

  void FromSimParticleStartPoint::setNtuplaInfo(art::Event& event, GenParticle& gen,TNtuple*& ntup) {

    float nt[11];
    nt[0] = event.id().event();
    nt[1] = gen.momentum().e();
    nt[2] = gen.position().x();
    nt[3] = gen.position().y();
    nt[4] = gen.position().z();
    nt[5] = gen.momentum().vect().mag();
    nt[6] = gen.momentum().cosTheta();
    nt[7] = gen.momentum().phi();
    nt[8] = gen.time();
    nt[9] = gen.pdgId();
    nt[10] = gen.generatorId().id();

    ntup->Fill(nt);
  }

  void FromSimParticleStartPoint::readVolumesToSelect(art::Event& event) {

      if (_inVolumes.size() == 0) return;
      art::Handle<PhysicalVolumeInfoCollection> volsHandle;
      event.getRun().getByLabel(_inModuleLabel,volsHandle);
      PhysicalVolumeInfoCollection const& vols(*volsHandle);


      for (size_t i=0; i<vols.size(); ++i) {
        PhysicalVolumeInfo const& theVol = vols.at(i);

        for (size_t j=0; j<_inVolumes.size(); ++j) {
          if (theVol.name().compare(0,_inVolumes[j].size(),_inVolumes[j]) == 0) {
            _selVolumes.push_back(i);
          }
        }
      }

      if (_diagLevel>-1) {
        cout << "Searching for all ";
        for (size_t i = 0; i < _inPdgId.size(); ++i ) {
          cout << _inPdgId[i] << " ";
        }
        cout << "start in ";
        for (size_t i = 0; i < _selVolumes.size(); ++i ) {
          cout << vols.at(_selVolumes[i]).name() << endl;
        }
        if (_selVolumes.size() == 0) {
          cout << "any volume" << endl;
        }
        cout << "with process code ";
        cout << _inProcessCodeToLook.name() << " ";
        cout << "to create a " << _outGenIdToCreate.name() << endl;
      }
  }

  void FromSimParticleStartPoint::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    std::unique_ptr<GenSimParticleLink> history(new GenSimParticleLink);
    art::ProductID gpc_pid = (event.getProductID<GenParticleCollection>());

    if (_firstevent) {

      readVolumesToSelect(event);
      _firstevent = false;
    }

    // The input collection
    art::Handle<mu2e::SimParticleCollection> inh;
    event.getByLabel(_inModuleLabel, inh);
    const SimParticleCollection& insims(*inh);

    if (_diagLevel>-1) {
      cout << "insims.size() " << insims.size() << endl;
      for ( SimParticleCollection::const_iterator i=inh->begin(); i!=inh->end(); ++i ){
        SimParticle const& sim = i->second;
          cout << "search simparticle id " << sim.id() << " pdg " << sim.pdgId() << endl;

        art::Ptr<SimParticle> sp_ptr = art::Ptr<SimParticle>(inh, sim.id().asInt());
        if (sp_ptr) {
          cout << "found simparticle id " << sp_ptr->id() << " pdg " << sp_ptr->pdgId()
            << " distance " << std::distance(inh->begin(), i)  << endl;
        }
        else {
          cout << "NOT found simparticle" << endl;
        }
      }
    }


    // Loop on SimParticle of the previous run
    if (_diagLevel>-1) cout << "Previous run SimParticle size " << insims.size() << endl;
    for(SimParticleCollection::const_iterator i=insims.begin(); i!=insims.end(); ++i) {

      SimParticle const& aParticle(i->second);
      //cout << "aParticle " << aParticle.id() << " pdg " << aParticle.pdgId() << endl;

      //Check the Pdg of the SimParticle
      Vpdg::iterator pdgFinder = find(_inPdgId.begin(),_inPdgId.end(), aParticle.pdgId());
      if (pdgFinder != _inPdgId.end()) {

        bool findvolume = false;
        if (_selVolumes.size() == 0) {
          findvolume = true;
        }
        else {

          for (vector<unsigned>::iterator volumeFinder = _selVolumes.begin();
            volumeFinder != _selVolumes.end(); ++volumeFinder) {
            if ( *volumeFinder == aParticle.endVolumeIndex()) {
              findvolume = true;
              break;
            }
          }
        }

        //Check the creation ProcessCode of the SimParticle
        if (findvolume) {

          if (_diagLevel>-1) {
            cout << "find particle " << aParticle.pdgId();
            cout << " in the selected volume " << endl;
          }
          if (_inProcessCodeToLook != aParticle.creationCode() ) {

            if (_diagLevel>-1) cout << "keep with creation code " << aParticle.creationCode().name() << endl;

            PDGCode::type pdgId = aParticle.pdgId();
            CLHEP::Hep3Vector pos = aParticle.startPosition();
            CLHEP::HepLorentzVector mom = aParticle.startMomentum();
            double time = aParticle.startGlobalTime();

            if (_diagLevel>-1) cout << "creating particle" << pdgId << endl;

            GenParticle outGen(pdgId, GenId::fromSimParticleStartPoint, pos, mom, time);

            output->push_back(outGen);

            //cout << "inh size " << (*inh).size() << endl;
            history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                               art::Ptr<SimParticle>(inh, aParticle.id().asInt()) );

            if (_doHistograms) {
              setNtuplaInfo(event, output->back(),_ntup);
            }
          }
          else {

            if (_diagLevel>-1) cout << "drop with creation code " << aParticle.creationCode().name() << endl;
          }
        }
      } // if (pdgFinder != _inPdgId.end())
    } // for

    event.put(std::move(output));
    event.put(std::move(history));

  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FromSimParticleStartPoint);
