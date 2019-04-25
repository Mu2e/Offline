// Read a SimParticle collection and create a GenParticleCollection from the end point of the former.
//
//
// $Id: FromSimParticleEndPoint_module.cc,v 1.8 2013/09/27 16:03:41 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/27 16:03:41 $
//
// Original author Gianni Onorato
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
#include "EventGenerator/inc/PiCaptureEffects.hh"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "EventGenerator/inc/GeneratorBase.hh"
#include "SeedService/inc/SeedService.hh"

using namespace std;

namespace mu2e {

  typedef vector<int> Vpdg;
  typedef vector<string> Vstring;

  class FromSimParticleEndPoint : public art::EDProducer {
  public:
    explicit FromSimParticleEndPoint(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& event) override;
    void setNtuplaInfo(art::Event& event, GenParticle& gen,TNtuple*& ntup);

    string _inModuleLabel;
    Vpdg _inPdgId;
    Vstring _inVolumes;
    string _inProcessCode;
    ProcessCode _inProcessCodeToLook;
    string _outGenId;
    GenId _outGenIdToCreate;
    double _piCaptureELow, _piCaptureEHi, _piCaptureNBins, _piCaptureProbInternalConversion;
    bool _doHistograms;
    int _diagLevel;
    bool _firstevent;
    vector<unsigned> _selVolumes;
    TNtuple* _ntup;
    PiCaptureEffects* _piCaptCreator;
    void readVolumesToSelect(art::Event& event);

 };


  FromSimParticleEndPoint::FromSimParticleEndPoint(fhicl::ParameterSet const& pset):
    EDProducer{pset},
    _inModuleLabel(pset.get<std::string>("inputG4ModuleLabel","g4run")),
    _inPdgId(pset.get<Vpdg>("inputPdgIds", Vpdg())),
    _inVolumes(pset.get<Vstring>("inputVolumes", Vstring())),
    _inProcessCode(pset.get<string>("inputProcessCode")),
    _outGenId(pset.get<string>("outputGenId")),
    _piCaptureELow(pset.get<double>("picaptureElow",38.2)),
    _piCaptureEHi(pset.get<double>("picaptureEhi",138.2)),
    _piCaptureNBins(pset.get<double>("picaptureNbins",1000)),
    _piCaptureProbInternalConversion(pset.get<double>("picaptureProbInternalConversion", 0.0069)),
    _doHistograms(pset.get<bool>("doHistograms",true)),
    _diagLevel(pset.get<int>("diagLevel", -1)),
    _firstevent(true)
  {
    produces<GenParticleCollection>();
    produces<GenSimParticleLink>();
    _inProcessCodeToLook = ProcessCode::findByName(_inProcessCode);
    _outGenIdToCreate    = GenId::findByName(_outGenId);

    if ( _doHistograms ){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "FromSPEndPoint" );

      _ntup = tfs->make<TNtuple>("fromspep","Generator information",
                                 "evt:E:x:y:z:p:costh:phi:t:pdgId:GenId");
    }

    if (_outGenIdToCreate == GenId::PiCaptureCombined) {
      auto& engine = createEngine(art::ServiceHandle<SeedService>()->getSeed());
      _piCaptCreator = new PiCaptureEffects(engine,
                                            _piCaptureProbInternalConversion,
                                            _piCaptureELow,
                                            _piCaptureEHi,
                                            _piCaptureNBins);
    }

    //Other possibilities to be implemented
    if (!(_outGenIdToCreate == GenId::PiCaptureCombined)) {
      throw cet::exception("MODEL") << "Option has been not implemented yet!\n";
    }

  }

  void FromSimParticleEndPoint::setNtuplaInfo(art::Event& event, GenParticle& gen,TNtuple*& ntup) {

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

  void FromSimParticleEndPoint::readVolumesToSelect(art::Event& event) {

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
        cout << "dead in ";
        for (size_t i = 0; i < _selVolumes.size(); ++i ) {
          cout << vols.at(_selVolumes[i]).name() << endl;
        }
        cout << "with process code ";
        cout << _inProcessCodeToLook.name() << " ";
        cout << "to create a " << _outGenIdToCreate.name() << endl;
      }
  }

  void FromSimParticleEndPoint::produce(art::Event& event) {


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

    // Loop on SimParticle of the previous run

    for(SimParticleCollection::const_iterator i=insims.begin(); i!=insims.end(); ++i) {

      SimParticle const& aParticle(i->second);

      //Check the Pdg of the SimParticle
      Vpdg::iterator pdgFinder = find(_inPdgId.begin(),_inPdgId.end(), aParticle.pdgId());
      if (pdgFinder != _inPdgId.end()) {

        //Check the end volume of the SimParticle
        bool findvolume = false;
        for (vector<unsigned>::iterator volumeFinder = _selVolumes.begin();
             volumeFinder != _selVolumes.end(); ++volumeFinder) {
          if ( *volumeFinder == aParticle.endVolumeIndex()) {
            findvolume = true;
            break;
          }
        }

        //Check the stopping ProcessCode of the SimParticle
        if (findvolume) {

          if (_diagLevel>-1) {
            cout << "find particle " << aParticle.pdgId();
            cout << " in the selected volume " << endl;
          }
          if (_inProcessCodeToLook == aParticle.stoppingCode() ) {

            if (_diagLevel>-1) cout << "with stopping code " << aParticle.stoppingCode().name() << endl;

            CLHEP::Hep3Vector pos = aParticle.endPosition();
            double time = aParticle.endGlobalTime();

            if (_outGenIdToCreate == GenId::PiCaptureCombined) {

              _piCaptCreator->defineOutput();

              if (_piCaptCreator->doPhoton()) {

                if (_diagLevel>1) cout << "creating photon" << endl;

                output->push_back(_piCaptCreator->outputGamma(pos, time));
                history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                                   art::Ptr<SimParticle>(inh, std::distance(insims.begin(), i)) );
                if (_doHistograms) {
                  setNtuplaInfo(event, output->back(),_ntup);
                }

              }

              if (_piCaptCreator->doElectron()) {

                if (_diagLevel>1) cout << "creating electron" << endl;

                output->push_back(_piCaptCreator->outputElec(pos, time));
                history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                                   art::Ptr<SimParticle>(inh, std::distance(insims.begin(), i)) );
                if (_doHistograms)
                  setNtuplaInfo(event, output->back(),_ntup);

              }

              if (_piCaptCreator->doPositron()) {

                if (_diagLevel>1) cout << "creating positron" << endl;

                output->push_back(_piCaptCreator->outputPosit(pos, time));
                history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                                   art::Ptr<SimParticle>(inh, std::distance(insims.begin(), i)) );
                if (_doHistograms)
                  setNtuplaInfo(event, output->back(),_ntup);

              }
            }
          }
        }
      }
    }

    event.put(std::move(output));
    event.put(std::move(history));

  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FromSimParticleEndPoint);
