//Author: S Middleton
//Date April 2020
//Purpose: Analyzer for stopped pions at Stage 3
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "Mu2eUtilities/inc/SimParticleGetTau.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TNtuple.h"
#include<TH1F.h>
#include<TH2F.h>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e{

   struct StopInfo {
      float sx;
      float sy;
      float sz;
      float spx;
      float spy;
      float spz;
      float spt;
      float st;
      float stau; // proper time, for stopped pion weights
      ProcessCode scode;
      PDGCode::type sid;
      StopInfo() : sx(), sy(), sz(), spx(), spy(), spz(), spt(), st(), stau(), scode(), sid() {}

      StopInfo(const art::Ptr<SimParticle>& p, float tt)
        : sx(p->endPosition().x())
        , sy(p->endPosition().y())
        , sz(p->endPosition().z())
        , spx(p->startMomentum().x())
        , spy(p->startMomentum().y())
        , spz(p->startMomentum().z())
        , spt(sqrt(p->startMomentum().x()*p->startMomentum().x()+p->startMomentum().y()*p->startMomentum().y()+p->startMomentum().z()*p->startMomentum().z()))
        , st(p->endGlobalTime())
        , stau(tt)
        , scode(p->stoppingCode())
        , sid(p->pdgId())
      {
       if(sqrt(p->startMomentum().x()*p->startMomentum().x()+p->startMomentum().y()*p->startMomentum().y()+p->startMomentum().z()*p->startMomentum().z()) >60){
          cout<<p->creationCode().name()<<endl;
        }
        if(!p->endDefined()) {
          throw cet::exception("BADINPUTS")
            <<"StoppedParticlesDumper: input SimParticle does not have end defined!\n";
        }
      }
    };

  double getCharge(PDGCode::type pdgId) {

    static GlobalConstantsHandle<ParticleDataTable> pdt;

    ParticleDataTable::maybe_ref info = pdt->particle(pdgId);

    if(!info.isValid()) {
      throw cet::exception("MISSINGINFO")<<"No valid PDG info for pdgId = "<<pdgId<<"\n";
    }

    return info.ref().charge();
  }

  //================================================================
  double getKineticEnergy(const StepPointMC& hit) {

    static GlobalConstantsHandle<ParticleDataTable> pdt;

    ParticleDataTable::maybe_ref info = pdt->particle(hit.simParticle()->pdgId());

    if(!info.isValid()) {
      throw cet::exception("MISSINGINFO")<<"No valid PDG info for hit = "<<hit<<"\n";
    }

    const double mass = info.ref().mass();
    return sqrt(hit.momentum().mag2() + std::pow(mass, 2)) - mass;
  }

    struct VDHit {
    float x;
    float y;
    float z;
    float time;

    float px;
    float py;
    float pz;
    float pmag;
    float ek;

    float charge;
    int   pdgId;
    unsigned particleId;
    unsigned volumeCopyNumber;

    VDHit() : x(std::numeric_limits<double>::quiet_NaN())
            , y(std::numeric_limits<double>::quiet_NaN())
            , z(std::numeric_limits<double>::quiet_NaN())

            , time(std::numeric_limits<double>::quiet_NaN())

            , px(std::numeric_limits<double>::quiet_NaN())
            , py(std::numeric_limits<double>::quiet_NaN())
            , pz(std::numeric_limits<double>::quiet_NaN())
            , pmag(std::numeric_limits<double>::quiet_NaN())
            , ek(std::numeric_limits<double>::quiet_NaN())

      , charge(std::numeric_limits<double>::quiet_NaN())
      , pdgId(0)
      , particleId(-1U)
      , volumeCopyNumber(-1U)
    {}

    //----------------------------------------------------------------
    VDHit(const StepPointMC& hit)
      : x(hit.position().x())
      , y(hit.position().y())
      , z(hit.position().z())

      , time(hit.time())

      , px(hit.momentum().x())
      , py(hit.momentum().y())
      , pz(hit.momentum().z())

      , pmag(hit.momentum().mag())
      , ek(getKineticEnergy(hit))

      , charge(getCharge(hit.simParticle()->pdgId()))

      , pdgId(hit.simParticle()->pdgId())
      , particleId(hit.simParticle()->id().asUint())
      , volumeCopyNumber(hit.volumeId())
    {}

  }; // struct VDHit




  class RPCAna : public art::EDAnalyzer {

     public:

     struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
      fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
      fhicl::Atom<art::InputTag> genTag{Name("GenParticleCollection"), Comment("gen particle info")};
      fhicl::Atom<art::InputTag> stepTag{Name("StepPointMCCollection"), Comment("step point info")};
      fhicl::Atom<art::InputTag> simTag{Name("SimParticleCollection"), Comment("sim particle info")};
      fhicl::Atom<bool> writeProperTime {Name("writeProperTime"),Comment("set to write tau."),false};
      fhicl::Sequence<int> decayOffPDGCodes {Name("decayOffPDGCodes"),Comment("pdg")};
      fhicl::Sequence<art::InputTag> hitCollections {Name("hitCollections"),Comment("A list of StepPointMCCollections")};

};
       typedef art::EDAnalyzer::Table<Config> Parameters;

       explicit RPCAna(const Parameters& conf);
       virtual ~RPCAna() {}


       virtual void beginJob();
       virtual void endJob();
	     virtual void analyze(const art::Event& e) override;

    private:
      Config _conf;
      int _diagLevel;
      int _nProcess;
      int _mcdiag;

      art::InputTag _genTag;
      art::InputTag _stepTag;
      art::InputTag _simTag; //VD
      bool writeProperTime_;
      std::vector<art::InputTag> hitColls_;
      std::vector<int> decayOffCodes_;
      const GenParticleCollection *_gencol;
      const SimParticleCollection *_simcol;
      const StepPointMCCollection *_stepcol;

      typedef std::vector<StepPointMCCollection> VspMC;
      VDHit hit1_, hit2_, hit3_;

      TTree* _Ntup;
      Int_t _evt, _run, _nSim;
      StopInfo data_;
      TH2F *pvz;
      TH2F *tvz;
      TH1F *weightedP;
      TH1F *weightedZ;
  
      void process(const art::Ptr<SimParticle>& p);
      bool is_leave(const SimParticle& p);
      bool findData(const art::Event& evt);

};


  RPCAna::RPCAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _mcdiag(conf().mcdiag()),
    _genTag(conf().genTag()),
    _stepTag(conf().stepTag()),
    _simTag(conf().simTag()),
    writeProperTime_(conf().writeProperTime())
  {
    if(writeProperTime_) {
      hitColls_ =  conf().hitCollections();
      decayOffCodes_ = conf().decayOffPDGCodes();
     }
    std::sort(decayOffCodes_.begin(), decayOffCodes_.end());
  }

  void RPCAna::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("RPCAna", "RPCAna");
    std::string simbranchDesc("sx/F:sy/F:sz/F:spx/F:spy/F:spz/F:spt/F:stime/F:scode/I:sid/I");
    if(writeProperTime_) {
      simbranchDesc += ":tauNormalized/F";
    } 

    static const char spbranchDesc1[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pmag/F:ek/F:charge/F:pdgId/I:particleId/i:volumeCopy/i";
    static const char spbranchDesc2[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pmag/F:ek/F:charge/F:pdgId/I:particleId/i:volumeCopy/i";
    static const char spbranchDesc3[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pmag/F:ek/F:charge/F:pdgId/I:particleId/i:volumeCopy/i";
    _Ntup->Branch("simParticles", &data_, simbranchDesc.c_str());
    _Ntup->Branch("vacuuas1",  &hit1_, spbranchDesc1);
    _Ntup->Branch("vacuuas2",  &hit2_, spbranchDesc2);
     _Ntup->Branch("virtualdetector",  &hit3_, spbranchDesc3);
     pvz = tfs->make<TH2F>("Total Mom v z [mm] ", "Total Mom v Z[mm]", 40, 5400, 6300, 100, 0, 100);
    tvz = tfs->make<TH2F>("Global Time v Z [mm] ", "Global Time v Z[mm]", 40, 5400,6300, 100, 100, 700);

  }


void RPCAna::analyze(const art::Event& event) {

  _evt = event.id().event();
  _run = event.run();

  if(!findData(event))
  throw cet::exception("RECO")<<"No data in  event"<< endl;

  ++_nProcess;
  if (_nProcess%1000000==0) std::cout<<"Processing event "<<_nProcess<<std::endl;

  _nSim = _simcol->size();

  const auto sc = event.getValidHandle<SimParticleCollection>(_simTag);
  for(const auto& p: *sc) {
    if(is_leave(p.second)) {
      art::Ptr<SimParticle> pp(sc, p.first.asUint());
      process(pp);
    }
  }
   VspMC spMCColls;
  for (const auto& iColl : hitColls_ ){
    auto spColl = event.getValidHandle<StepPointMCCollection>(iColl);
    spMCColls.push_back( *spColl );
  }
  for (const auto& i : spMCColls[0]){
    hit1_ = VDHit(i);
  }
  for (const auto& i : spMCColls[1]){
    hit2_ = VDHit(i);
  }
  for (const auto& i : spMCColls[2]){
    hit3_ = VDHit(i);
  }
  _nSim = _simcol->size();
  _Ntup->Fill();
}

void RPCAna::process(const art::Ptr<SimParticle>& p) {
  const float tau = -1;//writeProperTime_ ? SimParticleGetTau::calculate(p,spMCColls,decayOffCodes_) : -1;
  data_ = StopInfo(p, tau);
  //double weight = exp(-1*p->endProperTime()/tau);
  pvz->Fill(p->endPosition().z(), (sqrt(p->startMomentum().x()*p->startMomentum().x()+p->startMomentum().y()*p->startMomentum().y()+p->startMomentum().z()*p->startMomentum().z())));
  tvz->Fill(p->endPosition().z(),p->endGlobalTime());

}

bool RPCAna::is_leave(const SimParticle& p) {
    return p.daughters().empty();
}

bool RPCAna::findData(const art::Event& evt){
  _simcol=0;
  auto simpart = evt.getValidHandle<SimParticleCollection>(_simTag);
  _simcol = simpart.product();

	return _simcol!=0; 
}

 void RPCAna::endJob(){}

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCAna);
