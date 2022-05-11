// Write into an ntuple information about time, position, and
// (optionally) proper time of SimParticle end points.  There are two
// ways to specify the set of particles to process:
//
// 1) dumpSimParticleLeaves=false (default), inputCollection is a
//    SimParticlePtrCollection that explicitly lists what to dump.
//
// 2) dumpSimParticleLeaves=true, inputCollection is a SimParticle
//    collection.  The leaves of the SimParticle tree will be dumped.
//
// Andrei Gaponenko, 2013, 2015

#include <string>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"

#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

#include <algorithm>
#include <iterator>

namespace mu2e {

  namespace {

    typedef std::vector<StepPointMCCollection> VspMC;

    // This should be minimal info, we'll want to load this
    // in memory in consumer jobs.  This is NOT an analysis ntuple!
    struct StopInfo {
      float x;
      float y;
      float z;
      float px;
      float py;
      float pz;
      float pt;
      float t;
      float tau; // proper time, for stopped pion weights
      ProcessCode code;
      PDGCode::type id;
      float x0;
      float y0;
      float z0;
      float productionP;
      float angle;
      StopInfo() : x(), y(), z(), px(), py(), pz(), pt(), t(), tau(), code(), id(), x0(), y0(), z0(), productionP(), angle() {}

      StopInfo(const art::Ptr<SimParticle>& p, const VspMC& spMCcolls, float tt)
        : x(p->endPosition().x())
        , y(p->endPosition().y())
        , z(p->endPosition().z())
        , px(p->startMomentum().x())
        , py(p->startMomentum().y())
        , pz(p->startMomentum().z())
        , pt(sqrt(p->startMomentum().x()*p->startMomentum().x()+p->startMomentum().y()*p->startMomentum().y()+p->startMomentum().z()*p->startMomentum().z()))
        , t(p->endGlobalTime())
        , tau(tt)
        , code(p->stoppingCode())
        , id(p->pdgId())
        , x0(p->startPosition().x())
        , y0(p->startPosition().y())
        , z0(p->startPosition().z())
      {
        //std::cout<<sqrt(p->startMomentum().x()*p->startMomentum().x()+p->startMomentum().y()*p->startMomentum().y()+p->startMomentum().z()*p->startMomentum().z())<<"   "<<p->startPosition().x()<<" "<<p->startPosition().y()<<" "<<p->startPosition().z()<<std::endl;
        std::cout<<"PDG "<<p->pdgId()<<std::endl;
        /*int n=0;
        double pdg = 211;
        art::Ptr<SimParticle>  sp =p->parent();
			  while(abs(pdg) == 211){
			    art::Ptr<SimParticle> temppart = sp->parent();
			    pdg = temppart->pdgId();
			    if(abs(pdg) == 211){ 
			      sp = temppart;
			    }
			    std::cout<<n<<" "<<pdg<<std::endl;
			    n++;
			  }
			  
			  productionP = sqrt(sp->startMomentum().x()*sp->startMomentum().x()+sp->startMomentum().y()*sp->startMomentum().y()+sp->startMomentum().z()*sp->startMomentum().z());
			  std::cout<<productionP<<std::endl;*/

        if(!p->endDefined()) {
          throw cet::exception("BADINPUTS")
            <<"myStoppedParticleDumper: input SimParticle does not have end defined!\n";
        }
      }
    };

  }// namespace

  //================================================================
  class myStoppedParticleDumper : public art::EDAnalyzer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<bool> dumpSimParticleLeaves {
        Name("dumpSimParticleLeaves"),
          Comment("The mode selector.  Set it to true to dump info on all the leaves\n"
                  "of a SimParticleCollection.   If it is set to false, the input should\n"
                  "be of the SimParticlePtrCollection type, and its full content will be included.\n"
                  ),
          false
          };

      fhicl::Atom<art::InputTag> inputCollection {
        Name("inputCollection"),
          Comment("InputTag of the collection to use.")
          };

      fhicl::Atom<bool> writeProperTime {
        Name("writeProperTime"),
          Comment("Compute and write out normalized proper time tau of particles.\n"
                  "This quantity is defined in such a way that exp(-tau) is the survival\n"
                  "probability of the particle.  It is used in cases when the decay\n"
                  "process is disabled during the simulation."
                  ),
          false
          };

      fhicl::Sequence<int> decayOffPDGCodes {
        Name("decayOffPDGCodes"),
          Comment("A list of PDG IDs of particles that had their decay process turned off during\n"
                  "the simulation. This must be specified if and only if writeProperTime is requested."),
          [this](){ return writeProperTime(); }
      };

      fhicl::Sequence<art::InputTag> hitCollections {
        Name("hitCollections"),
          Comment("A list of StepPointMCCollection-s via which different stages of simulation\n"
                  "are connected.  Must be specified if and only if writeProperTime is requested."),
          [this](){ return writeProperTime(); }
      };

    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit myStoppedParticleDumper(const Parameters& conf);

    void beginJob() override;
    void analyze(const art::Event& evt) override;
  private:
    bool dumpSimParticleLeaves_;
    art::InputTag input_;
    bool writeProperTime_;
    std::vector<art::InputTag> hitColls_;

    std::vector<int> decayOffCodes_;

    TTree *nt_;
    StopInfo data_;
    Float_t pod, ang,pdg;
    TH2F *pvz;
    TH2F *tvz;
    TH1F *weightedP;
    TH1F *weightedZ;
    bool is_leave(const SimParticle& p);
    void process(const art::Ptr<SimParticle>& p, const VspMC& spMCcolls);
  };

  //================================================================
  myStoppedParticleDumper::myStoppedParticleDumper(const Parameters& conf) :
    art::EDAnalyzer(conf),
    dumpSimParticleLeaves_(conf().dumpSimParticleLeaves()),
    input_(conf().inputCollection()),
    writeProperTime_(conf().writeProperTime()),
    nt_()
  {
    if(writeProperTime_) {
      hitColls_ =  conf().hitCollections();
      decayOffCodes_ = conf().decayOffPDGCodes();

      // must sort to use binary_search in SimParticleGetTau
      std::sort(decayOffCodes_.begin(), decayOffCodes_.end());
    }
  }

  //================================================================
  void myStoppedParticleDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    std::string branchDesc("x/F:y/F:z/F:px/F:py/F:pz/F:pt/F:time/F:code/I:id/I:productionP/F:x0/F:y0/F:z0/F");
    if(writeProperTime_) {
      branchDesc += ":tauNormalized/F";
    }
    nt_ = tfs->make<TTree>( "Pions", "Particles ntuple");
    nt_->Branch("Pions", &data_, branchDesc.c_str());
    nt_->Branch("Production", &pod, "pod/F");
    nt_->Branch("ProductionAngle", &ang, "ang/F");
    nt_->Branch("PDG", &pdg, "pdg/F");
    pvz = tfs->make<TH2F>("Total Mom v z [mm] ", "Total Mom v Z[mm]", 40, 5400, 6300, 100, 0, 100);
    tvz = tfs->make<TH2F>("Global Time v Z [mm] ", "Global Time v Z[mm]", 40, 5400,6300, 100, 100, 700);
		weightedZ =  tfs->make<TH1F>("Weighted Z", "Weighted MomZ", 40, 5400, 6300);
    //weightedP =  tfs->make<TH1F>("Weighted Mom", "Weighted Mom", 100, 0, 100);
    weightedP =  tfs->make<TH1F>("production", "production", 50, 0, 200);
  }

  //================================================================
  void myStoppedParticleDumper::analyze(const art::Event& event) {

    VspMC spMCColls;

    for ( const auto& iColl : hitColls_ ){
      auto spColl = event.getValidHandle<StepPointMCCollection>(iColl);
      spMCColls.push_back( *spColl );
    }


    if(dumpSimParticleLeaves_) {
      auto ih = event.getValidHandle<SimParticleCollection>(input_);
      for(const auto& p : *ih) {
        if(is_leave(p.second)) {
          art::Ptr<SimParticle> pp(ih, p.first.asUint());
          process(pp, spMCColls);
        }
      }
    }
    else {
      auto ih = event.getValidHandle<SimParticlePtrCollection>(input_);
      for(const auto& p : *ih) {
        process(p, spMCColls);
      }
    }

  }

  //================================================================
  void myStoppedParticleDumper::process(const art::Ptr<SimParticle>& p, const VspMC& spMCColls) {
    const float tau = writeProperTime_ ? SimParticleGetTau::calculate(p,spMCColls,decayOffCodes_) : -1;
     pdg = static_cast<double>(p->pdgId());
     if(abs(p->pdgId()) == 13)data_ = StopInfo(p, spMCColls, tau);
     //----pions
       int n=0;
        double pdg = 211;
        art::Ptr<SimParticle>  sp =p->parent();
			  while(abs(pdg) == 211){
			    art::Ptr<SimParticle> temppart = sp->parent();
			    pdg = temppart->pdgId();
			    if(abs(pdg) == 211){ 
			      sp = temppart;
			    }
			    std::cout<<n<<" "<<pdg<<std::endl;
			    n++;
			  }

			 double productionP = sqrt(sp->startMomentum().x()*sp->startMomentum().x()+sp->startMomentum().y()*sp->startMomentum().y()+sp->startMomentum().z()*sp->startMomentum().z());
			 ang = sp->startMomentum().z()/productionP;
			 pod = productionP;
			 
			// -----pions
     double weight = exp(-1*p->endProperTime()/tau);
     pvz->Fill(p->endPosition().z(), (sqrt(p->startMomentum().x()*p->startMomentum().x()+p->startMomentum().y()*p->startMomentum().y()+p->startMomentum().z()*p->startMomentum().z())));
    tvz->Fill(p->endPosition().z(),p->endGlobalTime());
    //weightedP->Fill(productionP);
    weightedZ->Fill(p->endPosition().z(),weight);
    nt_->Fill();


  }

  //================================================================
  bool myStoppedParticleDumper::is_leave(const SimParticle& p) {
    return p.daughters().empty();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::myStoppedParticleDumper);
