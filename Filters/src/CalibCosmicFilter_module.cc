//
// This module will take a set of product and re-write them
// after updating the Simpartile and StepPointMC art::Ptr's.
// It would typically be run after FilterG4Out reduces and 
// rewrites the SimParticle and StepPointMC's.  This usually leaves
// other products pointed to the old collections.  This modules
// makes them point to the filtered collections
//   The algorithm assumes there are only unique Simparticle
// numbers (keys), even if there are multiple collecitons.
// This is the case for production-style multi-stage jobs, 
// but you can make a mess where this is not true.
//
// The input is the new Simparticle and StepPointMC collection 
// tags and the tags of the collections to read, repoint, and write.  
// Only certain products are implemented, since each product 
// requires custom code.
//

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SelectorBase.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Common/Ptr.h"
// the products that can be re-written
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"


namespace mu2e {

  // a set of SimParticle numbers
  typedef std::set<art::Ptr<SimParticle>::key_type> SPset;

  namespace {
    // Adapter for compressSimParticleCollection()
    class ParticleSelector {
    public:
      ParticleSelector(const SPset& m):spset(m) { }
      bool operator[]( cet::map_vector_key key ) const {
        return spset.find(key.asUint()) != spset.end();
      }

    private:
      SPset spset;
    };
  } // end namespace

  class CalibCosmicFilter : public art::EDFilter {

    typedef std::vector<art::InputTag> InputTags;
    typedef std::vector<std::string> VS;

    art::InputTag _SPtag; // the tracker StrawDigiMC collection
    art::InputTag _trkDMCtag; // the tracker StrawDigiMC collection
    art::InputTag _calCSStag;  // the CaloShowerStep coll
    art::InputTag _crvSPMCtag;  // the CaloShowerStep coll
    InputTags _copySPMC; // StepPointMCs to repoint
    VS _copySPMCnames; // output names
    InputTags _copyTraj; // Trajectories to repoint
    VS _copyTrajNames; // output names

  public:
    explicit CalibCosmicFilter(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  private:
    std::size_t _minTrkDigis;
    double _minCalEnergy;
    std::size_t _minCrvSteps;
    int _verbose;
  };

  //================================================================
  CalibCosmicFilter::CalibCosmicFilter(const fhicl::ParameterSet& pset)
  {

    _trkDMCtag = art::InputTag(pset.get<std::string>("strawDigiMCs","MakeSD"));
    _calCSStag = art::InputTag(pset.get<std::string>("caloShowerSteps",
					     "CaloShowerStepFromStepPt"));
    _crvSPMCtag = art::InputTag(pset.get<std::string>("crvSteps","g4run"));
    // will need to re-write these collections
    produces<mu2e::StrawDigiMCCollection>(_trkDMCtag.instance());
    produces<mu2e::CaloShowerStepCollection>(_calCSStag.instance());

    // new filtered SimParticle collection
    produces<mu2e::SimParticleCollection>(""); 
    
    VS slist; 

    // lists of StepPointMC collections to copy and re-point
    slist = pset.get<VS>("copyHits",VS());
    for(auto const& i: slist) _copySPMC.emplace_back(art::InputTag(i));
    _copySPMCnames = pset.get<VS>("copyHitsNames",VS());
    for(auto const& i:_copySPMCnames) produces<mu2e::StepPointMCCollection>(i);

    if (slist.size() != _copySPMCnames.size()) {
      throw cet::exception("BADCONFIG") << "the number of collections"
		<< " in copyHits and copyHitsNames must be the same\n";
    }

    // lists of MCTrajectory collections to copy and re-point
    slist = pset.get<VS>("copyTrajs",VS());
    for(auto const& i: slist) _copyTraj.emplace_back(art::InputTag(i));
    _copyTrajNames = pset.get<VS>("copyTrajsNames",VS());
    for(auto const& i:_copyTrajNames) produces<mu2e::MCTrajectoryCollection>(i);

    if (slist.size() != _copyTrajNames.size()) {
      throw cet::exception("BADCONFIG") << "the number of collections"
		<< " in copyTrajs and copyTrajsNames must be the same\n";
    }

    // filter requirements
    _minTrkDigis = pset.get<int>("minTrkDigis",0);
    _minCalEnergy = pset.get<int>("minCalEnergy",-1.0);
    _minCrvSteps = pset.get<int>("minCrvSteps",0);

    _verbose = pset.get<int>("verbose",0);

  }

  //================================================================
  bool CalibCosmicFilter::filter(art::Event& event) {

    bool passed = false;

    auto trkDMC = 
      event.getValidHandle<mu2e::StrawDigiMCCollection>(_trkDMCtag);
    auto calCSS  = 
      event.getValidHandle<CaloShowerStepCollection>(_calCSStag);
    auto crvSPMC = 
      event.getValidHandle<mu2e::StepPointMCCollection>(_crvSPMCtag);
    auto simPh  = 
      event.getValidHandle<SimParticleCollection>("g4run");
    auto const& simP = *simPh;

    // compute filtering result
    if(trkDMC->size()>=_minTrkDigis) passed = true;
    
    double energy = 0;
    for(auto const& s: *calCSS) energy += s.energyMC();
    if(energy>=_minCalEnergy) passed = true;

    if(crvSPMC->size()>=_minCrvSteps) passed = true;


    if(_verbose>3) std::cout 
	    << "CalibCosmicFilter results of filtering " 
	    << "  trk hits: " << std::setw(4) << trkDMC->size() 
	    << " cut:" << std::setw(4) << _minTrkDigis
	    << "  cal energy: " 
	    << std::setw(9) << std::setprecision(3) << energy 
	    << "  cut:" 
	    << std::setw(9) << std::setprecision(3) << _minCalEnergy
	    << "  crv steps: " << std::setw(4) << crvSPMC->size() 
	    << " cut:" << std::setw(4) << _minCrvSteps
            << "  passed: " << passed << std::endl;


    //************* make the list of saved SimParticles
    SPset SPsave;  // std::set of SP keys

    // save SimParticles from tracker StrawDigiMC
    for(auto const& sd : *trkDMC) { // loop over StrawDigiMC
      SPsave.insert(sd.stepPointMC(TrkTypes::cal)->simParticle().key());
      SPsave.insert(sd.stepPointMC(TrkTypes::hv )->simParticle().key());
      for(auto const& s: sd.stepPointMCs()) {
	SPsave.insert(s->simParticle().key());
      }
    }
    if(_verbose>5) std::cout 
	    << "CalibCosmicFilter SimParticles saved after tracker " 
            << SPsave.size() << std::endl;

    // save the SimParticles pointed to by the CaloShowerSteps
    for(auto const& s: *calCSS) SPsave.insert(s.simParticle().key());
    if(_verbose>5) std::cout 
	    << "CalibCosmicFilter SimParticles saved after cal " 
            << SPsave.size() << std::endl;

    // save the SimParticles pointed to by SPMcs
    for(auto const& tag: _copySPMC) { // loop over string of input tags
      auto h = event.getValidHandle<mu2e::StepPointMCCollection>(tag);
      for(auto const& i : *h) { // loop over the SPMC
	SPsave.insert(i.simParticle().key());
      }
    }
    if(_verbose>5) std::cout 
	    << "CalibCosmicFilter SimParticles saved after copyHits " 
            << SPsave.size() << std::endl;

    // save the SimParticles pointed to by MCTrajectory
    for(auto const& tag: _copyTraj) { // loop over string of input tags
      auto h = event.getValidHandle<mu2e::MCTrajectoryCollection>(tag);
      for(auto const& i : *h) { // loop over the SPMC
	SPsave.insert(i.first.key());
      }
    }
    if(_verbose>5) std::cout 
	    << "CalibCosmicFilter SimParticles saved after trajHits " 
            << SPsave.size() << std::endl;

    // save parents of saved SimParticles
    auto tempSave = SPsave;
    for(const auto& iid : tempSave) {
      auto ptr = simP.find(cet::map_vector_key(iid))->second.parent();
      while(ptr.isNonnull()) {
	SPsave.insert(ptr.key());
	ptr = ptr->parent();
      }
    }
    if(_verbose>5) std::cout 
	    << "CalibCosmicFilter SimParticles saved after parents " 
            << SPsave.size() << std::endl;

    //**********  write filtered output collection of SimParticles

    // these are needed to set the art::Ptrs
    auto SPpid = getProductID<SimParticleCollection>("");
    auto SPpg  = event.productGetter(SPpid);
    // make a set of new SimParticle art pointers for use later
    // save in a map that converts a key to an art::Ptr
    std::map<std::size_t,art::Ptr<SimParticle>> SPPtrmap;
    for(auto const& i : SPsave) { // list of saved SimParticle keys
      SPPtrmap[i] = art::Ptr<SimParticle>(SPpid,i,SPpg);
    }

    std::unique_ptr<mu2e::SimParticleCollection> 
                       outSPp(new mu2e::SimParticleCollection());
    auto& outSP = *outSPp;

    ParticleSelector sel(SPsave);    
    compressSimParticleCollection(SPpid, SPpg, simP, sel, outSP);

    if(_verbose>5) std::cout 
	    << "CalibCosmicFilter SimParticles written " 
            << outSP.size() << std::endl;

    // write out the new SimParticle collection
    event.put(std::move(outSPp), "");

    //***************** write the rest of the re-pointed collections

    //re-write StepPointMC collections
    // big translation map, will need it later
    std::map<art::Ptr<StepPointMC>,art::Ptr<StepPointMC>> SPMCmap;

    for(std::size_t i=0; i<_copySPMC.size(); i++) { // loop input tags
      auto pid = getProductID<StepPointMCCollection>(_copySPMCnames[i]);
      auto pg  = event.productGetter(pid);

      std::unique_ptr<mu2e::StepPointMCCollection> 
	                  outSPMCp(new mu2e::StepPointMCCollection());
      auto h = event.getValidHandle<mu2e::StepPointMCCollection>(_copySPMC[i]);
      size_t j = 0;
      for(auto const& hit : *h) { // loop over the SPMC
	outSPMCp->emplace_back(hit); // copy hit
	outSPMCp->back().simParticle() = 
	  SPPtrmap[hit.simParticle().key()]; // update Ptr
	// save Ptr translation
	SPMCmap[art::Ptr<StepPointMC>(h,j)] = art::Ptr<StepPointMC>(pid,j,pg); 
	j++;
      }

      if(_verbose>5) std::cout 
		       << "CalibCosmicFilter StepPointMCs written "
		       << _copySPMCnames[i] << " " 
		       << outSPMCp->size() << std::endl;

      event.put(std::move(outSPMCp), _copySPMCnames[i]);
    }

    //re-write MCTrajectory collections
    for(std::size_t i=0; i<_copyTraj.size(); i++) { // loop over input tags
      std::unique_ptr<mu2e::MCTrajectoryCollection> 
	                  outTrajp(new mu2e::MCTrajectoryCollection());
      auto& outTraj = *outTrajp;
      auto h = event.getValidHandle<mu2e::MCTrajectoryCollection>(_copyTraj[i]);
      for(auto const& trp : *h) { // loop over the SPMC
	// trajectory coll is a map of art ptr and traj
	outTraj[ SPPtrmap[trp.first.key()] ] = trp.second;
	outTraj[ SPPtrmap[trp.first.key()] ].sim() = SPPtrmap[trp.first.key()];
      }

      if(_verbose>5) std::cout 
		       << "CalibCosmicFilter MCTrajectory written "
		       << _copyTrajNames[i] << " " 
		       << outTrajp->size() << std::endl;

      event.put(std::move(outTrajp), _copyTrajNames[i]);
    }

    //re-write StrawDigiMC collection
    std::unique_ptr<mu2e::StrawDigiMCCollection> 
      outDMCp(new mu2e::StrawDigiMCCollection());
    auto& outDMC = *outDMCp;
    art::Ptr<StepPointMC> stepMC[2];
    std::vector<art::Ptr<StepPointMC> > stepMCs;
    for(auto d: *trkDMC) {
      stepMC[TrkTypes::cal] = SPMCmap[d.stepPointMC(TrkTypes::cal)];
      stepMC[TrkTypes::hv ] = SPMCmap[d.stepPointMC(TrkTypes::hv )];
      stepMCs.clear();
      for(auto i:d.stepPointMCs()) stepMCs.push_back(SPMCmap[i]);
      outDMC.push_back( mu2e::StrawDigiMC(d, stepMC, stepMCs) );
    } // loop over digis

    if(_verbose>5) std::cout 
		     << "CalibCosmicFilter StrawDigiMC written "
		     << _trkDMCtag.instance() << " " 
		     << outDMCp->size() << std::endl;

    event.put(std::move(outDMCp), _trkDMCtag.instance());


    //re-write CaloShowerStep collection - already have handle to old
    std::unique_ptr<mu2e::CaloShowerStepCollection> 
      outCSSp(new mu2e::CaloShowerStepCollection());
    auto& outCSS = *outCSSp;
    for(auto c: *calCSS) {
      outCSS.push_back(c);
      outCSS.back().setSimParticle(SPPtrmap[c.simParticle().key()]);
    } // loop over steps

    if(_verbose>5) std::cout 
		     << "CalibCosmicFilter CaloShowerStep written "
		     << _calCSStag.instance() << " " 
		     << outCSSp->size() << std::endl;

    event.put(std::move(outCSSp), _calCSStag.instance());


    return passed;

  }

  //================================================================
  void CalibCosmicFilter::endJob() {
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CalibCosmicFilter);
