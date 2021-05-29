//
// This module will check all the art::Ptr's in various products.
// Which types of products to check are hard-coded, and
// includes SimParticle, StepPointMC, and others.
// The code will throw if it finds a bad Ptr.
// You can tell it to skip certain products if you know they
// have a problem, or to skip the dereference check.
// verbose=1 by default which prints a summary for each product,
// or set verbose to 0 to run quietly.
//

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/ModuleMacros.h"
// the products that can be re-written
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CaloHitMC.hh"
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"


namespace mu2e {


  class PointerCheck : public art::EDAnalyzer {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<bool> skipDereference{Name("skipDereference"), 
	  Comment("skip dereferencing pointers check"),false};
      fhicl::Atom<int> verbose{Name("verbose"), 
	  Comment("verbose flag, 0 to 10"),1};

      fhicl::Sequence<art::InputTag> skipSimParticle{ 
	fhicl::Name("skipSimParticle"),
          fhicl::Comment("InputTag for collections to skip"),
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipSimParticlePtr{ 
	fhicl::Name("skipSimParticlePtr"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipStepPointMC{ 
	fhicl::Name("skipStepPointMC"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipMCTracjectory{ 
	fhicl::Name("skipMCTracjectory"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipStrawDigiMC{ 
	fhicl::Name("skipStrawDigiMC"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipCaloHitMC{ 
	fhicl::Name("skipCaloHitMC"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipCaloShowerStep{ 
	fhicl::Name("skipCaloShowerStep"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
      fhicl::Sequence<art::InputTag> skipCrvDigiMC{ 
	fhicl::Name("skipCrvDigiMC"),
          fhicl::Comment("InputTag for collections to skip"), 
	  std::vector<art::InputTag>()
	  };
    };

    // this line is required by art to allow the command line help print
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit PointerCheck(const Parameters& conf);
    void analyze(art::Event const&  event) override;
  private:

    typedef std::vector<int> pcounts;
    typedef std::vector<art::InputTag> InputTags;
    typedef std::vector<std::string> VS;

    bool _skipDereference; // do not do dereference check
    int _verbose; // 0 = none, 1 = print as you go

    InputTags _SPtags; // SimParticles to skip
    InputTags _SPPtrtags; // SimParticlePtr to skip
    InputTags _SPMCtags; // StepointMC to skip
    InputTags _trajtags; // MCTrajectories to skip
    InputTags _trkDMCtags; // StrawDigiMC to skip
    InputTags _calDMCtags;  // CaloHitMC to skip
    InputTags _calCSStags;  // CaloShowerStep to skip
    InputTags _crvDMCtags;  // CaloShowerStep to skip

    void printProvenance(art::Provenance const& p);
    bool excludedCollection(art::Provenance const& p, InputTags const& tags);

    bool checkSimParticle(SimParticleCollection const& coll);
    bool checkStepPointMC(StepPointMCCollection const& coll);
    bool checkMCTrajectory(MCTrajectoryCollection const& coll);
    bool checkSimParticlePtr(SimParticlePtrCollection const& coll);
    bool checkStrawDigiMC(StrawDigiMCCollection const& coll);
    bool checkCaloHitMC(CaloHitMCCollection const& coll);
    bool checkCaloShowerStep(CaloShowerStepCollection const& coll);
    bool checkCrvDigiMC(CrvDigiMCCollection const& coll);

  };

  //================================================================
  PointerCheck::PointerCheck(const Parameters& conf):
    art::EDAnalyzer(conf),
    _skipDereference(conf().skipDereference()),
    _verbose(conf().verbose()),
    _SPtags(conf().skipSimParticle()),
    _SPPtrtags(conf().skipSimParticlePtr()),
    _SPMCtags(conf().skipStepPointMC()),
    _trajtags(conf().skipMCTracjectory()),
    _trkDMCtags(conf().skipStrawDigiMC()),
    _calDMCtags(conf().skipCaloHitMC()),
    _calCSStags(conf().skipCaloShowerStep()),
    _crvDMCtags(conf().skipCrvDigiMC())
  {
  }

  //================================================================
  void PointerCheck::analyze(art::Event const& event) {

    std::vector< art::Handle<SimParticleCollection> > vah_sp = event.getMany<SimParticleCollection>();
    for (auto const & ah : vah_sp) {
      if(!excludedCollection(*ah.provenance(),_SPtags)) {
	printProvenance(*ah.provenance());
	checkSimParticle(*ah);
      } // endif excluded
    } // loop over handles

    std::vector< art::Handle<StepPointMCCollection> > vah_st = event.getMany<StepPointMCCollection>();
    for (auto const & ah : vah_st) {
      if(!excludedCollection(*ah.provenance(),_SPMCtags)) {
	printProvenance(*ah.provenance());
	checkStepPointMC(*ah);
      } // endif excluded
    } // loop over handles

    std::vector< art::Handle<MCTrajectoryCollection> > vah_tr = event.getMany<MCTrajectoryCollection>();
    for (auto const & ah : vah_tr) {
      if(!excludedCollection(*ah.provenance(),_trajtags)) {
	printProvenance(*ah.provenance());
	checkMCTrajectory(*ah);
      } // endif excluded
    } // loop over handles

    std::vector< art::Handle<SimParticlePtrCollection> > vah_pp = event.getMany<SimParticlePtrCollection>();
    for (auto const & ah : vah_pp) {
      if(!excludedCollection(*ah.provenance(),_trajtags)) {
	printProvenance(*ah.provenance());
	checkSimParticlePtr(*ah);
      } // endif excluded
    } // loop over handles

    std::vector< art::Handle<StrawDigiMCCollection> > vah_sd = event.getMany<StrawDigiMCCollection>();
    for (auto const & ah : vah_sd) {
      if(!excludedCollection(*ah.provenance(),_trkDMCtags)) {
	printProvenance(*ah.provenance());
	checkStrawDigiMC(*ah);
      } // endif excluded
    } // loop over handles

    std::vector< art::Handle<CaloHitMCCollection> > vah_cd = event.getMany<CaloHitMCCollection>();
    for (auto const & ah : vah_cd) {
      if(!excludedCollection(*ah.provenance(),_calDMCtags)) {
	printProvenance(*ah.provenance());
	checkCaloHitMC(*ah);
      } // endif excluded
    } // loop over handles
    std::cout << "Startgin shower steps" << std::endl;
    std::vector< art::Handle<CaloShowerStepCollection> > vah_cs = event.getMany<CaloShowerStepCollection>();
    for (auto const & ah : vah_cs) {
      if(!excludedCollection(*ah.provenance(),_calCSStags)) {
	printProvenance(*ah.provenance());
	checkCaloShowerStep(*ah);
      } // endif excluded
    } // loop over handles

    std::cout << "Starting crvMC" << std::endl;
    std::vector< art::Handle<CrvDigiMCCollection> > vah_cm = event.getMany<CrvDigiMCCollection>();
    for (auto const & ah : vah_cm) {
      if(!excludedCollection(*ah.provenance(),_crvDMCtags)) {
	printProvenance(*ah.provenance());
	checkCrvDigiMC(*ah);
      } // endif excluded
    } // loop over handles

  
  }

  void PointerCheck::printProvenance(art::Provenance const& p) {
    if(_verbose<1) return;
    std::cout 
      << " " << p.friendlyClassName() 
      << ":" << p.moduleLabel() 
      << ":" << p.productInstanceName() 
      << ":" << p.processName() 
      << std::endl;
    return;
  }

  bool PointerCheck::excludedCollection(art::Provenance const& p, 
					InputTags const& tags) {

    // a tag has four fields, Classname already determined
    // treat empty fields as wildcards that match
    bool exclude = true;
    for(auto const& tag : tags) {
      if( !tag.label().empty() && 
	  tag.label() != p.moduleLabel() ) exclude = false;
      if( !tag.instance().empty() && 
	  tag.instance() != p.productInstanceName() ) exclude = false;
      if( !tag.process().empty() &&
	  tag.process() != p.processName() ) exclude = false;
      if(exclude) return true;
    }
    return false;
  }

  bool PointerCheck::checkSimParticle(SimParticleCollection const& coll) {

    int n,nn,na,ni,ns,ne;
    n=nn=na=ni=ns=ne=0;
    int gn,ga,gi; // genparts
    gn=ga=gi=0;
    int nd,dn,da,di; // daughters
    nd=dn=da=di=0;
    // check internal and parents
    for(auto const& x: coll) { // loop over the collection
      n++;
      auto const& i = x.first;
      auto const& s = x.second;
      // check that the map index equals the SP internal index
      if(i==s.id()) {
	ne++;
      } else {
	throw cet::exception("BadSimParticleKey")
	  << "SimParticle map_vector index does not match id()";
      }
      // check parent pointers
      auto const& p = s.parent();
      if(p.isNonnull()) {
	nn++;
	if(p.isAvailable()) {
	  na++;
	  if(!_skipDereference && p.get()) {
	    ni++;
	  }
	}
      } else { // null
	if(s.isPrimary()) ns++;
      }
      
      auto const& g = s.genParticle();
      if(g.isNonnull()) {
	gn++;
	if(g.isAvailable()) {
	  ga++;
	  if(!_skipDereference && g.get()) {
	    gi++;
	  }
	}
      }
      
      // check daughters
      for(auto const& d: s.daughters()) {
	nd++;
	if(d.isNonnull()) {
	  dn++;
	  if(d.isAvailable()) {
	    da++;
	    if(!_skipDereference && d.get()) {
	      di++;
	    }
	  }
	}
      } // loop over daughters
      
    } // loop over SP in coll
    
    bool rc;
    rc = (ne==n && nn+ns==n && ni==nn && gi==gn && di==nd);

    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " SimParticles checked" << std::endl;
    std::cout << std::setw(8) << ne << " have matching internal index" << std::endl;
    std::cout << std::setw(8) << ns << " parents null and isPrimary()=true" << std::endl;
    std::cout << std::setw(8) << nn << " parents nonnull" << std::endl;
    std::cout << std::setw(8) << na << " parents available" << std::endl;
    std::cout << std::setw(8) << ni << " parents dereferenced" << std::endl;
    std::cout << std::setw(8) << gn << " gen particles nonnull" << std::endl;
    std::cout << std::setw(8) << ga << " gen particles available" << std::endl;
    std::cout << std::setw(8) << gi << " gen particles dereferenced" << std::endl;
    std::cout << std::setw(8) << nd << " daughters in arrays" << std::endl;
    std::cout << std::setw(8) << dn << " daughters nonnull" << std::endl;
    std::cout << std::setw(8) << da << " daughters available" << std::endl;
    std::cout << std::setw(8) << di << " daughters dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in SimParticle coll";

    return rc;
  }


  bool PointerCheck::checkStepPointMC(StepPointMCCollection const& coll) {

    int n,nn,na,ni;
    n=nn=na=ni=0;
    for(auto const& s: coll) { // loop over the collection
      n++;
      // check parent pointers
      auto const& p = s.simParticle();
      if(p.isNonnull()) {
	nn++;
	if(p.isAvailable()) {
	  na++;
	  if(!_skipDereference && p.get()) {
	    ni++;
	  }
	}
      }
    } // loop over SP in coll
    
    bool rc = (ni==n);
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " StepPointMC checked" << std::endl;
    std::cout << std::setw(8) << nn << " parents nonnull" << std::endl;
    std::cout << std::setw(8) << na << " parents available" << std::endl;
    std::cout << std::setw(8) << ni << " parents dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in StepPointMC coll";

    return rc;
  }

  bool PointerCheck::checkMCTrajectory(MCTrajectoryCollection const& coll) {

    int n,ne,nn,na,ni;
    n=ne=nn=na=ni=0;
    for(auto const& p: coll) { // loop over the collection
      auto const& a = p.first; //art::Ptr to SimParticle
      auto const& t = p.second; // tracjetory, with embedded ptr
      auto const& b = t.sim();
      n++;
      // check parent pointers
      if(a==b) {
	ne++;
      } else {
	throw cet::exception("BadSimParticlePtr")
	  << "MCTrajectory external and internal SimParticle pointers do not agree";
      }
      if(b.isNonnull()) {
	nn++;
	if(b.isAvailable()) {
	  na++;
	  if(!_skipDereference && b.get()) {
	    ni++;
	  }
	}
      }
    } // loop over SP in coll
    
    bool rc = (ni==n);
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " MCTrajectory checked" << std::endl;
    std::cout << std::setw(8) << ne << " external and internal Ptrs match" << std::endl;
    std::cout << std::setw(8) << nn << " parents nonnull" << std::endl;
    std::cout << std::setw(8) << na << " parents available" << std::endl;
    std::cout << std::setw(8) << ni << " parents dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in MCTrajectory coll";

    return rc;
  }

  bool PointerCheck::checkSimParticlePtr(SimParticlePtrCollection const& coll) {

    int n,nn,na,ni;
    n=nn=na=ni=0;
    for(auto const& p: coll) { // loop over the collection, a vector of Ptr
      n++;
      if(p.isNonnull()) {
	nn++;
	if(p.isAvailable()) {
	  na++;
	  if(!_skipDereference && p.get()) {
	    ni++;
	  }
	}
      }
    } // loop over SP in coll
    
    bool rc = (ni==n);
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " SimParticlePtr checked" << std::endl;
    std::cout << std::setw(8) << nn << " Ptrs nonnull" << std::endl;
    std::cout << std::setw(8) << na << " Ptrs available" << std::endl;
    std::cout << std::setw(8) << ni << " Ptrs dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in SimParticlePtr coll";

    return rc;
  }

  bool PointerCheck::checkStrawDigiMC(StrawDigiMCCollection const& coll) {

    int n,np,nn,na,ni;
    n=np=nn=na=ni=0;
    std::vector<art::Ptr<StrawGasStep> > ptrs;
    for(auto const& d: coll) { // loop over the collection
      n++;
      // assemble all the pointer in the object
      ptrs.push_back(d.strawGasStep(StrawEnd::cal));
      ptrs.push_back(d.strawGasStep(StrawEnd::hv ));
      // check them
      for(auto const& p: ptrs) {
	np++;
	if(p.isNonnull()) {
	  nn++;
	  if(p.isAvailable()) {
	    na++;
	    if(!_skipDereference && p.get()) {
	      ni++;
	    }
	  }
	}
      }
    } // loop over SDMC in coll
    
    bool rc = (ni==np);
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " StrawDigiMC checked" << std::endl;
    std::cout << std::setw(8) << np << " Ptrs checked" << std::endl;
    std::cout << std::setw(8) << nn << " Ptrs nonnull" << std::endl;
    std::cout << std::setw(8) << na << " Ptrs available" << std::endl;
    std::cout << std::setw(8) << ni << " Ptrs dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in StrawDigiMC coll";

    return rc;
  }

  bool PointerCheck::checkCaloHitMC(CaloHitMCCollection const& coll) {

    int n,n2,nn,na,ni;
    n=n2=nn=na=ni=0;
    for(auto const& s: coll) { // loop over the collection
      n++;
      for(unsigned i=0; i<s.nParticles(); i++) {
	n2++;
	auto const& p = s.energyDeposit(i).sim();
	if(p.isNonnull()) {
	  nn++;
	  if(p.isAvailable()) {
	    na++;
	    if(!_skipDereference && p.get()) {
	      ni++;
	    }
	  }
	}
      } // loop over Simparticles in a digi
    } // loop over SP in coll
    
    bool rc = (ni==n);
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " CaloHitMC checked" << std::endl;
    std::cout << std::setw(8) << n2 << " SimParticle Ptrs checked" << std::endl;
    std::cout << std::setw(8) << nn << " Ptrs nonnull" << std::endl;
    std::cout << std::setw(8) << na << " Ptrs available" << std::endl;
    std::cout << std::setw(8) << ni << " Ptrs dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in CaloHitMC coll";

    return rc;
  }

  bool PointerCheck::checkCaloShowerStep(CaloShowerStepCollection const& coll) {

    int n,nn,na,ni;
    n=nn=na=ni=0;
    for(auto const& s: coll) { // loop over the collection
      n++;
      // check parent pointers
      auto const& p = s.simParticle();
      if(p.isNonnull()) {
	nn++;
	if(p.isAvailable()) {
	  na++;
	  if(!_skipDereference && p.get()) {
	    ni++;
	  }
	}
      }
    } // loop over SS in coll
    
    bool rc = (ni==n);
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " CaloShowerStep checked" << std::endl;
    std::cout << std::setw(8) << nn << " parents nonnull" << std::endl;
    std::cout << std::setw(8) << na << " parents available" << std::endl;
    std::cout << std::setw(8) << ni << " parents dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in CaloShowerStep coll";

    return rc;
  }

  bool PointerCheck::checkCrvDigiMC(CrvDigiMCCollection const& coll) {

    int n,nn,na,ni,np;
    n=nn=na=ni=0;
    for(auto const& s: coll) { // loop over the collection
      n++;
      // check parent pointers
      auto const& p = s.GetSimParticle();
      if(p.isNonnull()) {
	nn++;
	if(p.isAvailable()) {
	  na++;
	  if(!_skipDereference && p.get()) {
	    ni++;
	  }
	}
      }
    } // loop over SS in coll
    
    bool rc = (n<10 || (n>10 && ni>=n/2));
    if(_verbose>0) {
    // report
      std::cout << "Note: about 5-10% of the digis are from noise and have empty SimParticle" << std::endl;
    std::cout << std::setw(8) << n  << " CrvDigiMC SimParticle Ptr checked" << std::endl;
    std::cout << std::setw(8) << nn << " parents nonnull" << std::endl;
    std::cout << std::setw(8) << na << " parents available" << std::endl;
    std::cout << std::setw(8) << ni << " parents dereferenced" << std::endl;
    }

    if(!rc) throw cet::exception("BadArtPtr") << " in CrvDigiMC coll";

    n=np=nn=na=ni=0;
    std::vector<art::Ptr<CrvStep> > ptrs;
    for(auto const& d: coll) { // loop over the collection
      n++;
      // assemble all the pointer in the object
      ptrs = d.GetCrvSteps();

      // check them
      for(auto const& p: ptrs) {
	np++;
	if(p.isNonnull()) {
	  nn++;
	  if(p.isAvailable()) {
	    na++;
	    if(!_skipDereference && p.get()) {
	      ni++;
	    }
	  }
	}
      }
    } // loop over SDMC in coll
    
    rc = rc && (n<10 || (n>10 && ni>=n/2));
    if(_verbose<1 && !rc) return rc;
    // report

    std::cout << std::setw(8) << n  << " CrvDigiMC CrvStep Ptr checked" << std::endl;
    std::cout << std::setw(8) << np << " Ptrs checked" << std::endl;
    std::cout << std::setw(8) << nn << " Ptrs nonnull" << std::endl;
    std::cout << std::setw(8) << na << " Ptrs available" << std::endl;
    std::cout << std::setw(8) << ni << " Ptrs dereferenced" << std::endl;

    if(!rc) throw cet::exception("BadArtPtr") << " in CrvDigiMC coll";

    return rc;
  }



} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PointerCheck)
