//
//  This module creates the CrvStep objects used in downstream digi simulation, using the
//  G4 StepPointMCs
//
//  Original author: Ralf Ehrlich (based on code by David Brown (LBNL), Krzysztof Genser)
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "HepPDT/ParticleData.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include <utility>
#include <algorithm>
// root
#include "TTree.h"

using namespace std;
using HepPDT::ParticleData;
namespace mu2e 
{

  class CrvStepsFromStepPointMCs : public art::EDProducer 
  {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config 
      {
	fhicl::Atom<int> debug{ Name("debugLevel"),
	  Comment("Debug Level"), 0};
	fhicl::Atom<int> diag{ Name("diagLevel"),
	  Comment("Diag Level"), 0};
	fhicl::Atom<unsigned> csize{ Name("OutputCollectionSize"),
	  Comment("Estimated size of output collection"), 2000};
	fhicl::Atom<string> stepPointsInstance { Name("stepPointsInstance"),
	  Comment("Crv StepPointMC product instance name"),"CRV"};
	fhicl::Sequence<string> stepPointsModuleLabels { Name("stepPointsModuleLabels"),
	  Comment("Crv StepPointMC module label names")};
	fhicl::Atom<unsigned> startSize { Name("StartSize"),
	  Comment("Starting size for crvIndex-particle vector"),4};
	fhicl::Atom<bool> removeNeutralParticles{ Name("removeNeutralParticles"),
	  Comment("Removes steps of neutral particles"),true};
	fhicl::Atom<bool> useTotalEDep{ Name("useTotalEDep"),
	  Comment("Use total energy deposition instead. Visible is a new default"),false};
	fhicl::Atom<bool> noPostPositionAvailable{ Name("noPostPositionAvailable"),
	  Comment("No postposition is available for CRY3/4 samples"),false};

      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit CrvStepsFromStepPointMCs(const Parameters& conf);

    private:
      void beginJob() override;
      void beginRun(art::Run& run) override;
      void produce(art::Event& e) override;
      typedef pair<CRSScintillatorBarIndex,cet::map_vector_key> BarIndexTrackPair; // key for pair of Crv barIndex, SimParticle
      typedef art::Ptr<StepPointMC> SPMCPtr;
      typedef vector<SPMCPtr> SPMCPtrV;
      typedef map<BarIndexTrackPair, SPMCPtrV> StepPointMap; // steps by Crv barIndex and SimParticle
      typedef art::Handle<StepPointMCCollection> SPMCCH;
      typedef vector< SPMCCH > SPMCCHV;
      void fillMap(SPMCCH const& spmcch, StepPointMap& stepPointMap, GlobalConstantsHandle<ParticleDataTable> &pdt);
      void fillSteps(SPMCPtrV const& spmcptrv, CRSScintillatorBarIndex const& barIndex,
	             ParticleData const* pdata, cet::map_vector_key trackId,
                     unique_ptr<CrvStepCollection> &crvSteps, vector<pair<size_t,size_t> > &spmcIndices);
      void fillStepDiag(CrvStep const& crvStep, SPMCPtrV const& spmcptrv, size_t firstIndex, size_t lastIndex);
      static bool compareStepTimes(SPMCPtr a, SPMCPtr b) {return(a->time()<b->time());}

      int      _debug, _diag;
      bool     _removeNeutralParticles;
      bool     _useTotalEDep;
      bool     _noPostPositionAvailable;
      unsigned _csize, _ssize;
      // StepPointMC selector
      // This selector will select only data products with the given instance name.
      // optionally exclude modules: this is a fix
      art::Selector   _selector;
      bool            _firstEvent;
      // diagnostic histograms
      TTree*          _csdiag;
      Float_t         _cslen, _csdist, _csedep, _csPartP;
      vector<Float_t> _slen, _sdist, _sedep, _sPartP;
      Int_t           _csPartPDG;
      vector<Int_t>   _sPartPDG;
  };

  CrvStepsFromStepPointMCs::CrvStepsFromStepPointMCs(const Parameters& config )  : 
    art::EDProducer{config},
    _debug(config().debug()),
    _diag(config().diag()),
    _removeNeutralParticles(config().removeNeutralParticles()),
    _useTotalEDep(config().useTotalEDep()),
    _noPostPositionAvailable(config().noPostPositionAvailable()),
    _csize(config().csize()),
    _ssize(config().startSize()),
    _selector(art::ProductInstanceNameSelector(config().stepPointsInstance())),
    _firstEvent(true)
  {
    consumesMany<StepPointMCCollection>();
    produces <CrvStepCollection>();
    const std::vector<std::string> &labels = config().stepPointsModuleLabels();
    if(labels.size()==0)
    {
      throw cet::exception("SIM")<<"mu2e::CrvStepsFromStepPointMCs: No StepPointMC Module Labels specificied" << endl;
    }
    if(std::find(labels.begin(),labels.end(),"*")==labels.end())
    {
      for(size_t i=0; i<labels.size(); ++i)
      {
        art::Selector tmpSelector{art::ProductInstanceNameSelector(config().stepPointsInstance()) &&
                                  art::ModuleLabelSelector(labels.at(i))};
        if(i==0) _selector=tmpSelector;
        else _selector=art::Selector{_selector || tmpSelector};
      }
    }
  }

  void CrvStepsFromStepPointMCs::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    if(_diag > 1)
    {
      _csdiag=tfs->make<TTree>("csdiag","CrvStep diagnostics");
      _csdiag->Branch("cslen",&_cslen,"cslen/F");
      _csdiag->Branch("csdist",&_csdist,"csdist/F");
      _csdiag->Branch("csedep",&_csedep,"csedep/F");
      _csdiag->Branch("csPartP",&_csPartP,"csPartP/F");
      _csdiag->Branch("csPartPDG",&_csPartPDG,"csPartPDG/I");
      _csdiag->Branch("slen",&_slen);
      _csdiag->Branch("sdist",&_sdist);
      _csdiag->Branch("sedep",&_sedep);
      _csdiag->Branch("sPartP",&_sPartP);
      _csdiag->Branch("sPartPDG",&_sPartPDG);
    }
  }

  void CrvStepsFromStepPointMCs::beginRun( art::Run& run)
  {
  }

  void CrvStepsFromStepPointMCs::produce(art::Event& event) 
  {
    GlobalConstantsHandle<ParticleDataTable> pdt;

    //CrvSteps
    unique_ptr<CrvStepCollection> crvSteps(new CrvStepCollection);
    crvSteps->reserve(_csize);
 
    // Get all of the tracker StepPointMC collections from the event:
    // This selector will select only data products with the given instance name.
    SPMCCHV stepsHandles = event.getMany<StepPointMCCollection>(_selector);

    // Informational message on the first event.
    if(_firstEvent)
    {
      mf::LogInfo log("COSMIC");
      log << "mu2e::CrvStepsFromStepPointMCs will use StepPointMCs from: \n";
      for(SPMCCHV::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i!=e; ++i)
      {
	art::Provenance const& prov(*(i->provenance()));
	log  << "   " << prov.branchName() << "\n";
      }
      _firstEvent = false;
    }

    if(stepsHandles.empty())
    {
      throw cet::exception("SIM")<<"mu2e::CrvStepsFromStepPointMCs: No StepPointMC collections found for CRV" << endl;
    }

    // diagnostic counters
    unsigned nspmcs(0);
    // Loop over StepPointMC collections
    for(auto const& handle : stepsHandles) 
    {
      StepPointMCCollection const& steps(*handle);
      nspmcs += steps.size();

      // Loop over the StepPointMCs in this collection and sort them by Crv counter and SimParticle
      StepPointMap stepPointMap; // map of step points by Crv counter and sim particle
      fillMap(handle, stepPointMap, pdt);

      // convert the CrvBarIndex/SimParticle pair steps into CrvStep objects and fill the collection.  
      for(auto istepPointMap = stepPointMap.begin(); istepPointMap != stepPointMap.end(); ++istepPointMap)
      {
	auto& spmcptrv = istepPointMap->second;  // step pointer vector
        stable_sort(spmcptrv.begin(), spmcptrv.end(), compareStepTimes); //TODO: can be removed, if all StepPointMCs are time ordered

	auto const& trackId  = istepPointMap->first.second; // track id
	auto const& barIndex = istepPointMap->first.first; //bar index
	auto const& simptr = spmcptrv.front()->simParticle(); //sim particle
	auto pref = pdt->particle(simptr->pdgId());
	if(!pref.isValid()) continue;
	ParticleData const* pdata = &pref.ref();

        // indices in the StepPointMC vector where a new CrvStep starts and ends
        vector<pair<size_t,size_t> > spmcIndices;

        // fill the CrvSteps
	fillSteps(spmcptrv,barIndex,pdata,trackId,crvSteps,spmcIndices);
        size_t nNewCrvSteps = spmcIndices.size(); //number of CrvSteps created

        for(size_t i=0; i<nNewCrvSteps; ++i)
        {
          int iCrvStep = (int)crvSteps->size()-(int)nNewCrvSteps+i;
          CrvStep const& crvStep = crvSteps->at(iCrvStep);
          if(_diag>1) fillStepDiag(crvStep,spmcptrv,spmcIndices[i].first,spmcIndices[i].second);
          if(_debug>1)
          {
	    cout << " CrvSteps with " << spmcIndices[i].second-spmcIndices[i].first+1 << " steps,"
                 << " barIndex = " << crvStep.barIndex()  << " SimParticle Key = " << crvStep.simParticle()->id()
                 << " edep = " << crvStep.visibleEDep() << " pathlen = " << crvStep.pathLength()
                 << " glen = " << sqrt((crvStep.endPos()-crvStep.startPos()).mag2())
                 << " startTime = " << crvStep.startTime() << " endTime = " << crvStep.endTime() << endl;
          }
        }
      } // end of pair loop
    } // end of collection loop
    if(_debug > 0)
    {
      cout << "Total number of CrvSteps " << crvSteps->size() << " , StepPointMCs = " << nspmcs << endl;
    }

    event.put(move(crvSteps));
  }

  void CrvStepsFromStepPointMCs::fillSteps(SPMCPtrV const& spmcptrv, CRSScintillatorBarIndex const& barIndex,
                                           ParticleData const* pdata, cet::map_vector_key trackId, 
                                           unique_ptr<CrvStepCollection> &crvSteps, vector<pair<size_t,size_t> > &spmcIndices)
  {
    //the sequence of steps needs to be devided into two (or more) CrvSteps, if the particle leaves the scintillator and returns later
    //this situtation can be identified by the end process code of Transportation
    bool newCrvStep(true);

    // variables we accumulate for all StepPoints in this pair
    double edep(0.0), pathlen(0.0);
    // keep track of the first and last step
    SPMCPtr first;
    SPMCPtr last;
    // loop over all StepPoints
    for(size_t i=0; i<spmcptrv.size(); ++i)
    {
      const auto& spmcptr = spmcptrv[i];
      if(newCrvStep)
      {
        newCrvStep=false;
        edep      =0.0;
        pathlen   =0.0;
        first     =spmcptr;
        spmcIndices.emplace_back(i,i);
      }

      if(_useTotalEDep)
        edep+=spmcptr->totalEDep();
      else
        edep+=spmcptr->visibleEDep();

      pathlen+=spmcptr->stepLength();
      last    =spmcptr;
      spmcIndices.back().second=i;

      if(spmcptr->endProcessCode()==ProcessCode::Transportation) newCrvStep=true;

      //analyze the accumulated steps at the end of the sequence or if the track leaves the scintillator
      if(i+1==spmcptrv.size() || newCrvStep)
      {
        if(first.isNull() || last.isNull())
          throw cet::exception("SIM")<<"mu2e::CrvStepsFromStepPointMCs: No first or last step" << endl;

        // Define the position at entrance and exit
        CLHEP::Hep3Vector startPos = first->position();
        CLHEP::Hep3Vector endPos   = last->postPosition();
        if(_noPostPositionAvailable) endPos = last->position() + last->momentum().unit()*last->stepLength();

        double startTime = first->time();
        double endTime   = last->time();

        CLHEP::Hep3Vector startMomV = first->momentum();
        CLHEP::Hep3Vector endMomV   = last->momentum();

        //endMomV is the start momentum of the last StepPointMC
        //need to calculate the end momentum based on the energy loss
        double mass            = pdata->mass();
        double endMom2         = endMomV.mag2();
        double endEnergyBefore = sqrt(endMom2 + mass*mass);
        double endEnergyAfter  = endEnergyBefore - last->totalEDep(); //TODO: does not take the energy of daughter particles into account

        endMom2          = endEnergyAfter*endEnergyAfter - mass*mass;
        double endMom    = 0.0;
        if(endMom2>0.0) endMom=sqrt(endMom2);

        //calculating the end time of the last step
        double avgEnergy = 0.5*(endEnergyBefore+endEnergyAfter);
        double avgGamma  = avgEnergy/mass;
        double avgBeta   = sqrt(1.0-1.0/(avgGamma*avgGamma));
        double velocity  = avgBeta*CLHEP::c_light;
        endTime         += last->stepLength()/velocity;

        // create the CrvStep and emplace it back into the vector of crvSteps
        crvSteps->emplace_back(first->barIndex(), edep, startTime, endTime, 
                               Geom::toXYZVec(startPos), Geom::toXYZVec(endPos),
                               Geom::toXYZVec(startMomV), endMom,
                               pathlen, first->simParticle());
      }
    }
  }

  void CrvStepsFromStepPointMCs::fillMap(SPMCCH const& spmcch, StepPointMap& stepPointMap, GlobalConstantsHandle<ParticleDataTable> &pdt)
  {
    StepPointMCCollection const& steps(*spmcch);
    for(size_t ispmc=0; ispmc<steps.size(); ++ispmc) 
    {
      const auto& step = steps[ispmc];

      if(_removeNeutralParticles)
      {
        auto pref = pdt->particle(step.simParticle()->pdgId());
        bool charged=true;
        if(pref.isValid()) charged=(pref.ref().charge()!=0.0);
        if(!charged) continue;  //skip neutral particles
      }

      // create key
      const CRSScintillatorBarIndex &barIndex = step.barIndex();
      cet::map_vector_key trackId = step.simParticle().get()->id();
      BarIndexTrackPair barIndexTrackPair(barIndex,trackId);

      // create ptr to this step
      SPMCPtr  spmcptr(spmcch,ispmc);

      // get vector of steps from map for this key, or create the vector, if the key doesn't exist yet
      vector<SPMCPtr> &spmcptrv = stepPointMap[barIndexTrackPair];
      if(spmcptrv.empty()) spmcptrv.reserve(_ssize);
      spmcptrv.push_back(spmcptr);
    }
  }

  void CrvStepsFromStepPointMCs::fillStepDiag(CrvStep const& crvStep, SPMCPtrV const& spmcptrv, size_t firstIndex, size_t lastIndex) 
  {
    _cslen = crvStep.pathLength();
    _csdist = sqrt((crvStep.endPos()-crvStep.startPos()).mag2());
    _csedep = crvStep.visibleEDep();
    auto const& spmcptr = spmcptrv.front();
    _csPartP = spmcptr->momentum().mag();
    _csPartPDG = spmcptr->simParticle()->pdgId();

    _slen.clear();
    _sdist.clear();
    _sedep.clear();
    _sPartP.clear();
    _sPartPDG.clear();
    for(size_t i=firstIndex; i<=lastIndex; ++i)
    {
      auto const& spmc = *spmcptrv.at(i);
      _slen.push_back(spmc.stepLength());
      _sdist.push_back(sqrt((spmc.postPosition()-spmc.position()).mag2())); 
      _sedep.push_back(spmc.visibleEDep()); 
      _sPartP.push_back(spmc.simParticle()->startMomentum().vect().mag());
      _sPartPDG.push_back(spmc.simParticle()->pdgId());
    }
    _csdiag->Fill();
  }

}

DEFINE_ART_MODULE(mu2e::CrvStepsFromStepPointMCs)
