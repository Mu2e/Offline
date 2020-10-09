//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
// Original author Bertrand Echenard
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"



namespace mu2e {


  class CaloInspector : public art::EDAnalyzer {

     public:

       explicit CaloInspector(fhicl::ParameterSet const& pset);
       virtual ~CaloInspector() { }

       virtual void beginJob();
       virtual void endJob() {};

       virtual void analyze(const art::Event& e);


     private:
       std::string caloCrystalModuleLabel_;
       std::string caloShowerSimModuleLabel_;
       std::string caloDigiTruthModuleLabel_;
       int diagLevel_;

       TH1F *hcryMatchE_,*hcryNoMatchE_,*hsimMatchT_,*hsimMatchE_,*hsimNoMatchE_,*hsimNoMatchT_,*hsimMatchELT_;
       TH1F *hSimNoMatchCE_;

  };


  CaloInspector::CaloInspector(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    caloCrystalModuleLabel_     (pset.get<std::string>("caloCrystalModuleLabel")),
    caloShowerSimModuleLabel_   (pset.get<std::string>("caloShowerSimModuleLabel")),
    caloDigiTruthModuleLabel_   (pset.get<std::string>("caloDigiTruthModuleLabel")),
    diagLevel_                  (pset.get<int>("diagLevel",0))
  {}

  void CaloInspector::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;
       
       hcryMatchE_      = tfs->make<TH1F>("hcryMatchE",    "crystal Edep match;Edep (MeV);Entries",      100,0,100);
       hcryNoMatchE_    = tfs->make<TH1F>("hcryNoMatchE",  "crystal Edep NO match;Edep (MeV);Entries",   100,0,100);
       hsimMatchT_      = tfs->make<TH1F>("hsimMatchDT",   "sim Time match;DT (ns);Entries",             100,0,500);
       hsimMatchELT_    = tfs->make<TH1F>("hsimMatchELT",  "sim Edep large T Match;Edep (MeV);Entries",  100,0,100);
       hsimMatchE_      = tfs->make<TH1F>("hsimMatchE",    "sim Edep match;Edep (MeV);Entries",          100,0,100);
       hsimNoMatchE_    = tfs->make<TH1F>("hsimNoMatchE",  "sim Edep NO match;Edep (MeV);Entries",       100,0,100);
       hsimNoMatchT_    = tfs->make<TH1F>("hsimNoMatchT",  "sim Time match;time (ns);Entries",           100,0,2000);
       hSimNoMatchCE_   = tfs->make<TH1F>("hSimNoMatchCE", "sim E no match Ce match;Edep (MeV);Entries", 100,0,100);
  }



  void CaloInspector::analyze(const art::Event& event) 
  {      
      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(caloCrystalModuleLabel_, caloCrystalHitsHandle);
      const CaloCrystalHitCollection& caloCrystalHits(*caloCrystalHitsHandle);

      //Calorimeter shower sims
      art::Handle<CaloShowerSimCollection> caloShowerSimHandle;
      event.getByLabel(caloShowerSimModuleLabel_, caloShowerSimHandle);
      const CaloShowerSimCollection& caloShowerSims(*caloShowerSimHandle);
      
      //Calorimeter digi truth assignment
      art::Handle<CaloDigiMCTruthAssn> caloDigiTruthHandle;
      event.getByLabel(caloDigiTruthModuleLabel_, caloDigiTruthHandle);
      const CaloDigiMCTruthAssn& caloDigiTruth(*caloDigiTruthHandle);


      for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
      {
          const CaloCrystalHit& hit = caloCrystalHits.at(ic);

          //Find the caloDigiMC in the truth map          
          auto itMC = caloDigiTruth.begin();
          while (itMC != caloDigiTruth.end()) {if (itMC->first.get() == &hit) break; ++itMC;}
          unsigned nCrySims = (itMC != caloDigiTruth.end()) ? itMC->second->nParticles() : 0;
          
          if (nCrySims>0) hcryMatchE_->Fill(hit.energyDep());
          else            hcryNoMatchE_->Fill(hit.energyDep());
          
          if (diagLevel_ > 0) std::cout<<"Crystal "<<hit.id()<<" "<<hit.energyDep()<<" "<<hit.time()<<" "<<nCrySims<<std::endl;

          for (unsigned i=0;i< nCrySims;++i)
	  {	                      
	      const auto& eDepMC = itMC->second->energyDeposit(i);
              hsimMatchT_->Fill(eDepMC.time()-hit.time());
              if (eDepMC.time()-hit.time() > 5 ) hsimMatchELT_->Fill(eDepMC.energyDep());
              if (diagLevel_ > 0) std::cout<<"Sim "<<eDepMC.sim()->id().asInt()<<" "<<eDepMC.sim()->pdgId()<<" "<<eDepMC.sim()->creationCode()
                                           <<" "<<eDepMC.time()<<" "<<eDepMC.energyDep()<<" "<<eDepMC.momentumIn()<<std::endl;

          }
          if (diagLevel_ > 0) std::cout<<std::endl;
      }

      // set of all sims matched to reco hits
      std::set<int> simMatch;
      for (const auto& kv :caloDigiTruth ) 
      {
         for (const auto& edep : kv.second->energyDeposits()) simMatch.insert(edep.sim()->id().asInt());
      }

      //check sims that are not matched to calo hit
      for (const auto& showerSim : caloShowerSims)
      {
         int simId = showerSim.sim()->id().asInt();
         bool isMatched = simMatch.find(simId) != simMatch.end();
         if (isMatched)  hsimMatchE_->Fill(showerSim.energyDep() );
         else            hsimNoMatchE_->Fill(showerSim.energyDep() ); 
         
         if (!isMatched) hsimNoMatchT_->Fill(showerSim.time());      
      }
 

      double sumEno(0);
      for (const auto& showerSim : caloShowerSims)
      {
         int simId = showerSim.sim()->id().asInt();
         bool isMatched = simMatch.find(simId) != simMatch.end();
         
         auto parent(showerSim.sim());
         while (parent->hasParent()) parent = parent->parent();                     
	 bool isConversion(false);
         if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) isConversion=true;
         
         if (!isConversion) continue;
         
         
         if (!isMatched) sumEno += showerSim.energyDep();
      }

      hSimNoMatchCE_->Fill(sumEno);
  }



}  

DEFINE_ART_MODULE(mu2e::CaloInspector);
