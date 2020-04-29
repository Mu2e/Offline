#ifndef TimeAndPhiClusterFinderDiag_hh
#define TimeAndPhiClusterFinderDiag_hh

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "TrkPatRec/inc/TimeAndPhiClusterFinder_types.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

#include "TTree.h"


#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"



namespace mu2e {
 
  using TimeAndPhiClusterFinderTypes::Data_t;
  using TimeAndPhiClusterFinderTypes::Config;
  
  class TimeAndPhiClusterFinderDiag: public mu2e::ModuleHistToolBase  {

    public:
        
        explicit TimeAndPhiClusterFinderDiag(const Config& config);
        explicit TimeAndPhiClusterFinderDiag(const fhicl::ParameterSet& PSet);
       ~TimeAndPhiClusterFinderDiag() = default;

       virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override;
       virtual int fillHistograms(void* Data, int Mode = -1) override ;


    private:
 
       void setupBranches();

       bool           treeInit_;
       TTree*         trkdiag_;
       Data_t*        data_;
       bool           mcdiag_;
       art::InputTag  mcdigisTag_;
  };
     
  
  TimeAndPhiClusterFinderDiag::TimeAndPhiClusterFinderDiag(const Config& config) : 
    treeInit_(  false),
    mcdiag_(    config.mcDiag()),
    mcdigisTag_(config.strawDigiMCCollection())
  {}
  
  TimeAndPhiClusterFinderDiag::TimeAndPhiClusterFinderDiag(const fhicl::ParameterSet& pset) : 
    treeInit_(  false),
    mcdiag_(    pset.get<bool>(         "MCDiag",                true)),
    mcdigisTag_(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD"))
  {}

   

  int TimeAndPhiClusterFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs)
  {
      trkdiag_ = Tfs->make<TTree>("tpcdiag","time and phi cluster diagnostics");
      return 0;
  }
  
  void TimeAndPhiClusterFinderDiag::setupBranches()
  {          
      trkdiag_->Branch("iev",       &data_->iev_,       "iev/I");
      trkdiag_->Branch("Nch",       &data_->Nch_,       "Nch/I");
      trkdiag_->Branch("chSel",     &data_->chSel_,     "chSel[Nch]/I");
      trkdiag_->Branch("chTime",    &data_->chTime_,    "chTime[Nch]/F");
      trkdiag_->Branch("Ncal",      &data_->Ncal_,      "Ncal/I");
      trkdiag_->Branch("calTime",   &data_->calTime_,   "calTime[Ncal]/F");      
      trkdiag_->Branch("nhit1",     &data_->nhit1_,     "nhit1/I");
      trkdiag_->Branch("hitIdx1",   &data_->hitIdx1_,   "hitIdx1[nhit1]/I");         
      trkdiag_->Branch("nclu1",     &data_->nclu1_,     "nclu1[nhit1]/I");
      trkdiag_->Branch("nhit2",     &data_->nhit2_,     "nhit2/I");
      trkdiag_->Branch("hitIdx2",   &data_->hitIdx2_,   "hitIdx2[nhit2]/I");
      trkdiag_->Branch("nclu2",     &data_->nclu2_,     "nclu2[nhit2]/I");
      trkdiag_->Branch("clu2Time",  &data_->clu2Time_,  "clu2Time[nhit2]/F");
      trkdiag_->Branch("clu2Phi",   &data_->clu2Phi_,   "clu2Phi[nhit2]/F");
      trkdiag_->Branch("clu2Rad",   &data_->clu2Rad_,   "clu2Rad[nhit2]/F");
      trkdiag_->Branch("nhit3",     &data_->nhit3_,     "nhit3/I");
      trkdiag_->Branch("hitIdx3",   &data_->hitIdx3_,   "hitIdx3[nhit3]/I");          
      trkdiag_->Branch("nclu3",     &data_->nclu3_,     "nclu3[nhit3]/I");
      trkdiag_->Branch("nhitMVA",   &data_->nhitMVA_,   "nhitMVA/I");
      trkdiag_->Branch("ncluMVA",   &data_->ncluMVA_,   "ncluMVA[nhitMVA]/I");
      trkdiag_->Branch("hitIdxMVA", &data_->hitIdxMVA_, "hitIdxMVA[nhitMVA]/I");
      trkdiag_->Branch("MVAsel",    &data_->MVAsel_,    "MVAsel[nhitMVA]/I");
      trkdiag_->Branch("MVAsig",    &data_->MVAsig_,    "MVAsig[nhitMVA]/I");
      trkdiag_->Branch("MVAvar1",   &data_->MVAvar1_,   "MVAvar1[nhitMVA]/F");
      trkdiag_->Branch("MVAvar2",   &data_->MVAvar2_,   "MVAvar2[nhitMVA]/F");
      trkdiag_->Branch("MVAvar3",   &data_->MVAvar3_,   "MVAvar3[nhitMVA]/F");
      trkdiag_->Branch("MVAvar4",   &data_->MVAvar4_,   "MVAvar4[nhitMVA]/F");
      trkdiag_->Branch("MVAvar5",   &data_->MVAvar5_,   "MVAvar5[nhitMVA]/F");
      
      treeInit_ = true;
  }

  
  int TimeAndPhiClusterFinderDiag::fillHistograms(void* Data, int Mode)
  {
     data_ = static_cast<Data_t*>(Data);

     // we need to set up the branches once the data_ structure has been allocated
     if (!treeInit_) setupBranches();
     
     // fix the MC truth for the MVA, and set MVAsel for conversion clusters
     if (mcdiag_)
     {
         auto& mcdigis = *data_->event_->getValidHandle<StrawDigiMCCollection>(mcdigisTag_);

         std::set<int> signalCluster;
         for (int i=0;i<data_->nhitMVA_;++i) 
         {
            std::vector<StrawDigiIndex> dids;
            data_->chcol_->fillStrawDigiIndices(*data_->event_,data_->hitIdxMVA_[i],dids);      
            const StrawDigiMC& mcdigi        = mcdigis.at(dids[0]);// taking 1st digi: is there a better idea??
            const art::Ptr<SimParticle>& spp = mcdigi.earlyStrawGasStep()->simParticle();
            int genId(-1);
            if (spp->genParticle().isNonnull()) genId = spp->genParticle()->generatorId().id();
            data_->MVAsig_[i]  = genId;  
            
            data_->MVAsel_[i]=0;
            if (genId==2) signalCluster.insert(data_->ncluMVA_[i]);
         }
         
         for (int i=0;i<data_->nhitMVA_;++i) 
            if (signalCluster.find(data_->ncluMVA_[i]) != signalCluster.end()) data_->MVAsel_[i]=1;
     }     
     
     trkdiag_->Fill();
     return 0; 
  }
   
}


DEFINE_ART_CLASS_TOOL(mu2e::TimeAndPhiClusterFinderDiag)
#endif
