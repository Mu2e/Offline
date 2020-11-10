#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <cmath>



namespace mu2e {


  class CaloCrystalHitFromHit : public art::EDProducer 
  {
    public:
       struct Config 
       {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;        
           fhicl::Atom<std::string> caloDigisModuleLabel{ Name("caloDigisModuleLabel"), Comment("Calo Digi module label")};
           fhicl::Atom<double>      time4Merge          { Name("time4Merge"),           Comment("Maximal time differnce to merge two RoID signals (ns)")};
           fhicl::Atom<int>         diagLevel           { Name("diagLevel"),            Comment("Diagnosis level")};
       };

       explicit CaloCrystalHitFromHit(const art::EDProducer::Table<Config>& config) :
         EDProducer{config},
         caloDigisToken_ {consumes<CaloRecoDigiCollection>(config().caloDigisModuleLabel())},
         time4Merge_     (config().time4Merge()),
         diagLevel_      (config().diagLevel())
       {
           produces<CaloCrystalHitCollection>();
       }

       void beginJob() override;
       void produce(art::Event& e) override;


    private:
       typedef art::Ptr<CaloRecoDigi> CaloRecoDigiPtr;

       void makeCaloHits(CaloCrystalHitCollection&,const art::ValidHandle<CaloRecoDigiCollection>&);
       void fillBuffer  (int, int, double, double, double, double, std::vector<CaloRecoDigiPtr>&, CaloCrystalHitCollection&);

       const art::ProductToken<CaloRecoDigiCollection> caloDigisToken_;      
       double                                          time4Merge_;
       int                                             diagLevel_;

       TH1F*  hEdep_;
       TH1F*  hTime_;
       TH1F*  hNRo_;
       TH1F*  hEdep_Cry_;
       TH1F*  hDelta_;
       TH2F*  hNRo2_;
       TH1F*  hEdep1_;
       TH1F*  hEdep2_;
  };


  //--------------------------------------------
  void CaloCrystalHitFromHit::beginJob()
  {
      if (diagLevel_ > 2) 
      {
          art::ServiceHandle<art::TFileService> tfs;
          hEdep_     = tfs->make<TH1F>("hEdep",   "Hit energy deposition",        200,   0.,  500);
          hTime_     = tfs->make<TH1F>("hTime",   "Hit time ",                  12000,   0., 2000);
          hNRo_      = tfs->make<TH1F>("hNRo",    "Number RO ",                    10,   0.,   10);
          hEdep_Cry_ = tfs->make<TH1F>("hEdepCry","Energy deposited per crystal",2000,   0., 2000);
          hDelta_    = tfs->make<TH1F>("hDelta",  "Hit time difference",          200, -20,    20);
          hNRo2_     = tfs->make<TH2F>("hNRo2",   "Number RO ",                    5,    0., 5, 50, 0, 50);
          hEdep1_    = tfs->make<TH1F>("hEdep1",  "Hit energy deposition",        200,   0.,  100);
          hEdep2_    = tfs->make<TH1F>("hEdep2",  "Hit energy deposition",        200,   0.,  100);
      }
  }


  //------------------------------------------------------------
  void CaloCrystalHitFromHit::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout<<"[CaloCrystalHitFromHit::produce] end"<<std::endl;
      
      const auto& recoCaloDigisHandle = event.getValidHandle(caloDigisToken_);
      auto caloHits                   = std::make_unique<CaloCrystalHitCollection>();
      
      makeCaloHits(*caloHits, recoCaloDigisHandle);
      
      event.put(std::move(caloHits));
      
      if (diagLevel_ > 0) std::cout<<"[CaloCrystalHitFromHit::produce] end"<<std::endl;
  }


  //--------------------------------------------------------------------------------------------------------------
  void CaloCrystalHitFromHit::makeCaloHits(CaloCrystalHitCollection& caloHits, const art::ValidHandle<CaloRecoDigiCollection>& recoCaloDigisHandle)
  {
    const Calorimeter& cal = *(GeomHandle<Calorimeter>());
    const auto& recoCaloDigis = *recoCaloDigisHandle;
    if (recoCaloDigis.empty()) return;
    const CaloRecoDigi* base = &recoCaloDigis.front(); 


    // fill the map that associate for each crystal the corresponding CaloRecoDigi indexes
    std::vector<std::vector<const CaloRecoDigi*>> hitMap(cal.nRO(),std::vector<const CaloRecoDigi*>());
    
    for (unsigned i=0; i< recoCaloDigis.size(); ++i)
    {
        int crystalId = cal.caloInfo().crystalByRO(recoCaloDigis[i].ROid());
        hitMap[crystalId].push_back(&recoCaloDigis[i]);
    }

    float totEnergyRec(0);
    for (unsigned crystalId=0;crystalId<hitMap.size();++crystalId)
    {
        std::vector<const CaloRecoDigi*> &hits = hitMap[crystalId];
        if (hits.empty()) continue;

        std::sort(hits.begin(),hits.end(),[](const auto& a, const auto& b){return a->time() < b->time();});

        auto startHit = hits.begin();
        auto endHit   = hits.begin();

        std::vector<CaloRecoDigiPtr> buffer;
        int nRoid(0);
        double timeW(0),timeWtot(0),eDepTot(0),eDepTotErr(0);

        while (endHit != hits.end())
        {
            double deltaTime = (*endHit)->time()-(*startHit)->time();
            if (diagLevel_ > 2) hDelta_->Fill(deltaTime);

            if (deltaTime > time4Merge_)
            {
                double time    = timeW/timeWtot;
                double timeErr = 1.0/sqrt(timeWtot);

                fillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, caloHits);
                totEnergyRec += eDepTot/float(nRoid);

                buffer.clear();
                timeW      = 0.0;
                timeWtot   = 0.0;
                eDepTot    = 0.0;
                eDepTotErr = 0.0;
                nRoid      = 0;
                startHit   = endHit;
             } 
             else
             {
                double wt  = 1.0/(*endHit)->timeErr()/(*endHit)->timeErr();
                timeWtot   += wt;
                timeW      += wt*(*endHit)->time();

                eDepTot    += (*endHit)->energyDep();
                eDepTotErr += (*endHit)->energyDepErr() * (*endHit)->energyDepErr();

                ++nRoid;

                size_t index = *endHit - base;
                buffer.push_back(art::Ptr<CaloRecoDigi>(recoCaloDigisHandle, index));

                ++endHit;
              }

          }

          //flush last buffer
          double time    = timeW/timeWtot;
          double timeErr = 1.0/sqrt(timeWtot);
          fillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, caloHits);
          totEnergyRec += eDepTot/float(nRoid);
      }

      if ( diagLevel_ > 1 ) std::cout<<"[CaloCrystalHitFromHit::produce] produced RecoCrystalHits with caloHits.size() = "<<caloHits.size()<<std::endl;
      if ( diagLevel_ > 1 ) std::cout<<"[CaloCrystalHitFromHit::produce] Total energy reconstructed = "<<totEnergyRec<<std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------
  void CaloCrystalHitFromHit::fillBuffer(int crystalId,int nRoid,double time,double timeErr,double eDep,double eDepErr,
                                         std::vector<CaloRecoDigiPtr>& buffer, CaloCrystalHitCollection& caloHits)
  {
      caloHits.emplace_back(CaloCrystalHit(crystalId, nRoid, time, timeErr, eDep, eDepErr, buffer));

      if (diagLevel_ > 2) std::cout<<"[CaloCrystalHitFromHit] created hit in crystal id="<<crystalId<<"\t with time="
                                   <<time<<"\t eDep="<<eDep<<"\t  from "<<nRoid<<" RO"<<std::endl;
                    
      if (diagLevel_ > 2)
      {
          hTime_->Fill(time);
          hEdep_->Fill(eDep);
          hNRo_->Fill(nRoid);
          hEdep_Cry_->Fill(crystalId,eDep);
          hNRo2_->Fill(nRoid,eDep);
          if (nRoid==1) hEdep1_->Fill(eDep);
          if (nRoid==2) hEdep2_->Fill(eDep);
      }
  }




}

DEFINE_ART_MODULE(mu2e::CaloCrystalHitFromHit);
