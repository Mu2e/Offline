//Author: S Middleton
//Date: Nov 2019
//Purpose: For the purpose of crystal hit making from online

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/NewCaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/NewCaloRecoDigiCollection.hh"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <cmath>

namespace mu2e {


  class NewCaloCrystalHitFromHit : public art::EDProducer {
  public:

    explicit NewCaloCrystalHitFromHit(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloDigisToken_{consumes<NewCaloRecoDigiCollection>(pset.get<std::string>("caloDigisModuleLabel"))},
      time4Merge_          (pset.get<double>     ("time4Merge")),
      diagLevel_           (pset.get<int>        ("diagLevel",0))
    {
      produces<NewCaloCrystalHitCollection>();
    }

    void beginJob() override;
    void produce(art::Event& e) override;

  private:

    typedef art::Ptr<NewCaloRecoDigi> CaloRecoDigiPtr;

    art::ProductToken<NewCaloRecoDigiCollection> const caloDigisToken_;
    double      time4Merge_;
    int         diagLevel_;

    std::vector<std::vector<const NewCaloRecoDigi*>> hitMap_;  


    void MakeCaloCrystalHits(NewCaloCrystalHitCollection& CaloCrystalHits, const art::ValidHandle<NewCaloRecoDigiCollection>& recoCaloDigisHandle);

    void FillBuffer(int crystalId, int nRoid, double time, double timeErr, double eDep, double eDepErr,
                    std::vector<CaloRecoDigiPtr>& buffer, NewCaloCrystalHitCollection& CaloCrystalHits);
  };


  //--------------------------------------------
  void NewCaloCrystalHitFromHit::beginJob()
  {}


  //------------------------------------------------------------
  void NewCaloCrystalHitFromHit::produce(art::Event& event)
  {
    auto const& recoCaloDigisHandle = event.getValidHandle(caloDigisToken_);
    auto CaloCrystalHits = std::make_unique<NewCaloCrystalHitCollection>();
    MakeCaloCrystalHits(*CaloCrystalHits, recoCaloDigisHandle);
    event.put(std::move(CaloCrystalHits));
  }


  //--------------------------------------------------------------------------------------------------------------
  void NewCaloCrystalHitFromHit::MakeCaloCrystalHits(NewCaloCrystalHitCollection& CaloCrystalHits, art::ValidHandle<NewCaloRecoDigiCollection> const& recoCaloDigisHandle)
  {
    Calorimeter const &cal = *(GeomHandle<Calorimeter>());
    auto const& recoCaloDigis = *recoCaloDigisHandle;
    NewCaloRecoDigi const* base = &recoCaloDigis.front(); // TODO What if recoCaloDigis is empty?


    if (cal.nRO() > int(hitMap_.size()))
      for (int i = hitMap_.size(); i<= cal.nRO(); ++i) hitMap_.push_back(std::vector<const NewCaloRecoDigi*>());
    for (size_t i=0; i<hitMap_.size(); ++i) hitMap_[i].clear();

    for (unsigned int i=0; i< recoCaloDigis.size(); ++i)
      {
        int crystalId = cal.caloInfo().crystalByRO(recoCaloDigis[i].ROid());
        hitMap_[crystalId].push_back(&recoCaloDigis[i]);
      }


    for (unsigned int crystalId=0;crystalId<hitMap_.size();++crystalId)
      {
        std::vector<const NewCaloRecoDigi*> &hits = hitMap_[crystalId];
        
        //check if empty:
        if (hits.empty()) continue;

        //sort hits in terms of time:
        std::sort(hits.begin(),hits.end(),[](const auto a, const auto b){return a->time() < b->time();});

        //find hit:
        auto startHit = hits.begin();
        auto endHit   = hits.begin();

        //create a buffer 
        std::vector<CaloRecoDigiPtr> buffer;
        double timeW(0);
        double eDepTot(0),eDepTotErr(0);
        int nRoid(0);

        //loop through hits:
        while (endHit != hits.end())
          {
            //time:
            double deltaTime = (*endHit)->time()-(*startHit)->time();
        
            //if > than set merge time:
            if (deltaTime > time4Merge_) 
              {
                double time = timeW/nRoid;
                double timeErr = 0;
                //fill that buffer: //TODO -->peakpos, flags etc.
                FillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, CaloCrystalHits);
                //clear:
                buffer.clear();
                timeW      = 0.0;
                eDepTot    = 0.0;
                eDepTotErr = 0.0;
                nRoid      = 0;
                startHit   = endHit;
              }
            else //if not then:
              {
                //add up times
                timeW      += (*endHit)->time();
                //add up energies
                eDepTot    += (*endHit)->energyDep();
                eDepTotErr += (*endHit)->energyDepErr() * (*endHit)->energyDepErr();
                //increment no. ROs
                ++nRoid;
                //get index:
                size_t index = *endHit - base;
                //add to end of hit buffer:
                buffer.push_back(art::Ptr<NewCaloRecoDigi>(recoCaloDigisHandle, index));
                //increment:
                ++endHit;
              }

          }
        //get time and error
        double time = timeW/nRoid;
        double timeErr = 0;
        //Finally fill buffer:
        FillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, CaloCrystalHits);
      }


    if ( diagLevel_ > 0 )
      {
        printf("[NewCaloCrystalHitFromHit::produce] produced RecoCrystalHits ");
        printf(": CaloCrystalHits.size()  = %i \n", int(CaloCrystalHits.size()));
      }

  }

  void NewCaloCrystalHitFromHit::FillBuffer(int const crystalId,
                                         int const nRoid,
                                         double const time,
                                         double const timeErr,
                                         double const eDep,
                                         double const eDepErr,
                                         std::vector<CaloRecoDigiPtr>& buffer,
                                         NewCaloCrystalHitCollection& CaloCrystalHits)
  {
 
    CaloCrystalHits.emplace_back(NewCaloCrystalHit(crystalId, nRoid, time, timeErr, eDep, eDepErr, buffer));

    if (diagLevel_ > 1)
      {
        std::cout<<"[NewCaloCrystalHitFromHit] created hit in crystal id="<<crystalId<<"\t with time="<<time<<"\t eDep="<<eDep<<"\t  from "<<nRoid<<" RO"<<std::endl;

        
      }
  }


}

DEFINE_ART_MODULE(mu2e::NewCaloCrystalHitFromHit);
