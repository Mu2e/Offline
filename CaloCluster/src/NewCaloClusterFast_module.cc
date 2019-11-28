#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/NewCaloDigi.hh"
#include "RecoDataProducts/inc/NewCaloDigiCollection.hh"
#include "RecoDataProducts/inc/NewCaloClusterCollection.hh"
#include "RecoDataProducts/inc/NewCaloCrystalHit.hh"
#include "RecoDataProducts/inc/NewCaloCrystalHitCollection.hh"

#include <iostream>
#include <string>
#include <queue>

namespace mu2e {


  class NewCaloClusterFast : public art::EDProducer {

  public:

    explicit NewCaloClusterFast(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloCrystalToken_{consumes<NewCaloCrystalHitCollection>(pset.get<std::string>("caloCrystalModuleLabel"))},
      digiSampling_(pset.get<double>("digiSampling")),
      windowPeak_(pset.get<unsigned>("windowPeak")),
      extendSecond_(pset.get<bool>("extendSecond")),
      blindTime_(pset.get<double>("blindTime")),
      endTimeBuffer_(pset.get<double>("endTimeBuffer")),
      minEnergy_( pset.get<double>("minEnergy")),
      timeCorrection_(pset.get<double>("timeCorrection")),
      adcToEnergy_(pset.get<double>("adcToEnergy")),
      diagLevel_(pset.get<int>("diagLevel",0)),
      includeCrystalHits_(pset.get<bool>("includeCrystalHits")),
      window_(2*windowPeak_+1)
    {
     
      produces<NewCaloClusterCollection>();

    }
    
    typedef std::list<const NewCaloCrystalHit*>    CaloCrystalList;

    virtual ~NewCaloClusterFast() {};
    virtual void produce(art::Event& e) override;

  private:
    art::ProductToken<NewCaloCrystalHitCollection> const caloCrystalToken_;
    double       digiSampling_;
    unsigned     windowPeak_;
    bool         extendSecond_;
    double       blindTime_;
    double       endTimeBuffer_;
    double       minEnergy_;
    double       timeCorrection_;
    double       adcToEnergy_;
    int          diagLevel_;
    bool         includeCrystalHits_;
    unsigned     window_;

    
    int mbtime_;

    art::ProductID                _crystalHitsPtrID;
    art::EDProductGetter const*  _crystalHitsPtrGetter;

    void MakeClusters(NewCaloClusterCollection& recoClusters, const art::Handle<NewCaloCrystalHitCollection> & recoCrystalHit);
    void FormCluster(const Calorimeter* cal, NewCaloCrystalHit *crystalSeed, std::vector<CaloCrystalList>& idHitVec)  ;
    void FillCluster(NewCaloClusterCollection& caloClustersColl, const CaloCrystalList& clusterPtrList, const art::Handle<NewCaloCrystalHitCollection>& CaloCrystalHitsHandle);
    
  };


  //-------------------------------------------------------
  void NewCaloClusterFast::produce(art::Event& event)
  {

       if (diagLevel_ > 0) std::cout<<"[CaloClusterFast::produce] begin"<<std::endl;

       ConditionsHandle<AcceleratorParams> accPar("ignored");
       mbtime_ = accPar->deBuncherPeriod; //TODO

       art::Handle<NewCaloCrystalHitCollection> caloCrystalHitsHandle;
       bool const success = event.getByToken(caloCrystalToken_, caloCrystalHitsHandle);
       if (!success) return;

       auto recoClustersColl = std::make_unique<NewCaloClusterCollection>();
      
       recoClustersColl->reserve(10);

       _crystalHitsPtrID     = event.getProductID<NewCaloCrystalHitCollection>();
       _crystalHitsPtrGetter = event.productGetter(_crystalHitsPtrID);

       MakeClusters( *recoClustersColl, caloCrystalHitsHandle);

       if ( diagLevel_ > 3 )
       {
           printf("[NewCaloClusterFast::produce] produced RecoCrystalHits ");
           printf(", recoClustersColl size  = %i \n", int(recoClustersColl->size()));
       }
    
       event.put(std::move(recoClustersColl));

       if (diagLevel_ > 0) std::cout<<"[NewCaloClusterFast::produce] end"<<std::endl;

       return;
  }



  void NewCaloClusterFast::MakeClusters(NewCaloClusterCollection& recoClusters, const art::Handle<NewCaloCrystalHitCollection> & CaloCrystalHitsHandle)
  {
      CaloCrystalList  clusterList;
      
      const NewCaloCrystalHitCollection& recoCrystalHits(*CaloCrystalHitsHandle);
      if (recoCrystalHits.empty()) return;

      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      const Calorimeter* cal = ch.get();
      
      unsigned offsetT0_ = unsigned(blindTime_/digiSampling_);
      
      if (recoCrystalHits.empty()) return;
      std::vector<CaloCrystalList>      clusterList, caloIdHitMap(cal->nCrystal());
      
      std::list<const NewCaloCrystalHit*>  seeds_;
    
      for (const auto& hit : recoCrystalHits)
      {
        if (hit.energyDep() < minEnergy_) continue;
        caloIdHitMap[hit.id()].push_back(&hit);
        seeds_.push_back(&hit);
       
      }

    seeds_.sort([](const NewCaloCrystalHit* a, const NewCaloCrystalHit* b) {return a->energyDep() > b->energyDep();});

    while( !seeds_.empty() )
      {
        const NewCaloCrystalHit* crystalSeed = *seeds_.begin();
        if (crystalSeed->energyDep() < minEnergy_) break;

        CaloCrystalList List = FormCluster(caloIdHitMap);

        clusterList.push_back(List);
      
        for (const auto& hit : List) { seeds_.remove(hit); } //TODO
      }
    for (auto cluster : clusterList) { FillCluster(recoClusters, cluster, CaloCrystalHitsHandle); }
}
    
CaloCrystalList   NewCaloFastCluster::FormCluster(const Calorimeter* cal_, NewCaloCrystalHit *crystalSeed_, std::vector<CaloCrystalList>& idHitVec)  
       { 

            CaloCrystalList clusterList_; 
            CaloCrystalList crystalsToVisit_;
    
            double seedTime_ crystalSeed->time();
            clusterList_.push_front(crystalSeed_);
            crystalToVisit_.push(crystalSeed_->id());  

            CaloCrystalList& cryL = idHitVec[crystalSeed_->id()];
            cryL.erase(std::find(cryL.begin(), cryL.end(), crystalSeed_));

            while (!crystalToVisit_.empty())
            {            
                 int visitId         = crystalToVisit_.front();
                 isVisited_[visitId] = 1;

                 std::vector<int> const& neighborsId = cal_->crystal(visitId).neighbors();
                 for (auto& iId : neighborsId)
                 {               
                     if (isVisited_[iId]) continue;
                     isVisited_[iId]=1;


                     CaloCrystalList& list = idHitVec[iId];
                     auto it=list.begin();
                     while(it != list.end())
                     {
                         NewCaloCrystalHit const* hit = *it;
                         if (std::abs(hit->time() - seedTime_) < deltaTime_).......//TODO
                         { 
                            if (hit->energyDep() < minEnergy_) { crystalToVisit_.push(iId); }
                            clusterList_.push_front(hit);
                            it = list.erase(it);   
                         } 
                         else {++it;}
                     } 
                     
                 }
                                       
                 crystalToVisit_.pop();                 
            }
            
           clusterList_.sort([] (NewCaloCrystalHit const* lhs, NewCaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();} );               
       } 
       return clusterList_;
}


//----------------------------------------------------------------------------------------------------------
  void NewCaloFastCluster::FillCluster(NewCaloClusterCollection& caloClustersColl, const CaloCrystalList& clusterPtrList, const art::Handle<NewCaloCrystalHitCollection>& CaloCrystalHitsHandle)
  {

    const NewCaloCrystalHitCollection& CaloCrystalHits(*CaloCrystalHitsHandle);
    const NewCaloCrystalHit* caloCrystalHitBase = &CaloCrystalHits.front();

    std::vector<art::Ptr<NewCaloCrystalHit>> caloCrystalHitsPtrVector;
    double totalEnergy(0),totalEnergyErr(0);
    
    for (auto clusterPrt : clusterPtrList)
      {
        
	    totalEnergy    += clusterPrt->energyDep();
        totalEnergyErr += clusterPrt->energyDepErr()*clusterPrt->energyDepErr();
        //TODO positions
        size_t idx = (clusterPrt - caloCrystalHitBase);
        caloCrystalHitsPtrVector.push_back( art::Ptr<NewCaloCrystalHit>(CaloCrystalHitsHandle,idx) );
      }

    totalEnergyErr = sqrt(totalEnergyErr);
    double time    = (*clusterPtrList.begin())->time();
    double timeErr = (*clusterPtrList.begin())->timeErr();
   

    caloClustersColl.emplace_back(NewCaloCluster(time,timeErr,totalEnergy,totalEnergyErr,caloCrystalHitsPtrVector,false)); //TODO


    if (diagLevel_ > 1)
      {
        std::cout<<"This cluster contains "<<clusterPtrList.size()<<" crystals, id= ";
        for (auto clusterPrt : clusterPtrList) std::cout<<clusterPrt->id()<<" ";
        std::cout<<" with energy="<<totalEnergy<<" and time="<<time<<std::endl;;
      }
  }

  


using mu2e::NewCaloClusterFast;
DEFINE_ART_MODULE(NewCaloClusterFast);


