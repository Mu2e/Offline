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
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include <iostream>
#include <string>
#include <queue>


namespace {

    struct FastHit
    {
       FastHit(int crId, int index, int val) : crId_(crId), index_(index), val_(val) {};
       int crId_;
       int index_;
       int val_;
    };
}




namespace mu2e {


  class CaloClusterFast : public art::EDProducer {

  public:

    explicit CaloClusterFast(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloDigiModuleLabel_(pset.get<std::string>("caloDigiModuleLabel")),
      digiSampling_(       pset.get<double>("digiSampling")),
      windowPeak_(         pset.get<unsigned>("windowPeak")),
      minAmp_(             pset.get<unsigned>("minAmplitude")),
      minSeedAmp_(         pset.get<unsigned>("minSeedAmplitude")),
      extendSecond_(       pset.get<bool>("extendSecond")),
      blindTime_(          pset.get<double>("blindTime")),
      endTimeBuffer_(      pset.get<double>("endTimeBuffer")),
      minEnergy_(          pset.get<double>("minEnergy")),
      timeCorrection_(     pset.get<double>("timeCorrection")),
      adcToEnergy_(        pset.get<double>("adcToEnergy")),
      diagLevel_(          pset.get<int>("diagLevel",0)),
      includeCrystalHits_( pset.get<bool>("includeCrystalHits")),
      window_(2*windowPeak_+1),
      hitList_()
    {
      produces<CaloCrystalHitCollection>();
      produces<CaloClusterCollection>();
      seeds_.reserve(100);
    }

    virtual ~CaloClusterFast() {};
    virtual void produce(art::Event& e);



  private:

    std::string  caloDigiModuleLabel_;
    double       digiSampling_;
    unsigned     windowPeak_;
    double       minAmp_;
    double       minSeedAmp_;
    bool         extendSecond_;
    double       blindTime_;
    double       endTimeBuffer_;
    double       minEnergy_;
    double       timeCorrection_;
    double       adcToEnergy_;
    int          diagLevel_;
    bool         includeCrystalHits_;
    unsigned     window_;

    std::deque<int> deque_;
    std::vector<std::list<FastHit>> hitList_;
    std::vector<FastHit*> seeds_;
    int mbtime_;

    art::ProductID                _crystalHitsPtrID;
    art::EDProductGetter const*  _crystalHitsPtrGetter;

    void extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle,
                         CaloClusterCollection&                 recoClusters,
                         CaloCrystalHitCollection&              recoCrystalHit);

  };


  //-------------------------------------------------------
  void CaloClusterFast::produce(art::Event& event)
  {

       if (diagLevel_ > 0) std::cout<<"[CaloClusterFast::produce] begin"<<std::endl;

       ConditionsHandle<AcceleratorParams> accPar("ignored");
       mbtime_ = accPar->deBuncherPeriod;


       //Get the calorimeter Digis
       art::Handle<CaloDigiCollection> caloDigisHandle;
       event.getByLabel(caloDigiModuleLabel_, caloDigisHandle);

       auto crystalHitColl   = std::make_unique<CaloCrystalHitCollection>();
       auto recoClustersColl = std::make_unique<CaloClusterCollection>();
       // std::unique_ptr<CaloCrystalHitCollection> crystalHitColl  (new CaloCrystalHitCollection);
       // std::unique_ptr<CaloClusterCollection>    recoClustersColl(new CaloClusterCollection);
       recoClustersColl->reserve(10);

       _crystalHitsPtrID     = event.getProductID<CaloCrystalHitCollection>();
       _crystalHitsPtrGetter = event.productGetter(_crystalHitsPtrID);

       extractRecoDigi( caloDigisHandle, *recoClustersColl, *crystalHitColl);

       if ( diagLevel_ > 3 )
       {
           printf("[CaloClusterFast::produce] produced RecoCrystalHits ");
           printf(", recoClustersColl size  = %i \n", int(recoClustersColl->size()));
       }
       event.put(std::move(crystalHitColl));
       event.put(std::move(recoClustersColl));

       if (diagLevel_ > 0) std::cout<<"[CaloClusterFast::produce] end"<<std::endl;

       return;
  }




  //--------------------------------------------------------------------------------------
  void CaloClusterFast::extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle,
                                        CaloClusterCollection&                 recoClusters,
                                        CaloCrystalHitCollection&              recoCrystalHits)
  {

      const CaloDigiCollection& caloDigis(*caloDigisHandle);

      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      const Calorimeter* cal = ch.get();
      int nro = cal->caloInfo().nROPerCrystal();

      unsigned offsetT0_ = unsigned(blindTime_/digiSampling_);
      unsigned nBinTime  = unsigned (mbtime_ - blindTime_ + endTimeBuffer_) / digiSampling_;
      unsigned winOffsetT0_ = windowPeak_+offsetT0_ ;

      if (hitList_.size() < nBinTime ) hitList_ = std::vector<std::list<FastHit>>(nBinTime,std::list<FastHit>());
      else for (auto& hl : hitList_) hl.clear();

      seeds_.clear();


      for (const auto& caloDigi : caloDigis)
      {
           if (caloDigi.roId() % nro) continue;

           double t0    = caloDigi.t0();
           int    roId  = caloDigi.roId();
           int    crId  = roId / nro; //cal->crystalByRO(roId); use short form for efficiency
           //double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
           const auto& waveform = caloDigi.waveform();

           auto it   = waveform.begin();
           auto tail = waveform.begin();
           unsigned countdown(window_ - 1);
           unsigned nCount(0);
           deque_.clear();

           for (; it != waveform.end(); ++it)
           {
                while (!deque_.empty() && *it > deque_.back()) deque_.pop_back();
                deque_.push_back(*it);

                if (countdown > 0) countdown--;
                else
                {
                   if (deque_.front()> minAmp_ && deque_.front()== *std::prev(it,windowPeak_) && *std::prev(it,windowPeak_) != *std::prev(it,windowPeak_-1))
                   {
                       int index = int(t0/digiSampling_) + nCount - winOffsetT0_;
                       hitList_[index].push_back(FastHit(crId,index,deque_.front()));
                       if (deque_.front()> minSeedAmp_) seeds_.emplace_back(&(hitList_[index].back()));
                   }

                   if (*tail == deque_.front()) deque_.pop_front();
                   ++tail;
                }
                ++nCount;
           }
      }

      if (seeds_.empty()) return;
      std::vector<art::Ptr<CaloRecoDigi> > caloRecoDigi;


      /* In this version, we do the usual seed + neighbors clustering. There is the option to look at next-to-neighbors crystals
         instead of simply neighbors crystals (extended = true).
         For the trigger, it might be useful to look at the energy in the first / second rings. In that case, it would be
         better to start the clustering considering only those crystals, and then continu egrowing the clusters by looking at
         neighbors.
         One could also implement a split-off cluyster recovery if it is fast enough.
      */
      for (auto& seed : seeds_)
      {
         if (seed->val_ < 1) continue;
         std::vector<art::Ptr<CaloCrystalHit>> hitsPtr;
         unsigned ncry(1);
         int      cluEnergy(seed->val_);
         double   xc = cal->crystal(seed->crId_).localPosition().x()*seed->val_;
         double   yc = cal->crystal(seed->crId_).localPosition().y()*seed->val_;
         //seed->val_ = 0;

         std::queue<int> crystalToVisit;
         for (const auto& nid : cal->neighbors(seed->crId_)) crystalToVisit.push(nid);

         double     hitTime = (seed->index_+offsetT0_)*digiSampling_-timeCorrection_;
         if (includeCrystalHits_) {
           recoCrystalHits.push_back(CaloCrystalHit(seed->crId_, 2, hitTime, 0, seed->val_*adcToEnergy_, 0., caloRecoDigi));
         }
         seed->val_ = 0;

         while (!crystalToVisit.empty())
         {
            int nid = crystalToVisit.front();

            for (auto& hit : hitList_[seed->index_])
            {
               if (hit.crId_ != nid || hit.val_ <0.5) continue;
               cluEnergy += hit.val_;
               xc += cal->crystal(nid).localPosition().x()*hit.val_;
               yc += cal->crystal(nid).localPosition().y()*hit.val_;
               ++ncry;
               if (includeCrystalHits_)  {
                 hitTime = (seed->index_+offsetT0_)*digiSampling_-timeCorrection_;
                 recoCrystalHits.push_back(CaloCrystalHit(hit.crId_, 2, hitTime, 0, hit.val_*adcToEnergy_, 0., caloRecoDigi));
               }
               hit.val_ = 0;
               for (const auto& neighbor : cal->neighbors(nid)) crystalToVisit.push(neighbor);
               if (extendSecond_) for (const auto& nneighbor : cal->nextNeighbors(nid)) crystalToVisit.push(nneighbor);
            }
            for (auto& hit : hitList_[seed->index_-1])
            {
               if (hit.crId_ != nid || hit.val_ <0.5) continue;
               cluEnergy += hit.val_;
               xc += cal->crystal(nid).localPosition().x()*hit.val_;
               yc += cal->crystal(nid).localPosition().y()*hit.val_;
               ++ncry;
               if (includeCrystalHits_) {
                 hitTime = (seed->index_-1 + offsetT0_)*digiSampling_-timeCorrection_;
                 recoCrystalHits.push_back(CaloCrystalHit(hit.crId_, 2, hitTime, 0, hit.val_*adcToEnergy_, 0., caloRecoDigi));
               }
               hit.val_ = 0;
               for (const auto& neighbor : cal->neighbors(nid)) crystalToVisit.push(neighbor);
               if (extendSecond_) for (const auto& nneighbor : cal->nextNeighbors(nid)) crystalToVisit.push(nneighbor);
           }
            for (auto& hit : hitList_[seed->index_+1])
            {
               if (hit.crId_ != nid || hit.val_ <0.5) continue;
               cluEnergy += hit.val_;
               xc += cal->crystal(nid).localPosition().x()*hit.val_;
               yc += cal->crystal(nid).localPosition().y()*hit.val_;
               ++ncry;
               if (includeCrystalHits_) {
                 hitTime = (seed->index_+1 + offsetT0_)*digiSampling_-timeCorrection_;
                 recoCrystalHits.push_back(CaloCrystalHit(hit.crId_, 2, hitTime, 0, hit.val_*adcToEnergy_, 0., caloRecoDigi));
               }
               hit.val_ = 0;
               for (const auto& neighbor : cal->neighbors(nid)) crystalToVisit.push(neighbor);
               if (extendSecond_) for (const auto& nneighbor : cal->nextNeighbors(nid)) crystalToVisit.push(nneighbor);
            }

            crystalToVisit.pop();
         }

         if (includeCrystalHits_){
           int   firstCrystal_index = recoCrystalHits.size() - ncry;
           int   lastCrystal_index  = recoCrystalHits.size();
           for (int k=firstCrystal_index; k<lastCrystal_index; ++k){
             art::Ptr<CaloCrystalHit> aPtr = art::Ptr<CaloCrystalHit>(_crystalHitsPtrID, k, _crystalHitsPtrGetter);
             hitsPtr.push_back(aPtr);
           }
         }

         double eDep = cluEnergy*adcToEnergy_;
         if (eDep > minEnergy_)
         {
             xc /= cluEnergy;
             yc /= cluEnergy;

             double time  = (seed->index_+offsetT0_)*digiSampling_-timeCorrection_;
             int iSection = cal->crystal(seed->crId_).diskId();

             CaloCluster cluster(iSection,time,0.0,eDep,0.0,hitsPtr,ncry,0.0);
             cluster.cog3Vector(CLHEP::Hep3Vector(xc,yc,0));

             recoClusters.emplace_back(std::move(cluster));
         }

      }








  }

}

using mu2e::CaloClusterFast;
DEFINE_ART_MODULE(CaloClusterFast);
