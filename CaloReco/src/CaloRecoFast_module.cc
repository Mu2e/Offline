#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"


#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigiFast.hh"
#include "RecoDataProducts/inc/CaloRecoDigiFastCollection.hh"

#include <iostream>
#include <string>

#include "/usr/include/valgrind/callgrind.h"

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


  class CaloRecoFast : public art::EDProducer {

    public:

        explicit CaloRecoFast(fhicl::ParameterSet const& pset) :
          caloDigiModuleLabel_ (pset.get<std::string>("caloDigiModuleLabel")),
          digiSampling_        (pset.get<double>     ("digiSampling")),
          windowPeak_          (pset.get<unsigned>   ("windowPeak")),
          minAmp_              (pset.get<unsigned>   ("minAmplitude")),
          minSeedAmp_          (pset.get<unsigned>   ("minSeedAmplitude")),
          blindTime_           (pset.get<double>     ("blindTime")),
          endTimeBuffer_       (pset.get<double>     ("endTimeBuffer")),  // ns
          minEnergy_           (pset.get<double>     ("minEnergy")),  // ns
          diagLevel_           (pset.get<int>        ("diagLevel",0)),
          window_              (2*windowPeak_+1),
          hitList_(300,std::list<FastHit>())
          {
             produces<CaloRecoDigiFastCollection>();            
             seeds_.reserve(100);            
          }

        virtual ~CaloRecoFast() {};
	virtual void produce(art::Event& e);



    private:

       std::string  caloDigiModuleLabel_;
       double       digiSampling_;
       unsigned     windowPeak_;
       double       minAmp_;
       double       minSeedAmp_;
       double       blindTime_;
       double       endTimeBuffer_; 
       double       minEnergy_; 
       int          diagLevel_;
       unsigned     window_;
       
       std::deque<int> deque_;
       std::vector<std::list<FastHit>> hitList_;
       std::vector<FastHit*> seeds_;
       int mbtime_;
 
       void extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle,
                            CaloRecoDigiFastCollection &recoCaloHits);
       
  };


  //-------------------------------------------------------
  void CaloRecoFast::produce(art::Event& event)
  {

       if (diagLevel_ > 0) std::cout<<"[CaloRecoFast::produce] begin"<<std::endl;
       
       ConditionsHandle<AcceleratorParams> accPar("ignored");
       mbtime_ = accPar->deBuncherPeriod;


       //Get the calorimeter Digis
       art::Handle<CaloDigiCollection> caloDigisHandle;
       event.getByLabel(caloDigiModuleLabel_, caloDigisHandle);

       std::unique_ptr<CaloRecoDigiFastCollection> recoCaloDigiColl(new CaloRecoDigiFastCollection);
       recoCaloDigiColl->reserve(10);
       extractRecoDigi( caloDigisHandle, *recoCaloDigiColl);

       if ( diagLevel_ > 3 )
       {
           printf("[CaloRecoFast::produce] produced RecoCrystalHits ");
           printf(", recoCaloDigiColl size  = %i \n", int(recoCaloDigiColl->size()));
       }
       event.put(std::move(recoCaloDigiColl));

       if (diagLevel_ > 0) std::cout<<"[CaloRecoFast::produce] end"<<std::endl;

       return;
  }




  //--------------------------------------------------------------------------------------
  void CaloRecoFast::extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle,
                                            CaloRecoDigiFastCollection &recoCaloHits)
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
      std::sort(seeds_.begin(),seeds_.end(),[](const FastHit* a, const FastHit* b) {return a->val_ > b->val_;});
      
      
      //loop over seeds
      for (auto& seed : seeds_)
      {
         if (seed->val_ < 1) continue; 

         int cluEnergy(seed->val_);
         double xc = cal->crystal(seed->crId_).position().x()+3904;
         double yc = cal->crystal(seed->crId_).position().y();
         seed->val_ = 0; 
         
         for (const auto& nid : cal->neighbors(seed->crId_))
         {
            for (auto& hit : hitList_[seed->index_])
            {
               if (hit.crId_ != nid) continue;
               cluEnergy += hit.val_;
               hit.val_=0;
            }
            for (auto& hit : hitList_[seed->index_-1])
            {
               if (hit.crId_ != nid) continue;
               cluEnergy += hit.val_;
               hit.val_=0;
            }
            for (auto& hit : hitList_[seed->index_+1])
            {
               if (hit.crId_ != nid) continue;
               cluEnergy += hit.val_;
               hit.val_=0;
            }
         } 
         
         for (const auto& nid : cal->nextNeighbors(seed->crId_))
         {
            for (auto& hit : hitList_[seed->index_])
            {
               if (hit.crId_ != nid) continue;
               cluEnergy += hit.val_;
               hit.val_=0;
            }
            for (auto& hit : hitList_[seed->index_-1])
            {
               if (hit.crId_ != nid) continue;
               cluEnergy += hit.val_;
               hit.val_=0;
            }
            for (auto& hit : hitList_[seed->index_+1])
            {
               if (hit.crId_ != nid) continue;
               cluEnergy += hit.val_;
               hit.val_=0;
            }
         } 
         
         double time = (seed->index_+offsetT0_)*digiSampling_-20.0;
         double eDep = cluEnergy*0.04844/1.05;
         if (eDep > minEnergy_) recoCaloHits.emplace_back(CaloRecoDigiFast(time, eDep,xc,yc));                  
      } 

  }

}

using mu2e::CaloRecoFast;
DEFINE_ART_MODULE(CaloRecoFast);

