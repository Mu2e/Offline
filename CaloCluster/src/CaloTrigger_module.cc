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
#include "RecoDataProducts/inc/CaloTrigSeed.hh"
#include "RecoDataProducts/inc/CaloTrigSeedCollection.hh"

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


  class CaloTrigger : public art::EDProducer {

  public:

    explicit CaloTrigger(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloDigiModuleLabel_(pset.get<std::string>("caloDigiModuleLabel")),
      digiSampling_(       pset.get<double>("digiSampling")),
      windowPeak_(         pset.get<unsigned>("windowPeak")),
      minAmp_(             pset.get<unsigned>("minAmplitude")),
      minSeedAmp_(         pset.get<unsigned>("minSeedAmplitude")),
      extendSecond_(       pset.get<bool>("extendSecond")),
      blindTime_(          pset.get<double>("blindTime")),
      endTimeBuffer_(      pset.get<double>("endTimeBuffer")),
      minEnergy_(          pset.get<float>("minEnergy")),
      timeCorrection_(     pset.get<double>("timeCorrection")),
      adcToEnergy_(        pset.get<double>("adcToEnergy")),	   
      diagLevel_(          pset.get<int>("diagLevel",0)),	   
      DNTBINs_(          pset.get<int>("dntbins",1)),
      window_(2*windowPeak_+1),
      hitList_()
    {
      produces<CaloTrigSeedCollection>();
      seeds_.reserve(100);            
    }

    virtual ~CaloTrigger() {};
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
    float       minEnergy_; 
    double       timeCorrection_; 
    double       adcToEnergy_; 
    int          diagLevel_;
    int DNTBINs_;
    unsigned     window_;
       
    std::vector<std::list<FastHit>> hitList_;
    std::vector<FastHit*> seeds_;
    std::vector<CaloTrigSeed> trigseeds_;
    int mbtime_;

    art::ProductID               _crystalHitsPtrID;
    art::EDProductGetter const*  _crystalHitsPtrGetter;
 
    void extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle,
			 CaloTrigSeedCollection&                 trigSeeds);
       
  };


  //-------------------------------------------------------
  void CaloTrigger::produce(art::Event& event) 
  {

    if (diagLevel_ > 0) std::cout<<"[CaloTrigSeed::produce] begin"<<std::endl;
       
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    mbtime_ = accPar->deBuncherPeriod;
    if (diagLevel_ > 1) std::cout<<"Debuncher time: " << mbtime_ << std::endl;

    //Get the calorimeter Digis
    art::Handle<CaloDigiCollection> caloDigisHandle;
    event.getByLabel(caloDigiModuleLabel_, caloDigisHandle);

    auto trigSeedsColl = std::make_unique<CaloTrigSeedCollection>();
    trigSeedsColl->reserve(10);
        
    extractRecoDigi( caloDigisHandle, *trigSeedsColl);

    if ( diagLevel_ > 3 ){
      printf("[CaloTrigSeed::produce] produced trigSeedsColl size  = %i \n", int(trigSeedsColl->size()));
    }
    event.put(std::move(trigSeedsColl));
    
    if (diagLevel_ > 0) std::cout<<"[CaloTrigSeed::produce] end"<<std::endl;

    return;
  }




  //--------------------------------------------------------------------------------------
  void CaloTrigger::extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle,
				    CaloTrigSeedCollection&                 trigSeeds)
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

    for (const auto& caloDigi : caloDigis){
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
  
      std::deque<int> deque_;
      for (; it != waveform.end(); ++it){		
	while (!deque_.empty() && *it > deque_.back()) deque_.pop_back();
	deque_.push_back(*it);
	if (countdown > 0) countdown--;
	else{		   
	  if (deque_.front()> minAmp_ && deque_.front()== *std::prev(it,windowPeak_) && *std::prev(it,windowPeak_) != *std::prev(it,windowPeak_-1)){
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
    trigSeeds.clear();
    if (diagLevel_ > 1) std::cout << "Seed SIZE=" << seeds_.size() << std::endl; 
    for (auto& seed : seeds_){
      if (diagLevel_ > 1) std::cout << "Seed val=" << seed->val_ << " cryId=" << seed->crId_ << std::endl;
      if (seed->val_ < 1) continue; 
      int      cluEnergy(0);
      unsigned int idpeak= (unsigned int) seed->crId_;
      float epeak=seed->val_*adcToEnergy_;
      float xpeak= cal->crystal(seed->crId_).localPosition().x();
      float ypeak= cal->crystal(seed->crId_).localPosition().y();
      float rpeak;
      rpeak= sqrt(xpeak*xpeak+ypeak*ypeak);
      double xc(0.);
      double yc(0.);
      int ring1max(0);
      int ring1max2(0);
      std::queue<int> crystalRing1;
      for (const auto& nid : cal->neighbors(seed->crId_)) crystalRing1.push(nid);
      while (!crystalRing1.empty()){
	int nid = crystalRing1.front();
	for (int itimebin= seed->index_-DNTBINs_;itimebin<=seed->index_+DNTBINs_;++itimebin){
	  for (auto& hit : hitList_[itimebin]){
	    if (hit.crId_ != nid) continue;
	    if (hit.val_>ring1max){
	      ring1max2=ring1max;
	      ring1max=hit.val_;
	    }
	    else{
	      if (hit.val_>ring1max2){
		ring1max2=hit.val_;
	      }
	    }
	  }
	}
	crystalRing1.pop();
      }
      //
      int ring2max(0);
      std::queue<int> crystalRing2;
      for (const auto& nid : cal->nextNeighbors(seed->crId_)) crystalRing2.push(nid);
      while (!crystalRing2.empty()){
	int nid = crystalRing2.front();
	for (int itimebin= seed->index_-DNTBINs_;itimebin<=seed->index_+DNTBINs_;++itimebin){
	  for (auto& hit : hitList_[itimebin]){
	    if (hit.crId_ != nid) continue;
	    if (hit.val_>ring2max){
	      ring2max=hit.val_;
	    }
	  }
	}
	crystalRing2.pop();
      }
      float ring1emax= ring1max*adcToEnergy_;
      float ring1emax2= ring1max2*adcToEnergy_;
      float ring2emax= ring2max*adcToEnergy_;
      
      std::queue<int> crystalToVisit;
      crystalToVisit.push(seed->crId_);
      for (const auto& nid : cal->neighbors(seed->crId_)) crystalToVisit.push(nid);
      if (diagLevel_ > 1) std::cout << " To visit: " << crystalToVisit.size() << std::endl;

      while (!crystalToVisit.empty()){
	int nid = crystalToVisit.front();
	if (diagLevel_ > 1)  std::cout << " visiting: " << nid << std::endl;
	
	for (int itimebin= seed->index_-DNTBINs_;itimebin<=seed->index_+DNTBINs_;++itimebin){
	  for (auto& hit : hitList_[itimebin]){
	    if (hit.crId_ != nid || hit.val_ <0.5) continue;
	    cluEnergy += hit.val_;
	    xc += cal->crystal(nid).localPosition().x()*hit.val_;
	    yc += cal->crystal(nid).localPosition().y()*hit.val_;
	    hit.val_=-abs(hit.val_);
	    if (diagLevel_ > 1)  std::cout << " index add close to: " << nid << std::endl;
	    for (const auto& neighbor : cal->neighbors(nid)) crystalToVisit.push(neighbor);                
	    if (extendSecond_) for (const auto& nneighbor : cal->nextNeighbors(nid)) crystalToVisit.push(nneighbor);                
	  }
	}
	crystalToVisit.pop();
      }
         
      for (int itimebin= seed->index_-DNTBINs_;itimebin<=seed->index_+DNTBINs_;++itimebin){
	for (auto& hit : hitList_[itimebin]) hit.val_=abs(hit.val_);
      }
      float eDep = cluEnergy*adcToEnergy_;
      if (eDep > minEnergy_){
	xc /= cluEnergy;
	yc /= cluEnergy;
	float clutime  = (seed->index_+offsetT0_)*digiSampling_-timeCorrection_;
	
	// calorimeter trigger variables
	float tpeak    = (seed->index_+offsetT0_)*digiSampling_-timeCorrection_;
	//	int   iSection = cal->crystal(seed->crId_).diskId(); 

	CaloTrigSeed trigseed(idpeak,epeak,tpeak,rpeak,ring1emax,ring1emax2,ring2emax,eDep,clutime,xc,yc);
	trigSeeds.emplace_back(std::move(trigseed));
      }
    } 
  }

}

using mu2e::CaloTrigger;
DEFINE_ART_MODULE(CaloTrigger);
