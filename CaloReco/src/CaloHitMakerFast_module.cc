#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "art/Framework/Principal/Handle.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloHit.hh"

   
namespace
{
    struct HitInfo {
       HitInfo(float time, float eDep) : nSiPM_(1), time_(time), eDep_(eDep) {}

       void add(float time, float eDep) {time_ = (nSiPM_*time_+time)/(nSiPM_+1); eDep_=(nSiPM_*eDep_+eDep)/(nSiPM_+1); ++nSiPM_;}

       unsigned nSiPM_;
       float    time_, eDep_;
    };
}



namespace mu2e {

  class CaloHitMakerFast : public art::EDProducer {

    public:
        struct Config 
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag> caloDigiCollection { Name("caloDigiCollection"), Comment("CaloDigi collection name") }; 
            fhicl::Atom<double>        digiSampling       { Name("digiSampling"),       Comment("Digitization time sampling") }; 
            fhicl::Atom<double>        deltaTPulses       { Name("deltaTPulses"),       Comment("Maximum time difference between two signals ") }; 
            fhicl::Atom<double>        pulseRatioMin      { Name("pulseRatioMin"),      Comment("Min value of energy ratio between pulses ") }; 
            fhicl::Atom<double>        pulseRatioMax      { Name("pulseRatioMax"),      Comment("Max value of energy ratio between pulses ") }; 
            fhicl::Sequence<double>    caphriCrystalID    { Name("caphriCrystalID"),    Comment("List of caphri crystal ID") }; 
            fhicl::Atom<int>           diagLevel          { Name("diagLevel"),          Comment("Diag Level"),0 };
        };

        explicit CaloHitMakerFast(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloDigisToken_    {consumes<CaloDigiCollection>(config().caloDigiCollection())},
           digiSampling_      (config().digiSampling()),
           deltaTPulses_      (config().deltaTPulses()),
           pulseRatioMin_     (config().pulseRatioMin()),
           pulseRatioMax_     (config().pulseRatioMax()),
           caphriCrystalID_   (config().caphriCrystalID()),
           diagLevel_         (config().diagLevel())
        {
            produces<CaloHitCollection>("calo");
            produces<CaloHitCollection>("caphri");
        }

        void produce(art::Event& e) override;


    private:
        using pulseMapType = std::unordered_map<unsigned, std::vector<HitInfo>>;

        void extractHits(const CaloDigiCollection& caloDigis, CaloHitCollection& caloHitsColl, CaloHitCollection& caphriHitsColl);
        void addPulse(pulseMapType& pulseMap, unsigned crystalID, float time, float eDep);

        art::ProductToken<CaloDigiCollection> caloDigisToken_;
        double              digiSampling_;
        double              deltaTPulses_;
        double              pulseRatioMin_;
        double              pulseRatioMax_;
        std::vector<double> caphriCrystalID_;
        int                 diagLevel_;
    };



    //--------------------------------------------------------------------------------------------------------------
    void CaloHitMakerFast::produce(art::Event& event)
    {
        if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi] begin"<<std::endl;

        const auto& caloDigisHandle = event.getValidHandle(caloDigisToken_);
        const auto& caloDigis(*caloDigisHandle);

        auto caloHitsColl   = std::make_unique<CaloHitCollection>();
        auto caphriHitsColl = std::make_unique<CaloHitCollection>();
        extractHits(caloDigis,*caloHitsColl,*caphriHitsColl);

        event.put(std::move(caloHitsColl),  "calo");
        event.put(std::move(caphriHitsColl),"caphri");

        if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi] end"<<std::endl;
        return;
    }


   //--------------------------------------------------------------------------------------------------------------
   void CaloHitMakerFast::extractHits(const CaloDigiCollection& caloDigis, CaloHitCollection& caloHitsColl, CaloHitCollection& caphriHitsColl)
   {
       const Calorimeter& cal = *(GeomHandle<Calorimeter>());
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

       pulseMapType pulseMap;
       for (const auto& caloDigi : caloDigis)
       {
           int crystalID   = cal.caloIDMapper().crystalIDFromSiPMID(caloDigi.SiPMID());

           double baseline = (caloDigi.waveform().at(0)+caloDigi.waveform().at(1)+caloDigi.waveform().at(2)+caloDigi.waveform().at(3))/4.0; 
           double eDep     = (caloDigi.waveform().at(caloDigi.peakpos())-baseline)*calorimeterCalibrations->ADC2MeV(caloDigi.SiPMID()); 
           double time     = caloDigi.t0() + caloDigi.peakpos()*digiSampling_;                      //Giani's definition
           //double time     = caloDigi.t0() + (caloDigi.peakpos()+0.5)*digiSampling_ - shiftTime_; //Bertrand's definition

           addPulse(pulseMap, crystalID, time, eDep);
           if (diagLevel_ > 2) std::cout<<"[FastRecoDigiFromDigi] extracted Digi with crystalID="<<crystalID<<" eDep="<<eDep<<"\t time=" <<time<<std::endl;
       }

       for (auto& crystal : pulseMap)
       {
	   int crID = crystal.first;
           bool isCaphri = std::find(caphriCrystalID_.begin(), caphriCrystalID_.end(), crID) != caphriCrystalID_.end();
	   for (auto& info : crystal.second)
           {
               if (diagLevel_ > 1) std::cout<<"[FastRecoDigiFromDigi] extracted Hit with crystalID="<<crID<<" eDep="<<info.eDep_<<"\t time=" <<info.time_<<"\t nSiPM= "<<info.nSiPM_<<std::endl;
               if (isCaphri) caphriHitsColl.emplace_back(CaloHit(crID, info.nSiPM_, info.time_,info.eDep_));
               else          caloHitsColl.emplace_back(  CaloHit(crID, info.nSiPM_, info.time_,info.eDep_)); 
	   }
       }

       if ( diagLevel_ > 0 ) std::cout<<"[CaloHitMakerFast] extracted "<<caloHitsColl.size()<<" CaloDigis"<<std::endl;
       if ( diagLevel_ > 0 ) std::cout<<"[CaloHitMakerFast] extracted "<<caphriHitsColl.size()<<" CapriDigis"<<std::endl;
   }



   //--------------------------------------------------------------------------------------------------------------
   void CaloHitMakerFast::addPulse(pulseMapType& pulseMap, unsigned crystalID, float time, float eDep)
   {
       bool addNewHit(true);
       for (auto& pulse : pulseMap[crystalID])
       {
	   double pulseRatio = eDep / pulse.eDep_;
	   if ( pulseRatio < pulseRatioMax_ && pulseRatio > pulseRatioMin_  && std::fabs(pulse.time_ - time) < deltaTPulses_)
	   {
              pulse.add(time,eDep); 
              addNewHit = false;
              break;
	   }
       }

       if (addNewHit) pulseMap[crystalID].push_back(HitInfo(time, eDep));  
   }
}

DEFINE_ART_MODULE(mu2e::CaloHitMakerFast);
