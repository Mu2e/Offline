#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/CalorimeterCalibrations.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"


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
            fhicl::Atom<art::InputTag> pbttoken           { Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer")};
            fhicl::Atom<double>        digiSampling       { Name("digiSampling"),       Comment("Digitization time sampling") };
            fhicl::Atom<double>        deltaTPulses       { Name("deltaTPulses"),       Comment("Maximum time difference between two signals") };
            fhicl::Atom<double>        nPEperMeV          { Name("nPEperMeV"),          Comment("number of photo-electrons per MeV") };
            fhicl::Atom<double>        noiseLevelMeV      { Name("noiseLevelMeV"),      Comment("Noise level in MeV") };
            fhicl::Atom<double>        nSigmaNoise        { Name("nSigmaNoise"),        Comment("Maxnumber of sigma Noise to combine digi") };
            fhicl::Sequence<int>       caphriCrystalID    { Name("caphriCrystalID"),    Comment("List of caphri crystal ID") };
            fhicl::Atom<int>           diagLevel          { Name("diagLevel"),          Comment("Diag Level"),0 };
        };

        explicit CaloHitMakerFast(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloDigisToken_    {consumes<CaloDigiCollection>(config().caloDigiCollection())},
           pbttoken_          {consumes<ProtonBunchTime>(config().pbttoken())},
           digiSampling_      (config().digiSampling()),
           deltaTPulses_      (config().deltaTPulses()),
           nPEperMeV_         (config().nPEperMeV()),
           noise2_            (config().noiseLevelMeV()*config().noiseLevelMeV()),
           nSigmaNoise_       (config().nSigmaNoise()),
           caphriCrystalID_   (config().caphriCrystalID()),
           diagLevel_         (config().diagLevel())
        {
            produces<CaloHitCollection>("calo");
            produces<CaloHitCollection>("caphri");
            produces<IntensityInfoCalo>();
        }

        void produce(art::Event& e) override;


    private:
        using pulseMapType = std::unordered_map<unsigned, std::vector<HitInfo>>;

    void extractHits(const CaloDigiCollection& caloDigis, CaloHitCollection& caloHitsColl, CaloHitCollection& caphriHitsColl, IntensityInfoCalo&intCalo, double pbtOffset);
        void addPulse(pulseMapType& pulseMap, unsigned crystalID, float time, float eDep);

        art::ProductToken<CaloDigiCollection> caloDigisToken_;
        const  art::ProductToken<ProtonBunchTime>    pbttoken_;
        double              digiSampling_;
        double              deltaTPulses_;
        double              nPEperMeV_;
        double              noise2_;
        double              nSigmaNoise_;
        std::vector<int>    caphriCrystalID_;
        int                 diagLevel_;
    };



    //--------------------------------------------------------------------------------------------------------------
    void CaloHitMakerFast::produce(art::Event& event)
    {
        if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi] begin"<<std::endl;

        auto pbtH = event.getValidHandle(pbttoken_);
        const ProtonBunchTime& pbt(*pbtH);
        double pbtOffset = pbt.pbtime_;

        const auto& caloDigisHandle = event.getValidHandle(caloDigisToken_);
        const auto& caloDigis(*caloDigisHandle);

        auto caloHitsColl   = std::make_unique<CaloHitCollection>();
        auto caphriHitsColl = std::make_unique<CaloHitCollection>();
        auto intInfo        = std::make_unique<IntensityInfoCalo>();
        extractHits(caloDigis,*caloHitsColl,*caphriHitsColl,*intInfo,pbtOffset);

        event.put(std::move(caloHitsColl),  "calo");
        event.put(std::move(caphriHitsColl),"caphri");
        event.put(std::move(intInfo));

        if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi] end"<<std::endl;
        return;
    }


   //--------------------------------------------------------------------------------------------------------------
   void CaloHitMakerFast::extractHits(const CaloDigiCollection& caloDigis, CaloHitCollection& caloHitsColl, CaloHitCollection& caphriHitsColl, IntensityInfoCalo& intInfo, double pbtOffset)
   {

       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

       pulseMapType pulseMap;
       for (const auto& caloDigi : caloDigis)
       {
           int crystalID   = CaloSiPMId(caloDigi.SiPMID()).crystal().id();

           size_t nSamPed  = caloDigi.peakpos() > 3 ? 4 : std::max(caloDigi.peakpos()-1, 1);
           double baseline(0);
           for (size_t i=0; i<nSamPed; ++i){ baseline += caloDigi.waveform().at(i);}
           baseline /= nSamPed;
           double eDep     = (caloDigi.waveform().at(caloDigi.peakpos())-baseline)*calorimeterCalibrations->ADC2MeV(caloDigi.SiPMID());//FIXME! we should use the function ::Peak2MeV, I also think that we should: (i) discard the hit if eDep is <0 (noise/stange pulse), (ii) require a minimum pulse length. gianipez
           double time     = caloDigi.t0() + caloDigi.peakpos()*digiSampling_ - pbtOffset;                      //Giani's definition
           //double time     = caloDigi.t0() + (caloDigi.peakpos()+0.5)*digiSampling_ - shiftTime_; //Bertrand's definition

           addPulse(pulseMap, crystalID, time, eDep);
           if (diagLevel_ > 2) std::cout<<"[FastRecoDigiFromDigi] extracted Digi with crystalID="<<crystalID<<" eDep="<<eDep<<"\t time=" <<time<<std::endl;
       }

       unsigned short evtEnergy(0);
       for (auto& crystal : pulseMap)
       {
           int crID = crystal.first;
           bool isCaphri = std::find(caphriCrystalID_.begin(), caphriCrystalID_.end(), crID) != caphriCrystalID_.end();
           for (auto& info : crystal.second)
           {
               if (diagLevel_ > 1) std::cout<<"[FastRecoDigiFromDigi] extracted Hit with crystalID="<<crID<<" eDep="<<info.eDep_<<"\t time=" <<info.time_<<"\t nSiPM= "<<info.nSiPM_<<std::endl;
               if (isCaphri) caphriHitsColl.emplace_back(CaloHit(crID, info.nSiPM_, info.time_,info.eDep_));
               else          caloHitsColl.emplace_back(  CaloHit(crID, info.nSiPM_, info.time_,info.eDep_));
               evtEnergy += info.eDep_;
           }
       }

       intInfo.setCaloEnergy(evtEnergy);
       intInfo.setNCaloHits(pulseMap.size());
       intInfo.setNCaphriHits(caphriHitsColl.size());

       if ( diagLevel_ > 0 ) std::cout<<"[CaloHitMakerFast] extracted "<<caloHitsColl.size()<<" CaloDigis"<<std::endl;
       if ( diagLevel_ > 0 ) std::cout<<"[CaloHitMakerFast] extracted "<<caphriHitsColl.size()<<" CapriDigis"<<std::endl;

   }



   //--------------------------------------------------------------------------------------------------------------
   void CaloHitMakerFast::addPulse(pulseMapType& pulseMap, unsigned crystalID, float time, float eDep)
   {

       bool addNewHit(true);
       for (auto& pulse : pulseMap[crystalID])
       {
           if (std::fabs(pulse.time_ - time) > deltaTPulses_) continue;

           float ratio  = (eDep-pulse.eDep_)/(eDep+pulse.eDep_);
           float eMean  = (eDep+pulse.eDep_)/2.0;
           float sigmaR = 0.707*sqrt(1.0/eMean/nPEperMeV_ + noise2_/eMean/eMean);

           if (abs(ratio) > nSigmaNoise_*sigmaR) continue;

           pulse.add(time,eDep);
           addNewHit = false;
           break;
       }

       if (addNewHit) pulseMap[crystalID].push_back(HitInfo(time, eDep));
   }
}

DEFINE_ART_MODULE(mu2e::CaloHitMakerFast)
