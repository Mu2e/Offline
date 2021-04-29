#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/ProtonBunchTime.hh"

#include <iostream>
#include <string>
   

namespace mu2e {

  class FastCaloRecoDigiFromDigi : public art::EDProducer {

      public:
          struct Config 
          {
              using Name    = fhicl::Name;
              using Comment = fhicl::Comment;
              fhicl::Atom<art::InputTag> caloShowerCollection  { Name("caloDigiCollection"),     Comment("CaloDigi collection name") }; 
              fhicl::Atom<art::InputTag> ProtonBunchTime       { Name("ProtonBunchTime"),        Comment("Proton bunch time collection name") }; 
              fhicl::Atom<double>        digiSampling          { Name("digiSampling"),           Comment("Digitization time sampling") }; 
              fhicl::Atom<double>        minDigiE              { Name("minDigiE"),               Comment("Minimum hit energy to digitize") }; 
              fhicl::Atom<double>        shiftTime             { Name("shiftTime"),              Comment("Time correction between peak and start time ") }; 
              fhicl::Atom<int>           diagLevel             { Name("diagLevel"),              Comment("Diag Level"),0 };
          };
          
          explicit FastCaloRecoDigiFromDigi(const art::EDProducer::Table<Config>& config) :
             EDProducer{config},
             caloDigisToken_    {consumes<CaloDigiCollection>(config().caloShowerCollection())},
             pbtToken_          {consumes<ProtonBunchTime>   (config().ProtonBunchTime())},
             digiSampling_      (config().digiSampling()),
             minDigiE_          (config().minDigiE()),
             shiftTime_         (config().shiftTime()),
             diagLevel_         (config().diagLevel())
          {
              produces<CaloRecoDigiCollection>();
          }


          void produce(art::Event& e) override;


      private:
          void extractRecoDigi(const art::ValidHandle<CaloDigiCollection>& caloDigis, CaloRecoDigiCollection& recoCaloHits, double pbtime);

          art::ProductToken<CaloDigiCollection> caloDigisToken_;
          art::ProductToken<ProtonBunchTime>    pbtToken_;
          double digiSampling_;
          double minDigiE_ ;
          double shiftTime_ ;
          int    diagLevel_;
      };



      void FastCaloRecoDigiFromDigi::produce(art::Event& event)
      {
          if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi::produce] begin"<<std::endl;

          const auto& caloDigisH = event.getValidHandle(caloDigisToken_);
          auto recoCaloDigiColl = std::make_unique<CaloRecoDigiCollection>();
          
          const auto& pbtH = event.getValidHandle(pbtToken_);
          const ProtonBunchTime& pbt(*pbtH);
          double pbtime = pbt.pbtime_;

          extractRecoDigi(caloDigisH, *recoCaloDigiColl, pbtime);
          if ( diagLevel_ > 0 )std::cout<<"[FastCaloRecoDigiFromDigi::produce] extracted "<<recoCaloDigiColl->size()<<" RecoDigis"<<std::endl;

          event.put(std::move(recoCaloDigiColl));
          if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi::produce] end"<<std::endl;

          return;
      }

      void FastCaloRecoDigiFromDigi::extractRecoDigi(const art::ValidHandle<CaloDigiCollection>& caloDigisHandle, CaloRecoDigiCollection &recoCaloHits, double protonBunchTime)
      {
          ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
          auto const& caloDigis = *caloDigisHandle;
          CaloDigi const* base  = &caloDigis.front(); 

          for (const auto& caloDigi : caloDigis)
          {
              int    roId     = caloDigi.SiPMID();
              double t0       = caloDigi.t0();
                      // TODO:+ calorimeterCalibrations->timeOffset(roId);
              double adc2MeV  = calorimeterCalibrations->Peak2MeV(roId);
              size_t index = &caloDigi - base;
              art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

              double time = t0 + (caloDigi.peakpos()+0.5)*digiSampling_ - protonBunchTime - shiftTime_; 
              double eDep = (caloDigi.waveform().at(caloDigi.peakpos()))*adc2MeV; 
              if (eDep <  minDigiE_) continue;
              double eDepErr = 0*adc2MeV;
              double timeErr = 0;

              if (diagLevel_ > 1)
                  std::cout<<"[FastRecoDigiFromDigi::extractAmplitude] extracted Digi with roId =  "<<roId<<"  eDep = "<<eDep<<" time = " <<time<<std::endl;

              recoCaloHits.emplace_back(CaloRecoDigi(caloDigiPtr, eDep,eDepErr,time,timeErr,0,1,false));
          }
      }
}

DEFINE_ART_MODULE(mu2e::FastCaloRecoDigiFromDigi);
