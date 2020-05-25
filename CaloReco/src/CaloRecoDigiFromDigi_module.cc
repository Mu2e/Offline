#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"
#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/LogNormalProcessor.hh"
#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/RawProcessor.hh"

#include <iostream>
#include <string>


namespace mu2e {

  class CaloRecoDigiFromDigi : public art::EDProducer 
  {
     public:
        enum processorStrategy {NoChoice, RawExtract, LogNormalFit, FixedFast};

        struct Config 
        {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;        

           fhicl::Table<mu2e::RawProcessor::Config>       proc_raw_conf  { Name("RawProcessor"),        Comment("Raw processor config") };
           fhicl::Table<mu2e::LogNormalProcessor::Config> proc_log_conf  { Name("LogNormalProcessor"),  Comment("Fixed fast processor config") };
           fhicl::Table<mu2e::FixedFastProcessor::Config> proc_fixed_conf{ Name("FixedFastProcessor"),  Comment("Log normal fit processor config") };                    
           fhicl::Atom<std::string> caloDigiModuleLabel                  { Name("caloDigiModuleLabel"), Comment("Calo Digi module label") };
           fhicl::Atom<std::string> processorStrategy                    { Name("processorStrategy"),   Comment("Digi reco processor name") };
           fhicl::Atom<double>      digiSampling                         { Name("digiSampling"),        Comment("Calo ADC sampling time (ns)") };
           fhicl::Atom<double>      maxChi2Cut                           { Name("maxChi2Cut"),          Comment("Chi2 cut for keeping reco digi") };
           fhicl::Atom<int>         diagLevel                            { Name("diagLevel"),           Comment("Diagnosis level") };
        };

        explicit CaloRecoDigiFromDigi(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloDigisToken_    {consumes<CaloDigiCollection>(config().caloDigiModuleLabel())},
           processorStrategy_ (config().processorStrategy()),
           digiSampling_      (config().digiSampling()),
           maxChi2Cut_        (config().maxChi2Cut()),
           diagLevel_         (config().diagLevel())
        {
            produces<CaloRecoDigiCollection>();

            std::map<std::string, processorStrategy> spmap;
            spmap["RawExtract"]   = RawExtract;
            spmap["LogNormalFit"] = LogNormalFit;
            spmap["FixedFast"]    = FixedFast;

            switch (spmap[processorStrategy_])
            {
                case RawExtract:
                {
                    waveformProcessor_ = std::make_unique<RawProcessor>(config().proc_raw_conf());
                    break;
                }
                case LogNormalFit:
                {
                    waveformProcessor_ = std::make_unique<LogNormalProcessor>(config().proc_log_conf());
                    break;
                }
                case FixedFast:
                {
                    waveformProcessor_ = std::make_unique<FixedFastProcessor>(config().proc_fixed_conf());
                    break;
                }
                default:
                {
                    throw cet::exception("CATEGORY")<< "Unrecognized processor in CaloHitsFromDigis module";
                }
            }
        }

        void beginRun(art::Run& aRun) override;
        void produce(art::Event& e) override;

     private:
        void extractRecoDigi(const art::ValidHandle<CaloDigiCollection>& caloDigis, CaloRecoDigiCollection& recoCaloHits);


        const  art::ProductToken<CaloDigiCollection> caloDigisToken_;
        const  std::string processorStrategy_;
        double digiSampling_;
        double maxChi2Cut_;
        int    diagLevel_;

        std::unique_ptr<WaveformProcessor> waveformProcessor_;
  };


  //-------------------------------------------------------
  void CaloRecoDigiFromDigi::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout<<"[CaloRecoDigiFromDigi::produce] begin"<<std::endl;

      const auto& caloDigisH = event.getValidHandle(caloDigisToken_);
      auto recoCaloDigiColl  = std::make_unique<CaloRecoDigiCollection>();

      extractRecoDigi(caloDigisH, *recoCaloDigiColl);

      event.put(std::move(recoCaloDigiColl));

      if (diagLevel_ > 0) std::cout<<"[CaloRecoDigiFromDigi::produce] end"<<std::endl;
      return;
  }


  //-----------------------------------------------------------------------------
  void CaloRecoDigiFromDigi::beginRun(art::Run& aRun)
  {
      waveformProcessor_->initialize();
  }


  //--------------------------------------------------------------------------------------
  void CaloRecoDigiFromDigi::extractRecoDigi(const art::ValidHandle<CaloDigiCollection>& caloDigisHandle,
                                             CaloRecoDigiCollection &recoCaloHits)
  {
      const auto& caloDigis = *caloDigisHandle;
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      std::vector<double> x,y;
      for (const auto& caloDigi : caloDigis)
      {
          int    roId     = caloDigi.roId();
          double t0       = caloDigi.t0();
          double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
          const std::vector<int>& waveform = caloDigi.waveform();

          size_t index = &caloDigi - &caloDigis.front();
          art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

          x.clear();y.clear();
          for (unsigned int i=0;i<waveform.size();++i)
          {
              x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
              y.push_back(waveform.at(i));
          }

          if (diagLevel_ > 3)
          {
              std::cout<<"[CaloRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
              for (const auto& val : waveform) {std::cout<< val<<" ";} std::cout<<std::endl;
          }

          waveformProcessor_->reset();
          waveformProcessor_->extract(x,y);

          for (int i=0;i<waveformProcessor_->nPeaks();++i)
          {
              double eDep      = waveformProcessor_->amplitude(i)*adc2MeV;
              double eDepErr   = waveformProcessor_->amplitudeErr(i)*adc2MeV;
              double time      = waveformProcessor_->time(i);
              double timeErr   = waveformProcessor_->timeErr(i);
              bool   isPileUp  = waveformProcessor_->isPileUp(i);
              double chi2      = waveformProcessor_->chi2();
              int    ndf       = waveformProcessor_->ndf();

              if (diagLevel_ > 1) std::cout<<"[CaloRecoDigiFromDigi::extractAmplitude] extract "<<roId<<"   i="<<i<<"  eDep="<<eDep
                                           <<" time="<<time<<"  chi2="<<chi2<<std::endl;

              if (chi2/ndf > maxChi2Cut_) continue;
              recoCaloHits.emplace_back(CaloRecoDigi(roId, caloDigiPtr, eDep, eDepErr, time, timeErr, chi2, ndf, isPileUp));
          }
      }     
  }



}

DEFINE_ART_MODULE(mu2e::CaloRecoDigiFromDigi);
