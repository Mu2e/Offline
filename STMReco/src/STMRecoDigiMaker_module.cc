#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/CalorimeterCalibrations.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"
#include "Offline/RecoDataProducts/inc/STMRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/STMReco/inc/STMWaveformProcessor.hh"
#include "Offline/STMReco/inc/STMTemplateWFProcessor.hh"
#include "Offline/STMReco/inc/STMRawWFProcessor.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"

#include <iostream>
#include <string>
#include <sstream>


namespace mu2e {

  class STMRecoDigiMaker : public art::EDProducer
  {
     public:
        enum processorStrategy {NoChoice, RawExtract, Template};

        struct Config
        {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Table<mu2e::STMRawWFProcessor::Config>      proc_raw_conf       { Name("RawProcessor"),        Comment("Raw processor config") };
           fhicl::Table<mu2e::STMTemplateWFProcessor::Config> proc_templ_conf     { Name("TemplateProcessor"),   Comment("Log normal fit processor config") };
           fhicl::Atom<art::InputTag>                          stmDigiCollection  { Name("stmDigiCollection"),  Comment("STM Digi module label") };
           fhicl::Atom<art::InputTag>                          pbttoken            { Name("ProtonBunchTimeTag"),  Comment("ProtonBunchTime producer")};
           fhicl::Atom<std::string>                            processorStrategy   { Name("processorStrategy"),   Comment("Digi reco processor name") };
           fhicl::Atom<double>                                 digiSampling        { Name("digiSampling"),        Comment("STM ADC sampling time (ns)") };
           fhicl::Atom<double>                                 maxChi2Cut          { Name("maxChi2Cut"),          Comment("Chi2 cut for keeping reco digi") };
           fhicl::Atom<int>                                    maxPlots            { Name("maxPlots"),            Comment("Maximum number of adcs plots") };
           fhicl::Atom<int>                                    diagLevel           { Name("diagLevel"),           Comment("Diagnosis level") };
        };

        explicit STMRecoDigiMaker(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           stmDigisToken_    {consumes<STMDigiCollection>(config().stmDigiCollection())},
           pbttoken_          {consumes<ProtonBunchTime>(config().pbttoken())},
           processorStrategy_ (config().processorStrategy()),
           digiSampling_      (config().digiSampling()),
           maxChi2Cut_        (config().maxChi2Cut()),
           maxPlots_          (config().maxPlots()),
           diagLevel_         (config().diagLevel())
        {
            produces<STMRecoDigiCollection>();

            std::map<std::string, processorStrategy> spmap;
            spmap["RawExtract"]  = RawExtract;
            spmap["TemplateFit"] = Template;

            switch (spmap[processorStrategy_])
            {
                case RawExtract:
                {
                    waveformProcessor_ = std::make_unique<STMRawWFProcessor>(config().proc_raw_conf());
                    break;
                }
                case Template:
                {
                    waveformProcessor_ = std::make_unique<STMTemplateWFProcessor>(config().proc_templ_conf());
                    break;
                }
                default:
                {
                    throw cet::exception("CATEGORY")<< "Unrecognized processor in STMHitsFromDigis module";
                }
            }
        }

        void beginRun(art::Run& aRun) override;
        void produce(art::Event& e) override;

     private:
        void extractRecoDigi(const art::ValidHandle<STMDigiCollection>&, STMRecoDigiCollection&, double );

        const  art::ProductToken<STMDigiCollection> stmDigisToken_;
        const  art::ProductToken<ProtonBunchTime>    pbttoken_;
        const  std::string                           processorStrategy_;
        double                                       digiSampling_;
        double                                       maxChi2Cut_;
        int                                          maxPlots_;
        int                                          diagLevel_;
        std::unique_ptr<STMWaveformProcessor>       waveformProcessor_;
  };


  //-------------------------------------------------------
  // The t0 time in stmDigi has been corrected to the DR time frame for backward compatibility, so we only need to
  // correct for the jitter between the DR marker and the nearest clock
  void STMRecoDigiMaker::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout<<"[STMRecoDigiMaker::produce] begin"<<std::endl;

      const auto& stmDigisH = event.getValidHandle(stmDigisToken_);
      auto recoSTMDigiColl  = std::make_unique<STMRecoDigiCollection>();

      auto pbtH = event.getValidHandle(pbttoken_);
      const ProtonBunchTime& pbt(*pbtH);
      double pbtOffset = pbt.pbtime_;

      extractRecoDigi(stmDigisH, *recoSTMDigiColl, pbtOffset);

      event.put(std::move(recoSTMDigiColl));

      if (diagLevel_ > 0) std::cout<<"[STMRecoDigiMaker::produce] end"<<std::endl;
  }


  //--------------------------------------------------
  void STMRecoDigiMaker::beginRun(art::Run& aRun)
  {
      waveformProcessor_->initialize();
  }


  //------------------------------------------------------------------------------------------------------------
  void STMRecoDigiMaker::extractRecoDigi(const art::ValidHandle<STMDigiCollection>& stmDigisHandle,
                                          STMRecoDigiCollection &recoSTMHits, double pbtOffset)
  {

      const auto& stmDigis = *stmDigisHandle;
      ConditionsHandle<STMrimeterCalibrations> stmrimeterCalibrations("ignored");

      double totEnergyReco(0);
      std::vector<double> x{},y{};
      for (const auto& stmDigi : stmDigis)
      {
          int    DetID   = stmDigi.DetID();
          double t0       = stmDigi.t0();
          double adc2MeV  = stmrimeterCalibrations->ADC2MeV(DetID);
          const std::vector<int>& adcs = stmDigi.adcs();

          size_t index = &stmDigi - &stmDigis.front();
          art::Ptr<STMDigi> stmDigiPtr(stmDigisHandle, index);

          x.clear();y.clear();
          for (unsigned int i=0;i<adcs.size();++i)
          {
              x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
              y.push_back(adcs.at(i));
          }

          waveformProcessor_->reset();
          waveformProcessor_->extract(x,y);
          if (diagLevel_ > 2) std::cout<<"STMRecoDigiMaker found "<<waveformProcessor_->nPeaks()<<" peaks for DetID="<<DetID<<std::endl;

          for (int i=0;i<waveformProcessor_->nPeaks();++i)
          {
              double eDep      = waveformProcessor_->amplitude(i)*adc2MeV;
              double eDepErr   = waveformProcessor_->amplitudeErr(i)*adc2MeV;
              double time      = waveformProcessor_->time(i) - pbtOffset; // correct to time since protons
              double timeErr   = waveformProcessor_->timeErr(i);
              bool   isPileUp  = waveformProcessor_->isPileUp(i);
              double chi2      = waveformProcessor_->chi2();
              int    ndf       = waveformProcessor_->ndf();

              if (diagLevel_ > 2) std::cout<<"Found reco digi hit with eDep="<<eDep<<"  time="<<time<<" chi2="<<chi2<<"  ndf="<<ndf<<std::endl;
              if (chi2/float(ndf) > maxChi2Cut_) continue;

              if (DetID%2==0) totEnergyReco += eDep;
              recoSTMHits.emplace_back(STMRecoDigi(stmDigiPtr, eDep, eDepErr, time, timeErr, chi2, ndf, isPileUp));
          }
      }

      if (diagLevel_ > 1) std::cout<<"[STMRecoDigiMaker] Total energy reco "<<totEnergyReco <<std::endl;
  }


}

DEFINE_ART_MODULE(mu2e::STMRecoDigiMaker)
