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
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/CaloReco/inc/CaloWaveformProcessor.hh"
#include "Offline/CaloReco/inc/CaloTemplateWFProcessor.hh"
#include "Offline/CaloReco/inc/CaloRawWFProcessor.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"

#include <iostream>
#include <string>
#include <sstream>


namespace mu2e {

  class CaloRecoDigiMaker : public art::EDProducer
  {
     public:
        enum processorStrategy {NoChoice, RawExtract, Template};

        struct Config
        {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Table<mu2e::CaloRawWFProcessor::Config>      proc_raw_conf       { Name("RawProcessor"),        Comment("Raw processor config") };
           fhicl::Table<mu2e::CaloTemplateWFProcessor::Config> proc_templ_conf     { Name("TemplateProcessor"),   Comment("Log normal fit processor config") };
           fhicl::Atom<art::InputTag>                          caloDigiCollection  { Name("caloDigiCollection"),  Comment("Calo Digi module label") };
           fhicl::Atom<art::InputTag>                          pbttoken            { Name("ProtonBunchTimeTag"),  Comment("ProtonBunchTime producer")};
           fhicl::Atom<std::string>                            processorStrategy   { Name("processorStrategy"),   Comment("Digi reco processor name") };
           fhicl::Atom<double>                                 digiSampling        { Name("digiSampling"),        Comment("Calo ADC sampling time (ns)") };
           fhicl::Atom<double>                                 maxChi2Cut          { Name("maxChi2Cut"),          Comment("Chi2 cut for keeping reco digi") };
           fhicl::Atom<int>                                    maxPlots            { Name("maxPlots"),            Comment("Maximum number of waveform plots") };
           fhicl::Atom<int>                                    diagLevel           { Name("diagLevel"),           Comment("Diagnosis level") };
        };

        explicit CaloRecoDigiMaker(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloDigisToken_    {consumes<CaloDigiCollection>(config().caloDigiCollection())},
           pbttoken_          {consumes<ProtonBunchTime>(config().pbttoken())},
           processorStrategy_ (config().processorStrategy()),
           digiSampling_      (config().digiSampling()),
           maxChi2Cut_        (config().maxChi2Cut()),
           maxPlots_          (config().maxPlots()),
           diagLevel_         (config().diagLevel())
        {
            produces<CaloRecoDigiCollection>();

            std::map<std::string, processorStrategy> spmap;
            spmap["RawExtract"]  = RawExtract;
            spmap["TemplateFit"] = Template;

            switch (spmap[processorStrategy_])
            {
                case RawExtract:
                {
                    waveformProcessor_ = std::make_unique<CaloRawWFProcessor>(config().proc_raw_conf());
                    break;
                }
                case Template:
                {
                    waveformProcessor_ = std::make_unique<CaloTemplateWFProcessor>(config().proc_templ_conf());
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
        void extractRecoDigi(const art::ValidHandle<CaloDigiCollection>&, CaloRecoDigiCollection&, double );

        const  art::ProductToken<CaloDigiCollection> caloDigisToken_;
        const  art::ProductToken<ProtonBunchTime>    pbttoken_;
        const  std::string                           processorStrategy_;
        double                                       digiSampling_;
        double                                       maxChi2Cut_;
        int                                          maxPlots_;
        int                                          diagLevel_;
        std::unique_ptr<CaloWaveformProcessor>       waveformProcessor_;
  };


  //-------------------------------------------------------
  // The t0 time in caloDigi has been corrected to the DR time frame for backward compatibility, so we only need to
  // correct for the jitter between the DR marker and the nearest clock
  void CaloRecoDigiMaker::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout<<"[CaloRecoDigiMaker::produce] begin"<<std::endl;

      auto recoCaloDigiColl  = std::make_unique<CaloRecoDigiCollection>();

      // try {
      //   const auto& caloDigisH = event.getValidHandle(caloDigisToken_);
      //   auto pbtH = event.getValidHandle(pbttoken_);
      //   const ProtonBunchTime& pbt(*pbtH);
      //   double pbtOffset = pbt.pbtime_;
        
      //   extractRecoDigi(caloDigisH, *recoCaloDigiColl, pbtOffset);
      // }
      // catch(...) {
      //   if (diagLevel_ > 0) std::cout<<"[CaloRecoDigiMaker::produce] ERROR: input coll not found" << std::endl;
      // }

      event.put(std::move(recoCaloDigiColl));

      if (diagLevel_ > 0) std::cout<<"[CaloRecoDigiMaker::produce] end"<<std::endl;
  }


  //--------------------------------------------------
  void CaloRecoDigiMaker::beginRun(art::Run& aRun)
  {
      waveformProcessor_->initialize();
  }


  //------------------------------------------------------------------------------------------------------------
  void CaloRecoDigiMaker::extractRecoDigi(const art::ValidHandle<CaloDigiCollection>& caloDigisHandle,
                                          CaloRecoDigiCollection &recoCaloHits, double pbtOffset)
  {

      const auto& caloDigis = *caloDigisHandle;
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      double totEnergyReco(0);
      std::vector<double> x{},y{};
      for (const auto& caloDigi : caloDigis)
      {
          int    SiPMID   = caloDigi.SiPMID();
          double t0       = caloDigi.t0();
          double adc2MeV  = calorimeterCalibrations->ADC2MeV(SiPMID);
          const std::vector<int>& waveform = caloDigi.waveform();

          size_t index = &caloDigi - &caloDigis.front();
          art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

          x.clear();y.clear();
          for (unsigned int i=0;i<waveform.size();++i)
          {
              x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
              y.push_back(waveform.at(i));
          }

          waveformProcessor_->reset();
          waveformProcessor_->extract(x,y);
          if (diagLevel_ > 2) std::cout<<"CaloRecoDigiMaker found "<<waveformProcessor_->nPeaks()<<" peaks for SiPMID="<<SiPMID<<std::endl;

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

              if (SiPMID%2==0) totEnergyReco += eDep;
              recoCaloHits.emplace_back(CaloRecoDigi(caloDigiPtr, eDep, eDepErr, time, timeErr, chi2, ndf, isPileUp));
          }
      }

      if (diagLevel_ > 1) std::cout<<"[CaloRecoDigiMaker] Total energy reco "<<totEnergyReco <<std::endl;
  }


}

DEFINE_ART_MODULE(mu2e::CaloRecoDigiMaker)
