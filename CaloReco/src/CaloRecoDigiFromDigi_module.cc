#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/LogNormalProcessor.hh"
#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/RawProcessor.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include <iostream>
#include <string>


namespace mu2e {



  class CaloRecoDigiFromDigi : public art::EDProducer {

  public:

    enum processorStrategy {NoChoice, RawExtract, LogNormalFit, FixedFast};

    explicit CaloRecoDigiFromDigi(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloDigisToken_{consumes<CaloDigiCollection>(pset.get<std::string>("caloDigiModuleLabel"))},
      processorStrategy_   (pset.get<std::string>("processorStrategy")),
      digiSampling_        (pset.get<double>     ("digiSampling")),
      maxChi2Cut_          (pset.get<double>     ("maxChi2Cut")),
      diagLevel_           (pset.get<int>        ("diagLevel",0)),
      nplot_(0)
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
            auto const& param = pset.get<fhicl::ParameterSet>("RawProcessor", {});
            waveformProcessor_ = std::make_unique<RawProcessor>(param);
            break;
          }

        case LogNormalFit:
          {
            auto const& param = pset.get<fhicl::ParameterSet>("LogNormalProcessor", {});
            waveformProcessor_ = std::make_unique<LogNormalProcessor>(param);
            break;
          }

        case FixedFast:
          {
            auto const& param = pset.get<fhicl::ParameterSet>("FixedFastProcessor", {});
            waveformProcessor_ = std::make_unique<FixedFastProcessor>(param);
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

    art::ProductToken<CaloDigiCollection> const caloDigisToken_;
    std::string const processorStrategy_;
    double const digiSampling_;
    double       maxChi2Cut_;
    int          diagLevel_;
    int          nplot_;

    std::unique_ptr<WaveformProcessor> waveformProcessor_;

    void extractRecoDigi(art::ValidHandle<CaloDigiCollection> const& caloDigis,
                         CaloRecoDigiCollection& recoCaloHits);

  };


  //-------------------------------------------------------
  void CaloRecoDigiFromDigi::produce(art::Event& event)
  {

    if (diagLevel_ > 0) std::cout<<"[CaloRecoDigiFromDigi::produce] begin"<<std::endl;

    //Get the calorimeter Digis
    auto const& caloDigisH = event.getValidHandle(caloDigisToken_);

    auto recoCaloDigiColl = std::make_unique<CaloRecoDigiCollection>();
    extractRecoDigi(caloDigisH, *recoCaloDigiColl);

    if ( diagLevel_ > 3 )
      {
        printf("[CaloRecoDigiFromDigi::produce] produced RecoCrystalHits ");
        printf(", recoCaloDigiColl size  = %i \n", int(recoCaloDigiColl->size()));
      }

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
  void CaloRecoDigiFromDigi::extractRecoDigi(art::ValidHandle<CaloDigiCollection> const& caloDigisHandle,
                                             CaloRecoDigiCollection &recoCaloHits)
  {

    std::vector<double> x,y;

    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

    auto const& caloDigis = *caloDigisHandle;
    CaloDigi const* base = &caloDigis.front(); // What if caloDigis is empty?


    for (const auto& caloDigi : caloDigis)
      {
        int    roId     = caloDigi.roId();
        double t0       = caloDigi.t0();
        double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
        const std::vector<int>& waveform = caloDigi.waveform();

        size_t index = &caloDigi - base;
        art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

        x.clear();
        y.clear();
        for (unsigned int i=0;i<waveform.size();++i)
          {
            x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
            y.push_back(waveform.at(i));
          }

        if (diagLevel_ > 3)
          {
            std::cout<<"[CaloRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
            for (auto const& val : waveform) {std::cout<< val<<" ";} std::cout<<std::endl;
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

            if (diagLevel_ > 1)
              {
                std::cout<<"[CaloRecoDigiFromDigi::extractAmplitude] extract "<<roId<<"   i="<<i<<"  eDep="<<eDep<<" time="<<time<<"  chi2="<<chi2<<std::endl;
              }

            if (chi2/ndf > maxChi2Cut_) continue;

            recoCaloHits.emplace_back(CaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,chi2,ndf,isPileUp));
          }

      }

  }

}

DEFINE_ART_MODULE(mu2e::CaloRecoDigiFromDigi);
