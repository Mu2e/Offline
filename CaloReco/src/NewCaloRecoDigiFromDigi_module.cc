#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/LogNormalProcessor.hh"
#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/RawProcessor.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/NewCaloDigi.hh"
#include "RecoDataProducts/inc/NewCaloDigiCollection.hh"
#include "RecoDataProducts/inc/NewCaloRecoDigi.hh"
#include "RecoDataProducts/inc/NewCaloRecoDigiCollection.hh"

#include <iostream>
#include <string>


#include <iostream>
#include <string>


namespace mu2e {

  class NewCaloRecoDigiFromDigi : public art::EDProducer {

  public:

  
    explicit NewCaloRecoDigiFromDigi(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloDigisToken_{consumes<NewCaloDigiCollection>(pset.get<std::string>("caloDigiModuleLabel"))},
      digiSampling_        (pset.get<double>     ("digiSampling")),
      maxChi2Cut_          (pset.get<double>     ("maxChi2Cut")),
      diagLevel_           (pset.get<int>        ("diagLevel",0)),
      nplot_(0)
    {
      produces<NewCaloRecoDigiCollection>();
      auto const& param = pset.get<fhicl::ParameterSet>("RawProcessor", {});
      waveformProcessor_ = std::make_unique<RawProcessor>(param);
       
    }

    void beginRun(art::Run& aRun) override;
    void produce(art::Event& e) override;

  private:

    art::ProductToken<NewCaloDigiCollection> const caloDigisToken_;
    double const digiSampling_;
    double       maxChi2Cut_;
    int          diagLevel_;
    int          nplot_;

    std::unique_ptr<WaveformProcessor> waveformProcessor_;

    void extractRecoDigi(art::ValidHandle<NewCaloDigiCollection> const& caloDigis,
                         NewCaloRecoDigiCollection& recoCaloHits);

  };


  //-------------------------------------------------------
  void NewCaloRecoDigiFromDigi::produce(art::Event& event)
  {

    if (diagLevel_ > 0) std::cout<<"[NewCaloRecoDigiFromDigi::produce] begin"<<std::endl;

    auto const& caloDigisH = event.getValidHandle(caloDigisToken_);

    auto recoCaloDigiColl = std::make_unique<NewCaloRecoDigiCollection>();
    extractRecoDigi(caloDigisH, *recoCaloDigiColl);

    if ( diagLevel_ > 3 )
      {
        printf("[NewCaloRecoDigiFromDigi::produce] produced RecoCrystalHits ");
        printf(", recoCaloDigiColl size  = %i \n", int(recoCaloDigiColl->size()));
      }

    event.put(std::move(recoCaloDigiColl));

    if (diagLevel_ > 0) std::cout<<"[NewCaloRecoDigiFromDigi::produce] end"<<std::endl;

    return;
  }

  //-----------------------------------------------------------------------------
  void NewCaloRecoDigiFromDigi::beginRun(art::Run& aRun)
  {
    waveformProcessor_->initialize();
  }

  //--------------------------------------------------------------------------------------
  void NewCaloRecoDigiFromDigi::extractRecoDigi(art::ValidHandle<NewCaloDigiCollection> const& caloDigisHandle, NewCaloRecoDigiCollection &recoCaloHits)
  {

    std::vector<double> x,y;
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

    auto const& caloDigis = *caloDigisHandle;
    NewCaloDigi const* base = &caloDigis.front(); // What if caloDigis is empty?
    
    for (const auto& caloDigi : caloDigis)
      {
	uint16_t errFlag = caloDigi.errorFlag();
        if(errFlag){ continue; } ///NOTE: Added for new format! Function: to skip if flagged!
        int    roId     = caloDigi.roId();
        double t0       = caloDigi.t0();
        //float peak	=caloDigi.peakpos(); //TODO
        //uint8_t eventMode = caloDigi.eventMode();
        double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
        ///TODO - use this event mode to chose calibration settings!

        const std::vector<int>& waveform = caloDigi.waveform();

        size_t index = &caloDigi - base;
        art::Ptr<NewCaloDigi> caloDigiPtr(caloDigisHandle, index);

        x.clear();
        y.clear();
        for (unsigned int i=0;i<waveform.size();++i)
          {
            x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
            y.push_back(waveform.at(i));
          }

        if (diagLevel_ > 3)
          {
            std::cout<<"[NewCaloRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
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
                std::cout<<"[NewCaloRecoDigiFromDigi::extractAmplitude] extract "<<roId<<"   i="<<i<<"  eDep="<<eDep<<" time="<<time<<"  chi2="<<chi2<<std::endl;
              }

            if (chi2/ndf > maxChi2Cut_) continue;

            recoCaloHits.emplace_back(NewCaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,chi2,ndf,isPileUp));
          }

      }

  }

}

DEFINE_ART_MODULE(mu2e::NewCaloRecoDigiFromDigi);
