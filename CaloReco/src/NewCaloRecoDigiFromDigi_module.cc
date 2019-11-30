//Author: S Middleton
//Date: Nov 2019
//Purpose: To make CaloRecoHits from the Digis using the new data product structure (for Online compatib.)

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
      windowPeak_         (pset.get<int>    ("windowPeak")),
      minPeakAmplitude_   (pset.get<double> ("minPeakAmplitude")),
      shiftTime_          (pset.get<double> ("shiftTime")),
      scaleFactor_        (pset.get<double> ("scaleFactor")),
      diagLevel_           (pset.get<int>        ("diagLevel",0))
     {
      produces<NewCaloRecoDigiCollection>();

    }

    void beginRun(art::Run& aRun);
    void produce(art::Event& e) override;

  private:

    art::ProductToken<NewCaloDigiCollection> const caloDigisToken_;
    double const digiSampling_;
    double       maxChi2Cut_;
    int windowPeak_ ;
    double minPeakAmplitude_ ;
    double shiftTime_ ;
    double       scaleFactor_;
    int          diagLevel_;
   
   
    void extractRecoDigi(art::ValidHandle<NewCaloDigiCollection> const& caloDigis, NewCaloRecoDigiCollection& recoCaloHits);

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
  {}

  //--------------------------------------------------------------------------------------
  void NewCaloRecoDigiFromDigi::extractRecoDigi(art::ValidHandle<NewCaloDigiCollection> const& caloDigisHandle, NewCaloRecoDigiCollection &recoCaloHits)
  {

    std::vector<double> x,y;
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");//TODO

    auto const& caloDigis = *caloDigisHandle;
    //if(caloDigis.size() == 0){ continue; } TODO
    NewCaloDigi const* base = &caloDigis.front(); 
    
    for (const auto& caloDigi : caloDigis)
      {
    
	uint16_t errFlag = caloDigi.errorFlag();
        if(errFlag){ continue; } 
        
        int    roId     = caloDigi.roId();
        double t0       = caloDigi.t0();
        
       
        //uint8_t eventMode = caloDigi.eventMode();
        //TODO: if(evenMode == is somethings do this){} --> change calib work
        double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
       
        const std::vector<int>& waveform = caloDigi.waveform();

        size_t index = &caloDigi - base;
        art::Ptr<NewCaloDigi> caloDigiPtr(caloDigisHandle, index);

        x.clear();
        y.clear();
        for (unsigned int i=0;i<waveform.size();++i)
          {
	   
            		x.push_back(t0 + (i+0.5)*digiSampling_); //-timeCorrection_ windowPeak?
            		y.push_back(waveform.at(i));
          }

        if (diagLevel_ > 3)
          {
            std::cout<<"[NewCaloRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
            for (auto const& val : waveform) {std::cout<< val<<" ";} std::cout<<std::endl;
          }

           unsigned int nPeaks_ = caloDigi.peakpos().size();
           double chi2   = 0;
           int ndf    = 1;
           bool isPileUp;
           if(nPeaks_ > 1) { isPileUp = true ; } 
	   else { isPileUp = false; }
           for (unsigned int i=0;i<nPeaks_;++i)
           { 
		 
                if(y[caloDigi.peakpos()[i]] <  minPeakAmplitude_) {continue;}
                double eDep= scaleFactor_*y[caloDigi.peakpos()[i]]*adc2MeV;
                double eDepErr =  0*adc2MeV;
                double time =  x[caloDigi.peakpos()[i]] - shiftTime_;
                double timeErr = 0;

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
