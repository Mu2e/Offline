//
// An EDProducer Module that reads Caloshower and produces CaloDigiPacked
//
// Fill the readout response from the deosited energy. Photo-statistic fluctuations are calculated independently
// for each sensor and must be done at this stage.
//
// The output is split between the different digitization boards
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CaloMC/inc/CaloPulseShape.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>




namespace mu2e {


  class CaloDigiFromShower : public art::EDProducer {
  public:

    explicit CaloDigiFromShower(fhicl::ParameterSet const& pset) :
      EDProducer{pset},
      // Parameters
      caloShowerToken_{consumes<CaloShowerStepROCollection>(pset.get<std::string>("caloShowerModuleLabel"))},
      blindTime_             (pset.get<double>     ("blindTime")),
      addNoise_              (pset.get<bool>       ("addNoise")),           //flag for adding or not Gaussian noise
      noise_                 (pset.get<double>     ("noise")),           // mV
      thresholdVoltage_      (pset.get<double>     ("thresholdVoltage")), // mV
      thresholdAmplitude_    (pset.get<double>     ("thresholdAmplitude")), //mV
      energyScale_           (pset.get<double>     ("energyScale")),       // mV/MeV
      digiSampling_          (pset.get<double>     ("digiSampling")),       // ns
      nBits_                 (pset.get<int>        ("nBits")),       // number of digitizer bits
      dynamicRange_          (pset.get<double>     ("dynamicRange")),       // mV
      endTimeBuffer_         (pset.get<double>     ("endTimeBuffer")),  // ns
      bufferDigi_            (pset.get<int>        ("bufferDigi")),  //# timestamps
      pulseIntegralSteps_    (pset.get<int>        ("pulseIntegralSteps")),         //# integral steps
      diagLevel_             (pset.get<int>        ("diagLevel",0)),
      engine_                (createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      randGauss_             (engine_),
      pulseShape_(CaloPulseShape(digiSampling_,pulseIntegralSteps_))
    {
      produces<CaloDigiCollection>();

      maxADCCounts_ =  1 << nBits_;
      ADCTomV_      = dynamicRange_/float(maxADCCounts_);
      mVToADC_      = float(maxADCCounts_)/dynamicRange_;
      nROperCard_   = 40;
      nBinsPeak_    = pset.get<size_t>("nBinsPeak");
    }

  private:

    void produce(art::Event& e) override;
    void beginRun(art::Run& aRun) override;

    art::ProductToken<CaloShowerStepROCollection> const caloShowerToken_;

    double                  blindTime_;
    double                  mbtime_;

    bool                    addNoise_;
    double                  noise_;

    double                  thresholdVoltage_;
    int                     thresholdAmplitude_;
    double                  energyScale_;
    double                  digiSampling_;
    int                     nBits_;
    double                  dynamicRange_;
    double                  endTimeBuffer_;
    int                     bufferDigi_;
    int                     pulseIntegralSteps_;

    int                     diagLevel_;
    CLHEP::HepRandomEngine& engine_;
    CLHEP::RandGaussQ       randGauss_;
    CaloPulseShape          pulseShape_;

    int                     maxADCCounts_;
    double                  ADCTomV_;
    double                  mVToADC_;
    const Calorimeter*      calorimeter_;
    int                     nROperCard_;
    
    size_t  		    nBinsPeak_;
    std::vector< std::vector<double> > pulseDigitized_;
    std::vector< std::vector<double> > waveforms_;


    void   resetWaveforms();
    void   makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs, CaloDigiCollection&);
    void   fillWaveforms(const CaloShowerStepROCollection& caloShowerStepROs);
    void   readoutResponse(int ROID, double energyCorr, double time);
    void   buildOutputDigi(CaloDigiCollection& caloDigiColl);
    void   diag0(int iRO,std::vector<double>& itWave );
    void   diag1(int iRO, double time, std::vector<int>& wf );
  };


  //-----------------------------------------------------------------------------
  void CaloDigiFromShower::beginRun(art::Run& aRun)
  {
    pulseShape_.buildShapes();
    if ( diagLevel_ > 3) pulseShape_.printShape();
  }



  //---------------------------------------------------------
  void CaloDigiFromShower::produce(art::Event& event)
  {
    if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiFromShower::produce] begin" << std::endl;

    //update condition cache
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    mbtime_ = accPar->deBuncherPeriod;

    auto const& caloShowerStepROs = *event.getValidHandle(caloShowerToken_);
    auto caloDigiColl = std::make_unique<CaloDigiCollection>();
    makeDigitization(caloShowerStepROs, *caloDigiColl);

    event.put(std::move(caloDigiColl));

    if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiFromShower::produce] end" << std::endl;
  }

  //-------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs,CaloDigiCollection& caloDigiColl)
  {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    calorimeter_ = ch.get();

    resetWaveforms();
    fillWaveforms(caloShowerStepROs);
    buildOutputDigi(caloDigiColl);
  }



  //-------------------------------------------------
  void CaloDigiFromShower::resetWaveforms()
  {
    unsigned int nWaveforms   = calorimeter_->nCrystal()*calorimeter_->caloInfo().nROPerCrystal();
    unsigned int waveformSize = (mbtime_ - blindTime_ + endTimeBuffer_) / digiSampling_;

    if (waveforms_.size() < nWaveforms || (waveforms_.size() && waveforms_[0].size() < waveformSize))
      {
        waveforms_.clear();
        for (unsigned int i=0;i< nWaveforms;++i) waveforms_.push_back(std::vector<double>(waveformSize,0.0));
      }


    //if we want noise, fill all waveforms with noise, noise_ is in mV, sample is in counts, noise must be positive
    //otherwise simply reset the waveform to zero

    if (addNoise_)
      for (unsigned int i=0; i<nWaveforms; ++i) std::generate(waveforms_[i].begin(),waveforms_[i].end(),
                                                              [&] {return  std::max(0.0,randGauss_.fire(0.0,noise_)*mVToADC_);});
    else
      for (unsigned int i=0; i<nWaveforms; ++i) std::fill(waveforms_[i].begin(),waveforms_[i].end(),0);

  }


  //----------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::fillWaveforms(const CaloShowerStepROCollection& caloShowerStepROs)
  {
    double totalEdep(0);

    for (const auto& caloShowerStepRO : caloShowerStepROs)
      {
        readoutResponse(caloShowerStepRO.ROID(),caloShowerStepRO.energy(),caloShowerStepRO.time());
        if (diagLevel_ > 0) totalEdep +=  caloShowerStepRO.energy();
      }

    if (diagLevel_ > 1) std::cout<<"[CaloDigiFromShower::fillWaveforms] total energy processed "
                                 <<totalEdep/calorimeter_->caloInfo().nROPerCrystal()<<std::endl;
  }


  //-----------------------------------------------------------------------------------------
  void CaloDigiFromShower::readoutResponse(int ROID, double energyCorr, double time)
  {
    double                     pulseAmp       = energyCorr*energyScale_;
    double                     timeCorr       = time - blindTime_;
    int                        startSample    = timeCorr/digiSampling_;
    int                        precisionIndex = (timeCorr/digiSampling_ - int(timeCorr/digiSampling_))*pulseIntegralSteps_;
    std::vector<double>&       waveform       = waveforms_.at(ROID);
    const std::vector<double>& pulse          = pulseShape_.pulseDigitized(precisionIndex);
    int                        stopSample     = std::min(startSample+pulse.size(), waveform.size());

    for (int timeSample = startSample; timeSample < stopSample; ++timeSample)
      {
        int     funcIdx   = timeSample - startSample;
        double  funcValue = pulse.at(funcIdx);
        double  wfAmp     = pulseAmp*funcValue;
        double  ADCCounts = wfAmp*mVToADC_;

        waveform.at(timeSample) += ADCCounts;
        if ( waveform.at(timeSample) > maxADCCounts_) waveform.at(timeSample) = maxADCCounts_;
      }

  }


  //----------------------------------------------------------------------------
  void CaloDigiFromShower::buildOutputDigi(CaloDigiCollection& caloDigiColl)
  {

    for (unsigned int iRO=0; iRO<waveforms_.size(); ++iRO)
      {

        if (diagLevel_ > 5) std::cout<<"wfContent content (timesample: waveContent, funcValue)"<<std::endl;
        if (diagLevel_ > 4) diag0(iRO,waveforms_.at(iRO));


        std::vector<double>& itWave = waveforms_.at(iRO);

        int waveSize = itWave.size();
        int timeSample(0);
        while (timeSample < waveSize)
          {
            double waveContent = itWave.at(timeSample);
            double funcValue   = waveContent*ADCTomV_;

            if (diagLevel_ > 5 && waveContent > 0) printf("wfContent (%4i:  %4i, %9.3f) \n", timeSample, int(waveContent), funcValue);
            if (funcValue < thresholdVoltage_) {++timeSample; continue;}


            // find the starting / stopping point of the peak
            // the stopping point is the first value below the threshold _and_ the buffer is also below the threshold

            int sampleStart = std::max(timeSample - bufferDigi_,0);
            int sampleStop  = timeSample;
            for (; sampleStop < waveSize; ++sampleStop)
              {
                int sampleCheck = std::min(sampleStop+bufferDigi_+1,waveSize-1);
                double waveOverBuffer = *std::max_element(&itWave.at(sampleStop),&itWave.at(sampleCheck));
                if (waveOverBuffer*ADCTomV_ < thresholdVoltage_) break;
              }
            sampleStop = std::min(sampleStop + bufferDigi_, waveSize-1);

            timeSample = sampleStop+1;  //forward the scanning time


            if (sampleStop == sampleStart) continue;  //check if peak is acceptable and digitize

		
		float wfInt(0);
		size_t peakP(0);
		for(size_t i =sampleStart; i<sampleStop-nBinsPeak_;++i){

			float sum(0);
			for(size_t j=0; j< nBinsPeak_;++j){
				sum+=itWave.at(i+j);
			}
			if(sum>wfInt){
				wfInt = sum;
				peakP = i+nBinsPeak_/2;
			}


		}
			
		double sampleMax = itWave.at(peakP);
		peakP = peakP - sampleStart;
            	if (sampleMax*ADCTomV_ < thresholdAmplitude_) continue;


            int t0 = int(sampleStart*digiSampling_+ blindTime_);
            std::vector<int> wf;
            for (int i=sampleStart; i<=sampleStop; ++i) wf.push_back(int(itWave.at(i)));

            caloDigiColl.emplace_back( CaloDigi(iRO,t0,wf, peakP) );

            if (diagLevel_ > 4) diag1(iRO,t0,wf);
          }
      }
  }




  void CaloDigiFromShower::diag0(int iRO,std::vector<double>& itWave )
  {
    if (*std::max_element(itWave.begin(),itWave.end())<1) return;
    std::cout<<"CaloDigiFromShower::fillOutoutRO] Waveform content for readout "<<iRO<<std::endl;
    for (const auto  &v : itWave) std::cout<<v<<" ";
    std::cout<<std::endl;
  }

  void CaloDigiFromShower::diag1(int iRO, double time, std::vector<int>& wf )
  {
    std::cout<<"Created caloDigi with roID = "<<iRO<<"  time="<<time<<" and content ";
    for (const auto  &v : wf) std::cout<<v<<" ";
    std::cout<<std::endl;
  }

}

DEFINE_ART_MODULE(mu2e::CaloDigiFromShower);
