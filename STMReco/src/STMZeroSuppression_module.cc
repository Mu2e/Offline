//
// Create zero-suppressed STMWaveformDigis from unsuppressed STMWaveformDigis
// Original authors: Claudia Alvarez-Garcia, Alex Keshavarzi, and Mark Lancaster (see DocDB-43057 for details)
// Adapted for Offline: Andy Edmonds
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include <utility>
#include <algorithm>

// root
#include "TH1F.h"
#include "TF1.h"

#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"


namespace mu2e {

  class STMZeroSuppression : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"), Comment("InputTag for STMWaveformDigiCollection")};
        fhicl::Atom<double> tbefore{ Name("tbefore"), Comment("Store this time before the peak [ns]")};
        fhicl::Atom<double> tafter{ Name("tafter"), Comment("Store this time after the peak [ns]")};
        fhicl::Atom<double> threshold{ Name("threshold"), Comment("Threshold to define the peak [ADC/ct]")};
        fhicl::Atom<double> window{ Name("window"), Comment("Calculate the gradient between ADC values this number of elements away from each other")};
        fhicl::Atom<double> naverage{ Name("naverage"), Comment("Number of ADC values to average the gradient over")};
        fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMZeroSuppression(const Parameters& conf);

    private:
    void beginJob() override;
    void produce(art::Event& e) override;

    void calculateGradient(const std::vector<int16_t>& adcs);
    void averageGradient();
    void findPeaks(const STMEnergyCalib& stmEnergyCalib);
    void chooseStartsAndEnds(); // taking into account any overlapping data

    int _verbosityLevel;
    art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
    STMChannel _channel;

    double _tbefore; // time before the peak [ns]
    double _tafter; // time after the peak [ns]
    double _threshold; // threshold

    unsigned long int _nadc; // number of ADC values in unsuppressed waveform
    unsigned long int _window; // distance between two ADC values to calculate the gradient for
    int _naverage; // number of ADC values to average the gradient over
    std::vector<int16_t> _gradient; // vector to store gradient (difference between consecutive ADC values)
    std::vector<double> _avgradient; // vector to store averaged gradient
    std::vector<double> _avtime; // vector to store times at averaged gradient points
    std::vector<unsigned long int> _peaks; // vector with each peak time in clock ticks
    std::vector<size_t> _starts; // start positions of each zero-suppressed waveform
    std::vector<size_t> _ends; // end positions of each zero-suppressed waveform
    std::vector<size_t> _finalstarts; // start positions of each zero-suppressed waveform after taking into account overlapping data
    std::vector<size_t> _finalends; // end positions of each zero-suppressed waveform after taking into account overlapping data
  };

  STMZeroSuppression::STMZeroSuppression(const Parameters& config )  :
    art::EDProducer{config}
    ,_verbosityLevel(config().verbosityLevel())
    ,_stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(config().stmWaveformDigisTag()))
    ,_channel(STMUtils::getChannel(config().stmWaveformDigisTag()))
    ,_tbefore(config().tbefore())
    ,_tafter(config().tafter())
    ,_threshold(config().threshold())
    ,_window(config().window())
    ,_naverage(config().naverage())
  {
    produces<STMWaveformDigiCollection>();
  }

  void STMZeroSuppression::beginJob() {
  }

  void STMZeroSuppression::produce(art::Event& event) {
    // create output
    auto waveformsHandle = event.getValidHandle(_stmWaveformDigisToken);
    std::unique_ptr<STMWaveformDigiCollection> outputSTMWaveformDigis(new STMWaveformDigiCollection());


    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get prodition

    if (_verbosityLevel > 0) {
      std::cout << "STM Zero-Suppression Algorithm Parameters:" << std::endl;
      std::cout << "\ttbefore = " << _tbefore << " ns" << std::endl;
      std::cout << "\ttafter = " << _tafter << " ns" << std::endl;
      std::cout << "\tthreshold = " << _threshold << std::endl;
    }

    for (const auto& waveform : *waveformsHandle) {
      const auto& adcs = waveform.adcs();
      _nadc = adcs.size();
      _gradient.clear();
      _avgradient.clear();
      _avtime.clear();
      _peaks.clear();

      _gradient.reserve(_nadc-_window);
      _avgradient.reserve(_nadc-_window);
      _avtime.reserve(_nadc-_window);
      _peaks.reserve(_nadc-_window);

      calculateGradient(adcs);
      averageGradient();
      if(_verbosityLevel){std::cout << _channel.name() << " " << _channel.id() << std::endl;}
      findPeaks(stmEnergyCalib); // pass it the prodition because it needs to get the sampling frequency
      chooseStartsAndEnds();

      const auto& n_zp_waveforms = _finalstarts.size();
      for (size_t i_zp_waveform = 0; i_zp_waveform < n_zp_waveforms; ++i_zp_waveform) {
        const auto& i_start = _finalstarts.at(i_zp_waveform);
        const auto& i_end = _finalends.at(i_zp_waveform);
        std::vector<int16_t> zp_adcs(waveform.adcs().begin()+i_start, waveform.adcs().begin()+i_end);
        STMWaveformDigi stm_waveform(waveform.trigTimeOffset()+i_start, zp_adcs);
        outputSTMWaveformDigis->push_back(stm_waveform);
      }
    }

    if (_verbosityLevel > 0) {
      std::cout << _channel.name() << ": " << outputSTMWaveformDigis->size() << " waveforms found" << std::endl;
    }
    event.put(std::move(outputSTMWaveformDigis));
  }

  void STMZeroSuppression::calculateGradient(const std::vector<int16_t>& adcs) {
    //Calculate the time and gradient vectors
    for(unsigned long int i=0;i<_nadc-_window;i++){
      _gradient.push_back(adcs[i+_window]-adcs[i]);
    }
  }

  void STMZeroSuppression::averageGradient() {
    //Calculate the average of the gradient each _naverage ADC values to avoid fluctuations
    unsigned long int j=0;
    double av_gradient=0;
    unsigned long int h=0; // will be the number of elements in the averaged gradient vector

    unsigned long int n_gradient_points = _gradient.size(); // = _nadc - _window
    while(j<n_gradient_points){
      //Each point of the gradient and ADCtime averaged with the (_naverage-1) following points, each point is the mean of _naverage points of the gradient
      if((j+_naverage)>(n_gradient_points)){_naverage= (n_gradient_points)-j;}
      av_gradient=0;
      for(int k=0;k<_naverage;k++){
        av_gradient= av_gradient+_gradient[j+k];
      }
      _avgradient.push_back(av_gradient/_naverage);
      _avtime.push_back(j);

      j=j+_naverage;
      h++;
    }
  }

  void STMZeroSuppression::findPeaks(const STMEnergyCalib& stmEnergyCalib) {
    //Initial values
    bool found_peak=false;
    unsigned int peak=0;
    unsigned int peakcounter=0;
    unsigned long int n_avgradient_points = _avgradient.size();
    _starts.clear();
    _ends.clear();
    unsigned int nadcBefore = STMUtils::convertToClockTicks(_tbefore, _channel, stmEnergyCalib); // number of samples before peak
    unsigned int nadcAfter = STMUtils::convertToClockTicks(_tafter, _channel, stmEnergyCalib); // number of samples after peak
    //Store positions in clock ticks for the peaks found, fill _peaks[]
    for(unsigned long int i=0;i<n_avgradient_points;i++){

      if(_avgradient[i]>_threshold){
        found_peak=false;
        continue;
      }
      //skip the rest indexes of the peak after the peak that has already been stored
      if((_avgradient[i]<_threshold)&&(found_peak==true)){
        continue;
      }

      if((_avgradient[i]<_threshold)&&(found_peak==false)){
        found_peak=true;
        peak=_avtime[i];
        _peaks[peakcounter]=peak;
        if (peak<nadcBefore) {
          _starts.push_back(0); // too close to the start of the waveform so can't go tbefore back
        }
        else {
          _starts.push_back(peak - nadcBefore);
        }
        if (peak>_nadc-nadcAfter) {
          _ends.push_back(_nadc); // too close to the end of the waveform so can't go tafter forwaed
        }
        else {
          _ends.push_back(peak + nadcAfter);
        }

        peakcounter++;
      }
    }
  }

  void STMZeroSuppression::chooseStartsAndEnds() {
    // Now go through and account for overlapped data
    if (_starts.size() != 0) { // need to be careful just in case there were no peaks found (i.e. just noise)
      unsigned int i_peak = 0;
      size_t current_start = _starts.at(i_peak);
      size_t current_end = _ends.at(i_peak);

      unsigned int peakcounter = _starts.size();
      _finalstarts.clear();
      _finalends.clear();
      _finalstarts.reserve(peakcounter);
      _finalends.reserve(peakcounter);
      while (i_peak < peakcounter) {
        if (i_peak == (peakcounter-1)) { // the final peak
          _finalstarts.push_back(current_start);
          _finalends.push_back(_ends.at(i_peak));
          break;
        }
        else if (_starts.at(i_peak+1) <= _ends.at(i_peak)) { // peaks are overlapping
          current_end = _ends.at(i_peak+1); // update so that we will go to the end of the second peak
          ++i_peak;
        }
        else if (_starts.at(i_peak+1) > _ends.at(i_peak)) { // peaks don't overlap
          _finalstarts.push_back(current_start);
          _finalends.push_back(current_end);

          // go to next peak
          current_start = _starts.at(i_peak+1);
          current_end = _ends.at(i_peak+1);
          ++i_peak;
        }
      }
    }
  }
}

DEFINE_ART_MODULE(mu2e::STMZeroSuppression)
