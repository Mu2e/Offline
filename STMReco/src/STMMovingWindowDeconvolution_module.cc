//
// Perform MWD algorithm on zero-suppressed digis
// Original authors: Claudia Alvarez-Garcia, Alex Keshavarzi, and Mark Lancaster (see DocDB-XXXXX for details)
// Adapted for Offline: Andy Edmonds
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include <utility>
#include <algorithm>

#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "TH1.h"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class STMMovingWindowDeconvolution : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"), Comment("InputTag for STMWaveformDigiCollection")};
        fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level")};
        fhicl::Atom<double> tau{Name("tau"), Comment("Decay constant of the adcs (used in the deconvolution step) [ns]")};
        fhicl::Atom<double> M{Name("M"), Comment("M parameter (number of samples to differentiate between)")};
        fhicl::Atom<double> L{Name("L"), Comment("L parameter (number of samples to average over)")};
        fhicl::Atom<double> nsigma_cut{Name("nsigma_cut"), Comment("Number of sigma away from baseline_mean to cut (for finding peaks)")};
        fhicl::Atom<double> thresholdgrad{Name("thresholdgrad"), Comment("Threshold on gradient to cut out peaks when calculating baseline")};
        fhicl::OptionalAtom<std::string> xAxis{ Name("xAxis"), Comment("Choice of x-axis unit for histograms if verbosity level >= 5: \"sample_number\", \"adcs_time\", or \"event_time\"") };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMMovingWindowDeconvolution(const Parameters& conf);

    private:
    void beginJob() override;
    void produce(art::Event& e) override;

    void deconvolve(const STMWaveformDigi& adcs, std::vector<double>& deconvolved_data, const STMEnergyCalib& stmEnergyCalib);
    void differentiate(const std::vector<double>& deconvolved_data, std::vector<double>& differentiated_data);
    void average(const std::vector<double>& differentiated_data, std::vector<double>& averaged_data);
    void calculate_baseline(const std::vector<double>& averaged_data, double& mean, double& stddev);
    void find_peaks(const std::vector<double>& averaged_data, std::vector<double>& peak_heights, std::vector<double>& peak_times, const double baseline_mean, const double baseline_stddev);

    void make_debug_histogram(const art::Event& event, int count, const STMWaveformDigi& adcs, const STMEnergyCalib& stmEnergyCalib, const std::vector<double>& deconvolved_data, const std::vector<double>& differentiated_data, const std::vector<double>& averaged_data, const double baseline_mean, const double baseline_stddev, const std::vector<double>& peak_heights, const std::vector<double>& peak_times);

    int _verbosityLevel;
    art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;
    STMChannel _channel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;

    double _tau; // decay time of adcs [ns] (used in deconvolution step)
    double _M; // M-parameter (used in differentiation step)
    double _L; // L-parameter (used in averaging step)
    double _nsigma_cut; // number of sigma away from baseline mean to cut (used in find_peaks)
    double _thresholdgrad; // threshold on gradient

    std::string _xAxis; // optional parameter for x-axis unit if plotting histograms
  };

  STMMovingWindowDeconvolution::STMMovingWindowDeconvolution(const Parameters& config ) :
    art::EDProducer{config}
    ,_verbosityLevel(config().verbosityLevel())
    ,_stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(config().stmWaveformDigisTag()))
    ,_channel(STMUtils::getChannel(config().stmWaveformDigisTag()))
    ,_tau(config().tau())
    ,_M(config().M())
    ,_L(config().L())
    ,_nsigma_cut(config().nsigma_cut())
    ,_thresholdgrad(config().thresholdgrad())
  {
    produces<STMMWDDigiCollection>();

    if (!config().xAxis(_xAxis)) {
      if (_verbosityLevel >= 5) {
        throw cet::exception("STMMovingWindowDecomposition") << "No xAxis scale defined despite requesting verbosity level >= 5" << std::endl;
      }
    }
  }

  void STMMovingWindowDeconvolution::beginJob() {
  }

  void STMMovingWindowDeconvolution::produce(art::Event& event) {
    // create output
    unique_ptr<STMMWDDigiCollection> outputMWDDigis(new STMMWDDigiCollection);
    auto adcsDigisHandle = event.getValidHandle(_stmWaveformDigisToken);

    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get prodition

    std::vector<double> deconvolved_data;
    std::vector<double> differentiated_data;
    std::vector<double> averaged_data;
    int count = 0;
    for (const auto& adcs : *adcsDigisHandle) {

      // clear out data from previous adcs
      deconvolved_data.clear();
      deconvolved_data.reserve(adcs.adcs().size());
      differentiated_data.clear();
      differentiated_data.reserve(adcs.adcs().size());
      averaged_data.clear();
      averaged_data.reserve(adcs.adcs().size());

      deconvolve(adcs, deconvolved_data, stmEnergyCalib);
      differentiate(deconvolved_data, differentiated_data);
      average(differentiated_data, averaged_data);

      double baseline_mean = 0;
      double baseline_stddev = 0;
      calculate_baseline(averaged_data, baseline_mean, baseline_stddev);

      std::vector<double> peak_heights;
      std::vector<double> peak_times;
      find_peaks(averaged_data, peak_heights, peak_times, baseline_mean, baseline_stddev);
      for (size_t i_peak = 0; i_peak < peak_heights.size(); ++i_peak) {
        STMMWDDigi mwd_digi(peak_times[i_peak], -1*peak_heights[i_peak]); // peak_heights are negative, make them positive here
        outputMWDDigis->push_back(mwd_digi);
      }

      if (_verbosityLevel >= 5) {
        make_debug_histogram(event, count, adcs, stmEnergyCalib, deconvolved_data, differentiated_data, averaged_data, baseline_mean, baseline_stddev, peak_heights, peak_times);
      }

      ++count;
    }
    if (_verbosityLevel > 0) {
      std::cout << _channel.name() << ": " << outputMWDDigis->size() << " MWD digis found" << std::endl;
    }
    event.put(std::move(outputMWDDigis));
  }

  void STMMovingWindowDeconvolution::deconvolve(const STMWaveformDigi& adcs, std::vector<double>& deconvolved_data, const STMEnergyCalib& stmEnergyCalib) {
    const auto pedestal = stmEnergyCalib.pedestal(_channel);
    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);
    const auto& input_data = adcs.adcs();
    deconvolved_data.push_back(input_data[0] - pedestal);
    for(size_t i=1; i<input_data.size(); i++){
      deconvolved_data.push_back((input_data[i]-pedestal)-(1-(nsPerCt/_tau))*(input_data[i-1]-pedestal) + deconvolved_data[i-1]);
    }
  }

  void STMMovingWindowDeconvolution::differentiate(const std::vector<double>& deconvolved_data, std::vector<double>& differentiated_data) {
    for (size_t i = 0; i < _M; ++i) {
      differentiated_data.push_back(deconvolved_data[i]);
    }
    for (size_t i = _M; i < deconvolved_data.size(); ++i) {
      differentiated_data.push_back(deconvolved_data[i] - deconvolved_data[i-_M]);
    }
  }

  void STMMovingWindowDeconvolution::average(const std::vector<double>& differentiated_data, std::vector<double>& averaged_data) {

    double sum = 0.;
    // sum the first L-1 elements of differentiated data
    // and set the first L-1 elements of averaged data
    for (size_t i = 0; i < _L-1; ++i) {
      sum += differentiated_data[i];
      averaged_data.push_back(differentiated_data[i]);
    }
    sum += differentiated_data[_L-1];
    averaged_data.push_back(sum/_L);

    for (size_t i = _L; i < differentiated_data.size(); ++i) {
      sum += differentiated_data[i]-differentiated_data[i-_L]; // move the sum across one sample
      averaged_data.push_back(sum/_L);
    }
  }

  void STMMovingWindowDeconvolution::calculate_baseline(const std::vector<double>& averaged_data, double& mean, double& stddev){

    int k = _M;
    int nadc = averaged_data.size();

    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::mean, tag::variance> > acc_data_without_peaks;

    // Remove peaks so that we can calculate the baseline of the averaged data
    while (k < nadc){
      double gradient = averaged_data[k+1] - averaged_data[k];
      if(gradient < _thresholdgrad){ // if the gradient is too sharp (i.e. we have hit a peak)
        k = k + (_M+2*_L); // jump ahead a little bit
        continue;
      }
      else {
        acc_data_without_peaks(averaged_data[k]);
        k++;
      }
    }

    mean = extract_result<tag::mean>(acc_data_without_peaks);
    double variance = extract_result<tag::variance>(acc_data_without_peaks);
    stddev = std::sqrt(variance);
  }

  void STMMovingWindowDeconvolution::find_peaks(const std::vector<double>& averaged_data, std::vector<double>& peak_heights, std::vector<double>& peak_times, const double baseline_mean, const double baseline_stddev) {

    double threshold_cut = baseline_mean - _nsigma_cut*baseline_stddev;

    int n = averaged_data.size();
    double lowest_height = 0;
    int lowest_height_time = -1; // in clock ticks

    for( int i = _M; i < n; i++){
      //      std::cout << "i = " << i << ", avg[i] = " << averaged_data[i] << ": ";
      if (averaged_data[i] < threshold_cut){ // the waveforms are negative so if we go below this threshold we have seen a peak
        //        std::cout << " below threshold, ";
        if ((averaged_data[i] < averaged_data[i-1]) && averaged_data[i] < lowest_height){ // if the current value is lower than the previous value and lower than the lowest value we've seen so far
          lowest_height = averaged_data[i]; // record the lowest height

          if (lowest_height_time == -1) {
            lowest_height_time = i; // record the time we cross the threshold
          }
        }
        else {
          //          std::cout << " not low enough." << std::endl;
          continue;
        }
      }

      if (lowest_height_time == -1) { // this will be true if we haven't seen a peak yet
        //        std::cout << " not found peak yet." << std::endl;
        continue;
      }
      else if (averaged_data[i] > threshold_cut){ // if we have seen a peak and go above the cut
        //        std::cout << " all done.";
        // record the height and time
        peak_heights.push_back(lowest_height - baseline_mean);
        peak_times.push_back(lowest_height_time); // ct

        lowest_height_time=-1; // reset to 0 so we can find a new peak
        lowest_height = 0;
      }
      //      std::cout << std::endl;
    }
  }

  void STMMovingWindowDeconvolution::make_debug_histogram(const art::Event& event, int count, const STMWaveformDigi& adcs, const STMEnergyCalib& stmEnergyCalib, const std::vector<double>& deconvolved_data, const std::vector<double>& differentiated_data, const std::vector<double>& averaged_data, const double baseline_mean, const double baseline_stddev, const std::vector<double>& peak_heights, const std::vector<double>& peak_times) {
    art::ServiceHandle<art::TFileService> tfs;
    std::stringstream histsuffix;
    histsuffix.str("");
    histsuffix << "_evt" << event.event() << "_wvf" << count;

    const auto pedestal = stmEnergyCalib.pedestal(_channel);
    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);
    Binning binning = STMUtils::getBinning(adcs, _xAxis, nsPerCt);
    TH1D* h_adcs = tfs->make<TH1D>(("h_adcs"+histsuffix.str()).c_str(), "Waveform", binning.nbins(),binning.low(),binning.high());
    TH1D* h_deconvolved = tfs->make<TH1D>(("h_deconvolved"+histsuffix.str()).c_str(), "Deconvolution", binning.nbins(),binning.low(),binning.high());
    TH1D* h_differentiated = tfs->make<TH1D>(("h_differentiated"+histsuffix.str()).c_str(), "Differentiated", binning.nbins(),binning.low(),binning.high());
    TH1D* h_averaged = tfs->make<TH1D>(("h_averaged"+histsuffix.str()).c_str(), "Averaged", binning.nbins(),binning.low(),binning.high());
    TH1D* h_baseline_mean = tfs->make<TH1D>(("h_baseline_mean"+histsuffix.str()).c_str(), "Baseline Mean", binning.nbins(),binning.low(),binning.high());
    TH1D* h_baseline_mean_plus_stddev = tfs->make<TH1D>(("h_baseline_mean_plus_stddev"+histsuffix.str()).c_str(), "Baseline Mean + StdDev", binning.nbins(),binning.low(),binning.high());
    TH1D* h_baseline_mean_minus_stddev = tfs->make<TH1D>(("h_baseline_mean_minus_stddev"+histsuffix.str()).c_str(), "Baseline Mean - StdDev", binning.nbins(),binning.low(),binning.high());
    TH1D* h_peak_threshold = tfs->make<TH1D>(("h_peak_threshold"+histsuffix.str()).c_str(), "Threshold", binning.nbins(),binning.low(),binning.high());

    for (size_t i = 0; i < deconvolved_data.size(); ++i) {
      h_adcs->SetBinContent(i+1, adcs.adcs()[i] - pedestal); // remove the pedestal
      h_deconvolved->SetBinContent(i+1, deconvolved_data[i]);
      h_differentiated->SetBinContent(i+1, differentiated_data[i]);
      h_averaged->SetBinContent(i+1, averaged_data[i]);
      h_baseline_mean->SetBinContent(i+1, baseline_mean);
      h_baseline_mean_plus_stddev->SetBinContent(i+1, baseline_mean + baseline_stddev);
      h_baseline_mean_minus_stddev->SetBinContent(i+1, baseline_mean - baseline_stddev);
      h_peak_threshold->SetBinContent(i+1, baseline_mean - _nsigma_cut*baseline_stddev);
    }
    TH1D* h_peaks = tfs->make<TH1D>(("h_peaks"+histsuffix.str()).c_str(), "Peaks", binning.nbins(),binning.low(),binning.high());
    for (size_t i_peak = 0; i_peak < peak_heights.size(); ++i_peak) {
      //      std::cout << "t = " << peak_times[i_peak] << ", E = " << peak_heights[i_peak] << std::endl;
      h_peaks->SetBinContent(peak_times[i_peak]+1, peak_heights[i_peak]);
    }
  }
}

DEFINE_ART_MODULE(mu2e::STMMovingWindowDeconvolution)
