// Perform MWD algorithm on zero-suppressed digis
// Original authors: Claudia Alvarez-Garcia, Alex Keshavarzi, and Mark Lancaster (see DocDB-XXXXX for details)
// Adapted for Offline: Andy Edmonds
// Adapted: Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <utility>

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

// boost includes
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TH1.h"
#include "TTree.h"


using CLHEP::Hep3Vector;
namespace mu2e {
  class STMMovingWindowDeconvolution : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"), Comment("InputTag for STMWaveformDigiCollection")};
        fhicl::Atom<double> tau{Name("tau"), Comment("Decay constant of the waveform (used in the deconvolution step) [ns]")};
        fhicl::Atom<double> M{Name("M"), Comment("M parameter (number of samples to differentiate between)")};
        fhicl::Atom<double> L{Name("L"), Comment("L parameter (number of samples to average over)")};
        fhicl::Atom<double> nsigma_cut{Name("nsigma_cut"), Comment("Number of sigma away from baseline_mean to cut (for finding peaks)")};
        fhicl::Atom<double> thresholdgrad{Name("thresholdgrad"), Comment("Threshold on gradient to cut out peaks when calculating baseline")};
        fhicl::Atom<double> defaultBaselineMean{Name("defaultBaselineMean"), Comment("Default mean to use for baseline when ZS removes all baseline data, in ADC values")};
        fhicl::Atom<double> defaultBaselineSD{Name("defaultBaselineSD"), Comment("Default standard deviation to use for baseline when ZS removes all baseline data, in ADC values")};
        fhicl::OptionalAtom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level")};
        fhicl::OptionalAtom<std::string> xAxis{ Name("xAxis"), Comment("Choice of x-axis unit for histograms if verbosity level >= 5: \"sample_number\", \"waveform_time\", or \"event_time\"") };
        fhicl::OptionalAtom<bool> makeTTreeMWD{ Name("makeTTreeMWD"), Comment("Controls whether to make the TTree with branches ADC, deconvoluted, differentiated, averaged")};
        fhicl::OptionalAtom<bool> makeTTreeEnergies{ Name("makeTTreeEnergies"), Comment("Controls whether to make the TTree with branches time, E")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMMovingWindowDeconvolution(const Parameters& conf);

    private:
    void beginJob() override;
    void produce(art::Event& e) override;

    void deconvolve();
    void differentiate();
    void average();
    void calculate_baseline();
    void find_peaks();
    void make_debug_histogram(const art::Event& event, int count, const STMWaveformDigi& waveform, const STMEnergyCalib& stmEnergyCalib, const std::vector<double>& deconvolved_data, const std::vector<double>& differentiated_data, const std::vector<double>& averaged_data, const double baseline_mean, const double baseline_stddev, const std::vector<double>& peak_heights, const std::vector<double>& peak_times);

    // fhicl variables
    art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;  // token of required data
    STMChannel channel;                                                   // interface to prodicitions service
    double tau = 0.0;                                                     // decay time of waveform [ns] (used in deconvolution step)
    double M = 0.0;                                                       // M-parameter (used in differentiation step)
    double L = 0.0;                                                       // L-parameter (used in averaging step)
    double nsigma_cut = 0.0;                                              // number of sigma away from baseline mean to cut (used in find_peaks)
    double thresholdgrad = 0.0;                                           // threshold on gradient
    double defaultBaselineMean = 0.0;
    double defaultBaselineSD = 0.0;
    int verbosityLevel = 0;                                               // std cout verbosity level
    std::string _xAxis = "";                                              // optional parameter for x-axis unit if plotting histograms
    bool makeTTreeMWD = false;                                            // controls whether to make MWD process TTree
    bool makeTTreeEnergies = false;                                       // controls whether to make results TTree

    // Proditions service
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;

    // MWD analysis variables
    std::vector<int16_t> ADCs;                // input waveform ADCs
    unsigned long int nADCs = 0;              // size of input data
    unsigned long int i = 0;                  // iterator
    std::vector<double> deconvolved_data;
    float pedestal = 0.0;                     // ADC pedestal
    float nsPerCt = 0.0;                      // ADC time step
    double timeFactor = 0.0;                  // variable part of deconvolve()
    std::vector<double> differentiated_data;
    std::vector<double> averaged_data;
    double sum = 0.0;                         // variable part of avarege()
    int count = 0;                            // counter variable

    // Peak finding variables
    double baseline_mean = 0.0;
    double baseline_stddev = 0.0;
    double gradient = 0.0;            // selects baseline vs peak data in calculate_baseline()
    double threshold_cut = 0.0;       // threshold for peak presence in find_peaks()
    double lowest_height = 0;
    int lowest_height_time = -1; // in clock ticks
    std::vector<double> peak_heights;
    std::vector<double> peak_times;
    size_t nPeaks = 0;                // number of peaks found
    bool foundBaselineData = false;

    // TTree variables
    TTree* ttree;
    int16_t ADC = 0, E = 0;
    uint32_t time = 0;
    double deconvoluted = 0.0, differentiated = 0.0, averaged = 0.0;
    uint eventId = 0;
  };

  STMMovingWindowDeconvolution::STMMovingWindowDeconvolution(const Parameters& conf) :
    art::EDProducer{conf},
    _stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(conf().stmWaveformDigisTag())),
    channel(STMUtils::getChannel(conf().stmWaveformDigisTag())),
    tau(conf().tau()),
    M(conf().M()),
    L(conf().L()),
    nsigma_cut(conf().nsigma_cut()),
    thresholdgrad(conf().thresholdgrad()),
    defaultBaselineMean(conf().defaultBaselineMean()),
    defaultBaselineSD(conf().defaultBaselineSD()) {
      produces<STMMWDDigiCollection>();
      verbosityLevel = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;
      if (verbosityLevel > 10)
        verbosityLevel = 10;
      _xAxis = conf().xAxis() ? *(conf().xAxis()) : "";
      makeTTreeMWD = conf().makeTTreeMWD() ? *(conf().makeTTreeMWD()) : false;
      makeTTreeEnergies = conf().makeTTreeEnergies() ? *(conf().makeTTreeEnergies()) : false;
      if (makeTTreeMWD) {
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "STMMovingWindowDeconvolution pulse finding ttree");
        ttree->Branch("ADC", &ADC, "ADC/S");
        ttree->Branch("deconvoluted", &deconvoluted, "deconvoluted/D");
        ttree->Branch("differentiated", &differentiated, "differentiated/D");
        ttree->Branch("averaged", &averaged, "averaged/D");
        ttree->Branch("eventId", &eventId, "eventId/i");
        ttree->Branch("time", &time, "time/i");
      };
      if (makeTTreeEnergies) {
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "STMMovingWindowDeconvolution pulse height ttree");
        ttree->Branch("time", &time, "time/I");
        ttree->Branch("E", &E, "E/S");
        ttree->Branch("eventId", &eventId, "eventId/i");
        ttree->Branch("time", &time, "time/i");
      };
      if (_xAxis != "") {
        if (verbosityLevel >= 5) {
          throw cet::exception("STMMovingWindowDecomposition") << "No xAxis scale defined despite requesting verbosity level >= 5" << std::endl;
        };
      };
    };

  void STMMovingWindowDeconvolution::beginJob() {
    if (verbosityLevel) {
      std::cout << "STM Moving Window Deconvolution" << std::endl;
      std::cout << "\tAlgorithm parameters" << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "tau"            << tau           << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "M"              << M             << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "L"              << L             << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "nsigma_cut"     << nsigma_cut    << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "thresholdgrad"  << thresholdgrad << std::endl;
      std::cout << "\tChannel: " << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "Name" << channel.name()                  << std::endl;
      std::cout << std::left << "\t\t" << std::setw(15) << "ID"   << static_cast<int>(channel.id())  << std::endl;
      std::cout << std::endl; // buffer line
    };
  };

  void STMMovingWindowDeconvolution::produce(art::Event& event) {
    // create output
    std::unique_ptr<STMMWDDigiCollection> outputMWDDigis(new STMMWDDigiCollection);
    auto waveformDigisHandle = event.getValidHandle(_stmWaveformDigisToken);
    // get prodition
    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id());
    pedestal = stmEnergyCalib.pedestal(channel);
    nsPerCt = stmEnergyCalib.nsPerCt(channel);
    timeFactor = 1 - (nsPerCt / tau);
    count = 0;
    eventId = event.id().event();
    for (STMWaveformDigi waveform : *waveformDigisHandle) {
      // clear out data from previous waveform
      ADCs.clear();
      deconvolved_data.clear();
      differentiated_data.clear();
      averaged_data.clear();
      peak_heights.clear();
      peak_times.clear();
      ADCs = waveform.adcs();
      nADCs = ADCs.size();

      if (M < L)
        throw cet::exception("Configuration", "L (" + std::to_string(L) + ") is greater than M (" + std::to_string(M) + "), reconfigure\n");
      if (M > nADCs)
        M = nADCs;
      if (L > nADCs)
        L = nADCs;

      deconvolve();
      differentiate();
      average();
      calculate_baseline();
      find_peaks();

      nPeaks = peak_heights.size();
      std::cout << "nPeaks: " << nPeaks << std::endl;
      for (i = 0; i < nPeaks; ++i) {
        STMMWDDigi mwd_digi(peak_times[i], -1 * peak_heights[i]); // peak_heights are negative, make them positive here
        outputMWDDigis->push_back(mwd_digi);
        if (makeTTreeEnergies) {
          time = mwd_digi.time();
          E = mwd_digi.energy();
          ttree->Fill();
        };
        if (verbosityLevel > 2)
          std::cout << "energy: " << mwd_digi.energy() << std::endl;
      };

      // Save data to TTree
      if (makeTTreeMWD) {
        time = waveform.trigTimeOffset();
        for (i = 0; i < nADCs; i++) {
          ADC = ADCs[i];
          deconvoluted =  deconvolved_data[i];
          differentiated = differentiated_data[i];
          averaged = averaged_data[i];
          ttree->Fill();
          time++;
        };
      };

      if (verbosityLevel >= 5)
        make_debug_histogram(event, count, waveform, stmEnergyCalib, deconvolved_data, differentiated_data, averaged_data, baseline_mean, baseline_stddev, peak_heights, peak_times);
      ++count;
    };

    if (verbosityLevel)
      std::cout << "MWD: " << channel.name() << ": " << outputMWDDigis->size() << " MWD digis found" << std::endl;
    event.put(std::move(outputMWDDigis));
  };

  void STMMovingWindowDeconvolution::deconvolve() {
    if (verbosityLevel > 2) {
      std::cout << "MWD: input ADCs (" << ADCs.size() << "): ";
      for (int16_t data : ADCs)
        std::cout << data << ", ";
      std::cout << "\n" << std::endl;
    };
    deconvolved_data.push_back(ADCs[0] - pedestal);
    for(i = 1; i < nADCs; i++)
      deconvolved_data.push_back((ADCs[i] - pedestal) - timeFactor * (ADCs[i - 1] - pedestal) + deconvolved_data[i - 1]);
    if (verbosityLevel > 2) {
      std::cout << "MWD: deconvoluted data (" << deconvolved_data.size() << "): ";
      for (double data : deconvolved_data)
        std::cout << data << ", ";
      std::cout << "\n" << std::endl;
    };
  };

  void STMMovingWindowDeconvolution::differentiate() {
    for (i = 0; i < M; i++)
      differentiated_data.push_back(deconvolved_data[i]);
    for (i = M; i < nADCs; i++)
      differentiated_data.push_back(deconvolved_data[i] - deconvolved_data[i - M]);
    if (verbosityLevel > 2) {
      std::cout << "MWD: differentiated data (" << differentiated_data.size() << "): ";
      for (double data : differentiated_data)
        std::cout << data << ", ";
      std::cout << "\n" << std::endl;
    };
  };

  void STMMovingWindowDeconvolution::average() {
    // sum the first L-1 elements of differentiated data
    // and set the first L-1 elements of averaged data
    sum = 0.0;
    for (i = 0; i < L - 1; i++) {
      sum += differentiated_data[i];
      averaged_data.push_back(differentiated_data[i]); // TODO: should this be divided by sum/i?
    };
    sum += differentiated_data[L - 1];
    averaged_data.push_back(sum/L);
    for (i = L; i < nADCs; ++i) {
      sum += differentiated_data[i] - differentiated_data[i - L]; // move the sum across one sample
      averaged_data.push_back(sum/L);
    };
    if (verbosityLevel > 2) {
      std::cout << "MWD: averaged data (" << averaged_data.size() << "): ";
      for (double data : averaged_data)
        std::cout << data << ", ";
      std::cout << "\n" << std::endl;
    };
  };

  void STMMovingWindowDeconvolution::calculate_baseline() {
    i = M;
    foundBaselineData = false;
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::mean, tag::variance> > acc_data_without_peaks;
    // Remove peaks so that we can calculate the baseline of the averaged data
    while (i < nADCs){
      gradient = averaged_data[i + 1] - averaged_data[i];
      if(gradient < thresholdgrad) // if the gradient is too sharp (i.e. we have hit a peak)
        i += (M + 2 * L); // jump ahead a little bit
      else {
        acc_data_without_peaks(averaged_data[i]);
        i++;
        foundBaselineData = true;
      };
    };
    baseline_mean   = foundBaselineData ? extract_result<tag::mean>(acc_data_without_peaks) : defaultBaselineMean;
    baseline_stddev = foundBaselineData ? std::sqrt(extract_result<tag::variance>(acc_data_without_peaks)) : defaultBaselineSD;
  };

  void STMMovingWindowDeconvolution::find_peaks() {
    threshold_cut = baseline_mean - nsigma_cut * baseline_stddev;
    std::cout << "Threshold cut: " << threshold_cut << std::endl;
    std::cout << "M: " << M << std::endl;
    std::cout << "nADCs: " << nADCs << std::endl;
    lowest_height = 0;
    lowest_height_time = -1; // in clock ticks

    for(i = M; i < nADCs; i++){
      if (averaged_data[i] < threshold_cut) { // the waveforms are negative so if we go below this threshold we have seen a peak
        if (averaged_data[i] < averaged_data[i - 1] && averaged_data[i] < lowest_height){ // if the current value is lower than the previous value and lower than the lowest value we've seen so far
          lowest_height = averaged_data[i]; // record the lowest height
          if (lowest_height_time == -1)
            lowest_height_time = i; // record the time we cross the threshold
        }
        else
          continue;
      };
      if (lowest_height_time == -1) // this will be true if we haven't seen a peak yet
        continue;
      else if (averaged_data[i] > threshold_cut) { // if we have seen a peak and go above the cut
        peak_heights.push_back(lowest_height - baseline_mean);
        peak_times.push_back(lowest_height_time); // ct
        lowest_height_time = -1; // reset to 0 so we can find a new peak
        lowest_height = 0;
      };
    };
  };

  void STMMovingWindowDeconvolution::make_debug_histogram(const art::Event& event, int count, const STMWaveformDigi& waveform, const STMEnergyCalib& stmEnergyCalib, const std::vector<double>& deconvolved_data, const std::vector<double>& differentiated_data, const std::vector<double>& averaged_data, const double baseline_mean, const double baseline_stddev, const std::vector<double>& peak_heights, const std::vector<double>& peak_times) {
    art::ServiceHandle<art::TFileService> tfs;
    std::stringstream histsuffix;
    histsuffix.str("");
    histsuffix << "_evt" << event.event() << "_wvf" << count;

    pedestal = stmEnergyCalib.pedestal(channel);
    nsPerCt = stmEnergyCalib.nsPerCt(channel);
    Binning binning = STMUtils::getBinning(waveform, _xAxis, nsPerCt);
    TH1D* h_waveform = tfs->make<TH1D>(("h_waveform"+histsuffix.str()).c_str(), "Waveform", binning.nbins(),binning.low(),binning.high());
    TH1D* h_deconvolved = tfs->make<TH1D>(("h_deconvolved"+histsuffix.str()).c_str(), "Deconvolution", binning.nbins(),binning.low(),binning.high());
    TH1D* h_differentiated = tfs->make<TH1D>(("h_differentiated"+histsuffix.str()).c_str(), "Differentiated", binning.nbins(),binning.low(),binning.high());
    TH1D* h_averaged = tfs->make<TH1D>(("h_averaged"+histsuffix.str()).c_str(), "Averaged", binning.nbins(),binning.low(),binning.high());
    TH1D* h_baseline_mean = tfs->make<TH1D>(("h_baseline_mean"+histsuffix.str()).c_str(), "Baseline Mean", binning.nbins(),binning.low(),binning.high());
    TH1D* h_baseline_mean_plus_stddev = tfs->make<TH1D>(("h_baseline_mean_plus_stddev"+histsuffix.str()).c_str(), "Baseline Mean + StdDev", binning.nbins(),binning.low(),binning.high());
    TH1D* h_baseline_mean_minus_stddev = tfs->make<TH1D>(("h_baseline_mean_minus_stddev"+histsuffix.str()).c_str(), "Baseline Mean - StdDev", binning.nbins(),binning.low(),binning.high());
    TH1D* h_peak_threshold = tfs->make<TH1D>(("h_peak_threshold"+histsuffix.str()).c_str(), "Threshold", binning.nbins(),binning.low(),binning.high());

    for (size_t i = 0; i < deconvolved_data.size(); ++i) {
      h_waveform->SetBinContent(i+1, waveform.adcs()[i] - pedestal); // remove the pedestal
      h_deconvolved->SetBinContent(i+1, deconvolved_data[i]);
      h_differentiated->SetBinContent(i+1, differentiated_data[i]);
      h_averaged->SetBinContent(i+1, averaged_data[i]);
      h_baseline_mean->SetBinContent(i+1, baseline_mean);
      h_baseline_mean_plus_stddev->SetBinContent(i+1, baseline_mean + baseline_stddev);
      h_baseline_mean_minus_stddev->SetBinContent(i+1, baseline_mean - baseline_stddev);
      h_peak_threshold->SetBinContent(i+1, baseline_mean - nsigma_cut * baseline_stddev);
    }
    TH1D* h_peaks = tfs->make<TH1D>(("h_peaks"+histsuffix.str()).c_str(), "Peaks", binning.nbins(),binning.low(),binning.high());
    for (size_t i_peak = 0; i_peak < peak_heights.size(); ++i_peak) {
      //      std::cout << "t = " << peak_times[i_peak] << ", E = " << peak_heights[i_peak] << std::endl;
      h_peaks->SetBinContent(peak_times[i_peak]+1, peak_heights[i_peak]);
    }
  }
}

DEFINE_ART_MODULE(mu2e::STMMovingWindowDeconvolution)
