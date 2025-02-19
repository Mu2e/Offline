// Create zero-suppressed STMWaveformDigis from unsuppressed STMWaveformDigis
// Original authors: Claudia Alvarez-Garcia, Alex Keshavarzi, and Mark Lancaster (see DocDB-43057 for details)
// Adapted for Offline: Andy Edmonds
// Adapted: Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"


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
        fhicl::Atom<unsigned long int> window{ Name("window"), Comment("Calculate the gradient between ADC values this number of elements away from each other")};
        fhicl::Atom<unsigned long int> naverage{ Name("naverage"), Comment("Number of ADC values to average the gradient over")};
        fhicl::OptionalAtom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level")};
        fhicl::OptionalAtom<bool> makeTTreeGradients{ Name("makeTTreeGradients"), Comment("Controls whether to make the TTree with branches ADC, gradient, averagedGradient")};
        fhicl::OptionalAtom<bool> makeTTreeWaveforms{ Name("makeTTreeWaveforms"), Comment("Controls whether to make the TTree with branches time, ADC")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMZeroSuppression(const Parameters& conf);
    private:
      void beginJob() override;
      void produce(art::Event& e) override;
      void calculateGradient();
      void averageGradient();
      void findPeaks();
      void chooseStartsAndEnds(); // taking into account any overlapping data

      // fhicl variables
      art::ProductToken<STMWaveformDigiCollection> stmWaveformDigisToken;
      STMChannel stmChannel;
      double tBefore = 0.0;
      double tAfter = 0.0;
      double threshold = 0.0;
      unsigned long int window = 0;
      unsigned long int nAverage = 0;
      int verbosityLevel;
      bool makeTTreeGradients = false;
      bool makeTTreeWaveforms = false;

      // Averaging variables
      unsigned long int nAverageInitial = 0;  // number of ADC values to average the gradient over
      double av_gradient = 0;
      unsigned long int i = 0, j = 0;         // iterator variables

      // Averaging variables
      std::vector<int16_t> ADCs;        // ADC values
      unsigned long int nADCs = 0;      // number of ADC values in unsuppressed waveform
      std::vector<int16_t> zeros;       // For checking if any data is present
      std::vector<int16_t> gradients;   // ADC gradient (difference between consecutive ADC values)
      unsigned long int nGradients = 0; // number of ADC values in unsuppressed waveform minus the window length
      std::vector<double> avGradients;  // ADC averaged averaged gradient
      std::vector<double> avPeakTimes;  // times at averaged gradient points, as multiplier of clock ticks

      // Peak finding variables
      bool found_peak = false;                  // controls whether a peak has been found (i.e. if recording)
      unsigned int nADCBefore = 0;              // number of samples before peak in number of ADC values
      unsigned int nADCAfter = 0;               // number of samples after peak in number of ADC values
      std::vector<unsigned long int> peakTimes; // vector with each peak time in clock ticks
      unsigned long int nPeakTimes = 0;         // number of peak times
      unsigned int peakTime = 0;                // peak start time in clock ticks
      std::vector<size_t> peakStartTimes;       // start positions of each zero-suppressed waveform
      size_t peakStartTime = 0;                 // start position of each zero-suppressed waveform
      std::vector<size_t> peakEndTimes;         // end positions of each zero-suppressed waveform
      size_t peakEndTime = 0;                   // end position of each zero-suppressed waveform
      std::vector<size_t> finalPeakStartTimes;  // start positions of each zero-suppressed waveform after taking into account overlapping data
      std::vector<size_t> finalPeakEndTimes;    // end positions of each zero-suppressed waveform after taking into account overlapping data

      // New data product generation variables
      size_t nZSwaveforms = 0;                  // number of peaks to include in output product
      std::vector<int16_t> ZSADCs;              // zero suppressed waveform
      size_t k = 0;                             // iterator

      // Proditions service
      ProditionsHandle<STMEnergyCalib> stmEnergyCalibHandle;

      // TTree variables
      TTree* ttree;
      int16_t ADC = 0, gradient = 0, averagedGradient = 0;
      uint32_t time = 0, timeStep = 0;
  };

  STMZeroSuppression::STMZeroSuppression(const Parameters& conf):
    art::EDProducer{conf},
    stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(conf().stmWaveformDigisTag())),
    stmChannel(STMUtils::getChannel(conf().stmWaveformDigisTag())),
    tBefore(conf().tbefore()),
    tAfter(conf().tafter()),
    threshold(conf().threshold()),
    window(conf().window()),
    nAverage(conf().naverage()) {
      produces<STMWaveformDigiCollection>();
      nAverageInitial = nAverage;
      verbosityLevel = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;
      if (verbosityLevel > 10)
        verbosityLevel = 10;
      makeTTreeGradients = conf().makeTTreeGradients() ? *(conf().makeTTreeGradients()) : false;
      makeTTreeWaveforms = conf().makeTTreeWaveforms() ? *(conf().makeTTreeWaveforms()) : false;
      if (makeTTreeGradients) {
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "STMZeroSuppression gradients ttree");
        ttree->Branch("ADC", &ADC, "ADC/S");
        ttree->Branch("gradient", &gradient, "gradient/S");
        ttree->Branch("averagedGradient", &averagedGradient, "averagedGradient/S");
      };
      if (makeTTreeWaveforms) {
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "STMZeroSuppression waveforms ttree");
        ttree->Branch("time", &time, "time/I");
        ttree->Branch("ADC", &ADC, "ADC/S");
      };
    };

  void STMZeroSuppression::beginJob() {
    if (verbosityLevel) {
      std::cout << "STM Zero-Suppression Algorithm Parameters:" << std::endl;
      std::cout << std::left << "\t" << std::setw(15) << "tbefore"   << tBefore   << " ns" << std::endl;
      std::cout << std::left << "\t" << std::setw(15) << "tafter"    << tAfter    << " ns" << std::endl;
      std::cout << std::left << "\t" << std::setw(15) << "threshold" << threshold << std::endl;
      std::cout << std::endl; // buffer line
      std::cout << "STM Channel" << std::endl;
      std::cout << std::left << "\t" << std::setw(15) << "Name" << stmChannel.name()                  << std::endl;
      std::cout << std::left << "\t" << std::setw(15) << "ID"   << static_cast<int>(stmChannel.id())  << std::endl;
      std::cout << std::endl; // buffer line
    };
    return;
  };

  void STMZeroSuppression::produce(art::Event& event) {
    // create output
    auto waveformsHandle = event.getValidHandle(stmWaveformDigisToken);
    std::unique_ptr<STMWaveformDigiCollection> outputSTMWaveformDigis(new STMWaveformDigiCollection());
    // Get the prodition
    const STMEnergyCalib& stmEnergyCalib = stmEnergyCalibHandle.get(event.id());
    nADCBefore = STMUtils::convertToClockTicks(tBefore, stmChannel, stmEnergyCalib);
    nADCAfter = STMUtils::convertToClockTicks(tAfter, stmChannel, stmEnergyCalib);
    timeStep = tBefore/nADCBefore;
    if (verbosityLevel > 9) {
      std::cout << "ZS findPeaks fitting parameters" << std::endl;
      std::cout << std::left << std::setw(15) << "nADCBefore" << nADCBefore << std::endl;
      std::cout << std::left << std::setw(15) << "nADCAfter"  << nADCAfter  << std::endl;
      std::cout << std::endl; // buffer line
    };

    for (const auto& waveform : *waveformsHandle) {
      // Get the data product
      ADCs = waveform.adcs();
      // Assign the correct amount of space to all the data vectors
      nADCs = ADCs.size();
      zeros.clear();
      zeros.assign(nADCs, 0);
      if (ADCs == zeros) // Needed for simulation work
        continue;
      nGradients = nADCs - window;
      gradients.clear();
      avGradients.clear();
      avPeakTimes.clear();
      peakTimes.clear();

      // Calculate the ADC gradients
      calculateGradient();
      // Average the ADC gradients
      averageGradient();
      // Find the peaks
      findPeaks();
      // Remove overlapping peaks
      chooseStartsAndEnds();
      // Generate the output waveforms
      nZSwaveforms = finalPeakStartTimes.size();
      for (k = 0; k < nZSwaveforms; ++k) {
        peakStartTime = finalPeakStartTimes[k];
        peakEndTime = finalPeakEndTimes[k];
        ZSADCs.clear();
        ZSADCs.assign(waveform.adcs().begin() + peakStartTime, waveform.adcs().begin() + peakEndTime);
        if (verbosityLevel > 5) {
          std::cout << "ZSADCs: ";
          for (auto i : ZSADCs)
            std::cout << i << ", ";
          std::cout << "\n" << std::endl;
        };
        STMWaveformDigi stm_waveform(waveform.trigTimeOffset() + peakStartTime, ZSADCs);
        outputSTMWaveformDigis->push_back(stm_waveform);
      };
      // Save data to TTree
      if (makeTTreeGradients) {
        for (i = 0; i < nADCs; i++) {
          ADC = ADCs[i];
          j = i > nGradients ? nGradients : i;
          gradient = gradients[j];
          averagedGradient = avGradients[j];
          ttree->Fill();
        };
      };
      if (makeTTreeWaveforms) {
        time = (uint32_t) waveform.trigTimeOffset() + peakStartTime - timeStep;
        for (i = 0; i < nADCs; i++) {
          ADC = ADCs[i];
          time += timeStep;
          ttree->Fill();
        };
      };
    };
    if (verbosityLevel)
      std::cout << stmChannel.name() << ": " << outputSTMWaveformDigis->size() << " waveforms found" << std::endl;
    event.put(std::move(outputSTMWaveformDigis));
    return;
  };

  void STMZeroSuppression::calculateGradient() {
    // Print
    if (verbosityLevel > 9) {
      std::cout << "ZS: ADCs (" << ADCs.size() << " entries): ";
      for (int16_t ADC : ADCs)
        std::cout << ADC << ", ";
      std::cout << "\n" << std::endl;
    };
    // Calculate the time and gradient vectors
    for(i = 0; i < nGradients; i++)
      gradients.push_back(ADCs[i + window] - ADCs[i]);
    // Print gradients
    if (verbosityLevel > 9) {
      std::cout << "ZS: Gradient (" << gradients.size() << " entries): ";
      for (int16_t gradient : gradients)
        std::cout << gradient << ", ";
      std::cout << "\n" << std::endl;
    };
    return;
  };

  void STMZeroSuppression::averageGradient() {
    // Calculate the average of the gradient each nAverage ADC values to avoid fluctuations
    nAverage = nAverageInitial;
    for (i = 0; i < nGradients; i++) {
      // Each point of the gradient and ADCtime averaged with the (nAverage - 1) following points, each point is the mean of nAverage points of the gradient
      if((i + nAverage) > nGradients)
        nAverage = nGradients - i; // Use the rest of the vector
      // Determine the sum
      for (j = 0; j < nAverage; j++)
        av_gradient = av_gradient + gradients[i + j];
      // Save the results
      avGradients.push_back(av_gradient / nAverage);
      avPeakTimes.push_back(i);
      // Reset for next loop
      av_gradient = 0;
    };
    // Print
    if (verbosityLevel > 9) {
      std::cout << "ZS: Avg Gradient (" << avGradients.size() << " entries): ";
      for (int16_t grad : avGradients)
        std::cout << grad << ", ";
      std::cout << "\n" << std::endl;
    };
    return;
  };

  void STMZeroSuppression::findPeaks() {
    // STMEnergyCalib required to get the sampling frequency
    found_peak = false;
    peakStartTimes.clear();
    peakEndTimes.clear();
    // Store positions in clock ticks for the peaks found
    for (i = 0; i < nGradients; ++i) {
      // Check if current ADC is above threshold
      if (avGradients[i] > threshold) {
        found_peak=false;
        continue;
      };
      // Skip the rest indexes of the peak after the peak that has already been stored
      if (avGradients[i] < threshold && found_peak)
        continue;
      // If there is no data above threshold yet
      if (avGradients[i] < threshold && !found_peak) {
        found_peak = true;
        peakTime = avPeakTimes[i];
        peakTimes.push_back(peakTime);
        // Assign the appropriate peak start and end time
        peakStartTimes.push_back(peakTime < nADCBefore ? 0 : peakTime - nADCBefore);
        peakEndTimes.push_back(peakTime > (nADCs - nADCAfter) ? nADCs : peakTime + nADCAfter);
      };
    };
    return;
  };

  void STMZeroSuppression::chooseStartsAndEnds() {
    // Setup the vectors
    finalPeakStartTimes.clear();
    finalPeakEndTimes.clear();
    nPeakTimes = peakTimes.size();
    // Now go through and account for overlapped data
    if (nPeakTimes == 0) // need to be careful just in case there were no peaks found (i.e. just noise)
      return;
    // Collect the peaks
    for (i = 0; i < nPeakTimes; ++i) {
      peakStartTime = peakStartTimes[i];
      peakEndTime = peakEndTimes[i];
      if (i == (nPeakTimes - 1)) { // the final peak
        finalPeakStartTimes.push_back(peakStartTime);
        finalPeakEndTimes.push_back(peakEndTime);
        break;
      }
      else if (peakStartTimes[i + 1] <= peakEndTime) { // peaks are overlapping
        peakEndTime = peakEndTimes.at(i + 1); // update so that we will go to the end of the second peak
      }
      else if (peakStartTimes[i + 1] > peakEndTime) { // peaks don't overlap
        finalPeakStartTimes.push_back(peakStartTime);
        finalPeakEndTimes.push_back(peakEndTime);
      };
    };
    return;
  };
};

DEFINE_ART_MODULE(mu2e::STMZeroSuppression)
