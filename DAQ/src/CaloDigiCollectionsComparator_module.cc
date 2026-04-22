// ======================================================================
//
// CaloDigiCollectionsComparator_module: Compare two CaloDigiCollections
//
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace mu2e {

class CaloDigiCollectionsComparator : public art::EDAnalyzer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> referenceTag{fhicl::Name("referenceTag"),
                                            fhicl::Comment("Reference CaloDigiCollection")};
    fhicl::Atom<art::InputTag> candidateTag{fhicl::Name("candidateTag"),
                                            fhicl::Comment("Candidate CaloDigiCollection")};
    fhicl::Atom<bool> failOnMismatch{fhicl::Name("failOnMismatch"),
                                     fhicl::Comment("Throw if any mismatch is found"),
                                     true};
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"),
                               fhicl::Comment("Diagnostic verbosity"),
                               0};
    fhicl::Atom<size_t> maxMismatchesToPrint{fhicl::Name("maxMismatchesToPrint"),
                                             fhicl::Comment("Maximum mismatch entries to print"),
                                             10};
    fhicl::Atom<size_t> maxWaveformSamplesToPrint{
      fhicl::Name("maxWaveformSamplesToPrint"),
      fhicl::Comment("Maximum waveform samples to print in mismatch diagnostics"),
      64};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CaloDigiCollectionsComparator(Parameters const& config);

  void analyze(art::Event const& event) override;
  void endJob() override;

private:
  struct FlatDigi {
    int sipmID;
    int t0;
    int peakpos;
    std::vector<int> waveform;

    bool operator<(FlatDigi const& other) const {
      return std::tie(sipmID, t0, peakpos, waveform) <
             std::tie(other.sipmID, other.t0, other.peakpos, other.waveform);
    }

    bool operator==(FlatDigi const& other) const {
      return sipmID == other.sipmID && t0 == other.t0 && peakpos == other.peakpos &&
             waveform == other.waveform;
    }
  };

  static std::vector<FlatDigi> normalize(mu2e::CaloDigiCollection const& digis);
  static std::string describe(FlatDigi const& digi, size_t maxWaveformSamplesToPrint);

  art::InputTag referenceTag_;
  art::InputTag candidateTag_;
  bool failOnMismatch_;
  int diagLevel_;
  size_t maxMismatchesToPrint_;
  size_t maxWaveformSamplesToPrint_;

  size_t totalEvents_{0};
  size_t matchingEvents_{0};
  size_t mismatchedEvents_{0};
};

CaloDigiCollectionsComparator::CaloDigiCollectionsComparator(Parameters const& config)
    : art::EDAnalyzer{config}
    , referenceTag_(config().referenceTag())
    , candidateTag_(config().candidateTag())
    , failOnMismatch_(config().failOnMismatch())
    , diagLevel_(config().diagLevel())
    , maxMismatchesToPrint_(config().maxMismatchesToPrint())
    , maxWaveformSamplesToPrint_(config().maxWaveformSamplesToPrint()) {}

std::vector<CaloDigiCollectionsComparator::FlatDigi>
CaloDigiCollectionsComparator::normalize(mu2e::CaloDigiCollection const& digis) {
  std::vector<FlatDigi> out;
  out.reserve(digis.size());
  for (auto const& digi : digis) {
    out.push_back(FlatDigi{digi.SiPMID(), digi.t0(), digi.peakpos(), digi.waveform()});
  }
  // Sort so comparison is robust against ordering differences between producers.
  std::sort(out.begin(), out.end());
  return out;
}

std::string CaloDigiCollectionsComparator::describe(FlatDigi const& digi,
                                                    size_t maxWaveformSamplesToPrint) {
  std::ostringstream os;
  os << "SiPMID=" << digi.sipmID << ", t0=" << digi.t0 << ", peakpos=" << digi.peakpos
     << ", waveformSize=" << digi.waveform.size();
  if (!digi.waveform.empty()) {
    os << ", waveform={";
    size_t const n = digi.waveform.size();
    if (n <= maxWaveformSamplesToPrint) {
      for (size_t i = 0; i < n; ++i) {
        if (i != 0) {
          os << ",";
        }
        os << digi.waveform[i];
      }
    } else {
      size_t const head = maxWaveformSamplesToPrint / 2;
      size_t const tail = maxWaveformSamplesToPrint - head;
      for (size_t i = 0; i < head; ++i) {
        if (i != 0) {
          os << ",";
        }
        os << digi.waveform[i];
      }
      os << ",...,";
      for (size_t i = n - tail; i < n; ++i) {
        if (i != n - tail) {
          os << ",";
        }
        os << digi.waveform[i];
      }
    }
    os << "}";
  }
  return os.str();
}

void CaloDigiCollectionsComparator::analyze(art::Event const& event) {
  ++totalEvents_;

  auto const& referenceHandle = event.getValidHandle<mu2e::CaloDigiCollection>(referenceTag_);
  auto const& candidateHandle = event.getValidHandle<mu2e::CaloDigiCollection>(candidateTag_);

  auto reference = normalize(*referenceHandle);
  auto candidate = normalize(*candidateHandle);

  std::vector<std::string> mismatches;
  if (reference.size() != candidate.size()) {
    std::ostringstream os;
    os << "Collection sizes differ: reference=" << reference.size()
       << ", candidate=" << candidate.size();
    mismatches.push_back(os.str());
  }

  size_t compareCount = std::min(reference.size(), candidate.size());
  for (size_t i = 0; i < compareCount && mismatches.size() < maxMismatchesToPrint_; ++i) {
    if (!(reference[i] == candidate[i])) {
      std::ostringstream os;
      os << "Digi mismatch at sorted index " << i << "\n"
        << "  reference: " << describe(reference[i], maxWaveformSamplesToPrint_) << "\n"
        << "  candidate: " << describe(candidate[i], maxWaveformSamplesToPrint_);
      mismatches.push_back(os.str());
    }
  }

  if (mismatches.empty()) {
    ++matchingEvents_;
    if (diagLevel_ > 1) {
      std::cout << "[CaloDigiCollectionsComparator] Event " << event.id()
                << " collections are equivalent (" << reference.size() << " digis)" << std::endl;
    }
    return;
  }

  ++mismatchedEvents_;

  std::ostringstream summary;
  summary << "[CaloDigiCollectionsComparator] Event " << event.id() << " mismatch between "
          << referenceTag_ << " and " << candidateTag_ << "\n";
  for (auto const& mismatch : mismatches) {
    summary << mismatch << "\n";
  }

  if (diagLevel_ > 0 || failOnMismatch_) {
    // Keep human-readable details even when failOnMismatch is false.
    std::cout << summary.str();
  }

  if (failOnMismatch_) {
    throw cet::exception("CALODIGI_COMPARE") << summary.str();
  }
}

void CaloDigiCollectionsComparator::endJob() {
  if (diagLevel_ > 0) {
    std::cout << "\n ----- [CaloDigiCollectionsComparator] Summary ----- " << std::endl;
    std::cout << "Total events: " << totalEvents_ << std::endl;
    std::cout << "Matching events: " << matchingEvents_ << std::endl;
    std::cout << "Mismatched events: " << mismatchedEvents_ << std::endl;
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloDigiCollectionsComparator)
