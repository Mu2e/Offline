// ======================================================================
//
// StrawDigiCollectionsComparator_module: Compare two StrawDigi collections
// (and optionally their associated StrawDigiADCWaveform collections)
//
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace mu2e {

class StrawDigiCollectionsComparator : public art::EDAnalyzer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> referenceDigiTag{fhicl::Name("referenceDigiTag"),
                                                fhicl::Comment("Reference StrawDigiCollection")};
    fhicl::Atom<art::InputTag> candidateDigiTag{fhicl::Name("candidateDigiTag"),
                                                fhicl::Comment("Candidate StrawDigiCollection")};
    fhicl::Atom<bool> compareWaveforms{fhicl::Name("compareWaveforms"),
                                       fhicl::Comment("Compare StrawDigiADCWaveformCollection"),
                                       true};
    fhicl::Atom<art::InputTag> referenceWaveformTag{
      fhicl::Name("referenceWaveformTag"),
      fhicl::Comment("Reference StrawDigiADCWaveformCollection"),
      art::InputTag("makeSD")};
    fhicl::Atom<art::InputTag> candidateWaveformTag{
      fhicl::Name("candidateWaveformTag"),
      fhicl::Comment("Candidate StrawDigiADCWaveformCollection"),
      art::InputTag("makeSD")};
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
    fhicl::Atom<bool> normalizeWaveformsTo10Bit{
      fhicl::Name("normalizeWaveformsTo10Bit"),
      fhicl::Comment("Apply 10-bit normalization (sample & 0x3FF) to both reference and candidate waveforms"),
      false};
    fhicl::Atom<bool> normalizePmpTo10Bit{
      fhicl::Name("normalizePmpTo10Bit"),
      fhicl::Comment("Apply 10-bit normalization (PMP & 0x3FF) to both reference and candidate digis"),
      false};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit StrawDigiCollectionsComparator(Parameters const& config);

  void analyze(art::Event const& event) override;
  void endJob() override;

private:
  struct FlatDigi {
    uint16_t strawID;
    uint32_t tdc0;
    uint32_t tdc1;
    uint32_t tot0;
    uint32_t tot1;
    uint32_t pmp;
    std::vector<int> waveform;

    bool operator<(FlatDigi const& other) const {
      return std::tie(strawID, tdc0, tdc1, tot0, tot1, pmp, waveform) <
             std::tie(other.strawID, other.tdc0, other.tdc1, other.tot0, other.tot1, other.pmp,
                      other.waveform);
    }

    bool operator==(FlatDigi const& other) const {
      return strawID == other.strawID && tdc0 == other.tdc0 && tdc1 == other.tdc1 &&
             tot0 == other.tot0 && tot1 == other.tot1 && pmp == other.pmp &&
             waveform == other.waveform;
    }
  };

  static std::vector<FlatDigi> normalize(mu2e::StrawDigiCollection const& digis,
                                         mu2e::StrawDigiADCWaveformCollection const* waveforms,
                                         bool compareWaveforms,
                                         bool normalizeWaveformsTo10Bit,
                                         bool normalizePmpTo10Bit);
  static std::string describe(FlatDigi const& digi, size_t maxWaveformSamplesToPrint);

  art::InputTag referenceDigiTag_;
  art::InputTag candidateDigiTag_;
  bool compareWaveforms_;
  art::InputTag referenceWaveformTag_;
  art::InputTag candidateWaveformTag_;
  bool failOnMismatch_;
  int diagLevel_;
  size_t maxMismatchesToPrint_;
  size_t maxWaveformSamplesToPrint_;
  bool normalizeWaveformsTo10Bit_;
  bool normalizePmpTo10Bit_;

  size_t totalEvents_{0};
  size_t matchingEvents_{0};
  size_t mismatchedEvents_{0};
};

StrawDigiCollectionsComparator::StrawDigiCollectionsComparator(Parameters const& config)
    : art::EDAnalyzer{config}
    , referenceDigiTag_(config().referenceDigiTag())
    , candidateDigiTag_(config().candidateDigiTag())
    , compareWaveforms_(config().compareWaveforms())
    , referenceWaveformTag_(config().referenceWaveformTag())
    , candidateWaveformTag_(config().candidateWaveformTag())
    , failOnMismatch_(config().failOnMismatch())
    , diagLevel_(config().diagLevel())
    , maxMismatchesToPrint_(config().maxMismatchesToPrint())
    , maxWaveformSamplesToPrint_(config().maxWaveformSamplesToPrint())
    , normalizeWaveformsTo10Bit_(config().normalizeWaveformsTo10Bit())
    , normalizePmpTo10Bit_(config().normalizePmpTo10Bit()) {}

std::vector<StrawDigiCollectionsComparator::FlatDigi> StrawDigiCollectionsComparator::normalize(
    mu2e::StrawDigiCollection const& digis, mu2e::StrawDigiADCWaveformCollection const* waveforms,
  bool compareWaveforms, bool normalizeWaveformsTo10Bit, bool normalizePmpTo10Bit) {
  if (compareWaveforms && waveforms != nullptr && waveforms->size() != digis.size()) {
    throw cet::exception("STRAWDIGI_COMPARE")
        << "StrawDigi and StrawDigiADCWaveform collection sizes differ: digis=" << digis.size()
        << ", waveforms=" << waveforms->size();
  }

  std::vector<FlatDigi> out;
  out.reserve(digis.size());

  for (size_t i = 0; i < digis.size(); ++i) {
    auto const& digi = digis[i];
    std::vector<int> waveform;

    if (compareWaveforms && waveforms != nullptr) {
      auto const& samples = waveforms->at(i).samples();
      waveform.reserve(samples.size());
      if (normalizeWaveformsTo10Bit) {
        for (auto const sample : samples) {
          waveform.push_back(static_cast<int>(sample & 0x3FF));
        }
      } else {
        waveform.assign(samples.begin(), samples.end());
      }
    }

    auto const normalizedPmp = normalizePmpTo10Bit
                   ? static_cast<uint32_t>(digi.PMP() & 0x3FF)
                   : static_cast<uint32_t>(digi.PMP());

    out.push_back(FlatDigi{static_cast<uint16_t>(digi.strawId().asUint16()),
                 static_cast<uint32_t>(digi.TDC(mu2e::StrawEnd::cal)),
                 static_cast<uint32_t>(digi.TDC(mu2e::StrawEnd::hv)),
                 static_cast<uint32_t>(digi.TOT(mu2e::StrawEnd::cal)),
                 static_cast<uint32_t>(digi.TOT(mu2e::StrawEnd::hv)),
                 normalizedPmp,
                 std::move(waveform)});
  }

  // Sort so comparison is robust against ordering differences between producers.
  std::sort(out.begin(), out.end());
  return out;
}

std::string StrawDigiCollectionsComparator::describe(FlatDigi const& digi,
                                                     size_t maxWaveformSamplesToPrint) {
  std::ostringstream os;
  os << "strawID=" << digi.strawID << ", tdc0=" << digi.tdc0 << ", tdc1=" << digi.tdc1
     << ", tot0=" << digi.tot0 << ", tot1=" << digi.tot1 << ", pmp=" << digi.pmp
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

void StrawDigiCollectionsComparator::analyze(art::Event const& event) {
  ++totalEvents_;

  auto const& referenceDigis = event.getValidHandle<mu2e::StrawDigiCollection>(referenceDigiTag_);
  auto const& candidateDigis = event.getValidHandle<mu2e::StrawDigiCollection>(candidateDigiTag_);

  mu2e::StrawDigiADCWaveformCollection const* referenceWaveforms = nullptr;
  mu2e::StrawDigiADCWaveformCollection const* candidateWaveforms = nullptr;
  if (compareWaveforms_) {
    referenceWaveforms =
        event.getValidHandle<mu2e::StrawDigiADCWaveformCollection>(referenceWaveformTag_).product();
    candidateWaveforms =
        event.getValidHandle<mu2e::StrawDigiADCWaveformCollection>(candidateWaveformTag_).product();
  }

  auto reference = normalize(*referenceDigis,
                             referenceWaveforms,
                             compareWaveforms_,
                             normalizeWaveformsTo10Bit_,
                             normalizePmpTo10Bit_);
  auto candidate = normalize(*candidateDigis,
                             candidateWaveforms,
                             compareWaveforms_,
                             normalizeWaveformsTo10Bit_,
                             normalizePmpTo10Bit_);

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
      std::cout << "[StrawDigiCollectionsComparator] Event " << event.id()
                << " collections are equivalent (" << reference.size() << " digis)" << std::endl;
    }
    return;
  }

  ++mismatchedEvents_;

  std::ostringstream summary;
  summary << "[StrawDigiCollectionsComparator] Event " << event.id() << " mismatch between "
          << referenceDigiTag_ << " and " << candidateDigiTag_ << "\n";
  for (auto const& mismatch : mismatches) {
    summary << mismatch << "\n";
  }

  if (diagLevel_ > 0 || failOnMismatch_) {
    // Keep human-readable details even when failOnMismatch is false.
    std::cout << summary.str();
  }

  if (failOnMismatch_) {
    throw cet::exception("STRAWDIGI_COMPARE") << summary.str();
  }
}

void StrawDigiCollectionsComparator::endJob() {
  if (diagLevel_ > 0) {
    std::cout << "\n ----- [StrawDigiCollectionsComparator] Summary ----- " << std::endl;
    std::cout << "Total events: " << totalEvents_ << std::endl;
    std::cout << "Matching events: " << matchingEvents_ << std::endl;
    std::cout << "Mismatched events: " << mismatchedEvents_ << std::endl;
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StrawDigiCollectionsComparator)
