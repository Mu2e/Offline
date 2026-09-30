// ====================================================================
//
// STMDigisFromFragments: create all types of STMDigis from STMFragments
//
// ======================================================================

// framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// offline includes
#include "Offline/RecoDataProducts/inc/STMFragmentSummary.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "Offline/RecoDataProducts/inc/STMEventHeader.hh"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

// c++
#include <string>
#include <sstream>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

namespace art
{
    class STMDigisFromFragments;
}
using art::STMDigisFromFragments;
class art::STMDigisFromFragments : public EDProducer
{
public:
    struct Config {
        //General Configs
        fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
        fhicl::OptionalAtom<int> verbosityLevel {fhicl::Name("verbosityLevel"),
            fhicl::Comment("Verbosity Level for debugging purposes")};
        fhicl::Atom<bool> saveSTMFragSummary{fhicl::Name("saveSTMFragSummary"),
            fhicl::Comment("Whether to save Fragment Summary"), false};

        //HPGe Configs
        fhicl::Atom<bool> saveRawWaveformsWithHeaderHPGe {fhicl::Name("saveRawWaveformsWithHeaderHPGe"),
            fhicl::Comment("Whether to save raw fragment waveforms with header for HPGe for debugging purposes"), false};
        fhicl::Atom<bool> saveRawWaveformsHPGe {fhicl::Name("saveRawWaveformsHPGe"),
            fhicl::Comment("Whether to save raw waveforms for HPGe for debugging purposes"), false};
        fhicl::Atom<bool> saveZSWaveformsHPGe {fhicl::Name("saveZSWaveformsHPGe"),
            fhicl::Comment("Whether to save zero-suppressed waveforms for HPGe for debugging purposes"), false};

        //LaBr Configs
        fhicl::Atom<bool> saveRawWaveformsWithHeaderLaBr {fhicl::Name("saveRawWaveformsWithHeaderLaBr"),
            fhicl::Comment("Whether to save raw fragment waveforms with header for LaBr for debugging purposes"), false};
        fhicl::Atom<bool> saveRawWaveformsLaBr {fhicl::Name("saveRawWaveformsLaBr"),
            fhicl::Comment("Whether to save raw waveforms for LaBr for debugging purposes"), false};
        fhicl::Atom<bool> saveZSWaveformsLaBr {fhicl::Name("saveZSWaveformsLaBr"),
            fhicl::Comment("Whether to save zero-suppressed waveforms for LaBr for debugging purposes"), false};
    };

    explicit STMDigisFromFragments (const art::EDProducer::Table<Config>& config);
    virtual void produce(Event &) override;
    void endJob() override;

private:
    // art input tags
    art::InputTag _stmFragmentsTag;

    // Fhicl parameters
    int _verbosityLevel{0};
    bool _saveSTMFragSummary{false};
    // HPGe fhicl parameters
    bool _saveRawWaveformsWithHeaderHPGe{false};
    bool _saveRawWaveformsHPGe{false};
    bool _saveZSWaveformsHPGe{false};
    // LaBr fhicl parameters
    bool _saveRawWaveformsWithHeaderLaBr{false};
    bool _saveRawWaveformsLaBr{false};
    bool _saveZSWaveformsLaBr{false};

    // General Fragment variables
    size_t _totalEvents{0};
    size_t _totalNonContainers{0};
    size_t _totalUnknownContainers{0};
    size_t _totalFragments{0};
    size_t _totalContainers{0};
    size_t _totalInnerFrags{0};
    size_t _totalUnreadInnerFrags{0};
    size_t _totalRawFragsSeen{0};
    size_t _totalZSFragsSeen{0};
    size_t _totalPHFragsSeen{0};
    size_t _totalGoodRawFrags{0};
    size_t _totalGoodZSFrags{0};
    size_t _totalGoodPHFrags{0};
    size_t _totalZeroRawFrags{0};
    size_t _totalZeroZSFrags{0};
    size_t _totalZeroPHFrags{0};
    size_t _totalEmptyRawFrags{0};
    size_t _totalEmptyZSFrags{0};
    size_t _totalEmptyPHFrags{0};

    // zs errors
    size_t _totalZSLengthMismatch{0};
    size_t _totalZSLengthMismatchHPGe{0};
    size_t _totalZSLengthMismatchLaBr{0};
    size_t _totalZSRegionMismatch{0};
    size_t _totalZSRegionMismatchHPGe{0};
    size_t _totalZSRegionMismatchLaBr{0};

    // Header-related variables
    size_t _totalRawFragsPrescaled{0};
    size_t _totalZSFragsPrescaled{0};
    size_t _totalRawFragsFlaggedBadOnly{0};
    size_t _totalRawFragsFlaggedMissingOnly{0};
    size_t _totalRawFragsFlaggedBadAndMissing{0};
    size_t _totalRawFragsWithInvalidHeaders{0};
    size_t _totalRawFragsWithInvalidAnchors{0};
    size_t _totalPHCountMismatch{0};

    // Detector Summary variables
    size_t _totalEventsWithBothDetectors{0};
    size_t _totalEventsWithOnlyHPGe{0};
    size_t _totalEventsWithOnlyLaBr{0};
    size_t _totalEventsWithNeitherDetector{0}; // Should sum to total art events

    // HPGe job-level variables
    size_t _totalContainersHPGe{0};
    size_t _totalInnerFragsHPGe{0};
    size_t _totalRawFragsSeenHPGe{0};
    size_t _totalZSFragsSeenHPGe{0};
    size_t _totalPHFragsSeenHPGe{0};
    size_t _totalGoodRawFragsHPGe{0};//Good means fragment has data
    size_t _totalGoodZSFragsHPGe{0};
    size_t _totalGoodPHFragsHPGe{0};
    size_t _totalZeroRawFragsHPGe{0};// Zero means fragment has data but all values are zero
    size_t _totalZeroZSFragsHPGe{0};
    size_t _totalZeroPHFragsHPGe{0};
    size_t _totalEmptyRawFragsHPGe{0};// Empty means fragment has no data
    size_t _totalEmptyZSFragsHPGe{0};
    size_t _totalEmptyPHFragsHPGe{0};
    size_t _totalRawFragsFlaggedBadOnlyHPGe{0};
    size_t _totalRawFragsFlaggedMissingOnlyHPGe{0};
    size_t _totalRawFragsFlaggedBadAndMissingHPGe{0};
    size_t _totalRawFragsWithInvalidHeadersHPGe{0};
    size_t _totalRawFragsWithInvalidAnchorsHPGe{0};
    size_t _totalZSFragsSkippedDueToRawFlagHPGe{0};// We want to skip ZS if raw is flagged bad/missing
    size_t _totalPHFragsSkippedDueToRawFlagHPGe{0};// We want to skip PH if raw is flagged bad/missing
    size_t _totalRawFragsPrescaledHPGe{0};// track how many raw fragments were prescaled
    size_t _totalZSFragsPrescaledHPGe{0};// track how many zs fragments were prescaled
    size_t _totalPHCountMismatchHPGe{0};

    // LaBr job-level variables
    size_t _totalContainersLaBr{0};
    size_t _totalInnerFragsLaBr{0};
    size_t _totalRawFragsSeenLaBr{0};
    size_t _totalZSFragsSeenLaBr{0};
    size_t _totalPHFragsSeenLaBr{0};
    size_t _totalGoodRawFragsLaBr{0};//Good means fragment has data
    size_t _totalGoodZSFragsLaBr{0};
    size_t _totalGoodPHFragsLaBr{0};
    size_t _totalZeroRawFragsLaBr{0};// Zero means fragment has data but all values are zero
    size_t _totalZeroZSFragsLaBr{0};
    size_t _totalZeroPHFragsLaBr{0};
    size_t _totalEmptyRawFragsLaBr{0};// Empty means fragment has no data
    size_t _totalEmptyZSFragsLaBr{0};
    size_t _totalEmptyPHFragsLaBr{0};
    size_t _totalRawFragsFlaggedBadOnlyLaBr{0};
    size_t _totalRawFragsFlaggedMissingOnlyLaBr{0};
    size_t _totalRawFragsFlaggedBadAndMissingLaBr{0};
    size_t _totalRawFragsWithInvalidHeadersLaBr{0};
    size_t _totalRawFragsWithInvalidAnchorsLaBr{0};
    size_t _totalZSFragsSkippedDueToRawFlagLaBr{0};// We want to skip ZS if raw is flagged bad/missing
    size_t _totalPHFragsSkippedDueToRawFlagLaBr{0};// We want to skip PH if raw is flagged bad/missing
    size_t _totalRawFragsPrescaledLaBr{0};// track how many raw fragments were prescaled
    size_t _totalZSFragsPrescaledLaBr{0};// track how many ZS fragments were prescaled
    size_t _totalPHCountMismatchLaBr{0};
    size_t _totalNoHitsPHFrags{0};
    size_t _totalNoHitsPHFragsHPGe{0};
    size_t _totalNoHitsPHFragsLaBr{0};

    // Used to save ZS information
    struct ZSRegion{
        uint32_t offset;
        std::vector<int16_t> adcs;
    };

    // Used to track the expected raw header information, event based
    struct RawHeaderState {
        uint16_t expectedZSLength{0};
        uint16_t expectedZSRegions{0};
        uint16_t expectedPHCount{0};

        bool containsZSInfo{false};
        bool containsPHInfo{false};

        bool rawPrescaled{false};
        bool zsPrescaled{false};

        uint16_t rawPrescaleValue{0};
        uint16_t zsPrescaleValue{0};

        bool skipCurrentSetDueToRawFlags{false};
        bool skipCurrentSetDueToInvalidHeader{false};
        // save the following now for mapping
        uint64_t eventWindowTag{0};
        uint8_t eventMode{0};
        uint64_t adcClock{0};
        uint64_t dtcClock{0};
    };

    struct FragmentCounters {
        size_t seen{0};
        size_t prescaled{0};
        size_t unread{0};
        size_t empty{0};
        size_t zero{0};
        size_t good{0};
        size_t noHits{0};
    };

    struct DetectorSpecificEventMetrics {
        FragmentCounters raw;
        FragmentCounters zs;
        FragmentCounters ph;

        size_t rawFragsFlaggedBadOnly{0};
        size_t rawFragsFlaggedMissingOnly{0};
        size_t rawFragsFlaggedBadAndMissing{0};
        size_t rawFragsWithInvalidHeaders{0};
        size_t rawFragsWithInvalidAnchors{0};
        size_t zsFragsSkippedDueToRawFlag{0};
        size_t phFragsSkippedDueToRawFlag{0};
        size_t setsSkippedDueToRawFlag{0};
        size_t setsSkippedDueToInvalidHeaders{0};
        size_t zsFragsSkippedDueToInvalidRawHeader{0};
        size_t phFragsSkippedDueToInvalidRawHeader{0};
        size_t zsFragsSkippedDueToNoPrecedingRawHeader{0};
        size_t phFragsSkippedDueToNoPrecedingRawHeader{0};
        size_t phCountMismatch{0};

        size_t zsLengthMismatch{0};
        size_t zsRegionMismatch{0};

    };

}; // STMDigisFromFragments

// Configuration section
// =====================
STMDigisFromFragments::STMDigisFromFragments(const art::EDProducer::Table<Config>& config)
    :art::EDProducer{config}
    ,_stmFragmentsTag(config().stmTag())
    ,_verbosityLevel(config().verbosityLevel() ? *(config().verbosityLevel()) : 0)
    ,_saveSTMFragSummary(config().saveSTMFragSummary())
    ,_saveRawWaveformsWithHeaderHPGe(config().saveRawWaveformsWithHeaderHPGe())
    ,_saveRawWaveformsHPGe(config().saveRawWaveformsHPGe())
    ,_saveZSWaveformsHPGe(config().saveZSWaveformsHPGe())
    ,_saveRawWaveformsWithHeaderLaBr(config().saveRawWaveformsWithHeaderLaBr())
    ,_saveRawWaveformsLaBr(config().saveRawWaveformsLaBr())
    ,_saveZSWaveformsLaBr(config().saveZSWaveformsLaBr())
    // PH digis are always produced
{
    // Products depending on configuration switches
    if (_saveSTMFragSummary){
        produces<mu2e::STMFragmentSummaryCollection>("stmFragSummaryHPGe");
        produces<mu2e::STMFragmentSummaryCollection>("stmFragSummaryLaBr");
    }

    //change products to be maps of event headers
    // HPGe
    if (_saveRawWaveformsWithHeaderHPGe){produces<mu2e::STMWaveformDigiCollectionMap>("rawWithHeaderHPGe"); }
    if (_saveRawWaveformsHPGe){produces<mu2e::STMWaveformDigiCollectionMap>("rawHPGe"); }
    if (_saveZSWaveformsHPGe){produces<mu2e::STMWaveformDigiCollectionMap>("zsHPGe"); }
    produces<mu2e::STMPHDigiCollectionMap>("phHPGe");
    // LaBr
    if (_saveRawWaveformsWithHeaderLaBr){ produces<mu2e::STMWaveformDigiCollectionMap>("rawWithHeaderLaBr"); }
    if (_saveRawWaveformsLaBr){ produces<mu2e::STMWaveformDigiCollectionMap>("rawLaBr"); }
    if (_saveZSWaveformsLaBr){ produces<mu2e::STMWaveformDigiCollectionMap>("zsLaBr"); }
    produces<mu2e::STMPHDigiCollectionMap>("phLaBr");
}

// Start of Event processing
// =========================
void STMDigisFromFragments::produce(Event& event)
{
    // Extra frag counters for this event
    size_t containerFragsThisEvent{0};
    size_t innerFragsThisEvent{0};
    size_t containerFragsHPGeThisEvent{0};
    size_t containerFragsLaBrThisEvent{0};
    size_t innerFragsHPGeThisEvent{0};
    size_t innerFragsLaBrThisEvent{0};
    // Quality of frags
    size_t unknownFragsThisEvent{0};
    size_t badHPGeFragsThisEvent{0};
    size_t missingHPGeFragsThisEvent{0};
    size_t badLaBrFragsThisEvent{0};
    size_t missingLaBrFragsThisEvent{0};

    // Unknown Container
    size_t unknownContainersThisEvent{0};
    size_t nonContainersThisEvent{0};

    ++_totalEvents; // Increments total event counter

    RawHeaderState rawHeaderHPGe;
    RawHeaderState rawHeaderLaBr;
    DetectorSpecificEventMetrics LaBrEventMetrics;
    DetectorSpecificEventMetrics HPGeEventMetrics;

    // Set product ID

    // Frag Summaries
    std::unique_ptr<mu2e::STMFragmentSummaryCollection> stmFragSummaryHPGe(new mu2e::STMFragmentSummaryCollection);
    std::unique_ptr<mu2e::STMFragmentSummaryCollection> stmFragSummaryLaBr(new mu2e::STMFragmentSummaryCollection);

    // map start - keep explicit for now
    // HPGe
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>> rawWaveformDigisWithHeaderHPGe(new std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>);
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>> rawWaveformDigisHPGe(new std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>);
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>> zsWaveformDigisHPGe(new std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>);
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMPHDigiCollection>> phDigisHPGe(new std::map<mu2e::STMEventHeader, mu2e::STMPHDigiCollection>);
    // LaBr
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>> rawWaveformDigisWithHeaderLaBr(new std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>);
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>> rawWaveformDigisLaBr(new std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>);
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>> zsWaveformDigisLaBr(new std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>);
    std::unique_ptr<std::map<mu2e::STMEventHeader, mu2e::STMPHDigiCollection>> phDigisLaBr(new std::map<mu2e::STMEventHeader, mu2e::STMPHDigiCollection>);

    // booleans for tracking what detectors are present in the event
    bool eventHasHPGe{false};
    bool eventHasLaBr{false};

    // Get STM Fragments from the event, via handle
    art::Handle<artdaq::Fragments> STMFragmentsHandle;
    event.getByLabel(_stmFragmentsTag, STMFragmentsHandle);
    const auto& STMFragments = STMFragmentsHandle.product();
    uint16_t outerFragID{0};

    // Loop over outer fragments
    for (const auto& frag : *STMFragments) {
        ++_totalFragments;
        outerFragID = frag.fragmentID();
        if (_verbosityLevel > 0) {
            mf::LogDebug("STMDigisFromFragments")
            << "Processing outer fragment with ID : " << outerFragID << "\n";
        }

        // Check if this is a container fragment
        if (frag.type() == artdaq::Fragment::ContainerFragmentType) {

            mu2e::STMFragment container_frag(frag);
            artdaq::ContainerFragment cont_frag(frag);
            size_t blocks = cont_frag.block_count();

            if (container_frag.isHPGeContainer()) {
                if (_verbosityLevel > 2) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "Processing HPGe container fragment with ID : " << outerFragID << "\n";
                }
                ++_totalContainersHPGe; eventHasHPGe = true; _totalInnerFragsHPGe += blocks;
                innerFragsHPGeThisEvent += blocks;
                ++containerFragsHPGeThisEvent;
            } else if (container_frag.isLaBrContainer()) {
                if (_verbosityLevel > 2) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "Processing LaBr container fragment with ID : " << outerFragID << "\n";
                }
                ++_totalContainersLaBr; eventHasLaBr = true; _totalInnerFragsLaBr += blocks;
                innerFragsLaBrThisEvent += blocks;
                ++containerFragsLaBrThisEvent;
            } else {
                mf::LogWarning("STMDigisFromFragments")
                << "Encountered an unknown STM Container Fragment \n"
                << "Frag ID : "  << outerFragID << "\n"
                << "Event   : " << _totalEvents << "\n";

                ++_totalUnknownContainers;
                ++unknownContainersThisEvent;
                continue; // Skips the rest of this container frag if it is unknown
            }
            ++containerFragsThisEvent;
            ++_totalContainers;
            _totalInnerFrags += blocks;
            innerFragsThisEvent += blocks;

            // loop over container blocks
            for (size_t i = 0 ; i < cont_frag.block_count(); ++i) {
                auto inner_frag = cont_frag.at(i);
                mu2e::STMFragment stm_frag(*inner_frag);
                mu2e::STMWaveformDigi stm_waveform;

                if (stm_frag.isRaw()){
                    ++_totalRawFragsSeen;
                    // Determine which detector this raw fragment is in
                    bool const isHPGe = stm_frag.isHPGe();
                    bool const isLaBr = stm_frag.isLaBr();

                    if (!isHPGe && !isLaBr) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "Encountered raw fragment that is neither HPGe nor LaBr\n";

                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    isHPGe ? ++_totalRawFragsSeenHPGe : ++_totalRawFragsSeenLaBr;

                    auto& headerState = isHPGe ? rawHeaderHPGe : rawHeaderLaBr;
                    auto& eventMetrics  = isHPGe ? HPGeEventMetrics : LaBrEventMetrics;

                    // reset current raw header state
                    headerState = RawHeaderState{};
                    ++eventMetrics.raw.seen;
                    // reset eventHeader related variables
                    headerState.eventWindowTag = 0;
                    headerState.eventMode = 0;
                    headerState.adcClock = 0;
                    headerState.dtcClock = 0;

                    // check header+adcs is not less than expected header length (22)
                    if(stm_frag.dataWords() < stm::RawHeader::WORDS){
                        mf::LogWarning("STMDigisFromFragments")
                        << "Raw fragment has fewer data words than expected header length\n";

                        ++_totalRawFragsWithInvalidHeaders;
                        ++_totalUnreadInnerFrags;
                        ++eventMetrics.rawFragsWithInvalidHeaders;
                        isHPGe ? ++_totalRawFragsWithInvalidHeadersHPGe : ++_totalRawFragsWithInvalidHeadersLaBr;
                        ++eventMetrics.setsSkippedDueToInvalidHeaders;
                        headerState.skipCurrentSetDueToInvalidHeader = true;
                        continue;
                    }

                    // Confirm Raw Header has Valid Anchors
                    // For now just use as a diagnostic tool
                    if (!stm_frag.hasValidAnchors()) {
                        if (_verbosityLevel > 1) {
                            std::ostringstream msg;
                            auto const* temp_dataPtr = stm_frag.dataBegin();
                            msg << std::hex // switch to hexadecimal
                            << "\nStart of raw header : 0x" << static_cast<uint16_t>(temp_dataPtr[stm::RawHeader::ANCHOR_START])
                            << "\nEnd of raw header   : 0x" << static_cast<uint16_t>(temp_dataPtr[stm::RawHeader::ANCHOR_END])
                            << std::dec << "\n"; // switch to decimal

                            mf::LogDebug("STMDigisFromFragments")
                            << msg.str();
                        }

                        ++eventMetrics.rawFragsWithInvalidAnchors;
                        ++_totalRawFragsWithInvalidAnchors;
                        isHPGe ? ++_totalRawFragsWithInvalidAnchorsHPGe : ++_totalRawFragsWithInvalidAnchorsLaBr;
                    }

                    // Check if bad or missing or both
                    if (stm_frag.badData() && stm_frag.missing()) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "Raw Fragment flagged as bad and missing\n";

                        ++eventMetrics.rawFragsFlaggedBadAndMissing;
                        ++_totalRawFragsFlaggedBadAndMissing;
                        isHPGe ? ++_totalRawFragsFlaggedBadAndMissingHPGe : ++_totalRawFragsFlaggedBadAndMissingLaBr;
                    } else if (stm_frag.badData()) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "Raw Fragment flagged as bad\n";

                        ++eventMetrics.rawFragsFlaggedBadOnly;
                        ++_totalRawFragsFlaggedBadOnly;
                        isHPGe ? ++_totalRawFragsFlaggedBadOnlyHPGe : ++_totalRawFragsFlaggedBadOnlyLaBr;
                    } else if (stm_frag.missing()) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "Raw Fragment flagged as missing\n";

                        ++eventMetrics.rawFragsFlaggedMissingOnly;
                        ++_totalRawFragsFlaggedMissingOnly;
                        isHPGe ? ++_totalRawFragsFlaggedMissingOnlyHPGe : ++_totalRawFragsFlaggedMissingOnlyLaBr;
                    }

                    bool const badOrMissing = stm_frag.badData() || stm_frag.missing();
                    if (badOrMissing) {
                        headerState.skipCurrentSetDueToRawFlags = true;
                        ++eventMetrics.setsSkippedDueToRawFlag;
                        continue;
                    }

                    // Check if raw fragment is prescaled to determine whether to skip checks
                    // If the fragment is prescaled we do not save the waveform
                    // If the fragment is prescaled we also do not want to check its payload
                    // We do want relevant infrormation for the zs and ph
                    headerState.rawPrescaled = stm_frag.rawPrescaled();
                    headerState.rawPrescaleValue = stm_frag.rawPrescaleValue();

                    // Extract rest of the raw header information
                    // Even if raw is prescaled we still want to save zs info just in case
                    // before checking whats in the payload
                    headerState.containsZSInfo = true;
                    headerState.containsPHInfo = true;
                    headerState.zsPrescaled = stm_frag.zsPrescaled();
                    headerState.zsPrescaleValue = stm_frag.zsPrescaleValue();
                    headerState.expectedZSLength = stm_frag.zsLength();
                    headerState.expectedZSRegions = stm_frag.zsRegions();
                    headerState.expectedPHCount = stm_frag.phCount();

                    // extract EWT, mode/spillFlags, adcClock, dtcClock
                    headerState.eventWindowTag = stm_frag.eventWindowTag();
                    headerState.eventMode = stm_frag.spillFlag();
                    headerState.adcClock = stm_frag.adcClock();
                    headerState.dtcClock = stm_frag.dtcClock();

                    // create stmEventHeader
                    mu2e::STMEventHeader stm_event_header(
                        headerState.eventWindowTag,
                        headerState.eventMode,
                        headerState.adcClock,
                        headerState.dtcClock);

                    if (headerState.rawPrescaled) {
                        if(_verbosityLevel > 2) {
                            mf::LogDebug("STMDigisFromFragments")
                            << "Raw Fragment is prescaled";
                        }

                        ++eventMetrics.raw.prescaled;
                        ++_totalRawFragsPrescaled;
                        isHPGe ? ++_totalRawFragsPrescaledHPGe : ++_totalRawFragsPrescaledLaBr;

                        continue;
                    }

                    auto payloadPtr = stm_frag.payloadBegin();
                    auto payloadWords = stm_frag.payloadWords();

                    // 22 Header Words + adcs < = dataWords
                    // This check comes after ensuring raw is not prescaled
                    if (stm::RawHeader::WORDS + payloadWords > stm_frag.dataWords()) {
                        // this is an invalid-header case
                        mf::LogWarning("STMDigisFromFragments")
                        << "Raw Fragment is reading more than what fragment stores";

                        ++_totalRawFragsWithInvalidHeaders;
                        isHPGe ? ++_totalRawFragsWithInvalidHeadersHPGe : ++_totalRawFragsWithInvalidHeadersLaBr;
                        ++eventMetrics.rawFragsWithInvalidHeaders;
                        ++eventMetrics.setsSkippedDueToInvalidHeaders;
                        ++_totalUnreadInnerFrags;
                        headerState.skipCurrentSetDueToInvalidHeader = true;
                        continue;
                    }

                    // Ideally by here we know the raw frag is good and its payload can be checked
                    // Check if raw frag is empty
                    bool allZeros = true;
                    if (payloadWords == 0) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nFound an empty raw fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--Raw Frag\n";

                        //increment
                        ++eventMetrics.raw.empty;
                        ++_totalEmptyRawFrags;
                        isHPGe ? ++_totalEmptyRawFragsHPGe : ++_totalEmptyRawFragsLaBr;
                        continue;
                    }
                    // Check if any payload words exist
                    for (size_t k = 0; k < payloadWords; ++k){
                        if (payloadPtr[k] != 0) {
                            allZeros = false;
                            break;
                        }
                    }
                    if (allZeros) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nFound an all-zero raw fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--Raw Frag\n";
                        //counter increment for all-zero raw fragments
                        ++eventMetrics.raw.zero;
                        ++_totalZeroRawFrags;
                        isHPGe? ++_totalZeroRawFragsHPGe: ++_totalZeroRawFragsLaBr;
                        continue;
                    }
                    // At this point the raw fragment has payload words and is non zero filled
                    // We continue processing this raw fragment
                    if (_verbosityLevel > 2) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "\nFound a good raw fragment at Event : " << _totalEvents << "\n"
                        << "Fragment Index : " << i << "\n"
                        << "--Raw Frag\n";
                    }

                    // Print first few payload words for inspection
                    if (_verbosityLevel > 3) {
                        std::ostringstream msg;
                        msg << "First few payload words(adcs) for inspection : " ;
                        for (size_t w = 0; w < std::min(payloadWords, static_cast<size_t>(10)); ++w){
                            msg << payloadPtr[w] << " ";
                        }
                        mf::LogDebug("STMDigisFromFragments")
                        << msg.str();
                    }
                    // Raw Header inspection
                    if (_verbosityLevel > 4) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "\nRaw Header for inspection \n"
                        << "Raw Length       : " << stm_frag.rawLength() << "\n"
                        << "ZS  Length       : " << stm_frag.zsLength() << "\n"
                        << "ZS  Regions      : " << stm_frag.zsRegions() << "\n"
                        << "Event            : " << _totalEvents << "\n"
                        << "Inner Frag Index : " << i << "\n";
                    }
                    // From here good raw frags end up and get saved
                    ++eventMetrics.raw.good;
                    ++_totalGoodRawFrags;
                    isHPGe ? ++_totalGoodRawFragsHPGe : ++_totalGoodRawFragsLaBr;


                    // Save waveforms based on fcl configurations
                    // Save Raw Waveform With Header Info - HPGe
                    if (_saveRawWaveformsWithHeaderHPGe && isHPGe){
                        auto dataPtr = stm_frag.dataBegin();
                        auto dataWords = stm_frag.dataWords();
                        // set stm_waveform
                        stm_waveform.set_data(dataWords,dataPtr);
                        // use map
                        (*rawWaveformDigisWithHeaderHPGe)[stm_event_header].push_back(stm_waveform);
                    }
                    // Save Raw Waveform Map - HPGe
                    if (_saveRawWaveformsHPGe && isHPGe){
                        // set stm_waveform
                        stm_waveform.set_data(payloadWords, payloadPtr);

                        // Give map a key and value, key is stmHeader,
                        // *rawwaveformDigisHPGe is the map
                        // stm_event_header is the collection stored under this header
                        // .push_back(stm_waveform) adds waveform to the collection
                        (*rawWaveformDigisHPGe)[stm_event_header].push_back(stm_waveform);
                    }
                    // Save Raw Waveform With Header Info - LaBr
                    if (_saveRawWaveformsWithHeaderLaBr && isLaBr){
                        auto dataPtr = stm_frag.dataBegin();
                        auto dataWords = stm_frag.dataWords();
                        // set stm_waveform
                        stm_waveform.set_data(dataWords,dataPtr);
                        // set map here
                        (*rawWaveformDigisWithHeaderLaBr)[stm_event_header].push_back(stm_waveform);
                    }
                    // Save Raw Waveform Map - LaBr
                    if (_saveRawWaveformsLaBr && isLaBr){
                        // set stm_waveform
                        stm_waveform.set_data(payloadWords, payloadPtr);
                        // set map here
                        (*rawWaveformDigisLaBr)[stm_event_header].push_back(stm_waveform);
                    }

                    if (_verbosityLevel > 2) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "Raw with frag index" << i
                        << ", detector = " << (isHPGe? "HPGe": "LaBr") << std::boolalpha
                        << ", prescaled = " << headerState.rawPrescaled
                        << "\n";
                    }
                }// End of Raw fragment check
                else if (stm_frag.isZS()){
                    ++_totalZSFragsSeen;
                    // Determine which detector this ZS fragment belongs to
                    bool const isHPGe = stm_frag.isHPGe();
                    bool const isLaBr = stm_frag.isLaBr();

                    // Double check that the fragment belongs to one of the expected detectors
                    if (!isHPGe && !isLaBr) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "Encountered a ZS fragment that is neither HPGe nor LaBr";

                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    auto& headerState = isHPGe ? rawHeaderHPGe : rawHeaderLaBr;
                    auto& eventMetrics  = isHPGe ? HPGeEventMetrics : LaBrEventMetrics;

                    isHPGe ? ++_totalZSFragsSeenHPGe : ++_totalZSFragsSeenLaBr;
                    ++eventMetrics.zs.seen;

                    // Check
                    if (_verbosityLevel > 2) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "ZS frag index = " << i
                        << ", detector = " << (isHPGe ? "HPGe" : "LaBr") << std::boolalpha
                        << "\n";
                    }

                    // Decide Here to skip based on previous raw fragment information
                    if (headerState.skipCurrentSetDueToInvalidHeader) {
                        ++eventMetrics.zsFragsSkippedDueToInvalidRawHeader;
                        continue;
                    }

                    if (headerState.skipCurrentSetDueToRawFlags) {
                        isHPGe ? ++_totalZSFragsSkippedDueToRawFlagHPGe : ++_totalZSFragsSkippedDueToRawFlagLaBr;
                        ++eventMetrics.zsFragsSkippedDueToRawFlag;
                        continue;
                    }

                    // Check that there was a raw header before this ZS fragment
                    if (!headerState.containsZSInfo) {
                        ++eventMetrics.zsFragsSkippedDueToNoPrecedingRawHeader;
                        continue;
                    }

                    // Extract zs variables from Raw Header
                    bool zsInfoWasExtracted = headerState.containsZSInfo;
                    bool zsPrescaled = headerState.zsPrescaled;
                    bool rawPrescaled = headerState.rawPrescaled;
                    uint16_t zsRegions = headerState.expectedZSRegions;
                    uint16_t zsLength = headerState.expectedZSLength;
                    uint16_t rawPrescaleValue = headerState.rawPrescaleValue;
                    uint16_t zsPrescaleValue = headerState.zsPrescaleValue;

                    // Extract for eventHeader (EWT, mode/spillFlag, adcClock, dtcClock)
                    uint64_t zsEWT = headerState.eventWindowTag;
                    uint8_t zsMode = headerState.eventMode;
                    uint64_t zsAdcClock = headerState.adcClock;
                    uint64_t zsdtcClock = headerState.dtcClock;

                    // create stmEventHeader
                    mu2e::STMEventHeader stm_event_header(
                        zsEWT,
                        zsMode,
                        zsAdcClock,
                        zsdtcClock);

                    // Decide Here to skip based on zs prescale information
                    if (zsPrescaled) {
                        if(_verbosityLevel > 2) {
                            mf::LogDebug("STMDigisFromFragments")
                            << "ZS Fragment is prescaled\n";
                        }

                        isHPGe ? ++_totalZSFragsPrescaledHPGe : ++_totalZSFragsPrescaledLaBr;
                        ++_totalZSFragsPrescaled;
                        ++eventMetrics.zs.prescaled;
                        continue;
                    }

                    auto payloadPtr = stm_frag.payloadBegin();
                    auto payloadWords = stm_frag.payloadWords();

                    // Check if Empty
                    if (payloadWords == 0) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nFound an empty zs fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--ZS Frag\n";

                        ++_totalEmptyZSFrags;
                        isHPGe ? ++_totalEmptyZSFragsHPGe : ++_totalEmptyZSFragsLaBr;
                        ++eventMetrics.zs.empty;
                        continue;
                    }

                    // At this point the zs fragmant is non-empty
                    // May contain valid data

                    // Print first 20 payload adcs for inspection
                    if (_verbosityLevel > 3) {
                        std::ostringstream msg;
                        msg << "\nFirst few payload words for inspection : ";
                        for (size_t w = 0; w < std::min(payloadWords, static_cast<size_t>(10)); ++w) {
                            msg << payloadPtr[w] << " ";
                        }
                        mf::LogDebug("STMDigisFromFragments")
                        << msg.str();
                    }

                    // Definitions for payload references
                    auto const* dataPtr = stm_frag.dataBegin();
                    auto const dataWords = stm_frag.dataWords();
                    auto dataEndPtr = dataPtr + dataWords;
                    size_t regionCounter {0};
                    size_t zsTotalLengthCalculated {0};
                    uint16_t lastZSIndexRecorded {0};
                    uint16_t lastZSLengthRecorded {0};
                    bool malformedZS {false};
                    std::vector<ZSRegion> regions;

                    // Check if data is zero filled using parsing for ZS
                    bool sawADC{false};
                    bool allADCSZero{true};
                    // May need to be adjusted, we would only read the first pulse

                    if (_verbosityLevel > 5) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "\nData Words for inspection : " << "\n"
                        << "Data Words % 4 : " << dataWords % 4 << "\n"
                        << "--ZS Frag\n";
                    }

                    while (dataPtr + 2 <= dataEndPtr){
                      if (zsInfoWasExtracted && regionCounter >= zsRegions) {
                        break;
                      }
                        uint16_t currentZSIndex = static_cast<uint16_t>(dataPtr[0]);
                        uint16_t currentZSLength = static_cast<uint16_t>(dataPtr[1]);
                        auto adc = dataPtr + 2;

                        if (currentZSLength > static_cast<size_t>(dataEndPtr - adc)) {
                            malformedZS = true;
                            break;
                        }

                        // Check if any is zero filled payload for this region
                        for (size_t sample = 0; sample < currentZSLength; ++sample){
                            sawADC = true;
                            if (adc[sample] !=0) {
                                allADCSZero = false;
                            }
                        }
                        // Record the region
                        bool const saveZS = isHPGe ? _saveZSWaveformsHPGe : _saveZSWaveformsLaBr;
                        if (saveZS){
                            regions.push_back({currentZSIndex, std::vector<int16_t>(adc, adc + currentZSLength)});
                        }

                        uint32_t trigTimeOffset = currentZSIndex;

                        // Print Check per segment
                        if (_verbosityLevel > 5) {
                            mf::LogDebug("STMDigisFromFragments")
                            << "\nZS Segment Check :" << "\n"
                            << " , ZS Index  = " << currentZSIndex
                            << " , ZS Length = " << currentZSLength
                            << " , ZS Offset = " << trigTimeOffset << "\n";
                        }
                        // Update Variables
                        lastZSIndexRecorded = currentZSIndex;
                        lastZSLengthRecorded = currentZSLength;
                        zsTotalLengthCalculated += lastZSLengthRecorded;
                        ++regionCounter;
                        dataPtr = adc + currentZSLength;
                    } // End of While loop

                    if (_verbosityLevel > 4){
                        //Summary
                        mf::LogDebug("STMDigisFromFragments")
                        << "ZS Summary : "
                        << " , last ZS Index = " << lastZSIndexRecorded
                        << " , last ZS Length = " << lastZSLengthRecorded
                        << " , ZS Total Length = " << zsTotalLengthCalculated
                        << " , Region Counter = " << regionCounter << "\n"
                        << ", Frag Index = " << i << "\n";
                    }

                    // In case ZS has some weird behavior
                    if (malformedZS) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nMalformed ZS detected at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n";

                        ++eventMetrics.zs.unread;
                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    // Add exceptions
                    if (zsInfoWasExtracted) {
                        if (zsLength != zsTotalLengthCalculated) {
                            mf::LogWarning("STMDigisFromFragments")
                            << "\n=== ZS LENGTH MISMATCH ===\n"
                            << "ZS Length Count from Raw Header : " << zsLength << "\n"
                            << "ZS Length Calculated : " << zsTotalLengthCalculated << "\n"
                            // General Information about where error was found
                            << "Found at Event : " << _totalEvents << "\n"
                            << "Found at Frag Index : " << i << "\n"
                            << "Raw Prescaled : " << (rawPrescaled ? "Yes" : "No") << "\n"
                            << "Raw Prescale Value : " << rawPrescaleValue << "\n"
                            << "ZS Prescaled : " << (zsPrescaled ? "Yes" : "No") << "\n"
                            << "ZS Prescale Value : " << zsPrescaleValue << "\n"
                            << "Found at HPGe Container Frag : " << (isHPGe ? "Yes" : "No") << "\n"
                            << "Found at LaBr Container Frag : " << (isLaBr ? "Yes" : "No") << "\n";

                            ++_totalZSLengthMismatch;
                            isHPGe? ++_totalZSLengthMismatchHPGe: ++_totalZSLengthMismatchLaBr;
                            ++eventMetrics.zsLengthMismatch;
                            continue;
                        }
                        if ( zsRegions != regionCounter) {
                            mf::LogWarning("STMDigisFromFragments")
                            << "\n=== ZS REGION COUNT MISMATCH ===\n"
                            << "ZS Region Count from Raw Header : " << zsRegions << "\n"
                            << "ZS Region Count Calculated : " << regionCounter << "\n"
                            // General Information about where error was found
                            << "Found at Event : " << _totalEvents << "\n"
                            << "Found at Frag Index : " << i << "\n"
                            << "Raw Prescaled : " << (rawPrescaled ? "Yes" : "No") << "\n"
                            << "Raw Prescale Value : " << rawPrescaleValue << "\n"
                            << "ZS Prescaled : " << (zsPrescaled ? "Yes" : "No") << "\n"
                            << "ZS Prescale Value : " << zsPrescaleValue << "\n"
                            << "Found at HPGe Container Frag : " << (isHPGe ? "Yes" : "No") << "\n"
                            << "Found at LaBr Container Frag : " << (isLaBr ? "Yes" : "No") << "\n";

                            ++_totalZSRegionMismatch;
                            isHPGe? ++_totalZSRegionMismatchHPGe: ++_totalZSRegionMismatchLaBr;
                            ++eventMetrics.zsRegionMismatch;
                            continue;
                        }
                    }

                    if (sawADC && allADCSZero) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nFound a zero-filled zs fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--ZS Frag\n";

                        ++eventMetrics.zs.zero;
                        ++_totalZeroZSFrags;
                        isHPGe ? ++_totalZeroZSFragsHPGe : ++_totalZeroZSFragsLaBr;
                        continue;
                    }

                    // At this point we know the fragment was good
                    if (_verbosityLevel > 2) {
                        mf::LogDebug("STMDigisFromFragments")
                        <<"\nFound a good ZS fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--ZS Frag\n";
                    }
                    for (auto& region : regions) {
                        // set waveform digi with offset and adcs
                        mu2e::STMWaveformDigi zsDigi(region.offset, region.adcs);

                        // Emplacing
                        if (isHPGe && _saveZSWaveformsHPGe) {
                          // set map for HPGe
                            (*zsWaveformDigisHPGe)[stm_event_header].push_back(zsDigi);
                        } else if (isLaBr && _saveZSWaveformsLaBr) {
                            // set map for LaBr
                           (*zsWaveformDigisLaBr)[stm_event_header].push_back(zsDigi);
                        }
                    }

                    // Increment remaining counters
                    ++_totalGoodZSFrags;
                    isHPGe ? ++_totalGoodZSFragsHPGe : ++_totalGoodZSFragsLaBr;
                    ++eventMetrics.zs.good;

                }// End of ZS fragmment check
                else if (stm_frag.isPH()){
                    ++_totalPHFragsSeen;
                    // Determine which detector this PH fragment belongs to
                    bool const isHPGe = stm_frag.isHPGe();
                    bool const isLaBr = stm_frag.isLaBr();

                    if (!isHPGe && !isLaBr) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "Encountered PH fragment that is neither HPGe nor LaBr";

                        ++_totalUnreadInnerFrags;
                        continue;
                    }
                    auto& headerState = isHPGe ? rawHeaderHPGe : rawHeaderLaBr;
                    auto& eventMetrics  = isHPGe ? HPGeEventMetrics : LaBrEventMetrics;

                    isHPGe ? ++_totalPHFragsSeenHPGe : ++_totalPHFragsSeenLaBr;
                    ++eventMetrics.ph.seen;

                    // Skip if header was malformed
                    if (headerState.skipCurrentSetDueToInvalidHeader) {
                        ++eventMetrics.phFragsSkippedDueToInvalidRawHeader;
                        continue;
                    }

                    // Skip if Raw Fragment was Bad or Missing
                    if (headerState.skipCurrentSetDueToRawFlags) {
                        ++eventMetrics.phFragsSkippedDueToRawFlag;
                        isHPGe ? ++_totalPHFragsSkippedDueToRawFlagHPGe : ++_totalPHFragsSkippedDueToRawFlagLaBr;
                        continue;
                    }

                    // Check if a raw header was extracted before this PH fragment
                    if (!headerState.containsPHInfo) {
                        ++eventMetrics.phFragsSkippedDueToNoPrecedingRawHeader;
                        continue;
                    }

                    // Extract ph varibales from Raw Header
                    uint16_t phCount = headerState.expectedPHCount;

                    // Extract for eventHeader (EWT, mode/spillFlag, adcClock, dtcClock)
                    uint64_t phEWT = headerState.eventWindowTag;
                    uint8_t phMode = headerState.eventMode;
                    uint64_t phAdcClock = headerState.adcClock;
                    uint64_t phdtcClock = headerState.dtcClock;

                    // Construct stmEventHeader for PH fragment
                    mu2e::STMEventHeader stm_event_header(
                        phEWT,
                        phMode,
                        phAdcClock,
                        phdtcClock);

                    // Check if PH fragment is empty
                    auto payloadPtr = stm_frag.payloadBegin();
                    auto payloadWords = stm_frag.payloadWords();
                    bool allPHAreZeros= true;

                    if (phCount == 0) {
                        // no hits reported by raw header
                        if (_verbosityLevel > 2) {
                            mf::LogDebug("STMDigisFromFragments")
                            << "\nNo hits reported for this fragment at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "--PH Frag\n";
                        }
                        ++eventMetrics.ph.noHits;
                        ++_totalNoHitsPHFrags;
                        isHPGe ? ++_totalNoHitsPHFragsHPGe : ++_totalNoHitsPHFragsLaBr;
                        continue;
                    }

                    if (payloadWords == 0) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nFound an empty PH fragment at Event : " <<  _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--PH Frag\n";
                        // increment
                        // fragment has no payload
                        ++eventMetrics.ph.empty;
                        ++_totalEmptyPHFrags;
                        isHPGe ? ++_totalEmptyPHFragsHPGe : ++_totalEmptyPHFragsLaBr;
                        continue;
                    }

                    // Check if we under count ph Pairs in fragment comapred to raw header count
                    size_t expectedWords = 2 * static_cast<size_t>(phCount);
                    if (payloadWords < expectedWords) {
                        // Note payloadWords can include data padding
                        // Ideally it would not be less than ph pair count
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nPH Payload is shorter than Raw Header Expects\n"
                        << "Expected Words : " << expectedWords << "\n"
                        << "Payload Words  : " << payloadWords << "\n";

                        // report
                        ++_totalPHCountMismatch;
                        isHPGe ? ++_totalPHCountMismatchHPGe : ++_totalPHCountMismatchLaBr;
                        ++eventMetrics.phCountMismatch;
                        ++eventMetrics.ph.unread;
                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    // phCount becomes our upper bound
                    size_t nPairsToRead = phCount;

                    // Check if zero filled
                    // At this stage we expect PHCount to be a valid estimate in frag
                    // But if all ph's second entry return 0 then this is zero filled
                    // Since its a (time,PH) pair we will check every second entry for the PH value
                    for (size_t k = 0; k < nPairsToRead; ++k) {
                        if (payloadPtr[k*2+1] !=0) {
                            allPHAreZeros = false;
                            break;
                        }
                    }

                    if (allPHAreZeros) {
                        mf::LogWarning("STMDigisFromFragments")
                        << "\nFound a zero-filled PH Fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "---PH Frag\n";

                        // counter
                        ++eventMetrics.ph.zero;
                        ++_totalZeroPHFrags;
                        isHPGe ? ++_totalZeroPHFragsHPGe : ++_totalZeroPHFragsLaBr;
                        continue;
                    }

                    if (_verbosityLevel > 4) {
                        std::ostringstream msg;
                        msg << "\nFirst few payload words: ";

                        for (size_t k = 0; k < std::min<size_t>(payloadWords,10); ++k) {
                            msg << payloadPtr[k] << " ";
                        }

                        mf::LogDebug("STMDigisFromFragments")
                        << msg.str();
                    }

                    for (size_t i_phPair = 0; i_phPair < nPairsToRead ; ++i_phPair){
                        size_t i_PH = 2 * i_phPair;
                        uint32_t time = static_cast<uint16_t>(payloadPtr[i_PH]);
                        int16_t const pulseHeight = payloadPtr[i_PH + 1];
                        mu2e::STMPHDigi PH_digi(time,pulseHeight);

                        // Emplace Back (Always On)
                        if (isHPGe) {
                          (*phDigisHPGe)[stm_event_header].emplace_back(PH_digi);
                        }
                        if (isLaBr) {
                          (*phDigisLaBr)[stm_event_header].emplace_back(PH_digi);
                        }
                    }

                    // At this point we have a good PH fragment
                    if (_verbosityLevel > 2) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "\nFound a good PH Fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "---PH Frag\n";
                    }

                    // Last counter increment
                    ++_totalGoodPHFrags;
                    isHPGe ? ++_totalGoodPHFragsHPGe : ++_totalGoodPHFragsLaBr;
                    ++eventMetrics.ph.good;

                }// End of PH fragment check
                else {
                    // fallback for unreadable fragment
                    ++_totalUnreadInnerFrags; //Job Summary Counter
                    ++unknownFragsThisEvent;

                    mf::LogWarning("STMDigisFromFragments")
                    << "Encountered an unreadable inner fragment " << "\n"
                    << "Frag Index  : " << i << "\n"
                    << "Frag ID     : " << inner_frag->fragmentID() << "\n"
                    << "Event       : " << _totalEvents << "\n";

                }// End of else non-raw/zs/pg frag
            }
        } else {
            // fallback for non-container fragment
            mf::LogWarning("STMDigisFromFragments")
            << "\nEncountered a non-container fragment " << "\n"
            << "Fragment ID : " << frag.fragmentID() << "\n"
            << "Event       : " << _totalEvents << "\n";

            ++_totalNonContainers;
            ++nonContainersThisEvent;
            continue;
        }
    } // End of frags loop

    // Update event type counters based on the detectors in event
    if (eventHasHPGe && eventHasLaBr) {
        ++_totalEventsWithBothDetectors;
    } else if (eventHasHPGe) {
        ++_totalEventsWithOnlyHPGe;
    } else if (eventHasLaBr) {
        ++_totalEventsWithOnlyLaBr;
    } else {
        ++_totalEventsWithNeitherDetector;
    }

    // HPGe Remaining Fragment Counts
    badHPGeFragsThisEvent = HPGeEventMetrics.rawFragsFlaggedBadOnly
    + HPGeEventMetrics.rawFragsFlaggedBadAndMissing;

    missingHPGeFragsThisEvent = HPGeEventMetrics.rawFragsFlaggedMissingOnly
    + HPGeEventMetrics.rawFragsFlaggedBadAndMissing;

    // LaBr Remaining Fragment Counts
    badLaBrFragsThisEvent = LaBrEventMetrics.rawFragsFlaggedBadOnly
    + LaBrEventMetrics.rawFragsFlaggedBadAndMissing;

    missingLaBrFragsThisEvent = LaBrEventMetrics.rawFragsFlaggedMissingOnly
    + LaBrEventMetrics.rawFragsFlaggedBadAndMissing;

    // Save STMFrag Summary here
    if (_saveSTMFragSummary) {
        if (_verbosityLevel > 2) {
            mf::LogDebug("STMDigisFromFragments")
            << "\nSaving STMFragSummary\n";
        }
        stmFragSummaryHPGe->emplace_back(
            containerFragsHPGeThisEvent,innerFragsHPGeThisEvent,
            badHPGeFragsThisEvent, missingHPGeFragsThisEvent,
            HPGeEventMetrics.zsFragsSkippedDueToRawFlag, HPGeEventMetrics.phFragsSkippedDueToRawFlag,
            HPGeEventMetrics.raw.prescaled, HPGeEventMetrics.zs.prescaled,
            HPGeEventMetrics.raw.good, HPGeEventMetrics.zs.good, HPGeEventMetrics.ph.good,
            HPGeEventMetrics.raw.zero, HPGeEventMetrics.zs.zero, HPGeEventMetrics.ph.zero,
            HPGeEventMetrics.raw.empty, HPGeEventMetrics.zs.empty, HPGeEventMetrics.ph.empty
        );

        stmFragSummaryLaBr->emplace_back(
            containerFragsLaBrThisEvent,innerFragsLaBrThisEvent,
            badLaBrFragsThisEvent, missingLaBrFragsThisEvent,
            LaBrEventMetrics.zsFragsSkippedDueToRawFlag, LaBrEventMetrics.phFragsSkippedDueToRawFlag,
            LaBrEventMetrics.raw.prescaled, LaBrEventMetrics.zs.prescaled,
            LaBrEventMetrics.raw.good, LaBrEventMetrics.zs.good, LaBrEventMetrics.ph.good,
            LaBrEventMetrics.raw.zero, LaBrEventMetrics.zs.zero, LaBrEventMetrics.ph.zero,
            LaBrEventMetrics.raw.empty, LaBrEventMetrics.zs.empty, LaBrEventMetrics.ph.empty
        );
    }

    // Get Number of Raw, ZS Waveforms saved and PH Digis saved for this event
    auto totalDigis = [](const auto& digiMap) {
        size_t total = 0;
        for (const auto& entry : digiMap) {
            total += entry.second.size();
        }
        return total;
    };

    // Final Move
    if (_verbosityLevel > 1) {
        // Event Summary -> tells us what happened per event
        mf::LogInfo("STMDigisFromFragments")
        << "\n========== STM EVENT SUMMARY - (Unpacking Module) ==========\n"
        << "\n--- Module Configuration For Job ---\n"
        << "Raw Waveforms with Header HPGe   : " << (_saveRawWaveformsWithHeaderHPGe ? "Yes" : "No") << "\n"
        << "Raw Waveforms HPGe               : " << (_saveRawWaveformsHPGe ? "Yes" : "No") << "\n"
        << "ZS Waveforms HPGe                : " << (_saveZSWaveformsHPGe ? "Yes" : "No") << "\n"
        << "PH Digis HPGe                    : Yes\n"
        << "Raw Waveforms with Header LaBr   : " << (_saveRawWaveformsWithHeaderLaBr ? "Yes" : "No") << "\n"
        << "Raw Waveforms LaBr               : " << (_saveRawWaveformsLaBr ? "Yes" : "No") << "\n"
        << "ZS Waveforms LaBr                : " << (_saveZSWaveformsLaBr ? "Yes" : "No") << "\n"
        << "PH Digis LaBr                    : Yes\n"

        << "\n--- Products Saved Per Event ---\n"
        << "Extracted Raw Waveforms With Header (HPGe)    : " << totalDigis(*rawWaveformDigisWithHeaderHPGe) << "\n"
        << "Extracted Raw Waveforms (HPGe)                : " << totalDigis(*rawWaveformDigisHPGe) << "\n"
        << "Extracted ZS  Waveforms (HPGe)                : " << totalDigis(*zsWaveformDigisHPGe) << "\n"
        << "Extracted PH  Digis     (HPGe)                : " << totalDigis(*phDigisHPGe) << "\n"
        << "\n"
        << "Extracted Raw Waveforms With Header (LaBr)    : " << totalDigis(*rawWaveformDigisWithHeaderLaBr) << "\n"
        << "Extracted Raw Waveforms (LaBr)                : " << totalDigis(*rawWaveformDigisLaBr) << "\n"
        << "Extracted ZS  Waveforms (LaBr)                : " << totalDigis(*zsWaveformDigisLaBr) << "\n"
        << "Extracted PH  Digis     (LaBr)                : " << totalDigis(*phDigisLaBr) << "\n"

        << "\n--- HPGe Summary Per Event---\n"
        << "Good  Raw Frags (HPGe)                        : " << HPGeEventMetrics.raw.good << "\n"
        << "Good  ZS  Frags (HPGe)                        : " << HPGeEventMetrics.zs.good << "\n"
        << "Good  PH  Frags (HPGe)                        : " << HPGeEventMetrics.ph.good << "\n"
        << "\n"
        << "Empty Raw Frags (HPGe)                        : " << HPGeEventMetrics.raw.empty << "\n"
        << "Empty ZS  Frags (HPGe)                        : " << HPGeEventMetrics.zs.empty << "\n"
        << "Empty PH  Frags (HPGe)                        : " << HPGeEventMetrics.ph.empty << "\n"
        << "\n"
        << "Zero  Raw Frags (HPGe)                        : " << HPGeEventMetrics.raw.zero << "\n"
        << "Zero  ZS  Frags (HPGe)                        : " << HPGeEventMetrics.zs.zero << "\n"
        << "Zero  PH  Frags (HPGe)                        : " << HPGeEventMetrics.ph.zero << "\n"
        << "\n"
        << "No Hits - PH Frags (HPGe)                     : " << HPGeEventMetrics.ph.noHits << "\n"
        << "\n"
        << "Bad Raw Frags (HPGe)                          : " << HPGeEventMetrics.rawFragsFlaggedBadOnly << "\n"
        << "Missing Raw Frags (HPGe)                      : " << HPGeEventMetrics.rawFragsFlaggedMissingOnly << "\n"
        << "Bad and Missing Raw Frags (HPGe)              : " << HPGeEventMetrics.rawFragsFlaggedBadAndMissing << "\n"
        << "\n"

        << "\n--- LaBr Summary Per Event ---\n"
        << "Good  Raw Frags (LaBr)                        : " << LaBrEventMetrics.raw.good << "\n"
        << "Good  ZS  Frags (LaBr)                        : " << LaBrEventMetrics.zs.good << "\n"
        << "Good  PH  Frags (LaBr)                        : " << LaBrEventMetrics.ph.good << "\n"
        << "\n"
        << "Empty Raw Frags (LaBr)                        : " << LaBrEventMetrics.raw.empty << "\n"
        << "Empty ZS  Frags (LaBr)                        : " << LaBrEventMetrics.zs.empty << "\n"
        << "Empty PH  Frags (LaBr)                        : " << LaBrEventMetrics.ph.empty << "\n"
        << "\n"
        << "Zero  Raw Frags (LaBr)                        : " << LaBrEventMetrics.raw.zero << "\n"
        << "Zero  ZS  Frags (LaBr)                        : " << LaBrEventMetrics.zs.zero << "\n"
        << "Zero  PH  Frags (LaBr)                        : " << LaBrEventMetrics.ph.zero << "\n"
        << "\n"
        << "No Hits - PH Frags (LaBr)                     : " << LaBrEventMetrics.ph.noHits << "\n"
        << "\n"
        << "Bad Raw Frags (LaBr)                          : " << LaBrEventMetrics.rawFragsFlaggedBadOnly << "\n"
        << "Missing Raw Frags (LaBr)                      : " << LaBrEventMetrics.rawFragsFlaggedMissingOnly << "\n"
        << "Bad and Missing Raw Frags (LaBr)              : " << LaBrEventMetrics.rawFragsFlaggedBadAndMissing << "\n"
        // Extra Filters
        << "\n=== Extra Filters Per Event ===\n"
        << "Container Frags                                           : " << containerFragsThisEvent << "\n"
        << "Inner Frags This Event                                    : " << innerFragsThisEvent << "\n"
        << "Unknown Frags                                             : " << unknownFragsThisEvent << "\n"
        << "Unknown Container Frags                                   : " << unknownContainersThisEvent << "\n"
        << "Non Container Frags                                       : " << nonContainersThisEvent << "\n"
        << "\n"
        << "Raw Frags With Invalid Headers (HPGe)                     : " << HPGeEventMetrics.rawFragsWithInvalidHeaders << "\n"
        << "Raw Frags With Invalid Anchors (HPGe)                     : " << HPGeEventMetrics.rawFragsWithInvalidAnchors << "\n"
        << "ZS Frags With Length Mismatch (HPGe)                      : " << HPGeEventMetrics.zsLengthMismatch << "\n"
        << "Zs Frags WIth Region Mismatch (HPGe)                      : " << HPGeEventMetrics.zsRegionMismatch << "\n"
        << "ZS Frags Skipped Due To Raw Flags (HPGe)                  : " << HPGeEventMetrics.zsFragsSkippedDueToRawFlag << "\n"
        << "ZS Frags Skipped Due To Invalid Raw Header (HPGe)         : " << HPGeEventMetrics.zsFragsSkippedDueToInvalidRawHeader << "\n"
        << "ZS Frags Skipped Due To No Preceding Raw Header (HPGe)    : " << HPGeEventMetrics.zsFragsSkippedDueToNoPrecedingRawHeader << "\n"
        << "PH Frags Skipped Due To Raw Flags (HPGe)                  : " << HPGeEventMetrics.phFragsSkippedDueToRawFlag << "\n"
        << "PH Frags Skipped Due To Invalid Raw Header (HPGe)         : " << HPGeEventMetrics.phFragsSkippedDueToInvalidRawHeader << "\n"
        << "PH Frags Skipped Due To No Preceding Raw Header (HPGe)    : " << HPGeEventMetrics.phFragsSkippedDueToNoPrecedingRawHeader << "\n"
        << "Raw/ZS/PH Sets Skipped Due To Raw Flags (HPGe)            : " << HPGeEventMetrics.setsSkippedDueToRawFlag << "\n"
        << "Raw/ZS/PH Sets Skipped Due To Invalid Raw Headers (HPGe)  : " << HPGeEventMetrics.setsSkippedDueToInvalidHeaders << "\n"
        << "\n"
        << "Raw Frags With Invalid Headers (LaBr)                     : " << LaBrEventMetrics.rawFragsWithInvalidHeaders << "\n"
        << "Raw Frags With Invalid Anchors (LaBr)                     : " << LaBrEventMetrics.rawFragsWithInvalidAnchors << "\n"
        << "ZS Frags With Length Mismatch (LaBr)                      : " << LaBrEventMetrics.zsLengthMismatch << "\n"
        << "ZS Frags With Region Mismatch (LaBr)                      : " << LaBrEventMetrics.zsRegionMismatch << "\n"
        << "ZS Frags Skipped Due To Raw Flags (LaBr)                  : " << LaBrEventMetrics.zsFragsSkippedDueToRawFlag << "\n"
        << "ZS Frags Skipped Due To Invalid Raw Header (LaBr)         : " << LaBrEventMetrics.zsFragsSkippedDueToInvalidRawHeader << "\n"
        << "ZS Frags Skipped Due To No Preceding Raw Header (LaBr)    : " << LaBrEventMetrics.zsFragsSkippedDueToNoPrecedingRawHeader << "\n"
        << "PH Frags Skipped Due To Raw Flags (LaBr)                  : " << LaBrEventMetrics.phFragsSkippedDueToRawFlag << "\n"
        << "PH Frags Skipped Due To Invalid Raw Header (LaBr)         : " << LaBrEventMetrics.phFragsSkippedDueToInvalidRawHeader << "\n"
        << "PH Frags Skipped Due To No Preceding Raw Header (LaBr)    : " << LaBrEventMetrics.phFragsSkippedDueToNoPrecedingRawHeader << "\n"
        << "Raw/ZS/PH Sets Skipped Due To Raw Flags (LaBr)            : " << LaBrEventMetrics.setsSkippedDueToRawFlag << "\n"
        << "Raw/ZS/PH Sets Skipped Due To Invalid Raw Headers (LaBr)  : " << LaBrEventMetrics.setsSkippedDueToInvalidHeaders << "\n"
        << "=================================\n";

    }
    // Frag Summary
    if (_saveSTMFragSummary) {
        event.put(std::move(stmFragSummaryHPGe), "stmFragSummaryHPGe");
        event.put(std::move(stmFragSummaryLaBr), "stmFragSummaryLaBr");
    }
    // HPGe
    if (_saveRawWaveformsWithHeaderHPGe) { event.put(std::move(rawWaveformDigisWithHeaderHPGe), "rawWithHeaderHPGe"); }
    if (_saveRawWaveformsHPGe) { event.put(std::move(rawWaveformDigisHPGe), "rawHPGe"); }
    if (_saveZSWaveformsHPGe) { event.put(std::move(zsWaveformDigisHPGe), "zsHPGe"); }
    //default to always save PH Digis
    event.put(std::move(phDigisHPGe), "phHPGe");
    // LaBr
    if (_saveRawWaveformsWithHeaderLaBr) { event.put(std::move(rawWaveformDigisWithHeaderLaBr), "rawWithHeaderLaBr"); }
    if (_saveRawWaveformsLaBr) { event.put(std::move(rawWaveformDigisLaBr), "rawLaBr"); }
    if (_saveZSWaveformsLaBr) { event.put(std::move(zsWaveformDigisLaBr), "zsLaBr"); }
    //default to always save PH Digis
    event.put(std::move(phDigisLaBr), "phLaBr");

}// End of produce()

void STMDigisFromFragments::endJob() {
    if (_verbosityLevel > 0 ){
        // Print job Summary
        mf::LogInfo("STMDigisFromFragments")
        << "\n========== STM JOB SUMMARY - (Unpacking Module) ==========\n"
        << "\n--- Module Configuration For Job ---\n"
        << "Raw Waveforms with Header HPGe   : " << (_saveRawWaveformsWithHeaderHPGe ? "Yes" : "No") << "\n"
        << "Raw Waveforms HPGe               : " << (_saveRawWaveformsHPGe ? "Yes" : "No") << "\n"
        << "ZS Waveforms HPGe                : " << (_saveZSWaveformsHPGe ? "Yes" : "No") << "\n"
        << "PH Digis HPGe                    : Yes\n"
        << "Raw Waveforms with Header LaBr   : " << (_saveRawWaveformsWithHeaderLaBr ? "Yes" : "No") << "\n"
        << "Raw Waveforms LaBr               : " << (_saveRawWaveformsLaBr ? "Yes" : "No") << "\n"
        << "ZS Waveforms LaBr                : " << (_saveZSWaveformsLaBr ? "Yes" : "No") << "\n"
        << "PH Digis LaBr                    : Yes\n"

        // Container and Inner Fragment Summary
        << "\n--- Container and Inner Fragment Job Summary ---\n"
        << "Total Art Events Processed                      : " << _totalEvents << "\n"
        << "Total Art Events w/ HPGe & LaBr                 : " << _totalEventsWithBothDetectors << "\n"
        << "Total Art Events w/ HPGe Only                   : " << _totalEventsWithOnlyHPGe << "\n"
        << "Total Art Events w/ LaBr Only                   : " << _totalEventsWithOnlyLaBr << "\n"
        << "Total Art Events w/ Neither HPGe nor LaBr       : " << _totalEventsWithNeitherDetector << "\n"

        << "Total HPGe Containers                           : " << _totalContainersHPGe << "\n"
        << "Total LaBr Containers                           : " << _totalContainersLaBr << "\n"
        << "Total Container Processed                       : " << _totalContainers << "\n"

        << "Total Inner Fragments Processed                 : " << _totalInnerFrags << "\n"
        << "Total Unreadable Inner Fragments                : " << _totalUnreadInnerFrags << "\n"

        << "Total Unknown Container Fragments               : " << _totalUnknownContainers << "\n"
        << "Total Non Container Fragments                   : " << _totalNonContainers << "\n"

        // Data Type Summary - pre-filtering
        << "\n--- Data Types Read ---\n"
        << "Total Raw Frags Seen                            : " << _totalRawFragsSeen << "\n"
        << "Total ZS  Frags Seen                            : " << _totalZSFragsSeen << "\n"
        << "Total PH  Frags Seen                            : " << _totalPHFragsSeen << "\n"
        << "Total Raw Prescaled Frags                       : " << _totalRawFragsPrescaled << "\n"
        << "Total ZS Prescaled Frags                        : " << _totalZSFragsPrescaled << "\n"
        << "Total Raw Frags Flagged Bad Only                : " << _totalRawFragsFlaggedBadOnly << "\n"
        << "Total Raw Frags Flagged Missing Only            : " << _totalRawFragsFlaggedMissingOnly << "\n"
        << "Total Raw Frags Flagged Bad And Missing         : " << _totalRawFragsFlaggedBadAndMissing << "\n"
        << "\n"
        << "Total RAW Frags seen (HPGe)                     : " << _totalRawFragsSeenHPGe << "\n"
        << "Total ZS  Frags seen (HPGe)                     : " << _totalZSFragsSeenHPGe << "\n"
        << "Total PH  Frags seen (HPGe)                     : " << _totalPHFragsSeenHPGe << "\n"
        << "\n"
        << "Total RAW Frags seen (LaBr)                     : " << _totalRawFragsSeenLaBr << "\n"
        << "Total ZS  Frags seen (LaBr)                     : " << _totalZSFragsSeenLaBr << "\n"
        << "Total PH  Frags seen (LaBr)                     : " << _totalPHFragsSeenLaBr << "\n"

        // Data Type Summary - post-filtering (General)
        << "\n--- Data Types Read Classified ---\n"
        << "Total Good  Raw Frags                           : " << _totalGoodRawFrags << "\n"
        << "Total Good  ZS  Frags                           : " << _totalGoodZSFrags << "\n"
        << "Total Good  PH  Frags                           : " << _totalGoodPHFrags << "\n"
        << "\n"
        << "Total Zero  Raw Frags                           : " << _totalZeroRawFrags << "\n"
        << "Total Zero  ZS  Frags                           : " << _totalZeroZSFrags << "\n"
        << "Total Zero  PH  Frags                           : " << _totalZeroPHFrags << "\n"
        << "\n"
        << "Total Empty Raw Frags                           : " << _totalEmptyRawFrags << "\n"
        << "Total Empty ZS  Frags                           : " << _totalEmptyZSFrags << "\n"
        << "Total Empty PH  Frags                           : " << _totalEmptyPHFrags << "\n"
        << "\n"
        << "Total Non Hits PH Frags                         : " << _totalNoHitsPHFrags << "\n"

        // Data Type Summary - post-filtering (HPGe)
        << "\n--- Data Types Read For HPGe ---\n"
        << "Total Good  Raw Frags (HPGe)                    : " << _totalGoodRawFragsHPGe << "\n"
        << "Total Good  ZS  Frags (HPGe)                    : " << _totalGoodZSFragsHPGe << "\n"
        << "Total Good  PH  Frags (HPGe)                    : " << _totalGoodPHFragsHPGe << "\n"
        << "\n"
        << "Total Zero  Raw Frags (HPGe)                    : " << _totalZeroRawFragsHPGe << "\n"
        << "Total Zero  ZS  Frags (HPGe)                    : " << _totalZeroZSFragsHPGe << "\n"
        << "Total Zero  PH  Frags (HPGe)                    : " << _totalZeroPHFragsHPGe << "\n"
        << "\n"
        << "Total Empty Raw Frags (HPGe)                    : " << _totalEmptyRawFragsHPGe << "\n"
        << "Total Empty ZS  Frags (HPGe)                    : " << _totalEmptyZSFragsHPGe << "\n"
        << "Total Empty PH  Frags (HPGe)                    : " << _totalEmptyPHFragsHPGe << "\n"
        << "\n"
        << "Total Non Hits PH Frgas (HPGe)                  : " << _totalNoHitsPHFragsHPGe << "\n"
        << "\n"
        << "Total Raw Prescaled Frags (HPGe)                : " << _totalRawFragsPrescaledHPGe << "\n"
        // Raw Prescaled Frags should be classified as their own data for now
        << "Total ZS Prescaled Frags (HPGe)                 : " << _totalZSFragsPrescaledHPGe << "\n"
        << "Total Raw Frags Flagged Bad Only (HPGe)         : " << _totalRawFragsFlaggedBadOnlyHPGe << "\n"
        << "Total Raw Frags Flagged Missing Only (HPGe)     : " << _totalRawFragsFlaggedMissingOnlyHPGe << "\n"
        << "Total Raw Frags Flagged Bad & Missing (HPGe)    : " << _totalRawFragsFlaggedBadAndMissingHPGe << "\n"

        // Data Type Summary - post-filtering (LaBr)
        << "\n--- Data Types Read For LaBr ---\n"
        << "Total Good  Raw Frags (LaBr)                    : " << _totalGoodRawFragsLaBr << "\n"
        << "Total Good  ZS  Frags (LaBr)                    : " << _totalGoodZSFragsLaBr << "\n"
        << "Total Good  PH  Frags (LaBr)                    : " << _totalGoodPHFragsLaBr << "\n"
        << "\n"
        << "Total Zero  Raw Frags (LaBr)                    : " << _totalZeroRawFragsLaBr << "\n"
        << "Total Zero  ZS  Frags (LaBr)                    : " << _totalZeroZSFragsLaBr << "\n"
        << "Total Zero  PH  Frags (LaBr)                    : " << _totalZeroPHFragsLaBr << "\n"
        << "\n"
        << "Total Empty Raw Frags (LaBr)                    : " << _totalEmptyRawFragsLaBr << "\n"
        << "Total Empty ZS  Frags (LaBr)                    : " << _totalEmptyZSFragsLaBr << "\n"
        << "Total Empty PH  Frags (LaBr)                    : " << _totalEmptyPHFragsLaBr << "\n"
        << "\n"
        << "Total Non Hits PH Frgas (LaBr)                  : " << _totalNoHitsPHFragsLaBr << "\n"
        << "\n"
        << "Total Raw Prescaled Frags (LaBr)                : " << _totalRawFragsPrescaledLaBr << "\n"
        << "Total ZS Prescaled Frags (LaBr)                 : " << _totalZSFragsPrescaledLaBr << "\n"
        << "Total Raw Frags Flagged Bad Only (LaBr)         : " << _totalRawFragsFlaggedBadOnlyLaBr << "\n"
        << "Total Raw Frags Flagged Missing Only (LaBr)     : " << _totalRawFragsFlaggedMissingOnlyLaBr << "\n"
        << "Total Raw Frags Flagged Bad & Missing (LaBr)    : " << _totalRawFragsFlaggedBadAndMissingLaBr << "\n"

        << "\n=== Extra Filters ===\n"
        << "Invalid Raw Headers                             : " << _totalRawFragsWithInvalidHeaders << "\n"
        << "Invalid Raw Headers (HPGe)                      : " << _totalRawFragsWithInvalidHeadersHPGe << "\n"
        << "Invalid Raw Headers (LaBr)                      : " << _totalRawFragsWithInvalidHeadersLaBr << "\n"
        << "Invalid Raw Anchors                             : " << _totalRawFragsWithInvalidAnchors << "\n"
        << "Invalid Raw Anchors (HPGe)                      : " << _totalRawFragsWithInvalidAnchorsHPGe << "\n"
        << "Invalid Raw Anchors (LaBr)                      : " << _totalRawFragsWithInvalidAnchorsLaBr << "\n"

        << "ZS Skipped Due To Raw Flags (HPGe)              : " << _totalZSFragsSkippedDueToRawFlagHPGe << "\n"
        << "ZS Skipped Due To Raw Flags (LaBr)              : " << _totalZSFragsSkippedDueToRawFlagLaBr << "\n"

        << "PH Skipped Due To Raw Flags (HPGe)              : " << _totalPHFragsSkippedDueToRawFlagHPGe << "\n"
        << "PH Skipped Due To Raw Flags (LaBr)              : " << _totalPHFragsSkippedDueToRawFlagLaBr << "\n"

        << "ZS Length Mismatches Encountered                : " << _totalZSLengthMismatch << "\n"
        << "ZS Length Mismatches Encountered (HPGe)         : " << _totalZSLengthMismatchHPGe << "\n"
        << "ZS Length Mismtaches Enocuntered (LaBr)         : " << _totalZSLengthMismatchLaBr << "\n"

        << "ZS Region Mismatches Encountered                : " << _totalZSRegionMismatch << "\n"
        << "ZS Region Mismatches Encountered (HPGe)         : " << _totalZSRegionMismatchHPGe << "\n"
        << "ZS Region Mismatches Encountered (LaBr)         : " << _totalZSRegionMismatchLaBr << "\n"

        << "PH Count Mismatches Encountered                 : " << _totalPHCountMismatch << "\n"
        << "PH Count Mismatches Encountered (HPGe)          : " << _totalPHCountMismatchHPGe << "\n"
        << "PH Count Mismatches Encountered (LaBr)          : " << _totalPHCountMismatchLaBr << "\n"

        << "\n===========================================================\n";

    }
}

DEFINE_ART_MODULE(STMDigisFromFragments)
