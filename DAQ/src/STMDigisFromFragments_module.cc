// ====================================================================
//
// STMDigisFromFragments: create all types of STMDigis from STMFragments
// Note that for the STM, a single art event contains many EWTs
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/STMFragmentSummary.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "Offline/RecoDataProducts/inc/STMEventHeader.hh"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>
#include "canvas/Persistency/Common/Ptr.h"

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdbool.h>
#include <map>

// This is version 2
// Meant to clean up and make the code more maintainable and readable

namespace art
{
    class STMDigisFromFragments;
}
using art::STMDigisFromFragments;
// can make easier but lets be explicit for now
//using STMWaveformDigiMap = std::map<mu2e::STMEventHeader, std::vector<mu2e::STMWaveformDigi>>;

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

    // Used to set and verify the state of the raw parent fragment for zs fragments, event based
    struct RawParentState {
        art::Ptr<mu2e::STMWaveformDigi> ptr;
        art::ProductID productID;
        size_t index{0};
        bool available{false};
    };

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

        bool rawHeaderIsValid{false};
        bool skipCurrentSet{false};

        // save the following now for mapping
        uint64_t eventWindowTag{0};
        uint8_t eventMode{0}; // Not sure what this is right now
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
        size_t phCountMismatch{0};

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
    /*
    // HPGe
    if (_saveRawWaveformsWithHeaderHPGe){ produces<mu2e::STMWaveformDigiCollection>("rawWithHeaderHPGe"); }
    if (_saveRawWaveformsHPGe){ produces<mu2e::STMWaveformDigiCollection>("rawHPGe"); }
    if (_saveZSWaveformsHPGe){ produces<mu2e::STMWaveformDigiCollection>("zsHPGe"); }
    produces<mu2e::STMPHDigiCollection>("phHPGe");
    // LaBr
    if (_saveRawWaveformsWithHeaderLaBr){ produces<mu2e::STMWaveformDigiCollection>("rawWithHeaderLaBr"); }
    if (_saveRawWaveformsLaBr){ produces<mu2e::STMWaveformDigiCollection>("rawLaBr"); }
    if (_saveZSWaveformsLaBr){ produces<mu2e::STMWaveformDigiCollection>("zsLaBr"); }
    produces<mu2e::STMPHDigiCollection>("phLaBr");
    */

    //change products to be maps of event headers
    // HPGe
    if (_saveRawWaveformsWithHeaderHPGe){produces<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>>("rawWithHeaderHPGe"); }
    if (_saveRawWaveformsHPGe){produces<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>>("rawHPGe"); }
    if (_saveZSWaveformsHPGe){produces<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>>("zsHPGe"); }
    produces<std::map<mu2e::STMEventHeader, mu2e::STMPHDigiCollection>>("phHPGe");
    // LaBr
    if (_saveRawWaveformsWithHeaderLaBr){ produces<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>>("rawWithHeaderLaBr"); }
    if (_saveRawWaveformsLaBr){ produces<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>>("rawLaBr"); }
    if (_saveZSWaveformsLaBr){ produces<std::map<mu2e::STMEventHeader, mu2e::STMWaveformDigiCollection>>("zsLaBr"); }
    produces<std::map<mu2e::STMEventHeader, mu2e::STMPHDigiCollection>>("phLaBr");
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

    ++_totalEvents; // Increments total event counter

    // Each RawParentState keeps track of the last Raw Fragment
    RawParentState rawParentHPGe;
    RawParentState rawParentLaBr;
    RawHeaderState rawHeaderHPGe;
    RawHeaderState rawHeaderLaBr;
    DetectorSpecificEventMetrics LaBrEventMetrics;
    DetectorSpecificEventMetrics HPGeEventMetrics;

    // Set product ID

    /*
    if (_saveRawWaveformsHPGe){
        rawParentHPGe.productID = event.getProductID<mu2e::STMWaveformDigiCollection>("rawHPGe");
    }
    if (_saveRawWaveformsLaBr){
        rawParentLaBr.productID = event.getProductID<mu2e::STMWaveformDigiCollection>("rawLaBr");
    }
    */

    // Frag Summaries
    std::unique_ptr<mu2e::STMFragmentSummaryCollection> stmFragSummaryHPGe(new mu2e::STMFragmentSummaryCollection);
    std::unique_ptr<mu2e::STMFragmentSummaryCollection> stmFragSummaryLaBr(new mu2e::STMFragmentSummaryCollection);

    /*
    // HPGe
    std::unique_ptr<mu2e::STMWaveformDigiCollection> rawWaveformDigisWithHeaderHPGe(new mu2e::STMWaveformDigiCollection);
    std::unique_ptr<mu2e::STMWaveformDigiCollection> rawWaveformDigisHPGe(new mu2e::STMWaveformDigiCollection);
    std::unique_ptr<mu2e::STMWaveformDigiCollection> zsWaveformDigisHPGe(new mu2e::STMWaveformDigiCollection);
    std::unique_ptr<mu2e::STMPHDigiCollection> phDigisHPGe(new mu2e::STMPHDigiCollection);

    // LaBr
    std::unique_ptr<mu2e::STMWaveformDigiCollection> rawWaveformDigisWithHeaderLaBr(new mu2e::STMWaveformDigiCollection);
    std::unique_ptr<mu2e::STMWaveformDigiCollection> rawWaveformDigisLaBr(new mu2e::STMWaveformDigiCollection);
    std::unique_ptr<mu2e::STMWaveformDigiCollection> zsWaveformDigisLaBr(new mu2e::STMWaveformDigiCollection);
    std::unique_ptr<mu2e::STMPHDigiCollection> phDigisLaBr(new mu2e::STMPHDigiCollection);
    */

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
        if  (_verbosityLevel > 2) { std::cout << "Processing outer fragment with ID: " << outerFragID << "\n"; }

        // Check if this is a container fragment
        if (frag.type() == artdaq::Fragment::ContainerFragmentType) {

            mu2e::STMFragment container_frag(frag);
            artdaq::ContainerFragment cont_frag(frag);
            size_t blocks = cont_frag.block_count();

            if (container_frag.isHPGeContainer()) {
                if  (_verbosityLevel > 2) { std::cout << "Processing HPGe container fragment with ID: " << outerFragID << "\n"; }
                ++_totalContainersHPGe; eventHasHPGe = true; _totalInnerFragsHPGe += blocks;
                innerFragsHPGeThisEvent += blocks;
                ++containerFragsHPGeThisEvent;
            } else if (container_frag.isLaBrContainer()) {
                if  (_verbosityLevel > 2) { std::cout << "Processing LaBr container fragment with ID: " << outerFragID << "\n"; }
                ++_totalContainersLaBr; eventHasLaBr = true; _totalInnerFragsLaBr += blocks;
                innerFragsLaBrThisEvent += blocks;
                ++containerFragsLaBrThisEvent;
            } else {
                if (_verbosityLevel > 0) {
                    std::cout << "Encountered an unknown STM Container Fragment \n"
                    << "Frag ID : "  << outerFragID << "\n"
                    << "Event   : " << _totalEvents << "\n";
                }
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
                        if (_verbosityLevel > 1) {
                            std::cout << "Encountered raw fragment that is neither HPGe nor LaBr\n";
                        }
                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    isHPGe ? ++_totalRawFragsSeenHPGe : ++_totalRawFragsSeenLaBr;

                    auto& headerState = isHPGe ? rawHeaderHPGe : rawHeaderLaBr;
                    auto& parentState = isHPGe ? rawParentHPGe : rawParentLaBr;
                    auto& eventMetrics  = isHPGe ? HPGeEventMetrics : LaBrEventMetrics;

                    // reset current raw header and parent state bool
                    headerState = RawHeaderState{};
                    parentState.ptr = {};
                    parentState.available = false;
                    parentState.index = 0;
                    ++eventMetrics.raw.seen;
                    // reset eventHeader related variables
                    headerState.eventWindowTag = 0;
                    headerState.eventMode = 0;
                    headerState.adcClock = 0;
                    headerState.dtcClock = 0;

                    // check header+adcs is not less than expected header length (22)
                    if(stm_frag.dataWords() < stm::RawHeader::WORDS){
                        if (_verbosityLevel > 1) {
                            std::cout << "Raw fragment has fewer data words than expected header length\n";
                        }
                        ++_totalRawFragsWithInvalidHeaders;
                        ++_totalUnreadInnerFrags;
                        ++eventMetrics.rawFragsWithInvalidHeaders;
                        isHPGe ? ++_totalRawFragsWithInvalidHeadersHPGe : ++_totalRawFragsWithInvalidHeadersLaBr;
                        ++eventMetrics.setsSkippedDueToInvalidHeaders;

                        headerState.skipCurrentSet = true;
                        parentState.ptr = {};
                        parentState.available = false;
                        continue;
                    }

                    // Confirm Raw Header has Valid Anchors
                    // For now just use as a diagnostic tool
                    if (!stm_frag.hasValidAnchors()) {
                        ++eventMetrics.rawFragsWithInvalidAnchors;
                        ++_totalRawFragsWithInvalidAnchors;
                        isHPGe ? ++_totalRawFragsWithInvalidAnchorsHPGe : ++_totalRawFragsWithInvalidAnchorsLaBr;
                    }
                    // Not needed for the moment
                    // headerState.rawHeaderIsValid = stm_frag.hasValidAnchors();

                    // Check if bad or missing or both
                    if (stm_frag.badData() && stm_frag.missing()) {
                        ++eventMetrics.rawFragsFlaggedBadAndMissing;
                        ++_totalRawFragsFlaggedBadAndMissing;
                        isHPGe ? ++_totalRawFragsFlaggedBadAndMissingHPGe : ++_totalRawFragsFlaggedBadAndMissingLaBr;
                    } else if (stm_frag.badData()) {
                        ++eventMetrics.rawFragsFlaggedBadOnly;
                        ++_totalRawFragsFlaggedBadOnly;
                        isHPGe ? ++_totalRawFragsFlaggedBadOnlyHPGe : ++_totalRawFragsFlaggedBadOnlyLaBr;
                    } else if (stm_frag.missing()) {
                        ++eventMetrics.rawFragsFlaggedMissingOnly;
                        ++_totalRawFragsFlaggedMissingOnly;
                        isHPGe ? ++_totalRawFragsFlaggedMissingOnlyHPGe : ++_totalRawFragsFlaggedMissingOnlyLaBr;
                    }

                    bool const badOrMissing = stm_frag.badData() || stm_frag.missing();
                    if (badOrMissing) {
                        headerState.skipCurrentSet = true;
                        parentState.ptr = {};
                        parentState.available = false;
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
                        ++eventMetrics.raw.prescaled;
                        ++_totalRawFragsPrescaled;
                        isHPGe ? ++_totalRawFragsPrescaledHPGe : ++_totalRawFragsPrescaledLaBr;

                        parentState.ptr = {};
                        parentState.available = false;
                        continue;
                    }

                    auto payloadPtr = stm_frag.payloadBegin();
                    auto payloadWords = stm_frag.payloadWords();

                    // 22 Header Words + adcs < = dataWords
                    // This check comes after ensuring raw is not prescaled
                    if (stm::RawHeader::WORDS + payloadWords > stm_frag.dataWords()) {
                        // this is an invalid-header case
                        ++_totalRawFragsWithInvalidHeaders;
                        isHPGe ? ++_totalRawFragsWithInvalidHeadersHPGe : ++_totalRawFragsWithInvalidHeadersLaBr;
                        ++eventMetrics.rawFragsWithInvalidHeaders;
                        ++eventMetrics.setsSkippedDueToInvalidHeaders;
                        ++_totalUnreadInnerFrags;
                        headerState.skipCurrentSet = true;
                        parentState.available = false;
                        continue;
                    }

                    // Ideally by here we know the raw frag is good and its payload can be checked
                    // Check if raw frag is empty
                    bool allZeros = true;
                    if (payloadWords == 0) {
                        if (_verbosityLevel > 2){
                            std::cout << "\nFound an empty raw fragment at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "--Raw Frag\n";
                        }
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
                        if (_verbosityLevel > 2){
                            std::cout << "\nFound an all-zero raw fragment at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "--Raw Frag\n";
                        }
                        //counter increment for all-zero raw fragments
                        ++eventMetrics.raw.zero;
                        ++_totalZeroRawFrags;
                        isHPGe? ++_totalZeroRawFragsHPGe: ++_totalZeroRawFragsLaBr;
                        continue;
                    }
                    // At this point the raw fragment has payload words and is non zero filled
                    // We continue processing this raw fragment
                    if (_verbosityLevel > 2) {
                        std::cout << "\nFound a good raw fragment at Event : " << _totalEvents << "\n"
                        << "Fragment Index : " << i << "\n"
                        << "--Raw Frag\n";
                    }

                    //Print first few payload words for inspection
                    if (_verbosityLevel > 3) {
                        std::cout << "First few payload words(adcs) for inspection : " ;
                        for (size_t w = 0; w < std::min(payloadWords, static_cast<size_t>(10)); ++w){
                            std::cout << payloadPtr[w] << " ,";
                        }
                        std::cout << "\n";
                    }
                    //Raw Header inspection
                    if (_verbosityLevel > 4) {
                        std::cout << "\nRaw Header for inspection \n"
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
                    // Save Raw Waveform and Save parent Ptr - HPGe
                    if (_saveRawWaveformsHPGe && isHPGe){
                        // set stm_waveform
                        stm_waveform.set_data(payloadWords, payloadPtr);

                        // Give map a key and value, key is stmHeader,
                        // *rawwaveformDigisHPGe is the map
                        // stm_event_header is the collection stored under this header
                        // .push_back(stm_waveform) adds waveform to the collection
                        (*rawWaveformDigisHPGe)[stm_event_header].push_back(stm_waveform);

                        //Save raw parent index
                        //parentState.index = rawWaveformDigisHPGe->size() - 1;
                        //parentState.ptr = art::Ptr<mu2e::STMWaveformDigi>(parentState.productID, parentState.index,
                        //    event.productGetter(parentState.productID));
                        //parentState.available = true;
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
                    // Save Raw Waveform and Save parent Ptr - LaBr
                    if (_saveRawWaveformsLaBr && isLaBr){
                        // set stm_waveform
                        stm_waveform.set_data(payloadWords, payloadPtr);
                        // set map here
                        (*rawWaveformDigisLaBr)[stm_event_header].push_back(stm_waveform);

                        //Save raw parent index
                        //parentState.index = rawWaveformDigisLaBr->size() - 1;
                        //parentState.ptr = art::Ptr<mu2e::STMWaveformDigi>(parentState.productID, parentState.index,
                        //event.productGetter(parentState.productID));
                        //parentState.available = true;
                    }
                    if (_verbosityLevel > 2) {
                      std::cout << "Raw with frag index" << i
                                << ", detector = " << (isHPGe? "HPGe": "LaBr") << std::boolalpha
                                << ", prescaled = " << headerState.rawPrescaled
                                << ", saved parent = " << parentState.available;

                        if (parentState.available) {
                          std::cout << ", raw index saved = " << parentState.index;
                            }
                      std::cout << "\n";
                    }
                }// End of Raw fragment check
                else if (stm_frag.isZS()){
                    ++_totalZSFragsSeen;
                    // Determine which detector this ZS fragment belongs to
                    bool const isHPGe = stm_frag.isHPGe();
                    bool const isLaBr = stm_frag.isLaBr();

                    // Double check that the fragment belongs to one of the expected detectors
                    if (!isHPGe && !isLaBr){
                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    auto& headerState = isHPGe ? rawHeaderHPGe : rawHeaderLaBr;
                    auto& parentState = isHPGe ? rawParentHPGe : rawParentLaBr;
                    auto& eventMetrics  = isHPGe ? HPGeEventMetrics : LaBrEventMetrics;

                    isHPGe ? ++_totalZSFragsSeenHPGe : ++_totalZSFragsSeenLaBr;
                    ++eventMetrics.zs.seen;

                    // Check
                    if (_verbosityLevel > 2) {
                      std::cout << "ZS frag index = " << i
                                << ", detector = " << (isHPGe ? "HPGe" : "LaBr") << std::boolalpha
                                << ", parent available = " << parentState.available;
                      if (parentState.available) {
                        std::cout << ", raw index saved = " << parentState.index;
                      }
                      std::cout << "\n";
                    }

                    // Extract zs variables from Raw Header
                    bool skipCurrentSet = headerState.skipCurrentSet;
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

                    // Decide Here to skip based on previous raw fragment information
                    if (skipCurrentSet){
                        isHPGe ? ++_totalZSFragsSkippedDueToRawFlagHPGe : ++_totalZSFragsSkippedDueToRawFlagLaBr;
                        ++eventMetrics.zsFragsSkippedDueToRawFlag;
                        continue;
                    }
                    // Decide Here to skip based on zs prescale information
                    if (zsPrescaled) {
                        isHPGe ? ++_totalZSFragsPrescaledHPGe : ++_totalZSFragsPrescaledLaBr;
                        ++_totalZSFragsPrescaled;
                        ++eventMetrics.zs.prescaled;
                        continue;
                    }

                    auto payloadPtr = stm_frag.payloadBegin();
                    auto payloadWords = stm_frag.payloadWords();

                    // Check if Empty
                    if (payloadWords == 0) {
                        if (_verbosityLevel > 2) {
                            std::cout << "\nFound an empty zs fragment at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "--ZS Frag\n";
                        }
                        ++_totalEmptyZSFrags;
                        isHPGe ? ++_totalEmptyZSFragsHPGe : ++_totalEmptyZSFragsLaBr;
                        ++eventMetrics.zs.empty;
                        continue;
                    }

                    // At this point the zs fragmant is non-empty
                    // May contain valid data

                    // Print first 20 payload adcs for inspection
                    if (_verbosityLevel > 3) {
                        std::cout << "\nFirst few payload words for inspection : ";
                        for (size_t w = 0; w < std::min(payloadWords, static_cast<size_t>(10)); ++w) {
                            std::cout << payloadPtr[w] << ",";
                        }
                        std::cout << "\n";
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
                        std::cout << "\nData Words for inspection : " << "\n"
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
                            std::cout << "\nZS Segment Check :" << "\n"
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
                        std::cout << "ZS Summary : "
                        << " , last ZS Index = " << lastZSIndexRecorded
                        << " , last ZS Length = " << lastZSLengthRecorded
                        << " , ZS Total Length = " << zsTotalLengthCalculated
                        << " , Region Counter = " << regionCounter << "\n"
                        << ", Frag Index = " << i << "\n";
                    }

                    // In case ZS has some weird behavior
                    if (malformedZS) {
                        if (_verbosityLevel > 2){
                            std::cout << "\nMalformed ZS detected at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n";
                        }
                        ++eventMetrics.zs.unread;
                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    // Add exceptions
                    if (zsInfoWasExtracted) {
                        if (zsLength != zsTotalLengthCalculated) {
                            throw cet::exception("STM_UNPACKING")
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

                        }
                        if ( zsRegions != regionCounter) {
                            throw cet::exception("STM_UNPACKING")
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
                        }
                    }

                    if (sawADC && allADCSZero) {
                        if (_verbosityLevel > 2) {
                            std::cout << "\nFound a zero-filled zs fragment at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "--ZS Frag\n";
                        }
                        ++eventMetrics.zs.zero;
                        ++_totalZeroZSFrags;
                        isHPGe ? ++_totalZeroZSFragsHPGe : ++_totalZeroZSFragsLaBr;
                        continue;
                    }

                    // At this point we know the fragment was good
                    if (_verbosityLevel > 2) {
                        std::cout << "\nFound a good ZS fragment at Event : " << _totalEvents << "\n"
                        << "Frag Index : " << i << "\n"
                        << "--ZS Frag\n";
                    }
                    for (auto& region : regions) {
                        // set waveform digi with offset and adcs
                        mu2e::STMWaveformDigi zsDigi(region.offset, region.adcs);
                        //if(parentState.available){
                            // Set parent pointer if available
                          //  zsDigi.setParent(parentState.ptr);
                        //}

                        // Emplacing
                        if (isHPGe && _saveZSWaveformsHPGe) {
                          // set map for HPGe
                            (*zsWaveformDigisHPGe)[stm_event_header].push_back(zsDigi);
                          //zsWaveformDigisHPGe->emplace_back(zsDigi);
                            if (_verbosityLevel > 2 && zsDigi.hasParent()) {
                                std::cout << "\nZS HPGe Parent Raw index  : " << zsDigi.parent().key()
                                << ", offset : " << zsDigi.trigTimeOffset() << "\n";
                            }

                        } else if (isLaBr && _saveZSWaveformsLaBr) {
                            // set map for LaBr
                           (*zsWaveformDigisLaBr)[stm_event_header].push_back(zsDigi);
                            //zsWaveformDigisLaBr->emplace_back(zsDigi);
                            if (_verbosityLevel > 2 && zsDigi.hasParent()) {
                                std::cout << "\nZS LaBr Parent Raw index  : " << zsDigi.parent().key()
                                << ", offset : " << zsDigi.trigTimeOffset() << "\n";
                            }
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
                        ++_totalUnreadInnerFrags;
                        continue;
                    }
                    auto& headerState = isHPGe ? rawHeaderHPGe : rawHeaderLaBr;
                    auto& eventMetrics  = isHPGe ? HPGeEventMetrics : LaBrEventMetrics;

                    isHPGe ? ++_totalPHFragsSeenHPGe : ++_totalPHFragsSeenLaBr;
                    ++eventMetrics.ph.seen;

                    // Extract ph varibales from Raw Header
                    bool skipCurrentSet = headerState.skipCurrentSet;
                    bool extractedPHInfo = headerState.containsPHInfo;
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

                    // Skip if Raw Fragment was Bad or Missing
                    if (skipCurrentSet) {
                        ++eventMetrics.phFragsSkippedDueToRawFlag;
                        isHPGe ? ++_totalPHFragsSkippedDueToRawFlagHPGe : ++_totalPHFragsSkippedDueToRawFlagLaBr;
                        continue;
                    }

                    // Check if PH fragment is empty
                    auto payloadPtr = stm_frag.payloadBegin();
                    auto payloadWords = stm_frag.payloadWords();
                    bool allPHAreZeros= true;

                    if (payloadWords == 0) {
                        if (_verbosityLevel > 2) {
                            std::cout << "\nFound an empty PH fragment at Event : " <<  _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "--PH Frag\n";
                        }
                        //increment
                        ++eventMetrics.ph.empty;
                        ++_totalEmptyPHFrags;
                        isHPGe ? ++_totalEmptyPHFragsHPGe : ++_totalEmptyPHFragsLaBr;
                        continue;
                    }

                    // Check that the payload has pairs (time, pulse height)
                    if (payloadWords % 2 != 0) {
                        if(_verbosityLevel > 2) {
                            std::cout << "\nFound a PH Fragment with odd number of payload words at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "---PH Frag\n";
                        }
                        ++eventMetrics.ph.unread;
                        ++_totalUnreadInnerFrags;
                        continue;
                    }

                    // Check if PH Pair count matches expected count from raw header
                    size_t const nPairsInFragment = payloadWords / 2;

                    // Add a print out to compare the two if there is a mismatch
                    // Maybe add an if (nPairsInfragment != phCount)
                    // Include a counter for this -> future decide whether to throw or not
                    if (extractedPHInfo && nPairsInFragment != phCount) {
                        if (_verbosityLevel > 2) {
                            std::cout << "\nPH Count Comparison\n"
                            << "Number of PH pairs we read in fragment : " << nPairsInFragment <<"\n"
                            << "Number of PH pars from raw header      : " << phCount << "\n";
                        }
                        ++_totalPHCountMismatch;
                        isHPGe ? ++_totalPHCountMismatchHPGe : ++_totalPHCountMismatchLaBr;
                        ++eventMetrics.phCountMismatch;
                    }

                    // Here you either read the lowest pair count from raw header or what is in fragment
                    // If not you read what what is in the fragment anyway
                    size_t const nPairsToRead = extractedPHInfo ? std::min<size_t>(nPairsInFragment, phCount) : nPairsInFragment;
                    if (_verbosityLevel > 2) {
                        std::cout << "\nPH Pair Count We Will Processes : " << nPairsToRead << "\n"
                        << "Pair Count coming from Fragment   : " << (nPairsInFragment < phCount ? "Yes" : "No") << "\n"
                        << "Pair Count coming from Raw Header : " << (phCount < nPairsInFragment ? "Yes" : "No") << "\n";
                    }

                    if (nPairsToRead ==0){
                        continue;
                    }

                    // Check if zero filled
                    // Since its a (time,PH) pair we will check every second entry for the PH value
                    for (size_t k = 0; k < nPairsToRead; ++k) {
                        if (payloadPtr[k*2+1] !=0) {
                            allPHAreZeros = false;
                            break;
                        }
                    }

                    if (allPHAreZeros) {
                        if (_verbosityLevel > 2) {
                            std:: cout << "\nFound a zero-filled PH Fragment at Event : " << _totalEvents << "\n"
                            << "Frag Index : " << i << "\n"
                            << "---PH Frag\n";
                        }
                        // counter
                        ++eventMetrics.ph.zero;
                        ++_totalZeroPHFrags;
                        isHPGe ? ++_totalZeroPHFragsHPGe : ++_totalZeroPHFragsLaBr;
                        continue;
                    }

                    if (_verbosityLevel > 4) {
                        std::cout << "\nFirst Few Payload Words : " ;
                        for (size_t k = 0; k < std::min<size_t>(payloadWords,20); ++k) {
                            std::cout << payloadPtr[k] << " ,";
                        }
                        std::cout << "Frag Index : " << i << "\n";
                    }

                    for (size_t i_phPair = 0; i_phPair < nPairsToRead ; ++i_phPair){
                        size_t i_PH = 2 * i_phPair;
                        uint32_t time = static_cast<uint16_t>(payloadPtr[i_PH]);
                        int16_t const pulseHeight = payloadPtr[i_PH + 1];
                        mu2e::STMPHDigi PH_digi(time,pulseHeight);

                        // Emplace Back (Always On)
                        if (isHPGe) {
                          //phDigisHPGe->emplace_back(PH_digi);
                          (*phDigisHPGe)[stm_event_header].emplace_back(PH_digi);
                        }
                        if (isLaBr) {
                          //phDigisLaBr->emplace_back(PH_digi);
                          (*phDigisLaBr)[stm_event_header].emplace_back(PH_digi);
                        }
                    }

                    // At this point we have a good PH fragment
                    if (_verbosityLevel > 2) {
                        std::cout << "\nFound a good PH Fragment at Event : " << _totalEvents << "\n"
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
                    if (_verbosityLevel > 0 ) {
                        std::cout << "Encountered an unreadable inner fragment " << "\n"
                        << "Frag Index  : " << i << "\n"
                        << "Frag ID     : " << inner_frag->fragmentID() << "\n"
                        << "Event       : " << _totalEvents << "\n";
                    }
                }// End of else non-raw/zs/pg frag
            }
        } else {
            // fallback for non-container fragment
            if (_verbosityLevel > 0) {
                std::cout << "\nEncountered a non-container fragment " << "\n"
                << "Fragment ID : " << frag.fragmentID() << "\n"
                << "Event       : " << _totalEvents << "\n";
            }
            ++_totalNonContainers;
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
            std::cout << "\nSaving STMFragSummary\n";
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

    // Final Move
    if (_verbosityLevel > 1) {
        // Event Summary -> tells us what happened per event
        std::cout << "\n========== STM EVENT SUMMARY - (Unpacking Module) ==========\n";
        std::cout << "\n--- Module Configuration ---\n";
        std::cout << "Raw Waveforms with Header HPGe   : " << (_saveRawWaveformsWithHeaderHPGe ? "Yes" : "No") << "\n";
        std::cout << "Raw Waveforms HPGe               : " << (_saveRawWaveformsHPGe ? "Yes" : "No") << "\n";
        std::cout << "ZS Waveforms HPGe                : " << (_saveZSWaveformsHPGe ? "Yes" : "No") << "\n";
        std::cout << "PH Digis HPGe                    : Yes\n";
        std::cout << "Raw Waveforms with Header LaBr   : " << (_saveRawWaveformsWithHeaderLaBr ? "Yes" : "No") << "\n";
        std::cout << "Raw Waveforms LaBr               : " << (_saveRawWaveformsLaBr ? "Yes" : "No") << "\n";
        std::cout << "ZS Waveforms LaBr                : " << (_saveZSWaveformsLaBr ? "Yes" : "No") << "\n";
        std::cout << "PH Digis LaBr                    : Yes\n";

        std::cout << "\n--- Products Saved ---\n";
        std::cout << "Extracted Raw Waveforms With Header (HPGe)    : " << rawWaveformDigisWithHeaderHPGe->size() << "\n";
        std::cout << "Extracted Raw Waveforms (HPGe)                : " << rawWaveformDigisHPGe->size() << "\n";
        std::cout << "Extracted ZS  Waveforms (HPGe)                : " << zsWaveformDigisHPGe->size() << "\n";
        std::cout << "Extracted PH  Digis     (HPGe)                : " << phDigisHPGe->size() << "\n";
        std::cout << "\n";
        std::cout << "Extracted Raw Waveforms With Header (LaBr)    : " << rawWaveformDigisWithHeaderLaBr->size() << "\n";
        std::cout << "Extracted Raw Waveforms (LaBr)                : " << rawWaveformDigisLaBr->size() << "\n";
        std::cout << "Extracted ZS  Waveforms (LaBr)                : " << zsWaveformDigisLaBr->size() << "\n";
        std::cout << "Extracted PH  Digis     (LaBr)                : " << phDigisLaBr->size() << "\n";

        std::cout << "\n--- Filtered Products (HPGe) ---\n";
        std::cout << "Good  Raw Frags (HPGe)                        : " << HPGeEventMetrics.raw.good << "\n";
        std::cout << "Good  ZS  Frags (HPGe)                        : " << HPGeEventMetrics.zs.good << "\n";
        std::cout << "Good  PH  Frags (HPGe)                        : " << HPGeEventMetrics.ph.good << "\n";
        std::cout << "\n";
        std::cout << "Empty Raw Frags (HPGe)                        : " << HPGeEventMetrics.raw.empty << "\n";
        std::cout << "Empty ZS  Frags (HPGe)                        : " << HPGeEventMetrics.zs.empty << "\n";
        std::cout << "Empty PH  Frags (HPGe)                        : " << HPGeEventMetrics.ph.empty << "\n";
        std::cout << "\n";
        std::cout << "Zero  Raw Frags (HPGe)                        : " << HPGeEventMetrics.raw.zero << "\n";
        std::cout << "Zero  ZS  Frags (HPGe)                        : " << HPGeEventMetrics.zs.zero << "\n";
        std::cout << "Zero  PH  Frags (HPGe)                        : " << HPGeEventMetrics.ph.zero << "\n";
        std::cout << "\n";
        std::cout << "Bad   Raw Frags (HPGe)                        : " << HPGeEventMetrics.rawFragsFlaggedBadOnly << "\n";
        std::cout << "Missing Raw Frags (HPGe)                      : " << HPGeEventMetrics.rawFragsFlaggedMissingOnly << "\n";
        std::cout << "Bad and Missing Raw Frags (HPGe)              : " << HPGeEventMetrics.rawFragsFlaggedBadAndMissing << "\n";
        std::cout << "\n";

        std::cout << "\n--- Filtered Products (LaBr) ---\n";
        std::cout << "Good  Raw Frags (LaBr)                        : " << LaBrEventMetrics.raw.good << "\n";
        std::cout << "Good  ZS  Frags (LaBr)                        : " << LaBrEventMetrics.zs.good << "\n";
        std::cout << "Good  PH  Frags (LaBr)                        : " << LaBrEventMetrics.ph.good << "\n";
        std::cout << "\n";
        std::cout << "Empty Raw Frags (LaBr)                        : " << LaBrEventMetrics.raw.empty << "\n";
        std::cout << "Empty ZS  Frags (LaBr)                        : " << LaBrEventMetrics.zs.empty << "\n";
        std::cout << "Empty PH  Frags (LaBr)                        : " << LaBrEventMetrics.ph.empty << "\n";
        std::cout << "\n";
        std::cout << "Zero  Raw Frags (LaBr)                        : " << LaBrEventMetrics.raw.zero << "\n";
        std::cout << "Zero  ZS  Frags (LaBr)                        : " << LaBrEventMetrics.zs.zero << "\n";
        std::cout << "Zero  PH  Frags (LaBr)                        : " << LaBrEventMetrics.ph.zero << "\n";
        std::cout << "\n";
        std::cout << "Bad   Raw Frags (LaBr)                        : " << LaBrEventMetrics.rawFragsFlaggedBadOnly << "\n";
        std::cout << "Missing Raw Frags (LaBr)                      : " << LaBrEventMetrics.rawFragsFlaggedMissingOnly << "\n";
        std::cout << "Bad and Missing Raw Frags (LaBr)              : " << LaBrEventMetrics.rawFragsFlaggedBadAndMissing << "\n";
        // Extra Filters
        std::cout << "\n=== Extra Filters ===\n";
        std::cout << "Container Frags                                           : " << containerFragsThisEvent << "\n";
        std::cout << "Inner Frags This Event                                    : " << innerFragsThisEvent << "\n";
        std::cout << "Unknown Frags                                             : " << unknownFragsThisEvent << "\n";
        std::cout << "Raw Frags With Invalid Headers (HPGe)                     : " << HPGeEventMetrics.rawFragsWithInvalidHeaders << "\n";
        std::cout << "Raw Frags With Invalid Anchors (HPGe)                     : " << HPGeEventMetrics.rawFragsWithInvalidAnchors << "\n";
        std::cout << "ZS Frags Skipped Due To Raw Flags (HPGe)                  : " << HPGeEventMetrics.zsFragsSkippedDueToRawFlag << "\n";
        std::cout << "PH Frags Skipped Due To Raw Flags (HPGe)                  : " << HPGeEventMetrics.phFragsSkippedDueToRawFlag << "\n";
        std::cout << "Raw/ZS/PH Sets Skipped Due To Raw Flags (HPGe)            : " << HPGeEventMetrics.setsSkippedDueToRawFlag << "\n";
        std::cout << "Raw/ZS/PH Sets Skipped Due To Invalid Raw Headers (HPGe)  : " << HPGeEventMetrics.setsSkippedDueToInvalidHeaders << "\n";
        std::cout << "\n";
        std::cout << "Raw Frags With Invalid Headers (LaBr)                     : " << LaBrEventMetrics.rawFragsWithInvalidHeaders << "\n";
        std::cout << "Raw Frags With Inavlid Anchors (LaBr)                     : " << LaBrEventMetrics.rawFragsWithInvalidAnchors << "\n";
        std::cout << "ZS Frags Skipped Due To Raw Flags (LaBr)                  : " << LaBrEventMetrics.zsFragsSkippedDueToRawFlag << "\n";
        std::cout << "PH Frags Skipped Due To Raw Flags (LaBr)                  : " << LaBrEventMetrics.phFragsSkippedDueToRawFlag << "\n";
        std::cout << "Raw/ZS/PH Sets Skipped Due To Raw Flags (LaBr)            : " << LaBrEventMetrics.setsSkippedDueToRawFlag << "\n";
        std::cout << "Raw/ZS/PH Sets Skipped Due To Invalid Raw Headers (LaBr)  : " << LaBrEventMetrics.setsSkippedDueToInvalidHeaders << "\n";
        std::cout << "=================================\n";

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
        std::cout << "\n========== STM JOB SUMMARY - (Unpacking Module) ==========\n";
        std::cout << "\n--- Module Configuration ---\n";
        std::cout << "Raw Waveforms with Header HPGe   : " << (_saveRawWaveformsWithHeaderHPGe ? "Yes" : "No") << "\n";
        std::cout << "Raw Waveforms HPGe               : " << (_saveRawWaveformsHPGe ? "Yes" : "No") << "\n";
        std::cout << "ZS Waveforms HPGe                : " << (_saveZSWaveformsHPGe ? "Yes" : "No") << "\n";
        std::cout << "PH Digis HPGe                    : Yes\n";
        std::cout << "Raw Waveforms with Header LaBr   : " << (_saveRawWaveformsWithHeaderLaBr ? "Yes" : "No") << "\n";
        std::cout << "Raw Waveforms LaBr               : " << (_saveRawWaveformsLaBr ? "Yes" : "No") << "\n";
        std::cout << "ZS Waveforms LaBr                : " << (_saveZSWaveformsLaBr ? "Yes" : "No") << "\n";
        std::cout << "PH Digis LaBr                    : Yes\n";

        // Container and Inner Fragment Summary
        std::cout << "\n--- Container and Inner Fragment Summary ---\n";
        std::cout << "Total Art Events Processed                      : " << _totalEvents << "\n";
        std::cout << "Total Art Events w/ HPGe & LaBr                 : " << _totalEventsWithBothDetectors << "\n";
        std::cout << "Total Art Events w/ HPGe Only                   : " << _totalEventsWithOnlyHPGe << "\n";
        std::cout << "Total Art Events w/ LaBr Only                   : " << _totalEventsWithOnlyLaBr << "\n";
        std::cout << "Total Art Events w/ Neither HPGe nor LaBr       : " << _totalEventsWithNeitherDetector << "\n";

        std::cout << "Total HPGe Containers                           : " << _totalContainersHPGe << "\n";
        std::cout << "Total LaBr Containers                           : " << _totalContainersLaBr << "\n";
        std::cout << "Total Container Processed                       : " << _totalContainers << "\n";
        std::cout << "Total Non Container Processed                   : " << _totalNonContainers << "\n";

        std::cout << "Total Inner Fragments Processed                 : " << _totalInnerFrags << "\n";
        std::cout << "Total Unreadable Inner Fragments                : " << _totalUnreadInnerFrags << "\n";

        // Data Type Summary - pre-filtering
        std::cout << "\n--- Data Types Read (Pre - Filtering) ---\n";
        std::cout << "Total Raw Frags Seen                            : " << _totalRawFragsSeen << "\n";
        std::cout << "Total ZS  Frags Seen                            : " << _totalZSFragsSeen << "\n";
        std::cout << "Total PH  Frags Seen                            : " << _totalPHFragsSeen << "\n";
        std::cout << "Total Raw Prescaled Frags                       : " << _totalRawFragsPrescaled << "\n";
        std::cout << "Total ZS Prescaled Frags                        : " << _totalZSFragsPrescaled << "\n";
        std::cout << "Total Raw Frags Flagged Bad Only                : " << _totalRawFragsFlaggedBadOnly << "\n";
        std::cout << "Total Raw Frags Flagged Missing Only            : " << _totalRawFragsFlaggedMissingOnly << "\n";
        std::cout << "Total Raw Frags Flagged Bad And Missing         : " << _totalRawFragsFlaggedBadAndMissing << "\n";
        std::cout << "\n";
        std::cout << "Total RAW Frags seen (HPGe)                     : " << _totalRawFragsSeenHPGe << "\n";
        std::cout << "Total ZS  Frags seen (HPGe)                     : " << _totalZSFragsSeenHPGe << "\n";
        std::cout << "Total PH  Frags seen (HPGe)                     : " << _totalPHFragsSeenHPGe << "\n";
        std::cout << "\n";
        std::cout << "Total RAW Frags seen (LaBr)                     : " << _totalRawFragsSeenLaBr << "\n";
        std::cout << "Total ZS  Frags seen (LaBr)                     : " << _totalZSFragsSeenLaBr << "\n";
        std::cout << "Total PH  Frags seen (LaBr)                     : " << _totalPHFragsSeenLaBr << "\n";

        // Data Type Summary - post-filtering (General)
        std::cout << "\n--- Data Types Read General (Post - Filtering) ---\n";
        std::cout << "Total Good  Raw Frags                           : " << _totalGoodRawFrags << "\n";
        std::cout << "Total Good  ZS  Frags                           : " << _totalGoodZSFrags << "\n";
        std::cout << "Total Good  PH  Frags                           : " << _totalGoodPHFrags << "\n";
        std::cout << "\n";
        std::cout << "Total Zero  Raw Frags                           : " << _totalZeroRawFrags << "\n";
        std::cout << "Total Zero  ZS  Frags                           : " << _totalZeroZSFrags << "\n";
        std::cout << "Total Zero  PH  Frags                           : " << _totalZeroPHFrags << "\n";
        std::cout << "\n";
        std::cout << "Total Empty Raw Frags                           : " << _totalEmptyRawFrags << "\n";
        std::cout << "Total Empty ZS  Frags                           : " << _totalEmptyZSFrags << "\n";
        std::cout << "Total Empty PH  Frags                           : " << _totalEmptyPHFrags << "\n";

        // Data Type Summary - post-filtering (HPGe)
        std::cout << "\n--- Data Types Read HPGe (Post - Filtering) ---\n";
        std::cout << "Total Good  Raw Frags (HPGe)                    : " << _totalGoodRawFragsHPGe << "\n";
        std::cout << "Total Good  ZS  Frags (HPGe)                    : " << _totalGoodZSFragsHPGe << "\n";
        std::cout << "Total Good  PH  Frags (HPGe)                    : " << _totalGoodPHFragsHPGe << "\n";
        std::cout << "\n";
        std::cout << "Total Zero  Raw Frags (HPGe)                    : " << _totalZeroRawFragsHPGe << "\n";
        std::cout << "Total Zero  ZS  Frags (HPGe)                    : " << _totalZeroZSFragsHPGe << "\n";
        std::cout << "Total Zero  PH  Frags (HPGe)                    : " << _totalZeroPHFragsHPGe << "\n";
        std::cout << "\n";
        std::cout << "Total Empty Raw Frags (HPGe)                    : " << _totalEmptyRawFragsHPGe << "\n";
        std::cout << "Total Empty ZS  Frags (HPGe)                    : " << _totalEmptyZSFragsHPGe << "\n";
        std::cout << "Total Empty PH  Frags (HPGe)                    : " << _totalEmptyPHFragsHPGe << "\n";
        std::cout << "\n";
        std::cout << "Total Raw Prescaled Frags (HPGe)                : " << _totalRawFragsPrescaledHPGe << "\n";
        // Raw Prescaled Frags should be classified as their own data for now
        std::cout << "Total ZS Prescaled Frags (HPGe)                 : " << _totalZSFragsPrescaledHPGe << "\n";
        std::cout << "Total Raw Frags Flagged Bad Only (HPGe)         : " << _totalRawFragsFlaggedBadOnlyHPGe << "\n";
        std::cout << "Total Raw Frags Flagged Missing Only (HPGe)     : " << _totalRawFragsFlaggedMissingOnlyHPGe << "\n";
        std::cout << "Total Raw Frags Flagged Bad & Missing (HPGe)    : " << _totalRawFragsFlaggedBadAndMissingHPGe << "\n";

        // Data Type Summary - post-filtering (LaBr)
        std::cout << "\n--- Data Types Read LaBr (Post - Filtering) ---\n";
        std::cout << "Total Good  Raw Frags (LaBr)                    : " << _totalGoodRawFragsLaBr << "\n";
        std::cout << "Total Good  ZS  Frags (LaBr)                    : " << _totalGoodZSFragsLaBr << "\n";
        std::cout << "Total Good  PH  Frags (LaBr)                    : " << _totalGoodPHFragsLaBr << "\n";
        std::cout << "\n";
        std::cout << "Total Zero  Raw Frags (LaBr)                    : " << _totalZeroRawFragsLaBr << "\n";
        std::cout << "Total Zero  ZS  Frags (LaBr)                    : " << _totalZeroZSFragsLaBr << "\n";
        std::cout << "Total Zero  PH  Frags (LaBr)                    : " << _totalZeroPHFragsLaBr << "\n";
        std::cout << "\n";
        std::cout << "Total Empty Raw Frags (LaBr)                    : " << _totalEmptyRawFragsLaBr << "\n";
        std::cout << "Total Empty ZS  Frags (LaBr)                    : " << _totalEmptyZSFragsLaBr << "\n";
        std::cout << "Total Empty PH  Frags (LaBr)                    : " << _totalEmptyPHFragsLaBr << "\n";
        std::cout << "\n";
        std::cout << "Total Raw Prescaled Frags (LaBr)                : " << _totalRawFragsPrescaledLaBr << "\n";
        std::cout << "Total ZS Prescaled Frags (LaBr)                 : " << _totalZSFragsPrescaledLaBr << "\n";
        std::cout << "Total Raw Frags Flagged Bad Only (LaBr)         : " << _totalRawFragsFlaggedBadOnlyLaBr << "\n";
        std::cout << "Total Raw Frags Flagged Missing Only (LaBr)     : " << _totalRawFragsFlaggedMissingOnlyLaBr << "\n";
        std::cout << "Total Raw Frags Flagged Bad & Missing (LaBr)    : " << _totalRawFragsFlaggedBadAndMissingLaBr << "\n";

        std::cout << "\n=== Extra Filters ===\n";
        std::cout << "Invalid Raw Headers                             : " << _totalRawFragsWithInvalidHeaders << "\n";
        std::cout << "Invalid Raw Headers (HPGe)                      : " << _totalRawFragsWithInvalidHeadersHPGe << "\n";
        std::cout << "Invalid Raw Headers (LaBr)                      : " << _totalRawFragsWithInvalidHeadersLaBr << "\n";
        std::cout << "Invalid Raw Anchors                             : " << _totalRawFragsWithInvalidAnchors << "\n";
        std::cout << "Invalid Raw Anchors (HPGe)                      : " << _totalRawFragsWithInvalidAnchorsHPGe << "\n";
        std::cout << "Invalid Raw Anchors (LaBr)                      : " << _totalRawFragsWithInvalidAnchorsLaBr << "\n";
        std::cout << "ZS Skipped Due To Raw Flags (HPGe)              : " << _totalZSFragsSkippedDueToRawFlagHPGe << "\n";
        std::cout << "ZS Skipped Due To Raw Flags (LaBr)              : " << _totalZSFragsSkippedDueToRawFlagLaBr << "\n";
        std::cout << "PH Skipped Due To Raw Flags (HPGe)              : " << _totalPHFragsSkippedDueToRawFlagHPGe << "\n";
        std::cout << "PH Skipped Due To Raw Flags (LaBr)              : " << _totalPHFragsSkippedDueToRawFlagLaBr << "\n";
        std::cout << "PH Count Mismatches Encountered                 : " << _totalPHCountMismatch << "\n";
        std::cout << "PH Count Mismatches Encountered (HPGe)          : " << _totalPHCountMismatchHPGe << "\n";
        std::cout << "PH Count Mismatches Encountered (LaBr)          : " << _totalPHCountMismatchLaBr << "\n";


        std::cout << "\n===========================================================\n";

    }
}

DEFINE_ART_MODULE(STMDigisFromFragments)
