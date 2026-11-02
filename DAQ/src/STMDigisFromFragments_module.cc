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
    void printJobSummary() const;
    void printDetailedJobSummary() const;
    void endJob() override;

private:
    // Art input tags
    // --------------
    art::InputTag _stmFragmentsTag;

    // Fhicl parameters
    // ----------------
    int _verbosityLevel{0};
    bool _saveSTMFragSummary{false};
    // HPGe
    bool _saveRawWaveformsWithHeaderHPGe{false};
    bool _saveRawWaveformsHPGe{false};
    bool _saveZSWaveformsHPGe{false};
    // LaBr
    bool _saveRawWaveformsWithHeaderLaBr{false};
    bool _saveRawWaveformsLaBr{false};
    bool _saveZSWaveformsLaBr{false};

    // General container frag / frag variables
    // ---------------------------------------
    size_t _totalEvents{0};                      // num of art evts
    size_t _totalContainers{0};                  // num of containers processed
    size_t _totalUnknownContainers{0};           // num of containers that are not hpge or labr
    size_t _totalNonContainers{0};               // num of non containers
    size_t _totalInnerFrags{0};                  // num of inner frags seen
    size_t _totalUnknownInnerFrags{0};           // num of unreadable inner frags
    size_t _totalEventsWithBothDetectors{0};     // num of art evts with both detectors
    size_t _totalEventsWithOnlyHPGe{0};          // num of art evts only with HPGe
    size_t _totalEventsWithOnlyLaBr{0};          // num of art evts only with LaBr
    size_t _totalEventsWithNeitherDetector{0};   // num of art evts with neither expected container but still interprted as

    // Counts for one fragment type for one detector
    struct fragCounts {
        size_t seen{0};        // frag seen and identified
        size_t prescaled{0};   // frag was prescaled (Raw & ZS ONLY)
        size_t empty{0};       // frag read an unexpected empty frag
        size_t zero{0};        // frag read an unexpected zero fill frag
        size_t good{0};        // frag had expected data
        size_t noRegions{0};   // frag reported no ZS regions
        size_t noHits{0};      // frag reported no PH Hits
        size_t unread{0};      // frag was unreadable (malformed)
    };

    // General counts that can be used for HPGe/LaBr & job/event
    struct detectorCounts {
        // set up frag
        fragCounts raw;
        fragCounts zs;
        fragCounts ph;

        // container & inner frag counts
        size_t container{0};
        size_t innerFrags{0};

        // Possible errors in raw frag
        size_t invalidHeaders{0};  // raw frag too short for its header
        size_t invalidAnchors{0};  // raw header start/end anchor check failed
        // Possible errors in Raw Header Flags read
        size_t badOnly{0};
        size_t missingOnly{0};
        size_t badAndMissing{0};
        // Possible errors related to zs
        size_t zsLengthMismatch{0};
        size_t zsRegionMismatch{0};
        size_t zsHeaderMismatch{0};
        // Possible errors related to ph
        size_t phCountMismatch{0};
        // Count how many ZS and PH frags were skipped due to InvalidHeader/Raw Flag
        size_t skippedZSRawFlags{0};
        size_t skippedPHRawFlags{0};
        size_t skippedZSInvalidHeader{0};
        size_t skippedPHInvalidHeader{0};
        size_t skippedZSNoRaw{0};
        size_t skippedPHNoRaw{0};

        // Count how many ZS,PH frags were skipped due to errors in raw
        size_t setsSkippedRawFlags{0};
        size_t setsSkippedInvalidHeaders{0};
        size_t setsSkippedInvalidAnchors{0};
    };

    // Used to save ZS information
    struct zsRegion{
        uint32_t offset;
        std::vector<int16_t> adcs;
    };

    // Used to track the raw header per set of (raw -> zs -> ph)
    struct rawHeaderState {
        //
        bool rawPrescaled{false};
        bool zsPrescaled{false};
        // Extract expected ZS and PH information
        uint16_t expectedZSLength{0};
        uint16_t expectedZSRegions{0};
        uint16_t expectedPHCount{0};
        // Used to determine whether raw header was read
        bool readZSInfo{false};
        bool readPHInfo{false};
        // Flag if raw header had an issue
        bool skipCurrentSetDueToRawFlags{false};
        bool skipCurrentSetDueToInvalidRawHeader{false};
        // Save for event header -> maybe only EWT?
        uint64_t eventWindowTag{0};
        uint8_t eventMode{0};
        uint64_t adcClock{0};
        uint64_t dtcClock{0};
        // Keep track of version type
        uint16_t headerVersion{0};
        bool isLegacyHeader{false};
        // Version dependent fields - only set for legacy header
        uint16_t rawPrescaleValue{0};
        uint16_t zsPrescaleValue{0};
    };

    // Used to track extract info after being added to map - used HPGe/ LaBr for job/event
    struct extractedInfo {
        size_t rawWithHeader{0};
        size_t rawWaveforms{0};
        size_t zsWaveforms{0};
        size_t phDigis{0};
    };

    //
    template <typename MapType>
    static size_t countDigis(const MapType& digiMap) {
        size_t total{0};
        for (const auto& entry : digiMap) {
            total += entry.second.size();
        }
        return total;
    }
    static mu2e::STMEventHeader makeEventHeader(const rawHeaderState& state) {
        return mu2e::STMEventHeader(state.eventWindowTag,
                                    state.eventMode,
                                    state.adcClock,
                                    state.dtcClock);
    }

    // Set up the job counters that will be used
    detectorCounts hpgeJobCounts;
    detectorCounts labrJobCounts;
    extractedInfo hpgeJobProducts;
    extractedInfo labrJobProducts;

    // Keeps track of raw header versions
    std::map<uint16_t, size_t> _rawHeaderVersions;

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
    // create products
    //----------------

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

    if (_saveSTMFragSummary){
        produces<mu2e::STMFragmentSummaryCollection>("stmFragSummaryHPGe");
        produces<mu2e::STMFragmentSummaryCollection>("stmFragSummaryLaBr");
    }
}

// Start of Event processing
// =========================
void STMDigisFromFragments::produce(Event& event)
{
    ++_totalEvents;

    // Event frag counters
    // -------------------
    size_t unknownFragsThisEvent{0};
    size_t unknownContainersThisEvent{0};
    size_t nonContainersThisEvent{0};

    // Set up
    rawHeaderState rawHeaderHPGe;
    rawHeaderState rawHeaderLaBr;

    // create null ptrs
    //-----------------
    // HPGe
    std::unique_ptr<mu2e::STMWaveformDigiCollectionMap>  rawWaveformDigisWithHeaderHPGe;
    std::unique_ptr<mu2e::STMWaveformDigiCollectionMap>  rawWaveformDigisHPGe;
    std::unique_ptr<mu2e::STMWaveformDigiCollectionMap>  zsWaveformDigisHPGe;
    std::unique_ptr<mu2e::STMPHDigiCollectionMap>        phDigisHPGe;
    std::unique_ptr<mu2e::STMFragmentSummaryCollection>  stmFragSummaryHPGe;

    // LaBr
    std::unique_ptr<mu2e::STMWaveformDigiCollectionMap>  rawWaveformDigisWithHeaderLaBr;
    std::unique_ptr<mu2e::STMWaveformDigiCollectionMap>  rawWaveformDigisLaBr;
    std::unique_ptr<mu2e::STMWaveformDigiCollectionMap>  zsWaveformDigisLaBr;
    std::unique_ptr<mu2e::STMPHDigiCollectionMap>        phDigisLaBr;
    std::unique_ptr<mu2e::STMFragmentSummaryCollection>  stmFragSummaryLaBr;

    // Assign ptrs to map (based on config)
    //------------------------------------
    // HPGe
    if (_saveRawWaveformsWithHeaderHPGe) {rawWaveformDigisWithHeaderHPGe = std::make_unique<mu2e::STMWaveformDigiCollectionMap>();}
    if (_saveRawWaveformsHPGe) {rawWaveformDigisHPGe = std::make_unique<mu2e::STMWaveformDigiCollectionMap>();}
    if (_saveZSWaveformsHPGe)  {zsWaveformDigisHPGe = std::make_unique<mu2e::STMWaveformDigiCollectionMap>();}
    phDigisHPGe = std::make_unique<mu2e::STMPHDigiCollectionMap>();

    // LaBr
    if (_saveRawWaveformsWithHeaderLaBr) {rawWaveformDigisWithHeaderLaBr = std::make_unique<mu2e::STMWaveformDigiCollectionMap>();}
    if (_saveRawWaveformsLaBr) {rawWaveformDigisLaBr = std::make_unique<mu2e::STMWaveformDigiCollectionMap>();}
    if (_saveZSWaveformsLaBr)  {zsWaveformDigisLaBr = std::make_unique<mu2e::STMWaveformDigiCollectionMap>();}
    phDigisLaBr = std::make_unique<mu2e::STMPHDigiCollectionMap>();

    // Frag Summaries
    if (_saveSTMFragSummary) {
        stmFragSummaryHPGe = std::make_unique<mu2e::STMFragmentSummaryCollection>();
        stmFragSummaryLaBr = std::make_unique<mu2e::STMFragmentSummaryCollection>();
    }

    // track what detectors are present per event
    // ------------------------------------------
    bool eventHasHPGe{false};
    bool eventHasLaBr{false};

    // Initialize event counter
    // ------------------------
    detectorCounts hpgeEventCounts;
    detectorCounts labrEventCounts;

    // Get STM Fragments from the event, via handle
    // ----------------------------------------------
    art::Handle<artdaq::Fragments> STMFragmentsHandle;
    event.getByLabel(_stmFragmentsTag, STMFragmentsHandle);
    const auto& STMFragments = STMFragmentsHandle.product();

    // Loop over outer fragments
    // -------------------------
    for (const auto& frag : *STMFragments) {
        const auto outerFragID = frag.fragmentID();
        if (_verbosityLevel >= 5) {
            mf::LogDebug("STMDigisFromFragments")
            << "Processing outer fragment with ID : "
            << " | Event : " << _totalEvents
            << " | Outer Frag ID : " << outerFragID;
        }

        // Check if this is a valid container fragment
        // -------------------------------------------
        if (frag.type() != artdaq::Fragment::ContainerFragmentType) {
            mf::LogWarning("STMDigisFromFragments")
            << "Non-container STM Fragment "
            << " | Event : " << _totalEvents
            << " | Outer Frag ID : " << outerFragID;
            ++_totalNonContainers;
            ++nonContainersThisEvent;
            continue;
        }

        const mu2e::STMFragment containerFrag(frag);
        artdaq::ContainerFragment container(frag);
        const auto blocks = container.block_count();

        // Check whether this is an expected container Frag
        // ------------------------------------------------
        if (containerFrag.isHPGeContainer()) {
            if (_verbosityLevel >= 6) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "Processing HPGe Container Frag"
                    << " | Event : " << _totalEvents
                    << " | Container Frag ID : " << outerFragID;
                }
                ++hpgeJobCounts.container;
                ++hpgeEventCounts.container;
                hpgeJobCounts.innerFrags += blocks;
                hpgeEventCounts.innerFrags += blocks;
                eventHasHPGe = true;
        } else if (containerFrag.isLaBrContainer()) {
            if (_verbosityLevel >= 6) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "Processing LaBr Container Frag"
                    << " | Event : " << _totalEvents
                    << " | Container Frag ID : " << outerFragID;
                }
                ++labrJobCounts.container;
                ++labrEventCounts.container;
                labrJobCounts.innerFrags += blocks;
                labrEventCounts.innerFrags += blocks;
                eventHasLaBr = true;
        } else {
            // Unknown Container Frag
            mf::LogWarning("STMDigisFromFragments")
                << "Unknown STM Container Frag"
                << " | Frag ID : " << outerFragID
                << " | Event : " << _totalEvents
                << " | Container Frag ID : " << outerFragID;

                ++_totalUnknownContainers;
                ++unknownContainersThisEvent;
                continue; // Skips the rest of this container frag if it is unknown
        }

        ++_totalContainers;
        _totalInnerFrags += blocks;

        // loop over container blocks
        // --------------------------
        for (size_t i = 0 ; i < container.block_count() ; ++i) {
            auto innerFrag = container.at(i);
            const mu2e::STMFragment stmFrag(*innerFrag);
            const bool isHPGe = stmFrag.isHPGe();
            const bool isLaBr = stmFrag.isLaBr();

            // Check for unknown inner frag
            if (!stmFrag.isRaw() && !stmFrag.isZS() && !stmFrag.isPH()){
                mf::LogWarning("STMDigisFromFragments")
                << "Unknown STM Inner Fragment"
                << " | Event : " << _totalEvents
                << " | Frag Index : " << i
                << " | Frag ID : " << innerFrag->fragmentID();

                ++_totalUnknownInnerFrags; //Job Summary Counter
                ++unknownFragsThisEvent;
                continue;
            }

            // Check if the inner frag falls under hpge or labr
            if (!isHPGe && !isLaBr) {
                mf::LogWarning("STMDigisFromFragments")
                << "STM Inner Fragment Is Neither HPGe nor LaBr"
                << " | Event : " << _totalEvents
                << " | Frag Index : " << i;

                ++_totalUnknownInnerFrags; ++unknownFragsThisEvent;
                continue;
            }

            // Create the relevant structs reference
            auto& eventCounts   = isHPGe ? hpgeEventCounts : labrEventCounts;
            auto& jobCounts     = isHPGe ? hpgeJobCounts   : labrJobCounts;
            auto& headerState   = isHPGe ? rawHeaderHPGe   : rawHeaderLaBr;

            // Note : Expected Inner Frag Processes (Raw -> ZS -> PH)
            if (stmFrag.isRaw()){
                // --------------
                // Check Raw Frag
                // --------------
                ++eventCounts.raw.seen; ++jobCounts.raw.seen;
                headerState = rawHeaderState{}; // resets raw header state upon entering raw

                // Check if header + adcs is not less than expected header length (version dependent)
                if(stmFrag.dataWords() < stmFrag.headerWords()) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "Raw Frag Has Fewer Data Words Than Expected"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i
                    << " | Version : " << stmFrag.headerVersion();

                    ++eventCounts.invalidHeaders; ++jobCounts.invalidHeaders;
                    ++eventCounts.setsSkippedInvalidHeaders; ++jobCounts.setsSkippedInvalidHeaders;
                    ++eventCounts.raw.unread; ++jobCounts.raw.unread;
                    headerState.skipCurrentSetDueToInvalidRawHeader = true;
                    continue;
                }
                ++_rawHeaderVersions[stmFrag.headerVersion()]; // Creates a map and count of headerVersion encountered (per job)
                // Confirm Raw Header has Valid Anchors
                // For now just use as a diagnostic tool
                if (!stmFrag.hasValidAnchors()) {
                    if (_verbosityLevel >= 6) {
                        const auto* dataPtr = stmFrag.dataBegin();
                        const size_t anchorEndIndex = stmFrag.isLegacyHeader() ? static_cast<size_t>(stm::RawHeaderLegacy::ANCHOR_END)
                                                                               : static_cast<size_t>(stm::RawHeader::ANCHOR_END);
                        std::ostringstream msg;
                        msg << std::hex // switch to hexadecimal
                        << "RAW ANCHORS : start = 0x" << static_cast<uint16_t>(dataPtr[stm::RawHeader::ANCHOR_START])
                        << " , end = 0x" << static_cast<uint16_t>(dataPtr[anchorEndIndex]) << std::dec // swithces back to decimal
                        << " | Version History     : " << stmFrag.headerVersion()
                        << " | LegacyHeader        : " << std::boolalpha << stmFrag.isLegacyHeader();
                        msg << "\n";

                        mf::LogDebug("STMDigisFromFragments")
                        << msg.str();
                    }
                    ++eventCounts.invalidAnchors; ++jobCounts.invalidAnchors;
                    //++eventCounts.setsSkippedInvalidAnchors; ++jobCounts.setsSkippedInvalidAnchors;

                    // first check if new files artdaq-core pr fixed this
                    //headerState.skipCurrentSetDueToInvalidHeader = true;
                    //continue; // debug first
                }

                // Check for bad &/or missing flags
                const bool bad = stmFrag.badData();
                const bool missing = stmFrag.missing();
                if (bad && missing) {
                    ++eventCounts.badAndMissing; ++jobCounts.badAndMissing;
                } else if (bad) {
                    ++eventCounts.badOnly; ++jobCounts.badOnly;
                } else if (missing) {
                    ++eventCounts.missingOnly; ++jobCounts.missingOnly;
                }

                if (bad || missing) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "RAW flagged " << (bad && missing ? "bad and missing" : bad ? "bad" : "missing")
                    << " | Event: " << _totalEvents
                    << " | Index: " << i;
                    ++eventCounts.setsSkippedRawFlags; ++jobCounts.setsSkippedRawFlags;
                    headerState.skipCurrentSetDueToRawFlags = true;
                    continue; // skip following frags if raw was malformed
                }

                // Save to rawHeader struct
                headerState.rawPrescaled = stmFrag.rawPrescaled();
                headerState.zsPrescaled = stmFrag.zsPrescaled();
                headerState.expectedZSLength = stmFrag.zsLength();
                headerState.expectedZSRegions = stmFrag.zsRegions();
                headerState.expectedPHCount = stmFrag.phCount();
                headerState.readZSInfo = true;
                headerState.readPHInfo = true;
                headerState.eventWindowTag = stmFrag.eventWindowTag();
                headerState.eventMode = stmFrag.spillFlag();
                headerState.adcClock = stmFrag.adcClock();
                headerState.dtcClock = stmFrag.dtcClock();
                headerState.headerVersion = stmFrag.headerVersion();
                headerState.isLegacyHeader = stmFrag.isLegacyHeader();

                // save only if reading legacy header
                if (headerState.isLegacyHeader) {
                    headerState.rawPrescaleValue = stmFrag.rawPrescaleValue();
                    headerState.zsPrescaleValue  = stmFrag.zsPrescaleValue();
                }

                if (_verbosityLevel >= 6) {
                    std::ostringstream msg;
                    msg << "RAW header | Event: " << _totalEvents << " | Index: " << i
                        << " | Detector: " << (isHPGe ? "HPGe" : "LaBr")
                        << " | Version: " << headerState.headerVersion
                        << " | Legacy: " << std::boolalpha << headerState.isLegacyHeader
                        << " | RAW prescaled: " << headerState.rawPrescaled
                        << " | ZS prescaled: " << headerState.zsPrescaled
                        << " | ZS regions: " << headerState.expectedZSRegions
                        << " | PH count: " << headerState.expectedPHCount;
                    if (headerState.isLegacyHeader) {
                        msg << " | RAW prescale value: " << headerState.rawPrescaleValue
                            << " | ZS prescale value: " << headerState.zsPrescaleValue;
                    }
                    mf::LogDebug("STMDigisFromFragments") << msg.str();
                }

                if (headerState.rawPrescaled) {
                    if(_verbosityLevel >= 5) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "Raw Fragment is prescaled"
                        << " | Event : " << _totalEvents
                        << " | Frag Index : " << i;
                    }
                        ++eventCounts.raw.prescaled; ++jobCounts.raw.prescaled;
                        continue; // Only count and report, prescaled data is still useful for zs and ph
                    }

                const auto payloadWords = stmFrag.payloadWords();
                // Should be header Words + adcs < = dataWords (version dependent)
                if (stmFrag.headerWords() + payloadWords > stmFrag.dataWords()) {
                    // this is an invalid-header case
                    mf::LogWarning("STMDigisFromFragments")
                    << "Raw Fragment Exceeds Fragment Data"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;

                    ++eventCounts.invalidHeaders; ++jobCounts.invalidHeaders;
                    ++eventCounts.setsSkippedInvalidHeaders; ++jobCounts.setsSkippedInvalidHeaders;
                    ++eventCounts.raw.unread; ++jobCounts.raw.unread;
                    headerState.readZSInfo = false;
                    headerState.readPHInfo = false;
                    headerState.skipCurrentSetDueToInvalidRawHeader = true;
                    continue;
                }

                // Ideally at this point we know the frag's payload can be inspected
                // payloadWords -> uses rawLength()
                if (payloadWords == 0) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "Unexpected Empty Raw Frag"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;
                    ++eventCounts.raw.empty; ++jobCounts.raw.empty;
                    continue;
                }

                // Check if any payload words exist
                const auto* payloadPtr = stmFrag.payloadBegin();
                const bool allZeros = std::all_of(payloadPtr, payloadPtr + payloadWords, [](int16_t adc) {return adc == 0;});

                if (allZeros) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "Unexpected Zero-Filled Raw Frag"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;
                    ++eventCounts.raw.zero; ++jobCounts.raw.zero;
                    continue;
                }

                // Print first few payload words for inspection
                if (_verbosityLevel >= 6) {
                    std::ostringstream msg;
                    msg << "First Few Payload Words(adcs) For Inspection : " ;
                    for (size_t w = 0 ; w < std::min(payloadWords, static_cast<size_t>(10)) ; ++w) {
                        msg << payloadPtr[w] << " ";
                    }
                    mf::LogDebug("STMDigisFromFragments")
                    << msg.str();
                }

                // At this point we know raw fragment has data
                if (_verbosityLevel >= 5) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "Good Raw Fragment"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;
                }
                // Good raw frags end up here and get saved
                ++eventCounts.raw.good; ++jobCounts.raw.good;
                const auto key = makeEventHeader(headerState);

                if (isHPGe) {
                    if (_saveRawWaveformsWithHeaderHPGe) {
                        mu2e::STMWaveformDigi waveformDigi;
                        waveformDigi.set_data(stmFrag.dataWords(), stmFrag.dataBegin());
                        (*rawWaveformDigisWithHeaderHPGe)[key].push_back(waveformDigi);
                    }
                    if (_saveRawWaveformsHPGe) {
                        mu2e::STMWaveformDigi waveformDigi;
                        waveformDigi.set_data(payloadWords, payloadPtr);
                        (*rawWaveformDigisHPGe)[key].push_back(waveformDigi);
                    }
                } else {
                    if (_saveRawWaveformsWithHeaderLaBr) {
                        mu2e::STMWaveformDigi waveformDigi;
                        waveformDigi.set_data(stmFrag.dataWords(), stmFrag.dataBegin());
                        (*rawWaveformDigisWithHeaderLaBr)[key].push_back(waveformDigi);
                    }
                    if (_saveRawWaveformsLaBr) {
                        mu2e::STMWaveformDigi waveformDigi;
                        waveformDigi.set_data(payloadWords, payloadPtr);
                        (*rawWaveformDigisLaBr)[key].push_back(waveformDigi);
                    }
                }

            } // end of raw
            else if (stmFrag.isZS()) {
                // -------------
                // Check ZS Frag
                // -------------
                ++eventCounts.zs.seen; ++jobCounts.zs.seen;

                // Skip corresponding zs frag?
                if (headerState.skipCurrentSetDueToInvalidRawHeader) {
                    ++eventCounts.skippedZSInvalidHeader; ++jobCounts.skippedZSInvalidHeader;
                    continue;
                }
                if (headerState.skipCurrentSetDueToRawFlags) {
                    ++eventCounts.skippedZSRawFlags; ++jobCounts.skippedZSRawFlags;
                    continue;
                }
                if (!headerState.readZSInfo) {
                    ++eventCounts.skippedZSNoRaw; ++jobCounts.skippedZSNoRaw;
                    continue;
                }
                if (headerState.zsPrescaled) {
                    if(_verbosityLevel >= 5) {
                            mf::LogDebug("STMDigisFromFragments")
                            << "ZS Frag is Prescaled"
                            << " | Event : " << _totalEvents
                            << " | Frag Index : " << i;
                    }
                    ++eventCounts.zs.prescaled; ++jobCounts.zs.prescaled;
                    continue;
                }

                const auto zsRegions = headerState.expectedZSRegions;
                const auto zsLength = headerState.expectedZSLength;
                if (zsRegions == 0 && zsLength == 0) {
                    if(_verbosityLevel >=5) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "ZS Regions & Length Were Zero"
                        << " | ZS Regions : " << zsRegions
                        << " | ZS Length  : " << zsLength
                        << " | Event : " << _totalEvents
                        << " | Frag Index : " << i;
                    }
                    ++eventCounts.zs.noRegions; ++jobCounts.zs.noRegions;
                    continue;
                }

                if ((zsRegions == 0 && zsLength > 0) || (zsRegions > 0 && zsLength == 0)) {
                    if(_verbosityLevel >= 5) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "ZS Header Mismatch"
                        << " | ZS Regions : " << zsRegions
                        << " | ZS Length : " << zsLength
                        << " | Event : " << _totalEvents << " | Frag Index : " << i;
                    }
                    ++eventCounts.zsHeaderMismatch; ++jobCounts.zsHeaderMismatch;
                    continue;
                }

                if (stmFrag.dataWords() <= stm::ZSHeader::WORDS) {
                    // We expect the frag to have at least two words at this point, index and length for zs header
                    mf::LogWarning("STMDigisFromFragments")
                    << "Unexpected Empty ZS Frag"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;
                    ++eventCounts.zs.empty; ++ jobCounts.zs.empty;
                    continue;
                }

                // Definitions
                const auto* dataPtr = stmFrag.dataBegin();
                const auto* endPtr = dataPtr + stmFrag.dataWords();
                size_t regionCounter = 0;
                size_t totalLength = 0;
                bool malformed = false;
                bool sawADC = false;
                bool allADCsZero = true;
                const bool saveZS = isHPGe ? (_saveZSWaveformsHPGe) : (_saveZSWaveformsLaBr);
                std::vector<zsRegion> regions;
                const auto payloadPtr = stmFrag.payloadBegin();


                if (_verbosityLevel >= 6) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "ZS Data Words for inspection : " << stmFrag.dataWords()
                    << "ZS Data Words % 4 : " << stmFrag.dataWords() % 4;
                }

                while (static_cast<size_t>(endPtr - dataPtr) >= stm::ZSHeader::WORDS &&
                       regionCounter < zsRegions) {
                    const uint16_t index = static_cast<uint16_t>(dataPtr[0]);
                    const uint16_t length = static_cast<uint16_t>(dataPtr[1]);
                    const auto* adcPtr = dataPtr + stm::ZSHeader::WORDS;
                    if (static_cast<size_t>(endPtr - adcPtr) < length) {
                        malformed = true;
                        break;
                    }
                    // print first couple adcs for inspection
                    if (_verbosityLevel >=6) {
                        std::ostringstream msg;
                        msg << "\nFirst Few Payload Words For ZS Inspection: ";
                        for (size_t w = 0; w < std::min<size_t>(length,10) ; ++w) {
                            msg << payloadPtr[w] << " ";
                        }
                        mf::LogDebug("STMDigisFromFragments")
                        << msg.str();
                    }

                    for (size_t j = 0; j < length; ++j) {
                        sawADC = true;
                        if (adcPtr[j] != 0) {allADCsZero = false;}
                    }
                    if (saveZS) {
                        regions.push_back({index, std::vector<int16_t>(adcPtr, adcPtr + length)});
                    }
                    if (_verbosityLevel >= 6) {
                        mf::LogDebug("STMDigisFromFragments")
                        << "ZS region"
                        << " | Event: " << _totalEvents
                        << " | Index: " << i
                        << " | Offset: " << index
                        << " | Length: " << length;
                    }
                    totalLength += length;
                    ++regionCounter;
                    dataPtr = adcPtr + length;
                }

                if (malformed) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "Malformed ZS Frag"
                    << " | Event: " << _totalEvents
                    << " | Index: " << i;
                    ++eventCounts.zs.unread; ++jobCounts.zs.unread;
                    continue;
                }
                if (totalLength != zsLength) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "ZS length mismatch"
                    << " | Event: " << _totalEvents
                    << " | Index: " << i
                    << " | Expected: " << zsLength
                    << " | Parsed: " << totalLength;
                    ++eventCounts.zsLengthMismatch; ++jobCounts.zsLengthMismatch;
                    continue;
                }
                if (regionCounter != zsRegions) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "ZS region mismatch:"
                    << " | Event: " << _totalEvents
                    << " | Index: " << i
                    << " | Expected: " << zsRegions
                    << " | Parsed: " << regionCounter;
                    ++eventCounts.zsRegionMismatch; ++jobCounts.zsRegionMismatch;
                    continue;
                }
                if (sawADC && allADCsZero) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "Unexpected Zero-filled ZS fragment"
                    << " | Event: " << _totalEvents
                    << " | Index: " << i;
                    ++eventCounts.zs.zero; ++jobCounts.zs.zero;
                    continue;
                }

                if (saveZS) {
                    const auto key = makeEventHeader(headerState);
                    auto& target = isHPGe ? *zsWaveformDigisHPGe : *zsWaveformDigisLaBr;
                    for (auto& region : regions) {
                        mu2e::STMWaveformDigi zsDigi(region.offset, region.adcs);
                        target[key].push_back(zsDigi);
                    }
                }

                // At this point we know the fragment was good
                if (_verbosityLevel >= 5) {
                    mf::LogDebug("STMDigisFromFragments")
                    << "Good ZS fragment at Event"
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;
                }
                ++eventCounts.zs.good; ++jobCounts.zs.good;
            } // end of zs
            else if (stmFrag.isPH()){
                // ---------
                // Check PH Frag
                // ---------
                ++eventCounts.ph.seen; ++jobCounts.ph.seen;
                if (headerState.skipCurrentSetDueToInvalidRawHeader) {
                    ++eventCounts.skippedPHInvalidHeader; ++jobCounts.skippedPHInvalidHeader;
                    continue;
                }
                if (headerState.skipCurrentSetDueToRawFlags) {
                    ++eventCounts.skippedPHRawFlags; ++jobCounts.skippedPHRawFlags;
                    continue;
                }
                if (!headerState.readPHInfo) {
                    ++eventCounts.skippedPHNoRaw; ++jobCounts.skippedPHNoRaw;
                    continue;
                }

                const auto phCount = headerState.expectedPHCount;
                if (phCount == 0) {
                    ++eventCounts.ph.noHits; ++jobCounts.ph.noHits;
                    continue;
                }
                const auto payloadWords = stmFrag.payloadWords();
                if (payloadWords == 0) {
                    mf::LogWarning("STMDigisFromFragments")
                    << "Unexpected Empty PH Frag "
                    << " | Event : " << _totalEvents
                    << " | Frag Index : " << i;
                    ++eventCounts.ph.empty; ++jobCounts.ph.empty;
                    ++eventCounts.ph.unread; ++jobCounts.ph.unread;
                    continue;
                }
                // Check if we under count ph Pairs in fragment compared to expected
                const size_t expectedWords = 2 * static_cast<size_t>(phCount);
                if (payloadWords < expectedWords) {
                    mf::LogWarning("STMDigisFromFragments")
                        << "PH Count Mismatch"
                        << " | Event : " << _totalEvents
                        << " | Index : " << i
                        << " | Expected Words : " << expectedWords
                        << " | Payload Words  : " << payloadWords;
                    ++eventCounts.phCountMismatch; ++jobCounts.phCountMismatch;
                    continue;
                }
                const auto* payloadPtr = stmFrag.payloadBegin();
                bool allPHZero = true;
                // Check if second entry has unexpected zero - note (time,ph)
                for (size_t k = 0; k < phCount; ++k) {
                    if (payloadPtr[k * 2 + 1] != 0) {
                        allPHZero = false;
                        break;
                    }
                }

                if (allPHZero) {
                    mf::LogWarning("STMDigisFromFragments")
                        << "Unexpected Zero-Filled PH Fragment"
                        << " | Event : " << _totalEvents
                        << " | Frag Index : " << i;
                    ++eventCounts.ph.zero; ++jobCounts.ph.zero;
                    continue;
                }

                const auto key = makeEventHeader(headerState);
                auto& target = isHPGe ? *phDigisHPGe : *phDigisLaBr;
                for (size_t l = 0 ; l < phCount ; ++l) {
                    const auto phIndex = 2 * l;
                    const uint32_t time = static_cast<uint16_t>(payloadPtr[phIndex]);
                    const int16_t ph = payloadPtr[phIndex + 1];
                    if (_verbosityLevel >= 6) {
                        std::ostringstream msg;
                        msg << "Time : " << time << " | PH : " << ph;
                        mf::LogDebug("STMDigisFromFragments") << msg.str();
                    }

                    mu2e::STMPHDigi phDigi(time,ph);
                    target[key].push_back(phDigi); // vs emplace_back(phDigi) ?
                }


                // At this point we have a good PH fragment
                if (_verbosityLevel >= 5) {
                    mf::LogDebug("STMDigisFromFragments")
                        << "Good PH Fragment at Event "
                        << " | Event : " << _totalEvents
                        << " | Frag Index : " << i;
                }
                // Good counter increment
                ++eventCounts.ph.good; ++jobCounts.ph.good;

                } // end of ph
            }
    } // end of frags loop

    // Update event type counters based on the detectors in event
    if (eventHasHPGe && eventHasLaBr) {++_totalEventsWithBothDetectors;}
    else if (eventHasHPGe)            {++_totalEventsWithOnlyHPGe;}
    else if (eventHasLaBr)            {++_totalEventsWithOnlyLaBr;}
    else                              {++_totalEventsWithNeitherDetector;}

    // Save STMFrag Summary here
    if (_saveSTMFragSummary) {
        if (_verbosityLevel >= 5) {
            mf::LogDebug("STMDigisFromFragments")
            << "Saving STMFragSummary";
        }
        auto appendSummary = [](mu2e::STMFragmentSummaryCollection& output, const detectorCounts& c) {
            output.emplace_back(
                c.container, c.innerFrags,
                c.badOnly + c.badAndMissing, c.missingOnly + c.badAndMissing,
                c.skippedZSRawFlags, c.skippedPHRawFlags,
                c.raw.prescaled, c.zs.prescaled,
                c.raw.good, c.zs.good, c.ph.good,
                c.raw.zero, c.zs.zero, c.ph.zero,
                c.raw.empty, c.zs.empty, c.ph.empty);
        };
        appendSummary(*stmFragSummaryHPGe, hpgeEventCounts);
        appendSummary(*stmFragSummaryLaBr, labrEventCounts);
    }

    // Extract info on how many waveforms and digi collections were extracted
    // ----------------------------------------------------------------------
    // Event Based
    extractedInfo hpgeEventProducts;
    extractedInfo labrEventProducts;
    // hpge
    if (rawWaveformDigisWithHeaderHPGe) {hpgeEventProducts.rawWithHeader = countDigis(*rawWaveformDigisWithHeaderHPGe);}
    if (rawWaveformDigisHPGe) {hpgeEventProducts.rawWaveforms = countDigis(*rawWaveformDigisHPGe);}
    if (zsWaveformDigisHPGe) {hpgeEventProducts.zsWaveforms = countDigis(*zsWaveformDigisHPGe);}
    hpgeEventProducts.phDigis = countDigis(*phDigisHPGe);
    // labr
    if (rawWaveformDigisWithHeaderLaBr) {labrEventProducts.rawWithHeader = countDigis(*rawWaveformDigisWithHeaderLaBr);}
    if (rawWaveformDigisLaBr) {labrEventProducts.rawWaveforms = countDigis(*rawWaveformDigisLaBr);}
    if (zsWaveformDigisLaBr) {labrEventProducts.zsWaveforms = countDigis(*zsWaveformDigisLaBr);}
    labrEventProducts.phDigis = countDigis(*phDigisLaBr);


    // Job Based - Accumulate over all events
    // hpge
    if (rawWaveformDigisWithHeaderHPGe) {hpgeJobProducts.rawWithHeader += countDigis(*rawWaveformDigisWithHeaderHPGe);}
    if (rawWaveformDigisHPGe) {hpgeJobProducts.rawWaveforms += countDigis(*rawWaveformDigisHPGe);}
    if (zsWaveformDigisHPGe) {hpgeJobProducts.zsWaveforms += countDigis(*zsWaveformDigisHPGe);}
    hpgeJobProducts.phDigis += countDigis(*phDigisHPGe);
    // labr
    if (rawWaveformDigisWithHeaderLaBr) {labrJobProducts.rawWithHeader += countDigis(*rawWaveformDigisWithHeaderLaBr);}
    if (rawWaveformDigisLaBr) {labrJobProducts.rawWaveforms += countDigis(*rawWaveformDigisLaBr);}
    if (zsWaveformDigisLaBr) {labrJobProducts.zsWaveforms += countDigis(*zsWaveformDigisLaBr);}
    labrJobProducts.phDigis += countDigis(*phDigisLaBr);

    // Compact Event Summary
    if (_verbosityLevel >= 3) {
        std::ostringstream out;
        auto row = [&](const std::string& name, size_t hpge, size_t labr) {
            out << std::left << std::setw(32) << name
                << std::right << std::setw(12) << hpge << std::setw(12) << labr
                << std::setw(12) << hpge + labr << '\n';
        };
        out << "\n========== STM EVENT SUMMARY (Unpacking Module) ==========\n"
            << "Art event: " << event.id() << " | Module event count: " << _totalEvents << '\n'
            << "--- Products Saved (event total) ---\n"
            << "Raw with header |  HPGe : " << hpgeEventProducts.rawWithHeader << " , LaBr : " << labrEventProducts.rawWithHeader << '\n'
            << "RAW waveforms   |  HPGe : " << hpgeEventProducts.rawWaveforms  << " , LaBr : " << labrEventProducts.rawWaveforms << '\n'
            << "ZS waveforms    |  HPGe : " << hpgeEventProducts.zsWaveforms   << " , LaBr : " << labrEventProducts.zsWaveforms << '\n'
            << "PH digis        |  HPGe : " << hpgeEventProducts.phDigis       << " , LaBr : " << labrEventProducts.phDigis << "\n\n";
        out << std::left << std::setw(32) << "Fragment metric"
            << std::right << std::setw(12) << "HPGe" << std::setw(12) << "LaBr"
            << std::setw(12) << "Total" << '\n';
        out << std::string(68, '-') << '\n';
        row("Containers", hpgeEventCounts.container, labrEventCounts.container);
        row("Inner fragments", hpgeEventCounts.innerFrags, labrEventCounts.innerFrags);
        row("Good RAW fragments", hpgeEventCounts.raw.good, labrEventCounts.raw.good);
        row("Good ZS fragments", hpgeEventCounts.zs.good, labrEventCounts.zs.good);
        row("Good PH fragments", hpgeEventCounts.ph.good, labrEventCounts.ph.good);
        row("No ZS regions", hpgeEventCounts.zs.noRegions, labrEventCounts.zs.noRegions);
        row("No PH hits", hpgeEventCounts.ph.noHits, labrEventCounts.ph.noHits);
        mf::LogInfo("STMDigisFromFragments") << out.str();
    }

    if (_verbosityLevel >= 4) {
        std::ostringstream out;
        auto row = [&](const std::string& name, size_t hpge, size_t labr) {
            out << std::left << std::setw(36) << name
                << std::right << std::setw(12) << hpge << std::setw(12) << labr
                << std::setw(12) << hpge + labr << '\n';
        };
        auto header = [&](const std::string& name) {
        out << std::left << std::setw(36) << name
            << std::right << std::setw(12) << "HPGe"
            << std::setw(12) << "LaBr"
            << std::setw(12) << "Total" << '\n'
            << std::string(72, '-') << '\n';
        };

        out << "\n======= STM Detailed Event Summary =======\n"
            << "Unknown inner fragments  : " << unknownFragsThisEvent << '\n'
            << "Unknown containers       : " << unknownContainersThisEvent << '\n'
            << "Non-container fragments  : " << nonContainersThisEvent << "\n\n";
        out << "\n";

        header("Raw Diagnostic");
        row("Invalid headers", hpgeEventCounts.invalidHeaders, labrEventCounts.invalidHeaders);
        row("Invalid anchors", hpgeEventCounts.invalidAnchors, labrEventCounts.invalidAnchors);
        row("Flagged bad only", hpgeEventCounts.badOnly, labrEventCounts.badOnly);
        row("Flagged missing only", hpgeEventCounts.missingOnly, labrEventCounts.missingOnly);
        row("Flagged bad and missing", hpgeEventCounts.badAndMissing, labrEventCounts.badAndMissing);
        out << "\n";

        header("ZS Diagnostic");
        row("Header mismatches", hpgeEventCounts.zsHeaderMismatch, labrEventCounts.zsHeaderMismatch);
        row("Length mismatches", hpgeEventCounts.zsLengthMismatch, labrEventCounts.zsLengthMismatch);
        row("Region mismatches", hpgeEventCounts.zsRegionMismatch, labrEventCounts.zsRegionMismatch);
        row("Skipped: RAW flags", hpgeEventCounts.skippedZSRawFlags, labrEventCounts.skippedZSRawFlags);
        row("Skipped: invalid RAW header", hpgeEventCounts.skippedZSInvalidHeader, labrEventCounts.skippedZSInvalidHeader);
        row("Skipped: no preceding RAW", hpgeEventCounts.skippedZSNoRaw, labrEventCounts.skippedZSNoRaw);
        out << "\n";

        header("PH Diagnostic");
        row("Count mismatches", hpgeEventCounts.phCountMismatch, labrEventCounts.phCountMismatch);
        row("Skipped: RAW flags", hpgeEventCounts.skippedPHRawFlags, labrEventCounts.skippedPHRawFlags);
        row("Skipped: invalid RAW header", hpgeEventCounts.skippedPHInvalidHeader, labrEventCounts.skippedPHInvalidHeader);
        row("Skipped: no preceding RAW", hpgeEventCounts.skippedPHNoRaw, labrEventCounts.skippedPHNoRaw);
        row("RAW sets skipped: flags", hpgeEventCounts.setsSkippedRawFlags, labrEventCounts.setsSkippedRawFlags);
        row("RAW sets skipped: invalid header", hpgeEventCounts.setsSkippedInvalidHeaders, labrEventCounts.setsSkippedInvalidHeaders);
        row("RAW sets skipped: invalid anchors", hpgeEventCounts.setsSkippedInvalidAnchors, labrEventCounts.setsSkippedInvalidAnchors);
        mf::LogInfo("STMDigisFromFragments") << out.str();
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

void STMDigisFromFragments::printJobSummary() const {
    std::ostringstream out;
    auto row = [&](const char* det, const char* type, const fragCounts& c) {
        out << std::left << std::setw(12) << det << std::setw(8) << type
            << std::right << std::setw(12) << c.seen << std::setw(12) << c.good
            << std::setw(12) << (std::string(type) == "PH" ? "-" : std::to_string(c.prescaled))
            << std::setw(12) << (std::string(type) == "PH" ? std::to_string(c.noHits) : "-")
            << std::setw(12) << (std::string(type) == "ZS" ? std::to_string(c.noRegions) : "-") << '\n';
    };
    out << "\n======= STM Job Summary =======\n"
        << "Art events processed: " << _totalEvents << '\n'
        << "Inner fragments: " << _totalInnerFrags << "\n\n"
        << std::left << std::setw(12) << "Detector" << std::setw(8) << "Type"
        << std::right << std::setw(12) << "Seen" << std::setw(12) << "Good"
        << std::setw(12) << "Prescaled" << std::setw(12) << "No hits"
        << std::setw(12) << "No regions" << '\n'
        << std::string(80, '-') << '\n';
    row("HPGe", "RAW", hpgeJobCounts.raw);
    row("HPGe", "ZS", hpgeJobCounts.zs);
    row("HPGe", "PH", hpgeJobCounts.ph);
    row("LaBr", "RAW", labrJobCounts.raw);
    row("LaBr", "ZS", labrJobCounts.zs);
    row("LaBr", "PH", labrJobCounts.ph);
    out << "\nProducts Saved (job total):\n"
        << "HPGe | rawWithHeader : " << hpgeJobProducts.rawWithHeader
        << " , raw : " << hpgeJobProducts.rawWaveforms
        << " , ZS : " << hpgeJobProducts.zsWaveforms
        << " , PH : " << hpgeJobProducts.phDigis << '\n'
        << "LaBr | rawWithHeader : " << labrJobProducts.rawWithHeader
        << " , raw : " << labrJobProducts.rawWaveforms
        << " , ZS : " << labrJobProducts.zsWaveforms
        << " , PH : " << labrJobProducts.phDigis << '\n';
    mf::LogInfo("STMDigisFromFragments") << out.str();
}

void STMDigisFromFragments::printDetailedJobSummary() const {
    std::ostringstream out;
    auto row = [&](const std::string& name, size_t hpge, size_t labr) {
        out << std::left << std::setw(36) << name
            << std::right << std::setw(12) << hpge << std::setw(12) << labr
            << std::setw(12) << hpge + labr << '\n';
    };
    auto header = [&](const std::string& name) {
        out << std::left << std::setw(36) << name
            << std::right << std::setw(12) << "HPGe"
            << std::setw(12) << "LaBr"
            << std::setw(12) << "Total" << '\n'
            << std::string(72, '-') << '\n';
    };
    out << "\n======= STM Detailed Job Summary =======\n"
        << "Events with both detectors : " << _totalEventsWithBothDetectors << '\n'
        << "Only HPGe                  : " << _totalEventsWithOnlyHPGe << '\n'
        << "Only LaBr                  : " << _totalEventsWithOnlyLaBr << '\n'
        << "Neither                    : " << _totalEventsWithNeitherDetector << '\n'
        << "Recognized STM containers  : " << _totalContainers << '\n'
        << "Unknown containers         : " << _totalUnknownContainers << '\n'
        << "Non-container fragments    : " << _totalNonContainers << '\n'
        << "Unknown inner fragments    : " << _totalUnknownInnerFrags << "\n\n";

    header("Frag Type");
    row("Containers", hpgeJobCounts.container, labrJobCounts.container);
    row("Inner fragments", hpgeJobCounts.innerFrags, labrJobCounts.innerFrags);
    out << "\n";

    header("Payload Classification");
    row("RAW empty", hpgeJobCounts.raw.empty, labrJobCounts.raw.empty);
    row("ZS empty", hpgeJobCounts.zs.empty, labrJobCounts.zs.empty);
    row("PH empty", hpgeJobCounts.ph.empty, labrJobCounts.ph.empty);
    row("RAW zero-filled", hpgeJobCounts.raw.zero, labrJobCounts.raw.zero);
    row("ZS zero-filled", hpgeJobCounts.zs.zero, labrJobCounts.zs.zero);
    row("PH zero-filled", hpgeJobCounts.ph.zero, labrJobCounts.ph.zero);
    row("RAW unreadable", hpgeJobCounts.raw.unread, labrJobCounts.raw.unread);
    row("ZS unreadable", hpgeJobCounts.zs.unread, labrJobCounts.zs.unread);
    row("PH unreadable", hpgeJobCounts.ph.unread, labrJobCounts.ph.unread);
    out << '\n';

    header("Raw Diagnostic");
    row("Invalid RAW headers", hpgeJobCounts.invalidHeaders, labrJobCounts.invalidHeaders);
    row("Invalid RAW anchors", hpgeJobCounts.invalidAnchors, labrJobCounts.invalidAnchors);
    row("RAW flagged bad only", hpgeJobCounts.badOnly, labrJobCounts.badOnly);
    row("RAW flagged missing only", hpgeJobCounts.missingOnly, labrJobCounts.missingOnly);
    row("RAW flagged bad and missing", hpgeJobCounts.badAndMissing, labrJobCounts.badAndMissing);
    out << '\n';

    header("ZS Diagnostic");
    row("ZS header mismatch", hpgeJobCounts.zsHeaderMismatch, labrJobCounts.zsHeaderMismatch);
    row("ZS length mismatch", hpgeJobCounts.zsLengthMismatch, labrJobCounts.zsLengthMismatch);
    row("ZS region mismatch", hpgeJobCounts.zsRegionMismatch, labrJobCounts.zsRegionMismatch);
    row("ZS skipped: RAW flags", hpgeJobCounts.skippedZSRawFlags, labrJobCounts.skippedZSRawFlags);
    row("ZS skipped: invalid RAW header", hpgeJobCounts.skippedZSInvalidHeader, labrJobCounts.skippedZSInvalidHeader);
    row("ZS skipped: no preceding RAW", hpgeJobCounts.skippedZSNoRaw, labrJobCounts.skippedZSNoRaw);
    out << '\n';

    header("PH Diagnostic");
    row("PH count mismatch", hpgeJobCounts.phCountMismatch, labrJobCounts.phCountMismatch);
    row("PH skipped: RAW flags", hpgeJobCounts.skippedPHRawFlags, labrJobCounts.skippedPHRawFlags);
    row("PH skipped: invalid RAW header", hpgeJobCounts.skippedPHInvalidHeader, labrJobCounts.skippedPHInvalidHeader);
    row("PH skipped: no preceding RAW", hpgeJobCounts.skippedPHNoRaw, labrJobCounts.skippedPHNoRaw);
    out << '\n';

    header("Raw|ZS|PH");
    row("RAW sets skipped: flags", hpgeJobCounts.setsSkippedRawFlags, labrJobCounts.setsSkippedRawFlags);
    row("RAW sets skipped: short/invalid", hpgeJobCounts.setsSkippedInvalidHeaders, labrJobCounts.setsSkippedInvalidHeaders);
    row("RAW sets skipped: anchors", hpgeJobCounts.setsSkippedInvalidAnchors, labrJobCounts.setsSkippedInvalidAnchors);
    mf::LogInfo("STMDigisFromFragments") << out.str();
}

void STMDigisFromFragments::endJob() {
    if (_verbosityLevel >= 1) {
        // Print job Summary
        std::ostringstream out;
        out << "\n========== STM JOB SUMMARY - (Unpacking Module) ==========\n"
        << "\n--- Module Configuration For Job ---\n"
        << "Raw Waveforms with Header HPGe   : " << (_saveRawWaveformsWithHeaderHPGe ? "Yes" : "No") << "\n"
        << "Raw Waveforms HPGe               : " << (_saveRawWaveformsHPGe ? "Yes" : "No") << "\n"
        << "ZS Waveforms HPGe                : " << (_saveZSWaveformsHPGe ? "Yes" : "No") << "\n"
        << "PH Digis HPGe                    : Yes\n"
        << "Raw Waveforms with Header LaBr   : " << (_saveRawWaveformsWithHeaderLaBr ? "Yes" : "No") << "\n"
        << "Raw Waveforms LaBr               : " << (_saveRawWaveformsLaBr ? "Yes" : "No") << "\n"
        << "ZS Waveforms LaBr                : " << (_saveZSWaveformsLaBr ? "Yes" : "No") << "\n"
        << "PH Digis LaBr                    : Yes" << "\n"
        << "Frag Summaries                   : " << (_saveSTMFragSummary ?"Yes" : "No" ) << "\n\n"
        << "--- RAW Header Versions Seen ---\n";

        if (_rawHeaderVersions.empty()) {
            out << "No readable format observed\n";
        } else {
            for (const auto& entry : _rawHeaderVersions) {
                out << "Version " << entry.first;
                if (entry.first == stm::RawHeader::VERSION_LEGACY) {
                    out << " (legacy)";
                }
                out << ": " << entry.second << " Raw fragments\n";
            }
        }

        mf::LogInfo("STMDigisFromFragments") << out.str();
        printJobSummary();
    }
    if (_verbosityLevel >= 2) {printDetailedJobSummary();}
}

DEFINE_ART_MODULE(STMDigisFromFragments)
