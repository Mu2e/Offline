//------------------------------------------------------------
// Make ComboHits from TrackerDataDecoders
//
//------------------------------------------------------------

// framework
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"

// conditions
#include "Offline/ConditionsBase/inc/TrackerCalibrationStructs.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/TrkHitReco/inc/PeakFit.hh"
#include "Offline/TrkHitReco/inc/StrawHitRecoUtils.hh"

#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <memory>

namespace art {
class ComboHitsFromDataDecoders;
}

class art::ComboHitsFromDataDecoders : public art::EDProducer {
public:
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  using ADCWFIter = mu2e::TrkTypes::ADCWaveform::const_iterator;
  struct Config {
    fhicl::Atom<int> diag{Name("diagLevel"), Comment("Diag level")};
    fhicl::Atom<int> print{Name("printLevel"), Comment("Print level")};
    fhicl::Atom<int> fittype{Name("FitType"), Comment("Waveform Fit Type")};
    fhicl::Atom<bool> usecc{Name("UseCalorimeter"), Comment("Use Calo cluster times to filter")};
    fhicl::Atom<float> clusterDt{Name("clusterDt"), Comment("Calo cluster time 1/2 window")};
    fhicl::Atom<float> minE{Name("minimumEnergy"), Comment("Minimum straw energy deposit (MeV)")};
    fhicl::Atom<float> maxE{Name("maximumEnergy"), Comment("Maximum straw energy deposit (MeV)")};
    fhicl::Atom<float> minR{Name("minimumRho"), Comment("Minimum transverse radius (mm)")};
    fhicl::Atom<float> maxR{Name("maximumRho"), Comment("Maximum transverse radius (mm)")};
    fhicl::Atom<float> ctE{Name("crossTalkEnergy"),
                           Comment("Energy to filter cross-talk in adajcent straws (MeV)")};
    fhicl::Atom<float> ctMinT{Name("crossTalkMinimumTime"),
                              Comment("Earliest time for cross-talk filter (nsec)")};
    fhicl::Atom<float> ctMaxT{Name("crossTalkMaximumTime"),
                              Comment("Latest time for cross-talk filter (nsec)")};
    fhicl::OptionalAtom<float> minTOff{Name("minimumTimeOffSpill"),
                            Comment("Earliest StrawDigi time to process (nsec)")};
    fhicl::OptionalAtom<float> maxTOff{Name("maximumTimeOffSpill"),
                            Comment("Latest StrawDigi time to process (nsec)")};
    fhicl::Atom<float> minT{Name("minimumTime"),
                            Comment("Earliest StrawDigi time to process (nsec)")};
    fhicl::Atom<float> maxT{Name("maximumTime"),
                            Comment("Latest StrawDigi time to process (nsec)")};
    fhicl::Atom<bool> filter{Name("FilterHits"),
                             Comment("Filter hits (alternative is to just flag)")};
    fhicl::Atom<bool> writesh{Name("WriteStrawHitCollection"), Comment("Save StrawHitCollection")};
    fhicl::Atom<bool> flagXT{Name("FlagCrossTalk"), Comment("Search for cross-talk")};
    fhicl::Atom<art::InputTag> cccTag{Name("CaloClusterCollectionTag"),
                                      Comment("CaloClusterCollection producer")};
    fhicl::Atom<art::InputTag> pbttoken{Name("ProtonBunchTimeTag"),
                                        Comment("ProtonBunchTime producer")};
    fhicl::Atom<art::InputTag> tfTag{Name("TrackerDataDecoderTag"),
                                     Comment("StrawDigiCollection producer")};
  };

  using Parameters = art::EDProducer::Table<Config>;
  explicit ComboHitsFromDataDecoders(Parameters const& config);
  void produce(art::Event& e) override;
  void beginJob() override;

private:
  float _minTOff, _maxTOff;
  bool _overrideminTOff, _overridemaxTOff;
  mu2e::StrawHitRecoUtils _shrUtils;
  bool _writesh; // write straw hits or not
  bool _flagXT;  // flag cross-talk
  bool _usecc;
  bool _useADCWF;
  int _printLevel;
  int _diagLevel;

  art::ProductToken<mu2e::CaloClusterCollection> const _ccctoken;
  art::ProductToken<mu2e::ProtonBunchTime> const
      _pbttoken; // name of the module that makes eventwindowmarkers
  art::InputTag _tfTag;
  std::unique_ptr<mu2e::TrkHitReco::PeakFit> _pfit; // peak fitting algorithm
  //
  // helper function
  //
  void analyze_tracker_(const mu2e::TrackerDataDecoder& cc,
                        std::unique_ptr<mu2e::StrawHitCollection> const& shCol,
                        std::unique_ptr<mu2e::ComboHitCollection> const& chCol, double pbtOffset,
                        mu2e::TrackerStatus const& trackerStatus, mu2e::StrawResponse const& srep,
                        const mu2e::CaloClusterCollection* caloClusters, mu2e::Tracker const& tt);
  mu2e::ProditionsHandle<mu2e::StrawResponse> _strawResponse_h;
  mu2e::ProditionsHandle<mu2e::TrackerStatus> _trackerStatus_h;
  mu2e::ProditionsHandle<mu2e::Tracker> _alignedTracker_h;
};

art::ComboHitsFromDataDecoders::ComboHitsFromDataDecoders(Parameters const& config) :
    art::EDProducer{config},
    _minTOff(config().minT()),
    _maxTOff(config().maxT()),
    _overrideminTOff(config().minTOff(_minTOff)),
    _overridemaxTOff(config().maxTOff(_maxTOff)),
    _shrUtils((mu2e::TrkHitReco::FitType)config().fittype(), config().diag(),
              mu2e::StrawIdMask::uniquestraw, // this module produces individual straw ComboHits
              config().writesh(), config().minT(), config().maxT(),
              _overrideminTOff, _minTOff, _overridemaxTOff, _maxTOff, config().minE(),
              config().maxE(), config().minR(), config().maxR(), config().filter(), config().ctE(),
              config().ctMinT(), config().ctMaxT(), config().usecc(), config().clusterDt()),
    _writesh(config().writesh()), _flagXT(config().flagXT()), _usecc(config().usecc()),
    _useADCWF(config().fittype() != mu2e::TrkHitReco::FitType::firmwarepmp),
    _printLevel(config().print()),
    _diagLevel(config().diag()), _ccctoken{mayConsume<mu2e::CaloClusterCollection>(
                                     config().cccTag())},
    _pbttoken{consumes<mu2e::ProtonBunchTime>(config().pbttoken())}, _tfTag(config().tfTag()) {
  produces<mu2e::ComboHitCollection>();
  if (_writesh)
    produces<mu2e::StrawHitCollection>();
  if (_printLevel > 0)
    std::cout << "In ComboHitsFromDataDecoders constructor " << std::endl;
}

//------------------------------------------------------------------------------------------
void art::ComboHitsFromDataDecoders::beginJob() {}

//------------------------------------------------------------------------------------------
void art::ComboHitsFromDataDecoders::produce(art::Event& event) {
  if (_printLevel > 0)
    std::cout << "In ComboHitsFromDataDecoders produce " << std::endl;
  const mu2e::Tracker& tt = _alignedTracker_h.get(event.id());
  auto const& srep = _strawResponse_h.get(event.id());
  //_tfTag = art::InputTag("test");

  size_t numTrkFrags = 0;
  auto fragmentHandle = event.getValidHandle<std::vector<mu2e::TrackerDataDecoder> >(_tfTag);

  for (auto frag : *fragmentHandle) {
    numTrkFrags++;
  }

  const mu2e::CaloClusterCollection* caloClusters(0);
  if (_usecc) {
    auto ccH = event.getValidHandle(_ccctoken);
    caloClusters = ccH.product();
  }

  art::Handle<mu2e::ProtonBunchTime> pbtHandle;
  auto pbtH = event.getValidHandle(_pbttoken);
  const mu2e::ProtonBunchTime& pbt(*pbtH);
  double pbtOffset = pbt.pbtime_;

  mu2e::TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());
  std::unique_ptr<mu2e::StrawHitCollection> shCol;
  if (_writesh) {
    shCol = std::unique_ptr<mu2e::StrawHitCollection>(new mu2e::StrawHitCollection);
  }
  std::unique_ptr<mu2e::ComboHitCollection> chCol(new mu2e::ComboHitCollection());

  if (numTrkFrags == 0) {
    std::cout << "[StrawAndCaloDigisFromFragments::produce] found no Tracker fragments!"
              << std::endl;
    if (_writesh)
      event.put(std::move(shCol));
    event.put(std::move(chCol));
    return;
  }

  if (_writesh) {
    shCol->reserve(numTrkFrags);
  }
  chCol->reserve(numTrkFrags);

  for (auto frag : *fragmentHandle) {
    analyze_tracker_(frag, shCol, chCol, pbtOffset, trackerStatus, srep, caloClusters, tt);
  }

  if (_writesh)
    event.put(std::move(shCol));
  event.put(std::move(chCol));
}

void art::ComboHitsFromDataDecoders::analyze_tracker_(
    const mu2e::TrackerDataDecoder& cc, std::unique_ptr<mu2e::StrawHitCollection> const& shCol,
    std::unique_ptr<mu2e::ComboHitCollection> const& chCol, double pbtOffset,
    mu2e::TrackerStatus const& trackerStatus, mu2e::StrawResponse const& srep,
    const mu2e::CaloClusterCollection* caloClusters, mu2e::Tracker const& tt) {

  if (_diagLevel > 1) {
    std::cout << std::endl;
    std::cout << "TrackerDataDecoder: ";
    std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
    std::cout << std::endl;
    std::cout << "\t"
              << "====== Example Block Sizes ======" << std::endl;
    for (size_t i = 0; i < 10; i++) {
      if (i < cc.block_count()) {
        std::cout << "\t" << i << "\t" << cc.blockSizeBytes(i) << std::endl;
      }
    }
    std::cout << "\t"
              << "=========================" << std::endl;
  }

  for (size_t curBlockIdx = 0; curBlockIdx < cc.block_count(); curBlockIdx++) {
    auto block = cc.dataAtBlockIndex(curBlockIdx);
    if (block == nullptr) {
      mf::LogError("StrawAndCaloDigisFromFragments")
          << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (_diagLevel > 1) {

      std::cout << "timestamp: "
                << static_cast<int>(hdr->GetEventWindowTag().GetEventWindowTag(true)) << std::endl;
      std::cout << "hdr->SubsystemID: " << static_cast<int>(hdr->GetSubsystemID()) << std::endl;
      std::cout << "dtcID: " << static_cast<int>(hdr->GetID()) << std::endl;
      std::cout << "rocID: " << static_cast<int>(hdr->GetLinkID()) << std::endl;
      std::cout << "packetCount: " << static_cast<int>(hdr->GetPacketCount()) << std::endl;
      std::cout << "EVB mode: " << static_cast<int>(hdr->GetEVBMode()) << std::endl;

      std::cout << std::endl;
    }

    // Parse phyiscs information from TRK packets
    if (hdr->GetPacketCount() > 0) {

      // Create the StrawDigi data products
      auto trkDataVec = cc.GetTrackerData(curBlockIdx, _useADCWF);
      if (trkDataVec.empty()) {
        mf::LogError("StrawAndCaloDigisFromFragments")
            << "Error retrieving Tracker data from DataBlock " << curBlockIdx
            << "! Aborting processing of this block!";
        continue;
      }

      double pmp(0.0);
      for (auto& trkDataPair : trkDataVec) {
        mu2e::StrawId sid(trkDataPair.first->StrawIndex);
        mu2e::TrkTypes::TDCValues tdc = {trkDataPair.first->TDC0(), trkDataPair.first->TDC1()};
        mu2e::TrkTypes::TOTValues tot = {trkDataPair.first->TOT0, trkDataPair.first->TOT1};
        if (!_useADCWF) {
          pmp = trkDataPair.first->PMP;
        } else {
          ADCWFIter maxiter;
          pmp = _shrUtils.peakMinusPedWF(trkDataPair.second, srep, maxiter);
        }
        // temporary hack, FIXME
        mu2e::EventWindowMarker ewm;
        ewm._spillType = mu2e::EventWindowMarker::onspill;
        ewm._eventLength = 1695.0;
        _shrUtils.createComboHit(ewm,-1,chCol, shCol, caloClusters, pbtOffset, sid, tdc, tot, pmp, trackerStatus,  srep, tt);
      }

      // flag straw and electronic cross-talk
      if (_flagXT) {
        _shrUtils.flagCrossTalk(shCol, chCol);
      }
    }
  }
  // cc.ClearUpgradedPackets();
}

DEFINE_ART_MODULE(art::ComboHitsFromDataDecoders)
