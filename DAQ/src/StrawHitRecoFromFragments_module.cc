//
// This module transforms StrawDigi objects into StrawHit objects
//
// Original author David Brown, LBNL
// Merged with flag and position creation B. Echenard, CalTech
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"

// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsBase/inc/TrackerCalibrationStructs.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

#include "Offline/TrkHitReco/inc/PeakFit.hh"
#include "Offline/TrkHitReco/inc/PeakFitRoot.hh"
#include "Offline/TrkHitReco/inc/PeakFitFunction.hh"
#include "Offline/TrkHitReco/inc/ComboPeakFitRoot.hh"
#include "Offline/TrkHitReco/inc/StrawHitRecoUtils.hh"

#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "art/Framework/Principal/Handle.h"
#include "mu2e-artdaq-core/Overlays/CalorimeterFragment.hh"
#include "mu2e-artdaq-core/Overlays/FragmentType.hh"
#include "mu2e-artdaq-core/Overlays/TrackerFragment.hh"
#include "mu2e-artdaq-core/Overlays/Mu2eEventFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "TH1F.h"

#include <memory>


namespace art {
class StrawHitRecoFromFragments;
}

class art::StrawHitRecoFromFragments : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag level")};
      fhicl::Atom<int> print{ Name("printLevel"), Comment("Print level")};
      fhicl::Atom<int> fittype { Name( "FitType"), Comment("Waveform Fit Type")};//, TrkHitReco::FitType::firmwarepmp;
      fhicl::Atom<bool> usecc{ Name("UseCalorimeter"), Comment("Use Calo cluster times to filter" )};
      fhicl::Atom<float>clusterDt{ Name("clusterDt"), Comment("Calo cluster time 1/2 window")};
      fhicl::Atom<float>minE{ Name("minimumEnergy"), Comment("Minimum straw energy deposit (MeV)")};
      fhicl::Atom<float>maxE{ Name("maximumEnergy"), Comment("Maximum straw energy deposit (MeV)")};
      fhicl::Atom<float>ctE{ Name("crossTalkEnergy"), Comment("Energy to filter cross-talk in adajcent straws (MeV)")};
      fhicl::Atom<float>ctMinT{ Name("crossTalkMinimumTime"), Comment("Earliest time for cross-talk filter (nsec)")};
      fhicl::Atom<float>ctMaxT{ Name("crossTalkMaximumTime"), Comment("Latest time for cross-talk filter (nsec)")};
      fhicl::Atom<float>minT{ Name("minimumTime"), Comment("Earliest StrawDigi time to process (nsec)")};
      fhicl::Atom<float>maxT{ Name("maximumTime"), Comment("Latest StrawDigi time to process (nsec)")};
      fhicl::Atom<bool>filter{ Name("FilterHits"), Comment("Filter hits (alternative is to just flag)") };
      fhicl::Atom<bool>writesh{ Name("WriteStrawHitCollection"), Comment("Save StrawHitCollection")};
      fhicl::Atom<bool>flagXT{ Name("FlagCrossTalk"), Comment("Search for cross-talk")};
      fhicl::Atom<art::InputTag> cccTag{ Name("CaloClusterCollectionTag"), Comment("CaloClusterCollection producer")};
      fhicl::Atom<art::InputTag> pbttoken{ Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer")};
      fhicl::Atom<art::InputTag> tfTag{ Name("TrackerFragmentTag"), Comment("StrawDigiCollection producer")};
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit StrawHitRecoFromFragments(Parameters const& config);
    void produce( art::Event& e) override;
    void beginRun( art::Run& run ) override;
    void beginJob() override;


  private:
    mu2e::TrkHitReco::FitType _fittype; // peak Fitter
    bool  _usecc;                   // use calorimeter cluster filtering
    float _clusterDt;               // maximum hit-calo lcuster time difference
    float _minE;             // energy range (MeV)
    float _maxE;             // energy range (MeV)
    float _ctE;                     // minimum charge to flag neighbors as cross talk
    float _ctMinT;                  // time relative to proton hit to flag cross talk (ns)
    float _ctMaxT;                  // time relative to proton hit to flag cross talk (ns)
    float _minT, _maxT;             // time range
    bool  _filter;                // filter the output, or just flag
    bool  _writesh;                // write straw hits or not
    bool  _flagXT; // flag cross-talk
    int   _printLevel;
    int   _diagLevel;
    mu2e::StrawIdMask _mask;
    mu2e::StrawEnd _end[2]; // helper
    float _invnpre; // cache
    float _invgainAvg; // cache
    float _invgain[96]; // cache
    unsigned _npre; //cache
    double pbtOffset;

    art::ProductToken<mu2e::CaloClusterCollection> const _ccctoken;
    art::ProductToken<mu2e::ProtonBunchTime> const _pbttoken; // name of the module that makes eventwindowmarkers
    art::InputTag _tfTag;
    std::unique_ptr<mu2e::TrkHitReco::PeakFit> _pfit; // peak fitting algorithm
    size_t npanels, nplanes;
    // diagnostic
    TH1F* _maxiter;
    // helper function
    //
    void analyze_tracker_(const mu2e::TrackerFragment& cc, std::unique_ptr<mu2e::StrawHitCollection> const& shCol,
        std::unique_ptr<mu2e::ComboHitCollection> const& chCol, mu2e::StrawHitRecoUtils &shrUtils, mu2e::TrackerStatus const& trackerStatus, mu2e::StrawResponse const& srep,
        const mu2e::CaloClusterCollection* caloClusters, mu2e::Tracker const& tt
        );
    float peakMinusPedAvg(mu2e::TrkTypes::ADCWaveform const& adcData) const;
    float peakMinusPed(mu2e::StrawId id, mu2e::TrkTypes::ADCWaveform const& adcData) const;
    float peakMinusPedFirmware(mu2e::StrawId id, mu2e::TrkTypes::ADCValue const& pmp) const;
    mu2e::ProditionsHandle<mu2e::StrawResponse> _strawResponse_h;
    mu2e::ProditionsHandle<mu2e::TrackerStatus> _trackerStatus_h;
    mu2e::ProditionsHandle<mu2e::Tracker> _alignedTracker_h;
};

art::StrawHitRecoFromFragments::StrawHitRecoFromFragments(Parameters const& config) :
  art::EDProducer{config},
  _fittype((mu2e::TrkHitReco::FitType) config().fittype()),
  _usecc(config().usecc()),
  _clusterDt(config().clusterDt()),
  _minE(config().minE()),
  _maxE(config().maxE()),
  _ctE(config().ctE()),
  _ctMinT(config().ctMinT()),
  _ctMaxT(config().ctMaxT()),
  _minT(config().minT()),
  _maxT(config().maxT()),
  _filter(config().filter()),
  _writesh(config().writesh()),
  _flagXT(config().flagXT()),
  _printLevel(config().print()),
  _diagLevel(config().diag()),
  _mask(mu2e::StrawIdMask::uniquestraw), // this module produces individual straw ComboHits
  _end{mu2e::StrawEnd::cal,mu2e::StrawEnd::hv}, // this should be in a general place, FIXME!
  _ccctoken{mayConsume<mu2e::CaloClusterCollection>(config().cccTag())},
  _pbttoken{consumes<mu2e::ProtonBunchTime>(config().pbttoken())},
  _tfTag(config().tfTag())
{
  produces<mu2e::ComboHitCollection>();
  if(_writesh)produces<mu2e::StrawHitCollection>();
  if (_printLevel > 0) std::cout << "In StrawHitRecoFromFragments constructor " << std::endl;
}

//------------------------------------------------------------------------------------------
void art::StrawHitRecoFromFragments::beginJob()
{
  if(_diagLevel > 0){
    art::ServiceHandle<art::TFileService> tfs;
    _maxiter   = tfs->make<TH1F>( "maxiter",  "ADC max",16,-0.5,15.5 );
  }
}

void art::StrawHitRecoFromFragments::beginRun(art::Run& run)
{
  auto const& srep = _strawResponse_h.get(run.id());
  // set cache for peak-ped calculation (default)
  _npre = srep.nADCPreSamples();
  _invnpre = 1.0/(float)_npre;
  _invgainAvg = srep.adcLSB()*srep.peakMinusPedestalEnergyScale()/srep.strawGain();
  for (int i=0;i<96;i++){
    mu2e::StrawId dummyId(0,0,i);
    _invgain[i] = srep.adcLSB()*srep.peakMinusPedestalEnergyScale(dummyId)/srep.strawGain();
  }

  // Detailed histogram-based waveform fits are no longer supported TODO!
  if (_fittype != mu2e::TrkHitReco::FitType::peakminusped && _fittype != mu2e::TrkHitReco::FitType::peakminuspedavg && _fittype != mu2e::TrkHitReco::FitType::firmwarepmp)
    throw cet::exception("RECO")<<"TrkHitReco: Peak fit " << _fittype << " not implemented " <<  std::endl;
}

//------------------------------------------------------------------------------------------
void art::StrawHitRecoFromFragments::produce(art::Event& event)
{
  if (_printLevel > 0) std::cout << "In StrawHitRecoFromFragments produce " << std::endl;

  const mu2e::Tracker& tt = _alignedTracker_h.get(event.id());

  nplanes = tt.nPlanes();
  npanels = tt.getPlane(0).nPanels();
  auto const& srep = _strawResponse_h.get(event.id());
  //_tfTag = art::InputTag("test");


  size_t numTrkFrags = 0;
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles =
      event.getMany<std::vector<artdaq::Fragment>>();

  for (const auto& handle : fragmentHandles) {
      if (!handle.isValid() || handle->empty()) {
          continue;
      }

      if (handle->front().type() == mu2e::detail::FragmentType::MU2EEVENT) {
          for (const auto& cont : *handle) {
              mu2e::Mu2eEventFragment mef(cont);
              for (size_t ii = 0; ii < mef.tracker_block_count(); ++ii) {
                  numTrkFrags++;
              }
          }
      }
      else {
          if (handle->front().type() == mu2e::detail::FragmentType::TRK) {
              for (auto frag : *handle) {
                  numTrkFrags++;
              }
          }
      }
  }

  //auto sdH = event.getValidHandle(_sdctoken);
  //const StrawDigiCollection& sdcol(*sdH);

  //const StrawDigiADCWaveformCollection *sdadccol(0);
  //if (_fittype != TrkHitReco::FitType::firmwarepmp) {
  //  auto sdawH = event.getValidHandle(_sdadctoken);
  //  sdadccol = sdawH.product();
  //}

  const mu2e::CaloClusterCollection* caloClusters(0);
  if(_usecc){
    auto ccH = event.getValidHandle(_ccctoken);
    caloClusters = ccH.product();
  }

  art::Handle<mu2e::ProtonBunchTime> pbtHandle;
  auto pbtH = event.getValidHandle(_pbttoken);
  const mu2e::ProtonBunchTime& pbt(*pbtH);
  pbtOffset = pbt.pbtime_;

  mu2e::TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());
  std::unique_ptr<mu2e::StrawHitCollection> shCol;
  if(_writesh){
    shCol = std::unique_ptr<mu2e::StrawHitCollection>(new mu2e::StrawHitCollection);
  }
  std::unique_ptr<mu2e::ComboHitCollection> chCol(new mu2e::ComboHitCollection());

  if (numTrkFrags == 0) {
    std::cout << "[StrawAndCaloDigisFromFragments::produce] found no Tracker fragments!"
      << std::endl;
    if(_writesh)event.put(std::move(shCol));
    event.put(std::move(chCol));
    return;
  }

  if(_writesh){
    shCol->reserve(numTrkFrags);
  }
  chCol->reserve(numTrkFrags);

  mu2e::StrawHitRecoUtils shrUtils(pbtOffset, _fittype, _npre, _invnpre, _invgainAvg, _invgain,
    _diagLevel, _maxiter, _mask, nplanes, npanels, _writesh, _minT, _maxT, _minE, _maxE, _filter, _flagXT,
    _ctE, _ctMinT, _ctMaxT, _usecc, _clusterDt, numTrkFrags);

  for (const auto& handle : fragmentHandles) {
      if (!handle.isValid() || handle->empty()) {
          continue;
      }

      if (handle->front().type() == mu2e::detail::FragmentType::MU2EEVENT) {
          for (const auto& cont : *handle) {
              mu2e::Mu2eEventFragment mef(cont);
              for (size_t ii = 0; ii < mef.tracker_block_count(); ++ii) {
                  auto pair = mef.trackerAtPtr(ii);
                  mu2e::TrackerFragment cc(pair);
                  analyze_tracker_(cc, shCol, chCol, shrUtils, trackerStatus, srep, caloClusters, tt);
              }
          }
      }
      else {
          if (handle->front().type() == mu2e::detail::FragmentType::TRK) {
              for (auto frag : *handle) {
                  mu2e::TrackerFragment cc(frag.dataBegin(), frag.dataSizeBytes());
                  analyze_tracker_(cc, shCol, chCol, shrUtils, trackerStatus, srep, caloClusters, tt);
              }
          }
      }
  }

  if(_writesh)event.put(std::move(shCol));
  event.put(std::move(chCol));
}


void art::StrawHitRecoFromFragments::analyze_tracker_(
    const mu2e::TrackerFragment& cc, std::unique_ptr<mu2e::StrawHitCollection> const& shCol,
    std::unique_ptr<mu2e::ComboHitCollection> const& chCol, mu2e::StrawHitRecoUtils &shrUtils,
    mu2e::TrackerStatus const& trackerStatus, mu2e::StrawResponse const& srep,
    const mu2e::CaloClusterCollection* caloClusters, mu2e::Tracker const& tt
    ) {

  if (_diagLevel > 1) {
    std::cout << std::endl;
    std::cout << "TrackerFragment: ";
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

      std::cout << "timestamp: " << static_cast<int>(hdr->GetEventWindowTag().GetEventWindowTag(true))
        << std::endl;
      std::cout << "hdr->SubsystemID: " << static_cast<int>(hdr->GetSubsystemID()) << std::endl;
      std::cout << "dtcID: " << static_cast<int>(hdr->GetID()) << std::endl;
      std::cout << "rocID: " << static_cast<int>(hdr->GetLinkID()) << std::endl;
      std::cout << "packetCount: " << static_cast<int>(hdr->GetPacketCount()) << std::endl;
      std::cout << "EVB mode: " << static_cast<int>(hdr->GetEVBMode()) << std::endl;

      std::cout << std::endl;
    }

    bool getADC = false;
    if (_fittype != mu2e::TrkHitReco::FitType::firmwarepmp)  getADC = true;

    // Parse phyiscs information from TRK packets
    if (hdr->GetPacketCount() > 0) {

      // Create the StrawDigi data products
      auto trkDataVec = cc.GetTrackerData(curBlockIdx,getADC);
      if (trkDataVec.empty()) {
        mf::LogError("StrawAndCaloDigisFromFragments")
          << "Error retrieving Tracker data from DataBlock " << curBlockIdx
          << "! Aborting processing of this block!";
        continue;
      }


      for (auto& trkDataPair : trkDataVec) {

        mu2e::StrawId sid(trkDataPair.first->StrawIndex);
        mu2e::TrkTypes::TDCValues tdc = {trkDataPair.first->TDC0(), trkDataPair.first->TDC1()};
        mu2e::TrkTypes::TOTValues tot = {trkDataPair.first->TOT0, trkDataPair.first->TOT1};
        mu2e::TrkTypes::ADCValue pmp = trkDataPair.first->PMP;

        shrUtils.createComboHit(-1,chCol, shCol, caloClusters, sid, tdc, tot, pmp, trkDataPair.second,
                trackerStatus,  srep, tt);

      }

      //flag straw and electronic cross-talk
      if(!_filter && _flagXT){
        shrUtils.flagCrossTalk(shCol, chCol);
      }
    }
  }
  //cc.ClearUpgradedPackets();
}


DEFINE_ART_MODULE(art::StrawHitRecoFromFragments);
