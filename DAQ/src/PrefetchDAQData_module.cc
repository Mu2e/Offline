// Prefetch DAQ data to avoid fetch times impacting module timing evaluations

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// data
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

// DAQ data
#include <artdaq-core/Data/Fragment.hh>


namespace mu2e {

  class PrefetchDAQData : public art::EDProducer {

  protected:

  public:
    struct Config {
      using  Name    = fhicl::Name;
      using  Comment = fhicl::Comment;

      fhicl::Atom<int>         debugLevel            {Name("debugLevel"            ), Comment("Debug output verbosity level"), 0};
      fhicl::Atom<bool>        fetchCaloDigis        {Name("fetchCaloDigis"        ), Comment("Prefetch calorimeter digi collections")};
      fhicl::Atom<bool>        fetchStrawDigis       {Name("fetchStrawDigis"       ), Comment("Prefetch straw tracker digi collections")};
      fhicl::Atom<bool>        fetchCaloFragments    {Name("fetchCaloFragments"    ), Comment("Prefetch calorimeter DAQ fragments")};
      fhicl::Atom<bool>        fetchTrkFragments     {Name("fetchTrkFragments"     ), Comment("Prefetch tracker DAQ fragments")};
      fhicl::Atom<bool>        fetchAllFragments     {Name("fetchAllFragments"     ), Comment("Prefetch all available DAQ fragments regardless of subdetector-specific fetch flags" )};
      fhicl::Atom<bool>        fetchEWM              {Name("fetchEWM"              ), Comment("Prefetch EventWindowMarker")};
      fhicl::Atom<bool>        fetchPBT              {Name("fetchPBT"              ), Comment("Prefetch ProtonBunchTime")};
      fhicl::Atom<std::string> caloDigiCollectionTag {Name("caloDigiCollectionTag" ), Comment("Input tag for the calorimeter digi collection")};
      fhicl::Atom<std::string> strawDigiCollectionTag{Name("strawDigiCollectionTag"), Comment("Input tag for the straw tracker digi collection")};
      fhicl::Atom<std::string> caloFragmentTag       {Name("caloFragmentTag"       ), Comment("Input tag for the calorimeter DAQ fragment collection")};
      fhicl::Atom<std::string> trkFragmentTag        {Name("trkFragmentTag"        ), Comment("Input tag for the tracker DAQ fragment collection")};
      fhicl::Atom<std::string> pbtTag                {Name("pbtTag"                ), Comment("Input tag for the ProtonBunchTime")};
      fhicl::Atom<std::string> ewmTag                {Name("ewmTag"                ), Comment("Input tag for the EventWindowMarker")};
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit PrefetchDAQData(const Parameters&);
    virtual ~PrefetchDAQData();
    virtual void beginJob();
    virtual void endJob();
    virtual void produce( art::Event& e);

    void   fake_access(const CaloDigi&    Digi);
    void   fake_access(const StrawDigi&   Digi);
    void   fake_access(const artdaq::Fragment&   Frag);

  private:

    bool findData(const art::Event& e);
    template <class T> const T* fetch(const art::Event& e, const art::InputTag& tag, unsigned long& nMissing);
                                        // control flags
    int   _debugLevel;
    bool  _fetchCaloDigis;
    bool  _fetchStrawDigis;
    bool  _fetchCaloFragments;
    bool  _fetchTrkFragments;
    bool  _fetchAllFragments;
    bool  _fetchEWM;
    bool  _fetchPBT;

    art::InputTag                  _cdTag;
    art::InputTag                  _sdTag;
    art::InputTag                  _cfTag;
    art::InputTag                  _tfTag;
    art::InputTag                  _ewmTag;
    art::InputTag                  _pbtTag;
                                        // cache of event objects
    const CaloDigiCollection*      _cdcol;
    const StrawDigiCollection*     _sdcol;
    const artdaq::Fragments*       _cfcol;
    const artdaq::Fragments*       _tfcol;
    const EventWindowMarker*       _ewm;
    const ProtonBunchTime*         _pbt;
    int                            _eventNum;
                                        // missing product counters
    unsigned long                  _nEvents;
    unsigned long                  _nMissingCD;
    unsigned long                  _nMissingSD;
    unsigned long                  _nMissingCF;
    unsigned long                  _nMissingTF;
    unsigned long                  _nMissingEWM;
    unsigned long                  _nMissingPBT;
  };

//-----------------------------------------------------------------------------
  PrefetchDAQData::PrefetchDAQData(const Parameters& config):
    art::EDProducer(config),
    _debugLevel        (config().debugLevel()),
    _fetchCaloDigis    (config().fetchCaloDigis()),
    _fetchStrawDigis   (config().fetchStrawDigis()),
    _fetchCaloFragments(config().fetchCaloFragments()),
    _fetchTrkFragments (config().fetchTrkFragments()),
    _fetchAllFragments (config().fetchAllFragments()),
    _fetchEWM          (config().fetchEWM()),
    _fetchPBT          (config().fetchPBT()),
    _cdTag             (config().caloDigiCollectionTag()),
    _sdTag             (config().strawDigiCollectionTag()),
    _cfTag             (config().caloFragmentTag()),
    _tfTag             (config().trkFragmentTag()),
    _ewmTag            (config().ewmTag()),
    _pbtTag            (config().pbtTag()),
    _nEvents(0), _nMissingCD(0), _nMissingSD(0), _nMissingCF(0), _nMissingTF(0), _nMissingEWM(0), _nMissingPBT(0)
  {}

  PrefetchDAQData::~PrefetchDAQData() {
  }

//-----------------------------------------------------------------------------
  void PrefetchDAQData::beginJob() {
  }

//-----------------------------------------------------------------------------
  void PrefetchDAQData::endJob() {
    const unsigned long nMissing = _nMissingCD + _nMissingSD + _nMissingCF + _nMissingTF + _nMissingEWM + _nMissingPBT;
    if (nMissing > 0) {
      mf::LogWarning("PrefetchDAQData") << "Missing products in " << _nEvents << " events:"
                                        << " CaloDigis: "      << _nMissingCD
                                        << " StrawDigis: "     << _nMissingSD
                                        << " calo fragments: " << _nMissingCF
                                        << " trk fragments: "  << _nMissingTF
                                        << " EWM: "            << _nMissingEWM
                                        << " PBT: "            << _nMissingPBT;
    }
  }

//-----------------------------------------------------------------------------
// Retrieve a product if present, counting and reporting when it is not. Repeated
// warnings are rate-limited by the message facility configuration for this category
//-----------------------------------------------------------------------------
  template <class T>
  const T* PrefetchDAQData::fetch(const art::Event& evt, const art::InputTag& tag, unsigned long& nMissing) {
    art::Handle<T> handle;
    if (evt.getByLabel(tag, handle) && handle.isValid()) return handle.product();
    ++nMissing;
    mf::LogWarning("PrefetchDAQData") << "Event " << evt.id() << ": product " << tag.encode()
                                      << " not found, skipping it";
    return nullptr;
  }

//-----------------------------------------------------------------------------
  void PrefetchDAQData::fake_access(const CaloDigi& Hit) {
  }

//-----------------------------------------------------------------------------
  void PrefetchDAQData::fake_access(const StrawDigi& X) {
  }
//-----------------------------------------------------------------------------
  void PrefetchDAQData::fake_access(const artdaq::Fragment& Frag) {
  }

//-----------------------------------------------------------------------------
  bool PrefetchDAQData::findData(const art::Event& evt){
    _cdcol   = 0;
    _sdcol   = 0;
    _cfcol   = 0;
    _tfcol   = 0;
    _ewm     = 0;
    _pbt     = 0;

//-----------------------------------------------------------------------------
// Retrieve the data from the event
//-----------------------------------------------------------------------------
    if (_fetchCaloDigis)     _cdcol = fetch<CaloDigiCollection> (evt, _cdTag, _nMissingCD);
    if (_fetchTrkFragments)  _tfcol = fetch<artdaq::Fragments>  (evt, _tfTag, _nMissingTF);
    if (_fetchCaloFragments) _cfcol = fetch<artdaq::Fragments>  (evt, _cfTag, _nMissingCF);
    if (_fetchStrawDigis)    _sdcol = fetch<StrawDigiCollection>(evt, _sdTag, _nMissingSD);

    if (_fetchAllFragments) {
      auto allH = evt.getMany<std::vector<artdaq::Fragment>>();
      for (auto& fs : allH) {
        for (auto& frag : *fs) {
          fake_access(frag);
        }
      }
    }

    if (_fetchEWM) _ewm = fetch<EventWindowMarker>(evt, _ewmTag, _nMissingEWM);
    if (_fetchPBT) _pbt = fetch<ProtonBunchTime>  (evt, _pbtTag, _nMissingPBT);

//-----------------------------------------------------------------------------
// prefetch data within the collections
//-----------------------------------------------------------------------------
    if (_cdcol) {
      int ncd = _cdcol->size();
      for(int i=0;i<ncd;++i){
        const CaloDigi& cd = _cdcol->at(i);
        fake_access(cd);
      }
    }

    if (_sdcol) {
      int nsd = _sdcol->size();
      for(int i=0;i<nsd;++i){
        const StrawDigi& sdigi = _sdcol->at(i);
        fake_access(sdigi);
      }
    }

    if (_cfcol) {
      int ncf = _cfcol->size();
      for(int i=0;i<ncf;++i){
        auto frag = _cfcol->at(i);
         fake_access(frag);
      }
    }
    if (_tfcol) {
      int ntf = _tfcol->size();
      for(int i=0;i<ntf;++i){
        auto frag = _tfcol->at(i);
         fake_access(frag);
      }
    }

    return true;
  }

//-----------------------------------------------------------------------------
  void PrefetchDAQData::produce(art::Event& Event) {

    _eventNum = Event.event();
    ++_nEvents;

    if (_debugLevel > 0) printf("[PrefetchDAQData::produce] event number: %10i\n",_eventNum);

    findData(Event);
  }

//------------------------------------------------------------------------------
// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(PrefetchDAQData)

}
