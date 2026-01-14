// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

// data
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

// DAQ data
//#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
//#include "artdaq-core-mu2e/Overlays/TrackerDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>


namespace mu2e {

  class PrefetchDAQData : public art::EDProducer {

  protected:

  public:
    struct Config {
      using  Name    = fhicl::Name;
      using  Comment = fhicl::Comment;
      fhicl::Atom<int>         debugLevel            {Name("debugLevel"            ), Comment("debugLevel"            ), 0};
      fhicl::Atom<int>         fetchCaloDigis        {Name("fetchCaloDigis"        ), Comment("fetchCaloDigis"        )};
      fhicl::Atom<int>         fetchStrawDigis       {Name("fetchStrawDigis"       ), Comment("fetchStrawDigis"       )};
      fhicl::Atom<int>         fetchCaloFragments    {Name("fetchCaloFragments"    ), Comment("fetchCaloFragments"    )};
      fhicl::Atom<int>         fetchTrkFragments     {Name("fetchTrkFragments"     ), Comment("fetchTrkFragments"     )};
      fhicl::Atom<int>         fetchAllFragments     {Name("fetchAllFragments"     ), Comment("fetchAllFragments"     )};
      fhicl::Atom<int>         fetchEWM              {Name("fetchEWM"              ), Comment("fetchEWM"              )};
      fhicl::Atom<int>         fetchPBT              {Name("fetchPBT"              ), Comment("fetchPBT"              )};
      fhicl::Atom<std::string> caloDigiCollectionTag {Name("caloDigiCollectionTag" ), Comment("caloDigiCollectionTag" )};
      fhicl::Atom<std::string> strawDigiCollectionTag{Name("strawDigiCollectionTag"), Comment("strawDigiCollectionTag")};
      fhicl::Atom<std::string> caloFragmentTag       {Name("caloFragmentTag"       ), Comment("caloFragmentTag"       )};
      fhicl::Atom<std::string> trkFragmentTag        {Name("trkFragmentTag"        ), Comment("trkFragmentTag"        )};
      fhicl::Atom<std::string> pbtTag                {Name("pbtTag"                ), Comment("pbtTag"                )};
      fhicl::Atom<std::string> ewmTag                {Name("ewmTag"                ), Comment("ewmTag"                )};
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit PrefetchDAQData(const Parameters&);
    virtual ~PrefetchDAQData();
    virtual void beginJob();
    virtual void produce( art::Event& e);

    void   fake_access(const CaloDigi&    Digi);
    void   fake_access(const StrawDigi&   Digi);
    void   fake_access(const artdaq::Fragment&   Frag);

  private:

    bool findData(const art::Event& e);
                                        // control flags
    int   _debugLevel;
    int   _fetchCaloDigis;
    int   _fetchStrawDigis;
    int   _fetchCaloFragments;
    int   _fetchTrkFragments;
    int   _fetchAllFragments;
    int   _fetchEWM;
    int   _fetchPBT;

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
    _pbtTag            (config().pbtTag())
  {}

  PrefetchDAQData::~PrefetchDAQData() {
  }

//-----------------------------------------------------------------------------
  void PrefetchDAQData::beginJob() {
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

    if (_fetchCaloDigis) {
      auto cdH = evt.getValidHandle<CaloDigiCollection>(_cdTag);
      _cdcol = cdH.product();
    }

    if (_fetchTrkFragments) {
      auto tfH = evt.getValidHandle<artdaq::Fragments>(_tfTag);
      _tfcol = tfH.product();
    }

    if (_fetchCaloFragments) {
      auto cfH = evt.getValidHandle<artdaq::Fragments>(_cfTag);
      _cfcol = cfH.product();
    }

    if (_fetchStrawDigis) {
      auto sdH = evt.getValidHandle<StrawDigiCollection>(_sdTag);
      _sdcol = sdH.product();
    }

    if (_fetchAllFragments) {
      auto allH = evt.getMany<std::vector<artdaq::Fragment>>();
      for (auto& fs : allH) {
        for (auto& frag : *fs) {
          fake_access(frag);
        }
      }
    }

    if (_fetchEWM) {
      auto handle = evt.getValidHandle<EventWindowMarker>(_ewmTag);
      _ewm = handle.product();
    }

    if (_fetchPBT) {
      auto handle = evt.getValidHandle<ProtonBunchTime>(_pbtTag);
      _pbt = handle.product();
    }
//-----------------------------------------------------------------------------
// prefetch data
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

    if (_debugLevel > 0) printf("[PrefetchDAQData::produce] event number: %10i\n",_eventNum);

    findData(Event);
  }

//------------------------------------------------------------------------------
// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(PrefetchDAQData)

}



