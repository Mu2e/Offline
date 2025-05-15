//
//  PM: straw hit filter, reject the noise
//
//  debugBit[0] : print all WFs
//  debugBit[1] : print good WFs
#include "TH1F.h"
//#include "TFolder.h"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
// mu2e
// data
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "TRACE/tracemf.h"

// c++
#include <iostream>
#include <memory>
#include <map>

namespace mu2e {
  class StrawHitFilter : public art::EDFilter {
  public:

    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool>            debugMode     {Name("debugMode"     ), Comment("Debug Mode, default:false")};
      fhicl::Sequence<std::string> debugBits     {Name("debugBits"     ), Comment("debug bits")};
      fhicl::Atom<std::string>     shCollTag     {Name("shCollTag"     ), Comment("StrawHitCollection tag") };
      fhicl::Atom<float>           maxDt         {Name("maxDt"         ), Comment("max abs(DT)")};
      fhicl::Atom<float>           minEDep       {Name("minEDep"       ), Comment("min EDep")};
      fhicl::Atom<int>             minNGoodHits  {Name("minNGoodHits"  ), Comment("minNStrawDigis")};
      fhicl::Atom<bool>            fillHistograms{Name("fillHistograms"), Comment("fill histogrms, default:false")};
    };

    struct Hist_t {
      TH1F* evt ;
      TH1F* nsht;
      TH1F* nshg;
      TH1F* dt  ;
      TH1F* edep;
    } _hist[2];

    using Parameters = art::EDFilter::Table<Config>;

    explicit StrawHitFilter(const Parameters& config);

    int book_histograms(int RunNumber);
    int fill_histograms(Hist_t* Hist);
    
    void      print_(const std::string&          Message,
                     const std::source_location& location = std::source_location::current());

  private:
    bool filter  (art::Event& ArtEvent) override;
    bool beginRun(art::Run&   Run     ) override;
    bool endRun  (art::Run&   Run     ) override;

    std::string              _shCollTag;

    bool                     _debugMode;
    std::vector<std::string> _debugBits;
    int                      _debugBit[100];
    
    float                    _maxDt;
    float                    _minEDep;
    int                      _minNGoodHits;
    int                      _fillHistograms;
    
    int                      _nevt;
    int                      _nevp;
    int                      _nsht;
    int                      _nshg;

    const mu2e::StrawHitCollection* _shc;

    const art::Event* _event;
    int               _rn;
  };

//-----------------------------------------------------------------------------
  StrawHitFilter::StrawHitFilter(const Parameters& conf)
    : art::EDFilter{conf},
      _shCollTag         (conf().shCollTag()),
      _debugMode         (conf().debugMode()),
      _debugBits         (conf().debugBits()),
      _maxDt             (conf().maxDt    ()),
      _minEDep           (conf().minEDep  ()),
      _minNGoodHits      (conf().minNGoodHits()),
      _fillHistograms    (conf().fillHistograms())
  {
//-----------------------------------------------------------------------------
// parse debug bits
//-----------------------------------------------------------------------------
    const char* key;
    // a flag is an integer!
    int nbits = _debugBits.size();
    for (int i=0; i<nbits; i++) {
      int index(0), value(0);
      key               = _debugBits[i].data();
      sscanf(key,"bit%i:%i",&index,&value);
      _debugBit[index]  = value;
    }
  }


//-----------------------------------------------------------------------------
  void StrawHitFilter::print_(const std::string& Message, const std::source_location& location) {
    std::cout << std::format(" event:{}:{}:{}",_event->run(),_event->subRun(),_event->event())
              << " " << location.file_name() << ":" << location.line()
      //            << location.function_name()
              << ": " << Message << std::endl;
  }
  
//-----------------------------------------------------------------------------
  int StrawHitFilter::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;

     TH1::AddDirectory(kFALSE);
   
     art::TFileDirectory d1 = tfs->mkdir("all");
     
    _hist[0].evt  = d1.make<TH1F>("evt" , Form("run:%06i Event number[0]", RunNumber), 10000,   0., 1.e7 );
    _hist[0].nsht = d1.make<TH1F>("nsht", Form("run:%06i N(shT)      [0]", RunNumber),   500, -0.5, 499.5);
    _hist[0].nshg = d1.make<TH1F>("nshg", Form("run:%06i N(shG)      [0]", RunNumber),   100, -0.5,  99.5);
    _hist[0].dt   = d1.make<TH1F>("dt"  , Form("run:%06i delta(T), ns[0]", RunNumber),   200, -100, 100. );
    _hist[0].edep = d1.make<TH1F>("edep", Form("run:%06i edep, keV   [0]", RunNumber),   100, -  2,   8. );
    
    art::TFileDirectory d2 = tfs->mkdir("passed");

    _hist[1].evt  = d2.make<TH1F>("evt" , Form("run:%06i Event number[1]", RunNumber), 10000,   0., 1.e7 );
    _hist[1].nsht = d2.make<TH1F>("nsht", Form("run:%06i N(shT)      [1]", RunNumber),   500, -0.5, 499.5);
    _hist[1].nshg = d2.make<TH1F>("nshg", Form("run:%06i N(shG)      [1]", RunNumber),   100, -0.5,  99.5);
    _hist[1].dt   = d2.make<TH1F>("dt"  , Form("run:%06i delta(T), ns[1]", RunNumber),   200, -100, 100. );
    _hist[1].edep = d2.make<TH1F>("edep", Form("run:%06i edep, keV   [1]", RunNumber),   100, -  2,   8. );

    return 0;
  }
  
  
//-----------------------------------------------------------------------------
  int StrawHitFilter::fill_histograms(Hist_t* Hist) {

    for (int i = 0; i<_nsht; ++i) {
      const mu2e::StrawHit* sh = &_shc->at(i);

      Hist->evt->Fill(_event->event());
      Hist->dt->Fill(sh->dt());
      Hist->edep->Fill(sh->energyDep()*1.e3);
    }
  
    Hist->nsht->Fill(_nsht);
    Hist->nshg->Fill(_nshg);
                     
    return 0;
  }

//-----------------------------------------------------------------------------
  bool StrawHitFilter::beginRun(art::Run& ArtRun) {
    _rn = ArtRun.run();

    if (_fillHistograms) {
                                        // make sure the job doesn't crash
      book_histograms(_rn);
    }
    return true;
  }
  
//-----------------------------------------------------------------------------
  bool StrawHitFilter::endRun(art::Run& run) {
    const float rate = (_nevt > 0) ? float(_nevp)/float(_nevt) : 0.f;
    TLOG(TLVL_DEBUG) << "passed:" << _nevp << " events out of:" << _nevt << " for a ratio of:" << rate;
    return true;
  }


  //-----------------------------------------------------------------------------
  bool StrawHitFilter::filter(art::Event& ArtEvent) {
    
    _event         = &ArtEvent;         // should always be the first line
    
    if (_debugMode) print_("-- START");
    
    ++_nevt;
    
    
    art::Handle<mu2e::StrawHitCollection> shch;
    if (!ArtEvent.getByLabel(_shCollTag, shch)) {
      TLOG(TLVL_ERROR) << "No straw hit collection tag:" << _shCollTag.data();
      _shc = nullptr;
    }
    else {
      _shc   =  shch.product();
      _nsht  = _shc->size();
    }

//-----------------------------------------------------------------------------
// process waveforms, count good hits
//-----------------------------------------------------------------------------
    _nshg          = 0;
    
    for (int i = 0; i<_nsht; ++i) {
      const mu2e::StrawHit* sh = &_shc->at(i);
      if (fabs(sh->dt())  > _maxDt  ) continue;
      if (sh->energyDep() < _minEDep) continue;
      _nshg++;
    }
  
    if (_fillHistograms) fill_histograms(&_hist[0]);
  
    if (_debugMode) print_(std::format("-- END, n good hits:{}",_nshg));
    
    bool passed = false;
    if (_nshg >= _minNGoodHits) {
      passed = true;
      _nevp++;
      if (_fillHistograms) fill_histograms(&_hist[1]);
    }
  
    return passed;
  }
}

DEFINE_ART_MODULE(mu2e::StrawHitFilter)
