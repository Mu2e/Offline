//
//  Filter for selecting events based on straw digis
//  Original author: Michael MacKenzie 03/2025
//  PM: make it real
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
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include "TRACE/tracemf.h"

// c++
#include <iostream>
#include <memory>
#include <map>

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"

// #include "Offline/DAQ/inc/TrkPanelMap_t.hh"

namespace mu2e {
  class StrawDigiFilter : public art::EDFilter {
  public:

    enum { kMaxNSamples = 100 };
    
    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string>     strawDigiCollTag   {Name("strawDigiCollTag")  , Comment("StrawDigiCollection tag"), "StrawDigisFromArtdaqFragments" };
      fhicl::Atom<unsigned>        minNStrawDigis     {Name("minNStrawDigis")    , Comment("minNStrawDigis")};
      fhicl::Atom<unsigned>        minNPlanes         {Name("minNPlanes")        , Comment("Minimum planes hit")};
      fhicl::Atom<int>             nSamplesBL         {Name("nSamplesBL")        , Comment("N (first) samples to determine the BL")};
      fhicl::Atom<float>           minPulseHeight     {Name("minPulseHeight")    , Comment("Minimum pulse height above the BL")};
      fhicl::Atom<float>           minGoodPulseHeight {Name("minGoodPulseHeight"), Comment("Minimum good pulse height above the BL")};
      fhicl::Atom<int>             minNGoodHits       {Name("minNGoodHits")      , Comment("Min numeber of good hits")};
      fhicl::Atom<int>             minNGoodPanels     {Name("minNGoodPanels")    , Comment("Min numebr of good panels")};
      fhicl::Atom<int>             debugLevel         {Name("debugLevel")        , Comment("Debug Level")};
      fhicl::Sequence<std::string> debugBits          {Name("debugBits"          )    , Comment("debug bits"                 ) };
      fhicl::Atom<bool  >          noFilter           {Name("noFilter")          , Comment("Don't filter anything")};
    };

    struct WfParam_t {
      int   fs;                         // first sample above _minPulseHeight
      float bl;                         // baseline
      float ph;                         // pulse height
      float q;                            // Q(positive)
      float qt;                           // Q(tail)
      float q_x_i;
      float ns;                           // nsamples in the charge integration
      float tm;                           // q_x_i/ns
    };

    struct Hist_t {
      TH1F* ph;
      TH1F* ngh;
      TH1F* ngp;
    } _hist;

    using Parameters = art::EDFilter::Table<Config>;

    explicit StrawDigiFilter(const Parameters& config);

    int book_histograms(int RunNumber);
    int fill_histograms();
    
    // takes unpacked waveform
    int       process_adc_waveform(float* Wf, int NSamples, WfParam_t* Wp);
    int       print_waveform      (const StrawDigiADCWaveform* Wf, WfParam_t* Wp);
    void      print_(const std::string&          Message,
                     const std::source_location& location = std::source_location::current());

    int                   offlineDtcID   (int DtcID);
    
  private:
    bool filter  (art::Event& ArtEvent) override;
    bool beginRun(art::Run&   Run     ) override;
    bool endRun  (art::Run&   Run     ) override;

    std::string   _sdTag;
    unsigned      _minndigis;
    unsigned      _minnplanes;
    int           _debugLevel;
    
    std::vector<std::string> _debugBits;
    int                      _debugBit[100];
    
    bool          _noFilter;
    int           _nSamplesBL;
    float         _minPulseHeight;
    float         _minGoodPulseHeight;
    int           _minNGoodHits;
    int           _minNGoodPanels;
    
    int           _ndigis;
    int           _n_good_hits[2][6];
    int           _n_good_panels;

    //     const TrkPanelMap_t* _panel_map[36][6];   // indexing: [plane][panel]
                                                     // counters
    unsigned      _nevt, _npass;
    const mu2e::StrawDigiCollection*            _digis;
    const mu2e::StrawDigiADCWaveformCollection* _sdwfc;

    std::vector<WfParam_t> _wfp;

    const art::Event* _event;
    int               _rn;

    ProditionsHandle<StrawElectronics> _stre_h;
    const StrawElectronics*            _strawElectronics;
    float                              _tdc_bin_ns;           // TDC bin, ns

    ProditionsHandle<TrackerPanelMap> _tpm_h;
    const TrackerPanelMap*            _trackerPanelMap;
  };

  StrawDigiFilter::StrawDigiFilter(const Parameters& conf)
    : art::EDFilter{conf},
      _sdTag             (conf().strawDigiCollTag()),
      _minndigis         (conf().minNStrawDigis()),
      _minnplanes        (conf().minNPlanes()),
      _debugLevel        (conf().debugLevel()),
      _debugBits         (conf().debugBits()),
      _noFilter          (conf().noFilter()),
      _nSamplesBL        (conf().nSamplesBL()),
      _minPulseHeight    (conf().minPulseHeight()),
      _minGoodPulseHeight(conf().minGoodPulseHeight()),
      _minNGoodHits      (conf().minNGoodHits()),
      _minNGoodPanels    (conf().minNGoodPanels()),
      _nevt              (0),
      _npass             (0)
  {
    // for (const TrkPanelMap_t* tpm = TrkPanelMap_data.begin(); tpm != TrkPanelMap_data.end(); ++tpm) {
    //   int plane = tpm->plane;
    //   int panel = tpm->panel;
    //   _panel_map[plane][panel] = tpm;
    // }
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
      
      TLOG(TLVL_DEBUG+1) << Form("... TrackerDQM: bit=%4i is set to %i\n",index,_debugBit[index]);
    }
  }


//-----------------------------------------------------------------------------
  void StrawDigiFilter::print_(const std::string& Message, const std::source_location& location) {
    std::cout << std::format(" event:{}:{}:{}",_event->run(),_event->subRun(),_event->event())
              << " " << location.file_name() << ":" << location.line()
      //            << location.function_name()
              << ": " << Message << std::endl;
  }
  
//-----------------------------------------------------------------------------
  int StrawDigiFilter::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
 
    _hist.ph            = tfs->make<TH1F>("ph"  , Form("run:%06i pulse height"   ,RunNumber),  100, 0., 1000.);
    _hist.ngh           = tfs->make<TH1F>("ngh" , Form("run:%06i N(good hits"    ,RunNumber),  100, 0.,  100.);
    _hist.ngp           = tfs->make<TH1F>("ngp" , Form("run:%06i N(good panels"  ,RunNumber),   20, 0.,   20.);
    
    return 0;
  }

  int StrawDigiFilter::fill_histograms() {

    for (int i=0; i<2; i++) {
      for (int panel=0; panel<6; panel++) {
        _hist.ngh->Fill(_n_good_hits[i][panel]);
      }
    }
    _hist.ngp->Fill(_n_good_panels);
    
    for(int i = 0; i < _ndigis; ++i) {
      _hist.ph->Fill(_wfp[i].ph);
    }
                     
    return 0;
  }


//-----------------------------------------------------------------------------
  int StrawDigiFilter::print_waveform(const StrawDigiADCWaveform* Wf, WfParam_t* Wp) {
    // print waveform
    std::string line;
    
    int ns = Wf->samples().size();
    int loc = 0;
    for (int i=0; i<ns; i++) {
      line += std::format(" {:4}",Wf->samples()[i]);
      loc++;
      if (loc >= 27) {
        printf("%s",line.data());
        line = "     ";
        loc  = 0;
      }
    }
    if (loc > 0) printf("%s",line.data());
    return 0;
  }
  
//-----------------------------------------------------------------------------
  int StrawDigiFilter::process_adc_waveform(float* Wf, int NSamples, WfParam_t* Wp) {
//-----------------------------------------------------------------------------
// straw man waveform processing
// 1. determine the baseline
//-----------------------------------------------------------------------------
    Wp->bl = 0;
    for (int i=0; i<_nSamplesBL; i++) {
      Wp->bl += Wf[i];
    }
    Wp->bl = Wp->bl/_nSamplesBL;
//-----------------------------------------------------------------------------
// 2. subtract the baseline and calculate the charge
//-----------------------------------------------------------------------------
    for (int i=0; i<NSamples; i++) {
      Wf[i] = Wf[i]-Wp->bl;
    }

    int tail =  0;
    Wp->fs   = -1;
    Wp->q    =  0;
    Wp->qt   =  0;
    Wp->ph   = -1;
    for (int i=_nSamplesBL; i<NSamples; i++) {
      if (Wf[i] > _minPulseHeight) {
        if (tail == 0) {
                                        // first sample above the threshold
          if (Wp->fs < 0) Wp->fs = i;

          Wp->q += Wf[i];
          if (Wf[i] > Wp->ph) {
            Wp->ph = Wf[i];
          }
        }
      }
      else if (Wf[i] < 0) {
        if (Wp->ph > 0) {
          tail  = 1;
        }
        if (tail == 1) Wp->qt -= Wf[i];
      }
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    if (Wp->q < 100) {
      TLOG(TLVL_DEBUG+1) << "event=" << _event->run() << ":"
                         << _event->subRun() << ":" << _event->event() 
                         << " Q=" << Wp->q;
    }
    return 0;
  }

  bool StrawDigiFilter::beginRun(art::Run& ArtRun) {
    _rn = ArtRun.run();
    
    art::EventID eid(_rn, 1, 1);
    _strawElectronics = &_stre_h.get(eid);
    _tdc_bin_ns = _strawElectronics->tdcLSB();      // 5./256 , ns
    // _tdc_bin             = (5./256.*1e-3);       // TDC bin width (Richie), in us
    ProditionsHandle<TrackerPanelMap> tpm_h;
    _trackerPanelMap = &tpm_h.get(eid);
      
    book_histograms(_rn);
    return true;
  }
  
//-----------------------------------------------------------------------------
  bool StrawDigiFilter::endRun(art::Run& run) {
    if(_debugLevel > 0){
      const float rate = (_nevt > 0) ? float(_npass)/float(_nevt) : 0.f;
      TLOG(TLVL_DEBUG) << "passed:" << _npass << " events out of:" << _nevt << " for a ratio of:" << rate;
    }
    return true;
  }
//-----------------------------------------------------------------------------
  bool StrawDigiFilter::filter(art::Event& ArtEvent) {

    _event         = &ArtEvent;         // should always be the first line
    
    if (_debugLevel > 0) print_("-- START");

    ++_nevt;

    _ndigis        = 0;
    _n_good_panels = 0;
    _wfp.clear();

    if (_debugLevel > 0) 
    
    for (int plane=0; plane<2; plane++) {
      for (int j=0; j<6; j++) {
        _n_good_hits[plane][j] = 0;
      }
    }

    // find the collection
    art::Handle<mu2e::StrawDigiCollection> digisH;
    if (!ArtEvent.getByLabel(_sdTag, digisH)) {
      TLOG(TLVL_ERROR) << "No straw digi collection found with tag:" << _sdTag.c_str();
      _digis = nullptr;
    } else {
      _digis  =  digisH.product();
      _ndigis = _digis->size();
    }

    art::Handle<mu2e::StrawDigiADCWaveformCollection> sdwfcH;
    if (!ArtEvent.getByLabel(_sdTag, sdwfcH)) {
      TLOG(TLVL_ERROR) << "No straw digi waveform collection found with tag:" << _sdTag.c_str();
      _sdwfc = nullptr;
    } else {
      _sdwfc = sdwfcH.product();
    }

    if(_debugLevel > 0) TLOG(TLVL_DEBUG) << "ndigis:" << _ndigis;
//-----------------------------------------------------------------------------
// process waveforms, count good hits
//-----------------------------------------------------------------------------
    WfParam_t wfpar;
    float     wf_data[kMaxNSamples];

    int header_printed(0);

    for(int i = 0; i<_ndigis; ++i) {
      const mu2e::StrawDigi* sd = &_digis->at(i);
      uint16_t     plane = sd->strawId().plane();  // this is the GEO ID of the panel
      uint16_t     panel = sd->strawId().panel();        // this is the GEO ID of the panel
      const TrkPanelMap::Row* tpm;
      try {
        tpm = _trackerPanelMap->panel_map_by_online_ind(plane,panel);
      }
      catch(...) {
//-----------------------------------------------------------------------------
// either DTC ID or link ID are corrupted. Haven't seen that so far, switch to the next ROC anyway
//-----------------------------------------------------------------------------
        print_(std::format("ERROR: either plane:{} or panel:{} is corrupted, skip ROC data\n",
                           plane,panel));
        continue;
      }
//-----------------------------------------------------------------------------
// unpacking: 
//-----------------------------------------------------------------------------
      const mu2e::StrawDigiADCWaveform* wf = &_sdwfc->at(i);
      int ns = wf->samples().size();

      if (ns > kMaxNSamples) {
        TLOG(TLVL_ERROR) << "ns:" << ns << " greater than kMaxNSamples:" << kMaxNSamples << ", truncate";
        ns = kMaxNSamples;
      }

      for (int i=0; i<ns; i++) {
        wf_data[i] = wf->samples().at(i);
      }
      
      process_adc_waveform(wf_data,ns,&wfpar);

      _wfp.push_back(wfpar);

      if (_debugLevel) {
        StrawId sid = sd->strawId();
        StrawDigiFlagDetail::mask_type /*uint8_t*/ sd_flag = *((StrawDigiFlagDetail::mask_type*) &sd->digiFlag());
        float t0        = sd->TDC()[0]*_tdc_bin_ns;  // in ns
        float t1        = sd->TDC()[1]*_tdc_bin_ns;

        int print_it = 0;
        if      (_debugLevel & _debugBit[0]) {                                     print_it = 1; }
        else if (_debugLevel & _debugBit[1]) { if (wfpar.ph > _minGoodPulseHeight) print_it = 1;}

        if (print_it) {
          if (header_printed == 0) {
            std::cout << std::format("  sid  station plane panel  name straw    T0       T1   FLAG  {:68s}     FS     PH       BL\n"," " );
            header_printed = 1;
          }
          std::cout << std::format("0x{:04x} {:5}   {:4}  {:4}    MN{:03d} {:3} {:8.1f} {:8.1f} 0x{:02x}",
                                   sid.asUint16(), sid.station(),sid.plane(),sid.panel(),tpm->mnid(),sid.straw(),t0,t1,sd_flag);
          print_waveform(wf,&wfpar);
          std::cout << std::format(" {:2d} {:8.3f} {:8.3f} ",wfpar.fs,wfpar.ph,wfpar.bl) << std::endl;
        }
      }

      if (wfpar.ph > _minGoodPulseHeight) {
        _n_good_hits[tpm->uniquePlane()][tpm->panel()]++;
      }
    }

    for (int plane=0; plane<2; plane++) {
      std::cout << "plane:" << plane << " ";
      for (int j=0; j<6; j++) {
        if (_n_good_hits[plane][j] >= _minNGoodHits) _n_good_panels++;
        std::cout << " j:" << j << " _n_good_hits:" << _n_good_hits[plane][j];
      }
      std::cout << std::endl;
    }
    
    fill_histograms();

    if (_debugLevel > 0) print_(std::format("-- END: n_good_panels:{}",_n_good_panels));

    return (_n_good_panels >= _minNGoodPanels);
  }

}

DEFINE_ART_MODULE(mu2e::StrawDigiFilter)
