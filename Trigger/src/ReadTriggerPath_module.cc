//
// An EDAnalyzer module that reads the Trigger Info to store the trigger path cut-flow
//
// Original author M. MacKenzie, based on ReadTriggerInfo
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "canvas/Utilities/InputTag.h"

//Services

//Data products
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

//Gen-level info
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

//Utilities
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"

//ROOT
#include "TH1.h"
#include "TH2.h"

#include <cmath>
#include <string>
#include <sstream>
#include <vector>

namespace mu2e {

  class ReadTriggerPath : public art::EDAnalyzer {

  public:
    enum {
      kMaxTriggers = 200,
      kNTrigInfo = 20 //histogram sets
    };

    // summary information about all triggers
    struct summaryInfoHist_ {
      TH1* _hTrigIndices; //trigger path indices that are firing
      TH1* _hTrigBits   ; //trigger path bits that are firing
      TH1* _hTrigPaths  ; //trigger path names that are firing
      TH2* _hTrig2D     ; //trigger path acceptance correlation
      TH1* _hNPOT       ; //N(POT) without cuts, matching TrigPOTEff binning
      TH1* _hTrigInfo   [kNTrigInfo  ]; //1D trigger information
      TH1* _hTrigModules[kMaxTriggers]; //cut-flow on the module chain for each trigger
      TH1* _hTrigPOTEff [kMaxTriggers]; //trigger efficiencies vs. N(POT)


      summaryInfoHist_() {
        // initialize each histogram to null
        _hTrigIndices = nullptr;
        _hTrigBits    = nullptr;
        _hTrigPaths   = nullptr;
        _hTrig2D      = nullptr;
        _hNPOT        = nullptr;
        for (int i = 0; i < kNTrigInfo; ++i) {
          _hTrigInfo   [i] = nullptr;
        }
        for (int i = 0; i < kMaxTriggers; ++i) {
          _hTrigModules[i] = nullptr;
          _hTrigPOTEff [i] = nullptr;
        }
      }
    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>               diagLevel          {Name("diagLevel"          ), Comment("Diagnostic printouts"                      ), 0 };
      fhicl::Sequence<std::string>   ignorePaths        {Name("ignorePaths"        ), Comment("Path tags to ignore in overlap/rejection"  ), {"Path"}}; //non-trigger paths typically contain "Path"
      fhicl::Atom<art::InputTag>     genCountTag        {Name("genCount"           ), Comment("GenEventCount label"                       ), "genCounter" };
      fhicl::Atom<art::InputTag>     PBITag             {Name("PBITag"             ), Comment("ProtonBunchIntensity label"                ), "PBISim" };
      fhicl::Atom<std::string>       processName        {Name("processName"        ), Comment("globalTrigger"                             ), ""  };
    };

    explicit ReadTriggerPath(const art::EDAnalyzer::Table<Config>& config);
    virtual ~ReadTriggerPath() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(const art::Run & run);
    virtual void beginSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);

    void     bookHistograms           ();
    void     bookTrigInfoHist         (art::ServiceHandle<art::TFileService>& Tfs, summaryInfoHist_& Hist);

    bool     isTrackFilter            (const std::string& module);
    bool     isHelixFilter            (const std::string& module);
    bool     isCaloFilter             (const std::string& module);
    bool     isGlobalFilter           (const std::string& module);
    bool     isTrackPath              (const std::string& path, TriggerResultsNavigator& trigNavig);
    bool     isHelixPath              (const std::string& path, TriggerResultsNavigator& trigNavig);
    bool     isCaloPath               (const std::string& path, TriggerResultsNavigator& trigNavig);
    bool     isGlobalPath             (const std::string& path, TriggerResultsNavigator& trigNavig);
    bool     isMinBiasPath            (const std::string& path);

    void     trigPathVal              (const int Index, TriggerResultsNavigator& trigNavig, summaryInfoHist_& Hist);
    void     evalTriggerRate          ();

  private:
    int                       _diagLevel;
    std::vector<std::string>  _ignorePaths;
    art::InputTag             _genCountTag;
    art::InputTag             _PBITag;
    std::string               _processName;

    size_t                    _nMaxTrig;
    std::vector<std::string>  _trigPaths;

    int                       _nProcess;
    double                    _nPOT;

    summaryInfoHist_          _sumHist;

    float  _minPOT, _maxPOT;

    bool _useNGen; //use gen event count for normalization if available
  };

  ReadTriggerPath::ReadTriggerPath(const art::EDAnalyzer::Table<Config>& config):
    art::EDAnalyzer{config},
    _diagLevel           (config() .diagLevel()      ),
    _ignorePaths         (config() .ignorePaths()    ),
    _genCountTag         (config() .genCountTag()    ),
    _PBITag              (config() .PBITag()         ),
    _processName         (config() .processName()    ),
    _nMaxTrig            (0                          ),
    _nProcess            (0                          ),
    _minPOT              (0.                         ),
    _maxPOT              (2.0e8                      )
  {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
  }

  //--------------------------------------------------------------------------------//
  // Initialize all output histograms
  void ReadTriggerPath::bookHistograms() {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    _nMaxTrig = _trigPaths.size();
    if(_nMaxTrig > kMaxTriggers) throw cet::exception("BADCONFIG") << __func__ << ": Input trigger path list size (" << _nMaxTrig << ") is larger than the maximum ("
                                                                   << kMaxTriggers << ")\n";
    art::ServiceHandle<art::TFileService> tfs;
    bookTrigInfoHist(tfs, _sumHist);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::bookTrigInfoHist(art::ServiceHandle<art::TFileService>& Tfs, summaryInfoHist_& Hist) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    art::TFileDirectory trigInfoDir = Tfs->mkdir("trigInfo");
    const int max_bits = 1000;
    Hist._hTrigBits      = trigInfoDir.make<TH1F>("hTrigBits"         , "Trigger bits"           ,  max_bits, -0.5, max_bits-0.5);
    Hist._hTrigPaths     = trigInfoDir.make<TH1F>("hTrigPaths"        , "Trigger paths"          , _nMaxTrig, -0.5, float(_nMaxTrig)-0.5);
    Hist._hTrigIndices   = trigInfoDir.make<TH1F>("hTrigIndices"      , "Trigger indices"        , _nMaxTrig, -0.5, float(_nMaxTrig)-0.5);

    Hist._hTrigInfo[0]  = trigInfoDir.make<TH1F>("hTrigEff"           , "Trigger efficiencies"   , _nMaxTrig, -0.5, float(_nMaxTrig)-0.5);
    Hist._hTrigInfo[1]  = trigInfoDir.make<TH1F>("hTrigRej"           , "Trigger rejections"     , _nMaxTrig, -0.5, float(_nMaxTrig)-0.5);
    Hist._hTrigInfo[2]  = trigInfoDir.make<TH1F>("hGlobalTrigEff"     , "Trigger efficiencies"   ,        10, -0.5, 9.5);
    Hist._hTrigInfo[3]  = trigInfoDir.make<TH1F>("hGlobalTrigRej"     , "Trigger rejections"     ,        10, -0.5, 9.5);
    Hist._hTrigInfo[4]  = trigInfoDir.make<TH1F>("hExclusiveTrigEff"  , "Exclusive trigger eff." , _nMaxTrig, -0.5, float(_nMaxTrig)-0.5);
    Hist._hTrigInfo[5]  = trigInfoDir.make<TH1F>("hExclusiveTrigRej"  , "Exclusive trigger rej." , _nMaxTrig, -0.5, float(_nMaxTrig)-0.5);

    const int npot_bins(200);
    Hist._hNPOT         = trigInfoDir.make<TH1F>("hNPOT"              , "N(POT);N(POT)/#mu-bunch", npot_bins, _minPOT, _maxPOT);

    Hist._hTrig2D       = trigInfoDir.make<TH2F>("hTrigOverlap"       , "Trigger overlap: N(x and y)/N(x)", _nMaxTrig, -0.5, _nMaxTrig-0.5, _nMaxTrig, -0.5, _nMaxTrig);

    art::TFileDirectory cutflowDir = Tfs->mkdir("cutflowDir");
    art::TFileDirectory potEffDir  = Tfs->mkdir("POTEffDir");
    //add a for loop for creating histograms showing the trigger path module filter path as well as the efficiency vs. N(POT)
    const int max_modules = 40;
    for(unsigned i = 0; i < _nMaxTrig; ++i) {
      Hist._hTrigModules[i] = cutflowDir.make<TH1F>(Form("hCutflow_%d",i), "" , max_modules, 0, max_modules);
      Hist._hTrigPOTEff [i] = potEffDir .make<TH1F>(Form("POTEff_%d"  , i), "", npot_bins, _minPOT, _maxPOT);
    }
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isTrackFilter(const std::string& module) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    return module.find("KSFilter") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isHelixFilter(const std::string& module) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    return module.find("HSFilter") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isCaloFilter(const std::string& module) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    return module.find("calo") != std::string::npos && module.find("Filter") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isGlobalFilter(const std::string& module) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    return (module.find("HSFilter")!= std::string::npos && module.find("CosmicHelix")!= std::string::npos) ||
      module.find("caloPhotonFilter")!= std::string::npos ||
      module.find("caloMVANNCEFilter")!= std::string::npos ||
      module.find("TSFilter")       != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isTrackPath(const std::string& path, TriggerResultsNavigator& trigNavig) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    bool passed = false;
    // check if the module labels are available for a more reliable check
    auto modules = trigNavig.triggerModules(path);
    if(modules.size() > 0) {
      for (auto module : modules) passed |= isTrackFilter(module);
    } else { // use the path name to determine if it's a track path
      passed |= path.find("tpr") != std::string::npos;
      passed |= path.find("cpr") != std::string::npos;
      passed |= path.find("apr") != std::string::npos;
      passed &= path.find("Helix") == std::string::npos;
    }
    return passed;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isHelixPath(const std::string& path, TriggerResultsNavigator& trigNavig) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    bool passed = false;
    // check if the module labels are available for a more reliable check
    auto modules = trigNavig.triggerModules(path);
    if(modules.size() > 0) {
      for (auto module : modules) passed |= isHelixFilter(module);
      passed &= !isTrackPath(path, trigNavig);
    } else { // use the path name to determine if it's a helix path
      passed |= path.find("tpr") != std::string::npos;
      passed |= path.find("cpr") != std::string::npos;
      passed |= path.find("apr") != std::string::npos;
      passed &= path.find("Helix") != std::string::npos;
    }
    return passed;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isCaloPath(const std::string& path, TriggerResultsNavigator& trigNavig) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    bool passed = false;
    // check if the module labels are available for a more reliable check
    auto modules = trigNavig.triggerModules(path);
    if(modules.size() > 0) {
      for (auto module : modules) passed |= isCaloFilter(module);
    } else { // use the path name to determine if it's a calo-only path
      passed |= path.find("calo") != std::string::npos;
      passed &= !(isTrackPath(path, trigNavig) || isHelixPath(path, trigNavig));
    }
    return passed;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isGlobalPath(const std::string& path, TriggerResultsNavigator& trigNavig) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    bool passed = false;
    // check if the module labels are available for a more reliable check
    auto modules = trigNavig.triggerModules(path);
    if(modules.size() > 0) {
      for (auto module : modules) passed |= isGlobalFilter(module);
    } else { // use the path name to determine if it's a global path
      passed = false; // FIXME: Determine global path names
    }
    return passed;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerPath::isMinBiasPath(const std::string& path) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    return path.find("minBias") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::trigPathVal(const int Index, TriggerResultsNavigator& trigNavig, summaryInfoHist_& Hist) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    const std::string path = trigNavig.getTrigPathNameByIndex(Index);
    const unsigned lastModule = trigNavig.indexLastModule(path);
    auto h = Hist._hTrigModules[Index];
    if(!h) throw cet::exception("BADCONFIG") << __func__ << ": Trigger path index " << Index << " is out of bounds for initialized histograms\n";
    if(h->GetEntries() == 0) h->SetTitle((path + " path cut-flow").c_str());

    // set the bin labels if needed
    if(h->Integral() <= 0.) {
      auto modules = trigNavig.triggerModules(path);
      auto size = modules.size();
      if(size > 0) {
        for(unsigned i = 0; i < size; ++i) {
          h->GetXaxis()->SetBinLabel(h->FindBin(i), (i < size) ? modules[i].c_str() : "Final");
        }
      }
    }

    // cut-flow along the trigger path
    for (unsigned i = 0; i < lastModule ; ++i) h->Fill(i);
  }


  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::beginJob() {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::endJob() {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    if(_nProcess <= 0) {
      if(_diagLevel > 0) printf("[ReadTriggerPath::%s] Setting N(processed) from %i to 1\n", __func__, _nProcess);
      _nProcess = 1;
    }
    if(_diagLevel > 0) printf("[ReadTriggerPath::%s] N(processed) = %i\n", __func__, _nProcess);

    //////////////////////////////////////////
    // Evaluate summary info

    if(!_sumHist._hTrigInfo[0]) return; // histograms were not booked

    // trigger efficiencies
    _sumHist._hTrigInfo[0]->Scale(1./_nProcess);
    _sumHist._hTrigInfo[2]->Scale(1./_nProcess);
    _sumHist._hTrigInfo[4]->Scale(1./_nProcess);

    // make trigger correlation efficiency relative to x-axis
    for(int xbin = 1; xbin <= _sumHist._hTrig2D->GetNbinsX(); ++xbin) {
      const double xval = _sumHist._hTrig2D->GetBinContent(xbin, xbin);
      if(xval <= 0.) continue;
      for(int ybin = 1; ybin <= _sumHist._hTrig2D->GetNbinsY(); ++ybin) {
        const double binc = _sumHist._hTrig2D->GetBinContent(xbin, ybin);
        _sumHist._hTrig2D->SetBinContent(xbin, ybin, binc/xval);
      }
    }

    // trigger rejections
    for(int bin = 0; bin <= _sumHist._hTrigInfo[1]->GetNbinsX()+1; ++bin) {
      const double binc = _sumHist._hTrigInfo[1]->GetBinContent(bin);
      if(binc > 0.) _sumHist._hTrigInfo[1]->SetBinContent(bin, _nProcess/binc);
    }
    for(int bin = 0; bin <= _sumHist._hTrigInfo[3]->GetNbinsX()+1; ++bin) {
      const double binc = _sumHist._hTrigInfo[3]->GetBinContent(bin);
      if(binc > 0.) _sumHist._hTrigInfo[3]->SetBinContent(bin, _nProcess/binc);
    }
    for(int bin = 0; bin <= _sumHist._hTrigInfo[5]->GetNbinsX()+1; ++bin) {
      const double binc = _sumHist._hTrigInfo[5]->GetBinContent(bin);
      if(binc > 0.) _sumHist._hTrigInfo[5]->SetBinContent(bin, _nProcess/binc);
    }

    // evaluate eff. vs. N(POT) plots
    const double npot_norm = _sumHist._hNPOT->Integral() / _nProcess;
    //manually divide to set binomial errors
    for(int bin = 1; bin <= _sumHist._hNPOT->GetNbinsX(); ++bin) {
      const double npot = _sumHist._hNPOT->GetBinContent(bin);
      if(npot <= 0.) continue;
      for(unsigned ihist = 0; ihist < _nMaxTrig; ++ihist) {
        auto h = _sumHist._hTrigPOTEff[ihist];
        const double val = h->GetBinContent(bin);
        const double rate = val/npot;
        const double err  = std::sqrt(rate*(1.-rate)/npot);
        h->SetBinContent(bin, rate*npot_norm);
        h->SetBinError  (bin, err *npot_norm);
      }
    }

    //////////////////////////////////////////
    // Print the cut-flow information

    if(_diagLevel > 0) {
      printf("[ReadTriggerPath::%s] Printing trigger path cut-flow information:\n", __func__);
      for(size_t trig = 0; trig < _nMaxTrig; ++trig) {
        auto h = _sumHist._hTrigModules[trig];
        //check if the trigger was defined
        if(std::string(h->GetTitle()).empty()) continue;
        printf(" %s:\n", h->GetTitle());
        for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
          const double binc = h->GetBinContent(ibin);
          if(binc <= 0.) continue;
          std::string label(h->GetXaxis()->GetBinLabel(ibin));
          if(label.empty()) printf("  %2i: %7.0f (%7.3f%%)\n", ibin, binc, binc*100./_nProcess);
          else              printf("  %40s: %7.0f (%7.3f%%)\n", label.c_str(), binc, binc*100./_nProcess);
        }
      }
    }
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::beginRun(const art::Run & run) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::beginSubRun(const art::SubRun& sr) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);
    art::Handle<mu2e::GenEventCount> genCountH;
    sr.getByLabel(_genCountTag, genCountH);
    _useNGen = genCountH.isValid();
    if(_useNGen) _nProcess += genCountH.product()->count();
    if(_diagLevel > 0) printf("[ReadTriggerPath::%s] use gen count product flag = %o\n", __func__, _useNGen);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerPath::analyze(const art::Event& event) {
    if(_diagLevel > 4) printf("[ReadTriggerPath::%s]\n", __func__);

    if(!_useNGen) _nProcess += 1;

    //get the number of POT
    art::Handle<mu2e::ProtonBunchIntensity> PBIH;
    event.getByLabel(_PBITag, PBIH);
    if(PBIH.isValid()) _nPOT = PBIH.product()->intensity();
    else _nPOT  = -1.;
    if(_diagLevel > 1) printf("[ReadTriggerPath::%s] N(POT) = %.2e\n", __func__, _nPOT);

    //get the TriggerResult
    std::ostringstream oss;
    oss << "TriggerResults::"<<_processName;
    art::InputTag const tag{oss.str()};
    if(_diagLevel > 1) printf("[ReadTriggerPath::%s] Trigger results tag = %s\n", __func__, oss.str().c_str());
    auto const trigResultsH   = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults*trigResults = trigResultsH.product();
    TriggerResultsNavigator   trigNavig(trigResults);
    if(_diagLevel > 1) trigNavig.print();

    //now set the value of _trigPaths
    _trigPaths = trigNavig.getTrigPaths();

    //check if the histograms have been initialized
    if(!_sumHist._hTrigInfo[0]) bookHistograms();

    if(_trigPaths.size() > _nMaxTrig) throw cet::exception("BADCONFIG") << "Number of triggers assumed (" << _nMaxTrig << ") is less than the observed (" << _trigPaths.size() << ")\n";

    //initialize bin labels
    if(_sumHist._hTrigInfo[0]->Integral() <= 0.) {
      for (unsigned int i=0; i< trigNavig.getTrigPaths().size(); ++i) {
        const std::string path = trigNavig.getTrigPathNameByIndex(i);
        _sumHist._hTrigInfo[0]->GetXaxis()->SetBinLabel(_sumHist._hTrigInfo[0]->FindBin(i), path.c_str());
        _sumHist._hTrigInfo[1]->GetXaxis()->SetBinLabel(_sumHist._hTrigInfo[1]->FindBin(i), path.c_str());
        _sumHist._hTrigInfo[4]->GetXaxis()->SetBinLabel(_sumHist._hTrigInfo[4]->FindBin(i), path.c_str());
        _sumHist._hTrigInfo[5]->GetXaxis()->SetBinLabel(_sumHist._hTrigInfo[5]->FindBin(i), path.c_str());
        _sumHist._hTrig2D     ->GetXaxis()->SetBinLabel(_sumHist._hTrig2D->GetXaxis()->FindBin(i), path.c_str());
        _sumHist._hTrig2D     ->GetYaxis()->SetBinLabel(_sumHist._hTrig2D->GetYaxis()->FindBin(i), path.c_str());
        _sumHist._hTrigPaths  ->Fill(path.c_str(), 0.);
      }
    }
    if(_sumHist._hTrigInfo[2]->Integral() <= 0.) {
      std::vector<std::string> labels = {"Total", "Tracks", "Helix", "Calo", "Calib", "MinBias", "Unknown"};
      for(unsigned index = 0; index < labels.size(); ++index) {
        _sumHist._hTrigInfo[2]->GetXaxis()->SetBinLabel(index+1, labels[index].c_str());
        _sumHist._hTrigInfo[3]->GetXaxis()->SetBinLabel(index+1, labels[index].c_str());
      }
    }

    //fill the histogram with the accepted trigger bits
    _sumHist._hNPOT->Fill(_nPOT);
    bool passed(false), track_passed(false), helix_passed(false), calo_passed(false), minbias_passed(false), unknown(false);
    unsigned naccept(0), exclusive_idx(-1);
    for (unsigned i=0; i< trigNavig.getTrigPaths().size(); ++i) {
      const std::string path = trigNavig.getTrigPathNameByIndex(i);
      if(_diagLevel > 4) {
        printf("[ReadTriggerPath::%s] Printing modules in path %s\n", __func__, path.c_str());
        auto modules = trigNavig.triggerModules(path);
        for(auto module : modules) printf(" %s\n", module.c_str());
      }

      if(trigNavig.accepted(path)) {
        _sumHist._hTrigBits->Fill(trigNavig.getTrigBitByName(path));
        _sumHist._hTrigIndices->Fill(trigNavig.getTrigPathIndex(path));
        _sumHist._hTrigPaths  ->Fill(path.c_str(), 1.);
        // _sumHist._hTrigPaths  ->Fill(i);
        _sumHist._hTrigInfo[0]->Fill(i);
        _sumHist._hTrigInfo[1]->Fill(i);
        if(_sumHist._hTrigPOTEff[i]->GetEntries() == 0) _sumHist._hTrigPOTEff[i]->SetTitle(Form("%s Eff. vs. N(POT);N(POT);#epsilon", path.c_str()));
        _sumHist._hTrigPOTEff[i]->Fill(_nPOT);
        bool skip_exclusive = false;
        for(auto tag : _ignorePaths) skip_exclusive |= path.find(tag) != std::string::npos;
        if(!skip_exclusive) {
          ++naccept;
          exclusive_idx = i;
        }
        const bool is_track   = isTrackPath  (path, trigNavig) && !skip_exclusive;
        const bool is_helix   = isHelixPath  (path, trigNavig) && !skip_exclusive;
        const bool is_calo    = isCaloPath   (path, trigNavig) && !skip_exclusive;
        const bool is_minbias = isMinBiasPath(path) && !skip_exclusive;
        passed         |= skip_exclusive;
        track_passed   |= is_track;
        helix_passed   |= is_helix;
        calo_passed    |= is_calo;
        minbias_passed |= is_minbias;
        unknown |= !(is_track || is_helix || is_calo || is_minbias || skip_exclusive);
        if(is_track + is_helix + is_calo + is_minbias > 1) {
          if(_diagLevel > -1) printf("[ReadTriggerPath::%s] Path %s identified as multiple types: track = %o helix = %o calo = %o minbias = %o\n",
                                     __func__, path.c_str(), is_track, is_helix, is_calo, is_minbias);
        }
        //Fill the correlation matrix
        for (unsigned j=0; j < trigNavig.getTrigPaths().size(); ++j) {
          if(trigNavig.accepted(trigNavig.getTrigPathNameByIndex(j)) && !skip_exclusive) _sumHist._hTrig2D->Fill(i,j);
        }
      }

      //fill the trigger path cut-flow
      trigPathVal(i, trigNavig, _sumHist);
    }

    // exclusive trigger acceptance tracking
    if(naccept == 1) {
      _sumHist._hTrigInfo[4]->Fill(exclusive_idx);
      _sumHist._hTrigInfo[5]->Fill(exclusive_idx);
    }

    // trigger rates by category
    if(passed) { //total trigger rates
      _sumHist._hTrigInfo[2]->Fill(0);
      _sumHist._hTrigInfo[3]->Fill(0);
    }
    if(track_passed) { //track trigger rates
      _sumHist._hTrigInfo[2]->Fill(1);
      _sumHist._hTrigInfo[3]->Fill(1);
    }
    if(track_passed) { //helix trigger rates
      _sumHist._hTrigInfo[2]->Fill(2);
      _sumHist._hTrigInfo[3]->Fill(2);
    }
    if(calo_passed) { //calo-only trigger rates
      _sumHist._hTrigInfo[2]->Fill(3);
      _sumHist._hTrigInfo[3]->Fill(3);
    }
    if(minbias_passed) { //min-bias trigger rates
      _sumHist._hTrigInfo[2]->Fill(5);
      _sumHist._hTrigInfo[3]->Fill(5);
    }
    if(unknown) { //unknown trigger rates
      _sumHist._hTrigInfo[2]->Fill(6);
      _sumHist._hTrigInfo[3]->Fill(6);
    }
  }
}

DEFINE_ART_MODULE(mu2e::ReadTriggerPath)
