#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "Offline/CalPatRec/inc/MergePatRec_types.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"

using namespace std;

namespace mu2e {
  using namespace MergePatRecTypes;
  using           CLHEP::Hep3Vector;

  class SimParticle;

  class MergePatRecDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets =  10,
      kNCprHistSets   =  10,
      kNTprHistSets   =  10
    };

    struct TrackHist_t {
      TH1F*  fMom;
      TH1F*  fNActive;
      TH1F*  fChi2Dof;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fNTprTracks;
      TH1F*  fNCprTracks;
    };

    struct Hist_t {
      EventHist_t* fEvent[kNEventHistSets];
      TrackHist_t*  fTpr [kNTprHistSets];
      TrackHist_t*  fCpr [kNCprHistSets];
    };
  protected:

    bool                                  _mcDiag;
    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    int                                   _eventNumber;

    Data_t*                               _data;                 // diag data, passed from the caller, cached
    Hist_t                                _hist;

  public:

    MergePatRecDiag(const fhicl::ParameterSet& PSet);
    ~MergePatRecDiag();

  private:

    void        bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir);
    void        bookTrackHistograms(TrackHist_t* Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms(EventHist_t* Hist);
//-----------------------------------------------------------------------------
// overriden virtual functions of the base class
//-----------------------------------------------------------------------------
  public:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1 ) override ;
    virtual int debug         (void* Data, int Mode = -1 ) override ;
  };


//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  MergePatRecDiag::MergePatRecDiag(const fhicl::ParameterSet& PSet):
    _mcDiag                (PSet.get<bool>         ("mcDiag"                       ))
  {
    printf(" MergePatRecDiag::MergePatRecDiag : HOORAY! \n");

    if (_mcDiag != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  MergePatRecDiag::~MergePatRecDiag() {
  }

  //-----------------------------------------------------------------------------
  void MergePatRecDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>("event" , "Event Number" , 100, 0., 1000.);
    Hist->fRunNumber       = Dir->make<TH1F>("run"   , "Run   Number" , 100, 0., 100000.);
    Hist->fNTprTracks      = Dir->make<TH1F>("ntpr"  , "N(TPR) tracks",  10, 0., 10.);
    Hist->fNCprTracks      = Dir->make<TH1F>("ncpr"  , "N(CPR) tracks",  10, 0., 10.);
  }

//-----------------------------------------------------------------------------
  void MergePatRecDiag::bookTrackHistograms(TrackHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fNActive = Dir->make<TH1F>("nactv" , "N(active hits)", 200, 0., 200.);
    Hist->fMom     = Dir->make<TH1F>("mom"   , "momentum"      , 200, 0., 200.);
    Hist->fChi2Dof = Dir->make<TH1F>("chi2d" , "Chi2/Ndof"     , 200, 0., 20.);
  }

//-----------------------------------------------------------------------------
// this routine is called once per job (likely, from beginJob)
// TH1::AddDirectory makes sure one can have histograms with the same name
// in different subdirectories
//-----------------------------------------------------------------------------
  int MergePatRecDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    TH1::AddDirectory(0);
    char folder_name[20];
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fEvent[i] = new EventHist_t;
        bookEventHistograms(_hist.fEvent[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book TrkPatRec histograms
//-----------------------------------------------------------------------------
    int book_tpr_histset[kNTprHistSets];
    for (int i=0; i<kNTprHistSets; i++) book_tpr_histset[i] = 0;

    book_tpr_histset[ 0] = 1;                // all tracks

    for (int i=0; i<kNTprHistSets; i++) {
      if (book_tpr_histset[i] != 0) {
        sprintf(folder_name,"tpr_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fTpr[i] = new TrackHist_t;
        bookTrackHistograms(_hist.fTpr[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book CalPatRec histograms
//-----------------------------------------------------------------------------
    int book_cpr_histset[kNCprHistSets];
    for (int i=0; i<kNCprHistSets; i++) book_cpr_histset[i] = 0;

    book_cpr_histset[ 0] = 1;                // all tracks

    for (int i=0; i<kNCprHistSets; i++) {
      if (book_cpr_histset[i] != 0) {
        sprintf(folder_name,"cpr_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fCpr[i] = new TrackHist_t;
        bookTrackHistograms(_hist.fCpr[i],&dir);
      }
    }

    return 0;
  }


//-----------------------------------------------------------------------------
  void  MergePatRecDiag::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();
    int run_number   = _data->event->run  ();

    Hist->fEventNumber->Fill(event_number);
    Hist->fRunNumber  ->Fill(run_number);

    int ntpr = _data->list_of_kseed_tpr->size();
    Hist->fNTprTracks->Fill(ntpr);

    int ncpr = _data->list_of_kseed_cpr->size();
    Hist->fNCprTracks->Fill(ncpr);
  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// main fill histograms function called once per event
// 'Mode' not used
//-----------------------------------------------------------------------------
  int MergePatRecDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;
//-----------------------------------------------------------------------------
// start from precalculating MC-specific info
//-----------------------------------------------------------------------------
    int en = _data->event->event();
    if (_mcDiag) {
      if (_eventNumber != en) {
          _eventNumber       = en;
      }
    }
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
    return 0;
  }
//-----------------------------------------------------------------------------
// fill TrkPatRec histograms
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// fill CalPatRec histograms
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
//-----------------------------------------------------------------------------
  int MergePatRecDiag::debug(void* Data, int Mode) {
//    _data = (Data_t*) Data;
    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::MergePatRecDiag)
