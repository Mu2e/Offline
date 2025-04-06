//
/*
#include "Offline/TrkPatRec/inc/KalSeedFit_types.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/TrkReco/inc/KalFitData.hh"

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1F.h"

namespace mu2e {

  using namespace KalSeedFitTypes;

  class KalSeedFitDiag : public mu2e::ModuleHistToolBase {
    public:

      enum {
        kNEventHistSets = 10,
        kNTrackHistSets = 10,
        kNHitHistSets   = 10
      };

      struct KalSeedData_t {
        int   nactive;
        int   ndeactivated;
        int   nrescued;
        float mom;

        KalSeedData_t(): nactive(0), ndeactivated(0), nrescued(0), mom(-1.) {}
      };

      struct EventHist_t {
        TH1F*  ntracks;                   //
      };

      struct TrackHist_t {
        TH1F*  nha;           // number of active hits
        TH1F*  nhd;           // number of deactivated hits
        TH1F*  nhr;           // number of rescued hits
        TH1F*  chi2dof;
        TH1F*  p;
      };

      struct HitHist_t {
        TH1F*  doca;            // distance of closest approach
      };

      struct Hist_t {
        EventHist_t* event[kNEventHistSets];
        TrackHist_t* track[kNTrackHistSets];
        HitHist_t*   hit  [kNHitHistSets  ];
      };

    protected:
      int                              _mcTruth;
      std::unique_ptr<McUtilsToolBase> _mcUtils;
      std::string                      _shDigiLabel;
      Hist_t                           _hist;            // owned
      Data_t*                          _data;            // cached

      //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    public:

      KalSeedFitDiag(const fhicl::ParameterSet& PSet);
      ~KalSeedFitDiag();

    private:

      int bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir);
      int bookTrackHistograms(TrackHist_t* Hist, art::TFileDirectory* Dir);
      int bookHitHistograms  (HitHist_t*   Hist, art::TFileDirectory* Dir);

      int fillEventHistograms(EventHist_t* Hist, Data_t* Data);
      int fillTrackHistograms(TrackHist_t* Hist, KalSeed* Track, KalSeedData_t* KsData);

      virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
      virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };


  //-----------------------------------------------------------------------------
  KalSeedFitDiag::KalSeedFitDiag(const fhicl::ParameterSet& PSet) {
    _mcTruth = PSet.get <int >("mcTruth");

    if (_mcTruth != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else               _mcUtils = std::make_unique<McUtilsToolBase>();
  }

  //-----------------------------------------------------------------------------
  KalSeedFitDiag::~KalSeedFitDiag() {
  }

  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->ntracks  = Dir->make<TH1F>("ntracks", "number of track candidates: all events", 21, -0.5, 20.5);
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::bookTrackHistograms(TrackHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->nha     = Dir->make<TH1F>("nha"  , "N(track hits)"            , 101,  -0.5, 100.5);
    Hist->nhd     = Dir->make<TH1F>("nhd"  , "N(deactivated track hits)", 101,  -0.5, 100.5);
    Hist->nhr     = Dir->make<TH1F>("nhr"  , "N(rescued track hits)"    , 101,  -0.5, 100.5);
    Hist->chi2dof = Dir->make<TH1F>("chi2d", "track chi2/dof"           , 100,   0. ,  10.);
    Hist->p       = Dir->make<TH1F>("p"    , "track momentum"           , 400,   0. , 200.);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::bookHitHistograms(HitHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->doca = Dir->make<TH1F>("doca","doca, [mm]", 1000, -20., 20);
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
    char folder_name[20];

    TH1::AddDirectory(0);
    //-----------------------------------------------------------------------------
    // book event-level histograms - fill them once per event
    //-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;   // all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.event[i] = new EventHist_t;
        bookEventHistograms(_hist.event[i],&dir);
      }
    }
    //-----------------------------------------------------------------------------
    // book track histograms
    //-----------------------------------------------------------------------------
    int book_track_histset[kNTrackHistSets];
    for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

    book_track_histset[ 0] = 1;   // all
    book_track_histset[ 1] = 1;   // nhits > 15

    for (int i=0; i<kNTrackHistSets; i++) {
      if (book_track_histset[i] != 0) {
        sprintf(folder_name,"trk_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.track[i] = new TrackHist_t;
        bookTrackHistograms(_hist.track[i],&dir);
      }
    }
    //-----------------------------------------------------------------------------
    // book hit histograms
    //-----------------------------------------------------------------------------
    int book_hit_histset[kNHitHistSets];
    for (int i=0; i<kNHitHistSets; i++) book_hit_histset[i] = 0;

    book_hit_histset[ 0] = 1;   // active hits
    book_hit_histset[ 1] = 1;   // non-active hits

    for (int i=0; i<kNHitHistSets; i++) {
      if (book_hit_histset[i] != 0) {
        sprintf(folder_name,"hit_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.hit[i] = new HitHist_t;
        bookHitHistograms(_hist.hit[i],&dir);
      }
    }

    return 0;
  }


  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::fillEventHistograms(EventHist_t* Hist, Data_t* Data) {
    int ntrk = Data->tracks->size();
    Hist->ntracks->Fill(ntrk);
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::fillTrackHistograms(TrackHist_t* Hist, KalSeed* Track, KalSeedData_t* KsData) {
    Hist->nha->Fill(KsData->nactive);
    Hist->nhd->Fill(KsData->ndeactivated);
    Hist->nhr->Fill(KsData->ndeactivated);

    double  chi2dof  = Track->chisquared()/(KsData->nactive-5);

    Hist->chi2dof->Fill(chi2dof);

    return 0;
  }

  //-----------------------------------------------------------------------------
  // Mode is not used
  //-----------------------------------------------------------------------------
  int KalSeedFitDiag::fillHistograms(void* Data, int Mode) {

    _data = (Data_t*) Data;

    //-----------------------------------------------------------------------------
    // fill event-level histograms
    //-----------------------------------------------------------------------------
    fillEventHistograms(_hist.event[0],_data);
    //-----------------------------------------------------------------------------
    // fill track-level histograms
    //-----------------------------------------------------------------------------
    int ntrk = _data->tracks->size();
    for (int i=0; i<ntrk; i++) {
      KalSeed* kseed = &_data->tracks->at(i);

      KalSeedData_t ksdata;

      vector<TrkStrawHitSeed> const& hot_l = kseed->hits();
      for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
        const TrkStrawHitSeed* hit = static_cast<const TrkStrawHitSeed*> (&(*it));
        StrawHitFlag flag = hit->flag();
        bool hitIsActive = flag.hasAllProperties(StrawHitFlagDetail::active);
        if (hitIsActive) ++ksdata.nactive;
        else             ++ksdata.ndeactivated;
      }

      ksdata.nrescued = _data->nrescued.at(i);
      ksdata.mom      = _data->mom.at(i);

      fillTrackHistograms(_hist.track[0],kseed,&ksdata);
      if (ksdata.nactive > 15) fillTrackHistograms(_hist.track[1],kseed,&ksdata);
    }
    return 0;
  }

  DEFINE_ART_CLASS_TOOL(KalSeedFitDiag)

}
*/
