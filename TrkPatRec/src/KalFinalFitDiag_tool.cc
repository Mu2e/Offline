///////////////////////////////////////////////////////////////////////////////
// diag mode: = 0 - most of the histograms
//            = 1 - doca histograms
///////////////////////////////////////////////////////////////////////////////
/*
#include "TH2.h"
#include "TH1.h"

#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/KalmanTrack/KalRep.hh"

#include "Offline/BTrkData/inc/TrkStrawHit.hh"

#include "Offline/BTrkData/inc/Doublet.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "Offline/TrkPatRec/inc/KalFinalFit_types.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/TrkReco/inc/KalFitData.hh"
#include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {
  using namespace KalFinalFitTypes;

  class KalFinalFitDiag : public mu2e::ModuleHistToolBase {
    public:
      enum {
        kNEventHistSets   = 10,
        kNTrackHistSets   = 20,
        kNDoubletHistSets = 10,
        kNHitHistSets     = 10
      };

      struct EventHist_t {
        TH1F*  ntracks;
      };

      struct TrackHist_t {
        TH1F*  nhits  ;
        TH1F*  chi2dof;
        TH1F*  p      ;
        TH1F*  tchDiskId[2];
        TH1F*  tchAdded    ;
        TH1F*  tchDepth [2];
        TH1F*  tchDoca  [2];
        TH1F*  tchDt    [2];
        TH1F*  trkPath  [2];
        TH1F*  tchEnergy[2];
        TH1F*  tchEp    [2];
      };

      struct DoubletHist_t {
        TH1F*  dSlope;
        TH1F*  chi2bc;   // coordinate part of the best chi2
        TH1F*  chi2bs;   // angular part of the best chi2
        TH1F*  chi2b ;   // best chi2
        TH1F*  chi2r ;   // chi2_best/chi2_next
      };

      struct HitData_t {
        float  doca;
        float  mcdoca;
        float  rdrift;
        float  dtCls;
      };

      struct HitHist_t {
        TH1F*  doca;
        TH1F*  xdoca;                // doca/hit_error
        TH2F*  doca_vs_mcdoca;
        TH2F*  rdrift_vs_mcdoca;
        TH1F*  dtCls;
      };

      struct Hist_t {
        EventHist_t*   _event  [kNEventHistSets  ];
        TrackHist_t*   _track  [kNTrackHistSets  ];
        DoubletHist_t* _doublet[kNDoubletHistSets];
        HitHist_t*     _hit    [kNHitHistSets    ];
      };

    protected:
      int                              _mcTruth;
      std::unique_ptr<McUtilsToolBase> _mcUtils;
      std::string                      _shDigiLabel;
      Hist_t                           _hist;              // owned
      Data_t*                          _data;

    public:

      KalFinalFitDiag(const fhicl::ParameterSet& PSet);
      ~KalFinalFitDiag();

    private:

      int  bookEventHistograms  (EventHist_t*   Hist, art::TFileDirectory* Dir);
      int  bookTrackHistograms  (TrackHist_t*   Hist, art::TFileDirectory* Dir);
      int  bookDoubletHistograms(DoubletHist_t* Hist, art::TFileDirectory* Dir);
      int  bookHitHistograms    (HitHist_t*     Hist, art::TFileDirectory* Dir);

      int  fillEventHistograms  (EventHist_t*   Hist, Data_t*    Data);
      int  fillDoubletHistograms(DoubletHist_t* Hist, Doublet*   D);
      int  fillHitHistograms    (HitHist_t*     Hist, HitData_t* Data);
      int  fillTrackHistograms  (TrackHist_t*   Hist, Data_t*    Data, const KalRep*krep);

      virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
      virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };

  KalFinalFitDiag::KalFinalFitDiag(const fhicl::ParameterSet& PSet) {
    _mcTruth = PSet.get <int >("mcTruth");

    if (_mcTruth != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else               _mcUtils = std::make_unique<McUtilsToolBase>();
  }

  //-----------------------------------------------------------------------------
  KalFinalFitDiag::~KalFinalFitDiag() {
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->ntracks  = Dir->make<TH1F>("ntracks", "number of track candidates: all events", 21, -0.5, 20.5);
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::bookDoubletHistograms(DoubletHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->dSlope = Dir->make<TH1F>("dslope", "Delta Slope"          , 200, -1,   1);
    Hist->chi2b  = Dir->make<TH1F>("chi2b" , "chi2(best)"           , 200,  0, 200);
    Hist->chi2bc = Dir->make<TH1F>("chi2bc", "chi2_coord(best)"     , 200,  0, 200);
    Hist->chi2bs = Dir->make<TH1F>("chi2bs", "chi2_slope(best)"     , 200,  0, 200);
    Hist->chi2r  = Dir->make<TH1F>("chi2r" , "chi2R (best/next)  OS", 200,  0,   1);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::bookTrackHistograms(TrackHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->nhits    = Dir->make<TH1F>("nhits", "N(track hits)"    , 101,  -0.5, 100.5);
    Hist->chi2dof  = Dir->make<TH1F>("chi2d", "track chi2/dof"   , 100,   0. ,  10.);
    Hist->p        = Dir->make<TH1F>("p"    , "track momentum"   , 400,   0. , 200.);
    //histogromas for the TrackCaloHit
    Hist->tchAdded = Dir->make<TH1F>("tchAdded "    , "trackCaloHit: added in the end: added-later"   , 2,   0. , 2.);

    Hist->tchDiskId[0] = Dir->make<TH1F>("tchDiskId0"    , "trackCaloHit: diskId; diskId"   , 3,   0. , 3.);
    Hist->tchDepth [0] = Dir->make<TH1F>("tchDepth0 "    , "trackCaloHit: hit-length"   , 500,   0. , 500.);
    Hist->tchDoca  [0] = Dir->make<TH1F>("tchDoca0  "    , "trackCaloHit: DOCA; DOCA [mm]"   , 400,   -200. , 200.);
    Hist->tchDt    [0] = Dir->make<TH1F>("tchDt0    "    , "trackCaloHit: dt @ calo; dt [0] = t_{calo} - t_{trk} [ns]"   , 400,   -10. , 10.);
    Hist->trkPath  [0] = Dir->make<TH1F>("trkPath0  "    , "trackCaloHit: track path length; trk-calo path [mm]"   , 500,   0. , 500.);
    Hist->tchEnergy[0] = Dir->make<TH1F>("tchEnergy0"    , "trackCaloHit: energy; E [MeV]"   , 240,   0. , 120.);
    Hist->tchEp    [0] = Dir->make<TH1F>("tchEp0"        , "trackCaloHit: E/p; E/p"          , 240,   0. , 1.2);

    Hist->tchDiskId[1] = Dir->make<TH1F>("tchDiskId1"    , "trackCaloHit added lastly: diskId; diskId"   , 3,   0. , 3.);
    Hist->tchDepth [1] = Dir->make<TH1F>("tchDepth1 "    , "trackCaloHit added lastly: hit-length"   , 500,   0. , 500.);
    Hist->tchDoca  [1] = Dir->make<TH1F>("tchDoca1  "    , "trackCaloHit added lastly: DOCA; DOCA [mm]"   , 400,   -200. , 200.);
    Hist->tchDt    [1] = Dir->make<TH1F>("tchDt1    "    , "trackCaloHit added lastly: dt @ calo; dt [1] = t_{calo} - t_{trk} [ns]"   , 400,   -10. , 10.);
    Hist->trkPath  [1] = Dir->make<TH1F>("trkPath1  "    , "trackCaloHit added lastly: track path length; trk-calo path [mm]"   , 500,   0. , 500.);
    Hist->tchEnergy[1] = Dir->make<TH1F>("tchEnergy1"    , "trackCaloHit added lastly: energy; E [MeV]"   , 240,   0. , 120.);
    Hist->tchEp    [1] = Dir->make<TH1F>("tchEp1"        , "trackCaloHit added lastly: E/p; E/p"          , 240,   0. , 1.2);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::bookHitHistograms(HitHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->doca  = Dir->make<TH1F>("doca"  ,"doca" , 1000, -20.,  20 );
    Hist->xdoca = Dir->make<TH1F>("xdoca" ,"xdoca", 1000, -20.,  20 );
    Hist->dtCls = Dir->make<TH1F>("dtCls" ,"dtCls; #Delta t = t_{calo-cluster}-t_{straw} - tof [ns]", 401, -200.5,  200.5 );
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
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

        _hist._event[i] = new EventHist_t;
        bookEventHistograms(_hist._event[i],&dir);
      }
    }
    //-----------------------------------------------------------------------------
    // book track histograms
    //-----------------------------------------------------------------------------
    int book_track_histset[kNTrackHistSets];
    for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

    book_track_histset[ 0] = 1;   // all
    book_track_histset[ 1] = 1;   // nhits > 20
    book_track_histset[ 2] = 1;   // nhits > 20 & p > 100

    for (int i=0; i<kNTrackHistSets; i++) {
      if (book_track_histset[i] != 0) {
        sprintf(folder_name,"trk_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist._track[i] = new TrackHist_t;
        bookTrackHistograms(_hist._track[i],&dir);
      }
    }
    //-----------------------------------------------------------------------------
    // book doublet histograms
    //-----------------------------------------------------------------------------
    int book_doublet_histset[kNDoubletHistSets];
    for (int i=0; i<kNDoubletHistSets; i++) book_doublet_histset[i] = 0;

    book_doublet_histset[ 0] = 1;   // OS
    book_doublet_histset[ 1] = 1;   // OS, correct assignment of the drift directions
    book_doublet_histset[ 2] = 1;   // SS
    book_doublet_histset[ 3] = 1;   // SS, correct assignment of the drift directions

    for (int i=0; i<kNDoubletHistSets; i++) {
      if (book_doublet_histset[i] != 0) {
        sprintf(folder_name,"dbl_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist._doublet[i] = new DoubletHist_t;
        bookDoubletHistograms(_hist._doublet[i],&dir);
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

        _hist._hit[i] = new HitHist_t;
        bookHitHistograms(_hist._hit[i],&dir);
      }
    }

    return 0;
  }


  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::fillEventHistograms(EventHist_t* Hist, Data_t* Data) {
    int ntrk = Data->tracks->size();
    Hist->ntracks->Fill(ntrk);
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::fillTrackHistograms(TrackHist_t* Hist, Data_t* Data, const KalRep*krep) {

    const TrkHitVector* thv = &(krep->hitVector());
    double  h1_fltlen, hn_fltlen, entlen;

    h1_fltlen      = krep->firstHit()->kalHit()->hit()->fltLen();
    hn_fltlen      = krep->lastHit ()->kalHit()->hit()->fltLen();
    entlen         = std::min(h1_fltlen,hn_fltlen);

    Hist->nhits     ->Fill( thv->size());
    Hist->chi2dof   ->Fill( krep->chisq()/(krep->nActive()-5.));
    Hist->p         ->Fill( krep->momentum(entlen).mag());

    if (Data->tchEnergy > 0){
      int  tchAdd = Data->tchAdded;
      Hist->tchAdded          ->Fill( tchAdd);
      Hist->tchDiskId[tchAdd] ->Fill( Data->tchDiskId);
      Hist->tchDepth [tchAdd] ->Fill( Data->tchDepth);
      Hist->tchDoca  [tchAdd] ->Fill( Data->tchDOCA);
      Hist->tchDt    [tchAdd] ->Fill( Data->tchDt);
      Hist->trkPath  [tchAdd] ->Fill( Data->tchTrkPath);
      Hist->tchEnergy[tchAdd] ->Fill( Data->tchEnergy);
      Hist->tchEp    [tchAdd] ->Fill( Data->tchEnergy/krep->momentum(entlen).mag());
    }
    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::fillDoubletHistograms(DoubletHist_t* Hist, Doublet* D) {

    float chi2b  = D->chi2Best      ();
    float chi2bc = D->chi2CoordBest();
    float chi2bs = D->chi2SlopeBest();
    float chi2r  = chi2b/D->fChi2[D->fINext];

    Hist->dSlope->Fill(D->bestDxDzRes());
    Hist->chi2b ->Fill(chi2b  );
    Hist->chi2bc->Fill(chi2bc);
    Hist->chi2bs->Fill(chi2bs);

    Hist->chi2r ->Fill(chi2r);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::fillHitHistograms(HitHist_t* Hist, HitData_t* HitData) {
    Hist->doca ->Fill(HitData->doca );
    Hist->dtCls->Fill(HitData->dtCls);
    return 0;
  }

  //-----------------------------------------------------------------------------
  // mode not used
  //-----------------------------------------------------------------------------
  int KalFinalFitDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;

    HitData_t hitData;
    //-----------------------------------------------------------------------------
    // fill event histograms
    //-----------------------------------------------------------------------------
    fillEventHistograms(_hist._event[0],_data);
    //-----------------------------------------------------------------------------
    // fill track histograms - per-track, can't have the loop here - some are
    // not reconstructed
    //-----------------------------------------------------------------------------
    int ntrk = _data->tracks->size();
    for (int i=0; i<ntrk; i++) {
      //-----------------------------------------------------------------------------
      // fill track-level histograms, a track doesn't necessarily have to be successfully
      // reconstructed
      //-----------------------------------------------------------------------------
      CLHEP::Hep3Vector        tdir;
      HepPoint                 tpos;

      const KalRep* krep = &_data->tracks->at(i);

      fillTrackHistograms(_hist._track[0], _data, krep);
      // hits on the track

      TrkHitVector const& hot_l = krep->hitVector();

      krep->traj().getInfo(0.0,tpos,tdir);
      // loop over track hits
      int nhits = krep->hitVector().size();

      if (nhits > 20)       fillTrackHistograms(_hist._track[1], _data, krep);
      double  h1_fltlen, hn_fltlen, entlen;

      h1_fltlen      = krep->firstHit()->kalHit()->hit()->fltLen();
      hn_fltlen      = krep->lastHit ()->kalHit()->hit()->fltLen();
      entlen         = std::min(h1_fltlen,hn_fltlen);
      double  mom    = krep->momentum(entlen).mag();

      if (nhits > 20 && mom>100.)    fillTrackHistograms(_hist._track[2], _data, krep);

      //-----------------------------------------------------------------------------
      // get the calorimeter cluster from the KalSeed
      //-----------------------------------------------------------------------------
      double             z_cls(-99999), time_cls(-9999);
      double             pitchAngle(0.67);//FIX ME! that should be parsed from the fcl
      double             meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity

      const CaloCluster* cluster = _data->kscol->at(i).caloCluster().get();

      if (cluster != 0)  {
        CLHEP::Hep3Vector gpos        = _data->calorimeter->geomUtil().diskToMu2e(cluster->diskID(),cluster->cog3Vector());
        CLHEP::Hep3Vector cog_cluster = _data->calorimeter->geomUtil().mu2eToTracker(gpos);

        z_cls    = cog_cluster.z(); // z-coordinate of the cluster in the tracker coordinate frame
        time_cls = cluster->time();
      }

      for (int i=0; i<nhits; ++i) {
        const mu2e::TrkStrawHit* hit = static_cast<TrkStrawHit*> (hot_l.at(i));
        int                hIndex    = hit->index();
        ComboHit const*    sh        = & _data->result->chcol->at(hIndex);
        Straw const&       straw     = _data->tracker->getStraw(sh->strawId());
        const CLHEP::Hep3Vector& hpos = straw.getMidPoint();
        const CLHEP::Hep3Vector& hdir = straw.getDirection();

        bool               found(false);

        for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
          hit = static_cast<const mu2e::TrkStrawHit*> (*it);
          if (!hit->isActive()) continue;
          int hit_index = hit->index();
          if (hIndex == hit_index) {
            found = true;
            break;
          }
        }
        //-----------------------------------------------------------------------------
        // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
        // estimate flight length along track.  This assumes a constant BField!!!
        //-----------------------------------------------------------------------------
        HepPoint          spt(hpos.x(),hpos.y(),hpos.z());
        TrkLineTraj       htraj(spt,hdir,-20,20);

        double  fltlen = (hpos.z()-tpos.z())/tdir.z();
        TrkPoca hitpoca(krep->traj(),fltlen,htraj,0.0);

        hitData.doca = hitpoca.doca();
        //-----------------------------------------------------------------------------
        // estimate the time residual between the straw hit and the calorimeter cluster
        // taking into account the hit TOF and mean drift time
        //-----------------------------------------------------------------------------
        double           dt(-9999.), z_straw, time, tof;
        if (cluster != 0){
          z_straw = hpos.z();
          time    = sh->time();
          tof     = (z_cls - z_straw)/sin(pitchAngle)/CLHEP::c_light;
          dt      = time_cls - (time + tof - meanDriftTime);
        }

        hitData.dtCls = dt;

        if (found) fillHitHistograms(_hist._hit[0],&hitData);
        else       fillHitHistograms(_hist._hit[1],&hitData);
      }
      //-----------------------------------------------------------------------------
      // doublet histograms
      //-----------------------------------------------------------------------------
      _data->dar->findDoublets(krep,_data->listOfDoublets);

      Doublet* d;
      int nd = _data->listOfDoublets->size();

      double mcdoca[2];
      //-----------------------------------------------------------------------------
      // doublet stores TrkStrawHit's
      //-----------------------------------------------------------------------------
      for (int i=0; i<nd; i++) {
        d = &_data->listOfDoublets->at(i);
        if (d->fNStrawHits == 2) {

          int same_sign = d->isSameSign();
          int a0        = d->fHit[0]->ambig();
          int a1        = d->fHit[1]->ambig();

          mcdoca[0]     = _mcUtils->mcDoca(_data->event,d->fHit[0]);
          mcdoca[1]     = _mcUtils->mcDoca(_data->event,d->fHit[1]);

          // bool h1_ok = ((a0 != 0) && (a0*d->fMcDoca[0] > 0));
          // bool h2_ok = ((a1 != 0) && (a1*d->fMcDoca[1] > 0));

          bool h1_ok = ((a0 != 0) && (a0*mcdoca[0] > 0));
          bool h2_ok = ((a1 != 0) && (a1*mcdoca[1] > 0));

          if   (! same_sign) {
            // OS doublet
            fillDoubletHistograms(_hist._doublet[0],d);
            if (h1_ok && h2_ok) fillDoubletHistograms(_hist._doublet[1],d);
          }
          else {
            // SS doublet
            fillDoubletHistograms(_hist._doublet[2],d);
            if (h1_ok && h2_ok) fillDoubletHistograms(_hist._doublet[3],d);
          }
        }
      }
    }
    return 0;
  }

  DEFINE_ART_CLASS_TOOL(KalFinalFitDiag)

}
*/
