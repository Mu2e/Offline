///////////////////////////////////////////////////////////////////////////////
// diag mode: = 0 - most of the histograms
//            = 1 - doca histograms
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/KalmanTrack/KalRep.hh"

#include "BTrkData/inc/TrkStrawHit.hh"

#include "RecoDataProducts/inc/Doublet.hh"

#include "CalPatRec/inc/CalTrkFit_types.hh"
#include "CalPatRec/inc/McUtilsToolBase.hh"
#include "CalPatRec/inc/ModuleHistToolBase.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {
  using namespace CalTrkFitTypes;
  
  class CalTrkFitDiag : public mu2e::ModuleHistToolBase {
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
    };

    struct DoubletHist_t {
      TH1F*  dSlope;
      TH1F*  chi2b ;
      TH1F*  chi2r ;
    };

    struct HitData_t {
      float  doca;
      float  dtCls;
    };
    
    struct HitHist_t {
      TH1F*  doca;
      TH1F*  dtCls;
    };

    struct Hist_t {
      EventHist_t*   _event  [kNEventHistSets];
      TrackHist_t*   _track  [kNTrackHistSets];
      DoubletHist_t* _doublet[kNDoubletHistSets];
      HitHist_t*     _hit    [kNHitHistSets];
    };

  protected:
    int                              _mcTruth;
    std::unique_ptr<McUtilsToolBase> _mcUtils;
    std::string                      _shDigiLabel;
    Hist_t                           _hist;              // owned 
    Data_t*                          _data;
  
  public:

    CalTrkFitDiag(const fhicl::ParameterSet& PSet);
    ~CalTrkFitDiag();

  private:

    int  bookEventHistograms  (EventHist_t*   Hist, art::TFileDirectory* Dir);
    int  bookTrackHistograms  (TrackHist_t*   Hist, art::TFileDirectory* Dir);
    int  bookDoubletHistograms(DoubletHist_t* Hist, art::TFileDirectory* Dir);
    int  bookHitHistograms    (HitHist_t*     Hist, art::TFileDirectory* Dir);
    
    int  fillEventHistograms  (EventHist_t*   Hist, Data_t*    Data);
    int  fillDoubletHistograms(DoubletHist_t* Hist, Doublet*   D);
    int  fillHitHistograms    (HitHist_t*     Hist, HitData_t* Data);
    int  fillTrackHistograms  (TrackHist_t*   Hist, KalRep*    KRep);
				  
    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };

  CalTrkFitDiag::CalTrkFitDiag(const fhicl::ParameterSet& PSet) {
    _mcTruth = PSet.get <int >("mcTruth"); 

    if (_mcTruth != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else               _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  CalTrkFitDiag::~CalTrkFitDiag() {
  }

//-----------------------------------------------------------------------------
  int CalTrkFitDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->ntracks  = Dir->make<TH1F>("ntracks", "number of track candidates: all events", 21, -0.5, 20.5);
    return 0;
  }
    
//-----------------------------------------------------------------------------
  int CalTrkFitDiag::bookDoubletHistograms(DoubletHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->dSlope = Dir->make<TH1F>("dslope", "Delta Slope",200, -1,   1);
    Hist->chi2b  = Dir->make<TH1F>("chi2b" , "chi2(best)",200,  0, 200);
    Hist->chi2r  = Dir->make<TH1F>("chi2r" , "chi2R (best/next)  OS",200,  0, 1);
    
    return 0;
  }
    
//-----------------------------------------------------------------------------
  int CalTrkFitDiag::bookTrackHistograms(TrackHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->nhits    = Dir->make<TH1F>("nhits", "N(track hits)"    , 101,  -0.5, 100.5);
    Hist->chi2dof  = Dir->make<TH1F>("chi2d", "track chi2/dof"   , 100,   0. ,  10.);
    Hist->p        = Dir->make<TH1F>("p"    , "track momentum"   , 400,   0. , 200.);
    
    return 0;
  }
    
//-----------------------------------------------------------------------------
  int CalTrkFitDiag::bookHitHistograms(HitHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->doca  = Dir->make<TH1F>("doca" ,"doca", 1000, -20.,  20 );
    Hist->dtCls = Dir->make<TH1F>("dtCls","dtCls; #Delta t = t_{calo-cluster}-t_{straw} - tof [ns]", 401, -200.5,  200.5 );
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTrkFitDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
    char folder_name[20];
    
    TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;		// all events

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

    book_track_histset[ 0] = 1;		// all 
    book_track_histset[ 1] = 1;		// nhits > 15

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

    book_doublet_histset[ 0] = 1;		// OS 
    book_doublet_histset[ 1] = 1;		// OS + MC truth
    book_doublet_histset[ 2] = 1;		// SS
    book_doublet_histset[ 3] = 1;		// SS + MC truth

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

    book_hit_histset[ 0] = 1;		// active hits
    book_hit_histset[ 1] = 1;		// non-active hits

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
  int CalTrkFitDiag::fillEventHistograms(EventHist_t* Hist, Data_t* Data) {
    int ntrk = Data->tracks->size();
    Hist->ntracks->Fill(ntrk);
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTrkFitDiag::fillTrackHistograms(TrackHist_t* Hist, KalRep* KRep) {
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTrkFitDiag::fillDoubletHistograms(DoubletHist_t* Hist, Doublet* D) {
    
    float chi2b = D->Chi2Best();
    float chi2r = chi2b/D->fChi2[D->fINext];

    Hist->dSlope->Fill(D->bestDxDzRes());
    Hist->chi2b ->Fill(chi2b);
    Hist->chi2r ->Fill(chi2r);
    
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTrkFitDiag::fillHitHistograms(HitHist_t* Hist, HitData_t* HitData) {
    Hist->doca ->Fill(HitData->doca );
    Hist->dtCls->Fill(HitData->dtCls);    
    return 0;
  }

//-----------------------------------------------------------------------------
// mode not used
//----------------------------------------------------------------------------- 
  int CalTrkFitDiag::fillHistograms(void* Data, int Mode) {
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
// reconstructions
//-----------------------------------------------------------------------------
      CLHEP::Hep3Vector        tdir;
      HepPoint                 tpos;

      const KalRep* krep = &_data->tracks->at(i);
					          // hits on the track
      
      TrkHitVector const& hot_l = krep->hitVector();

      krep->traj().getInfo(0.0,tpos,tdir);
					          // loop over track hits
      int nhits = krep->hitVector().size();

      //get the calorimeter cluster from the KalSeed
      double             z_cls(-99999), time_cls(-9999);
      double             pitchAngle(0.67);//FIX ME! that should be parsed from the fcl
      double             meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity
      const CaloCluster* cluster(0);
      cluster = _data->kscol->at(i).caloCluster().get();

      if (cluster != 0)  {
	CLHEP::Hep3Vector gpos        = _data->calorimeter->geomUtil().diskToMu2e(cluster->diskId(),cluster->cog3Vector());
	CLHEP::Hep3Vector cog_cluster = _data->calorimeter->geomUtil().mu2eToTracker(gpos);
	z_cls    = cog_cluster.z();//z-coordinate of the cluster in the tracker coordinate frame
	
	time_cls = cluster->time();
      }

      for (int i=0; i<nhits; ++i) {
	const mu2e::TrkStrawHit* hit = static_cast<TrkStrawHit*> (hot_l.at(i));
	int                hIndex    = hit->index();
	StrawHit const*    sh        = & _data->result->_shcol->at(hIndex);
	Straw const&       straw     = _data->tracker->getStraw(sh->strawIndex());
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
	  
	// convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	
	HepPoint          spt(hpos.x(),hpos.y(),hpos.z());
	TrkLineTraj       htraj(spt,hdir,-20,20);

	// estimate flightlength along track.  This assumes a constant BField!!!
	
	double           fltlen = (hpos.z()-tpos.z())/tdir.z();
	TrkPoca          hitpoca(krep->traj(),fltlen,htraj,0.0);

	hitData.doca = hitpoca.doca();
	
	//estiamte the time residual between the straw-hit and the calcorimeter cluster taking into account the tof and mean-rift time
	double           dt(-9999.);
	if (cluster != 0){
	  double           z_straw = hpos.z();
	  double           time    = sh->time();
	  double           tof     = (z_cls - z_straw)/sin(pitchAngle)/CLHEP::c_light;
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
      for (int i=0; i<nd; i++) {
	d = &_data->listOfDoublets->at(i);
	if (d->fNStrawHits == 2) {
	  
	  int same_sign = d->isSameSign();
	  int a0        = d->fHit[0]->ambig();
	  int a1        = d->fHit[1]->ambig();
	  
	  bool h1_ok = ((a0 != 0) && (a0*d->fMcDoca[0] > 0));
	  bool h2_ok = ((a1 != 0) && (a1*d->fMcDoca[1] > 0));
	  
	  if   (same_sign) {
	    fillDoubletHistograms(_hist._doublet[2],d);
	    if (h1_ok && h2_ok) fillDoubletHistograms(_hist._doublet[3],d);
	  }
	  else {
	    fillDoubletHistograms(_hist._doublet[0],d);
	    if (h1_ok && h2_ok) fillDoubletHistograms(_hist._doublet[0],d);
	  }
	}
      }
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(CalTrkFitDiag)

}

