///////////////////////////////////////////////////////////////////////////////
// diag mode: = 0 - most of the histograms
//            = 1 - doca histograms
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/BTrkData/inc/TrkStrawHit.hh"

#include "Offline/CalPatRec/inc/CalTimePeakFinder_types.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
// #include "CalPatRec/inc/KalFitResultNew.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "TH1.h"
#include "TH2.h"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

  using namespace CalTimePeakFinderTypes;

  class CalTimePeakFinderDiag : public mu2e::ModuleHistToolBase {
  public:

    enum {
      kNEventHistSets       = 10,
      kNClusterHistSets     = 10,
      kNTimeClusterHistSets = 10
    };

    struct EventHist_t {
      TH1F*   ncl;
      TH1F*   ncl50;
      TH1F*   nseeds;
      TH1F*   dt;
    };

    struct ClusterHist_t {
      TH1F*   energy;
      TH1F*   time;
      TH1F*   size;
      TH1F*   diskID;
    };

    struct TimeClusterHist_t {
      TH1F*  nhits;           // number of hits on a helix
      TH1F*  energy;
      TH1F*  time;
      TH2F*  time_vs_nhits;
      TH2F*  energy_vs_nhits;
    };

    struct Hist_t {
      EventHist_t*        _event  [kNEventHistSets];
      ClusterHist_t*      _cluster[kNClusterHistSets];
      TimeClusterHist_t*  _tcl    [kNTimeClusterHistSets];
    };

  protected:
    std::string  _shDigiLabel;

    std::unique_ptr<McUtilsToolBase> _mcUtils;

    Hist_t        _hist;
    Data_t*       _data;

  public:

    CalTimePeakFinderDiag(const fhicl::Table<mu2e::CalTimePeakFinderTypes::Config>& Conf);
    ~CalTimePeakFinderDiag();

  private:

    int  bookEventHistograms      (EventHist_t*       Hist, art::TFileDirectory* Dir);
    int  bookClusterHistograms    (ClusterHist_t*     Hist, art::TFileDirectory* Dir);
    int  bookTimeClusterHistograms(TimeClusterHist_t* Hist, art::TFileDirectory* Dir);

    int  fillEventHistograms      (EventHist_t*       Hist, Data_t* Data);
    int  fillClusterHistograms    (ClusterHist_t*     Hist, const CaloCluster* Cl);
    int  fillTimeClusterHistograms(TimeClusterHist_t* Hist, TimeCluster* Tcl);

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };


//-----------------------------------------------------------------------------
  CalTimePeakFinderDiag::CalTimePeakFinderDiag(const fhicl::Table<mu2e::CalTimePeakFinderTypes::Config>& Conf) {
  }

//-----------------------------------------------------------------------------
  CalTimePeakFinderDiag::~CalTimePeakFinderDiag() {
  }


//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->nseeds = Dir->make<TH1F>("nseeds0", "number of track candidates"  , 21, -0.5, 20.5);
    Hist->ncl    = Dir->make<TH1F>("ncl"    , "N(calorimeter clusters)"     , 500, 0, 500);
    Hist->ncl50  = Dir->make<TH1F>("ncl50"  , "N(calorimeter clusters) E>50", 100, 0, 100);
    Hist->dt     = Dir->make<TH1F>("dt"     , "hit dt; #Delta t = t_{cl} - (t_{ch}+TOF-<t_{drift}>) [ns]", 200, -100, 100);
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::bookClusterHistograms(ClusterHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->energy = Dir->make<TH1F>("energy", "Cluster energy", 200, 0,  200);
    Hist->time   = Dir->make<TH1F>("time"  , "Cluster time"  , 400, 0, 2000);
    Hist->size   = Dir->make<TH1F>("size"  , "Cluster size"  ,  50, 0,   50);
    Hist->diskID = Dir->make<TH1F>("disk"  , "Cluster diskID",   5, 0,    5);
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::bookTimeClusterHistograms(TimeClusterHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->nhits  = Dir->make<TH1F>("nhits" , "number of hits", 101, -0.5, 100.5);
    Hist->energy = Dir->make<TH1F>("energy", "cluster energy"       , 400, 0., 200.);
    Hist->time   = Dir->make<TH1F>("time0" , "cluster time"         , 1400, 300., 1700);

    Hist->time_vs_nhits   = Dir->make<TH2F>("time_vs_nhits"  ,"time vs nhits  ", 100, 0, 100, 2800, 300., 1700);
    Hist->energy_vs_nhits = Dir->make<TH2F>("energy_vs_nhits","energy vs nhits", 100, 0, 100, 400, 0, 200);
    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
    char folder_name[20];

    TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // events: all
    book_event_histset[ 1] = 1;                // events: with good seeds

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory tfdir = Tfs->mkdir(folder_name);

        _hist._event[i] = new EventHist_t;
        bookEventHistograms(_hist._event[i],&tfdir);
      }
    }
//-----------------------------------------------------------------------------
// book calorimeter cluster histograms
//-----------------------------------------------------------------------------
    int book_cluster_histset[kNClusterHistSets];
    for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

    book_cluster_histset[ 0] = 1;                // clusters: all
    book_cluster_histset[ 1] = 1;                // clusters: above time seeding threshold

    for (int i=0; i<kNClusterHistSets; i++) {
      if (book_cluster_histset[i] != 0) {
        sprintf(folder_name,"cl_%i",i);
        art::TFileDirectory tfdir = Tfs->mkdir(folder_name);

        _hist._cluster[i] = new ClusterHist_t;
        bookClusterHistograms(_hist._cluster[i],&tfdir);
      }
    }
//-----------------------------------------------------------------------------
// book time cluster histograms
//-----------------------------------------------------------------------------
    int book_timecluster_histset[kNTimeClusterHistSets];
    for (int i=0; i<kNTimeClusterHistSets; i++) book_timecluster_histset[i] = 0;

    book_timecluster_histset[ 0] = 1;                // all
    book_timecluster_histset[ 1] = 1;                // all

    for (int i=0; i<kNTimeClusterHistSets; i++) {
      if (book_timecluster_histset[i] != 0) {
        sprintf(folder_name,"tcl_%i",i);
        art::TFileDirectory tfdir = Tfs->mkdir(folder_name);

        _hist._tcl[i] = new TimeClusterHist_t;
        bookTimeClusterHistograms(_hist._tcl[i],&tfdir);
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::fillEventHistograms(EventHist_t* Hist, Data_t* Data) {
    int nseeds = Data->_tcColl->size();
    int ncl    = Data->ccCollection->size();
    int ncl50(0);

    for (int i=0; i<ncl; i++) {
      const CaloCluster* cl = &Data->ccCollection->at(i);
      if (cl->energyDep() > Data->minClusterEnergy) ncl50 += 1;
    }

    Hist->nseeds->Fill(nseeds);
    Hist->ncl->Fill(ncl);
    Hist->ncl50->Fill(ncl50);


    int nh = Data->dtvec.size();
    for (int i=0; i<nh; i++) {
      float dt = Data->dtvec[i];
      Hist->dt->Fill(dt);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::fillClusterHistograms(ClusterHist_t* Hist, const CaloCluster* Cl) {

    Hist->energy->Fill(Cl->energyDep());
    Hist->time->Fill(Cl->time());
    Hist->size->Fill(Cl->size());
    Hist->diskID->Fill(Cl->diskID());

    return 0;
  }

//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::fillTimeClusterHistograms(TimeClusterHist_t* Hist, TimeCluster* Tcl) {

    const CaloCluster* cl = Tcl->caloCluster().get();

    float  time   = cl->time();
    float  energy = cl->energyDep();
    int    nhits  = Tcl->hits().size();

    Hist->nhits ->Fill(nhits);
    Hist->energy->Fill(energy);
    Hist->time  ->Fill(time);

    Hist->time_vs_nhits  ->Fill(nhits, energy);
    Hist->energy_vs_nhits->Fill(nhits, time);
    return 0;
  }

//-----------------------------------------------------------------------------
// Mode is not used here
//-----------------------------------------------------------------------------
  int CalTimePeakFinderDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;
//-----------------------------------------------------------------------------
// Precalculate some parameters
//-----------------------------------------------------------------------------
    int nseeds       = _data->_tcColl->size();
    int n_good_seeds = 0;

    for (int i=0; i<nseeds; ++i) {
      TimeCluster* tcl = &_data->_tcColl->at(i);
      int nhits        = tcl->hits().size();
      if (nhits >= _data->minNHits) n_good_seeds += 1;
    }
//-----------------------------------------------------------------------------
// fill event histograms
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist._event[0],_data);
    if (n_good_seeds > 0) fillEventHistograms(_hist._event[1],_data);

//-----------------------------------------------------------------------------
// fill cluster histograms
//-----------------------------------------------------------------------------
    int ncl = _data->ccCollection->size();
    for (int i=0; i<ncl; ++i) {
      const CaloCluster* cl = &_data->ccCollection->at(i);
      fillClusterHistograms(_hist._cluster[0],cl);
      if (cl->energyDep() > _data->minClusterEnergy) fillClusterHistograms(_hist._cluster[1],cl);
    }
//-----------------------------------------------------------------------------
// fill time cluster histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<nseeds; ++i) {
      TimeCluster* tcl = &_data->_tcColl->at(i);
      fillTimeClusterHistograms(_hist._tcl[0],tcl);
      int nhits            = tcl->hits().size();
      if (nhits > _data->minNHits) {
        fillTimeClusterHistograms(_hist._tcl[1],tcl);
      }
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(CalTimePeakFinderDiag)

}

