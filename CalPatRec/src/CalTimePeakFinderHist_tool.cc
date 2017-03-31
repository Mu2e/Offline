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

#include "CalPatRec/inc/CalTimePeakFinder_types.hh"
#include "CalPatRec/inc/CprMcUtilsBase.hh"
#include "CalPatRec/inc/CprModuleHistBase.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"

#include "TTrackerGeom/inc/TTracker.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "TH1.h"
#include "TH2.h"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

  class CalTimePeakFinderHist : public mu2e::CprModuleHistBase {
  protected:
    int          _mcTruth;
    std::string  _shDigiLabel;
    //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;
  
    std::unique_ptr<CprMcUtilsBase> _mcUtils;

  public:

    CalTimePeakFinderHist(const fhicl::ParameterSet& PSet);
    ~CalTimePeakFinderHist();

  private:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) override ;
    virtual int fillHistograms(int Mode, const TObject* Data, TObject* Hist) override ;
  };


  CalTimePeakFinderHist::CalTimePeakFinderHist(const fhicl::ParameterSet& PSet) {
    _mcTruth = PSet.get <int >("mcTruth"); 

    if (_mcTruth != 0) _mcUtils = art::make_tool<CprMcUtilsBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else               _mcUtils = std::make_unique<CprMcUtilsBase>();
  }

//-----------------------------------------------------------------------------
  CalTimePeakFinderHist::~CalTimePeakFinderHist() {
  }


//-----------------------------------------------------------------------------
  int CalTimePeakFinderHist::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) {

    CalTimePeakFinder_Hist_t* hist = (CalTimePeakFinder_Hist_t*) Hist;

    hist->nseeds[0] = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
    hist->nseeds[1] = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15", 21, -0.5, 20.5);
    hist->nhits     = Tfs->make<TH1F>("nhits"  , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
    hist->energy[0] = Tfs->make<TH1F>("energy0", "cluster energy; E [MeV]"                   , 400, 0., 200.);
    hist->energy[1] = Tfs->make<TH1F>("energy1", "cluster energy, nhits > 15; E [MeV]"       , 400, 0., 200.);
    hist->time  [0] = Tfs->make<TH1F>("time0"  , "cluster time; t [ns]"                      , 2800, 300., 1700);
    hist->time  [1] = Tfs->make<TH1F>("time1"  , "cluster time, nhits > 15; t [ns]"          , 2800, 300., 1700);

    hist->time_vs_nhits   = Tfs->make<TH2F>("time_vs_nhits"  ,"time vs nhits  ; N [#]; t [ns]" , 100, 0, 100, 2800, 300., 1700);
    hist->energy_vs_nhits = Tfs->make<TH2F>("energy_vs_nhits","energy vs nhits; N [#]; E [MeV]", 100, 0, 100, 400, 0, 200);

    return 0;
  }


//-----------------------------------------------------------------------------
  int CalTimePeakFinderHist::fillHistograms(int Mode, const TObject* Data, TObject* Hist) {

    CalTimePeakFinder_Hist_t* hist = (CalTimePeakFinder_Hist_t*) Hist;
    CalTimePeakFinder_Data_t* data = (CalTimePeakFinder_Data_t*) Data;
//-----------------------------------------------------------------------------
// count number of MC hits produced by the MC particle of interest
// uncomment, when find a use case
//-----------------------------------------------------------------------------
    // int n_gen_hits = _NGenHits(data->event              ,
    // 			       data->timeOffsets        ,
    // 			       data->result->_shcol     ,
    // 			       data->result->shDigiLabel);
//-----------------------------------------------------------------------------
// fill histograms
//-----------------------------------------------------------------------------
    if (Mode == 0) {
      hist->nseeds[0]->Fill(data->nseeds[0]);
      hist->nseeds[1]->Fill(data->nseeds[1]);
    }
    else if (Mode == 1) {
//-----------------------------------------------------------------------------
// fill 'per time peak' histograms
//-----------------------------------------------------------------------------
      const CaloCluster* cl = data->cl;

      float  time   = cl->time();
      float  energy = cl->energyDep();
      int    nhits  = data->timeCluster->hits().size();

      hist->energy[0]->Fill(cl->energyDep());
      hist->time  [0]->Fill(cl->time());
      hist->nhits    ->Fill(nhits);

      if (nhits > data->minNHits) {
	hist->energy[1]->Fill(cl->energyDep());
	hist->time  [1]->Fill(cl->time());
      }

      hist->time_vs_nhits  ->Fill(nhits, energy);
      hist->energy_vs_nhits->Fill(nhits, time);
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(CalTimePeakFinderHist)

}

