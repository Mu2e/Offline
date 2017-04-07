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

#include "CalPatRec/inc/CalTrkFit_types.hh"
#include "CalPatRec/inc/CprMcUtilsBase.hh"
#include "CalPatRec/inc/CprModuleHistBase.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"

#include "TTrackerGeom/inc/TTracker.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

  class CalTrkFitHist : public mu2e::CprModuleHistBase {
  protected:
    int          _mcTruth;
    std::string  _shDigiLabel;
    //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;
  
    std::unique_ptr<CprMcUtilsBase> _mcUtils;

  public:

    CalTrkFitHist(const fhicl::ParameterSet& PSet);
    ~CalTrkFitHist();

  private:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) override ;
    virtual int fillHistograms(int Mode, const TObject* Data, TObject* Hist) override ;
  };


  CalTrkFitHist::CalTrkFitHist(const fhicl::ParameterSet& PSet) {
    _mcTruth = PSet.get <int >("mcTruth"); 

    if (_mcTruth != 0) _mcUtils = art::make_tool<CprMcUtilsBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else               _mcUtils = std::make_unique<CprMcUtilsBase>();
  }

//-----------------------------------------------------------------------------
  CalTrkFitHist::~CalTrkFitHist() {
  }


//-----------------------------------------------------------------------------
  int CalTrkFitHist::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) {

    CalTrkFit_Hist_t* hist = (CalTrkFit_Hist_t*) Hist;
  
    hist->nhits       = Tfs->make<TH1F>("nhits" , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
    hist->chi2[0]     = Tfs->make<TH1F>("chi20"  , "chi2 distribution: all tacks", 100, 0., 10.);
    hist->chi2[1]     = Tfs->make<TH1F>("chi21"  , "chi2 distribution:  tacks with nhits>15", 100, 0., 10.);
    hist->p   [0]     = Tfs->make<TH1F>("p0"  , "p distribution: all tacks", 400, 0., 200.);
    hist->p   [1]     = Tfs->make<TH1F>("p1"  , "p distribution:  tacks with nhits > 15", 400, 0., 200.);
  
    hist->kaldoca[0]  = Tfs->make<TH1F>("hkaldoca_0","doca kalfit active hits; doca [mm]", 1000, -20., 20);
    hist->kaldoca[1]  = Tfs->make<TH1F>("hkaldoca_1","doca kalfit non active hits; doca [mm]", 1000, -20., 20);
  
    hist->ntracks[0]  = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
    hist->ntracks[1]  = Tfs->make<TH1F>("nseeds1"  , "number of track candidates with nhits > 15 ", 21, -0.5, 20.5);

    for (int i=0; i<2; i++) {
      hist->dSlopeOS[i] = Tfs->make<TH1F>(Form("ds_os_%i"  ,i), Form("Delta Slope OS[%i]",i),200, -1,   1);
      hist->dSlopeSS[i] = Tfs->make<TH1F>(Form("ds_ss_%i"  ,i), Form("Delta Slope SS[%i]",i),200, -1,   1);
      hist->chi2bOS[i]  = Tfs->make<TH1F>(Form("chi2bos_%i",i), Form("chi2(best)  OS[%i]",i),200,  0, 200);
      hist->chi2bSS[i]  = Tfs->make<TH1F>(Form("chi2bss_%i",i), Form("chi2(best)  SS[%i]",i),200,  0, 200);

      hist->chi2rOS[i]  = Tfs->make<TH1F>(Form("chi2ros_%i",i), Form("chi2R (best/next)  OS[%i]",i),200,  0, 1);
      hist->chi2rSS[i]  = Tfs->make<TH1F>(Form("chi2rss_%i",i), Form("chi2R (best/next)  SS[%i]",i),200,  0, 1);
    }

    return 0;
  }


//-----------------------------------------------------------------------------
  int CalTrkFitHist::fillHistograms(int Mode, const TObject* Data, TObject* Hist) {

    //    static int                                   first_call  (1);
    //    static int                                   event_number(-1);

    CalTrkFit_Hist_t* hist = (CalTrkFit_Hist_t*) Hist;
    CalTrkFit_Data_t* data = (CalTrkFit_Data_t*) Data;
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
      hist->ntracks[0]->Fill(data->ntracks[0]);
      hist->ntracks[1]->Fill(data->ntracks[1]);
    }
    else if (Mode == 1) {
//-----------------------------------------------------------------------------
// fill reco doca histograms
//-----------------------------------------------------------------------------
      const mu2e::TrkStrawHit* hit;
      int                      hit_index;
      CLHEP::Hep3Vector        tdir;
      HepPoint                 tpos;

      KalRep* krep = data->result->_krep;

      TrkHitVector const& hot_l = krep->hitVector();

      krep->traj().getInfo(0.0,tpos,tdir);

					// loop over the track hits
      int nhits = krep->hitVector().size();

      for (int i=0; i<nhits; ++i) {
	hit = static_cast<TrkStrawHit*> (krep->hitVector().at(i));
	int               hIndex = hit->index();
	StrawHit const*   sh     = & data->result->_shcol->at(hIndex);
	Straw const&      straw  = data->tracker->getStraw(sh->strawIndex());
	CLHEP::Hep3Vector hpos   = straw.getMidPoint();
	CLHEP::Hep3Vector hdir   = straw.getDirection();
	bool              found  = false;

	// convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	HepPoint          spt(hpos.x(),hpos.y(),hpos.z());
	TrkLineTraj       htraj(spt,hdir,-20,20);
	// estimate flightlength along track.  This assumes a constant BField!!!
	double           fltlen = (hpos.z()-tpos.z())/tdir.z();
	TrkPoca          hitpoca(krep->traj(),fltlen,htraj,0.0);

	double           doca   = hitpoca.doca();
	for(auto it=hot_l.begin(); it<hot_l.end(); it++) {
	  hit = static_cast<const mu2e::TrkStrawHit*> (*it);
	  if (!hit->isActive()) continue;
	  hit_index = hit->index();
	  if (hIndex == hit_index) {
	    found = true;
	    break;
	  }
	}
	  
	if (found) hist->kaldoca[0]->Fill(doca);
	else       hist->kaldoca[1]->Fill(doca);
      }
//-----------------------------------------------------------------------------
// doublet histograms
//-----------------------------------------------------------------------------
      Doublet* d;
      int nd = data->listOfDoublets->size();
      for (int i=0; i<nd; i++) {
	d = &data->listOfDoublets->at(i);
	if (d->fNStrawHits == 2) {

	  int same_sign   = d->isSameSign();
	  float dsl       = d->bestDxDzRes();
	  float chi2_best = d->Chi2Best();
	  float chi2r     = chi2_best/d->fChi2[d->fINext];

	  int a0 = d->fHit[0]->ambig();
	  int a1 = d->fHit[0]->ambig();

	  bool h1_ok = ((a0 != 0) && (a0*d->fMcDoca[0] > 0));
	  bool h2_ok = ((a1 != 0) && (a1*d->fMcDoca[1] > 0));

	  if   (same_sign) {
	    hist->dSlopeSS[0]->Fill(dsl);
	    hist->chi2bSS[0]->Fill(chi2_best);
	    hist->chi2rSS[0]->Fill(chi2r);

	    if (h1_ok && h2_ok) {
	      hist->dSlopeSS[1]->Fill(dsl);
	      hist->chi2bSS[1]->Fill(chi2_best);
	      hist->chi2rSS[1]->Fill(chi2r);
	    }
	  }
	  else {
	    hist->dSlopeOS[0]->Fill(dsl);
	    hist->chi2bOS[0]->Fill(chi2_best);
	    hist->chi2rOS[0]->Fill(chi2r);
	    if (h1_ok && h2_ok) {
	      hist->dSlopeOS[1]->Fill(dsl);
	      hist->chi2bOS[1]->Fill(chi2_best);
	      hist->chi2rOS[1]->Fill(chi2r);
	    }
	  }
	}
      }
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(CalTrkFitHist)

}

