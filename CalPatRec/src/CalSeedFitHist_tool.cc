//

#include "CalPatRec/inc/CalSeedFit_types.hh"
#include "CalPatRec/inc/CprMcUtilsBase.hh"
#include "CalPatRec/inc/CprModuleHistBase.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

// #include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "TH1F.h"

namespace mu2e {


  class CalSeedFitHist : public mu2e::CprModuleHistBase {
  protected:
    int          _mcTruth;
    std::string  _shDigiLabel;
    //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;
    
    std::unique_ptr<CprMcUtilsBase> _mcUtils;

  public:

    CalSeedFitHist(const fhicl::ParameterSet& PSet);
    ~CalSeedFitHist();

  private:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) override ;
    virtual int fillHistograms(int Mode, const TObject* Data, TObject* Hist) override ;
  };


//-----------------------------------------------------------------------------
  CalSeedFitHist::CalSeedFitHist(const fhicl::ParameterSet& PSet) {
    _mcTruth = PSet.get <int >("mcTruth"); 

    if (_mcTruth != 0) _mcUtils = art::make_tool<CprMcUtilsBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else               _mcUtils = std::make_unique<CprMcUtilsBase>();
  }

//-----------------------------------------------------------------------------
  CalSeedFitHist::~CalSeedFitHist() {
  }


//-----------------------------------------------------------------------------
  int CalSeedFitHist::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) {

    CalSeedFit_Hist_t* hist = (CalSeedFit_Hist_t*) Hist;
 
    //    art::TFileDirectory sf_dir = tfs->mkdir("SeedFit");

    hist->nhits          = Tfs->make<TH1F>("nhits" , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
    hist->seeddoca   [0] = Tfs->make<TH1F>("hseeddoca_0","doca seedfit active hits; doca [mm]", 1000, -20., 20);
    hist->seeddoca   [1] = Tfs->make<TH1F>("hseeddoca_1","doca seedfit non active hits; doca [mm]", 1000, -20., 20);
    hist->seeddoca   [2] = Tfs->make<TH1F>("hseeddoca_2","doca seedfit all hits; doca [mm]", 1000, -20., 20);

    hist->chi2       [0] = Tfs->make<TH1F>("chi20"  , "chi2 distribution: all tacks", 100, 0., 10.);
    hist->chi2       [1] = Tfs->make<TH1F>("chi21"  , "chi2 distribution:  tacks with nhits>15", 100, 0., 10.);
    hist->p          [0] = Tfs->make<TH1F>("p0"  , "p distribution: all tacks", 400, 0., 200.);
    hist->p          [1] = Tfs->make<TH1F>("p1"  , "p distribution:  tacks with nhits > 15", 400, 0., 200.);
    hist->NpointsSeed[0] = Tfs->make<TH1F>("hnpseed_0","# of points de-activated by seedfit]", 50, 0., 50);
    hist->NpointsSeed[1] = Tfs->make<TH1F>("hnpseed_1","# of rescued points", 50, 0., 50);
           
    hist->ntracks    [0] = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
    hist->ntracks    [1] = Tfs->make<TH1F>("nseeds1"  , "number of track candidates with nhits > 15 ", 21, -0.5, 20.5);

    return 0;
  }


//-----------------------------------------------------------------------------
  int CalSeedFitHist::fillHistograms(int Mode, const TObject* Data, TObject* Hist) {

    CalSeedFit_Hist_t* hist = (CalSeedFit_Hist_t*) Hist;
    CalSeedFit_Data_t* data = (CalSeedFit_Data_t*) Data;


    if (Mode == 0) {
//-----------------------------------------------------------------------------
// event-level histograms
//-----------------------------------------------------------------------------
      hist->ntracks[0]->Fill(data->ntracks);
      hist->ntracks[1]->Fill(data->ntracks);
    }
    else if (Mode == 1) {
//-----------------------------------------------------------------------------
// track-level histograms
//-----------------------------------------------------------------------------
      int   nactive(0);
      int   ndeactivated(0);

      KalRep* krep = data->result->_krep;

      TrkHitVector const&  hot_l = krep->hitVector();

      for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
	const TrkStrawHit* hit = static_cast<const mu2e::TrkStrawHit*> (*it);
	if (hit->isActive()) ++nactive;
	else                 ++ndeactivated;
      }

      hist->NpointsSeed[0]->Fill(ndeactivated);
      hist->NpointsSeed[1]->Fill(data->nrescued[0]);

      double  chi2           = krep->chisq()/nactive;
      double  h1_fltlen      = krep->firstHit()->kalHit()->hit()->fltLen();
      double  hn_fltlen      = krep->lastHit ()->kalHit()->hit()->fltLen();
      double  entlen         = std::min(h1_fltlen, hn_fltlen);

      CLHEP::Hep3Vector fitmom = krep->momentum(entlen);

      hist->nhits  ->Fill(nactive);
      hist->chi2[0]->Fill(chi2);     
      hist->p   [0]->Fill(fitmom.mag());     

      if (nactive > 15){
	hist->chi2[1]->Fill(chi2);     
	hist->p   [1]->Fill(fitmom.mag());     
      }
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(CalSeedFitHist)

}
