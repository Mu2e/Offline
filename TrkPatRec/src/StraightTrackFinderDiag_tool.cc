#ifndef __TrkPatRec_StraightTrackFinderDiag_hh__
#define __TrkPatRec_StraightTrackFinderDiag_hh__

#include "TrkPatRec/inc/StraightTrackFinder_types.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "TrkReco/inc/StraightTrackFit.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1F.h"

namespace mu2e {

using namespace StraightTrackFinderTypes;
  
  class StraightTrackFinderDiag : public mu2e::ModuleHistToolBase {
  public:
      struct Hist_t {
          TH1F*  nhits;           // number of hits in track  
          TH1F*  nseeds;   	// number of seeds
	  TH1F*  niter;  		// number of fit iterations
          TH1F*  chi2XY;           // Chi^2 in XY of fit
         
      };
  protected:
      int     _mcTruth;
      Hist_t  _hist;            // owned
      Data_t* _data;            // cached
   public:
      StraightTrackFinderDiag(const fhicl::ParameterSet& PSet);
      ~StraightTrackFinderDiag();
   private:
      virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
      virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };

//-----------------------------------------------------------------------------
   StraightTrackFinderDiag::StraightTrackFinderDiag(const fhicl::ParameterSet& PSet) {
    printf(" StraightTrackFinderDiag::StraightTrackFinderDiag : Running! \n");
  }

//-----------------------------------------------------------------------------
  StraightTrackFinderDiag::~StraightTrackFinderDiag() {
  }

    
//-----------------------------------------------------------------------------
  int StraightTrackFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
      _hist.nseeds     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 20.5);
    //_hist.nseeds[1]     = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15;"    , 21, -0.5, 20.5);
    _hist.nhits       = Tfs->make<TH1F>("nhits"    , "number of hits on a track candidate"       , 401, -0.5, 800.5);
    _hist.niter     = Tfs->make<TH1F>("niter"     , "number of iterations in fit"   , 401, -0.5, 800.5);
    _hist.chi2XY   = Tfs->make<TH1F>("chisXY", "#chi^{2}/ndof from the xy fit "   , 2001, -0.05, 200.05);
      return 0;
  }
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Mode is not used
//-----------------------------------------------------------------------------
  int StraightTrackFinderDiag::fillHistograms(void* Data, int Mode) {
    std::cout<<"In Fill Hist in Diag Tool "<<std::endl;
    _data = (Data_t*) Data;
    _hist.nseeds->Fill(_data->nseeds[0]);
    //_hist.nseeds[1]->Fill(_data->nseeds[1]);
    /*
    for (int i=0; i<_data->nseeds[0]; i++) {
	_hist.nhits->Fill(_data->nhits[i]);
        _hist.chi2XY->Fill(_data->chi2XY[i]);
        _hist.niter->Fill(_data->niter[i]);
        
    }
    */
    return 0;
  }

  //}//End Class
}//End Mu2E NameSpace

DEFINE_ART_CLASS_TOOL(mu2e::StraightTrackFinderDiag)

#endif
