//S Middleton, Nov. 2018

#include "RecoDataProducts/inc/TimeCluster.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"


using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

//-----------------------------------------------------------------------------
// CosmicTrackFinderData
//-------------construct:----------------------------------
  CosmicTrackFinderData::CosmicTrackFinderData() {
    _chHitsToProcess.reserve(kNMaxChHits);
    _chHitsWPos     .reserve(kNMaxChHits);
  }



//--------------destruct:-----------------------------------------------
  CosmicTrackFinderData::~CosmicTrackFinderData() {
     
  }

//-----------------------------------------------------------------------------
  void CosmicTrackFinderData::orderID(CosmicTrackFinderData::ChannelID* X, CosmicTrackFinderData::ChannelID* O) {
/*
Orders the channels. Starts by finding the face of the X channel (i.e. the original start point channel). Sets "0" channel station and plane to X. Then if even station Face of "0" is 1-x and if odd face is x+2.
*/
    if (X->Panel % 2 == 0) X->Face = 0;
    else                   X->Face = 1; // define original face
    
    O->Station = X->Station; // stations already ordered
    O->Plane   = X->Plane;   // planes already ordered, but not necessary for ordered construct

    if (X->Station % 2 == 0) {
      if (X->Plane == 0) O->Face = 1 - X->Face;
      else               O->Face = X->Face + 2;
    }
    else {
      if (X->Plane == 0) O->Face = X->Face;
      else               O->Face = 3 - X->Face; // order face
    }
    
    O->Panel = int(X->Panel/2);                // order panel

    // int n = X->Station + X->Plane + X->Face;   // pattern has no intrinsic meaning, just works
    // if (n % 2 == 0) O->Layer = 1 - X->Layer;
    // else            O->Layer = X->Layer;       // order layer    
  }

//-----------------------------------------------------------------------------
// don't clear the diagnostics part.
//-----------------------------------------------------------------------------
  void CosmicTrackFinderData::clearTempVariables() {

    _timeCluster    = NULL;
    _timeClusterPtr = art::Ptr<TimeCluster>();
    
    _chHitsToProcess.clear();
    _chHitsWPos.clear();

    _nStrawHits = 0;
    _nComboHits = 0;
    
    
    
    _nFiltComboHits = 0;
    _nFiltStrawHits = 0;

    _nXYSh       = 0;


    _nStrawHits  = 0;
    _nComboHits  = 0;

    _nFiltComboHits = 0;
    _nFiltStrawHits = 0;
    //clear the panel-based structure
    for (int f=0; f<StrawId::_ntotalfaces; ++f) {
      FaceZ_t* facez = &_oTracker[f];
      facez->bestFaceHit = -1;
      facez->idChBegin   = -1;
      facez->idChEnd     = -1;
      for (int p=0; p<FaceZ_t::kNPanels; ++p){
	PanelZ_t* panelz  = &facez->panelZs[p];   
	panelz->idChBegin = -1;
	panelz->idChEnd   = -1;
      }
    }
    
    
  }

//-----------------------------------------------------------------------------
// don't clear the diagnostics part.
//-----------------------------------------------------------------------------
   void CosmicTrackFinderData::clearResults() {

    // _goodhits.clear();
    
    //_sxy.clear();
    
   
    _nXYSh       = 0;
    _nStrawHits  = 0;
    _nComboHits  = 0;
   }


//-----------------------------------------------------------------------------
  void CosmicTrackFinderData::print(const char* Title) {

    printf(" CosmicTrackFinderData::print: %s\n",Title);
    
  }

};
