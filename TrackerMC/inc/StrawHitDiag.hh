#ifndef StrawHitDiag_HH
#define StrawHitDiag_HH

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "TrackerMC/inc/SHInfo.hh"
#include "TrkChargeReco/inc/PeakFitParams.hh"
#include "TrackerGeom/inc/Straw.hh"

#include "TFile.h"
#include "TTree.h"


namespace mu2e {


  class StrawHitDiag  {
     

     public:
     
        StrawHitDiag() {};
        ~StrawHitDiag() {};
	
        void init();
        void fill(const Straw& straw, const StrawHit &newhit,
                  const TrkChargeReco::PeakFitParams& params, 
                  const StrawDigiMCCollection* mcdigis, int isd);
	


    private:
       
        void fillDiagMC(const Straw& straw, const StrawDigiMC& mcdigi);
        
        TTree* _shdiag;
        SHID   _shid; 
        TrkChargeReco::PeakFitParams _peakfit; 
        Float_t _edep, _time, _dt;
        SHMCInfo _shmcinfo; 

  };

}
#endif
