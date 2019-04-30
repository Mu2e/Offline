//
// Diagnostics for track-calo matching
// $Author: brownd $ 
// $Date: 2014/09/20 14:34:22 $
//
#ifndef TrkCaloDiag_HH
#define TrkCaloDiag_HH
// framework
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
// Track Calo matching
#include "RecoDataProducts/inc/TrackClusterMatch.hh"
#include "TrkDiag/inc/TrkCaloInfo.hh"
// particleId
#include "ParticleID/inc/PIDLogLRatio.hh"
#include "ParticleID/inc/PIDLogL1D.hh"
#include "ParticleID/inc/PIDLogLEp.hh"
// ROOT incldues
#include "TTree.h"
#include "Rtypes.h"
// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

namespace mu2e {

  class TrkCaloDiag {
    public:
    TrkCaloDiag(TrkParticle const& tpart, TrkFitDirection const& fdir, fhicl::ParameterSet const& pset);
    // find the required data products in the event
    void findData(const art::Event& event);
    // add calo info for a particular track
    void addCaloInfo(KalRep const* krep);
    // add calo info for a particular match
    void fillCaloInfo(TrackClusterMatch const& tcm, TrkCaloInfo& tcinfo);
    // add the calo info branches to the tree
    void addBranches(TTree* tree,const char* suffix="");
    art::Handle<TrackClusterMatchCollection> const& caloMatchHandle() { return _caloMatchHandle; }
    private:
// calorimeter matching labels 
    std::string _caloMatchModule;
    art::Handle<TrackClusterMatchCollection> _caloMatchHandle;
    // branch variables
    std::vector<TrkCaloInfo> _caloinfo;
    Int_t _ncalo;

// PID configuration
    typedef PIDLogLRatio<PIDLogL1D> PIDdt;
    typedef PIDLogLRatio<PIDLogLEp> PIDEp;
    // there is no default constructor for PID classes so the configuration MUST be
    // accurate to even instantiantiate these objects
    PIDdt _pid_dt;
    PIDEp _pid_ep;

  };

}
#endif

