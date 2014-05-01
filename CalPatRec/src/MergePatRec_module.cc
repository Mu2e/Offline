///////////////////////////////////////////////////////////////////////////////
// $Id: MergePatRec_module.cc,v 1.6 2014/05/01 19:31:09 murat Exp $
// $Author: murat $ 
// $Date: 2014/05/01 19:31:09 $
// takes inputs from two track finding algorithms, produces one track collection 
// on output to be used for analysis
///////////////////////////////////////////////////////////////////////////////
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"

#include "TROOT.h"
#include "TFolder.h"

#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#
#include "CalPatRec/inc/AlgorithmIDCollection.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
using namespace std; 

namespace mu2e {
  class MergePatRec : public art::EDProducer {
  public:
    explicit MergePatRec(fhicl::ParameterSet const&);
    virtual ~MergePatRec();
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    virtual void produce(art::Event& event ); 
    void endJob();
  private:
    unsigned         _iev;
					// configuration parameters
    int              _diag, _debug;
    int              _printfreq;
    bool             _addhits; 
    TrkParticle      _tpart;	        // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search
					// event object labels
    std::string      _trkPatRecLabel;
    std::string      _calPatRecLabel;
    std::string      _iname;	// data instance name
  };
  
  MergePatRec::MergePatRec(fhicl::ParameterSet const& pset) :
    _diag          (pset.get<int>("diagLevel",0)),
    _debug         (pset.get<int>("debugLevel",0)),
    _tpart         ((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir          ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _trkPatRecLabel(pset.get<std::string>("trkPatReclabel","TrkPatRec")),
    _calPatRecLabel(pset.get<std::string>("calPatReclabel","CalPatRec"))
  {
    // tag the data product instance by the direction and particle type found by this module
    _iname = _fdir.name() + _tpart.name();
    produces<AlgorithmIDCollection>("DownstreameMinus");
    produces<KalRepCollection>     ("DownstreameMinus");
  }

  MergePatRec::~MergePatRec(){
  }
  
  void MergePatRec::beginJob() {
  }
  
  void MergePatRec::beginRun(art::Run& ){}
  
//-----------------------------------------------------------------------------
  void MergePatRec::produce(art::Event& AnEvent) {
    
    art::Handle<mu2e::KalRepCollection> tpr_h, cpr_h;
    mu2e::KalRepCollection  *list_of_kreps_tpr(0), *list_of_kreps_cpr(0);

    unique_ptr<KalRepCollection>      tracks   (new KalRepCollection     );
    unique_ptr<AlgorithmIDCollection> algs     (new AlgorithmIDCollection);
    unique_ptr<KalRepPtrCollection>   trackPtrs(new KalRepPtrCollection  );

    AnEvent.getByLabel(_trkPatRecLabel,"DownstreameMinus",tpr_h);
    AnEvent.getByLabel(_calPatRecLabel,"DownstreameMinus",cpr_h);

    list_of_kreps_tpr = (mu2e::KalRepCollection*) &(*tpr_h);
    list_of_kreps_cpr = (mu2e::KalRepCollection*) &(*cpr_h);
    
//     art::ProductID tpr_id(getProductID<KalRepCollection>(AnEvent,_trkPatRecLabel,_iname));
//     art::ProductID cpr_id(getProductID<KalRepCollection>(AnEvent,_calPatRecLabel,_iname));

    // assume less 10 tracks

    int   tpr_flag[10], cpr_flag[10], ntpr, ncpr;

    ntpr = list_of_kreps_tpr->size();
    ncpr = list_of_kreps_cpr->size();

    for (int i=0; i<10; i++) {
      tpr_flag[i] = 1;
      cpr_flag[i] = 1;
    }

    KalRep        *tpr, *cpr, *new_trk;
    Hep3Vector    cpr_mom, tpr_mom;
    short          best, mask;
    AlgorithmID   alg_id;

    for (int i1=0; i1<ntpr; i1++) {
      tpr     = (KalRep*) list_of_kreps_tpr->get(i1);
      tpr_mom = tpr->momentum();
      mask    = 1 << AlgorithmID::TrkPatRecBit;

      for (int i2=0; i2<ncpr; i2++) {
	cpr = (KalRep*) list_of_kreps_cpr->get(i2);
	cpr_mom = cpr->momentum();
//-----------------------------------------------------------------------------
// primitive check if this is the same track - require delta(p) less than 5 MeV/c
// ultimately - check the number of common hits
//-----------------------------------------------------------------------------
	if (fabs (cpr_mom.mag()-tpr_mom.mag()) < 5.) {

	  mask = mask | (1 << AlgorithmID::CalPatRecBit);

	  if (tpr->nDof() >= cpr->nDof()) {
	    new_trk = tpr->clone();
	    best    = AlgorithmID::TrkPatRecBit;
	  }
	  else {	
	    new_trk = cpr->clone();
	    best    = AlgorithmID::CalPatRecBit;
	  }

	  tracks->push_back(new_trk);

	  tpr_flag[i1] = 0;
	  cpr_flag[i2] = 0;
	  break;
	}
      }

      if (tpr_flag[i1] == 1) {
	new_trk = tpr->clone();
	tracks->push_back(new_trk);
	best = AlgorithmID::TrkPatRecBit;
      }

      alg_id.Set(best,mask);
      algs->push_back(alg_id);
    }
//-----------------------------------------------------------------------------
// account for presence of multiple tracks
//-----------------------------------------------------------------------------
    for (int i=0; i<ncpr; i++) {
      if (cpr_flag[i] == 1) {
	cpr = (KalRep*) list_of_kreps_cpr->get(i);
	new_trk = cpr->clone();
	tracks->push_back(new_trk);

	best = AlgorithmID::CalPatRecBit;
	mask = 1 << AlgorithmID::CalPatRecBit;

	alg_id.Set(best,mask);
	algs->push_back(alg_id);
      }
    }


    AnEvent.put(std::move(tracks),_iname);
    AnEvent.put(std::move(algs  ),_iname);
  }

  void MergePatRec::endJob() {
  }
  
}

using mu2e::MergePatRec;
DEFINE_ART_MODULE(MergePatRec);
