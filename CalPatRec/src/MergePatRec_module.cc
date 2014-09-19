///////////////////////////////////////////////////////////////////////////////
// $Id: MergePatRec_module.cc,v 1.8 2014/09/19 20:49:45 murat Exp $
// $Author: murat $ 
// $Date: 2014/09/19 20:49:45 $
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
    _tpart         ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir          ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _trkPatRecLabel(pset.get<std::string>("trkPatReclabel","TrkPatRec")),
    _calPatRecLabel(pset.get<std::string>("calPatReclabel","CalPatRec"))
  {
    // tag the data product instance by the direction and particle type found by this module
    _iname = _fdir.name() + _tpart.name();

    produces<AlgorithmIDCollection>(_iname);
    produces<KalRepPtrCollection>  (_iname);
  }

  MergePatRec::~MergePatRec(){
  }
  
  void MergePatRec::beginJob() {
  }
  
  void MergePatRec::beginRun(art::Run& ){}
  
//-----------------------------------------------------------------------------
  void MergePatRec::produce(art::Event& AnEvent) {
    
    art::Handle<mu2e::KalRepPtrCollection> tpr_h, cpr_h;

    mu2e::KalRepPtrCollection  *list_of_kreps_tpr(0), *list_of_kreps_cpr(0);

    unique_ptr<AlgorithmIDCollection> algs     (new AlgorithmIDCollection);
    unique_ptr<KalRepPtrCollection>   trackPtrs(new KalRepPtrCollection  );

    AnEvent.getByLabel(_trkPatRecLabel,_iname,tpr_h);
    AnEvent.getByLabel(_calPatRecLabel,_iname,cpr_h);

    list_of_kreps_tpr = (mu2e::KalRepPtrCollection*) &(*tpr_h);
    list_of_kreps_cpr = (mu2e::KalRepPtrCollection*) &(*cpr_h);
    
    // assume less than 100 tracks

    int const max_ntrk(100);
    int   tpr_flag[max_ntrk], cpr_flag[max_ntrk], ntpr, ncpr;

    ntpr = list_of_kreps_tpr->size();
    ncpr = list_of_kreps_cpr->size();

    for (int i=0; i<max_ntrk; i++) {
      tpr_flag[i] = 1;
      cpr_flag[i] = 1;
    }

    art::Ptr<KalRep>  tpr, cpr;
    Hep3Vector        cpr_mom, tpr_mom;
    short             best, mask;
    AlgorithmID       alg_id;

    for (int i1=0; i1<ntpr; i1++) {
      //      tpr     = (KalRep*) list_of_kreps_tpr->get(i1);
      tpr     = list_of_kreps_tpr->at(i1);
      tpr_mom = tpr->momentum();
      mask    = 1 << AlgorithmID::TrkPatRecBit;

      for (int i2=0; i2<ncpr; i2++) {
	cpr = list_of_kreps_cpr->at(i2);
	cpr_mom = cpr->momentum();
//-----------------------------------------------------------------------------
// primitive check if this is the same track - require delta(p) less than 5 MeV/c
// ultimately - check the number of common hits
//-----------------------------------------------------------------------------
	if (fabs (cpr_mom.mag()-tpr_mom.mag()) < 5.) {

	  mask = mask | (1 << AlgorithmID::CalPatRecBit);

	  if (tpr->nDof() >= cpr->nDof()) {
	    trackPtrs->push_back(tpr);
	    best    = AlgorithmID::TrkPatRecBit;
	  }
	  else {	
	    trackPtrs->push_back(cpr);
	    best    = AlgorithmID::CalPatRecBit;
	  }

	  tpr_flag[i1] = 0;
	  cpr_flag[i2] = 0;
	  break;
	}
      }

      if (tpr_flag[i1] == 1) {
	trackPtrs->push_back(tpr);
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
	cpr = list_of_kreps_cpr->at(i);
	trackPtrs->push_back(cpr);

	best = AlgorithmID::CalPatRecBit;
	mask = 1 << AlgorithmID::CalPatRecBit;

	alg_id.Set(best,mask);
	algs->push_back(alg_id);
      }
    }

    AnEvent.put(std::move(trackPtrs),_iname);
    AnEvent.put(std::move(algs     ),_iname);
  }

  void MergePatRec::endJob() {
  }
  
}

using mu2e::MergePatRec;
DEFINE_ART_MODULE(MergePatRec);
