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
#include "BTrk/BaBar/BaBar.hh"
#include "CalPatRec/inc/TrkDefHack.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "TROOT.h"
#include "TFolder.h"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "TrkDiag/inc/KalDiag.hh"

#include "CalPatRec/inc/ObjectDumpUtils.hh"

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
using CLHEP::Hep3Vector;

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
    int              _diag;
    int              _debugLevel;
    int              _printfreq;
    bool             _addhits; 
    TrkParticle      _tpart;	        // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search
					// event object labels
    std::string      _trkPatRecModuleLabel;
    std::string      _calPatRecModuleLabel;

    float            _minTrkQualA;
    float            _minTrkQualB;

    std::string      _iname;	        // data instance name

    KalDiag*         _kalDiag;
  };
  
  MergePatRec::MergePatRec(fhicl::ParameterSet const& pset) :
    _diag                (pset.get<int>("diagLevel",0)),
    _debugLevel          (pset.get<int>("debugLevel",0)),
    _tpart               ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir                ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _trkPatRecModuleLabel(pset.get<std::string>("trkPatRecModuleLabel","TrkPatRec")),
    _calPatRecModuleLabel(pset.get<std::string>("calPatRecModuleLabel","CalPatRec")),
    _minTrkQualA         (pset.get<float>("minTrkQualA")),
    _minTrkQualB         (pset.get<float>("minTrkQualB"))
  {
    // tag the data product instance by the direction and particle type found by this module
    _iname = _fdir.name() + _tpart.name();

    produces<AlgorithmIDCollection>  (_iname);
    produces<KalRepPtrCollection>    (_iname);
    
    _kalDiag = new KalDiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet()));
  }

  MergePatRec::~MergePatRec(){
  }
  
  void MergePatRec::beginJob() {
  }
  
  void MergePatRec::beginRun(art::Run& ){
  }
  
//-----------------------------------------------------------------------------
  void MergePatRec::produce(art::Event& AnEvent) {

					// assume less than 100 tracks
    int const   max_ntrk(100);
    int         tpr_flag[max_ntrk], cpr_flag[max_ntrk], ntpr(0), ncpr(0);

    art::Handle<mu2e::KalRepPtrCollection>    tpr_h, cpr_h;

    mu2e::KalRepPtrCollection  *list_of_kreps_tpr(0), *list_of_kreps_cpr(0);

    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection );
    unique_ptr<KalRepPtrCollection>    trackPtrs(new KalRepPtrCollection   );

    if (_debugLevel > 0) ObjectDumpUtils::printEventHeader(&AnEvent,"MergePatRec::produce");

    AnEvent.getByLabel(_trkPatRecModuleLabel,_iname,tpr_h);
    AnEvent.getByLabel(_calPatRecModuleLabel,_iname,cpr_h);
    
    if (tpr_h.isValid()) { 
      list_of_kreps_tpr = (mu2e::KalRepPtrCollection*) &(*tpr_h);
      ntpr              = list_of_kreps_tpr->size();
    }

    if (cpr_h.isValid()) {
      list_of_kreps_cpr = (mu2e::KalRepPtrCollection*) &(*cpr_h);
      ncpr              = list_of_kreps_cpr->size();
    }

    for (int i=0; i<max_ntrk; i++) {
      tpr_flag[i] = 1;
      cpr_flag[i] = 1;
    }

    art::Ptr<KalRep>          *tpr, *cpr;
    const KalRep              *krep_tpr, *krep_cpr;
    Hep3Vector                cpr_mom, tpr_mom;
    short                     best(-1),  mask;
    AlgorithmID               alg_id;
    TrkHitVector              tlist, clist;
    int                       nat, nac, natc;
    const mu2e::TrkStrawHit  *hitt, *hitc;
    //    double                    tpr_chisq, cpr_chisq;

    TrkInfo                   tpr_trkinfo, cpr_trkinfo;

    for (int i1=0; i1<ntpr; i1++) {
      tpr       = &list_of_kreps_tpr->at(i1);
      tpr_mom   = (*tpr)->momentum();
      //      tpr_chisq = (*tpr)->chisq();
      mask      = 1 << AlgorithmID::TrkPatRecBit;
      tlist     = (*tpr)->hitVector();
      nat       = (*tpr)->nActive();
      natc      = 0;

      krep_tpr = tpr->get();
      _kalDiag->kalDiag(krep_tpr,false);

      tpr_trkinfo = _kalDiag->_trkinfo;

      for (int i2=0; i2<ncpr; i2++) {
	cpr       = &list_of_kreps_cpr->at(i2);
	cpr_mom   = (*cpr)->momentum();
	//	cpr_chisq = (*cpr)->chisq();
	clist     = (*cpr)->hitVector();
	nac       = (*cpr)->nActive();

	krep_cpr  = cpr->get();
	_kalDiag->kalDiag(krep_cpr,false);
	cpr_trkinfo = _kalDiag->_trkinfo;
//-----------------------------------------------------------------------------
// primitive check if this is the same track - require delta(p) less than 5 MeV/c
// ultimately - check the number of common hits
//-----------------------------------------------------------------------------
	for(auto itt=tlist.begin(); itt<tlist.end(); itt++) {
	  hitt = static_cast<const mu2e::TrkStrawHit*> (*itt);
	  if (hitt->isActive()) {
	    for(auto itc=clist.begin(); itc<clist.end(); itc++) {
	      hitc = static_cast<const mu2e::TrkStrawHit*> (*itc);
	      if (hitc->isActive()) {
		if (&hitt->strawHit() == &hitc->strawHit()) {
		  natc += 1;
		  break;
		}
	      }
	    }
	  }
	}
//-----------------------------------------------------------------------------
// if > 50% of all hits are common, consider cpr and tpr to be the same
// logic of the choice: 
// 1. take the track which has more active hits
// 2. if two tracks have the same number of active hits, choose the one with better chi2
//
//-----------------------------------------------------------------------------
	if (natc > (nac+nat)/4.) {

	  mask = mask | (1 << AlgorithmID::CalPatRecBit);

	  if (cpr_trkinfo._trkqual > _minTrkQualA) {
	    trackPtrs->push_back(*cpr);
	    best    = AlgorithmID::CalPatRecBit;
	  }
	  else if (tpr_trkinfo._trkqual > _minTrkQualA) {
	    trackPtrs->push_back(*tpr);
	    best    = AlgorithmID::TrkPatRecBit; 
	  }
	  else if (cpr_trkinfo._trkqual > _minTrkQualB) {
	    trackPtrs->push_back(*cpr);
	    best    = AlgorithmID::CalPatRecBit; 
	  }
	  else {
//-----------------------------------------------------------------------------
// final choice : pick track with larger trkqual
//-----------------------------------------------------------------------------
	    if (tpr_trkinfo._trkqual > cpr_trkinfo._trkqual) {
	      trackPtrs->push_back(*tpr);
	      best    = AlgorithmID::TrkPatRecBit; 
	    }
	    else {
	      trackPtrs->push_back(*cpr);
	      best    = AlgorithmID::CalPatRecBit; 
	    }
	  }

	  tpr_flag[i1] = 0;
	  cpr_flag[i2] = 0;
	  break;
	}
      }

      if (tpr_flag[i1] == 1) {
	trackPtrs->push_back(*tpr);
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
	cpr = &list_of_kreps_cpr->at(i);

	trackPtrs->push_back(*cpr);

	best = AlgorithmID::CalPatRecBit;
	mask = 1 << AlgorithmID::CalPatRecBit;

	alg_id.Set(best,mask);
	algs->push_back(alg_id);
      }
    }

    if (_debugLevel > 0) ObjectDumpUtils::printKalRepCollection(&AnEvent,trackPtrs.get(),1);

    AnEvent.put(std::move(trackPtrs),_iname);
    AnEvent.put(std::move(algs     ),_iname);
  }


//-----------------------------------------------------------------------------
// end job : 
//-----------------------------------------------------------------------------
  void MergePatRec::endJob() {
  }
  
}

using mu2e::MergePatRec;
DEFINE_ART_MODULE(MergePatRec);
