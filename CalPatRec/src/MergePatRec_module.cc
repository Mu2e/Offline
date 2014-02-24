///////////////////////////////////////////////////////////////////////////////
// $Id: MergePatRec_module.cc,v 1.2 2014/02/24 18:25:38 murat Exp $
// $Author: murat $ 
// $Date: 2014/02/24 18:25:38 $
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
    unsigned _iev;
					// configuration parameters
    int _diag,_debug;
    int _printfreq;
    bool _addhits; 
					// event object labels
    std::string _trkPatRecLabel;
    std::string _calPatRecLabel;
  };
  
  MergePatRec::MergePatRec(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _trkPatRecLabel(pset.get<std::string>("trkPatReclabel","TrkPatRec")),
    _calPatRecLabel(pset.get<std::string>("calPatReclabel","CalPatRec"))
  {
    // tag the data product instance by the direction and particle type found by this module
    produces<KalRepCollection>("DownstreameMinus");
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

    unique_ptr<KalRepCollection> tracks(new KalRepCollection );

    AnEvent.getByLabel(_trkPatRecLabel,"DownstreameMinus",tpr_h);
    AnEvent.getByLabel(_calPatRecLabel,"DownstreameMinus",cpr_h);

    list_of_kreps_tpr = (mu2e::KalRepCollection*) &(*tpr_h);
    list_of_kreps_cpr = (mu2e::KalRepCollection*) &(*cpr_h);
    
    // assume less 10 tracks

    int   tpr_flag[10], cpr_flag[10], ntpr, ncpr;

    ntpr = list_of_kreps_tpr->size();
    ncpr = list_of_kreps_cpr->size();

    for (int i=0; i<10; i++) {
      tpr_flag[i] = 1;
      cpr_flag[i] = 1;
    }

    KalRep  *tpr, *cpr;

    Hep3Vector cpr_mom, tpr_mom;

    for (int i1=0; i1<ntpr; i1++) {
      tpr  = (KalRep*) list_of_kreps_tpr->at(i1);
      tpr_mom = tpr->momentum();

      for (int i2=0; i2<ncpr; i2++) {
	cpr = (KalRep*) list_of_kreps_cpr->at(i2);
	cpr_mom = cpr->momentum();
//-----------------------------------------------------------------------------
// primitive check if this is the same track - require delta(p) less than 5 MeV/c
// ultimately - check the number of common hits
//-----------------------------------------------------------------------------
	if (fabs (cpr_mom.mag()-tpr_mom.mag()) < 5.) {

					// decide that track is the same, 
					// use the version which has more hits

	  if (tpr->nDof() > cpr->nDof()) {
	    tracks->push_back(tpr->clone());
	    cpr_flag[i2] = 0;
	  }
	  else {	
	    tracks->push_back(cpr->clone());
	    tpr_flag[i2] = 0;
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// account for presence of multiple tracks
//-----------------------------------------------------------------------------
    for (int i=0; i<ncpr; i++) {
      if (cpr_flag[i] == 1) {
	cpr = (KalRep*) list_of_kreps_cpr->at(i);
	tracks->push_back(cpr->clone());
      }
    }


    AnEvent.put(std::move(tracks),"DownstreameMinus");
  }

  void MergePatRec::endJob() {
  }
  
}

using mu2e::MergePatRec;
DEFINE_ART_MODULE(MergePatRec);
