///////////////////////////////////////////////////////////////////////////////
// $Id: $
// $Author: $ 
// $Date: $
// takes inputs from two helix finding algorithms, produces one helix collection 
// on output to be used for the track seed-fit
//
// Original author P. Murat
//
//
///////////////////////////////////////////////////////////////////////////////
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/HelixParams.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "TROOT.h"
#include "TFolder.h"
#include "TVector2.h"

#include "RecoDataProducts/inc/HelixSeed.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TrkDiag/inc/KalDiag.hh"
#include "RecoDataProducts/inc/Doublet.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"

// CalPatRec
// #include "CalPatRec/inc/TrkDefHack.hh"
#include "CalPatRec/inc/LsqSums4.hh"
#include "CalPatRec/inc/ObjectDumpUtils.hh"

#include "CalPatRec/inc/AlgorithmIDCollection.hh"
// Xerces XML Parser
#include <xercesc/dom/DOM.hpp>

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
  class MergeHelixFinder : public art::EDProducer {
  public:
    explicit MergeHelixFinder(fhicl::ParameterSet const&);
    virtual ~MergeHelixFinder();
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    virtual void produce(art::Event& event ); 
    void endJob();

    double calculateWeight(Hep3Vector& HitPos   ,
			   Hep3Vector& StrawDir ,
			   Hep3Vector& HelCenter, 
			   double      Radius   ,
			   int         WeightMode);
    double helix_momentum(const HelixSeed* Helix);
    double helix_chi2xy  (const HelixSeed* Helix);
    double helix_chi2phiz(const HelixSeed* Helix);

    
  private:
    unsigned         _iev;
					// configuration parameters
    int              _diag;
    int              _debugLevel;
    float            _minTprChi2;
    float            _minCprChi2;
    int              _printfreq;
    bool             _addhits; 
					// event object labels
    std::string      _trkHelixFinderModuleLabel;
    std::string      _calHelixFinderModuleLabel;
  };
  
  MergeHelixFinder::MergeHelixFinder(fhicl::ParameterSet const& pset) :
    _diag                        (pset.get<int>("diagLevel" )),
    _debugLevel                  (pset.get<int>("debugLevel")),
    _minTprChi2                  (pset.get<int>("minTprChi2")),   
    _minCprChi2                  (pset.get<int>("minCprChi2")),   
    _trkHelixFinderModuleLabel   (pset.get<std::string>("trkHelixFinderModuleLabel"   )),
    _calHelixFinderModuleLabel   (pset.get<std::string>("calHelixFinderModuleLabel"   ))
  {

    produces<AlgorithmIDCollection>  ();
    produces<HelixSeedCollection>    ();
    
  }

  MergeHelixFinder::~MergeHelixFinder() {
  }
  
  void MergeHelixFinder::beginJob() {
  }
  
  void MergeHelixFinder::beginRun(art::Run& ) {
  }


//--------------------------------------------------------------------------------
// evaluate the momentum of a Helix
//--------------------------------------------------------------------------------
  double MergeHelixFinder::helix_momentum(const HelixSeed* Helix){
    const RobustHelix*        robustHel = &Helix->helix();

    double      radius   = robustHel->radius();
    double      lambda   = robustHel->lambda();  
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;//FIX ME! that work assuming 1 Tesla axial magntic field
    double      pT       = radius*mm2MeV;
    double      p        = pT/std::cos( std::atan(tanDip));
    
    if (p>0) return p;
    else     return -1;
  }
  

//--------------------------------------------------------------------------------
// define the function used for projecting the strawhit error along the radial 
// direction of the helix-circle
//
// note: imported from HeliFitHack class
//--------------------------------------------------------------------------------

  double  MergeHelixFinder::calculateWeight(Hep3Vector& HitPos   ,
					    Hep3Vector& StrawDir ,
					    Hep3Vector& HelCenter, 
					    double             Radius   ,
					    int                WeightMode) {//WeightMode = 1 is for XY chi2 , WeightMode = 0 is for Phi-z chi2
  
  double    rs(2.5);   // straw radius, mm
  double    ew(30.0);  // assumed resolution along the wire, mm
  
  double x  = HitPos.x();
  double y  = HitPos.y();
  double dx = x-HelCenter.x();
  double dy = y-HelCenter.y();
  
  double costh  = (dx*StrawDir.x()+dy*StrawDir.y())/sqrt(dx*dx+dy*dy);
  double sinth2 = 1-costh*costh;
  
  double wt(0);
                                              //scale the weight for having chi2/ndof distribution peaking at 1
  if ( WeightMode == 1){
    double e2     = ew*ew*sinth2+rs*rs*costh*costh;
    wt  = 1./e2;
    wt *= 0.2941; //FIX ME! it should be get from CalPatRec/fcl/prolog.fcl
  } else if (WeightMode ==0 ){
    double e2     = ew*ew*costh*costh+rs*rs*sinth2;
    wt     = Radius*Radius/e2;
    wt    *=  0.174;//FIX ME! it should be get from CalPatRec/fcl/prolog.fcl
  }
  
  return wt;
}

//--------------------------------------------------------------------------------
// evaluate the momentum of a Helix
//--------------------------------------------------------------------------------
  double MergeHelixFinder::helix_chi2xy  (const HelixSeed* Helix){
    const RobustHelix*        robustHel = &Helix->helix();
    const HelixHitCollection* hits      = &Helix->hits();
    const HelixHit*           hit(0);
    const mu2e::CaloCluster*  cluster   = Helix->caloCluster().get();
    int                       nhits     = hits->size();
    double                    radius    = robustHel->radius();
    CLHEP::Hep3Vector         pos(0), helix_pos(0), wdir(0), sdir(0), helix_center(0);
    double                    chi2(0);

    static const CLHEP::Hep3Vector zdir(0.0,0.0,1.0);

    helix_center = robustHel->center();

    ::LsqSums4 sxy;
    //add the stopping target center as in CalHeliFinderAlg.cc
    sxy.addPoint(0., 0., 1./900.);

    for (int i=0; i<nhits; ++i){
      hit       = &hits->at(i);
      pos       = hit->pos();
      wdir      = hit->wdir();
      sdir      = zdir.cross(wdir);
      helix_pos = pos;
      robustHel->position(helix_pos);
      
      double    weight   = calculateWeight(pos, sdir, helix_center, radius, 1);
      sxy.addPoint(pos.x(), pos.y(), weight);
    }

    if (cluster != 0){
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      const mu2e::Calorimeter*  _calorimeter = ch.get();      
      CLHEP::Hep3Vector         gpos = _calorimeter->geomUtil().diskToMu2e(cluster->diskId(),cluster->cog3Vector());
      CLHEP::Hep3Vector         tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
      pos       = CLHEP::Hep3Vector(tpos.x(), tpos.y(), tpos.z());
      helix_pos = CLHEP::Hep3Vector(pos);
      robustHel->position(helix_pos);

      double     weight_cl_xy = 1./100.;//FIX ME!
      sxy.addPoint(pos.x(), pos.y(), weight_cl_xy);
	   
    }

    //normalize the chi2
    chi2  = sxy.chi2DofCircle();

    return chi2;
  }
  

//--------------------------------------------------------------------------------
// evaluate the momentum of a Helix
//--------------------------------------------------------------------------------
  double MergeHelixFinder::helix_chi2phiz(const HelixSeed* Helix){
    const RobustHelix*        robustHel = &Helix->helix();
    const HelixHitCollection* hits      = &Helix->hits();
    const HelixHit*           hit(0);
    const mu2e::CaloCluster*  cluster   = Helix->caloCluster().get();
    int                       nhits     = hits->size();
    CLHEP::Hep3Vector         pos(0), wdir(0), sdir(0), helix_center(0);
    double                    phi(0), helix_phi(0);
    double                    chi2(0);
    double                    radius    = robustHel->radius();

    ::LsqSums4 srphi;
    static const CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
    
    helix_center = robustHel->center();

    for (int i=0; i<nhits; ++i){
      hit       = &hits->at(i);
      pos       = hit->pos();
      wdir      = hit->wdir();
      sdir      = zdir.cross(wdir);
      phi       = hit->phi();
      helix_phi = robustHel->fz0() + pos.z()/robustHel->lambda();
      
      double    weight   = calculateWeight(pos, sdir, helix_center, radius, 0);
      double    dPhi     = helix_phi - phi- M_PI/2.;
      while (dPhi > M_PI){
	phi    += 2*M_PI;
        dPhi   = helix_phi - phi;
      }
      while (dPhi < -M_PI){
	phi   -= 2*M_PI; 
	dPhi  = helix_phi - phi;
      }
      srphi.addPoint(pos.z(), phi, weight);
    }

    if (cluster != 0){
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      const mu2e::Calorimeter*  _calorimeter = ch.get();      
      CLHEP::Hep3Vector         gpos = _calorimeter->geomUtil().diskToMu2e(cluster->diskId(),cluster->cog3Vector());
      CLHEP::Hep3Vector         tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
         
      pos       = CLHEP::Hep3Vector(tpos.x(), tpos.y(), tpos.z());
      phi       = CLHEP::Hep3Vector(pos - helix_center).phi();
      phi       = TVector2::Phi_0_2pi(phi);
      helix_phi = robustHel->fz0() + pos.z()/robustHel->lambda();
      double     dPhi        = helix_phi - phi;
      while (dPhi > M_PI){
        dPhi  -= 2*M_PI; 
      }
      while (dPhi < -M_PI){
	dPhi  += 2*M_PI; 
      }
      double     weight_cl_phiz = 10.;
      srphi.addPoint(pos.z(), phi, weight_cl_phiz);
    }

    chi2 = srphi.chi2DofLine();

    return chi2;
  }



//-----------------------------------------------------------------------------
  void MergeHelixFinder::produce(art::Event& AnEvent) {

					// assume less than 100 tracks
    int const   max_ntrk(100);
    int         tpr_flag[max_ntrk], cpr_flag[max_ntrk], ntpr(0), ncpr(0);

    art::Handle<mu2e::HelixSeedCollection>    tpr_h, cpr_h;

    mu2e::HelixSeedCollection  *list_of_helices_tpr(0), *list_of_helices_cpr(0);

    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection );
    unique_ptr<HelixSeedCollection>    helixPtrs(new HelixSeedCollection   );

    if (_debugLevel > 0) ObjectDumpUtils::printEventHeader(&AnEvent,"MergeHelixFinder::produce");

    AnEvent.getByLabel(_trkHelixFinderModuleLabel,tpr_h);
    AnEvent.getByLabel(_calHelixFinderModuleLabel,cpr_h);
    
    if (tpr_h.isValid()) { 
      list_of_helices_tpr = (mu2e::HelixSeedCollection*) &(*tpr_h);
      ntpr              = list_of_helices_tpr->size();
    }

    if (cpr_h.isValid()) {
      list_of_helices_cpr = (mu2e::HelixSeedCollection*) &(*cpr_h);
      ncpr              = list_of_helices_cpr->size();
    }

    for (int i=0; i<max_ntrk; i++) {
      tpr_flag[i] = 1;
      cpr_flag[i] = 1;
    }

    const HelixSeed          *helix_tpr, *helix_cpr;
    //    double                    cpr_mom(0), tpr_mom(0);
    short                     best(-1),  mask;
    AlgorithmID               alg_id;
    HelixHitCollection        tlist, clist;
    int                       nat, nac, natc;
    const mu2e::HelixHit     *hitt, *hitc;
    double                    tpr_chi2xy(0), tpr_chi2phiz(0), tpr_chi2(0);
    double                    cpr_chi2xy(0), cpr_chi2phiz(0), cpr_chi2(0);

    for (int i1=0; i1<ntpr; i1++) {
      helix_tpr    = &list_of_helices_tpr->at(i1);
      //      tpr_mom      = helix_momentum(helix_tpr);
      tpr_chi2xy   = helix_chi2xy  (helix_tpr);
      tpr_chi2phiz = helix_chi2phiz(helix_tpr);
      tpr_chi2     = tpr_chi2xy + tpr_chi2phiz;
      mask         = 1 << AlgorithmID::TrkPatRecBit;
      tlist        = helix_tpr->hits();
      nat          = tlist.size();
      natc         = 0;

      for (int i2=0; i2<ncpr; i2++) {
	helix_cpr    = &list_of_helices_cpr->at(i2);
	cpr_chi2xy   = helix_chi2xy  (helix_tpr);
	cpr_chi2phiz = helix_chi2phiz(helix_tpr);
	cpr_chi2     = cpr_chi2xy + cpr_chi2phiz;
	//	cpr_mom      = helix_momentum(helix_cpr);
	clist        = helix_cpr->hits();
	nac          = clist.size();

//-----------------------------------------------------------------------------
// check the number of common hits: do we need to check also if they have 
// close momentum?
//-----------------------------------------------------------------------------
	for(size_t k=0; k<tlist.size(); ++k){ 
	  hitt = &tlist.at(k);
	  for(size_t l=0; l<clist.size(); l++){ 
	    hitc = &clist.at(l);
	    if (hitt->index() == hitc->index()) {
	      natc += 1;
	      break;
	    }
	  }
	}
//-----------------------------------------------------------------------------
// if > 50% of all hits are common, consider cpr and tpr to be the same
// logic of the choice: 
// 1. take the track which has more hits
// 2. if two tracks have the same number of active hits, choose the one with 
//    best chi2
//-----------------------------------------------------------------------------
	if (natc > (nac+nat)/4.) {

	  mask = mask | (1 << AlgorithmID::CalPatRecBit);

	  if ((tpr_chi2 < _minTprChi2) && (cpr_chi2 < _minCprChi2)) {
//-----------------------------------------------------------------------------
// both tracks are "good", choose the one with best chi2
//-----------------------------------------------------------------------------
	    if (tpr_chi2 >= cpr_chi2) {
	      helixPtrs->push_back(*helix_cpr);
	      best    = AlgorithmID::CalPatRecBit;
	    }
	    else {
	      helixPtrs->push_back(*helix_tpr);
	      best    = AlgorithmID::TrkPatRecBit;
	    }
	  }
	  else if (tpr_chi2 < _minTprChi2) {
//-----------------------------------------------------------------------------
// only TrkHelixFinder track is "good", choose it
//-----------------------------------------------------------------------------
	    helixPtrs->push_back(*helix_tpr);
	    best    = AlgorithmID::TrkPatRecBit; 
	  }
	  else if (cpr_chi2 < _minCprChi2) {
//-----------------------------------------------------------------------------
// only CalHelixFinder track is "good", choose it
//-----------------------------------------------------------------------------
	    helixPtrs->push_back(*helix_cpr);
	    best    = AlgorithmID::CalPatRecBit; 
	  }
	  else {
//-----------------------------------------------------------------------------
// neither track will be selected for analysis, make a choice anyway
//-----------------------------------------------------------------------------
	    if (tpr_chi2 < cpr_chi2) {
	      helixPtrs->push_back(*helix_tpr);
	      best    = AlgorithmID::TrkPatRecBit; 
	    }
	    else {
	      helixPtrs->push_back(*helix_cpr);
	      best    = AlgorithmID::CalPatRecBit; 
	    }
	  }

	  tpr_flag[i1] = 0;
	  cpr_flag[i2] = 0;
	  break;
	}
      }

      if (tpr_flag[i1] == 1) {
	helixPtrs->push_back(*helix_tpr);
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
	helix_cpr = &list_of_helices_cpr->at(i);

	helixPtrs->push_back(*helix_cpr);

	best = AlgorithmID::CalPatRecBit;
	mask = 1 << AlgorithmID::CalPatRecBit;

	alg_id.Set(best,mask);
	algs->push_back(alg_id);
      }
    }

    //    if (_debugLevel > 0) ObjectDumpUtils::printKalRepCollection(&AnEvent,helixPtrs.get(),1);

    AnEvent.put(std::move(helixPtrs));
    AnEvent.put(std::move(algs     ));
  }


//-----------------------------------------------------------------------------
// end job : 
//-----------------------------------------------------------------------------
  void MergeHelixFinder::endJob() {
  }
  
}

using mu2e::MergeHelixFinder;
DEFINE_ART_MODULE(MergeHelixFinder);
