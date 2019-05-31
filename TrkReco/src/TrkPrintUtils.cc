//


#include "Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "TrkReco/inc/TrkPrintUtils.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrkData/inc/TrkCaloHit.hh"

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
//CLHEP

#include "CLHEP/Vector/ThreeVector.h"

//art
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include <iostream>
#include <fstream>

namespace mu2e {

//-----------------------------------------------------------------------------
  TrkPrintUtils::TrkPrintUtils(const fhicl::ParameterSet& PSet) {
    _mcTruth         = PSet.get<int>("mcTruth");
    _strawHitCollTag = PSet.get<std::string>("strawHitCollTag");

    if (_mcTruth != 0) {
      fhicl::ParameterSet ps = PSet.get<fhicl::ParameterSet>("mcUtils");
      _mcUtils = art::make_tool  <McUtilsToolBase>(ps);
    }
    else {
      _mcUtils = std::make_unique<McUtilsToolBase>();
    }
  }

//-----------------------------------------------------------------------------
  TrkPrintUtils::~TrkPrintUtils() {
  }

//-----------------------------------------------------------------------------
// '*' in front of the hit drift radius: the hit drift sign is not defined
//     and the drift ambiguity has been set to 0
// '?': the drift sign determined by the resolver is different from the MC truth
// Option = ""              : print banner + track parameters (but not hits)
//        includes "banner" : print banner
//        includes "data"   : print track parameters
//        includes "hits    : print hits
//-----------------------------------------------------------------------------
  void  TrkPrintUtils::printTrack(const art::Event* AnEvent, const KalRep* Trk, const char* Option, const char* Message) {

    const long unsigned int kNotFound = std::string::npos;

    if (Trk == NULL)  return;

    std::string opt(Option);

    if (Message[0] != 0) printf("[TrkPrintUtils::printTrack] BEGIN called from %s \n",Message);

    int    nhits = Trk->hitVector().size();

    if ((opt == "") || (opt.find("banner") != kNotFound)) {
      printf("-----------------------------------------------------------------------------------------");
      printf("-----------------------------------------------------------------\n");
      printf("  TrkID       Address    N  NA      P       pT     momerr  costh    T0      T0Err   Omega");
      printf("      D0       Z0       Phi0   TanDip    Chi2  MeanRes      FCons\n");
      printf("-----------------------------------------------------------------------------------------");
      printf("-----------------------------------------------------------------\n");
    }

    if ((opt == "") || (opt.find("data") != kNotFound)) {
      CLHEP::Hep3Vector trk_mom;

      double h1_fltlen = Trk->firstHit()->kalHit()->hit()->fltLen() - 10;
      trk_mom          = Trk->momentum(h1_fltlen);
      double mom       = trk_mom.mag();
      double pt        = trk_mom.perp();

      BbrVectorErr      merr   = Trk->momentumErr(h1_fltlen);
      CLHEP::Hep3Vector momdir = trk_mom.unit();

      CLHEP::HepVector momvec(3);
      for (int i=0; i<3; i++) momvec[i] = momdir[i];
    
      double momerr    = sqrt(merr.covMatrix().similarity(momvec));

      double costh     = trk_mom.cosTheta();
      double chi2      = Trk->chisq();
      int    nhits     = Trk->hitVector().size();
      int    nactive   = Trk->nActive();
      double t0        = Trk->t0().t0();
      double t0err     = Trk->t0().t0Err();
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker front plane
//-----------------------------------------------------------------------------
      double s1     = Trk->firstHit()->kalHit()->hit()->fltLen();
      double s2     = Trk->lastHit ()->kalHit()->hit()->fltLen();
      double s      = std::min(s1,s2);

      double d0     = Trk->helix(s).d0();
      double z0     = Trk->helix(s).z0();
      double phi0   = Trk->helix(s).phi0();
      double omega  = Trk->helix(s).omega();
      double tandip = Trk->helix(s).tanDip();

      double fit_consistency = Trk->chisqConsistency().consistency();
      int q         = Trk->charge();

      double sr2 (0);

      for (int i=0; i<nhits; ++i) {
	const TrkStrawHit* hit = dynamic_cast<TrkStrawHit*> (Trk->hitVector().at(i));
	if ((hit == nullptr) || (! hit->isActive())) continue;
//------------------------------------------------------------------------------
// this is an active track straw hit
//-----------------------------------------------------------------------------
	double res, sigres;
	bool hasres = hit->resid(res, sigres, true);
	if (hasres) {
	  sr2 += res*res;
	}
      }

      double mean_res = sqrt(sr2/nactive);

      printf("%5i %16p %3i %3i %8.3f %7.3f  %8.4f %8.4f %7.3f %7.4f",
	     -1,
	     Trk,
	     nhits,
	     nactive,
	     q*mom,pt,momerr,costh,t0,t0err
	     );
      
      printf(" %8.5f %8.3f %8.3f %8.4f %7.4f",omega,d0,z0,phi0,tandip
	     );
      printf(" %8.3f %7.4f %10.3e\n",chi2,mean_res,fit_consistency);
    }

    if (opt.find("hits") == kNotFound) return;
//-----------------------------------------------------------------------------
// print detailed information about the track hits
//-----------------------------------------------------------------------------
    printf("----------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("------------------------------------------------------------------------\n");
    printf("  ih  SId     Flag    A    len         x        y        z      HitT  ");
    printf(" Pl Pn L  W     T0       Xs      Ys        Zs      resid  sigres");
    printf("   Rdrift   mcdoca  totErr rdrErr  t0Err penErr extErr  vinst      simID\n");
    printf("----------------------------------------------------------------------");
    printf("----------------------------------------------------------------");
    printf("------------------------------------------------------------------------\n");

    TrkStrawHit          *hit;
    CLHEP::Hep3Vector     pos;
    const ComboHit        *sh;
    const Straw           *straw;
    int                   ihit;
    double                len;
    HepPoint              plen;


    auto shcH = AnEvent->getValidHandle<ComboHitCollection>(_strawHitCollTag);
    const ComboHitCollection* shcol = shcH.product();

    ihit = 0;
    for (int it=0; it<nhits; ++it) {
      hit   = dynamic_cast<TrkStrawHit*> (Trk->hitVector().at(it));
      if (hit != 0) {
//------------------------------------------------------------------------------
// this is a track hit
//-----------------------------------------------------------------------------
	sh    = &hit->comboHit();  // in reality, 'sh' is a single straw hit
	straw = &hit->straw();
//-----------------------------------------------------------------------------
// as the straw hit doesn't store its index in a collection, determine it 
// assuming that the straw hit collection in a vector
//-----------------------------------------------------------------------------
	int loc = sh-&shcol->at(0);

	hit->hitPosition(pos);

	len   = hit->fltLen();
	plen  = Trk->position(len);

	int    sim_id = _mcUtils->strawHitSimId(AnEvent,loc);
	double mcdoca = _mcUtils->mcDoca       (AnEvent,hit); // loc,straw);

	ihit += 1;
	printf("%4i %5i 0x%08x %1i %9.3f %8.3f %8.3f %9.3f %8.3f",
	    ihit,
	    straw->id().asUint16(),
	    hit->hitFlag(),
	    hit->isActive(),
	    len,
	    plen.x(),plen.y(),plen.z(),
	    sh->time()
	    );

	printf(" %2i %2i %1i %2i",
	    straw->id().getPlane(),
	    straw->id().getPanel(),
	    straw->id().getLayer(),
	    straw->id().getStraw()
	    );

	printf(" %8.3f",hit->hitT0().t0());

	double res, sigres;
	bool hasres = hit->resid(res, sigres, true);

	if (hasres) {
	  printf(" %8.3f %8.3f %9.3f %7.3f %7.3f",
	      pos.x(),
	      pos.y(),
	      pos.z(),
	      res,
	      sigres
	      );
	}
	else {
	  printf(" %8.3f %8.3f %9.3f %7s %7s",
	      pos.x(),
	      pos.y(),
	      pos.z(),
	      "   -   ",
	      "   -   "
	      );
	}

	if      (hit->ambig()       == 0) printf(" * %6.3f",hit->driftRadius());
	else if (hit->ambig()*mcdoca > 0) printf("   %6.3f",hit->driftRadius()*hit->ambig());
	else                              printf(" ? %6.3f",hit->driftRadius()*hit->ambig());

	float extErr(-1.);
	if (hit->isActive()) extErr = hit->temperature()*0.0625; // external error, assume same for all hits

	printf("  %7.3f  %6.3f %6.3f %6.3f %6.3f %6.3f %6.4f %10i\n",
	       mcdoca,
	       hit->totalErr(),
	       hit->driftRadiusErr(),
	       hit->t0Err(),
	       hit->penaltyErr(),
	       extErr,                     
	       hit->driftVelocity(),
	       sim_id
	       );
      }
      else {
	TrkCaloHit const* chit   = dynamic_cast<TrkCaloHit*> (Trk->hitVector().at(it));
	if (chit != 0) {
//-----------------------------------------------------------------------------
// calorimeter hit
//-----------------------------------------------------------------------------
	  double res, sigres;
	  bool hasres = chit->resid(res, sigres, true);
	  printf("TrkCaloHit, time = %10.3f hitT0 = %10.3g +/- %10.3f",
		 chit->time(), chit->hitT0().t0(),chit->hitT0().t0Err());

	  if (hasres) printf(" resid = %10.3f +- %10.3f\n",res,sigres);
	  else	      printf(" no residual, hit error = %10.3f\n",chit->hitErr());
	}
      }
    }
   }
  
}
