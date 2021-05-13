///////////////////////////////////////////////////////////////////////////////
// takes inputs from two track finding algorithms, produces one track collection
// on output to be used for analysis
//
// development history:
// --------------------
// calPatRecMVAType = 0 : use log10(fitCons)
//                  = 1 : use chi2/ndof
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

#include "TROOT.h"
#include "TFolder.h"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/KalSeed.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TrkDiag/inc/KalDiag.hh"
#include "BTrkData/inc/Doublet.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"

#include "CalPatRec/inc/ObjectDumpUtils.hh"

#include "CalPatRec/inc/MergePatRec_types.hh"
#include "RecoDataProducts/inc/AlgorithmIDCollection.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "art/Utilities/make_tool.h"

// Xerces XML Parser
#include <xercesc/dom/DOM.hpp>

#include "Mu2eUtilities/inc/MVATools.hh"
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
  using namespace MergePatRecTypes;

  class MergePatRec : public art::EDProducer {
  public:
    explicit MergePatRec(fhicl::ParameterSet const&);
    virtual ~MergePatRec();
    void beginJob() override;
    void beginRun(art::Run&) override;
    void produce(art::Event& event) override;

    double         s_at_given_z(const KalRep* Krep, double Z);

    class ProbDist {
    public:

      static int fgIndex;

      TH1F* _h;
      TH1F* _hprob;

      ProbDist();
      ProbDist(TH1F* Hist);
      ProbDist(const char* Fn);
      
      ~ProbDist();
      
      double prob(double X);
      int    readHistogram(const char* Fn, TH1F** Hist);
    };


  private:
    unsigned         _iev;
                                        // configuration parameters
    int              _diagLevel;
    int              _debugLevel;
    int              _printfreq;
    bool             _addhits;
    TrkParticle      _tpart;            // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search
                                        // event object tokens
    art::ProductToken<KalRepPtrCollection> const _tprToken;
    art::ProductToken<KalRepPtrCollection> const _cprToken;
    art::ProductToken<KalSeedCollection> const _stprToken;
    art::ProductToken<KalSeedCollection> const _scprToken;

    float            _minTprQual;
    std::string      _trkPatRecMVAHist ;  // in .tab format

    int              _calPatRecMVAType ;  //
    float            _minCprQual;
    std::string      _calPatRecMVAHist ;  // in .tab format

    KalDiag*         _kalDiag;
    MVATools*        _calPatRecQualMVA;

    DoubletAmbigResolver* _dar;

    ProbDist*        _trkPatRecProb;
    ProbDist*        _calPatRecProb;

    Data_t                              _data;              // all data used
    std::unique_ptr<ModuleHistToolBase> _hmanager;
  };

  MergePatRec::MergePatRec(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _diagLevel              (pset.get<int>("diagLevel" )),
    _debugLevel             (pset.get<int>("debugLevel")),
    _tprToken{consumes<KalRepPtrCollection>(pset.get<std::string>("trkPatRecModuleLabel"))},
    _cprToken{consumes<KalRepPtrCollection>(pset.get<std::string>("calPatRecModuleLabel"))},
    _stprToken{consumes<KalSeedCollection>(pset.get<std::string>("trkPatRecModuleLabel"))},
    _scprToken{consumes<KalSeedCollection>(pset.get<std::string>("calPatRecModuleLabel"))}
  {

    produces<AlgorithmIDCollection>  ();
    produces<KalRepPtrCollection>    ();
    produces<TrkQualCollection>      ();
    produces<KalSeedCollection>      ();

    _kalDiag = new KalDiag(pset.get<fhicl::ParameterSet>("KalDiag", {}));
    _dar     = new DoubletAmbigResolver (pset.get<fhicl::ParameterSet>("DoubletAmbigResolver"),0.,0,0);

    fhicl::ParameterSet pset_cpr = pset.get<fhicl::ParameterSet>("calPatRecMVA", {});

    _calPatRecQualMVA = new mu2e::MVATools(pset_cpr);
    _calPatRecQualMVA->initMVA();

    ConfigFileLookupPolicy configFile;

    _calPatRecMVAType = pset_cpr.get<int>        ("MVAType");
    _minCprQual       = pset_cpr.get<float>      ("minTrkQual");
    _calPatRecMVAHist = configFile(pset_cpr.get<std::string>("trkQualHist"));

    fhicl::ParameterSet pset_tpr = pset.get<fhicl::ParameterSet>("trkPatRecMVA", {});
    _minTprQual       = pset_tpr.get<float>      ("minTrkQual");
    _trkPatRecMVAHist = configFile(pset_tpr.get<std::string>("trkQualHist"));
//-----------------------------------------------------------------------------
// probability distributions for TrkPatRec and CalPatRec tracks
//-----------------------------------------------------------------------------
    _calPatRecProb    = new ProbDist(_calPatRecMVAHist.data());
    _trkPatRecProb    = new ProbDist(_trkPatRecMVAHist.data());

    if (_diagLevel != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();
  }

  MergePatRec::~MergePatRec() {
    delete _calPatRecQualMVA;
    delete _kalDiag;
  }

  void MergePatRec::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  void MergePatRec::beginRun(art::Run& ) {
  }

//-----------------------------------------------------------------------------
// extrapolate track to a given Z
//-----------------------------------------------------------------------------
  double MergePatRec::s_at_given_z(const KalRep* Krep, double Z) {

    double  ds(10.), s0, s1(1.e6), s2(1.e6), z0, z1, z2, dzds, sz, sz1, z01;

    const TrkHit  *first(nullptr), *last(nullptr); 

    const TrkHitVector* hots = &Krep->hitVector();
    int nh = hots->size();

    for (int ih=0; ih<nh; ++ih) {
      const TrkHit* hit =  dynamic_cast<const TrkHit*> (hots->at(ih));
      if (hit   != nullptr) {
	if (first == nullptr) first = hit;
	last = hit;
      }
    }

    if (first) s1     = first->fltLen();
    if (last ) s2     = last ->fltLen();
    z1     = Krep->position(s1).z();
    z2     = Krep->position(s2).z();

    dzds   = (z2-z1)/(s2-s1);
//-----------------------------------------------------------------------------
// iterate once, choose the closest point
//-----------------------------------------------------------------------------
    if (fabs(Z-z1) > fabs(Z-z2)) {
      z0 = z2;
      s0 = s2;
    }
    else {
      z0 = z1;
      s0 = s1;
    }

    sz    = s0+(Z-z0)/dzds;

    z0     = Krep->position(sz).z();     // z0 has to be close to Z(TT_FrontPA)
    z01    = Krep->position(sz+ds).z();

    dzds   = (z01-z0)/ds;
    sz1    = sz+(Z-z0)/dzds;              // should be good enough

    return sz1;
  }

//-----------------------------------------------------------------------------
  void MergePatRec::produce(art::Event& AnEvent) {

                                        // assume less than 100 tracks
    int const   max_ntrk(100);
    int         tpr_flag[max_ntrk], cpr_flag[max_ntrk], ntpr(0), ncpr(0);

    art::Handle<mu2e::KalRepPtrCollection>    htpr, hcpr;
    art::Handle<mu2e::KalSeedCollection>      hstpr, hscpr;

    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

    auto algs = std::make_unique<AlgorithmIDCollection>();
    auto trackPtrs = std::make_unique<KalRepPtrCollection>();
    auto tqcol = std::make_unique<TrkQualCollection>();
    auto kscol = std::make_unique<KalSeedCollection>();

    if (_debugLevel > 0) ObjectDumpUtils::printEventHeader(&AnEvent,"MergePatRec::produce");

    htpr = AnEvent.getHandle<mu2e::KalRepPtrCollection>(_tprToken);
    hcpr = AnEvent.getHandle<mu2e::KalRepPtrCollection>(_cprToken);

    hstpr = AnEvent.getHandle<KalSeedCollection>(_stprToken);
    hscpr = AnEvent.getHandle<KalSeedCollection>(_scprToken);

    _data.event = &AnEvent;
    _data.list_of_kreps_tpr = nullptr;
    _data.list_of_kreps_cpr = nullptr;
    _data.list_of_kseed_tpr = nullptr;
    _data.list_of_kseed_cpr = nullptr;

    if (htpr.isValid()) {
      _data.list_of_kreps_tpr = htpr.product();
      ntpr                    = _data.list_of_kreps_tpr->size();
    }

    if (hcpr.isValid()) {
      _data.list_of_kreps_cpr = hcpr.product();
      ncpr                    = _data.list_of_kreps_cpr->size();
    }

    if (hstpr.isValid()) {
      _data.list_of_kseed_tpr = hstpr.product();
    }

    if (hscpr.isValid()) {
      _data.list_of_kseed_cpr = hscpr.product();
    }

    for (int i=0; i<max_ntrk; i++) {
      tpr_flag[i] = 1;
      cpr_flag[i] = 1;
    }

    const art::Ptr<KalRep>    *tpr, *cpr;
    const KalRep              *krep_tpr, *krep_cpr;
    const KalSeed             *tpr_kseed, *cpr_kseed;
    Hep3Vector                cpr_mom, tpr_mom;
    short                     best(-1),  mask;
    double                    best_MVAQual(-1);
    AlgorithmID               alg_id;
    TrkHitVector              tlist, clist;
    int                       nat, nac, natc;
    const mu2e::TrkStrawHit  *hitt, *hitc;
    double                    tpr_qual, cpr_qual;

    TrkInfo                   tpr_trkinfo, cpr_trkinfo;

    for (int i1=0; i1<ntpr; i1++) {
      tpr_kseed = &_data.list_of_kseed_tpr->at(i1);
      tpr       = &_data.list_of_kreps_tpr->at(i1);
      tpr_mom   = (*tpr)->momentum();
      //      tpr_chisq = (*tpr)->chisq();
      mask      = 1 << AlgorithmID::TrkPatRecBit;
      tlist     = (*tpr)->hitVector();
      nat       = (*tpr)->nActive();
      natc      = 0;

      krep_tpr = tpr->get();
      _kalDiag->kalDiag(krep_tpr,false);

      tpr_trkinfo = _kalDiag->_trkinfo;
      tpr_qual    = tpr_trkinfo._trkqual;

      for (int i2=0; i2<ncpr; i2++) {
        cpr_kseed = &_data.list_of_kseed_cpr->at(i2);
        cpr       = &_data.list_of_kreps_cpr->at(i2);
        krep_cpr  = cpr->get();

        _kalDiag->kalDiag(krep_cpr,false);
        cpr_trkinfo = _kalDiag->_trkinfo;

        cpr_mom   = krep_cpr->momentum();
        clist     = krep_cpr->hitVector();
        nac       = krep_cpr->nActive();

        vector<double>  pmva;
        pmva.resize(10);

        float na = nac;

        pmva[ 0] = na;
        pmva[ 1] = na/krep_cpr->nHits();

        if      (_calPatRecMVAType == 0) pmva[2] = log10(krep_cpr->chisqConsistency().consistency());
        else                             pmva[2] = krep_cpr->chisq()/(na-5.);
//-----------------------------------------------------------------------------
// determine, approximately, 'sz0' - flight length corresponding to the
// virtual detector at the tracker entrance
// 2019-06-21 P.Murat: now, with the introduction of TrkCaloHit, not all hits 
//            have the fltlen() defined
//-----------------------------------------------------------------------------
        double  h1_fltlen(1.e6), hn_fltlen(1.e6), entlen;

	// const TrkStrawHit* h1 = dynamic_cast<const TrkStrawHit*> (krep_cpr->firstHit()->kalHit()->hit());
	// const TrkStrawHit* h2 = dynamic_cast<const TrkStrawHit*> (krep_cpr->lastHit ()->kalHit()->hit());

	const TrkHit  *first(nullptr), *last(nullptr); 

	const TrkHitVector* hots = &krep_cpr->hitVector();
	int nh = hots->size();

	for (int ih=0; ih<nh; ++ih) {
	  const TrkHit* hit =  dynamic_cast<const TrkHit*> (hots->at(ih));
	  if (hit   != nullptr) {
	    if (first == nullptr) first = hit;
	    last = hit;
	  }
	}

        if (first) h1_fltlen = first->fltLen();
        if (last ) hn_fltlen = last->fltLen();

        entlen         = std::min(h1_fltlen,hn_fltlen);

        CLHEP::Hep3Vector fitmom = krep_cpr->momentum(entlen);
//-----------------------------------------------------------------------------
// pmva[3]: momentum error in the first point
// for consistency, use helical parameterization in the first point
//-----------------------------------------------------------------------------
	CLHEP::Hep3Vector momdir = fitmom.unit();
	BbrVectorErr      momerr = krep_cpr->momentumErr(entlen);
	
	CLHEP::HepVector momvec(3);
	for (int i=0; i<3; i++) momvec[i] = momdir[i];
	
	pmva[ 3] = sqrt(momerr.covMatrix().similarity(momvec));
	pmva[ 4] = krep_cpr->t0().t0Err();

	Hep3Vector tfront = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_FrontPA));
	double     zfront = tfront.z();
	double     sz0    = s_at_given_z(krep_cpr,zfront);

	HelixParams helx  = krep_cpr->helix(sz0);

	pmva[ 5] = helx.d0();
	pmva[ 6] = helx.d0()+2/helx.omega();

					// calculate number of doublets 

	vector<mu2e::Doublet> list_of_doublets;
	_dar->findDoublets(krep_cpr,&list_of_doublets);
//-----------------------------------------------------------------------------
// counting only 2+ hit doublets
//-----------------------------------------------------------------------------
	mu2e::Doublet*                     d;
	
	int   nd_tot(0), nd_os(0), nd_ss(0), ns;

	int nad(0);                     // number of doublets with at least one hit active
	
	int nd = list_of_doublets.size();
	for (int i=0; i<nd; i++) {
	  d  = (mu2e::Doublet*) &list_of_doublets.at(i);
	  ns = d->fNStrawHits;
	  
	  if (ns > 1) { 
	    nd_tot += 1;
	    if (d->isSameSign()) nd_ss += 1;
	    else                 nd_os += 1;

	    int active = 1;
	    for (int is=0; is<ns; is++) {
	      if (!d->fHit[is]->isActive()) {
		active = 0;
		break;
	      }
	    }
	    
	    if (active == 1) {
	      nad += 1;
	    }
	  }
	}
	
	pmva[ 7] = nad/na;
	pmva[ 8] = cpr_trkinfo._nnullambig/na;
	pmva[ 9] = cpr_trkinfo._nmatactive/na;

	cpr_qual = _calPatRecQualMVA->evalMVA(pmva);
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
		if (&hitt->comboHit() == &hitc->comboHit()) {
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
// 2. if two tracks have the same number of active hits, choose the one with 
//    "higher probability"
//-----------------------------------------------------------------------------
	if (natc > (nac+nat)/4.) {

	  mask = mask | (1 << AlgorithmID::CalPatRecBit);

	  double tpr_prob = _trkPatRecProb->prob(tpr_qual);
	  double cpr_prob = _calPatRecProb->prob(cpr_qual);

          if ((tpr_trkinfo._trkqual > _minTprQual) && (cpr_qual > _minCprQual)) {
//-----------------------------------------------------------------------------
// both tracks are "good", choose the one with higher probability
//-----------------------------------------------------------------------------
            if (tpr_prob >= cpr_prob) {
              trackPtrs->push_back(*tpr);
              kscol    ->push_back(*tpr_kseed);
              best    = AlgorithmID::TrkPatRecBit;
              best_MVAQual = tpr_qual;
            }
            else {
              trackPtrs->push_back(*cpr);
              kscol    ->push_back(*cpr_kseed);
              best    = AlgorithmID::CalPatRecBit;
              best_MVAQual = cpr_qual;
            }
          }
          else if (tpr_trkinfo._trkqual > _minTprQual) {
//-----------------------------------------------------------------------------
// only TrkPatRec track is "good", choose it
//-----------------------------------------------------------------------------
            trackPtrs->push_back(*tpr);
            kscol    ->push_back(*tpr_kseed);
            best    = AlgorithmID::TrkPatRecBit;
            best_MVAQual = tpr_qual;
          }
          else if (cpr_qual > _minCprQual) {
//-----------------------------------------------------------------------------
// only CalPatRec track is "good", choose it
//-----------------------------------------------------------------------------
            trackPtrs->push_back(*cpr);
            kscol    ->push_back(*cpr_kseed);
            best    = AlgorithmID::CalPatRecBit;
            best_MVAQual = cpr_qual;
          }
          else {
//-----------------------------------------------------------------------------
// neither track will be selected for analysis, make a choice anyway
//-----------------------------------------------------------------------------
            if (tpr_prob >= cpr_prob) {
              trackPtrs->push_back(*tpr);
              kscol    ->push_back(*tpr_kseed);
              best    = AlgorithmID::TrkPatRecBit;
              best_MVAQual = tpr_qual;
            }
            else {
              trackPtrs->push_back(*cpr);
              kscol    ->push_back(*cpr_kseed);
              best    = AlgorithmID::CalPatRecBit;
              best_MVAQual = cpr_qual;
            }
          }

          tpr_flag[i1] = 0;
          cpr_flag[i2] = 0;
          break;
        }
      }

      if (tpr_flag[i1] == 1) {
        trackPtrs->push_back(*tpr);
        kscol    ->push_back(*tpr_kseed);
        best = AlgorithmID::TrkPatRecBit;
        best_MVAQual = tpr_qual;
      }

      alg_id.Set(best,mask);
      algs->push_back(alg_id);

      // compute TrkQual for this track and save it
      // DUMMY FILLING. FIXME!
      TrkQual trkqual;
      trkqual.setMVAValue(best_MVAQual);

      tqcol->push_back(trkqual);
    }
//-----------------------------------------------------------------------------
// account for presence of multiple tracks
//-----------------------------------------------------------------------------
    for (int i=0; i<ncpr; i++) {
      if (cpr_flag[i] == 1) {
        cpr = &_data.list_of_kreps_cpr->at(i);
        cpr_kseed = &_data.list_of_kseed_cpr->at(i);

        trackPtrs->push_back(*cpr);
        kscol    ->push_back(*cpr_kseed);

        best = AlgorithmID::CalPatRecBit;
        mask = 1 << AlgorithmID::CalPatRecBit;

        alg_id.Set(best,mask);
        algs->push_back(alg_id);

 // compute TrkQual for this track and save it
 // DUMMY FILLING. FIXME!
        krep_cpr  = cpr->get();	_kalDiag->kalDiag(krep_cpr,false);
        cpr_trkinfo = _kalDiag->_trkinfo;

        cpr_mom   = krep_cpr->momentum();
        clist     = krep_cpr->hitVector();
        nac       = krep_cpr->nActive();

        vector<double>  pmva;
        pmva.resize(10);

        float na = nac;

        pmva[ 0] = na;
        pmva[ 1] = na/krep_cpr->nHits();

        if      (_calPatRecMVAType == 0) pmva[2] = log10(krep_cpr->chisqConsistency().consistency());
        else                             pmva[2] = krep_cpr->chisq()/(na-5.);
//-----------------------------------------------------------------------------
// determine, approximately, 'sz0' - flight length corresponding to the
// virtual detector at the tracker entrance
//-----------------------------------------------------------------------------
	double  h1_fltlen(1.e6), hn_fltlen(1.e6), entlen;

	// const TrkStrawHit* h1 = dynamic_cast<const TrkStrawHit*> (krep_cpr->firstHit()->kalHit()->hit());
	// const TrkStrawHit* h2 = dynamic_cast<const TrkStrawHit*> (krep_cpr->lastHit ()->kalHit()->hit());

	// if (h1) h1_fltlen = h1->fltLen();
	// if (h2) hn_fltlen = h2->fltLen();

	const TrkHit  *first(nullptr), *last(nullptr); 

	const TrkHitVector* hots = &krep_cpr->hitVector();
	int nh = hots->size();

	for (int ih=0; ih<nh; ++ih) {
	  const TrkHit* hit =  dynamic_cast<const TrkHit*> (hots->at(ih));
	  if (hit   != nullptr) {
	    if (first == nullptr) first = hit;
	    last = hit;
	  }
	}

        if (first) h1_fltlen = first->fltLen();
        if (last ) hn_fltlen = last->fltLen();
	entlen         = std::min(h1_fltlen,hn_fltlen);

	CLHEP::Hep3Vector fitmom = krep_cpr->momentum(entlen);
//-----------------------------------------------------------------------------
// pmva[3]: momentum error in the first point
// for consistency, use helical parameterization in the first point
//-----------------------------------------------------------------------------
        CLHEP::Hep3Vector momdir = fitmom.unit();
        BbrVectorErr      momerr = krep_cpr->momentumErr(entlen);

        CLHEP::HepVector momvec(3);
        for (int i=0; i<3; i++) momvec[i] = momdir[i];

        pmva[ 3] = sqrt(momerr.covMatrix().similarity(momvec));
        pmva[ 4] = krep_cpr->t0().t0Err();

        Hep3Vector tfront = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_FrontPA));
        double     zfront = tfront.z();
        double     sz0    = s_at_given_z(krep_cpr,zfront);

        HelixParams helx  = krep_cpr->helix(sz0);

        pmva[ 5] = helx.d0();
        pmva[ 6] = helx.d0()+2/helx.omega();

                                        // calculate number of doublets

        vector<mu2e::Doublet> list_of_doublets;
        _dar->findDoublets(krep_cpr,&list_of_doublets);

//-----------------------------------------------------------------------------
// counting only 2+ hit doublets
//-----------------------------------------------------------------------------
        mu2e::Doublet*                     d;
        //	mu2e::DoubletAmbigResolver::Data_t r;

        int   nd_tot(0), nd_os(0), nd_ss(0), ns;

        int nad(0);                     // number of doublets with at least one hit active

        int nd = list_of_doublets.size();
        for (int i=0; i<nd; i++) {
          d  = (mu2e::Doublet*) &list_of_doublets.at(i);
          ns = d->fNStrawHits;

          if (ns > 1) {
            nd_tot += 1;
            if (d->isSameSign()) nd_ss += 1;
            else                 nd_os += 1;

            int active = 1;
            for (int is=0; is<ns; is++) {
              if (!d->fHit[is]->isActive()) {
                active = 0;
                break;
              }
            }

            if (active == 1) {
              nad += 1;
            }
          }
        }

        pmva[ 7] = nad/na;
        pmva[ 8] = cpr_trkinfo._nnullambig/na;
        pmva[ 9] = cpr_trkinfo._nmatactive/na;

        double  cpr_MVAQual = _calPatRecQualMVA->evalMVA(pmva);
        TrkQual trkqual;
        trkqual.setMVAValue(cpr_MVAQual);

        tqcol->push_back(trkqual);
      }
    }

    if (_debugLevel > 0) ObjectDumpUtils::printKalRepCollection(&AnEvent,trackPtrs.get(),1);

    AnEvent.put(std::move(trackPtrs));
    AnEvent.put(std::move(algs     ));
    AnEvent.put(std::move(tqcol    ));
    AnEvent.put(std::move(kscol    ));
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    if (_diagLevel  > 0) _hmanager->fillHistograms(&_data);
    if (_debugLevel > 0) _hmanager->debug(&_data);
  }

  int MergePatRec::ProbDist::fgIndex(0);

//-----------------------------------------------------------------------------
  MergePatRec::ProbDist::ProbDist() {
    _h     = nullptr;
    _hprob = nullptr;
  }


//-----------------------------------------------------------------------------
  MergePatRec::ProbDist::ProbDist(TH1F* Hist) {

    _h     = (TH1F*) Hist->Clone(Form("hprob_dist_h_%i"  ,fgIndex));

    _hprob = (TH1F*) _h->Clone(Form("hprob_dist_hprob_%i",fgIndex));
    fgIndex += 1;

    _hprob->Reset();

    int nb = _h->GetNbinsX();

    double anorm = _h->Integral(1,nb);

    for (int i=1; i<=nb; i++) {
      double prob = _h->Integral(1,i)/anorm;
      _hprob->SetBinContent(i,prob);
    }
  }


//-----------------------------------------------------------------------------
// read histogram from a text file
//-----------------------------------------------------------------------------
  MergePatRec::ProbDist::ProbDist(const char* Fn) {

    _h = nullptr;
    readHistogram(Fn,&_h);
    _hprob =  (TH1F*) _h->Clone(Form("h_ProbDist_hprob_%i"  ,fgIndex));
    fgIndex += 1;

    _hprob->Reset();

    int nb = _h->GetNbinsX();

    double anorm = _h->Integral(1,nb);

    for (int i=1; i<=nb; i++) {
      double prob = _h->Integral(1,i)/anorm;
      _hprob->SetBinContent(i,prob);
    }
  }

//-----------------------------------------------------------------------------
  int MergePatRec::ProbDist::readHistogram(const char* Fn, TH1F** Hist) {
    FILE  *f;
    int    done = 0, nbx, loc(0), ix, line(0);
    char   c[1000], title[200], name[200];
    float  val, xmin, xmax;

    f = fopen(Fn,"r");
    if (f == 0) {
      Error("TEmuLogLH::ReadHistogram",Form("missing file %s\n",Fn));
      return -2;
    }

    if ((*Hist) != nullptr) delete (*Hist);

    while ( ((c[0]=getc(f)) != EOF) && !done) {

                                        // check if it is a comment line
      if (c[0] != '#') {
        ungetc(c[0],f);

        if (line == 0) {
          fscanf(f,"title: %s" ,title);
          line++;
        }
        else if (line == 1) {
          fscanf(f,"name: %s"  ,name);
          line++;
        }
        else if (line ==2) {
          fscanf(f,"nbx,xmin,xmax: %i %f %f"  ,&nbx,&xmin,&xmax);
          *Hist = new TH1F(name,title,nbx,xmin,xmax);
          line++;
        }
        else {
          for (int i=0; i<10; i++) {
            fscanf(f,"%f" ,&val);
            ix = loc + 1;
            (*Hist)->SetBinContent(ix,val);
            loc++;
          }
          line++;
        }
      }
                                        // skip the rest of the line
      fgets(c,100,f);

    }

    fclose(f);
    return 0;
  }



//-----------------------------------------------------------------------------
  double MergePatRec::ProbDist::prob(double X) {
    double f(0);

    int nb = _h->GetNbinsX();
                                        // assume all bins are the same
    double binw = _h->GetBinWidth(1);

    if      (X <  _hprob->GetBinCenter( 1)-binw/2) return 0.;
    else if (X >= _hprob->GetBinCenter(nb)+binw/2) return 1.;

    for (int i=1; i<=nb; i++) {
      double x = _h->GetBinCenter(i);
      if (x+binw/2 > X) {
        f = _hprob->GetBinContent(i);
        break;
      }
    }

    return f;
  }

}

using mu2e::MergePatRec;
DEFINE_ART_MODULE(MergePatRec);
