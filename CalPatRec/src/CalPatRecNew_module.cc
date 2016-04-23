///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Pattern recognition only, passes results to CalTrkFit
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalPatRecNew_module.hh"
#include "CalPatRec/inc/Ref.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "KalmanTests/inc/KalFitResult.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "TVector2.h"

using namespace std;
using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {
//-----------------------------------------------------------------------------
// comparison functor for sorting by Z(wire)
//-----------------------------------------------------------------------------
  struct straw_zcomp : public binary_function<hitIndex,hitIndex,bool> {
    bool operator()(hitIndex const& h1, hitIndex const& h2) {

      mu2e::GeomHandle<mu2e::TTracker> handle;
      const TTracker* t = handle.get();
      const Straw* s1 = &t->getStraw(StrawIndex(h1._index));
      const Straw* s2 = &t->getStraw(StrawIndex(h2._index));

      return s1->getMidPoint().z() < s2->getMidPoint().z();
    }
  }; // a semicolumn here is required

//-----------------------------------------------------------------------------
// module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
//-----------------------------------------------------------------------------
  CalPatRecNew::CalPatRecNew(fhicl::ParameterSet const& pset) :
    _diagLevel   (pset.get<int>   ("diagLevel"                      )),
    _debugLevel  (pset.get<int>   ("debugLevel"                     )),
    _printfreq   (pset.get<int>   ("printFrequency"                 )),
    _shLabel     (pset.get<string>("StrawHitCollectionLabel"        )),
    _shpLabel    (pset.get<string>("StrawHitPositionCollectionLabel")),
    _shfLabel    (pset.get<string>("StrawHitFlagCollectionLabel"    )),
    _trkseedLabel(pset.get<string>("TrackSeedCollectionLabel"       )),
    _tpart       ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir        ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _hfit        (pset.get<fhicl::ParameterSet>("HelixFitHack",fhicl::ParameterSet()))
  {
    produces<TrackSeedCollection>();
    //    produces<CalTimePeakCollection>();

    fHackData = new THackData("HackData","Hack Data");
    gROOT->GetRootFolder()->Add(fHackData);
//-----------------------------------------------------------------------------
// provide for interactive disanostics
//-----------------------------------------------------------------------------
    _ref = new Ref("CalPatRecNewRef","Ref to CalPatRecNew",&_hfit);

    TFolder* f_mu2e;

    f_mu2e = (TFolder*) gROOT->GetRootFolder()->FindObject("Mu2e");
    if (f_mu2e == NULL) f_mu2e = gROOT->GetRootFolder()->AddFolder("Mu2e","Mu2e Folder");

    if (f_mu2e) {
      _folder = f_mu2e->AddFolder("CalPatRecNew","CalPatRecNew Folder");
      _folder->Add(_ref);
    }

    //    fgTimeOffsets     = new SimParticleTimeOffset(pset.get<fhicl::ParameterSet>("TimeOffsets"));

    _helTraj = 0;
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalPatRecNew::~CalPatRecNew() {
    delete _ref;
    if (_helTraj) delete _helTraj;
    //    delete fStopwatch;
  }

//-----------------------------------------------------------------------------
  void CalPatRecNew::beginJob(){

    if(_diagLevel > 0) bookHistograms();

  }

//-----------------------------------------------------------------------------
  void CalPatRecNew::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
					// calibrations

    mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    _trackerCalib = tcal.operator ->();

  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRecNew::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;

    art::TFileDirectory hf_dir = tfs->mkdir("HelixFit");
    
    _hist.nseeds[0]            = tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
    _hist.nseeds[1]            = tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15", 21, -0.5, 20.5);
    _hist.helixFit.nhits       = hf_dir.make<TH1F>("nhits" , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
    _hist.helixFit.radius[0]   = hf_dir.make<TH1F>("radius0", "helix radius; r [mm]"                  , 401, -0.5, 400.5);
    _hist.helixFit.radius[1]   = hf_dir.make<TH1F>("radius1", "helix radius nhits > 15; r [mm]"       , 401, -0.5, 400.5);
    _hist.helixFit.pT [0]      = hf_dir.make<TH1F>("pT0"    , "transverse momentum; pT [MeV/c]"       , 400, -0.5, 200.5);
    _hist.helixFit.p  [0]      = hf_dir.make<TH1F>("p0"     , "momentum; p [MeV/c]"                   , 400, -0.5, 200.5);
    _hist.helixFit.pT [1]      = hf_dir.make<TH1F>("pT1"    , "transverse momentum nhits > 15; pT [MeV/c]"       , 400, -0.5, 200.5);
    _hist.helixFit.p  [1]      = hf_dir.make<TH1F>("p1"     , "momentum nhits > 15; p [MeV/c]"                   , 400, -0.5, 200.5);
    _hist.helixFit.chi2XY[0]   = hf_dir.make<TH1F>("chi2XY0", "normalized chi2-XY"                   , 200, 0., 20.);
    _hist.helixFit.chi2XY[1]   = hf_dir.make<TH1F>("chi2XY1", "normalized chi2-XY: nhits>15"         , 200, 0., 20.);
    _hist.helixFit.chi2ZPhi[0] = hf_dir.make<TH1F>("chi2ZPhi0", "normalized chi2-ZPhi"             , 200, 0., 20.);
    _hist.helixFit.chi2ZPhi[1] = hf_dir.make<TH1F>("chi2ZPhi1", "normalized chi2-ZPhi: nhits>15"   , 200, 0., 20.);
    _hist.helixFit.nhitsvspT   = hf_dir.make<TH2F>("nhitsvspT","nhits vs pT", 100, 0, 100, 400, 0, 200);
    _hist.helixFit.nhitsvsp    = hf_dir.make<TH2F>("nhitsvsp" ,"nhits vs p" , 100, 0, 100, 400, 0, 200);

  }

//-----------------------------------------------------------------------------
// find the input data objects
//-----------------------------------------------------------------------------
  bool CalPatRecNew::findData(const art::Event& evt) {

    //    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel, _strawhitsH)) {
      _shcol = _strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalPatRecNew::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (evt.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CalPatRecNew::findData: StrawHitPositionCollection with label=%s not found.\n",
             _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalPatRecNew::findData: StrawHitFlagCollection with label=%s not found.\n",
             _shfLabel.data());
    }

    art::Handle<mu2e::TrackSeedCollection> trkseedsH;
    if (evt.getByLabel(_trkseedLabel, trkseedsH)) {
      _trkSeeds = trkseedsH.product();
    }
    else {
      _trkSeeds = 0;
      printf(" >>> ERROR in CalPatRecNew::findData: TrackSeedCollection with label=%s not found.\n",
             _trkseedLabel.data());
    }

    art::Handle<mu2e::CalTimePeakCollection> caltimepeaksH;
    if (evt.getByLabel(_trkseedLabel, caltimepeaksH)) {
      _tpeaks = caltimepeaksH.product();
    }
    else {
      _tpeaks = 0;
      printf(" >>> ERROR in CalPatRecNew::findData: CaltimepeakCollection with label=%s not found.\n",
             _trkseedLabel.data());
    }
    
// //-----------------------------------------------------------------------------
// // find list of MC hits - for debugging only
// //-----------------------------------------------------------------------------
//     art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
//     evt.getByLabel(_shLabel,"StrawHitMCPtr",mcptrHandle);
//     if (mcptrHandle.isValid()) {
//       _listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
//     }
//     else {
//       _listOfMCStrawHits = NULL;
//     }

 
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_shpcol != 0) && (_trkSeeds != 0) && (_tpeaks != 0);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  void CalPatRecNew::produce(art::Event& event ) {
    const char*               oname = "CalPatRecNew::produce";
    int                       npeaks;

    static TrkDef             dummydef;
    static HelixDefHack       dummyhdef;

    static HelixFitHackResult dummyhfit(dummyhdef);

    static StrawHitFlag       esel(StrawHitFlag::energysel), flag;

    //    ConditionsHandle<AcceleratorParams> accPar("ignored");
    //    _mbtime = accPar->deBuncherPeriod;
    //    fgTimeOffsets->updateMap(event);

                                        // event printout
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    unique_ptr<TrackSeedCollection>    outseeds(new TrackSeedCollection);
    //    unique_ptr<CalTimePeakCollection>  tpeaks(_tpeaks);
                                        // find the data
    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n", oname);
                                                            goto END;
    }

//-----------------------------------------------------------------------------
// loop over found time peaks - for us, - "eligible" calorimeter clusters
//-----------------------------------------------------------------------------
    npeaks = _tpeaks->size();

    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      const CalTimePeak* tp = &_tpeaks->at(ipeak);

//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information
//-----------------------------------------------------------------------------
      HelixDefHack helixdef(_shcol, _shpcol, _shfcol, tp->_index, _tpart, _fdir);

      TrkDef       seeddef (helixdef);

                                        // track fitting objects for this peak

      HelixFitHackResult hf_result(helixdef);

//-----------------------------------------------------------------------------
// Step 1: pattern recognition. Find initial helical approximation of a track
//-----------------------------------------------------------------------------
      int rc = _hfit.findHelix(hf_result, tp);

      if (rc) {
//-----------------------------------------------------------------------------
// pattern recognition succeeded, the seed fit starts
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// convert the result to standard helix parameters, and initialize the seed definition helix
//-----------------------------------------------------------------------------

        HepVector hpar;
        HepVector hparerr;
        _hfit.helixParams(hf_result,hpar,hparerr);

        HepSymMatrix hcov = vT_times_v(hparerr);
//----------------------------------------------------------------------------------------
// work around missing default constructor
//-----------------------------------------------------------------------------
        if (_helTraj == 0)  _helTraj = new HelixTraj(hpar,hcov);
        else               *_helTraj = HelixTraj(hpar,hcov);

        seeddef.setHelix(*_helTraj);
//-----------------------------------------------------------------------------
// P.Murat: here hits are ordered by index - WHY?
// the Kalman fitter needs them ordered in Z(straw)
//-----------------------------------------------------------------------------
        std::vector<hitIndex> goodhits;

        _index = _hfit._xyzp;

        std::sort(_index.begin(), _index.end(), [ ]( const XYZPHack& lhs,
                                                     const XYZPHack& rhs )
                  {
                    return lhs._ind < rhs._ind;
                  } );

        _nindex = _index.size();

        for (int i=0; i< _nindex; ++i){
          if (_index[i].isOutlier()) continue;
          goodhits.push_back(_index[i]._ind);
        }
        seeddef.setIndices (goodhits);
  
	//fill seed information
	TrackSeed      tmpseed;
	
	art::Ptr<CaloCluster> clusterPtr = _trkSeeds->at(ipeak)._caloCluster;

	initTrackSeed(tmpseed, seeddef, hf_result, tp, _strawhitsH, clusterPtr);

	outseeds->push_back(tmpseed);
      }

      
      
    }

//--------------------------------------------------------------------------------    
// fill diagnostic if needed
//--------------------------------------------------------------------------------
    if (_diagLevel > 0) {
      int   nseeds = outseeds->size();
      
      _hist.nseeds[0]->Fill(nseeds);
      
      double          radius(0), nhits(0), pT(0), p(0), chi2XY(0), chi2ZPhi(0);
      int             nseedsCut0(0), nhitsMin(15);
      TrackSeed      *tmpseed;
      double          mm2MeV = 3./10.;

      for (int i=0; i<nseeds; ++i){
	tmpseed  = &outseeds->at(i);

	radius   = 1./fabs(tmpseed->omega());
	nhits    = tmpseed->_selectedTrackerHits.size();

	pT       = mm2MeV*radius;
	p        = pT/std::cos( std::atan(tmpseed->tanDip()));

	chi2XY   = tmpseed->chi2XY  ();
	chi2ZPhi = tmpseed->chi2ZPhi();
	
	if (nhits >= nhitsMin) {
	  ++nseedsCut0;
	  _hist.helixFit.pT      [1] ->Fill(pT);
	  _hist.helixFit.p       [1] ->Fill(p);
	  _hist.helixFit.radius  [1] ->Fill(radius);
	  _hist.helixFit.chi2XY  [1] ->Fill(chi2XY);
	  _hist.helixFit.chi2ZPhi[1] ->Fill(chi2ZPhi);
	}
  	

	_hist.helixFit.radius  [0] ->Fill(radius);
	_hist.helixFit.nhits       ->Fill(nhits);
	_hist.helixFit.pT      [0] ->Fill(pT);
	_hist.helixFit.p       [0] ->Fill(p);
	_hist.helixFit.chi2XY  [0] ->Fill(chi2XY);
	_hist.helixFit.chi2ZPhi[0] ->Fill(chi2ZPhi);
	
	_hist.helixFit.nhitsvspT ->Fill(nhits, pT);
	_hist.helixFit.nhitsvsp  ->Fill(nhits, p );
	
      }

      _hist.nseeds[1]->Fill(nseedsCut0);

      
    }

//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record
//-----------------------------------------------------------------------------
  END:;
    event.put(std::move(outseeds));

  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRecNew::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

//--------------------------------------------------------------------------------
  void CalPatRecNew::initTrackSeed(TrackSeed                             &TrkSeed, 
				   TrkDef                                &SeedDef  , 
				   HelixFitHackResult                    &HfResult ,
				   const CalTimePeak                     *TPeak    , 
				   art::Handle<mu2e::StrawHitCollection> &StrawhitsH,
				   art::Ptr<CaloCluster>                  ClusterPtr){

    //set helix parameters
    TrkSeed._fullTrkSeed._d0     = SeedDef.helix().d0();
    TrkSeed._fullTrkSeed._phi0   = SeedDef.helix().phi0();
    TrkSeed._fullTrkSeed._omega  = SeedDef.helix().omega();
    TrkSeed._fullTrkSeed._z0     = SeedDef.helix().z0();
    TrkSeed._fullTrkSeed._tanDip = SeedDef.helix().tanDip();
    
    for(int i=0;i<5;i++) {
      for(int j=0;j<5;j++) {
	TrkSeed._fullTrkSeed._covMtrx[i][j] = SeedDef.helixCovMatr()[i][j];
      }
    }
    
    int             shIndices = SeedDef.strawHitIndices().size();
    const mu2e::hitIndex *hIndex;

    for (int i=0; i<shIndices; ++i){
      hIndex = &SeedDef.strawHitIndices().at(i);
      TrkSeed._fullTrkSeed._selectedTrackerHitsIdx.push_back( mu2e::hitIndex( hIndex->_index, hIndex->_ambig) );
      TrkSeed._selectedTrackerHits.push_back( StrawHitPtr (StrawhitsH, hIndex->_index) );
    }
    
    const mu2e::CaloCluster *cluster = TPeak->Cluster();
    
    TrkSeed._t0          = cluster->time(); 
    TrkSeed._errt0       = cluster->cog3Vector().z();
    TrkSeed._caloCluster = ClusterPtr;
    
    TrkSeed._chi2XY      = HfResult._sxy.chi2DofCircle() / HfResult._sxy.qn();
    TrkSeed._chi2ZPhi    = HfResult._srphi.chi2DofLine() / HfResult._srphi.qn();
  }


}

using mu2e::CalPatRecNew;
DEFINE_ART_MODULE(CalPatRecNew);
