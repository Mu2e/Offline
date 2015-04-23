// $Id: CalPatRec_module.cc,v 1.15 2014/09/19 20:49:45 murat Exp $
// $Author: murat $ 
// $Date: 2014/09/19 20:49:45 $
//
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalPatRec_module.hh"
#include "CalPatRec/inc/Ref.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

using namespace std; 
using namespace boost::accumulators;

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
// module constructor
//-----------------------------------------------------------------------------
  CalPatRec::CalPatRec(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>        ("diagLevel")),
    _debug       (pset.get<int>        ("debugLevel")),
    _printfreq   (pset.get<int>        ("printFrequency")),
    _addhits     (pset.get<bool>       ("addhits")),
    _shLabel     (pset.get<string>("StrawHitCollectionLabel"        )),
    _shpLabel    (pset.get<string>("StrawHitPositionCollectionLabel")),
    _shfLabel    (pset.get<string>("StrawHitFlagCollectionLabel"    )),    
    _ccmLabel    (pset.get<string>("caloClusterModuleLabel"         )),

    _dtspecpar   (pset.get<string>("DeltaTSpectrumParams","nobackgroundnomarkovgoff")),
    _tsel        (pset.get<vector<string> >("TimeSelectionBits")),
    _hsel        (pset.get<vector<string> >("HelixFitSelectionBits")),
    _addsel      (pset.get<vector<string> >("AddHitSelectionBits",vector<string>{} )),
    _ksel        (pset.get<vector<string> >("KalmanFitSelectionBits")),
    _bkgsel      (pset.get<vector<string> >("BackgroundSelectionBits")),
    _addbkg      (pset.get<vector<string> >("AddHitBackgroundBits",vector<string>{})),
    _maxedep     (pset.get<double>("MaxStrawEDep",0.005)),
    //    fUseDoublets (pset.get<int>("useDoublets")),
    _mindt       (pset.get<double>("DtMin",-70.0)),
    _maxdt       (pset.get<double>("DtMax", 20.0)),
    _maxdtmiss   (pset.get<double>("DtMaxMiss",55.0)),
    _fbf         (pset.get<double>("PhiEdgeBuffer",1.1)),
    _maxnpeak        (pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits        (pset.get<int>   ("MinNHits" ,20)),
    _tmin            (pset.get<double>("tmin")),
    _tmax            (pset.get<double>("tmax")),
    _tbin            (pset.get<double>("tbin"             ,20.0)),
    _minClusterEnergy(pset.get<double>("minClusterEnergy" )),
    _minClusterSize  (pset.get<int>("minClusterSize" )),
    _ymin            (pset.get<double>("ymin"             ,4)),
    _1dthresh        (pset.get<double>("OneDPeakThreshold",4.0)),
    _pitchAngle      (pset.get<double>("_pitchAngle"      ,0.67)),
    _maxseeddoca     (pset.get<double>("MaxSeedDoca"      ,10.0)),
    _maxhelixdoca    (pset.get<double>("MaxHelixDoca"     ,40.0)),
    _maxadddoca      (pset.get<double>("MaxAddDoca"       )),
    _maxaddchi       (pset.get<double>("MaxAddChi"        ,4.0)),
    _tpart           ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir            ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _hfit            (pset.get<fhicl::ParameterSet>("HelixFitHack",fhicl::ParameterSet())),
    _seedfit         (pset.get<fhicl::ParameterSet>("SeedFitHack",fhicl::ParameterSet())),
    _kfit            (pset.get<fhicl::ParameterSet>("KalFitHack",fhicl::ParameterSet()))
		      //    _payloadSaver(pset)
		      //    , _kfitmc      (pset.get<fhicl::ParameterSet>("KalFitMC",fhicl::ParameterSet()))
  {
    //    fStopwatch = new TStopwatch();

					// tag the data product instance by the direction 
					// and particle type found by this fitter

    _iname      = _fdir.name() + _tpart.name();
    _iname_seed = _iname + "seed";
    produces<KalRepCollection>      (_iname);
    produces<KalRepPtrCollection>   (_iname);

    produces<KalRepPtrCollection>   (_iname_seed);

    produces<StrawHitFlagCollection>(_iname);
    produces<CalTimePeakCollection> (_iname);

    //    produces<KalRepPayloadCollection>();

					// set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);

    fNminMChits = 25;
    
    fQualityTrack = 0;
    
    fHackData = new THackData("HackData","Hack Data");
    gROOT->GetRootFolder()->Add(fHackData);
//-----------------------------------------------------------------------------
// provide for interactive disanostics
//-----------------------------------------------------------------------------
    _ref = new Ref("CalPatRecRef","Ref to CalPatRec",&_hfit,&_seedfit,&_kfit);

    TFolder* f_mu2e;

    f_mu2e = (TFolder*) gROOT->GetRootFolder()->FindObject("Mu2e"); 
    if (f_mu2e == NULL) f_mu2e = gROOT->GetRootFolder()->AddFolder("Mu2e","Mu2e Folder");

    if (f_mu2e) {
      _folder = f_mu2e->AddFolder("CalPatRec","CalPatRec Folder");
      _folder->Add(_ref);
    }
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalPatRec::~CalPatRec() {
    delete _ref;
    //    delete fStopwatch;
  }

//-----------------------------------------------------------------------------
  void CalPatRec::beginJob(){
    // create diagnostics ntuple if requested
    if(_diag > 0)createDiagnostics();
    // create a histogram of throughput: this is 
    // a basic diagnostic that should ALWAYS be on

    art::ServiceHandle<art::TFileService> tfs;

    _hist._cutflow = tfs->make<TH1F>("cutflow","Cutflow",10,-0.5,9.5);
    _hist._cutflow->GetXaxis()->SetBinLabel(1,"nhits(CE) >= 25");
    _hist._cutflow->GetXaxis()->SetBinLabel(2,"Time Peak");	 
    _hist._cutflow->GetXaxis()->SetBinLabel(3,"Helix Fit");	 
    _hist._cutflow->GetXaxis()->SetBinLabel(4,"Seed Fit");	 
    _hist._cutflow->GetXaxis()->SetBinLabel(5,"Kalman Fit");     
    _hist._cutflow->GetXaxis()->SetBinLabel(6,"Cut set C & p>100MeV/c");     

    _hist._Tpeaks    = tfs->make<TH1F>("hTpeaks",
				       "Time peaks per event",100,0,100);
    _hist._NfitIter  = tfs->make<TH1F>("hNfitIter",
				       "Numebr of fit iteration on kalman::fiIteration",
				       100,0,100);
     
    //   _hist._Tfit[0]   = tfs->make<TH1F>("hTfit0",
    // 				 "Time per event spent on kalman maketrack",200,0,10);
    //     _hist._Tfit[1]   = tfs->make<TH1F>("hTfit1",
    // 				 "Time per event spent on Kalman addHits",200,0,10);

    //     _hist._Ttot      = tfs->make<TH1F>("hTtot",
    // 				 "Total time per event",200,0,10);
    
    _hist._dfdzmode  = tfs->make<TH1F>("hdfdzmode",
				       "index of the loop on dfdz modes",20,0,20);
    _hist._radius    = tfs->make<TH1F>("hradius",
				       "radius of the theretical helix",
				       1000,0.,700.);
    _hist._phi0      = tfs->make<TH1F>("hphi0",
				       "#phi_{0} of the theretical helix",
				       1000,-10.,10.);
    _hist._tanlambda = tfs->make<TH1F>("htanlambda",
				       "#tan(#lambda) of the theretical helix",
				       600,0.,6.);
    
    _hist._dphidz[0] = tfs->make<TH1F>("hdphidz0","dfdz from calculateDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
				       10000,-0.01,0.01);

    _hist._dphidz[1] = tfs->make<TH1F>("hdphidz1","dfdz from findDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
				       1000,-0.001,0.001);
    _hist._dphidz[2] = tfs->make<TH1F>("hdphidz2","dfdz from doLinearFitDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
				       1000,-0.001,0.001);
    _hist._kdphidz[0] = tfs->make<TH1F>("hkdphidz0","dfdz from calculateDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					1000,-0.001,0.001);
    
    _hist._kdphidz[1] = tfs->make<TH1F>("hkdphidz1","dfdz from findDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					1000,-0.001,0.001);
    _hist._kdphidz[2] = tfs->make<TH1F>("hkdphidz2","dfdz from doLinearFitDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					1000,-0.001,0.001);
    
    _hist._0mode     = tfs->make<TH1F>("h0mode",
				       "point rescued with dfdz recalculation",
				       20,0,20);
    _hist._dist      = tfs->make<TH1F>("hdist",
				       "distance between strahit points in the timepeak and the prediction; dist [mm]",
				       1000, 0.,1e3);
    _hist._dz        = tfs->make<TH1F>("hdz",
				       "distance along z between the two strahits used for the pattern-reco; dz [mm]",
				       200, 600.,600.);
    _hist._Npoints   = tfs->make<TH1F>("hNpoints",
				       "Number of points belong to the predictedtrajectory; N-points [#]",
				       100, 0., 100.);
    _hist._chi2      = tfs->make<TH1F>("hchi2",
				       "#chi^{2} distribution for track candidate; #chi^{2}",
				       50000, 0., 5000.);

    _hist._distvsdz  = tfs->make<TH2F>("hdistvsdz",
				       "Distance from prediction vs z-distance form the seed",
				       1400, -3500., 3500.,
				       500, 0, 500);

    _hist._kdfdzmode  = tfs->make<TH1F>("hkdfdzmode",
					"index of the loop on dfdz modes when Kalman filter converged",20,0,20);
    _hist._kradius[0]    = tfs->make<TH1F>("hkradius0",
					   "radius of the theretical helix when Kalman filter converged",
					   2000,-100.,100.);
    _hist._kradius[1]    = tfs->make<TH1F>("hkradius1",
					   "radius of the theretical helix when Kalman filter converged + cut set ''C'' and p>100 MeV/c",
					   2000,-100.,100.);
    _hist._kphi0      = tfs->make<TH1F>("hkphi0",
					"#phi_{0} of the theretical helix when Kalman filter converged",
					1000,-10.,10.);
    _hist._ktanlambda = tfs->make<TH1F>("hktanlambda",
					"#tan(#lambda) of the theretical helix when Kalman filter converged",
					600,0.,6.);
    _hist._kdfdz[0]      = tfs->make<TH1F>("hkdfdz0","dfdz from pattern recognition when Kalman filter converged; (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					   1000,-0.001,0.001);
    _hist._kdfdz[1]      = tfs->make<TH1F>("hkdfdz1","dfdz from pattern recognition when Kalman filter converged + cut set ''C'' and p>100 MeV/c; (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					   1000,-0.001,0.001);

    _hist._seeddfdz[0]      = tfs->make<TH1F>("hseeddfdz0",
					      "dfdz from seedFit when Kalman filter converged; (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					      1000,-0.001,0.001);
    _hist._seeddfdz[1]      = tfs->make<TH1F>("hseeddfdz1",
					      "dfdz from seedFit when Kalman filter converged + cut set ''C'' and p>100 MeV/c; (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
					      1000,-0.001,0.001);
    _hist._k0mode  = tfs->make<TH1F>("hk0mode", "point rescued with dfdz recalculation", 20,0,20);
    _hist._kdist      = tfs->make<TH1F>("hkdist",
					"distance between strahit points in the timepeak and the prediction; dist [mm]",
					1000, 0.,1e3);
    _hist._kdz        = tfs->make<TH1F>("hkdz",
					"distance along z between the two strahits used for the pattern-reco; dz [mm]",
					1200, -600.,600.);
    _hist._kNpoints   = tfs->make<TH1F>("hkNpoints",
					"Number of points belong to the predictedtrajectory; N-points [#]",
					100, 0., 100.);

    _hist._kchi2      = tfs->make<TH1F>("hkchi2",
					"#chi^{2} distribution for track candidates with converged KF; #chi^{2}",
					50000, 0., 5000.);

    _hist._PhiResid[0]= tfs->make<TH1F>("hPhiResid0",
					"#phi residual 0; #phi_{straw} - #phi_{fit} [rad]",
					200, -1, 1.);

    _hist._PhiResid[1]= tfs->make<TH1F>("hPhiResid1",
					"#phi residual 1; #phi_{straw} - #phi_{fit} [rad]",
					200, -1, 1.);

    _hist._kdistvsdz[0]  = tfs->make<TH2F>("hkdistvsdz0",
					   "Distance from prediction versus z-distance form the seed in case also the kalman fit converged; z-distance from the seed [mm]; Distance from prediction [mm]",
					   1400, -3500., 3500.,
					   500, 0, 500);
    _hist._kdistvsdz[1]  = tfs->make<TH2F>("hkdistvsdz1",
					   "Distance from prediction versus z-distance form the seed in case also the kalman fit converged + cut set ''C'' and p>100 MeV/c; z-distance from the seed [mm]; Distance from prediction [mm]",
					   1400, -3500., 3500.,
					   500, 0, 500);

    _hist._drw[0]    = tfs->make<TH1F>("hdrw_0","r(fitw) - r(track)[0]", 200,-100.,100.);
    _hist._drw[1]    = tfs->make<TH1F>("hdrw_1","r(fitw) - r(track)[1]", 200,-100.,100.);
 
    _hist._seeddr[0]    = tfs->make<TH1F>("hseeddr_0","r(seedfit) - r(track)[0]", 200,-100.,100.);
    _hist._seeddr[1]    = tfs->make<TH1F>("hseeddr_1","r(seedfit) - r(track)[1]", 200,-100.,100.);

    _hist._chi2w[0]  = tfs->make<TH1F>("hchi2w_0","chi2(fitw)[0]", 500,0.,10.);
    _hist._chi2w[1]  = tfs->make<TH1F>("hchi2w_1","chi2(fitw)[1]", 500,0.,10.);
    
    _hist._chi2zphi[0]  = tfs->make<TH1F>("hchi2zphi_0","chi2(fitzphi)[0]", 500,0.,10);
    _hist._chi2zphi[1]  = tfs->make<TH1F>("hchi2zphi_1","chi2(fitzphi)[1]", 500,0.,10);

    _hist._seeddoca[0] = tfs->make<TH1F>("hseeddoca_0","doca seedfit active hits; doca [mm]", 1000, -20., 20);
    _hist._seeddoca[1] = tfs->make<TH1F>("hseeddoca_1","doca seedfit non active hits; doca [mm]", 1000, -20., 20);
    _hist._seeddoca[2] = tfs->make<TH1F>("hseeddoca_2","doca seedfit all hits; doca [mm]", 1000, -20., 20);

    _hist._doca[0]     = tfs->make<TH1F>("hdoca_0","doca helixfit active hits; doca [mm]"      , 3000, -60., 60);
    _hist._doca[1]     = tfs->make<TH1F>("hdoca_1","doca helixfit non active hits; doca [mm]"  , 3000, -60., 60);
    _hist._doca[2]     = tfs->make<TH1F>("hdoca_2","doca helixfit 1 active hits; doca [mm]"    , 3000, -60., 60);
    _hist._doca[3]     = tfs->make<TH1F>("hdoca_3","doca helixfit 1 non active hits; doca [mm]", 3000, -60., 60);
    //-----------------------------------------------------------------------------
    // doublet histograms
    //-----------------------------------------------------------------------------
    _hist._ndoublets[0]      = tfs->make<TH1F>("nd_0","N(doublets)"            ,   50,    0,  50);
    _hist._ndoublets[1]      = tfs->make<TH1F>("nd_1","N(OS doublets)"         ,   50,    0,  50);
    _hist._ndoublets[2]      = tfs->make<TH1F>("nd_2","N(SS doublets)"         ,   50,    0,  50);

    _hist._ndstraws[0]      = tfs->make<TH1F>("nds_0","N(straws/doublet)"      ,   10,    0,  10);
    _hist._ndstraws[1]      = tfs->make<TH1F>("nds_1","N(straws/doublet) OS "  ,   10,    0,  10);
    _hist._ndstraws[2]      = tfs->make<TH1F>("nds_2","N(straws/doublet) SS "  ,   10,    0,  10);

    _hist._tslope[0]        = tfs->make<TH1F>("tslope_0","Track Slope [0]"    ,  500, -2.5, 2.5);
    _hist._tslope[1]        = tfs->make<TH1F>("tslope_1","Track Slope [1]"    ,  500, -2.5, 2.5);
    _hist._tslope[2]        = tfs->make<TH1F>("tslope_2","Track Slope [2]"    ,  500, -2.5, 2.5);

    _hist._dsbest[0]        = tfs->make<TH1F>("dsbest_0","St(s)-Sl(t) best[0]", 1000, -2.5, 2.5);
    _hist._dsbest[1]        = tfs->make<TH1F>("dsbest_1","St(s)-Sl(t) best[1]", 1000, -2.5, 2.5);
    _hist._dsbest[2]        = tfs->make<TH1F>("dsbest_2","St(s)-Sl(t) best[2]", 1000, -2.5, 2.5);

    _hist._dsnext[0]        = tfs->make<TH1F>("dsnext_0","St(s)-Sl(t) next[0]", 1000, -2.5, 2.5);
    _hist._dsnext[1]        = tfs->make<TH1F>("dsnext_1","St(s)-Sl(t) next[1]", 1000, -2.5, 2.5);
    _hist._dsnext[2]        = tfs->make<TH1F>("dsnext_2","St(s)-Sl(t) next[2]", 1000, -2.5, 2.5);

    _hist._NpointsSeed [0]  = tfs->make<TH1F>("hnpseed_0","# of points de-activatedby seedfit; N points de-activated [#]", 50, 0., 50);
    _hist._NpointsSeed [1]  = tfs->make<TH1F>("hnpseed_1","# of points not found in the pattern-reco and rescued ; N points rescued [#]", 50, 0., 50);
 
    _hist._kaldoca    [0]   = tfs->make<TH1F>("hkaldoca_0","doca kalfit active hits; doca [mm]", 1000, -20., 20);
    _hist._kaldoca    [1]   = tfs->make<TH1F>("hkaldoca_1","doca kalfit non active hits; doca [mm]", 1000, -20., 20);

    _hist._NpointsRescued [0]  = tfs->make<TH1F>("hnprescued_0","fraction of points rescued [0]; np_{rescued}/np_{tot}", 100, 0., 1);
    _hist._NpointsRescued [1]  = tfs->make<TH1F>("hnprescued_1","fraction of points rescued [1]; np_{rescued}/np_{tot}", 100, 0., 1);

    _hist._ntracks           = tfs->make<TH1F>("ntracks","N(reconstructed tracks)", 10, 0, 10);
    _eventid = 0;
  }

//-----------------------------------------------------------------------------
  void CalPatRec::beginRun(art::Run& ){
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
  }

//-----------------------------------------------------------------------------
// find the input data objects 
//-----------------------------------------------------------------------------
  bool CalPatRec::findData(const art::Event& evt) {

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalPatRec::findData: StrawHitCollection with label=%s not found.\n",
	     _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (evt.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CalPatRec::findData: StrawHitPositionCollection with label=%s not found.\n",
	     _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalPatRec::findData: StrawHitFlagCollection with label=%s not found.\n",
	     _shfLabel.data());
    }

    art::Handle<CaloClusterCollection> ccH;
    if (evt.getByLabel(_ccmLabel, ccH)) {
      _ccCollection = ccH.product();
    }
    else {
      _ccCollection = 0;
      printf(" >>> ERROR in CalPatRec::findData: CaloClusterCollection with label=%s not found.\n",
	     _ccmLabel.data());
    }
//-----------------------------------------------------------------------------
// find list of MC hits - for debugging only
//-----------------------------------------------------------------------------
    art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_shLabel,"StrawHitMCPtr",mcptrHandle);
    if (mcptrHandle.isValid()) {
      _listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      _listOfMCStrawHits = NULL;
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_shpcol != 0) && (_ccCollection != 0);
  }


//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  void CalPatRec::produce(art::Event& event ) {
    const char* oname = "CalPatRec::produce";
    char        message[200];
    bool                      findhelix (false), findseed (false), findkal (false);
    int                       nhits;
    int                       npeaks;
    int                       gen_index, nhits_from_gen(0), sim_id;
 
    ::KalRep*                 krep;
					// dummy objects
    static TrkDef             dummydef;
    static HelixDefHack       dummyhdef;

    static HelixFitHackResult dummyhfit(dummyhdef);
    static KalFitResult       dummykfit(dummydef);

    static StrawHitFlag       esel(StrawHitFlag::energysel), flag;

    _ntracks = 0;
					// reset the fit iteration counter
    _kfit.setNIter(0);
//     t1 = fStopwatch->RealTime();
//     fStopwatch->Continue();
					// event printout
    _eventid = event.event();
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) cout<<"CalPatRec: event="<<_iev<<endl;
    //    _hist._cutflow->Fill(0.0);
					// create output

    _tpeaks = new CalTimePeakCollection;
    unique_ptr<KalRepCollection>       tracks(new KalRepCollection      );
    unique_ptr<KalRepPtrCollection>    trackPtrs(new KalRepPtrCollection);
    unique_ptr<KalRepPtrCollection>    seedTrkPtrs(new KalRepPtrCollection);
    unique_ptr<CalTimePeakCollection>  tpeaks(_tpeaks);
    
    _flags = new StrawHitFlagCollection();
    unique_ptr<StrawHitFlagCollection> flags (_flags);

    art::ProductID kalRepsID(getProductID<KalRepCollection>(event,_iname));

					// find the data
    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n",oname); 
                                                            goto END;
    }

                                // count the number of straw hits generated by the CE
    nhits = _shcol->size();
     
    for (int i=0; i<nhits; i++) {
      mu2e::PtrStepPointMCVector const& mcptr( _listOfMCStrawHits->at(i ) );
      const mu2e::StepPointMC* Step = mcptr[0].operator ->();
    
      gen_index = -1;
      if (Step) {
	art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle(); 
	
	if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
	else                         gen_index = -1;
	
	sim_id        = simptr->id().asInt();
      }
      if (gen_index >0 && sim_id == 1) ++nhits_from_gen;
    }
    
    if (nhits_from_gen >= fNminMChits)  _hist._cutflow->Fill(0.0);
    
    _kfit.setStepPointMCVectorCollection(_listOfMCStrawHits);
    _seedfit.setStepPointMCVectorCollection(_listOfMCStrawHits);
//-----------------------------------------------------------------------------
// all needed pieces of data have been found, 
// tighten the energy cut and copy flags, clear 
//-----------------------------------------------------------------------------
    nhits = _shcol->size();
    for (int i=0; i<nhits; i++) {
      flag = _shfcol->at(i);
      if (_shcol->at(i).energyDep() > _maxedep && flag.hasAllProperties(esel)) {
	flag.clear(esel);
      }
      _flags->push_back(flag);
    }
//-----------------------------------------------------------------------------
// find the time peaks in the time spectrum of selected hits.  
//-----------------------------------------------------------------------------
    findTimePeaks(_tpeaks);
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    _hist._Tpeaks->Fill(_tpeaks->size());

    if (_tpeaks->size()>0 && (nhits_from_gen >= fNminMChits)) _hist._cutflow->Fill(1.0);
//-----------------------------------------------------------------------------
// loop over found time peaks - for us, - "eligible" calorimeter clusters 
//-----------------------------------------------------------------------------
    npeaks = _tpeaks->size();
 
    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      CalTimePeak* tp = &_tpeaks->at(ipeak);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
      if (_debug > 0) {
	const StrawHit*     hit;
	int nhits = tp->NHits();
	printf(" peak # ipeak = %2i; nhits = %5i\n",ipeak,nhits);
	if (_debug > 1) {
	  for (int ih=0; ih<nhits; ih++) {
	    hitIndex ind = tp->_index[ih];
	    hit = &_shcol->at(ind._index);
	    printf("index = %5i time=%10.3f energy = %10.3f\n",
		   hit->strawIndex().asInt(),hit->time(),hit->energyDep());
	  }
	}
      }
//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information 
//-----------------------------------------------------------------------------
      HelixDefHack helixdef(_shcol,_shpcol,_flags,tp->_index,_tpart,_fdir);

					// set some identifiers
      helixdef.setEventId(_eventid);
      helixdef.setTrackId(ipeak);
					// copy this for the other fits

      TrkDef             seeddef(helixdef);
      TrkDef             kaldef (helixdef);

					// track fitting objects for this peak

      HelixFitHackResult hf_result(helixdef);
      KalFitResult       sf_result (seeddef);
      KalFitResult       kf_result  (kaldef);
//-----------------------------------------------------------------------------
// pattern recognition step - find initial approximation for Kalman fitter
//-----------------------------------------------------------------------------
      int rc = _hfit.findHelix(hf_result,tp);

      if (_debug > 0) {
	printf("[CalPatRec::produce] helixFit status = %i\n", rc);
	_hfit.printInfo(hf_result);
      }
      
      if (rc) {
//-----------------------------------------------------------------------------
// pattern recognition succeeded
//-----------------------------------------------------------------------------
	findhelix = true;
//-----------------------------------------------------------------------------
// convert the result to standard helix parameters, and initialize the seed definition helix
//-----------------------------------------------------------------------------
	double dz, dist;
	for (int i=0; i< fHackData->goodPoints(); ++i){
	  dz   = fHackData->fDz[i];
	  dist = fHackData->fDist[i]; 
	  _hist._distvsdz->Fill(dz, dist);
	}

	HepVector hpar;
	HepVector hparerr;
	_hfit.helixParams(hf_result,hpar,hparerr);
	
	HepSymMatrix hcov = vT_times_v(hparerr);
	HelixTraj helTraj(hpar,hcov);
	seeddef.setHelix(HelixTraj(hpar,hcov));
//-----------------------------------------------------------------------------
// P.Murat: here hits are ordered by index - WHY?  
// - Kalman fit needs them ordered in Z(straw)
//-----------------------------------------------------------------------------
	std::vector<hitIndex> goodhits;

	_index = _hfit.fxyzp;

	std::sort(_index.begin(), _index.end(), [ ]( const XYZPHack& lhs,
						     const XYZPHack& rhs )
		  {
		    return lhs._ind < rhs._ind;
		  } );

	_nindex = _index.size();

	Hep3Vector             shPos;
	int                    loc;
	const mu2e::StrawHit*  hit;
	const mu2e::Straw*     straw;

	if (_debug > 0) {
	  printf("[CalPatRec::printGoodHits]   Index   Straw    status       X         Y         Z         Zw\n");
	}

	for (int i=0; i< _nindex; ++i){
	  if (_index[i].isOutlier()) continue;

	  if (_debug > 0) {
	    shPos = _index[i]._pos;
	    loc   = _index[i]._ind;
	    hit   = &_shcol->at(loc);
	    straw = &_tracker->getStraw(hit->strawIndex());
	    printf("[CalPatRec::printGoodHits]  %6i  %6i    active   %8.3f  %8.3f  %9.3f %9.3f\n", 
		   loc, hit->strawIndex().asInt(), shPos.x(), shPos.y(), shPos.z(),
		   straw->getMidPoint().z()
		   );
	  }
	  goodhits.push_back(_index[i]._ind);
	}
//-----------------------------------------------------------------------------
//  Trajectory info
//-----------------------------------------------------------------------------
	Hep3Vector tdir;
	HepPoint   tpos;
	double     doca;
	helTraj.getInfo(0.0,tpos,tdir);

	for (int i=0; i< _nindex; ++i){
	  StrawHit const*   sh    = _index[i]._strawhit;
	  Straw const&      straw = _tracker->getStraw(sh->strawIndex());
	  CLHEP::Hep3Vector wpos  = straw.getMidPoint();
	  CLHEP::Hep3Vector wdir  = straw.getDirection();

	  // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	  HepPoint      wpt(wpos.x(),wpos.y(),wpos.z());
	  TrkLineTraj  wtraj(wpt,wdir,-20,20);
	  // estimate flightlength along track.  This assumes a constant BField!!!
	  double fltlen = (wpos.z()-tpos.z())/tdir.z();
	  TrkPoca   wpoca(helTraj,fltlen,wtraj,0.0);
	    
	  doca      = wpoca.doca();

	  if (_index[i].isOutlier()) _hist._doca[1]->Fill(doca);
	  else                       _hist._doca[0]->Fill(doca);
	}

	if (_debug > 0) {
	  printf("[CalPatRec::seeddef] goodhits = %lu over nIndex = %i\n", goodhits.size(), _nindex); 
	}

	seeddef.setIndices (goodhits);

	if (_debug > 0) {
	  int shIndeces = seeddef.strawHitIndices().size();
	  int nSh       = seeddef.strawHitCollection()->size();
	  printf("[CalPatRec::seedfit] START N-straws = %i N-indices = %i \n", nSh, shIndeces);
	  printf("[CalPatRec::seedfit]");
	  for(int i=0; i<5; ++i) {
	    printf(" hpar[%i] =  %5.5f", i, hpar[i]);
	  }
	  printf("\n");
	}
//-----------------------------------------------------------------------------
// seed fit - fit through the wires of found hits, not using the drift times
//-----------------------------------------------------------------------------
	_seedfit.makeTrack(sf_result, tp);
//--------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following diagnnostic
//--------------------------------------------------------------------------------
	if (_debug > 0) {
	  sprintf(message,
		  "CalPatRec::produce seedfit::makeTrack : fit_success = %i\n",
		  sf_result._fit.success());
	  _seedfit.printHits(sf_result,message);
	}

//--------------------------------------------------------------------------------
// 2015-03-23 G. Pezzu: fill info about the doca
//--------------------------------------------------------------------------------
	if (sf_result._fit.success()) {
	  findseed = true;
//-----------------------------------------------------------------------------
// use helix parameters by the seed fit to initialize the full Kalman fit 
//-----------------------------------------------------------------------------
	  double           locflt;
	  const HelixTraj* shelix;
	  shelix = dynamic_cast<const HelixTraj*>(sf_result._krep->localTrajectory(sf_result._krep->flt0(),locflt));
	  kaldef.setHelix(*shelix);
//----------------------------------------------------------------------
//2015-02-07 G. Pezzu added new selection using seedFit results	
//----------------------------------------------------------------------
	  goodhits.clear();
	  const mu2e::TrkStrawHit* hit;
	  int                      hit_index;
	  const TrkHotList*        hot_list = sf_result._krep->hotList();

					//  Trajectory info
	  Hep3Vector               tdir;
	  HepPoint                 tpos;
	  double                   doca, rdrift, fltlen;
	  bool                     found(false), active;
	  int                      npointsrescued(0), npointsdeactivated(0);
	  int                      banner_11_printed(0);
	  
	  sf_result._krep->traj().getInfo(0.0,tpos,tdir);

	  for (int i=0; i< _nindex; ++i) {
	    StrawHit const*   sh    = _index[i]._strawhit;
	    Straw const&      straw = _tracker->getStraw(sh->strawIndex());
	    CLHEP::Hep3Vector wpos  = straw.getMidPoint();
	    CLHEP::Hep3Vector wdir  = straw.getDirection();
	    rdrift = -9990;
	    found  = false;
	    active = false;

	    HepPoint      wpt(wpos.x(),wpos.y(),wpos.z());
	    TrkLineTraj   wtraj(wpt,wdir,-20,20);
					// estimate flightlength along track. This assumes a constant BField!!!
	    fltlen = (wpos.z()-tpos.z())/tdir.z();
	    TrkPoca   wpoca(sf_result._krep->traj(),fltlen,wtraj,0.0);
	    
	    doca      = wpoca.doca();
	    hit_index = _index[i]._ind;

	    for (auto it=hot_list->begin(); it<hot_list->end(); it++) {
	      hit    = (const mu2e::TrkStrawHit*) &(*it);
	      rdrift = hit->driftRadius();
	      int shIndex = int(hit->index());
	      if (hit_index == shIndex) {
		found = true;
		break;
	      }
	    }
	    
	    if ( std::fabs(doca) < _maxadddoca) {
	      active = true;
	      goodhits.push_back(hit_index);
	    }
	    
	    if (_debug > 0) {
	      if (banner_11_printed == 0) { 
		banner_11_printed = 1;
		printf("[CalPatRec::produce] -------------------------------------\n");
		printf("[CalPatRec::produce]  ih  A   Sind      Rdrift        doca\n");
		printf("[CalPatRec::produce] -------------------------------------\n");
	      }

	      printf("[CalPatRec::produce]  %2i  %1i  %5i  %10.3f  %10.3f \n", 
		     i, active? 1:0, straw.index().asInt(), rdrift, doca );
	    }
	    
	    _hist._seeddoca[2]->Fill(doca);

	      if (active) _hist._seeddoca[0]->Fill(doca);
	      else        _hist._seeddoca[1]->Fill(doca);

	    if (!found && active){
	      ++npointsrescued;
	    }
	  }
//-----------------------------------------------------------------------------
// fill histograms: after the SEED fit, the DOCA is the distance to the wire
//-----------------------------------------------------------------------------
	  for (int i=0; i< _nindex; ++i){
	    StrawHit const*   sh    = _index[i]._strawhit;
	    Straw const&      straw = _tracker->getStraw(sh->strawIndex());
	    CLHEP::Hep3Vector wpos  = straw.getMidPoint();
	    CLHEP::Hep3Vector wdir  = straw.getDirection();
	    rdrift = -9990;
	    found  = false;
	    active = false;
	    
	    HepPoint      wpt  (wpos.x(),wpos.y(),wpos.z());
	    TrkLineTraj   wtraj(wpt,wdir,-20,20);
	    fltlen = (wpos.z()-tpos.z())/tdir.z();
	    TrkPoca       wpoca(helTraj,fltlen,wtraj,0.0);
	    
	    doca   = wpoca.doca();
	    
	    if (_index[i].isOutlier()) _hist._doca[3]->Fill(doca);
	    else                       _hist._doca[2]->Fill(doca);
	  }

	  _hist._NpointsSeed[1]->Fill(npointsrescued);

	  for (auto it=hot_list->begin(); it<hot_list->end(); it++) {
	    hit = (const mu2e::TrkStrawHit*) &(*it);
	    if (!hit->isActive()) 
	      ++npointsdeactivated;
	  }
	  _hist._NpointsSeed[0]->Fill(npointsdeactivated);
//-----------------------------------------------------------------------------
// at this point the full kalman fit starts
//-----------------------------------------------------------------------------
	  kaldef.setIndices(goodhits);
	  if (_debug > 0) printf("CalPatRec::produce] calling _kfit.makeTrack\n");
	  //	  _kfit.setDecisionMode(0);
	  _kfit.makeTrack(kf_result,tp);

	  if (_debug > 0) {
	    printf("[CalPatRec::produce] kalfit status = %i\n", kf_result._fit.success());
	    _kfit.printHits(kf_result,"CalPatRec::produce kalfit_001");
	  }

	  if (kf_result._fit.success()) {
	    findkal = true;
	    
	    if (_addhits) {
//-----------------------------------------------------------------------------
// this is the default. First, add back the hits on this track
// if successfull, try to add missing hits, at this point external errors were 
// set to zero
// assume this is the last iteration
//-----------------------------------------------------------------------------
	      int final(1);
	      int last_iteration = _kfit.maxIteration();

	      _kfit.unweedHits(kf_result,_maxaddchi);
	      if (_debug > 0) _kfit.printHits(kf_result,"CalPatRec::produce after unweedHits");

	      std::vector<hitIndex> misshits;
	      findMissingHits(kf_result,misshits);

					// force drift sign decision to be made
	      //	      _kfit.setDecisionMode(1);
//-----------------------------------------------------------------------------
// if new hits have been added, add then and refit the track. 
// Otherwise - just refit the track one last time
// in both cases 
//-----------------------------------------------------------------------------
	      if (misshits.size() > 0) {
		_kfit.addHits(kf_result,_shcol,misshits, _maxaddchi, tp);
	      }
	      else {
		_kfit.fitIteration(kf_result,-1,tp,final);
	      }

	      if (_debug > 0) _kfit.printHits(kf_result,"CalPatRec::produce after addHits");
//-----------------------------------------------------------------------------
// and weed hits again to insure that addHits doesn't add junk
//-----------------------------------------------------------------------------
	      _kfit.weedHits(kf_result,last_iteration,final);
	    }
//-----------------------------------------------------------------------------
// now evaluate the T0 ad its error using the straw hits
//-----------------------------------------------------------------------------
	    _kfit.updateT0(kf_result);

	    if (_debug > 0) {
	      _kfit.printHits(kf_result,"CalPatRec::produce : final, after weedHits");
	    }

//-----------------------------------------------------------------------------
// done, fill debug histograms
//-----------------------------------------------------------------------------
	   const TrkHotList* hot_l = kf_result._krep->hotList();
	    
	    for (int i=0; i< _nindex; ++i){
	      StrawHit const*     sh = _index[i]._strawhit;
	      Straw const&     straw = _tracker->getStraw(sh->strawIndex());
	      CLHEP::Hep3Vector hpos = straw.getMidPoint();
	      CLHEP::Hep3Vector hdir = straw.getDirection();
	      found = false;
	      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	      HepPoint      spt(hpos.x(),hpos.y(),hpos.z());
	      TrkLineTraj htraj(spt,hdir,-20,20);
	      // estimate flightlength along track.  This assumes a constant BField!!!
	      double fltlen = (hpos.z()-tpos.z())/tdir.z();
	      TrkPoca hitpoca(kf_result._krep->traj(),fltlen,htraj,0.0);
	    
	      doca = hitpoca.doca();
	      for(TrkHotList::hot_iterator it=hot_l->begin(); it<hot_l->end(); it++) {
		hit   = (const mu2e::TrkStrawHit*) &(*it);
		if (!hit->isActive()) continue;
		hit_index = hit->index();
		if (int(_index[i]._ind) == hit_index){
		  found = true;
		  break;
		}
	      }
	      if (found) _hist._kaldoca[0]->Fill(doca);
	      else       _hist._kaldoca[1]->Fill(doca);
	    }
	  }
	}
      }

      if (kf_result._fit.success()) {
	_ntracks += 1;
//-----------------------------------------------------------------------------
// fit succeeded
// flag hits used in this track.  This should use the track id, FIXME!!! (in the BaBar code)
//-----------------------------------------------------------------------------
	if (ipeak<16) {
	  for (size_t ihit=0; ihit<kf_result._hits.size(); ++ihit){
	    const TrkStrawHit* tsh = kf_result._hits[ihit];
	    if (tsh->isActive()) {
	      _flags->at(tsh->index()).merge(StrawHitFlag::trackBit(ipeak));
	      _flags->at(tsh->index()).merge(StrawHitFlag::calosel);
	    }
	  }
	}
//-----------------------------------------------------------------------------
//  fill fit diagnostics histograms if requested
//-----------------------------------------------------------------------------
	fillFitDiag(event,ipeak,hf_result,sf_result,kf_result);
//-----------------------------------------------------------------------------
// save successful kalman fits in the event
//-----------------------------------------------------------------------------
	krep = kf_result.stealTrack();
	tracks->push_back(krep);
        int index = tracks->size()-1;
        trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	tp->SetCprIndex(tracks->size());
      } 
      else {
//-----------------------------------------------------------------------------
// fit failed, just delete the track
//-----------------------------------------------------------------------------
	kf_result.deleteTrack();
      }
//-----------------------------------------------------------------------------
// cleanup the seed fit
//-----------------------------------------------------------------------------
      sf_result.deleteTrack();
    }
//-----------------------------------------------------------------------------
// diagnostics in the end
//-----------------------------------------------------------------------------
    if (findhelix && (nhits_from_gen >= fNminMChits)) _hist._cutflow->Fill(2.0);
    if (findseed  && (nhits_from_gen >= fNminMChits)) _hist._cutflow->Fill(3.0);
    if (findkal   && (nhits_from_gen >= fNminMChits)) {
      _hist._cutflow->Fill(4.0);
      if (fQualityTrack > 0) {
	_hist._cutflow->Fill(5.0);
      }
    }

    if (_debug > 0) {
      if ((nhits_from_gen >= fNminMChits) && _tpeaks->size() > 0){
	if (_tpeaks->at(0)._tmin > 400.){
	  if (!findhelix){
	    printf("[CalPatRec::produce] LOOK AT: more than 25 MC hits and findHelix not converged! event = %i\n", _iev);
	  }if (findhelix && !findseed){
	    printf("[CalPatRec::produce] LOOK AT: findhelix converged and findseed not! event = %i\n", _iev);
	  }
	  if (findseed && !findkal){
	    printf("[CalPatRec::produce] LOOK AT: findseed converged and findkal not! event = %i\n", _iev);
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// fill event-level histograms
//-----------------------------------------------------------------------------
    if(_diag > 0) {
      _hist._NfitIter->Fill(_kfit.nIter());
      _hist._ntracks->Fill(_ntracks);
    }
//-----------------------------------------------------------------------------
// put tracks into the event
//-----------------------------------------------------------------------------
  END:;
    //    art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
    //    _payloadSaver.put(*tracks, tracksID, event);
    event.put(std::move(tracks),   _iname);
    event.put(std::move(trackPtrs),_iname);
    event.put(std::move(flags ),   _iname);
    event.put(std::move(tpeaks),   _iname);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  void CalPatRec::findTimePeaks(CalTimePeakCollection* TimePeakColl) {

    int                 ncl, nsh;
    double              time, dt, tof, zstraw, cl_time, stime; 
    double              xcl, ycl, zcl/*, dz_cl*/;
    const CaloCluster*  cl;
    const StrawHit*     hit;
    const Straw*        straw;
    Hep3Vector          gpos, tpos;
//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
    nsh = _shcol->size();
    ncl = _ccCollection->size();

    for (int ic=0; ic<ncl; ic++) {
      cl      = &_ccCollection->at(ic);

      if ( (cl->energyDep() > _minClusterEnergy) && 
	   (int(cl->size()) > _minClusterSize) ) {
	cl_time = cl->time();
//-----------------------------------------------------------------------------
// convert cluster coordinates defined in the disk frame to the detector
// coordinate system
//-----------------------------------------------------------------------------
	gpos = _calorimeter->fromSectionFrameFF(cl->sectionId(),cl->cog3Vector());
	tpos = _calorimeter->toTrackerFrame(gpos);

	xcl     = tpos.x();
	ycl     = tpos.y();
	zcl     = tpos.z();

	//	dz_cl   = zcl; // -_tracker->z0();
					// create time peak
	CalTimePeak tpeak(cl,xcl,ycl,zcl);

	tpeak._shcol  = _shcol;
	tpeak._shfcol = _shfcol;
	tpeak._tmin   = cl_time+_mindt;
	tpeak._tmax   = cl_time+_maxdt;
//-----------------------------------------------------------------------------
// record hits in time with each peak, and accept them if they have a minimum # of hits
//-----------------------------------------------------------------------------
	stime = 0;
	mu2e::StrawHitFlag flag;

	for(int istr=0; istr<nsh;++istr) {
	  flag = _flags->at(istr);
	  int hit_has_all_properties = flag.hasAllProperties(_hsel);
	  int bgr_hit                = flag.hasAnyProperty(_bkgsel);

	  //	  if (_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)) {
	  if (hit_has_all_properties && !bgr_hit) {
	    hit    = &_shcol->at(istr);
	    time   = hit->time();
	    straw  = &_tracker->getStraw(hit->strawIndex());
	    zstraw = straw->getMidPoint().z();
//-----------------------------------------------------------------------------
// estimate time-of-flight and calculate residual between the predicted and the hit times
//-----------------------------------------------------------------------------
	    tof = (zcl-zstraw)/sin(_pitchAngle)/CLHEP::c_light;
	    dt  = cl_time-(time+tof);

	    if ((dt < _maxdt) && (dt >= _mindt)) {
	      tpeak._index.push_back(istr);
	      stime += time;
	    }
	  }
	}

	tpeak._tpeak = stime/(tpeak.NHits()+1.e-12);
	if (tpeak.NHits() > _minnhits) TimePeakColl->push_back(tpeak);
      }
    }
  }

//-----------------------------------------------------------------------------
// 2014-04-08 P.M.: I don't think this function is called any more
//-----------------------------------------------------------------------------
  void CalPatRec::createTimePeak(CalTimePeakCollection* TimePeakColl) {
// find the median time
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc;
    unsigned nstrs = _shcol->size();
    double   time;

    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)) {
	time = _shcol->at(istr).time();
	tacc(time);
      }
    }

    int np = boost::accumulators::extract::count(tacc);  
    if(np >= _minnhits){
      double mtime  = median(tacc);
      // create a time peak from the full subset of selected hits
      CalTimePeak tpeak(0, 0., 0., 0.);
      for(unsigned istr=0; istr<nstrs;++istr){
	if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)){
	  tpeak._index.push_back(istr);
	}
      }
      tpeak._tpeak = mtime;
      TimePeakColl->push_back(tpeak);
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::filterOutliers(TrkDef&                    mytrk  , 
				 Trajectory const&          traj   , 
				 double                     maxdoca, 
				 std::vector<TrkHitFilter>& thfvec ) {
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;
    traj.getInfo(0.0,tpos,tdir);

    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    const StrawHitCollection* hits = mytrk.strawHitCollection();
    const std::vector<hitIndex>& indices = mytrk.strawHitIndices();
    std::vector<hitIndex> goodhits;

    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = hits->at(indices[ihit]._index);
      Straw const& straw = _tracker->getStraw(sh.strawIndex());
      CLHEP::Hep3Vector hpos = straw.getMidPoint();
      CLHEP::Hep3Vector hdir = straw.getDirection();
      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
      TrkLineTraj htraj(spt,hdir,-20,20);
      // estimate flightlength along track.  This assumes a constant BField!!!
      double fltlen = (hpos.z()-tpos.z())/tdir.z();
      TrkPoca hitpoca(traj,fltlen,htraj,0.0);
      // flag hits with small residuals
      if(fabs(hitpoca.doca()) < maxdoca){
	goodhits.push_back(indices[ihit]);
      }
      // optional diagnostics
      if(_diag > 0){
	// summarize the MC truth for this strawhit
	TrkHitFilter thfilter;
	HepPoint tpos =  traj.position(hitpoca.flt1());
	thfilter._pos = CLHEP::Hep3Vector(tpos.x(),tpos.y(),tpos.z());
	thfilter._doca = hitpoca.doca();
	thfvec.push_back(thfilter);
      }
    }
    // update track
    mytrk.setIndices(goodhits);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::findMissingHits(KalFitResult& kalfit,std::vector<hitIndex>& misshits) {
					//  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;

    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    if (_debug > 0) {
      printf("[CalPatRec::findMissingHits]      shId    sec     panel       doca   \n");
    }

    for(unsigned istr=0; istr<nstrs;++istr){
//----------------------------------------------------------------------//
// 2015-02-11 gianipez  and P. Murat changed the selection bit for searching missed   //
// hits                                                                 //
//----------------------------------------------------------------------//
//      if(_flags->at(istr).hasAllProperties(_addsel)&& !_flags->at(istr).hasAnyProperty(_addbkg)){

	StrawHit const& sh = _shcol->at(istr);
	if(fabs(_shcol->at(istr).time()-kalfit._krep->t0()._t0) < _maxdtmiss) {
	  // make sure we haven't already used this hit
	  std::vector<TrkStrawHit*>::iterator ifnd = find_if(kalfit._hits.begin(),kalfit._hits.end(),FindTrkStrawHit(sh));
	  if(ifnd == kalfit._hits.end()){
	    // good in-time hit.  Compute DOCA of the wire to the trajectory
	    Straw const& straw = _tracker->getStraw(sh.strawIndex());
	    CLHEP::Hep3Vector hpos = straw.getMidPoint();
	    CLHEP::Hep3Vector hdir = straw.getDirection();
	    // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	    HepPoint spt(hpos.x(),hpos.y(),hpos.z());
	    TrkLineTraj htraj(spt,hdir,-20,20);
	    // estimate flightlength along track.  This assumes a constant BField!!!
	    double fltlen = (hpos.z()-tpos.z())/tdir.z();
	    TrkPoca hitpoca(kalfit._krep->pieceTraj(),fltlen,htraj,0.0);

	    if (_debug>0){
	      printf("[CalPatRec::findMissingHits] %8i  %6i  %8i  %10.3f \n",
		     straw.index().asInt(),
		     straw.id().getDevice(),
		     straw.id().getSector(),
		     hitpoca.doca());
	    }


	    // flag hits with small residuals
	    if(fabs(hitpoca.doca()) < _maxadddoca){
	      misshits.push_back(istr);
	    }
	  }
	}
      }
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;

  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::fillStrawDiag() {

  }

  void CalPatRec::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;

    char rsname[100];
    char csname[100];
    char tsname[100];
    char lsname[100];
    char tdsname[100];
    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(csname,100,"convtspectrum%i",_iev);
    snprintf(tsname,100,"tighttspectrum%i",_iev);
    snprintf(lsname,100,"loosetspectrum%i",_iev);
    snprintf(tdsname,100,"tightnodeltatspectrum%i",_iev);
    _hist.ttsp = tfs->make<TH1F>(tsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    _hist.ttsp->SetLineColor(kCyan);
    _hist.ltsp = tfs->make<TH1F>(lsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    _hist.ltsp->SetLineColor(kGreen);
    _hist.rtsp = tfs->make<TH1F>(rsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    _hist.rtsp->SetLineColor(kBlue);
    _hist.ctsp = tfs->make<TH1F>(csname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    _hist.ctsp->SetLineColor(kRed);
    _hist.tdtsp = tfs->make<TH1F>(tdsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    _hist.tdtsp->SetLineColor(kOrange);

    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = _shcol->at(istr).time();
      bool conversion(false);
      // summarize the MC truth for this strawhit
//       if(_kfitmc.mcData()._mcsteps != 0) {
// 	const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr); 
// 	conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
//       }
      // fill plots
      _hist.rtsp->Fill(time);
      if(_flags->at(istr).hasAllProperties(_tsel)){
	_hist.ttsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)){
	_hist.tdtsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_ksel) && !_flags->at(istr).hasAnyProperty(_bkgsel)){
	_hist.ltsp->Fill(time);
      }
      if(conversion){
	_hist.ctsp->Fill(time);
      }
    }

  }

//-----------------------------------------------------------------------------
// this routine is called only if a track has been found
//-----------------------------------------------------------------------------
  void CalPatRec::fillFitDiag (art::Event&         event    ,
			       int                 ipeak    , 
			       HelixFitHackResult& hf_result,
			       KalFitResult&       sf_result, 
			       KalFitResult&       kf_result) {
    Hep3Vector tdir, seedMom;
    HepPoint   tpos;
    bool       found;
    //    int        nhits_peak;
    double     seedPt, seedPz, seedTanL, seedRadius, seeddfdz;    

    sf_result._krep->traj().getInfo(0.0,tpos,tdir);
    
    fQualityTrack = 0;
    
    KalRep*    krep = kf_result._krep;

    _ipeak = ipeak;
//-----------------------------------------------------------------------------
// fit status 
//-----------------------------------------------------------------------------
    _helixfail = hf_result._fit.failure();
    _seedfail  = sf_result._fit.failure();
    _kalfail   = kf_result._fit.failure();
//-----------------------------------------------------------------------------
// histograms 
//-----------------------------------------------------------------------------
    const TrkStrawHit*     hit;
    int                    hit_index;
    double                 fltlen;

    const TrkHotList* hotl = krep->hotList();
	    
    for (int i=0; i<_nindex; ++i){
      const StrawHit*   sh    = _index[i]._strawhit;
      Straw const&      straw = _tracker->getStraw(sh->strawIndex());
      CLHEP::Hep3Vector hpos  = straw.getMidPoint();
      CLHEP::Hep3Vector hdir  = straw.getDirection();

      fltlen        = (hpos.z()-tpos.z())/tdir.z();
//-----------------------------------------------------------------------------
// histogram distance of closest approach for hits on the track
//-----------------------------------------------------------------------------
      HepPoint      spt (hpos.x(),hpos.y(),hpos.z());

      TrkLineTraj   htraj(spt,hdir,-20,20);
      TrkPoca       hitpoca(krep->traj(),fltlen,htraj,0.0);
      double doca = hitpoca.doca();
      
      found = false;
      for(TrkHotList::hot_iterator it=hotl->begin();it<hotl->end(); it++) {
	hit   = (const mu2e::TrkStrawHit*) &(*it);
	if (!hit->isActive()) continue;
	hit_index = hit->index();
	if (int(_index[i]._ind) == hit_index) {
	  found = true;
	  break;
	}
      }
      if (found) _hist._kaldoca[0]->Fill(doca);
      else       _hist._kaldoca[1]->Fill(doca);
    }
//----------------------------------------------------------------------
// 2014-11-02 gianipez added some diagnostic
//----------------------------------------------------------------------
    _hist._kdfdzmode ->Fill(fHackData->TheoImode());
    _hist._kphi0     ->Fill(fHackData->TheoPhi0());
    _hist._ktanlambda->Fill(fHackData->TheoTanL());

    Hep3Vector mom = krep->momentum(0);
    double pt      = sqrt(mom.x()*mom.x() + mom.y()*mom.y());
    double pz      = mom.z();
    double tanL    = pt/pz;
    double radius  = fabs(1./krep->helix(0).omega());//mom.mag()*10./3.;//convert MeV to mm
    double kdfdz   = tanL/radius;
	
    _hist._kdfdz  [0]->Fill(fHackData->dfdz() - kdfdz);
    _hist._kradius[0]->Fill(fHackData->TheoRadius() - radius);
    _hist._k0mode    ->Fill(fHackData->mode0Points());
//------------------------------------------------------------------------------------------
//take info from the seed fit for filling diagnostic histograms
//------------------------------------------------------------------------------------------
    seedMom      = sf_result._krep->momentum(0);
    seedPt       = sqrt(seedMom.x()*seedMom.x() + seedMom.y()*seedMom.y());
    seedPz       = seedMom.z();
    seedTanL     = seedPt/seedPz;
    seedRadius   = fabs(1./sf_result._krep->helix(0).omega());//mom.mag()*10./3.;//convert MeV to mm
    seeddfdz     = seedTanL/seedRadius;
	  
    _hist._seeddfdz[0]->Fill(seeddfdz - kdfdz);
    _hist._seeddr  [0]->Fill(seedRadius - radius);
    _hist._drw     [0]->Fill(fHackData->fData[14]-radius);
    _hist._chi2w   [0]->Fill(fHackData->fData[15]);
    _hist._chi2zphi[0]->Fill(fHackData->fData[13]);
    
    _hist._kdz     ->Fill(fHackData->shDz());
    _hist._kNpoints->Fill(fHackData->goodPoints());
    _hist._kchi2   ->Fill(fHackData->chi2());
    
    _hist._NpointsRescued[0]->Fill(fHackData->rescuedPoints()/double(fHackData->goodPoints()));

    _hist._dphidz[0]->Fill(fHackData->fData[17]- kdfdz);
    _hist._dphidz[1]->Fill(fHackData->fData[18]- kdfdz);
    _hist._dphidz[2]->Fill(fHackData->fData[19]- kdfdz);
    
    double dz, dist, dphi;
    for (int i=0; i< fHackData->goodPoints(); ++i){
      dz   = fHackData->fDz[i];
      dist = fHackData->fDist[i]; 
      dphi = fHackData->fResid[i]; 
      _hist._kdistvsdz[0]->Fill(dz, dist);
      _hist._PhiResid [0]->Fill(dphi);
    }
//-----------------------------------------------------------------------------
// histograms for doublets
//-----------------------------------------------------------------------------
    double             dsl;
    Doublet*           d;
    const TrkStrawHit* dhit [2];
    int                layer[2], nd, nd_tot(0), nd_os(0), nd_ss(0), ns;
    
    std::vector<Doublet>* list_of_doublets = &kf_result._listOfDoublets;
    nd = list_of_doublets->size();

    for (int i=0; i<nd; i++) {
      d  = &list_of_doublets->at(i);
      ns = d->fNstrawHits;
					
      _hist._ndstraws[0]->Fill(ns);
      if (d->fOs == 0) _hist._ndstraws[1]->Fill(ns);
      else             _hist._ndstraws[2]->Fill(ns);

      if (ns > 1) { 
	nd_tot += 1;
	if (d->fOs == 0) nd_os += 1;
	else             nd_ss += 1;
      }

      if (ns != 2)                                          continue;

      for (int i=0; i<2; i++) {
	dhit  [i] = d->fHit[i];
	layer [i] = dhit[i]->straw().id().getLayer();
      }
      if (layer[0] == layer[1])                             continue;
//-----------------------------------------------------------------------------
// to simplify, consider only 2-hit doublets, determine best and next combinations
// some "doublets" may have 1 hit, some - more than two
//-----------------------------------------------------------------------------
      double dsl_best(1.e6), dsl_next(1.e6);
      //      int    ibest(-1)/*, inext(-1)*/;
      
      for (int is=0; is<4; is++) {
	dsl = d->fTrkDxDz-d->fDxDz[is];
	if (fabs(dsl) < fabs(dsl_best)) { 
	  dsl_next = dsl_best;
	  //	  inext    = ibest;
	  dsl_best = dsl;
	  // ibest    = is;
	}
	else if (fabs(dsl) < fabs(dsl_next)) {
	  dsl_next = dsl;
	  //	  inext    = is;
	}
      }
					// all doublets
      _hist._tslope[0]->Fill(d->fTrkDxDz);
      _hist._dsbest[0]->Fill(dsl_best);
      _hist._dsnext[0]->Fill(dsl_next);
      
      if (d->fOs == 0) {
					// OS doublets
	_hist._tslope[1]->Fill(d->fTrkDxDz);
	_hist._dsbest[1]->Fill(dsl_best);
	_hist._dsnext[1]->Fill(dsl_next);
      }
      else {
					// SS doublets
	_hist._tslope[2]->Fill(d->fTrkDxDz);
	_hist._dsbest[2]->Fill(dsl_best);
	_hist._dsnext[2]->Fill(dsl_next);

	if (fabs(d->fTrkDxDz) < 0.05) {
	  printf("Event : %10i id : %2i dxdz(trk) : %10.5f ntrk: %2i\n",
		 event.event(),i,d->fTrkDxDz,_ntracks);
	}
      }
    }

    _hist._ndoublets[0]->Fill(nd_tot);
    _hist._ndoublets[1]->Fill(nd_os);
    _hist._ndoublets[2]->Fill(nd_ss);
//----------------------------------------------------------------------
// 2015-01-17 G. Pezzu: tracks passing cut set "C", except for the timing cut!
//----------------------------------------------------------------------
    double  fMinFitCons      = 2.e-3;
    double  fMinNActive      = 25;
    double  fMaxT0Err        = 0.9;  		// ns
    double  fMaxFitMomErr    = 0.25;  		// MeV
    double  fMinTanDip       = tan(M_PI/6.);	// 0.5773
    double  fMaxTanDip       = 1.0;  
    double  fMinD1           = -80.;		//  mm
    double  fMaxD1           = 105.;
    double  minMom           = 100.;                // MeV/c

    BbrVectorErr      momerr = krep->momentumErr(0);

    CLHEP::Hep3Vector momdir =  krep->momentum(0).unit();
    HepVector momvec(3);
    for (int i=0; i<3; i++) momvec[i] = momdir[i];
    
    double fitmom_err = sqrt(momerr.covMatrix().similarity(momvec));
    
    if ( (krep->chisqConsistency().consistency() > fMinFitCons  ) &&
	 (krep->nActive()                        > fMinNActive  ) &&
	 (krep->t0().t0Err()                     < fMaxT0Err    ) &&
	 (fitmom_err                             < fMaxFitMomErr) &&
	 (krep->helix(0).tanDip()                > fMinTanDip   ) &&
	 (krep->helix(0).tanDip()                < fMaxTanDip   ) &&
	 (krep->helix(0).d0()                    < fMaxD1       ) &&
	 (krep->helix(0).d0()                    > fMinD1       ) &&
	 (mom.mag()                              > minMom       )    ) {

      fQualityTrack = 1;

      _hist._kradius [1]->Fill(fHackData->TheoRadius() - radius);
      _hist._kdfdz   [1]->Fill(fHackData->dfdz() - kdfdz);
      _hist._seeddfdz[1]->Fill(seeddfdz - kdfdz);
      _hist._seeddr  [1]->Fill(seedRadius - radius);
      _hist._drw     [1]->Fill(fHackData->fData[14]-radius);
      _hist._chi2w   [1]->Fill(fHackData->fData[15]);
      _hist._chi2zphi[1]->Fill(fHackData->fData[13]);
      
      _hist._kdphidz [0]->Fill(fHackData->fData[17]- kdfdz);
      _hist._kdphidz [1]->Fill(fHackData->fData[18]- kdfdz);
      _hist._kdphidz [2]->Fill(fHackData->fData[19]- kdfdz);
      
      for (int i=0; i< fHackData->goodPoints(); ++i){
	dz   = fHackData->fDz[i];
	dist = fHackData->fDist[i]; 
	dphi = fHackData->fResid[i]; 
	_hist._PhiResid[1]->Fill(dphi);
	_hist._kdistvsdz[1]->Fill(dz, dist);
      }
    }
  }
}

using mu2e::CalPatRec;
DEFINE_ART_MODULE(CalPatRec);
