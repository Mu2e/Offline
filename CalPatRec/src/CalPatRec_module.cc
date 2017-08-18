///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalPatRec_module.hh"
#include "CalPatRec/inc/Ref.hh"
#include "CalPatRec/inc/AlgorithmIDCollection.hh"

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
#include "CalPatRec/inc/KalFitResult.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"

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
  struct straw_zcomp : public binary_function<StrawHitIndex,StrawHitIndex,bool> {
    bool operator()(StrawHitIndex const& h1, StrawHitIndex const& h2) {

      mu2e::GeomHandle<mu2e::TTracker> handle;
      const TTracker* t = handle.get();
      const Straw* s1 = &t->getStraw(StrawIndex(h1));
      const Straw* s2 = &t->getStraw(StrawIndex(h2));

      return s1->getMidPoint().z() < s2->getMidPoint().z();
    }
  }; // a semicolumn here is required

//-----------------------------------------------------------------------------
// module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
//-----------------------------------------------------------------------------
  CalPatRec::CalPatRec(fhicl::ParameterSet const& pset) :
    _diagLevel   (pset.get<int>        ("diagLevel")),
    _debugLevel  (pset.get<int>        ("debugLevel")),
    _printfreq   (pset.get<int>        ("printFrequency")),
    _addhits     (pset.get<bool>       ("addhits")),
    _shLabel     (pset.get<string>("StrawHitCollectionLabel"        )),
    _shDigiLabel (pset.get<string>("StrawDigiCollectionLabel"       )),
    _shpLabel    (pset.get<string>("StrawHitPositionCollectionLabel")),
    _shfLabel    (pset.get<string>("StrawHitFlagCollectionLabel"    )),
    _ccmLabel    (pset.get<string>("caloClusterModuleLabel"         )),
    _crmLabel    (pset.get<string>("caloReadoutModuleLabel"         )),
    _chmccpLabel (pset.get<string>("calorimeterHitMCCrystalPtr"     )),

    //    _dtspecpar   (pset.get<string>("DeltaTSpectrumParams","nobackgroundnomarkovgoff")),
    _tsel        (pset.get<vector<string> >("TimeSelectionBits")),
    _hsel        (pset.get<vector<string> >("HelixFitSelectionBits")),
    _addsel      (pset.get<vector<string> >("AddHitSelectionBits",vector<string>{} )),
    _ksel        (pset.get<vector<string> >("KalmanFitSelectionBits")),
    _bkgsel      (pset.get<vector<string> >("BackgroundSelectionBits")),
    _addbkg      (pset.get<vector<string> >("AddHitBackgroundBits",vector<string>{})),
    _maxedep     (pset.get<double>("MaxStrawEDep",0.005)),
    _mindt       (pset.get<double>("DtMin")),
    _maxdt       (pset.get<double>("DtMax")),
    _maxdtmiss   (pset.get<double>("MaxDtMiss")),
    _final       (pset.get<int>   ("final")),
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
    _maxadddoca      (pset.get<double>("MaxAddDoca"      )),
    _maxaddchi       (pset.get<double>("MaxAddChi"        ,4.0)),
    _tpart           ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir            ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _hfit            (pset.get<fhicl::ParameterSet>("HelixFitHack",fhicl::ParameterSet())),
    _seedfit         (pset.get<fhicl::ParameterSet>("SeedFitHack",fhicl::ParameterSet())),
    _kfit            (pset.get<fhicl::ParameterSet>("KalFitHack",fhicl::ParameterSet())),
    _sfresult(0),
    _kfresult(0)
  {
    //    fStopwatch = new TStopwatch();

                                        // tag the data product instance by the direction
                                        // and particle type found by this fitter

    produces<KalRepCollection>      ();
    produces<KalRepPtrCollection>   ();
    produces<AlgorithmIDCollection> ();

    produces<StrawHitFlagCollection>();
    produces<CalTimePeakCollection> ();

                                        // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);

    _minNMCHits = 25;

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

    fNCaloEnergyCut   = 0;
    fNCaloSizeCut     = 0;
    fNHitsTimePeakCut = 0;
    fNTimeWindow      = 0;

    fgTimeOffsets     = new SimParticleTimeOffset(pset.get<fhicl::ParameterSet>("TimeOffsets"));

    _helTraj = 0;
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalPatRec::~CalPatRec() {
    delete _ref;
    if (_helTraj) delete _helTraj;
    //    delete fStopwatch;
  }

//-----------------------------------------------------------------------------
  void CalPatRec::beginJob(){

    if(_diagLevel > 0) bookHistograms();

    _eventid = 0;
  }

//-----------------------------------------------------------------------------
  void CalPatRec::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
					// calibrations

    mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    _trackerCalib = tcal.operator ->();

    _kfit.setTracker(_tracker);
    _kfit.setTrackerCalib(_trackerCalib);

    _seedfit.setTracker(_tracker);
    _seedfit.setTrackerCalib(_trackerCalib);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;

    art::TFileDirectory hf_dir = tfs->mkdir("HelixFit");
    art::TFileDirectory sf_dir = tfs->mkdir("SeedFit");
    art::TFileDirectory kf_dir = tfs->mkdir("KalFit");

    _hist.helixFit.nhits = hf_dir.make<TH1F>( "nhits", "N(Helix Hits)", 100, 0, 100);

    _hist._cutflow[0] = tfs->make<TH1F>("cutflow0","Cutflow",10,-0.5,9.5);
    _hist._cutflow[0]->GetXaxis()->SetBinLabel(1,"nhits(CE) >= 25");
    _hist._cutflow[0]->GetXaxis()->SetBinLabel(2,"Time Peak");
    _hist._cutflow[0]->GetXaxis()->SetBinLabel(3,"Helix Fit");
    _hist._cutflow[0]->GetXaxis()->SetBinLabel(4,"Seed Fit");
    _hist._cutflow[0]->GetXaxis()->SetBinLabel(5,"Kalman Fit");
    _hist._cutflow[0]->GetXaxis()->SetBinLabel(6,"Cut set C & p>100MeV/c");

    _hist._cutflow[1] = tfs->make<TH1F>("cutflow1","Cutflow",12,-0.5,11.5);
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(1, "Nhits(CE) >= 25 && p_{MC} > 100 MeV/c && SH-time >500 ns");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(2, Form("E_{calo} > %4.3f", _minClusterEnergy));
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(3, Form("Size calo-cluster > %i", _minClusterSize  ));
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(4, "Time window ");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(9, "SH time sel");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(5, "SH radius sel");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(6, "SH energy sel");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(7, "SH Delta ray sel");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(8, "SH isolation sel");
    _hist._cutflow[1]->GetXaxis()->SetBinLabel(10, "Helix search");

    _hist._dt[0]     = tfs->make<TH1F>("hdt0",
                                       "#Delta t = t_{calo} - #left( t_{straw} - <t_{drift}> + <tof>#right) [ns]",400,-200,200);
    _hist._dt[1]     = tfs->make<TH1F>("hdt1",
                                       "#Delta t = t_{calo} - #left( t_{straw} - <t_{drift}> + <tof>#right) [ns]",400,-200,200);

    _hist._Tpeaks    = tfs->make<TH1F>("hTpeaks",
                                       "Time peaks per event",100,0,100);
    _hist._NfitIter  = tfs->make<TH1F>("hNfitIter",
                                       "Number of fit iteration on kalman::fiIteration",
                                       100,0,100);

    _hist._dphi0 [0]     = tfs->make<TH1F>("hdphi0_0",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad]",
                                           1000,-10.,10.);

    _hist._dphi0 [1]     = tfs->make<TH1F>("hdphi0_1",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad]",
                                           1000,-1.,1.);

    _hist._dphi0 [2]     = tfs->make<TH1F>("hdphi0_2",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad]",
                                           1000,-1.,1.);

    _hist._dphi0 [3]     = tfs->make<TH1F>("hdphi0_3",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad]",
                                           1000,-1.,1.);

    _hist._dphidz[0] = tfs->make<TH1F>("hdphidz0","dfdz from calculateDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                       500,-0.01,0.01);

    _hist._dphidz[1] = tfs->make<TH1F>("hdphidz1","dfdz from findDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                       1000,-0.001,0.001);

    _hist._dphidz[2] = tfs->make<TH1F>("hdphidz2","dfdz from doLinearFitDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                       1000,-0.001,0.001);

    _hist._dphidz[3] = tfs->make<TH1F>("hdphidz3","dfdz from doLinearFitDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                       1000,-0.001,0.001);

    _hist._kdphidz[0] = tfs->make<TH1F>("hkdphidz0","dfdz from calculateDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                        500,-0.01,0.01);

    _hist._kdphidz[1] = tfs->make<TH1F>("hkdphidz1",
                                        "dfdz from findDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                        1000,-0.001,0.001);

    _hist._kdphidz[2] = tfs->make<TH1F>("hkdphidz2",
                                        "dfdz from doLinearFitDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                        1000,-0.001,0.001);

    _hist._kdphidz[3] = tfs->make<TH1F>("hkdphidz3",
                                        "dfdz from doLinearFitDfDz(); (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                        1000,-0.001,0.001);

    _hist._distvsdz  = tfs->make<TH2F>("hdistvsdz",
                                       "Distance from prediction vs z-distance form the seed",
                                       1400, -3500., 3500.,
                                       500, 0, 500);

    _hist._kradius[0]    = tfs->make<TH1F>("hkradius0",
                                           "radius of the theretical helix when Kalman filter converged",
                                           2000,-100.,100.);
    _hist._kradius[1]    = tfs->make<TH1F>("hkradius1",
                                           "radius of the theretical helix when Kalman filter converged + cut set ''C'' and p>100 MeV/c",
                                           2000,-100.,100.);
    _hist._kdphi0 [0]    = tfs->make<TH1F>("hkdphi0_0",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad] when Kalman filter converged",
                                           1000,-10.,10.);

    _hist._kdphi0 [1]    = tfs->make<TH1F>("hkdphi0_1",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad] when Kalman filter converged",
                                           1000,-10.,10.);

    _hist._kdphi0 [2]    = tfs->make<TH1F>("hkdphi0_2",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad] when Kalman filter converged",
                                           1000,-10.,10.);

    _hist._kdphi0 [3]    = tfs->make<TH1F>("hkdphi0_3",
                                           "#phi_{0}_{hel} - #phi_{0}_{trk} [rad] when Kalman filter converged",
                                           1000,-10.,10.);

    _hist._seeddfdz[0]      = tfs->make<TH1F>("hseeddfdz0",
                                              "dfdz from seedFit when Kalman filter converged; (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                              1000,-0.001,0.001);
    _hist._seeddfdz[1]      = tfs->make<TH1F>("hseeddfdz1",
                                              "dfdz from seedFit when Kalman filter converged + cut set ''C'' and p>100 MeV/c; (d#phi/dz)_{hel} - (d#phi/dz)_{trk} [rad/mm]",
                                              1000,-0.001,0.001);
    _hist._kdz        = tfs->make<TH1F>("hkdz",
                                        "distance along z between the two strahits used for the pattern-reco; dz [mm]",
                                        1200, -600.,600.);
    _hist._kNpoints   = tfs->make<TH1F>("hkNpoints",
                                        "Number of points belong to the predictedtrajectory; N-points [#]",
                                        100, 0., 100.);

    _hist._PhiResid[0]= tfs->make<TH1F>("hPhiResid0",
                                        "#phi residual 0; #phi_{straw} - #phi_{pat-rec} [rad]",
                                        200, -1, 1.);

    _hist._PhiResid[1]= tfs->make<TH1F>("hPhiResid1",
                                        "#phi residual when track passes cut set C and p > 100 MeV/c; #phi_{straw} - #phi_{pat-rec} [rad]",
                                        200, -1, 1.);
    _hist._PhiResid[2]= tfs->make<TH1F>("hPhiResid2",
                                        "#phi residual 1; #phi_{straw} - #phi_{kal-fit} [rad]",
                                        200, -1, 1.);

    _hist._kdistvsdz[0]  = tfs->make<TH2F>("hkdistvsdz0",
                                           "Dist from prediction[mm] versus z-dist form the seed[mm], kalman fit converged",
                                           1400, -3500., 3500.,
                                           500, 0, 500);
    _hist._kdistvsdz[1]  = tfs->make<TH2F>("hkdistvsdz1",
                                           "Dist from prediction[mm] versus z-dist form the seed[mm], kalman fit converged + setC + and p>100",
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

    _hist._ntracks           = tfs->make<TH1F>("ntracks","N(reconstructed tracks)", 10, 0,   10);
    _hist._nhits[0]          = tfs->make<TH1F>("nhits_0","N(straw hits)[0]"       , 500,0,  500);
    _hist._nhits[1]          = tfs->make<TH1F>("nhits_1","N(straw hits)[1]"       , 500,0,10000);

    _hist.rtsp               = tfs->make<TH1F>("rtsp" ,"rtsp"      , 200,0,2000);
    _hist.ttsp               = tfs->make<TH1F>("ttsp" ,"ttsp"      , 200,0,2000);
    _hist.tdtsp              = tfs->make<TH1F>("tdtsp","tdtsp"     , 200,0,2000);
    _hist.ltsp               = tfs->make<TH1F>("ltsp" ,"ltsp"      , 200,0,2000);
    _hist.ctsp               = tfs->make<TH1F>("ctsp" ,"ctsp"      , 200,0,2000);
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
    //    evt.getByLabel(_shDigiLabel,"StrawHitMCPtr",mcptrHandle);
    evt.getByLabel(_shDigiLabel,mcptrHandle);
    if (mcptrHandle.isValid()) {
      _listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      _listOfMCStrawHits = NULL;
    }

//------------------------------------------------------------------------------------------
// Utility to match  cloHits with MCtruth, simParticles and StepPoints
//------------------------------------------------------------------------------------------

    //Get calorimeter readout hits (2 readout / crystal as of today)
    art::Handle<CaloHitCollection> caloHitsHandle;
    if (evt.getByLabel(_crmLabel, caloHitsHandle)){
      _chcol = caloHitsHandle.product();
    }else{
      _chcol = 0;
    }

    //Get calorimeter readout hits MC level - energy/time/type
    art::Handle<CaloHitMCTruthCollection> caloHitMCTruthHandle;
    if (evt.getByLabel(_crmLabel, caloHitMCTruthHandle)){
      _chmccol = caloHitMCTruthHandle.product();
    }else {
      _chmccol = 0;
    }

    //Get stepPointMC for crystal readout hits
    art::Handle<PtrStepPointMCVectorCollection> mccaloptrHandle;
    if (evt.getByLabel(_crmLabel,_chmccpLabel,mccaloptrHandle)){
      _listOfMCCrystals = mccaloptrHandle.product();
    }else {
      _listOfMCCrystals = 0;
    }

    //Get simParticles and stepPointMC summary for crystal readout hits
    art::Handle<CaloHitSimPartMCCollection> caloHitSimMCHandle;
    if (evt.getByLabel(_crmLabel, caloHitSimMCHandle)){
      _chsmccol = caloHitSimMCHandle.product();
    }else {
      _chsmccol = 0;
    }

    if ((_chcol != 0) && (_chmccol != 0) && (_listOfMCCrystals != 0) && (_chsmccol != 0)){
      //      _caloHitNavigator = new CaloHitMCNavigator(*_chcol, *_chmccol, *_listOfMCCrystals, *_chsmccol);
      _caloHitNavigator = new CaloHitMCNavigator(*_chcol, *_chmccol, *_chsmccol);
    }else {
      _caloHitNavigator = 0;
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
    const char*               oname = "CalPatRec::produce";
    char                      message[200];
    bool                      findhelix (false), findseed (false), findkal (false);
    int                       nhits;
    int                       npeaks;
    int                       gen_index, sim_id;

    ::KalRep*                 krep;
                                        // dummy objects
    static TrkDefHack             dummydef;
    static HelixDefHack       dummyhdef;

    static HelixFitHackResult dummyhfit(dummyhdef);
    static KalFitResult       dummykfit(&dummydef);

    static StrawHitFlag       esel(StrawHitFlag::energysel), flag;

    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    fgTimeOffsets->updateMap(event);

    _ntracks = 0;
    _nhits_from_gen   = 0;

    fNCaloEnergyCut   = 0;
    fNCaloSizeCut     = 0;
    fNHitsTimePeakCut = 0;
    fNTimeWindow      = 0;

    for (int i=0; i<3; ++i){
      fSHSel[i] = 0;
      fSHBkg[i] = 0;
    }
                                        // reset the fit iteration counter
    _kfit.setNIter(0);
//     t1 = fStopwatch->RealTime();
//     fStopwatch->Continue();
                                        // event printout
    _eventid = event.event();
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : [START] : run = %10i subrun = %10i event = %10i\n",
				       oname,event.run(),event.subRun(),_eventid);

    unique_ptr<KalRepCollection>       tracks   (new KalRepCollection     );
    unique_ptr<KalRepPtrCollection>    trackPtrs(new KalRepPtrCollection  );
    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection);

    _tpeaks = new CalTimePeakCollection;
    unique_ptr<CalTimePeakCollection>  tpeaks(_tpeaks);

    _flags = new StrawHitFlagCollection();
    unique_ptr<StrawHitFlagCollection> flags (_flags);

    art::ProductID kalRepsID(getProductID<KalRepCollection>(event));

    double pEntrance(.0), step_time(-9999.);
    double time_threshold(500.);
                                        // find the data
    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n",oname);
                                                            goto END;
    }
//-----------------------------------------------------------------------------
// count the number of MC straw hits generated by the CE
//-----------------------------------------------------------------------------
    if (_listOfMCStrawHits == 0){
      nhits = 0;
    } else{
      nhits = _shcol->size();
    }

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
      if (gen_index >0 && sim_id == 1) {
        step_time = fgTimeOffsets->timeWithOffsetsApplied(*Step);
        step_time = fmod(step_time,_mbtime);
        if (step_time > time_threshold) {
          ++_nhits_from_gen;
          if (Step->momentum().mag() > pEntrance) {
            pEntrance = Step->momentum().mag();
          }
        }
      }
    }
    if (pEntrance < 100. ) _nhits_from_gen = 0;

    if (_diagLevel > 0) {
      if (_nhits_from_gen >= _minNMCHits)  _hist._cutflow[0]->Fill(0.0);
    }

    _kfit.setStepPointMCVectorCollection(_listOfMCStrawHits);
    _seedfit.setStepPointMCVectorCollection(_listOfMCStrawHits);
//-----------------------------------------------------------------------------
// all needed pieces of data have been found,
// tighten the energy cut and copy flags, clear
//-----------------------------------------------------------------------------
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
    if (_diagLevel > 0) {
      fillTimeDiag();
    }
//-----------------------------------------------------------------------------
// loop over found time peaks - for us, - "eligible" calorimeter clusters
//-----------------------------------------------------------------------------
    npeaks = _tpeaks->size();
    _clCE  = false;

    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      CalTimePeak* tp = &_tpeaks->at(ipeak);

      int nhits = tp->NHits();
      if (nhits > 400) {
	printf ("[%s] ERROR: Nhits(timepeak) = %i. SKIP, continue with the next time peak\n",
		oname,nhits);
	continue;
      }
//-----------------------------------------------------------------------------
// this is debug-only
//-----------------------------------------------------------------------------
      if (_diagLevel > 0) {
	CaloContentMC clutil(*_caloHitNavigator, *tp->Cluster());
	_clCE    = clutil.hasConversion();
      }

      if (_debugLevel > 0) {
        const StrawHit*     hit;
        int nhits = tp->NHits();
        printf(" peak # ipeak = %2i; nhits = %5i\n",ipeak,nhits);
        if (_debugLevel > 1) {
          for (int ih=0; ih<nhits; ih++) {
            StrawHitIndex ind = tp->_index[ih];
            hit  = &_shcol->at(ind);
	    flag = _shfcol->at(ind);
            printf("index = %5i flag: %10s time=%10.3f energy = %10.3f\n",
                   hit->strawIndex().asInt(),flag.hex().data(), hit->time(),hit->energyDep());
          }
        }
      }
//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information
//-----------------------------------------------------------------------------
      HelixDefHack helixdef(_shcol,_shpcol,_flags,tp->_index,_tpart,_fdir);

      TrkDefHack             seeddef(helixdef);
      TrkDefHack             kaldef (helixdef);

                                        // track fitting objects for this peak

      HelixFitHackResult hf_result(helixdef);

      init(_sfresult, &seeddef);
      init(_kfresult, &kaldef );
//-----------------------------------------------------------------------------
// Step 1: pattern recognition. Find initial helical approximation of a track
//-----------------------------------------------------------------------------
      int rc = _hfit.findHelix(hf_result,tp);

      if (_debugLevel > 0) {
        printf("[CalPatRec::produce] helixFit status = %i\n", rc);
        _hfit.printInfo(hf_result);
      }

      if (rc) {
//-----------------------------------------------------------------------------
// pattern recognition succeeded, the seed fit starts
//-----------------------------------------------------------------------------
        findhelix = true;
//-----------------------------------------------------------------------------
// convert the result to standard helix parameters, and initialize the seed definition helix
//-----------------------------------------------------------------------------
        if (_diagLevel > 0) {
          double dz, dist;
          for (int i=0; i< fHackData->goodPoints(); ++i){
            dz   = fHackData->fDz[i];
            dist = fHackData->fDist[i];
            _hist._distvsdz->Fill(dz, dist);
          }
        }

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
        std::vector<StrawHitIndex> goodhits;

        _index = _hfit._xyzp;

        std::sort(_index.begin(), _index.end(), [ ]( const XYZPHack& lhs,
                                                     const XYZPHack& rhs )
                  {
                    return lhs._ind < rhs._ind;
                  } );

        _nindex = _index.size();

        if (_debugLevel > 0) {

	  Hep3Vector*            shPos;
	  int                    loc;
	  const mu2e::StrawHit*  hit;
	  const mu2e::Straw*     straw;

          printf("[CalPatRec::printGoodHits]   Index       Flag  Straw    status       X         Y         Z         Zw\n");
	  for (int i=0; i< _nindex; ++i) {
	    //	    if (_index[i].isOutlier()) continue;

            shPos = &_index[i]._pos;
            loc   = _index[i]._ind;
            hit   = &_shcol->at(loc);
            straw = &_tracker->getStraw(hit->strawIndex());
            printf("[CalPatRec::printGoodHits]  %6i  %10s %6i    active   %8.3f  %8.3f  %9.3f %9.3f\n",
                   loc, 
		   _index[i]._flag.hex().data(),
		   hit->strawIndex().asInt(), 
		   shPos->x(), shPos->y(), shPos->z(),
                   straw->getMidPoint().z()
                   );
	  }
        }

        for (int i=0; i< _nindex; ++i){
          if (_index[i].isOutlier()) continue;
          goodhits.push_back(_index[i]._ind);
        }
        seeddef.setIndices (goodhits);
//-----------------------------------------------------------------------------
//  Trajectory info
//-----------------------------------------------------------------------------
        if (_diagLevel > 0) {
          Hep3Vector tdir;
          HepPoint   tpos;
          double     doca;
          _helTraj->getInfo(0.0,tpos,tdir);

          for (int i=0; i< _nindex; ++i){
            StrawHit const*   sh    = _index[i]._strawhit;
            Straw const&      straw = _tracker->getStraw(sh->strawIndex());
            CLHEP::Hep3Vector wpos  = straw.getMidPoint();
            CLHEP::Hep3Vector wdir  = straw.getDirection();

            // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
            HepPoint      wpt  (wpos.x(),wpos.y(),wpos.z());
            TrkLineTraj   wtraj(wpt,wdir,-20,20);

            // estimate flightlength along track.  This assumes a constant BField!!!
            double fltlen = (wpos.z()-tpos.z())/tdir.z();

            TrkPoca   wpoca(*_helTraj,fltlen,wtraj,0.0);

            doca      = wpoca.doca();

            if (_index[i].isOutlier()) _hist._doca[1]->Fill(doca);
            else                       _hist._doca[0]->Fill(doca);
          }
        }

        if (_debugLevel > 0) {
          printf("[CalPatRec::seeddef] goodhits = %lu over nIndex = %i\n", goodhits.size(), _nindex);
          int shIndices = seeddef.strawHitIndices().size();
          int nSh       = seeddef.strawHitCollection()->size();
          printf("[CalPatRec::seedfit] START Nstraws: = %i Nindices: %i", nSh, shIndices);
	  printf(" hpar: ");
          for(int i=0; i<5; ++i) printf("%10.5f", hpar[i]);
          printf("\n");
        }
//-----------------------------------------------------------------------------
// seed fit - fit through the wires of found hits, not using the drift times
//-----------------------------------------------------------------------------
        _seedfit.makeTrack(*_sfresult, tp->Cluster());
//--------------------------------------------------------------------------------
// 2014-11-24 gianipez added the following diagnnostic
//--------------------------------------------------------------------------------
        if (_debugLevel > 0) {
          sprintf(message,
                  "CalPatRec::produce seedfit::makeTrack : fit_success = %i\n",
                  _sfresult->_fit.success());
          _seedfit.printHits(*_sfresult,message);
        }
//--------------------------------------------------------------------------------
// 2015-03-23 G. Pezzu: fill info about the doca
//--------------------------------------------------------------------------------
        if (_sfresult->_fit.success()) {
          if (_clCE)  findseed = true;
//-----------------------------------------------------------------------------
// use helix parameters by the seed fit to initialize the full Kalman fit
//-----------------------------------------------------------------------------
          double           locflt;
          const HelixTraj* shelix;
          shelix = (const HelixTraj*) _sfresult->_krep->localTrajectory(_sfresult->_krep->flt0(),locflt);
          kaldef.setHelix(*shelix);
//----------------------------------------------------------------------
//2015-02-07 G. Pezzu added new selection using seedFit results
//----------------------------------------------------------------------
          goodhits.clear();
          const mu2e::TrkStrawHit* hit;
          int                      hit_index;
	  TrkHitVector const&  hot_l = _sfresult->_krep->hitVector();
          const StrawHit*          sh;
          const Straw*             straw;
          const CLHEP::Hep3Vector  *wpos, *wdir;

                                        //  Trajectory info
          Hep3Vector               tdir;
          HepPoint                 tpos;
          double                   doca, rdrift, fltlen;
          bool                     found(false), active;
          int                      banner_11_printed(0);

          _sfresult->_krep->traj().getInfo(0.0,tpos,tdir);

          _nrescued = 0;

          for (int i=0; i< _nindex; ++i) {
            hit_index = _index[i]._ind;
            sh        = _index[i]._strawhit;
            straw     = &_tracker->getStraw(sh->strawIndex());
            wpos      = &straw->getMidPoint();
            wdir      = &straw->getDirection();
            rdrift    = -9990;
            found     = false;
            active    = false;

            HepPoint      wpt  (wpos->x(),wpos->y(),wpos->z());
            TrkLineTraj   wtraj(wpt,*wdir,-20,20);
//-----------------------------------------------------------------------------
// estimate flightlength along track. This assumes a constant BField!!!
// in principle, this should work well enough, however, may want to check
//-----------------------------------------------------------------------------
            fltlen = (wpos->z()-tpos.z())/tdir.z();

            TrkPoca   wpoca(_sfresult->_krep->traj(),fltlen,wtraj,0.0);

            doca = wpoca.doca();

            for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
              hit    = static_cast<const mu2e::TrkStrawHit*> (*it);
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

            if (_debugLevel > 0) {
              if (banner_11_printed == 0) {
                banner_11_printed = 1;
                printf("[CalPatRec::produce] -------------------------------------\n");
                printf("[CalPatRec::produce]  ih  A   Sind      Rdrift        doca\n");
                printf("[CalPatRec::produce] -------------------------------------\n");
              }

              printf("[CalPatRec::produce]  %2i  %1i  %5i  %10.3f  %10.3f \n",
                     i, active? 1:0, straw->index().asInt(), rdrift, doca );
            }

            if (_diagLevel > 0) {
              _hist._seeddoca[2]->Fill(doca);

              if (found)  _hist._seeddoca[0]->Fill(doca);
              else        _hist._seeddoca[1]->Fill(doca);

              if (active) _hist._seeddoca[0]->Fill(doca);
              else        _hist._seeddoca[1]->Fill(doca);
            }

            if (!found && active) ++_nrescued;
          }

          if (_diagLevel > 0) fillSeedFitHistograms(*_sfresult);
//-----------------------------------------------------------------------------
// at this point the full kalman fit starts
//-----------------------------------------------------------------------------
          kaldef.setIndices(goodhits);
          if (_debugLevel > 0) printf("CalPatRec::produce] calling _kfit.makeTrack\n");

          _kfit.makeTrack(*_kfresult,tp->Cluster());

          if (_debugLevel > 0) {
            printf("[CalPatRec::produce] kalfit status = %i\n", _kfresult->_fit.success());
            _kfit.printHits(*_kfresult,"CalPatRec::produce kalfit_001");
          }

          if (_kfresult->_fit.success()) {
	    if (_clCE)  findkal = true;

	    if (_addhits) {
//-----------------------------------------------------------------------------
// this is the default. First, add back the hits on this track
// if successfull, try to add missing hits, at this point external errors were
// set to zero
// assume this is the last iteration
//-----------------------------------------------------------------------------
//	      int last_iteration = _kfit.maxIteration();
	      int last_iteration = -1;
	      
	      _kfit.unweedHits(*_kfresult,_maxaddchi);
	      if (_debugLevel > 0) _kfit.printHits(*_kfresult,"CalPatRec::produce after unweedHits");
	      
	      std::vector<StrawHitIndex> misshits;
	      findMissingHits(*_kfresult,misshits);
//-----------------------------------------------------------------------------
// if new hits have been added, add then and refit the track.
// Otherwise - just refit the track one last time
// in both cases
//-----------------------------------------------------------------------------
	      if (misshits.size() > 0) {
		_kfit.addHits(*_kfresult,_shcol,misshits, _maxaddchi, tp->Cluster());
	      }
	      else {
		_kfit.fitIteration(*_kfresult,last_iteration,tp->Cluster());
	      }

	      if (_debugLevel > 0) _kfit.printHits(*_kfresult,"CalPatRec::produce after addHits");
//-----------------------------------------------------------------------------
// and weed hits again to insure that addHits doesn't add junk
//-----------------------------------------------------------------------------
	      _kfit.weedHits(*_kfresult,last_iteration);
	    }
//-----------------------------------------------------------------------------
// now evaluate the T0 and its error using the straw hits
//-----------------------------------------------------------------------------
	    _kfit.updateT0(*_kfresult);
	    
	    if (_debugLevel > 0) {
	      _kfit.printHits(*_kfresult,"CalPatRec::produce : final, after weedHits");
	    }
//-----------------------------------------------------------------------------
// done, fill debug histograms
//-----------------------------------------------------------------------------
	    TrkHitVector const& hot_l = _sfresult->_krep->hitVector();

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
	      TrkPoca hitpoca(_kfresult->_krep->traj(),fltlen,htraj,0.0);

	      doca = hitpoca.doca();
	      for(auto it=hot_l.begin(); it<hot_l.end(); it++) {
		hit = static_cast<const mu2e::TrkStrawHit*> (*it);
		if (!hit->isActive()) continue;
		hit_index = hit->index();
		if (int(_index[i]._ind) == hit_index){
		  found = true;
		  break;
		}
	      }
	    
	      if (_diagLevel > 0) {
		if (found) _hist._kaldoca[0]->Fill(doca);
		else       _hist._kaldoca[1]->Fill(doca);
	      }
	    }
	  }
	}
      }

      if (_kfresult->_fit.success()) {
        _ntracks += 1;
//-----------------------------------------------------------------------------
// fit succeeded
// flag hits used in this track.  This should use the track id, FIXME!!! (in the BaBar code)
//-----------------------------------------------------------------------------
        if (ipeak<16) {
          for (size_t ihit=0; ihit<_kfresult->_hits.size(); ++ihit){
            TrkStrawHit* tsh = dynamic_cast<TrkStrawHit*>(_kfresult->_hits[ihit]);
            if (tsh == 0)      continue;
	    if (tsh->isActive()) {
              _flags->at(tsh->index()).merge(StrawHitFlag::trackBit(ipeak));
              _flags->at(tsh->index()).merge(StrawHitFlag::calosel);
            }
          }
        }
//-----------------------------------------------------------------------------
//  fill fit diagnostics histograms if requested
//-----------------------------------------------------------------------------
        if (_diagLevel > 0) fillFitDiag(event,ipeak,hf_result,*_sfresult,*_kfresult);
//-----------------------------------------------------------------------------
// save successful kalman fits in the event.
// start from _kfresult, as stealTrack clears the hit pointers
// _kfresult doesn't own anything
//-----------------------------------------------------------------------------
        krep = _kfresult->stealTrack();
        tracks->push_back(krep);
        int index = tracks->size()-1;
        trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
        tp->SetCprIndex(tracks->size());

	int best = AlgorithmID::CalPatRecBit;
	int mask = 1 << AlgorithmID::CalPatRecBit;

	algs->push_back(AlgorithmID(best,mask));
	
      }
      else {
//-----------------------------------------------------------------------------
// fit failed, just delete the track
//-----------------------------------------------------------------------------
        _kfresult->deleteTrack();
      }
//-----------------------------------------------------------------------------
// cleanup the seed fit - why it is not being done ?
//-----------------------------------------------------------------------------
//       _sfresult->deleteTrack();

      if (_debugLevel > 0) {
	if (_nhits_from_gen >= _minNMCHits) {
	  if (tp->_tmin > 400.){
	    if (!findhelix) {
	      printf("[CalPatRec::produce] LOOK AT: more than 25 MC hits and findHelix not converged! event = %i\n", _iev);
	    }
	    if (findhelix && !findseed){
	      printf("[CalPatRec::produce] LOOK AT: findhelix converged and findseed not! event = %i\n", _iev);
	    }
	    if (findseed && !findkal){
	      printf("[CalPatRec::produce] LOOK AT: findseed converged and findkal not! event = %i\n", _iev);
	    }
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// diagnostics in the end
//-----------------------------------------------------------------------------
    if (_diagLevel > 0) {
      if (findhelix && (_nhits_from_gen >= _minNMCHits)) {
	_hist._cutflow[0]->Fill(2.0);
	_hist._cutflow[1]->Fill(9.0);
      }

      if (findseed  && (_nhits_from_gen >= _minNMCHits)) _hist._cutflow[0]->Fill(3.0);
      if (findkal   && (_nhits_from_gen >= _minNMCHits)) {
	_hist._cutflow[0]->Fill(4.0);
	if (fQualityTrack > 0) {
	  _hist._cutflow[0]->Fill(5.0);
	}
	_hist._NfitIter->Fill(_kfit.nIter());
	_hist._ntracks->Fill(_ntracks);
      }
    }
//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record
//-----------------------------------------------------------------------------
  END:;
    event.put(std::move(tracks)   );
    event.put(std::move(trackPtrs));
    event.put(std::move(algs     ));
    event.put(std::move(flags )   );
    event.put(std::move(tpeaks)   );
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

    StrawHitFlag        energyFlag(StrawHitFlag::energysel);
    StrawHitFlag        timeFlag  (StrawHitFlag::timesel);
    StrawHitFlag        radiusFlag(StrawHitFlag::radsel);
    StrawHitFlag        deltaRayFlag(StrawHitFlag::delta);
    StrawHitFlag        isolatedFlag(StrawHitFlag::isolated);
//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
    nsh   = _shcol->size();
    ncl   = _ccCollection->size();
    _clCE = 0;
    for (int ic=0; ic<ncl; ic++) {
      cl      = &_ccCollection->at(ic);

      if (_diagLevel > 0) {
	CaloContentMC clutil(*_caloHitNavigator, *cl);
	_clCE    = clutil.hasConversion();
      }

      if ( cl->energyDep() > _minClusterEnergy) {

        if ((_diagLevel > 0) && _clCE) ++fNCaloEnergyCut;

        if ( (int(cl->size()) >= _minClusterSize) ) {
          if ((_diagLevel > 0) && _clCE) ++fNCaloSizeCut;

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

          //    dz_cl   = zcl; // -_tracker->z0();
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
          int   nhitsTimeWindow(0), nhitsHasTime(0), nhitsHasEnergy(0), nhitsHasRadius(0),
                nhitsNoDelta(0), nhitsNoIsolated(0);

          double meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity
          int    gen_index, sim_id, vol_id;

          for(int istr=0; istr<nsh;++istr) {
            flag = _flags->at(istr);

            int hit_has_all_properties = flag.hasAllProperties(_hsel);
            int bgr_hit                = flag.hasAnyProperty(_bkgsel);

            int hit_has_energy         = flag.hasAllProperties(energyFlag);
            int hit_has_time           = flag.hasAllProperties(timeFlag);
            int hit_has_radius         = flag.hasAllProperties(radiusFlag);

            int deltaRay_hit           = flag.hasAnyProperty(deltaRayFlag);
            int isolated_hit           = flag.hasAnyProperty(isolatedFlag);

            hit    = &_shcol->at(istr);
            time   = hit->time();
            straw  = &_tracker->getStraw(hit->strawIndex());
            zstraw = straw->getMidPoint().z();
//-----------------------------------------------------------------------------
// estimate time-of-flight and calculate residual between the predicted and the hit times
//-----------------------------------------------------------------------------
            tof = (zcl-zstraw)/sin(_pitchAngle)/CLHEP::c_light;
            dt  = cl_time-(time+tof-meanDriftTime);
//-----------------------------------------------------------------------------
// fill some diag histograms
//-----------------------------------------------------------------------------
            if ((_diagLevel > 0) && _clCE) {
					// is the straw hit from CE?
              const mu2e::StepPointMC* step(0);
              int nstraws = _listOfMCStrawHits->size();
              for (int i=0; i<nstraws; i++) {
                mu2e::PtrStepPointMCVector  const& mcptr(_listOfMCStrawHits->at(i));
                step = &(*mcptr.at(0));
                vol_id = step->volumeId();
                if (vol_id == straw->index().asInt()) {
					// step found - use the first one in the straw
                  break;
                }
              }

              gen_index = -1;
              if (step) {
                art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();

                if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
                else                         gen_index = -1;

                sim_id        = simptr->id().asInt();
              }

              if (gen_index >0 && sim_id == 1) {
                // hit from CE
                if ((dt < -70.) && (step->momentum().z()>0.) && (step->momentum().mag()>80.)){
                  if (_debugLevel > 0) {
                    printf("Event : %10i dt = %5.3f\n", _eventid, dt);
                  }
                }
                if ( (step->momentum().mag()>80.)&& (step->momentum().z()>0.) ) _hist._dt[1]->Fill(dt);
              }

              if (hit_has_all_properties && !bgr_hit) {
                _hist._dt[0]->Fill(dt);
              }
            }

            if ((dt < _maxdt) && (dt >= _mindt)) {
              if ((_diagLevel > 0) && _clCE){
                ++nhitsTimeWindow;

                if ( hit_has_time   ) ++nhitsHasTime;
                if ( hit_has_radius ) ++nhitsHasRadius;
                if ( hit_has_energy ) ++nhitsHasEnergy;
                if ( !deltaRay_hit  ) ++nhitsNoDelta;
                if ( !isolated_hit  ) ++nhitsNoIsolated;
              }

              if (hit_has_all_properties && !bgr_hit) {
                tpeak._index.push_back(istr);
                stime += time;
              }
            }
          }

          tpeak._tpeak = stime/(tpeak.NHits()+1.e-12);

          if ((_diagLevel > 0) && _clCE) {
            if (nhitsTimeWindow>_minnhits) ++fNTimeWindow;
            if (nhitsHasTime   >_minnhits) fSHSel[0] = 1;
            if (nhitsHasRadius >_minnhits) fSHSel[1] = 1;
            if (nhitsHasEnergy >_minnhits) fSHSel[2] = 1;
            if (nhitsNoDelta   >_minnhits) fSHBkg[0] = 1;
            if (nhitsNoIsolated>_minnhits) fSHBkg[1] = 1;
          }

          if (tpeak.NHits() > _minnhits) {
           if ((_diagLevel > 0) && _clCE) {
             ++fNHitsTimePeakCut;
           }

           TimePeakColl->push_back(tpeak);
          }
        }
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
  void CalPatRec::filterOutliers(TrkDefHack&                    mytrk  ,
                                 Trajectory const&          traj   ,
                                 double                     maxdoca,
                                 std::vector<TrkHitFilter>& thfvec ) {
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;
    traj.getInfo(0.0,tpos,tdir);

    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    const StrawHitCollection* hits = mytrk.strawHitCollection();
    const std::vector<StrawHitIndex>& indices = mytrk.strawHitIndices();
    std::vector<StrawHitIndex> goodhits;

    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = hits->at(indices[ihit]);
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
      if(_diagLevel > 0){
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
  void CalPatRec::findMissingHits(KalFitResult& kalfit,std::vector<StrawHitIndex>& misshits) {
                                        //  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;
    int        radius_ok;
    double     dt;

    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    if (_debugLevel > 0) {
      printf("[CalPatRec::findMissingHits]      shId    sec     panel       doca   \n");
    }

    for (unsigned istr=0; istr<nstrs;++istr) {
//----------------------------------------------------------------------
// 2015-02-11 gianipez and P. Murat changed the selection bit
//            for searching for missed hits
//----------------------------------------------------------------------
      StrawHit const& sh = _shcol->at(istr);
//-----------------------------------------------------------------------------
// I think, we want to check the radial bit: if it is set, than at least one of
// the two measured times is wrong...
//-----------------------------------------------------------------------------
      radius_ok = _shfcol->at(istr).hasAllProperties(StrawHitFlag::radsel);
      dt        = _shcol->at(istr).time()-kalfit._krep->t0()._t0;

      if (radius_ok && (fabs(dt) < _maxdtmiss)) {
        // make sure we haven't already used this hit
	TrkStrawHitVector tshv;
	convert(kalfit._hits, tshv);
        TrkStrawHitVector::iterator ifnd = find_if(tshv.begin(), tshv.end(),FindTrkStrawHit(sh));
        if(ifnd == tshv.end()){
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

          if (_debugLevel > 0) {
            printf("[CalPatRec::findMissingHits] %8i  %6i  %8i  %10.3f \n",
                   straw.index().asInt(),
                   straw.id().getPlane(),
                   straw.id().getPanel(),
                   hitpoca.doca());
          }
//-----------------------------------------------------------------------------
// flag hits with small residuals
//-----------------------------------------------------------------------------
          if (fabs(hitpoca.doca()) < _maxadddoca) {
            misshits.push_back(istr);
          }
        }
      }
    }
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::fillStrawDiag() {

  }

//-----------------------------------------------------------------------------
  void CalPatRec::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;

    bool conversion(false);
    StrawHitFlag  flag;

    int nhits = _shcol->size();

    _hist._nhits[0]->Fill(nhits);
    _hist._nhits[1]->Fill(nhits);

    for (int i=0; i<nhits; i++) {
      double time = _shcol->at(i).time();

      _hist.rtsp->Fill(time);

      flag = _flags->at(i);

      if (flag.hasAllProperties(_tsel)) _hist.ttsp->Fill(time);

      if (flag.hasAllProperties(_tsel) && ! flag.hasAnyProperty(_bkgsel)) _hist.tdtsp->Fill(time);
      if (flag.hasAllProperties(_ksel) && ! flag.hasAnyProperty(_bkgsel)) _hist.ltsp->Fill(time);

      if (conversion) _hist.ctsp->Fill(time);
    }

    _hist._Tpeaks->Fill(_tpeaks->size());

    if (_nhits_from_gen >= _minNMCHits) {
      _hist._cutflow[1]->Fill(.0);
      if (fNCaloEnergyCut>0){
	_hist._cutflow[1]->Fill(1.0);
	if (fNCaloSizeCut>0){
	  _hist._cutflow[1]->Fill(2.0);
	  if (fNTimeWindow>0){
	    _hist._cutflow[1]->Fill(3.0);
	    if (fSHSel[1] > 0){
	      _hist._cutflow[1]->Fill(4.0);
	      if (fSHSel[2] > 0){
		_hist._cutflow[1]->Fill(5.0);
		if (fSHBkg[0] > 0){
		  _hist._cutflow[1]->Fill(6.0);
		  if (fSHBkg[1] > 0){
		    _hist._cutflow[1]->Fill(7.0);
		    if (fSHSel[0] > 0){
		      _hist._cutflow[1]->Fill(8.0);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      if (_tpeaks->size() > 0 ) _hist._cutflow[0]->Fill(1.0);
    }
  }

//-----------------------------------------------------------------------------
// seed fit diagnostic histograms, DOCA is the distance to the wire
//-----------------------------------------------------------------------------
  void CalPatRec::fillSeedFitHistograms(KalFitResult& SFResult) {

    int                      ndeactivated(0);
    double                   doca, fltlen;
    Hep3Vector               tdir;
    HepPoint                 tpos;
    TrkHitVector const& hot_l = SFResult._krep->hitVector();
    const mu2e::TrkStrawHit* hit;
    KalRep*                  krep;
    //    const HelixTraj*         shelix;

    krep     = SFResult._krep;

    //    shelix   = (const HelixTraj*) krep->localTrajectory(krep->flt0(),locflt);
    krep->traj().getInfo(0.0,tpos,tdir);

    for (int i=0; i<_nindex; ++i){
      StrawHit const*    sh    = _index[i]._strawhit;
      Straw const&       straw = _tracker->getStraw(sh->strawIndex());
      const CLHEP::Hep3Vector* wpos  = &straw.getMidPoint();
      const CLHEP::Hep3Vector* wdir  = &straw.getDirection();

      HepPoint      wpt  (wpos->x(),wpos->y(),wpos->z());
      TrkLineTraj   wtraj(wpt,*wdir,-20,20);

      fltlen  = (wpos->z()-tpos.z())/tdir.z();

      TrkPoca    wpoca(*_helTraj,fltlen,wtraj,0.0);

      doca   = wpoca.doca();
      if (_index[i].isOutlier()) _hist._doca[1]->Fill(doca);
      else                       _hist._doca[0]->Fill(doca);
    }

    for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
      hit = static_cast<const mu2e::TrkStrawHit*> (*it);
      if (!hit->isActive()) ++ndeactivated;
    }

    _hist._NpointsSeed[0]->Fill(ndeactivated);
    _hist._NpointsSeed[1]->Fill(_nrescued);
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

    _sfresult->_krep->traj().getInfo(0.0,tpos,tdir);

    fQualityTrack = 0;
//-----------------------------------------------------------------------------
// HelixFit histograms
//-----------------------------------------------------------------------------
    _hist.helixFit.nhits->Fill(hf_result._nGoodPoints);
//-----------------------------------------------------------------------------
// KalFit   histograms
//-----------------------------------------------------------------------------
    const TrkStrawHit*     hit;
    int                    hit_index;
    double                 fltlen;

    KalRep*    krep = _kfresult->_krep;

    TrkHitVector const& hotl = krep->hitVector();

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
      for(auto it=hotl.begin();it<hotl.end(); it++) {
        hit   = static_cast<const mu2e::TrkStrawHit*> (*it);
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
    Hep3Vector mom = krep->momentum(0);
    double pt      = sqrt(mom.x()*mom.x() + mom.y()*mom.y());
    double pz      = mom.z();
    double tanL    = pt/pz;
    double radius  = fabs(1./krep->helix(0).omega());//mom.mag()*10./3.;//convert MeV to mm
    double kdfdz   = tanL/radius;

    TrkDifTraj const &traj = krep->traj();

    double    trkPhi, hitPhi, dphi;
    HepPoint  trkPos;
    double    tmpFltz0(-9999.);
    double    dz = 1e10;

    for(auto it=hotl.begin();it<hotl.end(); it++) {
      hit   = static_cast<const mu2e::TrkStrawHit*> (*it);
      if (!hit->isActive())            continue;
      fltlen = hit->fltLen();
      trkPos = traj.position(fltlen);

      if ( fabs(trkPos.z()) < dz ){
        tmpFltz0 = fltlen;
        dz       =  fabs(trkPos.z());
      }
    }

    HelixTraj trkHel(krep->helix(tmpFltz0).params(),krep->helix(tmpFltz0).covariance());
    trkHel.getInfo(tmpFltz0, tpos,tdir);

    double   fltz0 = tmpFltz0 - traj.position(tmpFltz0).z()/tdir.z();

    double d0      = krep->helix(fltz0).d0();
    double phiC    = krep->helix(fltz0).phi0();
    double x0      = -(radius + d0)*std::sin(phiC);
    double y0      =  (radius + d0)*std::cos(phiC);

    HepPoint   trkPosz0 = traj.position(fltz0);
    CLHEP::Hep3Vector trkHelPosZ0 (trkPosz0.x(), trkPosz0.y(), trkPosz0.z());
    CLHEP::Hep3Vector trkHelCenter(x0, y0, 0.);

    double kphi0   = CLHEP::Hep3Vector(trkHelPosZ0-trkHelCenter).phi() - trkPosz0.z()*kdfdz;

    while (kphi0 >  M_PI) kphi0 -= 2*M_PI;
    while (kphi0 < -M_PI) kphi0 += 2*M_PI;

    double ptrPhi0[4], dPhi0[4];
    ptrPhi0[0] = (fHackData->fData[22]);
    ptrPhi0[1] = (fHackData->fData[23]);
    ptrPhi0[2] = (fHackData->fData[24]);
    ptrPhi0[3] = (fHackData->phi0());

    for (int i=0; i<4; ++i){
      while (ptrPhi0[i] >  M_PI) ptrPhi0[i] -= 2*M_PI;
      while (ptrPhi0[i] < -M_PI) ptrPhi0[i] += 2*M_PI;
      dPhi0[i] = (ptrPhi0[i]) - (kphi0);
    }
//--------------------------------------------------------------------------------
// calculate the phi residual of the straw hit position using the track result
//--------------------------------------------------------------------------------
    for (auto it=hotl.begin();it<hotl.end(); it++) {
      hit   = static_cast<const mu2e::TrkStrawHit*> (*it);
      if (!hit->isActive())            continue;

      //calculate the hit phi
      hit_index = hit->index();

      StrawHitPosition const& shp = _shpcol->at(hit_index);
      CLHEP::Hep3Vector    hitPos = shp.pos();
      hitPhi = CLHEP::Hep3Vector(hitPos - trkHelCenter).phi();
      while (hitPhi >  M_PI) hitPhi  -= 2*M_PI;
      while (hitPhi < -M_PI) hitPhi  += 2*M_PI;

      //calculate the track phi at a given z
      fltlen = hit->fltLen();
      trkPos = traj.position(fltlen);
      CLHEP::Hep3Vector trkHelPos (trkPos.x(), trkPos.y(), trkPos.z());

      trkPhi  = CLHEP::Hep3Vector(trkHelPos-trkHelCenter).phi();
      while (trkPhi >  M_PI) trkPhi  -= 2*M_PI;
      while (trkPhi < -M_PI) trkPhi  += 2*M_PI;

      dphi = hitPhi - trkPhi;
      _hist._PhiResid [2]->Fill(dphi);
    }

    // _hist._kdfdz  [0]->Fill(fHackData->dfdz() - kdfdz);
    _hist._kradius[0]->Fill(fHackData->fData[20] - radius);
    //    _hist._k0mode    ->Fill(fHackData->mode0Points());
//------------------------------------------------------------------------------------------
//take info from the seed fit for filling diagnostic histograms
//------------------------------------------------------------------------------------------
    seedMom      = _sfresult->_krep->momentum(0);
    seedPt       = sqrt(seedMom.x()*seedMom.x() + seedMom.y()*seedMom.y());
    seedPz       = seedMom.z();
    seedTanL     = seedPt/seedPz;
    seedRadius   = fabs(1./_sfresult->_krep->helix(0).omega());//mom.mag()*10./3.;//convert MeV to mm
    seeddfdz     = seedTanL/seedRadius;

    _hist._seeddfdz[0]->Fill(seeddfdz - kdfdz);
    _hist._seeddr  [0]->Fill(seedRadius - radius);
    _hist._drw     [0]->Fill(fHackData->fData[14]-radius);
    _hist._chi2w   [0]->Fill(fHackData->fData[15]);
    _hist._chi2zphi[0]->Fill(fHackData->fData[13]);

    _hist._kdz     ->Fill(fHackData->shDz());
    _hist._kNpoints->Fill(fHackData->goodPoints());
    //    _hist._kchi2   ->Fill(fHackData->chi2());

    _hist._NpointsRescued[0]->Fill(fHackData->rescuedPoints()/double(fHackData->goodPoints()));

    _hist._dphidz[0]->Fill(fHackData->fData[17]- kdfdz);
    _hist._dphidz[1]->Fill(fHackData->fData[18]- kdfdz);
    _hist._dphidz[2]->Fill(fHackData->fData[19]- kdfdz);
    _hist._dphidz[3]->Fill(fHackData->dfdz()   - kdfdz);


    _hist._dphi0[0]->Fill(dPhi0[0]);
    _hist._dphi0[1]->Fill(dPhi0[1]);
    _hist._dphi0[2]->Fill(dPhi0[2]);
    _hist._dphi0[3]->Fill(dPhi0[3]);

    double dist;
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

    std::vector<Doublet>* list_of_doublets = &_kfresult->_listOfDoublets;
    nd = list_of_doublets->size();

    for (int i=0; i<nd; i++) {
      d  = &list_of_doublets->at(i);
      ns = d->fNStrawHits;

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
          //      inext    = ibest;
          dsl_best = dsl;
          // ibest    = is;
        }
        else if (fabs(dsl) < fabs(dsl_next)) {
          dsl_next = dsl;
          //      inext    = is;
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
    double  fMaxT0Err        = 0.9;             // ns
    double  fMaxFitMomErr    = 0.25;            // MeV
    double  fMinTanDip       = tan(M_PI/6.);    // 0.5773
    double  fMaxTanDip       = 1.0;
    double  fMinD1           = -80.;            //  mm
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

      _hist._kradius [1]->Fill(fHackData->fData[20] - radius);
      _hist._seeddfdz[1]->Fill(seeddfdz - kdfdz);
      _hist._seeddr  [1]->Fill(seedRadius - radius);
      _hist._drw     [1]->Fill(fHackData->fData[14]-radius);
      _hist._chi2w   [1]->Fill(fHackData->fData[15]);
      _hist._chi2zphi[1]->Fill(fHackData->fData[13]);

      _hist._kdphidz [0]->Fill(fHackData->fData[17]- kdfdz);
      _hist._kdphidz [1]->Fill(fHackData->fData[18]- kdfdz);
      _hist._kdphidz [2]->Fill(fHackData->fData[19]- kdfdz);
      _hist._kdphidz [3]->Fill(fHackData->dfdz()   - kdfdz);

      _hist._kdphi0[0]->Fill(dPhi0[0]);
      _hist._kdphi0[1]->Fill(dPhi0[1]);
      _hist._kdphi0[2]->Fill(dPhi0[2]);
      _hist._kdphi0[3]->Fill(dPhi0[3]);


      _hist._NpointsRescued[1]->Fill(fHackData->rescuedPoints()/double(fHackData->goodPoints()));
      for (int i=0; i< fHackData->goodPoints(); ++i){
        dz   = fHackData->fDz[i];
        dist = fHackData->fDist[i];
        dphi = fHackData->fResid[i];
        _hist._PhiResid[1]->Fill(dphi);
        _hist._kdistvsdz[1]->Fill(dz, dist);
      }
    }
  }


//-----------------------------------------------------------------------------
  void CalPatRec::init(KalFitResult*& KRes, TrkDefHack* TDef) {

    if (KRes != 0) {
      KRes->deleteTrack();
      delete KRes;
    }
    KRes = new KalFitResult();
    KRes->_tdef = TDef;
  }
}

using mu2e::CalPatRec;
DEFINE_ART_MODULE(CalPatRec);
