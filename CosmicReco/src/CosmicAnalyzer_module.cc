#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinder_types.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
//Mu2e Data Prods:
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
// MC Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
// Mu2e diagnostics
#include "TrkDiag/inc/ComboHitInfo.hh"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// ROOT incldues
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"
using namespace std; 

namespace mu2e 
{
  class CosmicAnalyzer : public art::EDAnalyzer {
    public:
      explicit CosmicAnalyzer(fhicl::ParameterSet const& pset);
      virtual ~CosmicAnalyzer();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private: 
      int  _diag;
      bool _mcdiag;
      
      //TTree Info:
      TTree* _cosmic_analysis;
      //seed Info:
      TH1F* _seed_a1;
      TH1F* _seed_b1;
      TH1F* _seed_a0;
      TH1F* _seed_b0;
      TH1F* _seed_a1XYZ;
      TH1F* _seed_b1XYZ;
      TH1F* _seed_a0XYZ;
      TH1F* _seed_b0XYZ;
      TH1F* _niters;
     //Look for geometric bias in seed fit:
      TH2F* _chiY_v_true_trackposY;
      TH2F* _chiX_v_true_trackposX;
      TH2F* _chi_v_true_theta;
     //Looking for hidden correlations between seed paramters
      TH2F* _seedA0_v_seedA1;
      TH2F* _seedB0_v_seedB1;
      TH2F* _seedA0_v_seedB1;
      TH2F* _seedB1_v_seedA1;
      TH2F* _seedA0_v_seedB0;
    //Looking for correlations between residuals and N hits:
      TH2F* _seedDeltaA0_v_N;
      TH2F* _seedDeltaA1_v_N;
      TH2F* _seedDeltaB0_v_N;
      TH2F* _seedDeltaB1_v_N;
    //Looking for correlations between seed residuals and true momentum:
      TH2F* _seedDeltaA0_v_MOM;
      TH2F* _seedDeltaA1_v_MOM;
      TH2F* _seedDeltaB0_v_MOM;
      TH2F* _seedDeltaB1_v_MOM;
    // True MC parameters in local coords
      TH1F* _truea1;
      TH1F* _trueb1;
      TH1F* _truea0;
      TH1F* _trueb0;
    //True MC paramets in global coords:
      TH1F* _truea1XYZ;
      TH1F* _trueb1XYZ;
      TH1F* _truea0XYZ;
      TH1F* _trueb0XYZ;
      TH1F* _SPtruea1XYZ;
      TH1F* _SPtrueb1XYZ;
      TH1F* _SPtruea0XYZ;
      TH1F* _SPtrueb0XYZ;
 
      //Angles Info:
      TH1F* _mc_phi_angle;
      TH1F* _reco_phi_angle;
      TH1F* _mc_theta_angle;
      TH1F* _reco_theta_angle;

      //Looking for correvariation between seed parameter and true parameters:
      TH1F* _A0SeedMCDiff;
      TH1F* _A1SeedMCDiff;
      TH1F* _B0SeedMCDiff;
      TH1F* _B1SeedMCDiff;
     //Look at distributions of covarience of each seed parameter:
      TH1F* _Seed_Cov_A0;
      TH1F* _Seed_Cov_A1;
      TH1F* _Seed_Cov_B0;
      TH1F* _Seed_Cov_B1;
      TH1F* _Seed_Cov_B0B1;
      TH1F* _Seed_Cov_A0A1;
     //look for coreelataions between true angle and paramets pulls
      TProfile* _A0pull_v_theta_true;
      TProfile* _B0pull_v_theta_true;
      TProfile* _A1pull_v_theta_true;
      TProfile* _B1pull_v_theta_true;
      //Seed Diagnostics:
      TH1F* _total_residualsX_init;
      TH1F* _total_pullsX_init;
      TH1F* _total_residualsY_init;
      TH1F* _total_pullsY_init;
      TH1F* _chisq_ndf_plot_final;
      TH1F* _chisq_ndf_plot_true_fit;
      TH1F* _chisq_ndf_plot_finalX;
      TH1F* _chisq_ndf_plot_finalY;
      TH1F* _chisq_ndf_plot_initX;
      TH1F* _chisq_ndf_plot_initY;
      TH1F* _chisq_ndf_plot_init;
      TH1F* _change_chisq_ndf_plot_X;
      TH1F* _change_chisq_ndf_plot_Y;
      TH1F* _change_chisq_ndf_plot_Total;
      TH1F* _total_residualsX_final;
      TH1F* _total_pullsX_final;
      TH1F* _total_residualsY_final;
      TH1F* _total_pullsY_final;
      TH1F* _FinalErrX;
      TH1F* _FinalErrY;
      TH1F* _InitErrX;
      TH1F* _InitErrY;
      TH1F* _FinalErrTot;
      TH1F* _InitErrTot;

      //Drift Parameters
      TH1F* _total_residualsX_Minuit;
      TH1F* _total_residualsY_Minuit;
      TH1F* _A0MinuitFitDiff;//difference between seed and minuit fits
      TH1F* _A1MinuitFitDiff;
      TH1F* _B0MinuitFitDiff;
      TH1F* _B1MinuitFitDiff;
      TH1F* _A0Minuit;
      TH1F* _A1Minuit;
      TH1F* _B0Minuit;
      TH1F* _B1Minuit;
      TH1F* _A0MinuitMCDiff;
      TH1F* _A1MinuitMCDiff;
      TH1F* _B0MinuitMCDiff;
      TH1F* _B1MinuitMCDiff;
      //Looking for correlations berween residuals and number of outliers
      TH2F* _DeltaA0_v_OUT;
      TH2F* _DeltaA1_v_OUT;
      TH2F* _DeltaB0_v_OUT;
      TH2F* _DeltaB1_v_OUT;

      //Drift diags:
      TH1F* _FullFitEndDOCAs;
      TH1F* _TrueDOCAs;
      TH1F* _FullFitEndTimeResiduals;
      TH1F* _StartTimeResiduals;
      TH1F* _TrueTimeResiduals;
      TH1F* _StartDOCAs;
      TH1F* _GaussianEndDOCAs;
      TH1F* _GaussianEndTimeResiduals;
      TH1F* _minuit_pullsX_final;
      TH1F* _minuit_pullsY_final;
      TH1F* _NLL;
      TH2F* _LL_v_MOM;
      TH2F* _LL_v_Nhits;
      TH2F* _AMBIG;
      //look for coreelataions between true angle and paramets pulls
      TProfile* _MinuitA0pull_v_theta_true;
      TProfile* _MinuitB0pull_v_theta_true;
      TProfile* _MinuitA1pull_v_theta_true;
      TProfile* _MinuitB1pull_v_theta_true;

      TH2F* _seedDeltaA0_v_minuitA0;
      TH2F* _seedDeltaA1_v_minuitA1;
      TH2F* _seedDeltaB0_v_minuitB0;
      TH2F* _seedDeltaB1_v_minuitB1;

      TH2F* _seedDeltaA0_v_LL;
      TH2F* _seedDeltaA1_v_LL;
      TH2F* _seedDeltaB0_v_LL;
      TH2F* _seedDeltaB1_v_LL;

     //Triggger AnalysisL:
      TH1F* _TrueMomDOCACut;
      TH1F* _TrueMomNoCuts;
      // add event id
      Int_t _evt; 

      // Event object Tags
      art::InputTag   _chtag;//combo
      art::InputTag   _tctag;//timeclusters
      art::InputTag   _costag;//Striaght tracks
      art::InputTag   _mcdigistag; //MC digis
      // Event data cache
      const ComboHitCollection* _chcol;
      const CosmicTrackSeedCollection* _coscol;
      const TimeClusterCollection* _tccol;
      const StrawDigiMCCollection* _mcdigis;
      //Numbers:
      Int_t _nsh, _nch; // # associated straw hits / event
      Int_t     _ntc; // # clusters/event
      Int_t _nhits, _nused; // # hits used
      Int_t _n_panels; // # panels
      Int_t _n_stations; // # stations
      Int_t _n_planes; // # stations
      int n_analyze =0;
      Float_t _hit_time, _hit_drift_time, _cluster_time, _dt;
	
      //Flags:
	Bool_t _StraightTrackInit, _StraightTrackConverged, _StraightTrackOK, _hitsOK;
      //offsets for MC
      SimParticleTimeOffset _toff;

      Int_t _strawid; // strawid info

      vector<ComboHitInfoMC> _chinfomc;

      bool findData(const art::Event& evt);
    };

    CosmicAnalyzer::CosmicAnalyzer(fhicl::ParameterSet const& pset) :
	art::EDAnalyzer(pset),
	_diag		(pset.get<int>("diagLevel",1)),
	_mcdiag		(pset.get<bool>("mcdiag",true)),
	_chtag		(pset.get<art::InputTag>("ComboHitCollection")),
	_tctag		(pset.get<art::InputTag>("TimeClusterCollection")),
	_costag		(pset.get<art::InputTag>("CosmicTrackSeedCollection")),
	_mcdigistag	(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")),
	_toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
        {}

    CosmicAnalyzer::~CosmicAnalyzer(){}

    void CosmicAnalyzer::beginJob() {
      // create diagnostics if requested...
      if(_diag > 0){
	art::ServiceHandle<art::TFileService> tfs;
	//Tree for detailed diagnostics
	_cosmic_analysis=tfs->make<TTree>("cosmic_analysis"," Diagnostics for Cosmic Track Fitting");

        //Create branches:
        _cosmic_analysis->Branch("evt",&_evt,"evt/I");  // add event id
        _cosmic_analysis->Branch("nhits",&_nhits,"nhits/I");
        _cosmic_analysis->Branch("StrawHitsInEvent", &_nsh, "StrawHitsInEvent/I");
	_cosmic_analysis->Branch("ComboHitsInEvent", &_nch, "ComboHitsInEvent/I");
        _cosmic_analysis->Branch("PanelsCrossedInEvent", &_n_panels, "PanelsCrossedInEvent/I");
        _cosmic_analysis->Branch("PlanesCrossedInEvent", &_n_planes, "PlanesCrossedInEvent/I");
        _cosmic_analysis->Branch("StatonsCrossedInEvent", &_n_stations, "StationsCrossedInEvent/I");
        _cosmic_analysis->Branch("TimeClustersInEvent", &_ntc, "TimeClusterInEvent/I"); 
        _cosmic_analysis->Branch("hit_time", &_hit_time, "hit_time/F");
        _cosmic_analysis->Branch("hit_drit_time", &_hit_drift_time, "hit_drift_time/F");
        _cosmic_analysis->Branch("cluster_time", &_cluster_time, "cluster_time/F");
        _cosmic_analysis->Branch("dt", &_dt, "dt/F");
        _cosmic_analysis->Branch("hitsOK",&_hitsOK,"hitsOK/B");
        _cosmic_analysis->Branch("StraightTrackInit",&_StraightTrackInit,"StraightTrackInit/B");
        _cosmic_analysis->Branch("StraightTrackOK",&_StraightTrackOK,"StraightTrackOK/B");
        _cosmic_analysis->Branch("StraightTrackConverged",&_StraightTrackConverged,"StraightTrackConverged/B");
        
	
	//--------------------------------Truth----------------------------------------//
	if(_mcdiag == true){
		_truea0 = tfs->make<TH1F>("True Track Parameter A0 ","True Track Parameter A0 " ,50,-2000, 2000);
		_truea0->GetXaxis()->SetTitle("Track Parameter A0");
		_truea0->SetStats();
		
		_trueb0 = tfs->make<TH1F>("True Track Parameter B0 "," True Track Parameter B0  " ,50,-2000, 2000);
		_trueb0->GetXaxis()->SetTitle("Track Parameter B0");
		_trueb0->SetStats();
		
		_truea1 = tfs->make<TH1F>("True Track Parameter A1  ","True Track Parameter A1 " ,50,-5,5);
		_truea1->GetXaxis()->SetTitle("Track Parameter A1");
		_truea1->SetStats();
		
		_trueb1 = tfs->make<TH1F>("True Track Parameter B1  ","True Track Parameter B1  " ,50,-5,5);
		_trueb1->GetXaxis()->SetTitle("Track Parameter B1");
		_trueb1->SetStats();

		_truea0XYZ = tfs->make<TH1F>("True Track Parameter A0 XYZ","True Track Parameter A0 XYZ " ,50,-10000, 10000);
		_truea0XYZ->GetXaxis()->SetTitle("Track Parameter A0 XYZ");
		_truea0XYZ->SetStats();
		
		_trueb0XYZ = tfs->make<TH1F>("True Track Parameter B0 XYZ","True Track Parameter B0 XYZ " ,50,-5000, 5000);
		_trueb0XYZ->GetXaxis()->SetTitle("Track Parameter B0 XYZ");
		_trueb0XYZ->SetStats();
		
		_truea1XYZ = tfs->make<TH1F>("True Track Parameter A1 XYZ ","True Track Parameter A1 XYZ" ,50,-5,5);
		_truea1XYZ->GetXaxis()->SetTitle("Track Parameter A1 XYZ");
		_truea1XYZ->SetStats();
		
		_trueb1XYZ = tfs->make<TH1F>("True Track Parameter B1 XYZ ","True Track Parameter B1 XYZ " ,50,-5,5);
		_trueb1XYZ->GetXaxis()->SetTitle("Track Parameter B1");
		_trueb1XYZ->SetStats();

		_SPtruea0XYZ = tfs->make<TH1F>("Step Point True Track Parameter A0 XYZ","Step Point True Track Parameter A0 XYZ " ,50,-10000, 10000);
		_SPtruea0XYZ->GetXaxis()->SetTitle("Track Parameter A0 XYZ");
		_SPtruea0XYZ->SetStats();
		
		_SPtrueb0XYZ = tfs->make<TH1F>("Step Point rue Track Parameter B0 XYZ","Step Point True Track Parameter B0 XYZ " ,50,-5000, 5000);
		_SPtrueb0XYZ->GetXaxis()->SetTitle("Track Parameter B0 XYZ");
		_SPtrueb0XYZ->SetStats();
		
		_SPtruea1XYZ = tfs->make<TH1F>("Step Point True Track Parameter A1 XYZ ","Step Point True Track Parameter A1 XYZ" ,50,-5,5);
		_SPtruea1XYZ->GetXaxis()->SetTitle("Track Parameter A1 XYZ");
		_SPtruea1XYZ->SetStats();
		
		_SPtrueb1XYZ = tfs->make<TH1F>("Step Point True Track Parameter B1 XYZ ","Step Point True Track Parameter B1 XYZ " ,50,-5,5);
		_SPtrueb1XYZ->GetXaxis()->SetTitle("Track Parameter B1");
		_SPtrueb1XYZ->SetStats();

		_mc_phi_angle = tfs->make<TH1F>("#phi_{true, fit}","#phi_{true, fit}" ,100,-3.141529,3.141529);
		_mc_phi_angle->GetXaxis()->SetTitle("#phi_{true, fit} [rad]");

		_mc_theta_angle = tfs->make<TH1F>("#theta_{true, fit}","#theta_{true, fit}" ,20,0,3.141529);
		_mc_theta_angle->GetXaxis()->SetTitle("#theta_{true, fit} [rad]");
	}
	_reco_phi_angle = tfs->make<TH1F>("#phi Angle of Reconstructred Track","#phi_{reco} Angle of Tracl " ,20,0,3.141529);
	_reco_phi_angle->GetXaxis()->SetTitle("#phi_{reco} [rad]");
	_reco_phi_angle->SetStats();
	

        //------------------------------------------Seed---------------------------- --------------//

	_seed_a0 = tfs->make<TH1F>("Seed Track Parameter A0 ","Seed Track Parameter A0 " ,50,-2000, 2000);
	_seed_a0->GetXaxis()->SetTitle("Seed Track Parameter A0");
	_seed_a0->SetStats();
	
	_seed_b0 = tfs->make<TH1F>("Seed  Track Parameter B0 ","Seed Track Parameter B0  " ,50,-2000, 2000);
	_seed_b0->GetXaxis()->SetTitle("Seed Track Parameter B0");
	_seed_b0->SetStats();
	
	_seed_a1 = tfs->make<TH1F>("Seed Track Parameter A1  ","Seed vTrack Parameter A1 " ,50,-1,1);
	_seed_a1->GetXaxis()->SetTitle("Seed Track Parameter A1");
	_seed_a1->SetStats();
	
	_seed_b1 = tfs->make<TH1F>("Seed Track Parameter B1  ","Seed Track Parameter B1  " ,50,-1,1);
	_seed_b1->GetXaxis()->SetTitle("Seed Track Parameter B1");
	_seed_b1->SetStats();

	_seed_a0XYZ = tfs->make<TH1F>("Seed Track Parameter A0 XYZ","Seed Track Parameter A0 XYZ " ,50,-10000, 10000);
	_seed_a0XYZ->GetXaxis()->SetTitle("Seed Track Parameter A0 XYZ");
	_seed_a0XYZ->SetStats();
	
	_seed_b0XYZ = tfs->make<TH1F>("Seedf Track Parameter B0 XYZ","Seed Track Parameter B0 XYZ " ,50,-5000, 5000);
	_seed_b0XYZ->GetXaxis()->SetTitle("Seed Track Parameter B0 XYZ");
	_seed_b0XYZ->SetStats();
	
	_seed_a1XYZ = tfs->make<TH1F>("Seedf Track Parameter A1 XYZ ","Seed Track Parameter A1 XYZ" ,50,-5,5);
	_seed_a1XYZ->GetXaxis()->SetTitle("Seed Track Parameter A1 XYZ");
	_seed_a1XYZ->SetStats();
	
	_seed_b1XYZ = tfs->make<TH1F>("Seed Track Parameter B1 XYZ ","Seed Track Parameter B1 XYZ " ,50,-5,5);
	_seed_b1XYZ->GetXaxis()->SetTitle("Seed Track Parameter B1");
	_seed_b1XYZ->SetStats();

	_niters  = tfs->make<TH1F>("Seed Number of Iterations Until Seed Fit Converged  ","Number of Iterations Until Seed Fit Converged " ,100,0, 100);
	_niters->GetXaxis()->SetTitle("Seed Number of Iterations Until Seed Fit Converged");
	_niters->SetStats();
	
	_seedDeltaA0_v_N  = tfs->make<TH2F>("#Delta A^{seed-MC}_{0} v #N_{hits}", "#N_{hits}", 50, -100,100, 50, 0,50);
	_seedDeltaA0_v_N->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{0}");
	_seedDeltaA0_v_N->GetYaxis()->SetTitle("#N_{hits}");

	_seedDeltaA1_v_N  = tfs->make<TH2F>("#Delta A^{seed-MC}_{1} v #N_{hits}", "#N_{hits}", 50, -0.5, 0.5, 50, 0, 50);
	_seedDeltaA1_v_N->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{1}");
	_seedDeltaA1_v_N->GetYaxis()->SetTitle("#N_{hits}");

	_seedDeltaB0_v_N = tfs->make<TH2F>("#Delta B^{seed-MC}_{0} v #N_{hits}", "#N_{hits}", 50, -100,100, 50, 0,50);
	_seedDeltaB0_v_N->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{0}");
	_seedDeltaB0_v_N->GetYaxis()->SetTitle("#N_{hits}");

	_seedDeltaB1_v_N  = tfs->make<TH2F>("#Delta B^{seed-MC}_{1} v #N_{hits}", "#N_{hits}", 50, -0.5, 0.5, 50, 0,50);
	_seedDeltaB1_v_N->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{1}");
	_seedDeltaB1_v_N->GetYaxis()->SetTitle("#N_{hits}");

	_seedDeltaA0_v_MOM  = tfs->make<TH2F>("#Delta A^{seed-MC}_{0} v #p_{true}", "#p_{true}", 50, -100,100, 50, 0,10);
	_seedDeltaA0_v_MOM->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{0}");
	_seedDeltaA0_v_MOM->GetYaxis()->SetTitle("#p_{true}");

	_seedDeltaA1_v_MOM  = tfs->make<TH2F>("#Delta A^{seed-MC}_{1} v #p_{true}", "#p_{true}", 50, -0.5, 0.5, 50, 0, 10);
	_seedDeltaA1_v_MOM->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{1}");
	_seedDeltaA1_v_MOM->GetYaxis()->SetTitle("#p_{true}");

	_seedDeltaB0_v_MOM = tfs->make<TH2F>("#Delta B^{seed-MC}_{0} v #p_{true}", "#p_{true}", 50, -100,100, 50, 0,10);
	_seedDeltaB0_v_MOM->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{0}");
	_seedDeltaB0_v_MOM->GetYaxis()->SetTitle("#p_{true}");

	_seedDeltaB1_v_MOM  = tfs->make<TH2F>("#Delta B^{seed-MC}_{1} v #p_{true}", "#p_{true}", 50, -0.5, 0.5, 50, 0,10);
	_seedDeltaB1_v_MOM->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{1}");
	_seedDeltaB1_v_MOM->GetYaxis()->SetTitle("#p_{true}");

	_DeltaA0_v_OUT = tfs->make<TH2F>("#Delta A_{0} v Outliers", "#Outliers", 50, -100,100, 20, 0,20);
	_DeltaA0_v_OUT->GetXaxis()->SetTitle("#Delta A_{0}");
	_DeltaA0_v_OUT->GetYaxis()->SetTitle("Outliers");

	_DeltaA1_v_OUT  = tfs->make<TH2F>("#Delta A_{1} v Outliers", "#Outliers", 50, -0.5, 0.5, 20, 0, 20);
	_DeltaA1_v_OUT->GetXaxis()->SetTitle("#Delta A_{1}");
	_DeltaA1_v_OUT->GetYaxis()->SetTitle("Oiutliers");

	_DeltaB0_v_OUT = tfs->make<TH2F>("#Delta B_{0} v Outliers", "Outliers", 50, -100,100, 20, 0,20);
	_DeltaB0_v_OUT->GetXaxis()->SetTitle("#Delta B_{0}");
	_DeltaB0_v_OUT->GetYaxis()->SetTitle("Outliers");

	_DeltaB1_v_OUT  = tfs->make<TH2F>("#Delta B_{1} v Outliers", "Outliers", 50, -0.5, 0.5, 20, 0,20);
	_DeltaB1_v_OUT->GetXaxis()->SetTitle("#Delta B_{1}");
	_DeltaB1_v_OUT->GetYaxis()->SetTitle("Outliers");

	_seedA0_v_seedA1 = tfs->make<TH2F>("Seed A_{0} v A_{1}", "Seed A_{0} v A_{1}", 50, -0.2, 0.2, 50, -1500, 1500);
	_seedA0_v_seedA1->GetXaxis()->SetTitle("A_{1}");
	_seedA0_v_seedA1->GetYaxis()->SetTitle("A_{0}");

	_seedB0_v_seedB1 = tfs->make<TH2F>("Seed B_{0} v B_{1}", "Seed B_{0} v B_{1}", 50, -0.2, 0.2, 50, -1500, 1500);
	_seedB0_v_seedB1->GetXaxis()->SetTitle("B_{1}");
	_seedB0_v_seedB1->GetYaxis()->SetTitle("B_{0}");
	
	_seedA0_v_seedB1 = tfs->make<TH2F>("Seed A_{0} v B_{1}", "Seed A_{0} v B_{1}", 50, -0.2, 0.2, 50, -1500, 1500);
	_seedA0_v_seedB1->GetXaxis()->SetTitle("B_{1}");
	_seedA0_v_seedB1->GetYaxis()->SetTitle("A_{0");
	
	_seedA0_v_seedB0 = tfs->make<TH2F>("Seed B_{0} v A_{0}", "Seed B_{0} v A_{0}", 50, -1500,1500, 50, -2000, 2000);
	_seedA0_v_seedB0->GetXaxis()->SetTitle("B_{0}");
	_seedA0_v_seedB0->GetYaxis()->SetTitle("A_{0} ");
	
	_seedB1_v_seedA1 = tfs->make<TH2F>("Seed B_{1} v A_{1}", "Seed B_{1} v A_{1}", 100, -0.2,0.2, 100, -0.2,0.2);
	_seedB1_v_seedA1->GetXaxis()->SetTitle("B_{1}");
	_seedB1_v_seedA1->GetYaxis()->SetTitle("A_{1} ");
	
	_A0SeedMCDiff = tfs->make<TH1F>("A^{seed}_{0}- A^{MC}_{0} ","A^{seed}_{0}- A^{MC}_{0}" ,100,-100, 100);
	_A0SeedMCDiff->GetXaxis()->SetTitle("#Delta #A_{0}");
	_A0SeedMCDiff->SetStats();
	
	_A1SeedMCDiff = tfs->make<TH1F>("A^{seed}_{1} - A^{MC}_{1}","A^{seed}_{1} - A^{MC}_{1}" ,100,-0.5 , 0.5);
	_A1SeedMCDiff->GetXaxis()->SetTitle("#Delta A_{1}");
	_A1SeedMCDiff->SetStats();
	
	_B0SeedMCDiff = tfs->make<TH1F>("B^{seed}_{0} - B^{seed}_{0}" ,"B^{seed}_{0} - B^{MC}_{0}",100,-100, 100);
	_B0SeedMCDiff->GetXaxis()->SetTitle("#Delta B^{seed}_{0}");
	_B0SeedMCDiff->SetStats();
	
	_B1SeedMCDiff = tfs->make<TH1F>("B^{seed}_{1} - B^{seed}_{1}", "B^{seed}_{1} - B^{seed}_{1}" ,100,-0.5, 0.5);
	_B1SeedMCDiff->GetXaxis()->SetTitle("#Delta B^{seed}_{1}");
	_B1SeedMCDiff->SetStats();
	
	_chiX_v_true_trackposX = tfs->make<TH2F>("Seed #chi^{2}/ndf v A_{0,true}", "Seed #chi^{2}/ndf v A_{0,true}", 100, -2000,2000, 100, 0,5);
	_chiX_v_true_trackposX->GetXaxis()->SetTitle("A_{0}");
	_chiX_v_true_trackposX->GetYaxis()->SetTitle("#chi{2}/ndf");
	
	_chiY_v_true_trackposY = tfs->make<TH2F>("Seed #chi^{2}/ndf v B_{0,true}", "Seed #chi^{2}/ndf v B_{0,true}", 100, -2000,2000, 100, 0,5);
	_chiY_v_true_trackposY->GetXaxis()->SetTitle("B_{0}");
	_chiY_v_true_trackposY->GetYaxis()->SetTitle("#chi{2}/ndf");
	
	_chi_v_true_theta = tfs->make<TH2F>("Seed #chi^{2}/ndf v #theta_{true}", "Seed #chi^{2}/ndf v #theta_{true}", 100, 0, 3.141529, 100, 0,2);
	_chi_v_true_theta->GetXaxis()->SetTitle("B_{0}");
	_chi_v_true_theta->GetYaxis()->SetTitle("#chi{2}/ndf");
	
	_Seed_Cov_A0 = tfs->make<TH1F>("#sigma_{A_{0}}", "#sigma_{A_{0}}" ,100,0, 100);
	_Seed_Cov_A0->GetXaxis()->SetTitle("#sigma_{A_{0}}");

	_Seed_Cov_A0A1 = tfs->make<TH1F>("Seed #sigma_{A_{0}A_{1}}", "Seed #sigma_{A_{0}A_{1}}" ,100,0, 100);
	_Seed_Cov_A0A1->GetXaxis()->SetTitle("#sigma_{A_{0}A_{1}}");

	_Seed_Cov_A1 = tfs->make<TH1F>("Seed #sigma_{A_{1}}", "Seed #sigma_{A_{1}}" ,100,0, 0.5);
	_Seed_Cov_A1->GetXaxis()->SetTitle("#sigma_{A_{1}}");
	
	_Seed_Cov_B0 = tfs->make<TH1F>("Seed #sigma_{B_{0}}", "Seed #sigma_{B_{0}}" ,100,0, 100);
	_Seed_Cov_B0->GetXaxis()->SetTitle("#sigma_{B_{0}}");
		
	_Seed_Cov_B1 = tfs->make<TH1F>("Seed #sigma_{B_{1}}", "Seed #sigma_{B_{1}}" ,100,0, 0.5);
	_Seed_Cov_B1->GetXaxis()->SetTitle("#sigma_{B_{1}}");

        _Seed_Cov_B0B1 = tfs->make<TH1F>("Seed #sigma_{B_{0}B_{1}}", "Seed #sigma_{B_{0}B_{1}}" ,100,0, 100);
	_Seed_Cov_B0B1->GetXaxis()->SetTitle("#sigma_{B_{0}B_{1}}");

	_A0pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull A_{0} v #theta_{true,fit}  ","Seed Parameter Pull A_{0} v #theta_{true,fit} ",100,0,3.141529,-200,200);
	_A0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_A0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{0}");
	
	_B0pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull B_{0} v #theta_{true,fit}  ","Seed Parameter Pull B_{0} v #theta_{true,fit} ",100,0,3.141529,-200,200);
	_B0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_B0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{0}");
	
	_A1pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull A_{1} v #theta_{true,fit}  ","Seed Parameter Pull A_{1} v #theta_{true,fit} ",100,0,3.141529,-10,10);
	_A1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_A1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{1}");
	
	_B1pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull B_{1} v #theta_{true,fit}  ","Seed Parameter Pull B_{1} v #theta_{true,fit} ",100,0,3.141529,-10,10);
	_B1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_B1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{1}");
	
	_chisq_ndf_plot_init = tfs->make<TH1F>("Seed Init ChiSq/NDF","Seed Init ChiSq/NDF" ,100,0, 10);
	_chisq_ndf_plot_init->GetXaxis()->SetTitle("Init. #Chi^{2}/N");
	
	_chisq_ndf_plot_final = tfs->make<TH1F>("Seed Final ChiSq/NDF","Seed Final ChiSq/NDF" ,100,0,10);
	_chisq_ndf_plot_final->GetXaxis()->SetTitle("Final #Chi^{2}/N");
	
	_chisq_ndf_plot_finalX = tfs->make<TH1F>("Seed Final ChiSq/NDF X''","Seed Final ChiSq/NDF X''" ,100,0,10);
	_chisq_ndf_plot_finalX->GetXaxis()->SetTitle("Final X '' #Chi^{2}/N");
	
        _chisq_ndf_plot_finalY = tfs->make<TH1F>("Seed Final ChiSq/NDF Y''","Seed Final ChiSq/NDF Y''" ,100,0,10);
	_chisq_ndf_plot_finalY->GetXaxis()->SetTitle("Final Y '' #Chi^{2}/N");
	
	_chisq_ndf_plot_initX = tfs->make<TH1F>("Seed Init Chi2/NDF  X''","Seed Init Chi2/NDF X''" ,100,0,10);
	_chisq_ndf_plot_initX->GetXaxis()->SetTitle("Initial X '' #Chi^{2}/N");
	
        _chisq_ndf_plot_initY = tfs->make<TH1F>("Seed Init Chi2/NDF Y''","initial chisq_ndf_plot Y''" ,100,0,10);
	_chisq_ndf_plot_initY->GetXaxis()->SetTitle("Initial Y '' #Chi^{2}/N");
	
	_change_chisq_ndf_plot_X = tfs->make<TH1F>("Seed Total Change in Chisq/ ndf X''","Total Change in Chisq/ ndf X''" ,50,-10,20);
	_change_chisq_ndf_plot_X->GetXaxis()->SetTitle("Initial X '' #Chi^{2}/N");
	
	_change_chisq_ndf_plot_Y = tfs->make<TH1F>("Seed Total Change in Chisq/ ndf Y''","Seed Total Change in Chisq/ ndf Y''" ,50,-10,20);
	_change_chisq_ndf_plot_Y->GetXaxis()->SetTitle("Initial Y '' #Chi^{2}/N");
	
	_change_chisq_ndf_plot_Total = tfs->make<TH1F>("Seed Total Change in #Chi^{2}/ ndf","Seed Total Change in #Chi^{2}/ ndf" ,50,-5,20);
	_change_chisq_ndf_plot_Total->GetXaxis()->SetTitle("Total Change in #Chi^{2}/ndf");
	
        _total_residualsX_init = tfs->make<TH1F>("Seedf Initial Residuals X'' ","Initial Residuals X'' " ,50,-1000,1000);
	_total_residualsX_init->GetXaxis()->SetTitle("Initial Residual X'' [mm]");
	_total_residualsX_init->SetStats();
	
	_total_residualsY_init = tfs->make<TH1F>("Seed Initial Residuals Y''","Seed Initial Residuals Y'' " ,50,-1000,1000);
	_total_residualsY_init->GetXaxis()->SetTitle("Initial Residual Y'' [mm]");
	_total_residualsY_init->SetStats();

	_total_pullsX_init = tfs->make<TH1F>("Seed Initial PullSeed  X","Initial Pull X''" ,200,-50, 50);
	_total_pullsX_init->GetXaxis()->SetTitle("Pull X");
	_total_pullsX_init->SetStats();
	
	_total_pullsY_init = tfs->make<TH1F>("Seed Initial Pull Y''","Seed Initial Pull Y''" ,200,-50, 50);
	_total_pullsY_init->GetXaxis()->SetTitle("Initial Pull Y");
	_total_pullsY_init->SetStats();
       
	_total_residualsX_final = tfs->make<TH1F>("Seed Final Residuals X''  ","Seed Final Residuals X''#chi^{2}/dof all  " ,100,-200,200);
	_total_residualsX_final->GetXaxis()->SetTitle("Residual X'' [mm]");
	_total_residualsX_final->SetStats();
	
	_total_residualsY_final = tfs->make<TH1F>("Seed Final Residuals Y'' ","Seed Final Residuals Y'' #chi^{2}/dof all " ,100,-200,200);
	_total_residualsY_final->GetXaxis()->SetTitle("Final Residual Y'' [mm]");
	_total_residualsY_final->SetStats();
	
	_total_pullsX_final = tfs->make<TH1F>("Seed Final Pull X'' #chi^{2}/dof alSeed l ","Seed Final Pull X''#chi^{2}/dof all " ,50,-2, 2);
	_total_pullsX_final->GetXaxis()->SetTitle("Final Pull X''");
	
	
	_total_pullsY_final = tfs->make<TH1F>("Seed Final Pull Y''#chi^{2}/dof all ","Seed Final Pull Y''#chi^{2}/dof all " ,50,-2, 2);
	_total_pullsY_final->GetXaxis()->SetTitle("Final Pull Y''");
	
	_FinalErrX = tfs->make<TH1F>("Seed Final Errors in X''","Seed Final Errors in X'' " ,50,0, 100);
	_FinalErrX->GetXaxis()->SetTitle("#sigma_{X''} [mm]");
	_FinalErrX->SetStats();
	
	_FinalErrY = tfs->make<TH1F>("Seed Final Errors in Y''","Seed Final Errors in Y'' " ,50,0, 100);
	_FinalErrY->GetXaxis()->SetTitle("#sigma_{Y''} [mm]");
	_FinalErrY->SetStats();
	
	_InitErrX = tfs->make<TH1F>("Seed Initial Errors in X''","Seed Initial Errors in X'' " ,50,0, 100);
	_InitErrX->GetXaxis()->SetTitle("#sigma_{X''} [mm]");
	_InitErrX->SetStats();
	
	_InitErrY = tfs->make<TH1F>("Seed Initial Errors in Y''","Seed Initial Errors in Y'' " ,50,0, 100);
	_InitErrY->GetXaxis()->SetTitle("#sigma_{Y''} [mm]");
	_InitErrY->SetStats();
	
	_FinalErrTot = tfs->make<TH1F>("Seed Final Errors Total","Seed Final Errors Total " ,50,0, 100);
	_FinalErrTot->GetXaxis()->SetTitle("#sigma_{Y''} [mm]");
	_FinalErrTot->SetStats();
	
	_InitErrTot = tfs->make<TH1F>("Seed Initial Errors Total","Seed Initial Errors Total " ,50,0, 100);
	_InitErrTot->GetXaxis()->SetTitle("#sigma_{X''} [mm]");
	_InitErrTot->SetStats();

	//------------------------Straw Fit-------------------------------//
	
	_total_residualsX_Minuit = tfs->make<TH1F>("Minuit Residuals X''  ","Minuit Final Residuals X'' " ,100,-200,200);
	_total_residualsX_Minuit->GetXaxis()->SetTitle("Minuit Res. X'' [mm]");
	_total_residualsX_Minuit->SetStats();
	
	_total_residualsY_Minuit = tfs->make<TH1F>("Minuit Residuals Y'' ","Minuit Residuals Y''  " ,100,-200,200);
	_total_residualsY_Minuit->GetXaxis()->SetTitle("Minuit Res. Y'' [mm]");
	_total_residualsY_Minuit->SetStats();
	
	_FullFitEndDOCAs = tfs->make<TH1F>("Final DOCAs ","Final DOCAs" ,50,0,10);
	_FullFitEndDOCAs->GetXaxis()->SetTitle("Final DOCAs [mm]");
	_FullFitEndDOCAs->SetStats();

	_GaussianEndDOCAs = tfs->make<TH1F>("Gaussian DOCAs ","Gaussian DOCAs" ,50,0,10);
	_GaussianEndDOCAs->GetXaxis()->SetTitle("Gaussian DOCAs [mm]");
	_GaussianEndDOCAs->SetStats();
	
	_StartDOCAs = tfs->make<TH1F>("Initial DOCAs ","InitialDOCAs" ,50, 0,10);
	_StartDOCAs->GetXaxis()->SetTitle("Initial DOCAs [mm]");
	_StartDOCAs->SetStats();

	_TrueDOCAs = tfs->make<TH1F>("True DOCAs ","TrueDOCAs" ,50, 0,10);
	_TrueDOCAs->GetXaxis()->SetTitle("True DOCAs [mm]");
	_TrueDOCAs->SetStats();

        _FullFitEndTimeResiduals = tfs->make<TH1F>("Final Time Res ","Final Time Resid" ,50,0,500);
	_FullFitEndTimeResiduals->GetXaxis()->SetTitle("Final Time Res [ns]");
	_FullFitEndTimeResiduals->SetStats();

	_GaussianEndTimeResiduals = tfs->make<TH1F>("Gaussian Time Res ","Gaussian Time Resid" ,50,0,500);
	_GaussianEndTimeResiduals->GetXaxis()->SetTitle("Gaussian Time Res [ns]");
	_GaussianEndTimeResiduals->SetStats();

	_StartTimeResiduals = tfs->make<TH1F>("Start Time Res ","Start Time Resid" ,50,0,500);
	_StartTimeResiduals->GetXaxis()->SetTitle("Start Time Res [ns]");
	_StartTimeResiduals->SetStats();

	_TrueTimeResiduals = tfs->make<TH1F>("True Time Res ","True Time Resid" ,50,0,500);
	_TrueTimeResiduals->GetXaxis()->SetTitle("True Time Res [ns]");
	_TrueTimeResiduals->SetStats();
	
	_A0MinuitFitDiff = tfs->make<TH1F>(" A0 Minuit Diff "," A0 Minuit Diff" ,100,-100,100);
	_A0MinuitFitDiff->GetXaxis()->SetTitle("#Delta A0 [mm]");
	_A0MinuitFitDiff->SetStats();
	
	_A1MinuitFitDiff = tfs->make<TH1F>(" A1 Minuit Diff "," A1 Minuit Diff" ,100,-3,3);
	_A1MinuitFitDiff->GetXaxis()->SetTitle("#Delta A1 ");
	_A1MinuitFitDiff->SetStats();
	
	_B0MinuitFitDiff = tfs->make<TH1F>(" B0 Minuit Diff "," B0 Minuit Diff" ,100,-100,100);
	_B0MinuitFitDiff->GetXaxis()->SetTitle("#Delta B0 [mm] ");
	_B0MinuitFitDiff->SetStats();
	
	_B1MinuitFitDiff = tfs->make<TH1F>(" B1 Minuit Diff "," B1 Minuit Diff" ,100,-3,3);
	_B1MinuitFitDiff->GetXaxis()->SetTitle("#Delta B1 ");
	_B1MinuitFitDiff->SetStats();
	
	_minuit_pullsX_final = tfs->make<TH1F>("Miuit Pull X'' ","Minuit Pull X'' " ,50,-2, 2);
	_minuit_pullsX_final->GetXaxis()->SetTitle("Minuit Pull X''");
	
	_minuit_pullsY_final = tfs->make<TH1F>("Minuit Pull Y'' ","Minuit Pull Y'' " ,50,-2, 2);
	_minuit_pullsY_final->GetXaxis()->SetTitle("Minuit Pull Y''");
	
	_NLL = tfs->make<TH1F>(" Log(L) "," Log(L)" ,200,0,1000);
	_NLL->GetXaxis()->SetTitle("Log(L) ");
	_NLL->SetStats();
	
	_A0Minuit = tfs->make<TH1F>(" A0_{Minuit} "," A0_{Minuit} " ,50,-10000,10000);
	_A0Minuit->GetXaxis()->SetTitle(" [mm]");
	_A0Minuit->SetStats();
	
	_A1Minuit = tfs->make<TH1F>(" A1_{Minuit}  "," A1_{Minuit}" ,50,-5,5);
	_A1Minuit->GetXaxis()->SetTitle("A1 ");
	_A1Minuit->SetStats();
	
	_B0Minuit = tfs->make<TH1F>(" B0_{Minuit}"," B0_{Minuit} " ,50,-5000,5000);
	_B0Minuit->GetXaxis()->SetTitle(" [mm]");
	_B0Minuit->SetStats();
	
	_B1Minuit = tfs->make<TH1F>(" B1_{Minuit}  "," B1_{Minuit} " ,50,-5,5);
	_B1Minuit->GetXaxis()->SetTitle(" B1 ");
	_B1Minuit->SetStats();
	
	_A0MinuitMCDiff = tfs->make<TH1F>(" A0_{Minuit} - A0_{MC} "," A0_{Minuit} - A0_{Minuit}" ,100,-100,100);
	_A0MinuitMCDiff->GetXaxis()->SetTitle("#Delta A0 [mm]");
	_A0MinuitMCDiff->SetStats();
	
	_A1MinuitMCDiff = tfs->make<TH1F>(" A1_{Minuit} - A1_{MC} "," A1_{Minuit} - A1_{MC}" ,100,-0.5,0.5);
	_A1MinuitMCDiff->GetXaxis()->SetTitle("#Delta A1 ");
	_A1MinuitMCDiff->SetStats();
	
	_B0MinuitMCDiff = tfs->make<TH1F>(" B0_{Minuit} - B0_{MC} "," B0_{Minuit} - B0_{MC}" ,100,-100,100);
	_B0MinuitMCDiff->GetXaxis()->SetTitle("#Delta B0 [mm]");
	_B0MinuitMCDiff->SetStats();
	
	_B1MinuitMCDiff = tfs->make<TH1F>(" B1_{Minuit} - B1_{MC}  "," B1_{Minuit} - B1_{MC}" ,100,-0.5,0.5);
	_B1MinuitMCDiff->GetXaxis()->SetTitle("#Delta B1 ");
	_B1MinuitMCDiff->SetStats();

        _seedDeltaA0_v_minuitA0  = tfs->make<TH2F>("#Delta A^{seed-MC}_{0} v A^{straw}_{0}", "A^{straw}_{0}", 50, -100,100, 50, -5000, 5000);
	_seedDeltaA0_v_minuitA0->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{0}");
	_seedDeltaA0_v_minuitA0->GetYaxis()->SetTitle("A^{straw}_{0}");

	_seedDeltaA1_v_minuitA1  = tfs->make<TH2F>("#Delta A^{seed-MC}_{1} v A^{straw}_{1}", "A^{straw}_{1}", 50, -0.5, 0.5, 50, -5, 5);
	_seedDeltaA1_v_minuitA1->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{1}");
	_seedDeltaA1_v_minuitA1->GetYaxis()->SetTitle("A^{straw}_{1}");

	_seedDeltaB0_v_minuitB0  = tfs->make<TH2F>("#Delta B^{seed-MC}_{0} v B^{straw}_{0}", "B^{straw}_{0}", 50, -100,100, 50, -5000, 5000);
	_seedDeltaB0_v_minuitB0->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{0}");
	_seedDeltaB0_v_minuitB0->GetYaxis()->SetTitle("B^{straw}_{0}");

	_seedDeltaB1_v_minuitB1  = tfs->make<TH2F>("#Delta B^{seed-MC}_{1} v B^{straw}_{1}", "B^{straw}_{1}", 50, -0.5, 0.5, 50, -5, 5);
	_seedDeltaB1_v_minuitB1->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{1}");
	_seedDeltaB1_v_minuitB1->GetYaxis()->SetTitle("B^{straw}_{1}");

	_seedDeltaA0_v_LL  = tfs->make<TH2F>("#Delta A^{seed-MC}_{0} v LogL", "End Log L", 50, -100,100, 50, 0,200);
	_seedDeltaA0_v_LL->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{0}");
	_seedDeltaA0_v_LL->GetYaxis()->SetTitle("End LogL");

	_seedDeltaA1_v_LL  = tfs->make<TH2F>("#Delta A^{seed-MC}_{1} v LogL", "End LogL", 50, -0.5, 0.5, 50, 0, 200);
	_seedDeltaA1_v_LL->GetXaxis()->SetTitle("#Delta A^{seed-MC}_{1}");
	_seedDeltaA1_v_LL->GetYaxis()->SetTitle("End LogL");

	_seedDeltaB0_v_LL = tfs->make<TH2F>("#Delta B^{seed-MC}_{0} v LogL", "End LogL", 50, -100,100, 50, 0,200);
	_seedDeltaB0_v_LL->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{0}");
	_seedDeltaB0_v_LL->GetYaxis()->SetTitle("B^{straw}_{0}");

	_seedDeltaB1_v_LL  = tfs->make<TH2F>("#Delta B^{seed-MC}_{1} v LogL", "End LogL", 50, -0.5, 0.5, 50, 0,200);
	_seedDeltaB1_v_LL->GetXaxis()->SetTitle("#Delta B^{seed-MC}_{1}");
	_seedDeltaB1_v_LL->GetYaxis()->SetTitle("End LogL");

	_MinuitA0pull_v_theta_true= tfs->make<TProfile>("Minuit Fit Parameter Pull A_{0} v #theta_{true,fit}  ","Minuit Fit Parameter Pull A_{0} v #theta_{true,fit} ",100,0,3.141529,-200,200);
	_MinuitA0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_MinuitA0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{0}");
	
	_MinuitB0pull_v_theta_true= tfs->make<TProfile>("Minuit Fit Parameter Pull B_{0} v #theta_{true,fit}  ","Minuit Fit Parameter Pull B_{0} v #theta_{true,fit} ",100,0,3.141529,-200,200);
	_MinuitB0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_MinuitB0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{0}");
	
	_MinuitA1pull_v_theta_true= tfs->make<TProfile>("Minuit Fit Parameter Pull A_{1} v #theta_{true,fit}  ","Minuit Fit Parameter Pull A_{1} v #theta_{true,fit} ",100,0,3.141529,-10,10);
	_MinuitA1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_MinuitA1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{1}");
	
	_MinuitB1pull_v_theta_true= tfs->make<TProfile>("Minuit Fit Parameter Pull B_{1} v #theta_{true,fit}  ","Minuit Fit Parameter Pull B_{1} v #theta_{true,fit} ",100,0,3.141529,-10,10);
	_MinuitB1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_MinuitB1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{1}");

	_LL_v_MOM  = tfs->make<TH2F>("Log(L) v p_{true}", "p_{true}", 50, 0,1000, 50, 0,10);
	_LL_v_MOM->GetXaxis()->SetTitle("Log(L)");
	_LL_v_MOM->GetYaxis()->SetTitle("p_{true}");

	_LL_v_Nhits = tfs->make<TH2F>("Log(L) v N_{hits}", "N_{Hits}", 100, 0, 1000, 30, 0,30);
	_LL_v_Nhits->GetXaxis()->SetTitle("Log(L)");
	_LL_v_Nhits->GetYaxis()->SetTitle("#N_{hits}");

	_AMBIG = tfs->make<TH2F>("True v Reco", "True v Reco", 2,-2,2,2,-2,2);
	_AMBIG->GetXaxis()->SetTitle("True");
	_AMBIG->GetYaxis()->SetTitle("Reco");
	if(_mcdiag==true){
		_TrueMomDOCACut = tfs->make<TH1F>("True Mom Tot After DOCA Cut", "True Mom After DOCA cut", 100,0,100);
		_TrueMomDOCACut->GetXaxis()->SetTitle("#p_{true} [MeV/c]");
		
		_TrueMomNoCuts = tfs->make<TH1F>("True Mom Tot No Cut", "True Mom No cut", 100,0,100);
		_TrueMomNoCuts->GetXaxis()->SetTitle("#p_{true} [MeV/c]");
	}
        }
      }
      void CosmicAnalyzer::analyze(const art::Event& event) {
       
        _evt = event.id().event();  // add event id
        if(!findData(event)) // find data
      		throw cet::exception("RECO")<<"No Time Clusters in event"<< endl; 
       
        //find time clusters:
    	_ntc = _tccol->size();
        _nch = _chcol->size();
        
        for(size_t itc=0; itc<_tccol->size();++itc){
		TimeCluster tc = (*_tccol)[itc];
        	_cluster_time =  tc._t0._t0;
                //terr  = tc._t0._t0err;
	}
	
        //loop over tracks
        for(size_t ist = 0;ist < _coscol->size(); ++ist){
        	n_analyze+=1;
        	
        	CosmicTrackSeed sts =(*_coscol)[ist];
		CosmicTrack st = sts._track;
		CosmicTrkFitFlag const& status = sts._status;
        	if (!status.hasAllProperties(CosmicTrkFitFlag::StraightTrackOK) ){continue;}
		if(st.converged == false or st.minuit_converged  == false) { continue;}
		std::vector<int> panels, planes, stations;


                if(_mcdiag ){
		        _chiX_v_true_trackposX->Fill(st.FitParams.A0, st.Diag.FinalChiX);
		        _chiY_v_true_trackposY->Fill(st.SeedTrueParams.B0, st.Diag.FinalChiY);
               
                	_mc_phi_angle->Fill(st.get_true_phi());
	                _mc_theta_angle->Fill(st.get_true_theta());                  
	                _chi_v_true_theta->Fill(st.get_true_theta(),st.get_true_chisq());
			
			_truea1->Fill(st.SeedTrueParams.A1);
		        _trueb1->Fill(st.SeedTrueParams.B1);
		        _truea0->Fill(st.SeedTrueParams.A0);
		        _trueb0->Fill(st.SeedTrueParams.B0);
	               
		        _truea1XYZ->Fill(st.TrueFitEquation.Dir.X());
		        _trueb1XYZ->Fill(st.TrueFitEquation.Dir.Y());
		        _truea0XYZ->Fill(st.TrueFitEquation.Pos.X());
		        _trueb0XYZ->Fill(st.TrueFitEquation.Pos.Y());

			_SPtruea1XYZ->Fill(st.TrueTrueTrackDirection.X());
		        _SPtrueb1XYZ->Fill(st.TrueTrueTrackDirection.Y());
		        _SPtruea0XYZ->Fill(st.TrueTrueTrackPosition.X());
		        _SPtrueb0XYZ->Fill(st.TrueTrueTrackPosition.Y());

			
		        _MinuitA0pull_v_theta_true->Fill(st.get_true_theta(), (st.TrueFitEquation.Pos.X()- st.MinuitFitParams.A0));
		        _MinuitA1pull_v_theta_true->Fill(st.get_true_theta(), (st.TrueFitEquation.Dir.X()- st.MinuitFitParams.A1 ));
		        _MinuitB0pull_v_theta_true->Fill(st.get_true_theta(), (st.TrueFitEquation.Pos.Y()- st.MinuitFitParams.B0));
		        _MinuitB1pull_v_theta_true->Fill(st.get_true_theta(), (st.TrueFitEquation.Dir.Y()- st.MinuitFitParams.B1));
			if(st.MinuitFitParams.Covarience.sigA0 !=0){
				_A0MinuitMCDiff->Fill((st.TrueFitEquation.Pos.X()- st.MinuitFitParams.A0));///(st.MinuitFitParams.Covarience.sigA0));
			}if(st.MinuitFitParams.Covarience.sigA1 !=0){
				_A1MinuitMCDiff->Fill((st.TrueFitEquation.Dir.X() - st.MinuitFitParams.A1));///(st.MinuitFitParams.Covarience.sigA1));
		 	}if(st.MinuitFitParams.Covarience.sigB0 !=0){
				_B0MinuitMCDiff->Fill((st.TrueFitEquation.Pos.Y()- st.MinuitFitParams.B0));///(st.MinuitFitParams.Covarience.sigB0));
				
			}if(st.MinuitFitParams.Covarience.sigB1 !=0){
				_B1MinuitMCDiff->Fill((st.TrueFitEquation.Dir.Y() - st.MinuitFitParams.B1));//(st.MinuitFitParams.Covarience.sigB1));
		       }
                }
		
		_reco_phi_angle->Fill(st.get_fit_phi()); 
		_niters->Fill(st.get_iter()); 
		
                _seed_a1->Fill(st.FitParams.A1);
                _seed_b1->Fill(st.FitParams.B1);
                _seed_a0->Fill(st.FitParams.A0);
                _seed_b0->Fill(st.FitParams.B0);

		_seed_a1XYZ->Fill(st.FitEquationXYZ.Dir.X());
                _seed_b1XYZ->Fill(st.FitEquationXYZ.Dir.Y());
                _seed_a0XYZ->Fill(st.FitEquationXYZ.Pos.X());
                _seed_b0XYZ->Fill(st.FitEquationXYZ.Pos.Y());

                _Seed_Cov_A0->Fill(st.FitParams.Covarience.sigA0);
                _Seed_Cov_A1->Fill(st.FitParams.Covarience.sigA1);
                _Seed_Cov_B0->Fill(st.FitParams.Covarience.sigB0);
                _Seed_Cov_B1->Fill(st.FitParams.Covarience.sigB1);
                _Seed_Cov_B0B1->Fill(st.FitParams.Covarience.sigB1B0);
		_Seed_Cov_A0A1->Fill(st.FitParams.Covarience.sigA1A0);
	
                _seedA0_v_seedA1->Fill(st.FitParams.A1,st.FitParams.A0);
                _seedB0_v_seedB1->Fill(st.FitParams.B1,st.FitParams.B0);
                _seedA0_v_seedB1->Fill(st.FitParams.B1,st.FitParams.A0);
	        _seedA0_v_seedB0->Fill(st.FitParams.A0,st.FitParams.B0);
	        _seedB1_v_seedA1->Fill(st.FitParams.B1,st.FitParams.A1);
	 	
                
	
		_A0MinuitFitDiff->Fill(st.MinuitFitParams.A0-st.FitEquationXYZ.Pos.X());
		_A1MinuitFitDiff->Fill(st.MinuitFitParams.A1-st.FitEquationXYZ.Dir.X());
		_B0MinuitFitDiff->Fill(st.MinuitFitParams.B0- st.FitEquationXYZ.Pos.Y());
		_B1MinuitFitDiff->Fill(st.MinuitFitParams.B1-st.FitEquationXYZ.Dir.Y());
		
 
	      _A0SeedMCDiff->Fill(st.TrueFitEquation.Pos.X()- st.FitEquationXYZ.Pos.X());
              _A1SeedMCDiff->Fill(st.TrueFitEquation.Dir.X() -st.FitEquationXYZ.Dir.X());
       	     _B0SeedMCDiff->Fill(st.TrueFitEquation.Pos.Y() - st.FitEquationXYZ.Pos.Y());
             _B1SeedMCDiff->Fill(st.TrueFitEquation.Dir.Y()  - st.FitEquationXYZ.Dir.Y());
	      
	        _A0Minuit->Fill(st.MinuitFitParams.A0);
	        _A1Minuit->Fill(st.MinuitFitParams.A1);
	        _B0Minuit->Fill(st.MinuitFitParams.B0);
	        _B1Minuit->Fill(st.MinuitFitParams.B1);

	        _seedDeltaA0_v_minuitA0->Fill(st.TrueFitEquation.Pos.X()- st.FitEquationXYZ.Pos.X(), st.MinuitFitParams.A0);
		_seedDeltaB0_v_minuitB0->Fill(st.TrueFitEquation.Pos.Y()- st.FitEquationXYZ.Pos.Y(), st.MinuitFitParams.B0);
                _seedDeltaA1_v_minuitA1->Fill(st.TrueFitEquation.Dir.X()- st.FitEquationXYZ.Dir.X(), st.MinuitFitParams.A1);
		_seedDeltaB1_v_minuitB1->Fill(st.TrueFitEquation.Dir.Y()- st.FitEquationXYZ.Dir.Y(), st.MinuitFitParams.B1);
	  
		_seedDeltaA0_v_LL->Fill(st.TrueFitEquation.Pos.X()- st.FitEquationXYZ.Pos.X(), st.DriftDiag.NLL);
		_seedDeltaA1_v_LL->Fill(st.TrueFitEquation.Dir.X()- st.FitEquationXYZ.Dir.X(), st.DriftDiag.NLL);
		_seedDeltaB0_v_LL->Fill(st.TrueFitEquation.Pos.Y()- st.FitEquationXYZ.Pos.Y(), st.DriftDiag.NLL);
		_seedDeltaB1_v_LL->Fill(st.TrueFitEquation.Dir.Y()- st.FitEquationXYZ.Dir.Y(), st.DriftDiag.NLL);

		_seedDeltaA0_v_N->Fill(st.TrueFitEquation.Pos.X()- st.FitEquationXYZ.Pos.X(), st.get_N());
		_seedDeltaA1_v_N->Fill(st.TrueFitEquation.Dir.X()- st.FitEquationXYZ.Dir.X(), st.get_N());
		_seedDeltaB0_v_N->Fill(st.TrueFitEquation.Pos.Y()- st.FitEquationXYZ.Pos.Y(), st.get_N());
		_seedDeltaB1_v_N->Fill(st.TrueFitEquation.Dir.Y()- st.FitEquationXYZ.Dir.Y(), st.get_N());

		double ptrue = sqrt(st.TrueTrueTrackDirection.Mag2());
		_seedDeltaA0_v_MOM->Fill(st.TrueFitEquation.Pos.X()- st.FitEquationXYZ.Pos.X(), ptrue);
		_seedDeltaB0_v_MOM->Fill(st.TrueFitEquation.Pos.Y()- st.FitEquationXYZ.Pos.Y(), ptrue);
		_seedDeltaA1_v_MOM->Fill(st.TrueFitEquation.Dir.X()- st.FitEquationXYZ.Dir.X(), ptrue);
		_seedDeltaB1_v_MOM->Fill(st.TrueFitEquation.Dir.Y()- st.FitEquationXYZ.Dir.Y(), ptrue);
		
		_chisq_ndf_plot_init->Fill(st.Diag.InitialChiTot);
                _chisq_ndf_plot_final->Fill(st.Diag.FinalChiTot);
                
                _chisq_ndf_plot_finalX->Fill(st.Diag.FinalChiX);
                _chisq_ndf_plot_finalY->Fill(st.Diag.FinalChiY);
                
                _chisq_ndf_plot_initX->Fill(st.Diag.InitialChiX);
                _chisq_ndf_plot_initY->Fill(st.Diag.InitialChiY);
                
                _change_chisq_ndf_plot_X->Fill(st.Diag.InitialChiX-st.Diag.FinalChiX);
                _change_chisq_ndf_plot_Y->Fill(st.Diag.InitialChiY-st.Diag.FinalChiY);
                _change_chisq_ndf_plot_Total->Fill(st.Diag.InitialChiTot-st.Diag.FinalChiTot);

		_DeltaA0_v_OUT->Fill(st.TrueFitEquation.Pos.X()- st.FitEquationXYZ.Pos.X(), st.n_outliers);
		_DeltaB0_v_OUT->Fill(st.TrueFitEquation.Pos.Y()- st.FitEquationXYZ.Pos.Y(), st.n_outliers);
		_DeltaA1_v_OUT->Fill(st.TrueFitEquation.Dir.X()- st.FitEquationXYZ.Dir.X(), st.n_outliers);
		_DeltaB1_v_OUT->Fill(st.TrueFitEquation.Dir.Y()- st.FitEquationXYZ.Dir.Y(), st.n_outliers);
		
		_NLL->Fill(st.DriftDiag.NLL);
		_LL_v_MOM->Fill(st.DriftDiag.NLL, ptrue);
		_LL_v_Nhits->Fill(st.DriftDiag.NLL, st.get_N());

                for(size_t i=0; i< st.Diag.InitErrTot.size();i++){
		    _InitErrTot->Fill(st.Diag.InitErrTot[i]); 
                    _total_residualsX_init->Fill(st.Diag.InitialResidualsX[i]);             
	            _total_pullsX_init->Fill(st.Diag.InitialResidualsX[i]/st.Diag.InitErrX[i]);
	            _InitErrX->Fill(st.Diag.InitErrX[i]);      
                    _total_residualsY_init->Fill(st.Diag.InitialResidualsY[i]);             
	            _total_pullsY_init->Fill(st.Diag.InitialResidualsY[i]/st.Diag.InitErrY[i]);
	            _InitErrY->Fill(st.Diag.InitErrY[i]); 
	        }

		for(size_t i=0; i< st.Diag.FinalErrTot.size();i++){
		    _FinalErrTot->Fill(st.Diag.FinalErrTot[i]);
                    _total_residualsX_final->Fill(st.Diag.FinalResidualsX[i]);          
	            _total_pullsX_final->Fill(st.Diag.FinalResidualsX[i]/st.Diag.FinalErrX[i]);
	            _FinalErrX->Fill(st.Diag.FinalErrX[i]);
                    _total_residualsY_final->Fill(st.Diag.FinalResidualsY[i]);            
	            _total_pullsY_final->Fill(st.Diag.FinalResidualsY[i]/st.Diag.FinalErrY[i]);
	            _FinalErrY->Fill(st.Diag.FinalErrTot[i]);
	            
	        }   
                for(size_t i=0; i<st.DriftDiag.GaussianEndDOCAs.size();i++){
		    _GaussianEndDOCAs->Fill(st.DriftDiag.GaussianEndDOCAs[i]);
		    _GaussianEndTimeResiduals->Fill(st.DriftDiag.GaussianEndTimeResiduals[i]);
		    _TrueMomNoCuts->Fill(ptrue); 
			
		}
		
		for(size_t i=0; i<st.DriftDiag.FullFitEndDOCAs.size();i++){
	            _StartDOCAs->Fill(st.DriftDiag.StartDOCAs[i]);
	            _StartTimeResiduals->Fill(st.DriftDiag.StartTimeResiduals[i]);
	            _FullFitEndDOCAs->Fill(st.DriftDiag.FullFitEndDOCAs[i]);
		    _FullFitEndTimeResiduals->Fill(st.DriftDiag.FullFitEndTimeResiduals[i]);
		    _TrueDOCAs->Fill(st.DriftDiag.TrueDOCAs[i]);
		    _TrueMomDOCACut->Fill(ptrue); 
		    _TrueTimeResiduals->Fill(st.DriftDiag.TrueTimeResiduals[i]);	      
	        }

		for(size_t i=0; i< st.DriftDiag.FullFitEndDOCAs.size();i++){
                    _total_residualsX_Minuit->Fill(st.DriftDiag.FinalResidualsX[i]);          
                    _total_residualsY_Minuit->Fill(st.DriftDiag.FinalResidualsY[i]);   
		    _minuit_pullsX_final->Fill(st.DriftDiag.FinalResidualsX[i]/st.DriftDiag.FinalErrX[i]);          
                    _minuit_pullsY_final->Fill(st.DriftDiag.FinalResidualsY[i]/st.DriftDiag.FinalErrY[i]);           
   		   _AMBIG->Fill(st.DriftDiag.RecoAmbigs[i], st.DriftDiag.TrueAmbigs[i]);
		   
	        }   
	        for(auto const& tseed : *_coscol) {   
                	CosmicTrkFitFlag const& status = tseed._status;
                	_hitsOK = status.hasAllProperties(CosmicTrkFitFlag::hitsOK);
                	_StraightTrackOK = status.hasAllProperties(CosmicTrkFitFlag::StraightTrackOK);
                	_StraightTrackConverged = status.hasAllProperties(CosmicTrkFitFlag::StraightTrackConverged);
                	_StraightTrackInit = status.hasAllProperties(CosmicTrkFitFlag::StraightTrackInit);
        	}
		for(size_t ich = 0;ich < _chcol->size(); ++ich){
                        ComboHit const& chit =(*_chcol)[ich];
			
                //-----------Fill diag details:----------//
                        _nhits = chit.nStrawHits();
                        _nsh = chit.nStrawHits(); 
                        panels.push_back(chit.strawId().panel());
		        planes.push_back(chit.strawId().plane());
			stations.push_back(chit.strawId().station());
		//-----------Hit details:---------------//
		        _hit_time = chit.time();
			_hit_drift_time = chit.driftTime();
                        _dt =  _hit_time - _cluster_time;
			}
                //----------------Get panels/planes/stations per track:------------------//
                _n_panels = std::set<float>( panels.begin(), panels.end() ).size();
		_n_planes = std::set<float>( planes.begin(), planes.end() ).size();
		_n_stations = std::set<float>( stations.begin(), stations.end() ).size();
	 
      }//end analyze
      _cosmic_analysis->Fill();
      
      cout<<"true + "<< _AMBIG->GetBinContent(2,2)/_AMBIG->Integral()<<endl;
      cout<<"true - "<< _AMBIG->GetBinContent(1,1)/_AMBIG->Integral()<<endl;	
      cout<<"false + "<< _AMBIG->GetBinContent(2,1)/_AMBIG->Integral()<<endl;
      cout<<"false - "<< _AMBIG->GetBinContent(1,2)/_AMBIG->Integral()<<endl;

     }
      bool CosmicAnalyzer::findData(const art::Event& evt){
	_chcol = 0; 
        _tccol = 0;
        _coscol = 0; 
	auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
	_chcol = chH.product();
	auto tcH = evt.getValidHandle<TimeClusterCollection>(_tctag);
	_tccol =tcH.product();
	auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_costag);
	_coscol =stH.product();
        if(_mcdiag){
           auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
           _mcdigis = mcdH.product();
           _toff.updateMap(evt);
        }
	return _chcol != 0 && _tccol!=0 && _coscol !=0 && (_mcdigis != 0 || !_mcdiag);
       }

}

using mu2e::CosmicAnalyzer;
DEFINE_ART_MODULE(CosmicAnalyzer);
    

    


