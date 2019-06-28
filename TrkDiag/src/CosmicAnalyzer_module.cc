#define _USE_MATH_DEFINES
//S. Middleton, Feb 2019
// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
//Mu2e Data Prods:
#include "DataProducts/inc/threevec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

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

      //Some Diag histograms:
      TH1F* _chisq_quant;
      TH2F* _chiY_v_true_trackposY;
      TH2F* _chiX_v_true_trackposX;
      TH2F* _chi_v_true_theta;
    
      TH2F* _A0_v_A1;
      TH2F* _B0_v_B1;
      TH2F* _A0_v_B1;
      TH2F* _B1_v_A1;
      TH2F* _A0_v_B0;
      TH1F* _chisq_ndf_plot_init;
      TH1F* _total_residualsX_init;
      TH1F* _total_pullsX_init;
      TH1F* _total_residualsY_init;
      TH1F* _total_pullsY_init;
      TH1F* _hiterrs_init;
      TH1F* _a1;
      TH1F* _b1;
      TH1F* _a0;
      TH1F* _b0;
      TH1F* _truea1;
      TH1F* _trueb1;
      TH1F* _truea0;
      TH1F* _trueb0;
      TH1F* _niters;
      TH1F* _chisq_ndf_plot_final;
      TH1F* _chisq_ndf_plot_true_fit;
      TH1F* _chisq_ndf_plot_finalX;
      TH1F* _chisq_ndf_plot_finalY;
      TH1F* _chisq_ndf_plot_initX;
      TH1F* _chisq_ndf_plot_initY;
      TH1F* _change_chisq_ndf_plot_X;
      TH1F* _change_chisq_ndf_plot_Y;
      TH1F* _change_chisq_ndf_plot_Total;
      TH1F* _total_residualsX_final;
      TH1F* _total_pullsX_final;
      TH1F* _total_residualsY_final;
      TH1F* _total_pullsY_final;
      TH1F* _hiterrs_final;
      TH1F* _ErrX;
      TH1F* _ErrY;
      TH1F* _InitErrX;
      TH1F* _InitErrY;
      TH1F* _mc_phi_angle;
      TH1F* _reco_phi_angle;
      TH1F* _mc_theta_angle;
      TH1F* _reco_theta_angle;
      TH1F* _parameter_pulls_A0;
      TH1F* _parameter_pulls_A1;
      TH1F* _parameter_pulls_B0;
      TH1F* _parameter_pulls_B1;
      TH1F* _Cov_Fit_A0;
      TH1F* _Cov_Fit_A1;
      TH1F* _Cov_Fit_B0;
      TH1F* _Cov_Fit_B1;
      TProfile* _A0pull_v_theta_true;
      TProfile* _B0pull_v_theta_true;
      TProfile* _A1pull_v_theta_true;
      TProfile* _B1pull_v_theta_true;
      TGraph* EFF;
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
      Float_t _hit_time, _hit_drift_time, _cluster_time;
	
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
       _cosmic_analysis->Branch("hitsOK",&_hitsOK,"hitsOK/B");
       _cosmic_analysis->Branch("StraightTrackInit",&_StraightTrackInit,"StraightTrackInit/B");
       _cosmic_analysis->Branch("StraightTrackOK",&_StraightTrackOK,"StraightTrackOK/B");
       _cosmic_analysis->Branch("StraightTrackConverged",&_StraightTrackConverged,"StraightTrackConverged/B");
        //Extra histograms for Fit Diags:
        _chisq_quant = tfs->make<TH1F>("PDF","PDF" ,50,0, 1);
	_chisq_quant->GetXaxis()->SetTitle("PDF");
	
	_chisq_ndf_plot_init = tfs->make<TH1F>("init chisq_ndf_plot","init chisq_ndf_plot" ,100,0, 10);
	_chisq_ndf_plot_init->GetXaxis()->SetTitle("Init. #Chi^{2}/N");
	
	_chisq_ndf_plot_final = tfs->make<TH1F>("final chisq_ndf_plot","final chisq_ndf_plot" ,100,0,10);
	_chisq_ndf_plot_final->GetXaxis()->SetTitle("Final #Chi^{2}/N");
	
	_chisq_ndf_plot_finalX = tfs->make<TH1F>("final chisq_ndf_plot X''","final chisq_ndf_plot X''" ,100,0,10);
	_chisq_ndf_plot_finalX->GetXaxis()->SetTitle("Final X '' #Chi^{2}/N");
	
        _chisq_ndf_plot_finalY = tfs->make<TH1F>("final chisq_ndf_plot Y''","final chisq_ndf_plot Y''" ,100,0,10);
	_chisq_ndf_plot_finalY->GetXaxis()->SetTitle("Final Y '' #Chi^{2}/N");
	
	_chisq_ndf_plot_true_fit = tfs->make<TH1F>("#chi^{2}/ndf StrawDigi Fit","#chi^{2}/ndf StrawDigi Fit" ,100,0,10);
	_chisq_ndf_plot_true_fit->GetXaxis()->SetTitle("#chi^{2}/ndf StrawDigi Fit");
	
	_chisq_ndf_plot_initX = tfs->make<TH1F>("initial chisq_ndf_plot X''","initial chisq_ndf_plot X''" ,100,0,10);
	_chisq_ndf_plot_initX->GetXaxis()->SetTitle("Initial X '' #Chi^{2}/N");
	
        _chisq_ndf_plot_initY = tfs->make<TH1F>("initial chisq_ndf_plot Y''","initial chisq_ndf_plot Y''" ,100,0,10);
	_chisq_ndf_plot_initY->GetXaxis()->SetTitle("Initial Y '' #Chi^{2}/N");
	
	_change_chisq_ndf_plot_X = tfs->make<TH1F>("Total Change in Chisq/ ndf X''","Total Change in Chisq/ ndf X''" ,50,-10,20);
	_change_chisq_ndf_plot_X->GetXaxis()->SetTitle("Initial X '' #Chi^{2}/N");
	
	_change_chisq_ndf_plot_Y = tfs->make<TH1F>("Total Change in Chisq/ ndf Y''","Total Change in Chisq/ ndf Y''" ,50,-10,20);
	_change_chisq_ndf_plot_Y->GetXaxis()->SetTitle("Initial Y '' #Chi^{2}/N");
	
	_change_chisq_ndf_plot_Total = tfs->make<TH1F>("Total Change in #Chi^{2}/ ndf","Total Change in #Chi^{2}/ ndf" ,50,-5,20);
	_change_chisq_ndf_plot_Total->GetXaxis()->SetTitle("Total Change in #Chi^{2}/ndf");
	
        _total_residualsX_init = tfs->make<TH1F>("Initial Residuals X'' ","Initial Residuals X'' " ,200,-500,500);
	_total_residualsX_init->GetXaxis()->SetTitle("Initial Residual X'' [mm]");
	_total_residualsX_init->SetStats();
	
	_total_residualsY_init = tfs->make<TH1F>("Initial Residuals Y''","Initial Residuals Y'' " ,200,-500,500);
	_total_residualsY_init->GetXaxis()->SetTitle("Initial Residual Y'' [mm]");
	_total_residualsY_init->SetStats();

	_total_pullsX_init = tfs->make<TH1F>("Initial Pull X","Initial Pull X''" ,200,-50, 50);
	_total_pullsX_init->GetXaxis()->SetTitle("Pull X");
	_total_pullsX_init->SetStats();
	
	_total_pullsY_init = tfs->make<TH1F>("Initial Pull Y''","Initial Pull Y''" ,200,-50, 50);
	_total_pullsY_init->GetXaxis()->SetTitle("Initial Pull Y");
	_total_pullsY_init->SetStats();
       
        _hiterrs_init = tfs->make<TH1F>("Initial Total Hit Error","Initial Total Hit Error" ,50,0, 100);
	_hiterrs_init->GetXaxis()->SetTitle("Initial Hit Error in Track Frame [mm]");
	_hiterrs_init->SetStats();
	
	 _total_residualsX_final = tfs->make<TH1F>("Final Residuals X''  ","Final Residuals X''#chi^{2}/dof all  " ,500,-100,100);
	_total_residualsX_final->GetXaxis()->SetTitle("Residual X'' [mm]");
	_total_residualsX_final->SetStats();
	
	_total_residualsY_final = tfs->make<TH1F>("Final Residuals Y'' ","FinalResiduals Y'' #chi^{2}/dof all " ,500,-100,100);
	_total_residualsY_final->GetXaxis()->SetTitle("Final Residual Y'' [mm]");
	_total_residualsY_final->SetStats();
	
	_total_pullsX_final = tfs->make<TH1F>("Final Pull X'' #chi^{2}/dof all ","Final Pull X''#chi^{2}/dof all " ,50,-2, 2);
	_total_pullsX_final->GetXaxis()->SetTitle("Final Pull X''");
	//_total_pullsX_final->SetStats(0);
	
	_total_pullsY_final = tfs->make<TH1F>("Final Pull Y''#chi^{2}/dof all ","Final Pull Y''#chi^{2}/dof all " ,50,-2, 2);
	_total_pullsY_final->GetXaxis()->SetTitle("Final Pull Y''");
	//_total_pullsY_final->SetStats(0);
       
        _hiterrs_final = tfs->make<TH1F>("Final Total Hit Error #chi^{2}/dof all ","Final Total Hit Error #chi^{2}/dof all " ,50,0, 100);
	_hiterrs_final->GetXaxis()->SetTitle("Hit Error in Track Frame [mm]");
	_hiterrs_final->SetStats();
	
	_a0 = tfs->make<TH1F>("Track Parameter A0 ","Track Parameter A0 " ,50,-2000, 2000);
	_a0->GetXaxis()->SetTitle("Track Parameter A0");
	_a0->SetStats();
	
	_b0 = tfs->make<TH1F>("Track Parameter B0 ","Track Parameter B0  " ,50,-2000, 2000);
	_b0->GetXaxis()->SetTitle("Track Parameter B0");
	_b0->SetStats();
	
	_a1 = tfs->make<TH1F>("Track Parameter A1  ","Track Parameter A1 " ,50,-1,1);
	_a1->GetXaxis()->SetTitle("Track Parameter A1");
	_a1->SetStats();
	
	_b1 = tfs->make<TH1F>("Track Parameter B1  ","Track Parameter B1  " ,50,-1,1);
	_b1->GetXaxis()->SetTitle("Track Parameter B1");
	_b1->SetStats();
	
	_truea0 = tfs->make<TH1F>("True Track Parameter A0 ","True Track Parameter A0 " ,50,-2000, 2000);
	_truea0->GetXaxis()->SetTitle("Track Parameter A0");
	_truea0->SetStats();
	
	_trueb0 = tfs->make<TH1F>("True Track Parameter B0 "," True Track Parameter B0  " ,50,-2000, 2000);
	_trueb0->GetXaxis()->SetTitle("Track Parameter B0");
	_trueb0->SetStats();
	
	_truea1 = tfs->make<TH1F>("True Track Parameter A1  ","True Track Parameter A1 " ,50,-1,1);
	_truea1->GetXaxis()->SetTitle("Track Parameter A1");
	_truea1->SetStats();
	
	_trueb1 = tfs->make<TH1F>("True Track Parameter B1  ","True Track Parameter B1  " ,50,-1,1);
	_trueb1->GetXaxis()->SetTitle("Track Parameter B1");
	_trueb1->SetStats();
	
	_niters  = tfs->make<TH1F>("Number of Iterations Unitl Converged  ","Number of Iterations Unitl Converged " ,100,0, 100);
	_niters->GetXaxis()->SetTitle("Number of Iterations Unitl Converged");
	_niters->SetStats();
	
	_ErrX = tfs->make<TH1F>("Errors in X''","Errors in X'' " ,50,0, 100);
	_ErrX->GetXaxis()->SetTitle("#sigma_{X''} [mm]");
	_ErrX->SetStats();
	
	_ErrY = tfs->make<TH1F>("Errors in Y''","Errors in Y'' " ,50,0, 100);
	_ErrY->GetXaxis()->SetTitle("#sigma_{Y''} [mm]");
	_ErrY->SetStats();
	
	_InitErrX = tfs->make<TH1F>("Initial Errors in X''","Initial Errors in X'' " ,50,0, 100);
	_InitErrX->GetXaxis()->SetTitle("#sigma_{X''} [mm]");
	_InitErrX->SetStats();
	
	_InitErrY = tfs->make<TH1F>("Initial Errors in Y''","Initial Errors in Y'' " ,50,0, 100);
	_InitErrY->GetXaxis()->SetTitle("#sigma_{Y''} [mm]");
	_InitErrY->SetStats();
	
	_mc_phi_angle = tfs->make<TH1F>("#phi_{true, fit}","#phi_{true, fit}" ,100,-3.141529,3.141529);
	_mc_phi_angle->GetXaxis()->SetTitle("#phi_{true, fit} [rad]");
	
	
	_reco_phi_angle = tfs->make<TH1F>("#phi Angle of Reconstructred Track","#phi_{reco} Angle of Tracl " ,20,0,3.141529);
	_reco_phi_angle->GetXaxis()->SetTitle("#phi_{reco} [rad]");
	_reco_phi_angle->SetStats();
	
	_mc_theta_angle = tfs->make<TH1F>("#theta_{true, fit}","#theta_{true, fit}" ,20,0,3.141529);
	_mc_theta_angle->GetXaxis()->SetTitle("#theta_{true, fit} [rad]");
	
	EFF = tfs->make<TGraph>(100);//"Reco Eff","Reco EFF " ,100,0,3.141529, 0, 10);
	
	_A0_v_A1 = tfs->make<TH2F>("A_{0} v A_{1}", "A_{0} v A_{1}", 50, -0.2, 0.2, 50, -1500, 1500);
	_A0_v_A1->GetXaxis()->SetTitle("A_{1}");
	_A0_v_A1->GetYaxis()->SetTitle("A_{0}");
	
	
	_B0_v_B1 = tfs->make<TH2F>("B_{0} v B_{1}", "B_{0} v B_{1}", 50, -0.2, 0.2, 50, -1500, 1500);
	_B0_v_B1->GetXaxis()->SetTitle("B_{1}");
	_B0_v_B1->GetYaxis()->SetTitle("B_{0}");
	
	_A0_v_B1 = tfs->make<TH2F>("A_{0} v B_{1}", "A_{0} v B_{1}", 50, -0.2, 0.2, 50, -1500, 1500);
	_A0_v_B1->GetXaxis()->SetTitle("B_{1}");
	_A0_v_B1->GetYaxis()->SetTitle("A_{0");
	
	_A0_v_B0 = tfs->make<TH2F>("B_{0} v A_{0}", "B_{0} v A_{0}", 50, -1500,1500, 50, -2000, 2000);
	_A0_v_B0->GetXaxis()->SetTitle("B_{0}");
	_A0_v_B0->GetYaxis()->SetTitle("A_{0} ");
	
	_B1_v_A1 = tfs->make<TH2F>("B_{1} v A_{1}", "B_{1} v A_{1}", 100, -0.2,0.2, 100, -0.2,0.2);
	_B1_v_A1->GetXaxis()->SetTitle("B_{1}");
	_B1_v_A1->GetYaxis()->SetTitle("A_{1} ");
	
	_parameter_pulls_A0 = tfs->make<TH1F>("Parameter Pull A_{0} ","Parameter Pull A_{0}" ,100,-200, 200);
	_parameter_pulls_A0->GetXaxis()->SetTitle("Parameter Pull #A_{0}");
	_parameter_pulls_A0->SetStats();
	
	_parameter_pulls_A1 = tfs->make<TH1F>("Parameter Pull A_{1}","Parameter Pull A_{1}" ,50,-10, 10);
	_parameter_pulls_A1->GetXaxis()->SetTitle("Parameter Pull A_{1}");
	_parameter_pulls_A1->SetStats();
	
	_parameter_pulls_B0 = tfs->make<TH1F>("Parameter Pull B_{0}" ,"Parameter Pull B_{0}",100,-200, 200);
	_parameter_pulls_A0->GetXaxis()->SetTitle("Parameter Pull B_{0}");
	_parameter_pulls_A0->SetStats();
	
	_parameter_pulls_B1 = tfs->make<TH1F>("Parameter Pull B_{1}", "Parameter Pull B_{1}" ,50,-10, 10);
	_parameter_pulls_B1->GetXaxis()->SetTitle("Parameter Pull B_{1}");
	_parameter_pulls_B1->SetStats();
	
	_chiX_v_true_trackposX = tfs->make<TH2F>("#chi^{2}/ndf v A_{0,true}", "#chi^{2}/ndf v A_{0,true}", 100, -2000,2000, 100, 0,5);
	_chiX_v_true_trackposX->GetXaxis()->SetTitle("A_{0}");
	_chiX_v_true_trackposX->GetYaxis()->SetTitle("#chi{2}/ndf");
	
	_chiY_v_true_trackposY = tfs->make<TH2F>("#chi^{2}/ndf v B_{0,true}", "#chi^{2}/ndf v B_{0,true}", 100, -2000,2000, 100, 0,5);
	_chiY_v_true_trackposY->GetXaxis()->SetTitle("B_{0}");
	_chiY_v_true_trackposY->GetYaxis()->SetTitle("#chi{2}/ndf");
	
	_chi_v_true_theta = tfs->make<TH2F>("#chi^{2}/ndf v #theta_{true}", "#chi^{2}/ndf v #theta_{true}", 100, 0, 3.141529, 100, 0,2);
	_chi_v_true_theta->GetXaxis()->SetTitle("B_{0}");
	_chi_v_true_theta->GetYaxis()->SetTitle("#chi{2}/ndf");
	
	_Cov_Fit_A0 = tfs->make<TH1F>("#sigma_{A_{0}}", "#sigma_{A_{0}}" ,100,0, 100);
	_Cov_Fit_A0->GetXaxis()->SetTitle("#sigma_{A_{0}}");
	
	
	_Cov_Fit_A1 = tfs->make<TH1F>("#sigma_{A_{1}}", "#sigma_{A_{1}}" ,100,0, 0.5);
	_Cov_Fit_A1->GetXaxis()->SetTitle("#sigma_{A_{1}}");
	
	
	_Cov_Fit_B0 = tfs->make<TH1F>("#sigma_{B_{0}}", "#sigma_{B_{0}}" ,100,0, 100);
	_Cov_Fit_B0->GetXaxis()->SetTitle("#sigma_{B_{0}}");
	
	
	_Cov_Fit_B1 = tfs->make<TH1F>("#sigma_{B_{1}}", "#sigma_{B_{1}}" ,100,0, 0.5);
	_Cov_Fit_B1->GetXaxis()->SetTitle("#sigma_{B_{1}}");
	
	_A0pull_v_theta_true= tfs->make<TProfile>("Parameter Pull A_{0} v #theta_{true,fit}  ","Parameter Pull A_{0} v #theta_{true,fit} ",100,0,3.141529,-200,200);
	_A0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_A0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{0}");
	
	_B0pull_v_theta_true= tfs->make<TProfile>("Parameter Pull B_{0} v #theta_{true,fit}  ","Parameter Pull B_{0} v #theta_{true,fit} ",100,0,3.141529,-200,200);
	_B0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_B0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{0}");
	
	_A1pull_v_theta_true= tfs->make<TProfile>("Parameter Pull A_{1} v #theta_{true,fit}  ","Parameter Pull A_{1} v #theta_{true,fit} ",100,0,3.141529,-10,10);
	_A1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_A1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{1}");
	
	_B1pull_v_theta_true= tfs->make<TProfile>("Parameter Pull B_{1} v #theta_{true,fit}  ","Parameter Pull B_{1} v #theta_{true,fit} ",100,0,3.141529,-10,10);
	_B1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
	_B1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{1}");
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
		
		//TrkFitFlag const& status = sts._status;
		std::vector<int> panels, planes, stations;
                _chisq_ndf_plot_init->Fill(st.get_initchisq_dof());
                _chisq_ndf_plot_final->Fill(st.get_finalchisq_dof());
                
                _chisq_ndf_plot_finalX->Fill(st.get_finalchisq_dofX());
                _chisq_ndf_plot_finalY->Fill(st.get_finalchisq_dofY());
              
                
                _chisq_ndf_plot_initX->Fill(st.get_initchisq_dofX());
                _chisq_ndf_plot_initY->Fill(st.get_initchisq_dofY());
                
                _change_chisq_ndf_plot_X->Fill(st.get_initchisq_dofX()-st.get_finalchisq_dofX());
                _change_chisq_ndf_plot_Y->Fill(st.get_initchisq_dofY()-st.get_finalchisq_dofY());
                _change_chisq_ndf_plot_Total->Fill(st.get_initchisq_dof()-st.get_finalchisq_dof());
                //_chisq_ndf_plot_true_fit->Fill(st.get_true_finalchisq_dof());
                
                _a1->Fill(st.get_track_parameters()[1]);
                _b1->Fill(st.get_track_parameters()[3]);
                _a0->Fill(st.get_track_parameters()[0]);
                _b0->Fill(st.get_track_parameters()[2]);
                
                _chiX_v_true_trackposX->Fill(st.get_track_parameters()[0], st.get_initchisq_dofX());
                _chiY_v_true_trackposY->Fill(st.get_track_parameters()[2], st.get_initchisq_dofY());
               
                if(_mcdiag ){
                	_mc_phi_angle->Fill(st.get_true_phi());//atan2(st.get_true_track_direction().y(),st.get_true_track_direction().x()));//st.get_true_phi());  
	                //float val = sqrt(st.get_true_track_direction().Perp2()) ;
	                _mc_theta_angle->Fill(st.get_true_theta());//atan2(val,st.get_true_track_direction().z()));                   
	                _chi_v_true_theta->Fill(st.get_true_theta(),st.get_true_finalchisq_dof());
	                _mc_theta_angle->SetStats(0);
			_mc_phi_angle->SetStats(0);
		        _truea1->Fill(st.get_true_track_parameters()[1]);
		        _trueb1->Fill(st.get_true_track_parameters()[3]);
		        _truea0->Fill(st.get_true_track_parameters()[0]);
		        _trueb0->Fill(st.get_true_track_parameters()[2]);
                        _parameter_pulls_A0->Fill((st.get_true_track_parameters()[0] - st.get_track_parameters()[0] )/(st.get_cov()[0]));
		        _parameter_pulls_A1->Fill((st.get_true_track_parameters()[1] - st.get_track_parameters()[1] )/(st.get_cov()[1]));
		        _parameter_pulls_B0->Fill(st.get_true_theta(), (st.get_true_track_parameters()[2] - st.get_track_parameters()[2] )/(st.get_cov()[2]));
		        _parameter_pulls_B1->Fill(st.get_true_theta(), (st.get_true_track_parameters()[3] - st.get_track_parameters()[3] )/(st.get_cov()[3]));
		        _A0pull_v_theta_true->Fill(st.get_true_theta(), (st.get_true_track_parameters()[0] - st.get_track_parameters()[0] )/(st.get_cov()[0]));
		        _A1pull_v_theta_true->Fill(st.get_true_theta(), (st.get_true_track_parameters()[1] - st.get_track_parameters()[1] )/(st.get_cov()[1]));
		        _B0pull_v_theta_true->Fill(st.get_true_theta(), (st.get_true_track_parameters()[2] - st.get_track_parameters()[2] )/(st.get_cov()[2]));
		        _B1pull_v_theta_true->Fill(st.get_true_theta(), (st.get_true_track_parameters()[3] - st.get_track_parameters()[3] )/(st.get_cov()[3]));
		       
                }
                _Cov_Fit_A0->Fill(st.get_cov()[0]);
                _Cov_Fit_A1->Fill(st.get_cov()[1]);
                _Cov_Fit_B0->Fill(st.get_cov()[2]);
                _Cov_Fit_B1->Fill(st.get_cov()[3]);
               
                _A0_v_A1->Fill(st.get_track_parameters()[1],st.get_track_parameters()[0]);
                _B0_v_B1->Fill(st.get_track_parameters()[3],st.get_track_parameters()[2]);
                _A0_v_B1->Fill(st.get_track_parameters()[3],st.get_track_parameters()[0]);
	        _A0_v_B0->Fill(st.get_track_parameters()[0],st.get_track_parameters()[2]);
	        _B1_v_A1->Fill(st.get_track_parameters()[3],st.get_track_parameters()[1]);
	 	if(st.get_chi2_quant()>0){//TODO: Fix some weird root issue here.....
                	_chisq_quant->Fill(st.get_chi2_quant());
                	std::cout<<"pdf "<<st.get_chi2_quant()<<std::endl;
                }
                
                _reco_phi_angle->Fill(st.get_fit_phi());
                
		_niters->Fill(st.get_iter()); 
                //-----------Fill Hist Details:----------//
		for(size_t i=0; i< st.get_init_hit_errorsTotal().size();i++){
		    _hiterrs_init->Fill(st.get_init_hit_errorsTotal()[i]); 
                    _total_residualsX_init->Fill(st.get_init_fit_residualsX()[i]);             
	            _total_pullsX_init->Fill(st.get_init_fit_pullsX()[i]);
	            _InitErrX->Fill(st.get_init_fit_residual_errorsX()[i]);      
                    _total_residualsY_init->Fill(st.get_init_fit_residualsY()[i]);             
	            _total_pullsY_init->Fill(st.get_init_fit_pullsY()[i]);
	            _InitErrY->Fill(st.get_init_fit_residual_errorsY()[i]); 
	        }
                //-----------Fill Hist Details:----------//
		for(size_t i=0; i< st.get_final_hit_errorsTotal().size();i++){
		    _hiterrs_final->Fill(st.get_final_hit_errorsTotal()[i]);
                    _total_residualsX_final->Fill(st.get_final_fit_residualsX()[i]);          
	            _total_pullsX_final->Fill( st.get_final_fit_pullsX()[i]);
	            _ErrX->Fill(st.get_final_fit_residual_errorsX()[i]);
                    _total_residualsY_final->Fill(st.get_final_fit_residualsY()[i]);            
	            _total_pullsY_final->Fill(st.get_final_fit_pullsY()[i]);
	            _ErrY->Fill(st.get_final_fit_residual_errorsY()[i]);
	            
	        }   
	        for(auto const& tseed : *_coscol) {   
                	TrkFitFlag const& status = tseed._status;
                	_hitsOK = status.hasAllProperties(TrkFitFlag::hitsOK);
                	_StraightTrackOK = status.hasAllProperties(TrkFitFlag::StraightTrackOK);
                	_StraightTrackConverged = status.hasAllProperties(TrkFitFlag::StraightTrackConverged);
                	_StraightTrackInit = status.hasAllProperties(TrkFitFlag::StraightTrackInit);
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
			}
                //----------------Get panels/planes/stations per track:------------------//
                _n_panels = std::set<float>( panels.begin(), panels.end() ).size();
		_n_planes = std::set<float>( planes.begin(), planes.end() ).size();
		_n_stations = std::set<float>( stations.begin(), stations.end() ).size();
	 
		//}//endS ST
	        
      }//end analyze
      
      //_mc_phi_angle->SaveAs("phiMCfit.root");
      
      _cosmic_analysis->Fill();
      /*
      if(_mcdiag > 0){
      
      for ( auto const& digimc : *_mcdigis ){
        
          //const CLHEP::Hep3Vector p = digimc.stepPointMC(mu2e::StrawEnd::cal)->momentum();
          art::Ptr<mu2e::SimParticle>& sim = digimc.stepPointMC(mu2e::StrawEnd::cal)->simParticle();
          double mag = sqrt((sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x()));
      	  const double theta_start = acos(sim.startMomentum().vect().z()/mag);
          const double phi_start = atan(sim.startMomentum().vect().y()/sim.startMomentum().vect().x());
          std::cout<<theta_start<<" "<<phi_start<<" "<<mag<<std::endl;
         }
      }*/
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
	   // update time offsets
           _toff.updateMap(evt);
        }
	return _chcol != 0 && _tccol!=0 && _coscol !=0 && (_mcdigis != 0 || !_mcdiag);
       }

}//End Namespace Mu2e

using mu2e::CosmicAnalyzer;
DEFINE_ART_MODULE(CosmicAnalyzer);
    

    


