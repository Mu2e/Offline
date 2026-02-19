//Author: S Middleton
//Date: April 2019
//Purpose: (DEPRECATED - Almost!) Not very nice analyzer for Cosmics. Would be better to use the CosmicTrackDetails. This will be removed soon but might be useful for comparison with old code.

#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "Offline/CosmicReco/inc/CosmicTrackFit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CosmicReco/inc/ComboHitInfoMC.hh"
//Mu2e Data Prods:
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/GenVector.hh"

//Utilities
#include "Offline/CosmicReco/inc/DriftFitUtils.hh"
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"

// Mu2e diagnostics
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

//ROOT
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
        struct Config{
              using Name=fhicl::Name;
              using Comment=fhicl::Comment;
              fhicl::Atom<int> diag{Name("diagLevel"), Comment("set to 1 for info"),1};
              fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),true};
              fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
              fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
              fhicl::Atom<art::InputTag> costag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmci track seed collection")};
              fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag"),"makeSD"};
            };
        typedef art::EDAnalyzer::Table<Config> Parameters;

        explicit CosmicAnalyzer(const Parameters& conf);
        virtual ~CosmicAnalyzer();
        virtual void beginJob() override;
        virtual void beginRun(const art::Run& r) override;
        virtual void analyze(const art::Event& e) override;
        virtual void endJob() override;
    private:

        Config _conf;

        int  _diag;
        bool _mcdiag;
        std::ofstream outputfile;
        art::InputTag   _chtag;//combo
        art::InputTag   _tctag;//timeclusters
        art::InputTag   _costag;//Striaght tracks
        art::InputTag   _mcdigistag; //MC digis
        const ComboHitCollection* _chcol;
        const CosmicTrackSeedCollection* _coscol;
        const TimeClusterCollection* _tccol;
        const StrawDigiMCCollection* _mcdigis;
        CosmicTrackMCInfo trueinfo;


        //TTree Info:
        TTree* _cosmic_analysis;

        //Seed Info:
        TH1F* _seed_a1;
        TH1F* _seed_b1;
        TH1F* _seed_a0;
        TH1F* _seed_b0;

        //Seed Info in Det Frame:
        TH1F* _seed_a1XYZ;
        TH1F* _seed_b1XYZ;
        TH1F* _seed_a0XYZ;
        TH1F* _seed_b0XYZ;
        TH1F* _niters;

        //Looking for geo bias:
        TH2F* _chi_v_true_theta;

        //Looking for hidden correlations between seed paramters
        TH2F* _seedA0_v_seedA1;
        TH2F* _seedB0_v_seedB1;
        TH2F* _seedA0_v_seedB1;
        TH2F* _seedB1_v_seedA1;
        TH2F* _seedA0_v_seedB0;

        //True MC paramets in global coords:
        TH1F* _truea1XYZ;
        TH1F* _trueb1XYZ;
        TH1F* _truea0XYZ;
        TH1F* _trueb0XYZ;

        //Angles Info:
        TH1F* _mc_phi_angle;
        TH1F* _reco_phi_angle;
        TH1F* _mc_theta_angle;
        TH1F* _reco_theta_angle;

        //Looking for difference between seed parameter and true parameters:
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

        //look for correlataions between true angle and parameters pulls
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

        //difference between seed and minuit fits
        TH1F* _A0MinuitFitDiff;
        TH1F* _A1MinuitFitDiff;
        TH1F* _B0MinuitFitDiff;
        TH1F* _B1MinuitFitDiff;

        //Track Parameters from end of minuit minimzation rotuine:
        TH1F* _A0Minuit;
        TH1F* _A1Minuit;
        TH1F* _B0Minuit;
        TH1F* _B1Minuit;
        //Difference in minuit derived parameters and MC true
        TH1F* _A0MinuitMCDiff;
        TH1F* _A1MinuitMCDiff;
        TH1F* _B0MinuitMCDiff;
        TH1F* _B1MinuitMCDiff;


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
        TH2F* _AMBIG;

         //Compare residual with end result :
        TH2F* _seedDeltaA0_v_minuitA0;
        TH2F* _seedDeltaA1_v_minuitA1;
        TH2F* _seedDeltaB0_v_minuitB0;
        TH2F* _seedDeltaB1_v_minuitB1;

        // add event id
        Int_t _evt;

        //Numbers:
        Int_t _nsh, _nch; // # associated straw hits / event
        Int_t _ntc; // # clusters/event
        Int_t _nhits, _nused; // # hits used
        Int_t _n_panels; // # panels
        Int_t _n_stations; // # stations
        Int_t _n_planes; // # stations

        Float_t _hit_time, _hit_drift_time, _cluster_time, _dt;

        ProditionsHandle<Tracker> _alignedTracker_h;

        Bool_t _StraightTrackInit, _StraightTrackConverged, _StraightTrackOK, _hitsOK;
        Int_t _strawid;
        vector<ComboHitInfoMC> _chinfomc;

        CosmicTrackMCInfo FitMC(const StrawDigiMCCollection*& _mcdigis);
        CosmicTrackMCInfo FillDriftMC(ComboHit const& chit, double reco_ambig, CosmicTrackMCInfo info, const Tracker* tracker);
        bool findData(const art::Event& evt);
    };

    CosmicAnalyzer::CosmicAnalyzer(const Parameters& conf) :
        art::EDAnalyzer(conf),
        _diag (conf().diag()),
        _mcdiag (conf().mcdiag()),
        _chtag (conf().chtag()),
        _tctag (conf().tctag()),
        _costag (conf().costag()),
        _mcdigistag (conf().mcdigistag())
    {
    }

    CosmicAnalyzer::~CosmicAnalyzer(){}


    void CosmicAnalyzer::beginJob() {

      if(_diag > 0){
        if(_mcdiag) {
                outputfile.open("CosmicAnalysis.csv");
        }
        art::ServiceHandle<art::TFileService> tfs;

        _cosmic_analysis=tfs->make<TTree>("cosmic_analysis"," Diagnostics for Cosmic Track Fitting");

        _cosmic_analysis->Branch("evt",&_evt,"evt/I");
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

        if(_mcdiag ){

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


                _mc_phi_angle = tfs->make<TH1F>("#phi_{true, fit}","#phi_{true, fit}" ,100,-M_PI,M_PI);
                _mc_phi_angle->GetXaxis()->SetTitle("#phi_{true, fit} [rad]");

                _mc_theta_angle = tfs->make<TH1F>("#theta_{true, fit}","#theta_{true, fit}" ,20,0,M_PI);
                _mc_theta_angle->GetXaxis()->SetTitle("#theta_{true, fit} [rad]");

                _AMBIG = tfs->make<TH2F>("True v Reco", "True v Reco", 2,-2,2,2,-2,2);
                _AMBIG->GetXaxis()->SetTitle("True");
                _AMBIG->GetYaxis()->SetTitle("Reco");

                _TrueDOCAs = tfs->make<TH1F>("True DOCAs ","TrueDOCAs" ,50, 0,10);
                _TrueDOCAs->GetXaxis()->SetTitle("True DOCAs [mm]");
                _TrueDOCAs->SetStats();

                _TrueTimeResiduals = tfs->make<TH1F>("True Time Res ","True Time Resid" ,50,0,500);
                _TrueTimeResiduals->GetXaxis()->SetTitle("True Time Res [ns]");
                _TrueTimeResiduals->SetStats();

                _A0pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull A_{0} v #theta_{true,fit}  ","Seed Parameter Pull A_{0} v #theta_{true,fit} ",100,0,M_PI,-200,200);
                _A0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
                _A0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{0}");

                _B0pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull B_{0} v #theta_{true,fit}  ","Seed Parameter Pull B_{0} v #theta_{true,fit} ",100,0,M_PI,-200,200);
                _B0pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
                _B0pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{0}");

                _A1pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull A_{1} v #theta_{true,fit}  ","Seed Parameter Pull A_{1} v #theta_{true,fit} ",100,0,M_PI,-10,10);
                _A1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
                _A1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull A_{1}");

                _B1pull_v_theta_true= tfs->make<TProfile>("Seed Parameter Pull B_{1} v #theta_{true,fit}  ","Seed Parameter Pull B_{1} v #theta_{true,fit} ",100,0,M_PI,-10,10);
                _B1pull_v_theta_true->GetXaxis()->SetTitle("#theta_{true}");
                _B1pull_v_theta_true->GetYaxis()->SetTitle("Parameter Pull B_{1}");

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
        }
        _reco_phi_angle = tfs->make<TH1F>("#phi Angle of Reconstructred Track","#phi_{reco} Angle of Tracl " ,20,0,M_PI);
        _reco_phi_angle->GetXaxis()->SetTitle("#phi_{reco} [rad]");
        _reco_phi_angle->SetStats();

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

        _FullFitEndTimeResiduals = tfs->make<TH1F>("Final Time Res ","Final Time Resid" ,50,0,500);
        _FullFitEndTimeResiduals->GetXaxis()->SetTitle("Final Time Res [ns]");
        _FullFitEndTimeResiduals->SetStats();

        _GaussianEndTimeResiduals = tfs->make<TH1F>("Gaussian Time Res ","Gaussian Time Resid" ,50,0,500);
        _GaussianEndTimeResiduals->GetXaxis()->SetTitle("Gaussian Time Res [ns]");
        _GaussianEndTimeResiduals->SetStats();

        _StartTimeResiduals = tfs->make<TH1F>("Start Time Res ","Start Time Resid" ,50,0,500);
        _StartTimeResiduals->GetXaxis()->SetTitle("Start Time Res [ns]");
        _StartTimeResiduals->SetStats();

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


        }
      }

    void CosmicAnalyzer::beginRun(const art::Run& run){}

    void CosmicAnalyzer::analyze(const art::Event& event) {
        Tracker const* tracker = _alignedTracker_h.getPtr(event.id()).get();
        _evt = event.id().event();
        if(!findData(event))
                throw cet::exception("RECO")<<"No Time Clusters in event"<< endl;

        _ntc = _tccol->size();
        _nch = _chcol->size();

        for(size_t itc=0; itc<_tccol->size();++itc){
                TimeCluster tc = (*_tccol)[itc];
                _cluster_time =  tc._t0._t0;

        }

        for(size_t ist = 0;ist < _coscol->size(); ++ist){

                CosmicTrackSeed sts =(*_coscol)[ist];
                CosmicTrack st = sts._track;
                TrkFitFlag const& status = sts._status;
                if (!status.hasAllProperties(TrkFitFlag::helixOK) ){ continue;}
                if(st.converged == false or st.minuit_converged  == false) { continue;}
                std::vector<int> panels, planes, stations;

                _reco_phi_angle->Fill(acos(st.FitEquation.Dir.x()/st.FitEquation.Dir.Mag2()));
                _niters->Fill(st.get_iter());

                _seed_a1->Fill(st.FitParams.A1);
                _seed_b1->Fill(st.FitParams.B1);
                _seed_a0->Fill(st.FitParams.A0);
                _seed_b0->Fill(st.FitParams.B0);

                _seed_a1XYZ->Fill(st.FitEquation.Dir.X());
                _seed_b1XYZ->Fill(st.FitEquation.Dir.Y());
                _seed_a0XYZ->Fill(st.FitEquation.Pos.X());
                _seed_b0XYZ->Fill(st.FitEquation.Pos.Y());

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

                _A0MinuitFitDiff->Fill(st.MinuitParams.A0-st.FitEquation.Pos.X());
                _A1MinuitFitDiff->Fill(st.MinuitParams.A1-st.FitEquation.Dir.X());
                _B0MinuitFitDiff->Fill(st.MinuitParams.B0- st.FitEquation.Pos.Y());
                _B1MinuitFitDiff->Fill(st.MinuitParams.B1-st.FitEquation.Dir.Y());

                _A0Minuit->Fill(st.MinuitParams.A0);
                _A1Minuit->Fill(st.MinuitParams.A1);
                _B0Minuit->Fill(st.MinuitParams.B0);
                _B1Minuit->Fill(st.MinuitParams.B1);

               if(_mcdiag){

                        trueinfo = FitMC(_mcdigis);

                        _mc_phi_angle->Fill(trueinfo.TruePhi);
                        _mc_theta_angle->Fill(trueinfo.TrueTheta);

                        _truea1XYZ->Fill(trueinfo.TrueFitEquation.Dir.X());
                        _trueb1XYZ->Fill(trueinfo.TrueFitEquation.Dir.Y());
                        _truea0XYZ->Fill(trueinfo.TrueFitEquation.Pos.X());
                        _trueb0XYZ->Fill(trueinfo.TrueFitEquation.Pos.Y());

                        if(st.MinuitParams.Covarience.sigA0 !=0){
                                _A0MinuitMCDiff->Fill((trueinfo.TrueFitEquation.Pos.X()- st.MinuitParams.A0));
                        }if(st.MinuitParams.Covarience.sigA1 !=0){
                                _A1MinuitMCDiff->Fill((trueinfo.TrueFitEquation.Dir.X() - st.MinuitParams.A1));
                                _B0MinuitMCDiff->Fill((trueinfo.TrueFitEquation.Pos.Y()- st.MinuitParams.B0));

                        }if(st.MinuitParams.Covarience.sigB1 !=0){
                                _B1MinuitMCDiff->Fill((trueinfo.TrueFitEquation.Dir.Y() - st.MinuitParams.B1));
                       }

                              _A0SeedMCDiff->Fill(trueinfo.TrueFitEquation.Pos.X()- st.FitEquation.Pos.X());
                              _A1SeedMCDiff->Fill(trueinfo.TrueFitEquation.Dir.X() -st.FitEquation.Dir.X());
                                    _B0SeedMCDiff->Fill(trueinfo.TrueFitEquation.Pos.Y() - st.FitEquation.Pos.Y());
                             _B1SeedMCDiff->Fill(trueinfo.TrueFitEquation.Dir.Y()  - st.FitEquation.Dir.Y());

                        _seedDeltaA0_v_minuitA0->Fill(trueinfo.TrueFitEquation.Pos.X()- st.FitEquation.Pos.X(), st.MinuitParams.A0);
                        _seedDeltaB0_v_minuitB0->Fill(trueinfo.TrueFitEquation.Pos.Y()- st.FitEquation.Pos.Y(), st.MinuitParams.B0);
                        _seedDeltaA1_v_minuitA1->Fill(trueinfo.TrueFitEquation.Dir.X()- st.FitEquation.Dir.X(), st.MinuitParams.A1);
                        _seedDeltaB1_v_minuitB1->Fill(trueinfo.TrueFitEquation.Dir.Y()- st.FitEquation.Dir.Y(), st.MinuitParams.B1);


            }

                _chisq_ndf_plot_init->Fill(st.Diag.InitialChiTot);
                _chisq_ndf_plot_final->Fill(st.Diag.FinalChiTot);

                _chisq_ndf_plot_finalX->Fill(st.Diag.FinalChiX);
                _chisq_ndf_plot_finalY->Fill(st.Diag.FinalChiY);

                _chisq_ndf_plot_initX->Fill(st.Diag.InitialChiX);
                _chisq_ndf_plot_initY->Fill(st.Diag.InitialChiY);

                _change_chisq_ndf_plot_X->Fill(st.Diag.InitialChiX-st.Diag.FinalChiX);
                _change_chisq_ndf_plot_Y->Fill(st.Diag.InitialChiY-st.Diag.FinalChiY);
                _change_chisq_ndf_plot_Total->Fill(st.Diag.InitialChiTot-st.Diag.FinalChiTot);

                for(size_t i=0; i<sts._straw_chits.size();i++){
                            ComboHit const& chit = sts._straw_chits[i];

                            double StartDOCA = DriftFitUtils::GetTestDOCA(chit, st.FitParams.A0,st.FitParams.A1, st.FitParams.B0, st.FitParams.B1,tracker);


                            double DOCA = DriftFitUtils::GetTestDOCA(chit, st.MinuitParams.A0,st.MinuitParams.A1, st.MinuitParams.B0, st.MinuitParams.B1,  tracker);

                           int RecoAmbig = DriftFitUtils::GetAmbig(chit, st.MinuitParams.A0,st.MinuitParams.A1, st.MinuitParams.B0, st.MinuitParams.B1, tracker);

                            _FullFitEndDOCAs->Fill(DOCA);
                            _StartDOCAs->Fill(StartDOCA);
                           if(_mcdiag){
                                trueinfo = FitMC(_mcdigis);
                                trueinfo = FillDriftMC(chit, RecoAmbig, trueinfo, tracker);
                                if(DriftFitUtils::GetTestDOCA(chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(), tracker)<2.5 and DOCA<2.5 and abs(trueinfo.TrueFitEquation.Pos.X() ) < 5000 and abs(trueinfo.TrueFitEquation.Pos.Y())<5000 and abs(trueinfo.TrueFitEquation.Dir.X())<5 and abs(trueinfo.TrueFitEquation.Dir.Y())<5){

                                        outputfile<<DriftFitUtils::GetTestDOCA(chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(), tracker)<<","<<DOCA<<","<<RecoAmbig/DriftFitUtils::GetAmbig(chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(),  tracker)<<","<<trueinfo.TrueTheta<<","<<st.MinuitParams.A0<<","<<st.MinuitParams.A1
<<","<<st.MinuitParams.B0<<","<<st.MinuitParams.B1<<","<<trueinfo.TrueFitEquation.Pos.X()<<","<<trueinfo.TrueFitEquation.Dir.X()
<<","<<trueinfo.TrueFitEquation.Pos.Y()<<","<<trueinfo.TrueFitEquation.Dir.Y()<<endl;
                        }
                        }
                }


                for(auto const& tseed : *_coscol) {
                        TrkFitFlag const& status = tseed._status;
                        _hitsOK = status.hasAllProperties(TrkFitFlag::hitsOK);
                        _StraightTrackOK = status.hasAllProperties(TrkFitFlag::helixOK);
                        _StraightTrackConverged = status.hasAllProperties(TrkFitFlag::helixConverged);
                        _StraightTrackInit = status.hasAllProperties(TrkFitFlag::circleInit);
                }

                for(size_t ich = 0;ich < _chcol->size(); ++ich){
                        ComboHit const& chit =(*_chcol)[ich];

                        _nhits = chit.nStrawHits();
                        _nsh = chit.nStrawHits();
                        panels.push_back(chit.strawId().panel());
                        planes.push_back(chit.strawId().plane());
                        stations.push_back(chit.strawId().station());

                        _hit_time = chit.time();
                        _hit_drift_time = chit.driftTime();
                        _dt =  _hit_time - _cluster_time;
                        }

                _n_panels = std::set<float>( panels.begin(), panels.end() ).size();
                _n_planes = std::set<float>( planes.begin(), planes.end() ).size();
                _n_stations = std::set<float>( stations.begin(), stations.end() ).size();

      }
      _cosmic_analysis->Fill();
      if(_mcdiag){
              cout<<"true + "<< _AMBIG->GetBinContent(2,2)/_AMBIG->Integral()<<endl;
              cout<<"true - "<< _AMBIG->GetBinContent(1,1)/_AMBIG->Integral()<<endl;
              cout<<"false + "<< _AMBIG->GetBinContent(2,1)/_AMBIG->Integral()<<endl;
              cout<<"false - "<< _AMBIG->GetBinContent(1,2)/_AMBIG->Integral()<<endl;

      }

     }

void CosmicAnalyzer::endJob() {
        if(_diag and _mcdiag) outputfile.close();
}

CosmicTrackMCInfo CosmicAnalyzer::FitMC(const StrawDigiMCCollection*& _mcdigis){
        ::BuildLinearFitMatrixSums S;
        CosmicTrackMCInfo TrackTrueInfo;

        StrawDigiMC hitP1;
        StrawDigiMC first = (*_mcdigis)[0];

        //Get StepPointMC:
        auto const& spmcp0= first.strawGasStep(StrawEnd::cal);
        XYZVectorF pos0(spmcp0->position().x(), spmcp0->position().y(), spmcp0->position().z());
        XYZVectorF dir0(spmcp0->momentum().x(), spmcp0->momentum().y(), spmcp0->momentum().z());

        for(size_t ich = 0;ich < _mcdigis->size(); ++ich){
            hitP1 = (*_mcdigis)[ich];

            //Get StepPointMC:
            auto const& spmcp = hitP1.strawGasStep(StrawEnd::cal);
            XYZVectorF posN(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());

            //Use Step Point MC direction as the True Axes:
            XYZVectorF ZPrime = spmcp->momentum().unit();

            //Store True Track details:
            TrackAxes TrueAxes = ParametricFit::GetTrackAxes(ZPrime);
            TrackTrueInfo.TrueTrackCoordSystem = (TrueAxes);

            //Apply routine to the True Tracks (for validation):
            XYZVectorF point(posN.x(), posN.y(), posN.z());
            XYZVectorF X(1,0,0);
            XYZVectorF Y(0,1,0);
            XYZVectorF Z(0,0,1);
            S.addPoint( point, X,Y,Z, 1,1);

        }

        TrackParams RawTrueParams(S.GetAlphaX()[0][0], S.GetAlphaX()[1][0], S.GetAlphaY()[0][0], S.GetAlphaY()[1][0]);

        XYZVectorF TruePos(S.GetAlphaX()[0][0], S.GetAlphaY()[0][0], 0);

        XYZVectorF TrueDir(S.GetAlphaX()[1][0], S.GetAlphaY()[1][0], 1);
        TrueDir = TrueDir.Unit();
        TrueDir = TrueDir/TrueDir.Z();

        pos0.SetX(pos0.X()-(dir0.X()*pos0.Z()/dir0.Z()));
        pos0.SetY(pos0.Y()-(dir0.Y()*pos0.Z()/dir0.Z()));
        pos0.SetZ(pos0.Z()-(dir0.Z()*pos0.Z()/dir0.Z()));
        dir0 = dir0/dir0.Z();

        TrackEquation TrueTrack(pos0, dir0);

        TrackTrueInfo.TrueFitEquation = (TrueTrack);
        TrackTrueInfo.TruePhi =(atan(TrueDir.y()/TrueDir.x()));
        TrackTrueInfo.TrueTheta = (acos(TrueDir.x()/sqrt(TrueDir.Mag2())));

        return TrackTrueInfo;
     }

  CosmicTrackMCInfo CosmicAnalyzer::FillDriftMC(ComboHit const& chit, double RecoAmbig, CosmicTrackMCInfo info, const Tracker* tracker){

        double true_doca = DriftFitUtils::GetTestDOCA(chit, info.TrueFitEquation.Pos.X(), info.TrueFitEquation.Dir.X(), info.TrueFitEquation.Pos.Y(),info.TrueFitEquation.Dir.Y(),  tracker);
        double trueambig = DriftFitUtils::GetAmbig(chit, info.TrueFitEquation.Pos.X(), info.TrueFitEquation.Dir.X(), info.TrueFitEquation.Pos.Y(),info.TrueFitEquation.Dir.Y(),  tracker);
        //double true_time_residual = DriftFitUtils::TimeResidual(true_doca, _srep, sts._t0.t0(), chit);
        info.Ambig.push_back(trueambig);
        info.TrueDOCA.push_back(true_doca);
        //info.TrueTimeResiduals.push_back(true_time_residual);
        _TrueDOCAs->Fill(true_doca);
        //_TrueTimeResiduals->Fill(true_time_residual);
        _AMBIG->Fill(RecoAmbig, trueambig);
        return info;
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

                _mcdigis=0;
                auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
                _mcdigis = mcdH.product();
        }
        return _chcol != 0 && _tccol!=0 && _coscol !=0 && (_mcdigis != 0 || !_mcdiag);
       }

}

using mu2e::CosmicAnalyzer;
DEFINE_ART_MODULE(CosmicAnalyzer)
