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
      TH1F* _niters;
      
      TH1F* _chisq_ndf_plot_final;
      TH1F* _chisq_ndf_plot_finalX;
      TH1F* _chisq_ndf_plot_finalY;
      TH1F* _total_residualsX_final;
      TH1F* _total_pullsX_final;
      TH1F* _total_residualsY_final;
      TH1F* _total_pullsY_final;
      TH1F* _hiterrs_final;
      TH2F* _chi2_v_angle;
      TH2F* _chi2_v_NHits;
      TH2F* _chi2_v_length;
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
       
	_chisq_ndf_plot_init = tfs->make<TH1F>("init chisq_ndf_plot","init chisq_ndf_plot" ,200,0, 1000);
	_chisq_ndf_plot_init->GetXaxis()->SetTitle("Init. #Chi^{2}/N");
	
	_chisq_ndf_plot_final = tfs->make<TH1F>("final chisq_ndf_plot","final chisq_ndf_plot" ,200,0,1000);
	_chisq_ndf_plot_final->GetXaxis()->SetTitle("Final #Chi^{2}/N");
	
	_chisq_ndf_plot_finalX = tfs->make<TH1F>("final chisq_ndf_plot X''","final chisq_ndf_plot X''" ,200,0,1000);
	_chisq_ndf_plot_finalX->GetXaxis()->SetTitle("Final X '' #Chi^{2}/N");
	
        _chisq_ndf_plot_finalY = tfs->make<TH1F>("final chisq_ndf_plot Y''","final chisq_ndf_plot Y''" ,200,0,1000);
	_chisq_ndf_plot_finalY->GetXaxis()->SetTitle("Final Y '' #Chi^{2}/N");
	
        _total_residualsX_init = tfs->make<TH1F>("Initial Residuals X'' ","Initial Residuals X'' " ,200,-1000,1000);
	_total_residualsX_init->GetXaxis()->SetTitle("Initial Residual X'' [mm]");
	_total_residualsX_init->SetStats();
	
	_total_residualsY_init = tfs->make<TH1F>("Initial Residuals Y''","Initial Residuals Y'' " ,200,-1000,1000);
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
	
	 _total_residualsX_final = tfs->make<TH1F>("Final Residuals X'' #chi^{2}/dof all ","Final Residuals X''#chi^{2}/dof all  " ,200,-1000,1000);
	_total_residualsX_final->GetXaxis()->SetTitle("Residual X'' [mm]");
	_total_residualsX_final->SetStats();
	
	_total_residualsY_final = tfs->make<TH1F>("Final Residuals Y''#chi^{2}/dof all ","Fi nalResiduals Y'' #chi^{2}/dof all " ,200,-1000,1000);
	_total_residualsY_final->GetXaxis()->SetTitle("Final Residual Y'' [mm]");
	_total_residualsY_final->SetStats();

	_total_pullsX_final = tfs->make<TH1F>("Final Pull X'' #chi^{2}/dof all ","Final Pull X''#chi^{2}/dof all " ,200,-50, 50);
	_total_pullsX_final->GetXaxis()->SetTitle("Final Pull X''");
	_total_pullsX_final->SetStats();
	
	_total_pullsY_final = tfs->make<TH1F>("Final Pull Y''#chi^{2}/dof all ","Finla Pull Y''#chi^{2}/dof all " ,200,-50, 50);
	_total_pullsY_final->GetXaxis()->SetTitle("Final Pull Y''");
	_total_pullsY_final->SetStats();
       
        _hiterrs_final = tfs->make<TH1F>("Final Total Hit Error #chi^{2}/dof all ","Final Total Hit Error #chi^{2}/dof all " ,50,0, 100);
	_hiterrs_final->GetXaxis()->SetTitle("Hit Error in Track Frame [mm]");
	_hiterrs_final->SetStats();
	
	_a0 = tfs->make<TH1F>("Track Parameter A0 #chi^{2}/dof all ","Track Parameter A0 #chi^{2}/dof all " ,50,-1000, 1000);
	_a0->GetXaxis()->SetTitle("Track Parameter A0");
	_a0->SetStats();
	
	_b0 = tfs->make<TH1F>("Track Parameter B0 #chi^{2}/dof all ","Track Parameter B0 #chi^{2}/dof all " ,50,-1000, 1000);
	_b0->GetXaxis()->SetTitle("Track Parameter B0");
	_b0->SetStats();
	
	
	_a1 = tfs->make<TH1F>("Track Parameter A1 #chi^{2}/dof all ","Track Parameter A1 #chi^{2}/dof all " ,100,-10, 10);
	_a1->GetXaxis()->SetTitle("Track Parameter A1");
	_a1->SetStats();
	
	_b1 = tfs->make<TH1F>("Track Parameter B1 #chi^{2}/dof all ","Track Parameter B1 #chi^{2}/dof all " ,100,-10, 10);
	_b1->GetXaxis()->SetTitle("Track Parameter B1");
	_b1->SetStats();
	
	_niters  = tfs->make<TH1F>("Number of Iterations Unitl Converged #chi^{2}/dof all ","Number of Iterations Unitl Converged #chi^{2}/dof all " ,100,0, 100);
	_niters->GetXaxis()->SetTitle("Number of Iterations Unitl Converged");
	_niters->SetStats();
	
	_chi2_v_angle  = tfs->make<TH2F>("chi2 v angle #chi^{2}/dof all ","chi2 v angle #chi^{2}/dof all " ,1000, 0, 1000, 100,0, 3.1415);
	_chi2_v_angle->GetYaxis()->SetTitle("Angle");
	_chi2_v_angle->SetStats();
	
	_chi2_v_NHits  = tfs->make<TH2F>("chi2 v nhits #chi^{2}/dof all ","chi2 v nits #chi^{2}/dof all " ,50, 0, 10, 50,0, 50);
	_chi2_v_NHits ->GetYaxis()->SetTitle("Nhits");
	_chi2_v_NHits ->SetStats();
	
	_chi2_v_length = tfs->make<TH2F>("Res Y'' v hit number  ","chi2 v hit number  " ,100, 0, 100, 100, 0, 20);
	_chi2_v_length ->GetYaxis()->SetTitle("res Y'' [mm]");
	_chi2_v_length ->SetStats();
	
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
        //loop over tracks"
        for(size_t ist = 0;ist < _coscol->size(); ++ist){
        	n_analyze+=1;
        	std::cout<<"number analyzed "<<n_analyze<<std::endl;
        	CosmicTrackSeed sts =(*_coscol)[ist];
		CosmicTrack st = sts._track;
		//TrkFitFlag const& status = sts._status;
		std::vector<int> panels, planes, stations;
                _chisq_ndf_plot_init->Fill(st.get_initchisq_dof());
                _chisq_ndf_plot_final->Fill(st.get_finalchisq_dof());
                
                _chisq_ndf_plot_finalX->Fill(st.get_finalchisq_dofX());
                _chisq_ndf_plot_finalY->Fill(st.get_finalchisq_dofY());
                
                double angle =st.getZPrime().Theta();
                _chi2_v_angle->Fill(st.get_finalchisq_dof(),angle);
                _chi2_v_NHits ->Fill( st.get_finalchisq_dof(), st.get_N());
                
               
                _a1->Fill(st.get_track_parameters()[1]);
                _b1->Fill(st.get_track_parameters()[3]);
                _a0->Fill(st.get_track_parameters()[0]);
                _b0->Fill(st.get_track_parameters()[2]);
               
                for(size_t i=0; i< st.get_iter().size();i++){
                    _niters->Fill(st.get_iter()[i]);
                    
                }
                //-----------Fill Hist Details:----------//
		for(size_t i=0; i< st.get_init_hit_errorsTotal().size();i++){
		    _hiterrs_init->Fill(st.get_init_hit_errorsTotal()[i]);
		}
		for(size_t i=0; i< st.get_init_fit_residualsX().size();i++){
		    double pullX = st.get_init_fit_residualsX()[i]/st.get_init_fit_residual_errorsX()[i];
                    _total_residualsX_init->Fill(st.get_init_fit_residualsX()[i]);             
	            _total_pullsX_init->Fill(pullX);
	            
	        }
	        for(size_t i=0; i< st.get_init_fit_residualsY().size();i++){
	            double pullY = st.get_init_fit_residualsY()[i]/st.get_init_fit_residual_errorsY()[i];
                    _total_residualsY_init->Fill(st.get_init_fit_residualsY()[i]);             
	            _total_pullsY_init->Fill(pullY);
	        }   
	         
	        
                //-----------Fill Hist Details:----------//
		for(size_t i=0; i< st.get_final_hit_errorsTotal().size();i++){
		    _hiterrs_final->Fill(st.get_final_hit_errorsTotal()[i]);
		}
		for(size_t i=0; i< st.get_final_fit_residualsX().size();i++){
		    double pullX = st.get_final_fit_residualsX()[i]/st.get_final_fit_residual_errorsX()[i];
                    _total_residualsX_final->Fill(st.get_final_fit_residualsX()[i]);             
	            _total_pullsX_final->Fill(pullX);
	        }
	        for(size_t i=0; i< st.get_final_fit_residualsY().size();i++){
	            double pullY = st.get_final_fit_residualsY()[i]/st.get_final_fit_residual_errorsY()[i];
                    _total_residualsY_final->Fill(st.get_final_fit_residualsY()[i]);             
	            _total_pullsY_final->Fill(pullY);
	            _chi2_v_length ->Fill(i,st.get_final_fit_residualsY()[i]);
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
			
		    
                //----------------Get panels/planes/stations per track:------------------//
                _n_panels = std::set<double>( panels.begin(), panels.end() ).size();
		_n_planes = std::set<double>( planes.begin(), planes.end() ).size();
		_n_stations = std::set<double>( stations.begin(), stations.end() ).size();
	 
		}//endS ST
	        
      }//end analyze
     
      _cosmic_analysis->Fill();
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
    

    


