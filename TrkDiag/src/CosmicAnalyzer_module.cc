#define _USE_MATH_DEFINES

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "TrkReco/inc/StraightTrackFit.hh"
#include "TrkPatRec/inc/StraightTrackFinder_types.hh"
#include "TrkReco/inc/StraightTrackFinderData.hh"

//Mu2e Data Prods:
#include "DataProducts/inc/threevec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"
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
      TH1F* _chisq_plot;
      TH1F* _chisq_ndf_plot;
      TH1F* _total_residuals_XY;
      TH1F* _total_pulls_XY;
      TH1F* _hiterrs_XY;
      TH1F* _fit_phi_angle;
      TH1F* _fit_theta_angle;
      TH1F* _reco_eff_v_angle;
      TH1F* _mc_theta_angle;
      TH1F* _mc_phi_angle;
      TH2F* _pull_v_cosine_XY;
      TH1F* _c;
      TH1F* _m0error;
      TH1F* _m1error;
      TH1F* _cerror;
      
      //Fit Parameters:
      Float_t _chisq, _chisq_ndf, _c0, _m0,_m1, _c_err, _m_0_err,_m_1_err;
      Float_t _resid;
      // add event id
      Int_t _evt; 

      // Event object Tags
      art::InputTag   _chtag;//combo
      art::InputTag   _tctag;//timeclusters
      art::InputTag   _sttag;//Striaght tracks
      art::InputTag   _mcdigistag; //MC digis
      // Event data cache
      const ComboHitCollection* _chcol;
      const StraightTrackSeedCollection* _stcol;
      const TimeClusterCollection* _tccol;
      const StrawDigiMCCollection* _mcdigis;
      //Numbers:
      Int_t _nsh, _nch; // # associated straw hits / event
      Int_t     _ntc; // # clusters/event
      Int_t _nhits, _nused; // # hits used
      Int_t _n_panels; // # panels
      Int_t _n_stations; // # stations
      Int_t _n_planes; // # stations

      Float_t _hit_time, _hit_drift_time, _cluster_time;
	
      //Flags:

      //offsets for MC
      SimParticleTimeOffset _toff;
      
      //Add in details of the tracker segments positions. This will be used for alignment later on:

      Int_t _strawid; // strawid info

      // per-hit diagnostics
      //vector<ComboHitInfo> _chinfo;
      vector<ComboHitInfoMC> _chinfomc;

      bool findData(const art::Event& evt);
    };

    CosmicAnalyzer::CosmicAnalyzer(fhicl::ParameterSet const& pset) :
	art::EDAnalyzer(pset),
	_diag		(pset.get<int>("diagLevel",1)),
	_mcdiag		(pset.get<bool>("mcdiag",true)),
	_chtag		(pset.get<art::InputTag>("ComboHitCollection")),
	_tctag		(pset.get<art::InputTag>("TimeClusterCollection")),
	_sttag		(pset.get<art::InputTag>("StraightTrackSeedCollection")),
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
       _cosmic_analysis->Branch("c0", &_c, "c/F");
       _cosmic_analysis->Branch("c0_error", &_c_err, "c0_error/F");
       _cosmic_analysis->Branch("m0_error", &_m_0_err, "m0_error/F");
       _cosmic_analysis->Branch("m1_error", &_m_1_err, "m1_error/F");
       _cosmic_analysis->Branch("m0", &_m0, "m0/F");
       _cosmic_analysis->Branch("m1", &_m1, "m1/F");
       _cosmic_analysis->Branch("hit_time", &_hit_time, "hit_time/F");
       _cosmic_analysis->Branch("hit_drit_time", &_hit_drift_time, "hit_drift_time/F");
       _cosmic_analysis->Branch("cluster_time", &_cluster_time, "cluster_time/F");
        
        //Extra histograms for Fit Diags:
        _chisq_plot = tfs->make<TH1F>("chisq_plot","chisq_plot" ,50,0, 10);
	_chisq_ndf_plot = tfs->make<TH1F>("chisq_ndf_plot","chisq_ndf_plot" ,20,0, 0.5);
	_chisq_plot->GetXaxis()->SetTitle("#Chi^{2}");
	_chisq_ndf_plot->GetXaxis()->SetTitle("#Chi^{2}/N");
	
        
        _total_residuals_XY = tfs->make<TH1F>("Residuals ","Residuals " ,100,-100,100);
	_total_residuals_XY->GetXaxis()->SetTitle("Residual  [mm]");
	_total_residuals_XY->SetStats();

	_total_pulls_XY = tfs->make<TH1F>("Pull ","Pull " ,100,-1, 1);
	_total_pulls_XY->GetXaxis()->SetTitle("Pull ");
	_total_pulls_XY->SetStats();

        _hiterrs_XY = tfs->make<TH1F>("hiterr","hiterr" ,100,0, 1000);
	_hiterrs_XY->GetXaxis()->SetTitle("Total Hit Error [mm]");
	_hiterrs_XY->SetStats();
        
        _pull_v_cosine_XY = tfs->make<TH2F>("pull v cos XY","pull v cos" ,100,-1,1, 100, -4,4);
	_pull_v_cosine_XY->GetXaxis()->SetTitle("cosine");
	_pull_v_cosine_XY->GetYaxis()->SetTitle("pull XY");
	_pull_v_cosine_XY->SetStats();

	_fit_phi_angle=tfs->make<TH1F>("#phi_angle","#phi_angle" ,20, 0,M_PI);
	_fit_phi_angle->GetXaxis()->SetTitle("#phi_{fit} of Track [rad]");
	_fit_phi_angle->SetStats();

	_fit_theta_angle=tfs->make<TH1F>("#theta_angle","polar_angle" ,20, 0,M_PI);
	_fit_theta_angle->GetXaxis()->SetTitle("#theta_{fit} of track [rad]");
	_fit_theta_angle->SetStats();

	_mc_theta_angle=tfs->make<TH1F>("mc_true_theta","mc true theta" ,20, 0,M_PI);
	_mc_theta_angle->GetXaxis()->SetTitle("#theta_{true} of muon [rad]");
	_mc_theta_angle->SetStats();

	_mc_phi_angle=tfs->make<TH1F>("mc_true_phi","mc true phi" ,20, 0,M_PI);
	_mc_phi_angle->GetXaxis()->SetTitle("#phi_{true} of muon [rad]");
	_mc_phi_angle->SetStats();

	_reco_eff_v_angle=tfs->make<TH1F>("Recon_theta","Recon_theta" ,20, 0,M_PI);//,20, 0,10);
	_reco_eff_v_angle->GetXaxis()->SetTitle("#theta_{true}");
	_reco_eff_v_angle->GetYaxis()->SetTitle("#N_{reco}/#N_{true}");
	_reco_eff_v_angle->SetStats(0);

	_m0error=tfs->make<TH1F>("m0error","m0error" ,100, 0, 10);
 	_m0error->GetXaxis()->SetTitle("Error on m0");
	_m0error->SetStats();

	_m1error=tfs->make<TH1F>("m1error","m1error" ,100, 0, 10);
 	_m1error->GetXaxis()->SetTitle("Error on m1");
	_m1error->SetStats();

	_cerror=tfs->make<TH1F>("cerror","cerror" ,100, 0, 1000);
 	_cerror->GetXaxis()->SetTitle("Error on c");
	_cerror->SetStats();

	_c=tfs->make<TH1F>("c0","c0" ,100, -10000, 10000);
 	_c->GetXaxis()->SetTitle("c0");
	_c->SetStats();
        
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
        for(size_t ist = 0;ist < _stcol->size(); ++ist){
        	StraightTrackSeed sts =(*_stcol)[ist];
		StraightTrack st = sts._track;
		//TrkFitFlag const& status = sts._status;
		std::vector<int> panels, planes, stations;
                //--------Fill Fit Paramters ----------//
	        _c0 = st.get_c_0();
		_c_err = st.get_c_0_err();
		_m_0_err = st.get_m_0_err();
		
 		_m_1_err = st.get_m_1_err();
		_m0 = st.get_m_0();
		_m1 = st.get_m_1();
		_c->Fill(st.get_c_0());
		_m0error->Fill(st.get_m_0_err());
                _m1error->Fill(st.get_m_1_err());
		_cerror->Fill(st.get_c_0_err());
		if(st.get_chisq_dof()>0){
			_chisq_ndf_plot->Fill(st.get_chisq_dof());
			_chisq_plot->Fill(st.get_chisq());
		}


		//--------Get Track Direction------------//
		 XYZVec track_dir( st.get_m_0(), 1, st.get_m_1());
                //XYZVec track_dir( st.get_m_0(), st.get_m_1(), st.get_m_2());
                //-----------------ANGLES-----------------//

		double phi = track_dir.Phi();//TVector3 def
		double theta = track_dir.Theta(); //TVector3 def
		double mc_phi_sum = 0;//For MC
		double mc_theta_sum = 0;//For MC
		

		for(size_t ich = 0;ich < _chcol->size(); ++ich){
                        ComboHit const& chit =(*_chcol)[ich];

                //-----------Fill diag details:----------//
                        _nhits = chit.nStrawHits();
                        _nsh = chit.nStrawHits();
                        panels.push_back(chit.strawId().panel());
		        planes.push_back(chit.strawId().plane());
			stations.push_back(chit.strawId().station());

		//-----------Hit details:---------------//
			//XYZVec cpos = chit.pos() - chit.wireDist()*chit.wdir();
		        _hit_time = chit.time();
			_hit_drift_time = chit.driftTime();
			

			double werr_mag = chit.wireRes(); //hit major error axis 
      			double terr_mag = chit.transRes(); //hit minor error axis 
			XYZVec const& wdir = chit.wdir();//direction along wire
      			XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire
			XYZVec major_axis = werr_mag*wdir;
      			XYZVec minor_axis = terr_mag*wtdir;
	                
      			
			 //double hit_error = sqrt((major_axis.Dot(track_dir)*major_axis.Dot(track_dir)))+(minor_axis.Dot(track_dir)*minor_axis.Dot(track_dir));
			double hit_error = sqrt(major_axis.Cross(track_dir).mag2()+minor_axis.Cross(track_dir).mag2());
			
	        
		//-----------Fill Hist Details:----------//
			
			
			for(size_t i=0; i< st.get_fit_residuals().size();i++)			{
				    double pull = st.get_fit_residuals()[i]/st.get_fit_residual_errors()[i];
		                    _total_residuals_XY->Fill(st.get_fit_residuals()[i]);
			            _total_pulls_XY->Fill(pull);
				    _pull_v_cosine_XY->Fill( wdir.Dot(track_dir), pull);
			            
                          }
		                _hiterrs_XY->Fill(hit_error);
            
			
	
		//------------MC INFO------------------//
		
		if(_mcdiag){
                   
				std::vector<StrawDigiIndex> shids;
	 			_chcol->fillStrawDigiIndices(event,ich,shids);
				if(shids.size() != chit.nStrawHits()){
	  				throw cet::exception("DIAG")<<"mu2e::ComboHitDiag: invalid ComboHit" << std::endl;
		                }
				StrawDigiMC const& mcd1 = _mcdigis->at(shids[0]);
				art::Ptr<StepPointMC> const& spmcp = mcd1.stepPointMC(StrawEnd::cal);
		                double mcposx = spmcp->position().x();
				double mcposy = spmcp->position().y();
				double mcposz = spmcp->position().z();
				//XYZVec mc_track_dir( mcposy/mcposx, 1, mcposy/mcposz);
					 
				 XYZVec mcpos = XYZVec(mcposx, mcposy,mcposz);
				// double mcpos_perp = sqrt(mcpos.Perp2());	
				 double mcpos_mag = sqrt(mcpos.Mag2());    
		                 //_mc_polar_angle->Fill(acos(mcposz/mcpos_perp));//causes infinitiy but right shape
				mc_phi_sum += acos(mcposx/mcpos_mag);
				mc_theta_sum += acos(mcposz/mcpos_mag);
				//_mc_phi_angle->Fill(acos(mcposx/mcpos_mag));
				//_mc_theta_angle->Fill(acos(mcposz/mcpos_mag));
			         
                                   
		}//end MCdiag
				
		}//end CH
		_mc_phi_angle->Fill(mc_phi_sum/_chcol->size());
		_mc_theta_angle->Fill(mc_theta_sum/_chcol->size());
			         	
		_fit_phi_angle->Fill(phi); 
		_fit_theta_angle->Fill(theta); 
                //----------------Get panels/planes/stations per track:------------------//
                _n_panels = std::set<double>( panels.begin(), panels.end() ).size();
		_n_planes = std::set<double>( planes.begin(), planes.end() ).size();
		_n_stations = std::set<double>( stations.begin(), stations.end() ).size();

                
	 		 
		}//endS ST
	        //----------------------------Get RecoEff:-------------------------------//
		if(_fit_phi_angle->GetNbinsX()!= _mc_phi_angle->GetNbinsX()){
  				throw cet::exception("DIAG")<<"WARNING: cant calc reco eff as bin numbers do not match!" << std::endl;
                }
                else{
			int N_polar_bins = _mc_phi_angle->GetNbinsX();
			for(int n=1;n<N_polar_bins; n++){
				
				double frac_recon = _fit_phi_angle->GetBinContent(n)/ _mc_phi_angle->GetBinContent(n);
				_reco_eff_v_angle->SetBinContent(n,frac_recon);
				
		}//end reco eff calc
      }//end analyze

      _cosmic_analysis->Fill();
      }

      bool CosmicAnalyzer::findData(const art::Event& evt){
	_chcol = 0; 
        _tccol = 0;
        _stcol = 0; 
	auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
	_chcol = chH.product();
	auto tcH = evt.getValidHandle<TimeClusterCollection>(_tctag);
	_tccol =tcH.product();
	auto stH = evt.getValidHandle<StraightTrackSeedCollection>(_sttag);
	_stcol =stH.product();
        if(_mcdiag){
           auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
           _mcdigis = mcdH.product();
	   // update time offsets
           _toff.updateMap(evt);
        }
	return _chcol != 0 && _tccol!=0 && _stcol !=0 && (_mcdigis != 0 || !_mcdiag);
       }



}//End Namespace Mu2e

using mu2e::CosmicAnalyzer;
DEFINE_ART_MODULE(CosmicAnalyzer);
    

    


