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
      TH1F* _chisq_plot;
      TH1F* _chisq_ndf_plot;
      TH1F* _total_residuals;
      TH1F* _total_pulls;
      TH1F* _hiterrs;
      
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

      Float_t _hit_time, _hit_drift_time, _cluster_time;
	
      //Flags:

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
        
        //Extra histograms for Fit Diags:
        _chisq_plot = tfs->make<TH1F>("chisq_plot","chisq_plot" ,50,0, 10);
	_chisq_ndf_plot = tfs->make<TH1F>("chisq_ndf_plot","chisq_ndf_plot" ,20,0, 0.5);
	_chisq_plot->GetXaxis()->SetTitle("#Chi^{2}");
	_chisq_ndf_plot->GetXaxis()->SetTitle("#Chi^{2}/N");
	
        
        _total_residuals = tfs->make<TH1F>("Residuals ","Residuals " ,20,-50,50);
	_total_residuals->GetXaxis()->SetTitle("Residual  [mm]");
	_total_residuals->SetStats();

	_total_pulls = tfs->make<TH1F>("Pull ","Pull " ,25,-1, 1);
	_total_pulls->GetXaxis()->SetTitle("Pull ");
	_total_pulls->SetStats();
       
        _hiterrs = tfs->make<TH1F>("hiterr","hiterr" ,10,0, 100);
	_hiterrs->GetXaxis()->SetTitle("Total Hit Error [mm]");
	_hiterrs->SetStats();
        
        
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
        	CosmicTrackSeed sts =(*_coscol)[ist];
		CosmicTrack st = sts._track;
		//TrkFitFlag const& status = sts._status;
		std::vector<int> panels, planes, stations;
                
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
			XYZVec track_dir = st.get_track_direction();
			
			double werr_mag = chit.wireRes(); //hit major error axis 
      			double terr_mag = chit.transRes(); //hit minor error axis 
			XYZVec const& wdir = chit.wdir();//direction along wire
      			XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire
			XYZVec major_axis = werr_mag*wdir;
      			XYZVec minor_axis = terr_mag*wtdir;
	               
			double hit_error = sqrt(major_axis.Dot(track_dir)*major_axis.Dot(track_dir)+minor_axis.Dot(track_dir)*minor_axis.Dot(track_dir));
			
	        
		//-----------Fill Hist Details:----------//
			
			
			for(size_t i=0; i< st.get_fit_residuals().size();i++)			{
				   
				    //double pull = st.get_fit_residuals()[ich]/st.get_fit_residual_errors()[ich];
		                    _total_residuals->Fill(st.get_fit_residuals()[ich]);
			           // _total_pulls->Fill(pull);
				    
			            
                          }
		                _hiterrs->Fill(hit_error);
            		
			
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
    

    


