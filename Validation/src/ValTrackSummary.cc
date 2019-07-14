
#include "art_root_io/TFileDirectory.h"
#include "Validation/inc/ValTrackSummary.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TMath.h"

int mu2e::ValTrackSummary::declare(art::TFileDirectory tfs) {

  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);

  _hNTr = tfs.make<TH1D>( "NTracks", "N Tracks", 11, -0.5, 10.5);
  _hNState = tfs.make<TH1D>( "NState", "N States", 21, -0.5, 20.0);
  _hp = tfs.make<TH1D>( "p", "p", 100, 0., 110.);
  _hpce = tfs.make<TH1D>( "pce", "p CE", 100, 95.0, 110.);
  _hD0 = tfs.make<TH1D>( "d0", "d0", 100, -100., 100.);
  _hPhi0 = tfs.make<TH1D>( "phi0", "phi0", 100, -TMath::Pi(), TMath::Pi());
  _hOmega = tfs.make<TH1D>( "omega", "omega", 100, 0.0, 0.02);
  _hZ0 = tfs.make<TH1D>( "Z0", "Z0", 100, -1000.0, 1000.0);
  _hTan = tfs.make<TH1D>( "tanDip", "tanDip", 100, 0.4, 1.8);
  _hT0 = tfs.make<TH1D>( "t0", "t0", 100, 400.0, 1800.0);

  _hNActive = tfs.make<TH1D>( "NActive", "N Active", 101, -0.5, 100.0);
  _hNDof = tfs.make<TH1D>( "NDof", "N DOF", 101, -0.5, 100.0);
  _hChi2N = tfs.make<TH1D>( "Chi2N", "Chi2/DOF", 100, 0.0, 10.0);
  _hCL = tfs.make<TH1D>( "FitConn", "Fit CL", 100, 0.0, 1.0);
  _hStatus = tfs.make<TH1D>( "Status", "Fit Status", 11, -0.5, 10.0);

  _hCuts = tfs.make<TH1D>( "Cuts", "Cut series", 8, 0.5, 8.5);
  int ibin=1;
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"All CE"); // bin 1, first visible
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"MC Selection");
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"KF Track fit"); // 3
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"Livegate");   // 5
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection"); 
  _hCuts->GetXaxis()->SetBinLabel(ibin++,"Momentum window");  // 8
  _hCuts->SetMinimum(0.0);
  _hPRes = tfs.make<TH1D>( "PRes", "R resolution", 200, -5.0, 3.0);

  return 0;
}

int mu2e::ValTrackSummary::fill(const mu2e::TrackSummaryCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(2.0);

  // p of highest momentum electron SimParticle with good tanDip
  double p_mc = mcTrkP(event);

  // the first of the cut series, number of events
  _hCuts->Fill(1.0);
  // MC CE found
  if(p_mc>100.0) _hCuts->Fill(2.0);

  _hNTr->Fill(coll.size());
  for(auto trk : coll) {
    //    TrackSummary const& trk = t;

    const std::vector<mu2e::TrackSummary::TrackStateAtPoint>& statev 
                                         = trk.states();
    _hNState->Fill(statev.size());

    // histogram first track state only
    if (!statev.empty()) {
      auto state = statev.begin();
      auto helix = state->helix();
      _hD0->Fill(helix.d0());
      double d0 = helix.d0();
      _hPhi0->Fill(helix.phi0());
      _hOmega->Fill(helix.omega());
      double om = helix.omega();
      _hZ0->Fill(helix.z0());
      _hTan->Fill(helix.tanDip());
      double td = helix.tanDip();
      double p = state->momentum().mag();
      _hp->Fill(p);
      _hpce->Fill(p);
      bool d0cut = (d0<105 && d0>-80 && (d0+2/om)>450 && (d0+2/om)<680);

      if(p_mc>100.0) {
	if(trk.fitstatus()>0) {
	  _hCuts->Fill(3.0);
	  if(trk.fitcon()>0.002) {
	    _hCuts->Fill(4.0);
	    if(trk.t0()>700.0 && trk.t0()<1695.0) {
	      _hCuts->Fill(5.0);
	      if(td>0.577&&td<1.0) {
		_hCuts->Fill(6.0);
		if(d0cut) {
		  _hCuts->Fill(7.0);
		  if(p>103.75 && p< 105.0) {
		    _hCuts->Fill(8.0);
		    _hPRes->Fill(p-p_mc);
		  }
		}
	      }
	    }
	  }
	}
      } // if(p_mc>100.0) 

    }  // statev not empty

    _hT0->Fill(trk.t0());
    _hNActive->Fill(trk.nactive());
    _hNDof->Fill(trk.ndof());
    double cn = (trk.ndof()>0 ? trk.chi2()/trk.ndof() : -0.1);
    _hChi2N->Fill(cn);
    _hCL->Fill(trk.fitcon());
    _hStatus->Fill(trk.fitstatus());

  }

  return 0;
}

double mu2e::ValTrackSummary::mcTrkP(art::Event const& event) {

  art::Handle<SimParticleCollection> simsHandle;
  event.getByLabel("g4run",simsHandle);
  if (!simsHandle.isValid()) return -1.0;

  SimParticleCollection const& simpcoll = *simsHandle;


  art::Handle<StepPointMCCollection> mcVDstepsHandle;
  event.getByLabel("g4run","virtualdetector",mcVDstepsHandle);
  if (!mcVDstepsHandle.isValid()) return -1.0;
  StepPointMCCollection const& spmccoll = *mcVDstepsHandle;


  /*

	if(vids.size() == 0 ||  (imcs->trackId() == trkid && find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
	  steps.push_back(imcs);
	}
      }
      // sort these in time
      sort(steps.begin(),steps.end(),timecomp());


  }

  for( MCStepItr imcs =mcsteps->begin();imcs!= mcsteps->end();imcs++){
	if(vids.size() == 0 ||  (imcs->trackId() == trkid && find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
	  steps.push_back(imcs);
	}
      }
      // sort these in time
      sort(steps.begin(),steps.end(),timecomp());



const StepPointMCCollection *_mcsteps, *_mcvdsteps

     _entvids.push_back(VirtualDetectorId::TT_FrontHollow);
    _entvids.push_back(VirtualDetectorId::TT_FrontPA); 

     if(entsteps.size() > 0 && vdg->exist(entsteps[0]->volumeId()))
	fillTrkInfoMCStep(entsteps.front(),_mcentinfo);




  KalDiag::findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, vector<int> const& vids,
  vector<MCStepItr>& steps) {
    steps.clear();
    if(mcsteps != 0){
      // Loop over the step points, and find the one corresponding to the given detector
      for( MCStepItr imcs =mcsteps->begin();imcs!= mcsteps->end();imcs++){
	if(vids.size() == 0 ||  (imcs->trackId() == trkid && find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
	  steps.push_back(imcs);
	}
      }
      // sort these in time
      sort(steps.begin(),steps.end(),timecomp());
    }
  }

  */


  double p = -1.0;
  double td,p0;
  for(auto sp : simpcoll) {
    const mu2e::SimParticle& part = sp.second;
    p0 = 0.0;
    td = 0.0;
    if(part.pdgId()==11) {

      bool found = false;
      double t0 = 1e10;
      CLHEP::Hep3Vector pv;

      for(auto step : spmccoll) {
	if( step.trackId() == part.id() && 
	    (step.volumeId()==VirtualDetectorId::TT_FrontHollow ||
	     step.volumeId()==VirtualDetectorId::TT_FrontPA) ) {
	  if(step.time()<t0) {
	    t0 = step.time();
	    found = true;
	    pv = step.momentum();
	  }
	}
      }
      if(found) {
	p0 = pv.mag();
	double pz = pv.z();
	double pt = pv.perp();
	if(pt!=0.0) {
	  td = pz/pt;
	} else {
	  td = 1e14;
	}

	if(p0>p && p0>100.0 && td>0.527 && td<1.2) {
	  p = p0;
	}
      } // if found
    } // if electron
  }
  return p;

}
