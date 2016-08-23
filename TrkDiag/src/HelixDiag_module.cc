//
// Helix Fit diagnostics
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// mu2e
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
// TrkDiag
#include "TrkDiag/inc/TrkMCTools.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixCollection.hh"
#include "DataProducts/inc/threevec.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// root
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TArc.h"
#include "TTree.h"
#include "TMarker.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
// C++
#include <functional>
#include <algorithm>
using namespace std; 
using namespace boost::accumulators;
using CLHEP::Hep3Vector;

namespace mu2e {
  class HelixDiag : public art::EDAnalyzer {
    public:

      explicit HelixDiag(fhicl::ParameterSet const& pset);
      virtual ~HelixDiag();

      virtual void beginJob();

      // This is called for each event.
      virtual void analyze(art::Event const& e);

    private:
// config parameters
    int _diag;
// diagnostic histograms
    TH1F *_rdiff, *_fdiff;
    TH1F *_rpull, *_fpull;
// diagnostics
    void plotXY(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitResult const& myhel) const;
    void plotZ(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitResult const& myhel) const;
 
  };

  HelixDiag::HelixDiag(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0))
  {

    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _rdiff = tfs->make<TH1F>("rdiff","Radial difference",100,-100,100.0);
      _rpull = tfs->make<TH1F>("rpull","Radial pull",100,-20.0,20.0);
      _fdiff = tfs->make<TH1F>("fdiff","phi difference",100,-5,5.0);
      _fpull = tfs->make<TH1F>("fpull","Phi pull",100,-20.0,20.0);
    }
  }

  HelixDiag::analyze(art::Event const& evt) {

    if(_diag>1 && plothelix){
      // fill graphs for display if requested
      plotXY(mytrk,xyzp,myhel);
      plotZ(mytrk,xyzp,myhel);
    }
  }

  void
  HelixDiag::plotXY(HelixDef const& mytrk, XYZPVector const& xyzp,
    HelixFitResult const& myhel) const {

    // Check that there are more than 10 CE hits first
    int n_ce_hits = 0;
    for (XYZPVector::const_iterator i_hit = xyzp.begin(); i_hit != xyzp.end(); ++i_hit) {
      if ((*i_hit).conversion()) {
	++n_ce_hits;
      }
    }

    if (n_ce_hits > 10) {
      static unsigned igraph = 0;
      igraph++;
      art::ServiceHandle<art::TFileService> tfs;
      char ce_stereo_used_name[100];
      snprintf(ce_stereo_used_name,100,"ce_stereo_used_shxy%i",igraph);
      char ce_stereo_notused_name[100];
      snprintf(ce_stereo_notused_name,100,"ce_stereo_notused_shxy%i",igraph);
      char ce_notstereo_used_name[100];
      snprintf(ce_notstereo_used_name,100,"ce_notstereo_used_shxy%i",igraph);
      char ce_notstereo_notused_name[100];
      snprintf(ce_notstereo_notused_name,100,"ce_notstereo_notused_shxy%i",igraph);
      char bkg_stereo_used_name[100];
      snprintf(bkg_stereo_used_name,100,"bkg_stereo_used_shxy%i",igraph);
      char bkg_stereo_notused_name[100];
      snprintf(bkg_stereo_notused_name,100,"bkg_stereo_notused_shxy%i",igraph);
      char bkg_notstereo_used_name[100];
      snprintf(bkg_notstereo_used_name,100,"bkg_notstereo_used_shxy%i",igraph);
      char bkg_notstereo_notused_name[100];
      snprintf(bkg_notstereo_notused_name,100,"bkg_notstereo_notused_shxy%i",igraph);
      char title[100];
      snprintf(title,100,"StrawHit XY trk %i;mm;rad",igraph);
      TH2F* ce_stereo_used = tfs->make<TH2F>(ce_stereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* ce_stereo_notused = tfs->make<TH2F>(ce_stereo_notused_name,title,100,-500,500,100,-500,500);
      TH2F* ce_notstereo_used = tfs->make<TH2F>(ce_notstereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* ce_notstereo_notused = tfs->make<TH2F>(ce_notstereo_notused_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_stereo_used = tfs->make<TH2F>(bkg_stereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_stereo_notused = tfs->make<TH2F>(bkg_stereo_notused_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_notstereo_used = tfs->make<TH2F>(bkg_notstereo_used_name,title,100,-500,500,100,-500,500);
      TH2F* bkg_notstereo_notused = tfs->make<TH2F>(bkg_notstereo_notused_name,title,100,-500,500,100,-500,500);

      ce_stereo_used->SetMarkerStyle(kFullTriangleUp);
      ce_stereo_used->SetMarkerColor(kRed);
      ce_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      ce_stereo_notused->SetMarkerColor(kRed);
      ce_notstereo_used->SetMarkerStyle(kFullCircle);
      ce_notstereo_used->SetMarkerColor(kRed);
      ce_notstereo_notused->SetMarkerStyle(kOpenCircle);
      ce_notstereo_notused->SetMarkerColor(kRed);
      bkg_stereo_used->SetMarkerStyle(kFullTriangleUp);
      bkg_stereo_used->SetMarkerColor(kGreen);
      bkg_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      bkg_stereo_notused->SetMarkerColor(kGreen);
      bkg_notstereo_used->SetMarkerStyle(kFullCircle);
      bkg_notstereo_used->SetMarkerColor(kGreen);
      bkg_notstereo_notused->SetMarkerStyle(kOpenCircle);
      bkg_notstereo_notused->SetMarkerColor(kGreen);

      for(unsigned ih=0;ih<xyzp.size();++ih){
	if(xyzp[ih].conversion()){
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      ce_notstereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      ce_notstereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	  }
	}
	else {
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      bkg_notstereo_used->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	    else {
	      bkg_notstereo_notused->Fill(xyzp[ih]._pos.x()-myhel._center.x(),xyzp[ih]._pos.y()-myhel._center.y());
	    }
	  }
	}
      }


      TArc* fitarc = new TArc(0.0,0.0,myhel._radius);
      fitarc->SetLineColor(kRed);
      fitarc->SetLineWidth(2);
      fitarc->SetFillStyle(0);
      // draw the detector boundaries
      static double innerrad(380.0);
      static double outerrad(680.0);
      TArc* indet = new TArc(-myhel._center.x(),-myhel._center.y(),innerrad);
      TArc* outdet = new TArc(-myhel._center.x(),-myhel._center.y(),outerrad);
      indet->SetLineColor(kBlue);
      indet->SetFillStyle(0);
      outdet->SetLineColor(kBlue);
      outdet->SetFillStyle(0);

      TArc* target = new TArc(-myhel._center.x(),-myhel._center.y(),_targetradius);
      target->SetLineColor(kBlack);
      target->SetFillStyle(0);
      // add these to the plot
      TList* flist = ce_stereo_used->GetListOfFunctions();
      flist->Add(fitarc);
      flist->Add(indet);
      flist->Add(outdet);
      flist->Add(target);

      if (mytrk.strawDigiMCCollection() != 0) {
	// Plot the MC true CE hits
	char mctruth_name[100];
	snprintf(mctruth_name,100,"mctshxy%i",igraph);
	TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-500,500,100,-500,500);
	mct->SetMarkerStyle(5);
	mct->SetMarkerColor(kMagenta);
	
	for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	    istr != mytrk.strawHitIndices().end(); ++istr){

	  StrawDigiMC const& mcdigi = mytrk.strawDigiMCCollection()->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if (conversion) {
	    mct->Fill(spmcp->position().x()-myhel._center.x(),spmcp->position().y()-myhel._center.y());
	  }
	}
      }
    }
  }

  void
  HelixDiag::plotZ(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitResult const& myhel) const {
    // Check that there are more than 10 CE hits first
    int n_ce_hits = 0;
    for (XYZPVector::const_iterator i_hit = xyzp.begin(); i_hit != xyzp.end(); ++i_hit) {
      if ((*i_hit).conversion()) {
	++n_ce_hits;
      }
    }

    if (n_ce_hits > 10) {
      static unsigned igraph = 0;
      igraph++;
      art::ServiceHandle<art::TFileService> tfs;

      char ce_stereo_used_name[100];
      snprintf(ce_stereo_used_name,100,"ce_stereo_used_shphiz%i",igraph);
      char ce_stereo_notused_name[100];
      snprintf(ce_stereo_notused_name,100,"ce_stereo_notused_shphiz%i",igraph);
      char ce_notstereo_used_name[100];
      snprintf(ce_notstereo_used_name,100,"ce_notstereo_used_shphiz%i",igraph);
      char ce_notstereo_notused_name[100];
      snprintf(ce_notstereo_notused_name,100,"ce_notstereo_notused_shphiz%i",igraph);
      char bkg_stereo_used_name[100];
      snprintf(bkg_stereo_used_name,100,"bkg_stereo_used_shphiz%i",igraph);
      char bkg_stereo_notused_name[100];
      snprintf(bkg_stereo_notused_name,100,"bkg_stereo_notused_shphiz%i",igraph);
      char bkg_notstereo_used_name[100];
      snprintf(bkg_notstereo_used_name,100,"bkg_notstereo_used_shphiz%i",igraph);
      char bkg_notstereo_notused_name[100];
      snprintf(bkg_notstereo_notused_name,100,"bkg_notstereo_notused_shphiz%i",igraph);
      char title[100];
      snprintf(title,100,"StrawHit #phi Z trk %i;mm;rad",igraph);
      TH2F* ce_stereo_used = tfs->make<TH2F>(ce_stereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* ce_stereo_notused = tfs->make<TH2F>(ce_stereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* ce_notstereo_used = tfs->make<TH2F>(ce_notstereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* ce_notstereo_notused = tfs->make<TH2F>(ce_notstereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_stereo_used = tfs->make<TH2F>(bkg_stereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_stereo_notused = tfs->make<TH2F>(bkg_stereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_notstereo_used = tfs->make<TH2F>(bkg_notstereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
      TH2F* bkg_notstereo_notused = tfs->make<TH2F>(bkg_notstereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);

      ce_stereo_used->SetMarkerStyle(kFullTriangleUp);
      ce_stereo_used->SetMarkerColor(kRed);
      ce_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      ce_stereo_notused->SetMarkerColor(kRed);
      ce_notstereo_used->SetMarkerStyle(kFullCircle);
      ce_notstereo_used->SetMarkerColor(kRed);
      ce_notstereo_notused->SetMarkerStyle(kOpenCircle);
      ce_notstereo_notused->SetMarkerColor(kRed);
      bkg_stereo_used->SetMarkerStyle(kFullTriangleUp);
      bkg_stereo_used->SetMarkerColor(kGreen);
      bkg_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
      bkg_stereo_notused->SetMarkerColor(kGreen);
      bkg_notstereo_used->SetMarkerStyle(kFullCircle);
      bkg_notstereo_used->SetMarkerColor(kGreen);
      bkg_notstereo_notused->SetMarkerStyle(kOpenCircle);
      bkg_notstereo_notused->SetMarkerColor(kGreen);

      for(unsigned ih=0;ih<xyzp.size();++ih){
	if(xyzp[ih].conversion()){
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      ce_notstereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      ce_stereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      ce_notstereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	  }
	}
	else {
	  if (xyzp[ih].use()) {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      bkg_notstereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }	      
	  }
	  else {
	    if (xyzp[ih].stereo()) {
	      bkg_stereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	    else {
	      bkg_notstereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	    }
	  }
	}
      }
      TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
      line->SetParameter(0,myhel._fz0);
      line->SetParameter(1,myhel._dfdz);
      line->SetLineColor(kRed);
      TList* flist = ce_stereo_used->GetListOfFunctions();
      flist->Add(line);

      if (mytrk.strawDigiMCCollection() != 0) {
	// Plot the MC true CE hits
	char mctruth_name[100];
	snprintf(mctruth_name,100,"mctshphiz%i",igraph);
	TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-1500,1500,100,-12.5,12.5);
	mct->SetMarkerStyle(5);
	mct->SetMarkerColor(kMagenta);
	
	for(std::vector<hitIndex>::const_iterator istr=mytrk.strawHitIndices().begin();
	    istr != mytrk.strawHitIndices().end(); ++istr){

	  StrawDigiMC const& mcdigi = mytrk.strawDigiMCCollection()->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if (conversion) {
	    mct->Fill(spmcp->position().z(),spmcp->position().phi());
	  }
	}
      }
    }
  }

    if(_diag > 0){
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	VALERR rad;
	xyzp[ixyzp].rinfo(center,rad);
	_rdiff->Fill(rad._val - rmed);
	_rpull->Fill((rad._val - rmed)/rad._err);
      }
    }

  if(_diag > 0){
	  _fdiff->Fill(xyzp[ixyzp]._phi-phiex);
	  _fpull->Fill((xyzp[ixyzp]._phi-phiex)/fz._phi._err);
	}
 
	// Is this from a conversion hit?
	if(mytrk.strawDigiMCCollection() != 0) {
	  StrawDigiMC const& mcdigi = mytrk.strawDigiMCCollection()->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if (conversion) {
	    pos.setConversion(true);
	  }
	}

}
