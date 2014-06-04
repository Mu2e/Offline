///////////////////////////////////////////////////////////////////////////////
// $Id: TrackCaloMatching_module.cc,v 1.3 2014/06/04 16:46:27 murat Exp $
// $Author: murat $
// $Date: 2014/06/04 16:46:27 $
//
// Original author G. Pezzullo
//
// 2014-06-03 P.Murat: will not work with the vanes-based geometry
///////////////////////////////////////////////////////////////////////////////

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"

// conditions
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"

#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterLink.hh"
#include "TrackCaloMatching/inc/TrackClusterMatch.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
#include "CaloCluster/inc/CaloClusterUtilities.hh"

// Other includes.
#include "cetlib/exception.h"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include "cetlib/pow.h"
#include "BaBar/BaBar/include/Constants.hh"


using namespace std;
using cet::square;
using cet::sum_of_squares;

namespace mu2e {

  class TrackCaloMatching : public art::EDProducer {
  protected:

    std::string     _fitterModuleLabel;
    TrkParticle     _tpart;
    
    TrkFitDirection _fdir;
    std::string     _iname;
					// Diagnostic level
    int             _diagLevel;
    int             _addEnergyToChiSquare;

    double          _emcEnergyThreshold = 10.0;

					// this is a threshold on the time difference between impact time 
					// of the reco-trk and the EMC.This might help the case of mismatching
					// caused by cosmics
    double          _deltaTimeTollerance;
					// Label of the calo clusters  maker
    string          _caloClusterModuleLabel;
    string          _caloClusterAlgorithm;
    string          _caloClusterSeeding;
    string          _caloClusterCollName;
					// Label of the extrapolated impact points

    string          _trkToCaloExtrapolModuleLabel;

    double          _timeChiSquare, _energyChiSquare, _posVChiSquare, _posWChiSquare;

					// no offset in Y ?
    double   _solenoidOffSetX;
    double   _solenoidOffSetZ;

    bool     _skipEvent;
    bool     _firstEvent;

  public:

    explicit TrackCaloMatching(fhicl::ParameterSet const& pset):
      _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
      _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
      _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _addEnergyToChiSquare(pset.get<int>("addEnergyToChiSquare",1)),
      _deltaTimeTollerance(pset.get<double>("deltaTimeTollerance", 50.0)),
      _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
      _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
      _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
      _caloClusterCollName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
      _trkToCaloExtrapolModuleLabel(pset.get<std::string>("trkToCaloExtrapolModuleLabel", "TrkExtrapol")),
      _firstEvent(true)
    {
//-----------------------------------------------------------------------------
// Tell the framework what we make.
//-----------------------------------------------------------------------------
      _iname = _fdir.name() + _tpart.name();
      produces<TrackClusterMatchCollection>();
    }

    virtual ~TrackCaloMatching() {
    }

    void beginJob();
    void endJob  () {}
    void beginRun(art::Run& aRun);
    void produce (art::Event& e);

    void fromTrkToMu2eFrame(HepPoint& vec, CLHEP::Hep3Vector& res);
    
    void doMatching(art::Event & evt, bool skip);
  };


//-----------------------------------------------------------------------------
  void TrackCaloMatching::fromTrkToMu2eFrame(HepPoint&   vec, 
					     CLHEP::Hep3Vector& res) 
  {
    res.setX(vec.x() - _solenoidOffSetX);
    res.setZ(vec.z() - _solenoidOffSetZ);
    res.setY(vec.y());
  }

//-----------------------------------------------------------------------------
  void TrackCaloMatching::beginJob() {
  }



//-----------------------------------------------------------------------------
  void TrackCaloMatching::beginRun(art::Run& aRun) {
    art::ServiceHandle<GeometryService> geom;

					// calculate DS offsets

    _solenoidOffSetX =  geom->config().getDouble("mu2e.solenoidOffset");
    _solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");
  }



//-----------------------------------------------------------------------------
  void TrackCaloMatching::produce(art::Event & evt) {
    doMatching(evt, _skipEvent);
  }


//-----------------------------------------------------------------------------
  void TrackCaloMatching::doMatching(art::Event & evt, bool skip){
    const char* oname = "TrackCaloMatching::doMatching";

    int        nclusters, ntracks, nex, nmatches, vane_id;

    double     cl_v, cl_w, cl_time, cl_energy;
    double     trk_v, trk_w, trk_mom, trk_time;
    double     sigmaV, sigmaW, sigmaT, sigmaE, chiQ;
    double     s1, s2, smean;
    double     nx, ny, dv, dw, dvv, dww, xv, xw, xt, xe;

    double                   chi2_max(1.e12);

    int                      iex, icl, ltrk;
    double                   chi2_best[100][4];
    int                      iex_best [100][4];
    int                      icl_best [100][4];

    const TrkToCaloExtrapol  *extrk;
    const KalRep             *krep;
    const CaloCluster        *cl;

    CLHEP::Hep3Vector        tmpV, cogVaneFrame, tmpPosVaneFrame;
    CLHEP::Hep3Vector        mom, pos; 
    HepPoint                 point;
//-----------------------------------------------------------------------------
// Get handle to calorimeter
//-----------------------------------------------------------------------------
    art::ServiceHandle<GeometryService> geom;
    GeomHandle<Calorimeter> cg;
  
    art::Handle<KalRepCollection> trksHandle;
    evt.getByLabel(_fitterModuleLabel,_iname,trksHandle);
    KalRepCollection const& trks = *trksHandle;
    ntracks = trks.size();
  
    art::Handle<CaloClusterCollection> caloClusters;
    evt.getByLabel(_caloClusterModuleLabel,_caloClusterCollName, caloClusters );
    nclusters = caloClusters->size();

    art::Handle<TrkToCaloExtrapolCollection>  trjExtrapols;
    evt.getByLabel(_trkToCaloExtrapolModuleLabel, trjExtrapols);
    nex = trjExtrapols->size();

    if (_diagLevel > 2) {
      printf(" %s: Event Number : %8i ntracks: %2i nclusters: %4i nex: %2i\n",
	     oname,evt.event(), ntracks, nclusters, nex);
    }

    unique_ptr<TrackClusterMatchCollection> tcmcoll(new TrackClusterMatchCollection);
    tcmcoll->reserve(nex);
//-----------------------------------------------------------------------------
// here the execution really starts
//-----------------------------------------------------------------------------
    if (ntracks == 0)                                         goto END;

    for (int it=0; it<ntracks; it++) {
      for (int iv=0; iv<4; iv++) {
	chi2_best[it][iv] = chi2_max+1;
	iex_best [it][iv] = -1;
	icl_best [it][iv] = -1;
      }
    }
//-----------------------------------------------------------------------------
// loop over the extrapolated tracks (TrkToCaloExtrapol's), 
// each of them corresponds to a track-to-disk intersection, so one track may 
// have two corresponding TrkToCaloExtrapol things
//-----------------------------------------------------------------------------
    for (int jex=0; jex<nex; jex++) {

      if(_diagLevel > 2) {
	printf(" %s : jex = %5i\n", oname, jex);
      }
    
      extrk = &trjExtrapols->at(jex);
      krep  = *extrk->trk();
//-----------------------------------------------------------------------------
// track index, important: store one, the best, intersection per track per vane
//-----------------------------------------------------------------------------
      ltrk  = extrk->trackNumber();
      if (ltrk > 100) {
	printf(">>> ERROR in %s: more than 100 tracks, ltrk = %i, skip the rest\n",
	       oname,ltrk);
                                                            goto NEXT_INTERSECTION;
      }
      else if (ltrk < 0) {
	printf(">>> ERROR in %s: ltrk < 0:, ltrk = %i, skip the rest\n", 
	       oname,ltrk);
                                                            goto NEXT_INTERSECTION;
      }
//-----------------------------------------------------------------------------
// apparently, ntupling stuff was mixed in, almost removed
//-----------------------------------------------------------------------------
      vane_id    = extrk->vaneId();

      trk_time   = extrk->time();
					// assume Z(middle of the disk)

      s1      = extrk->pathLengthEntrance();
      s2      = extrk->pathLengthExit    ();
      smean   = (s1+s2)/2;
					// 

      mom     = krep->momentum(smean);

      smean   = smean-60.*mom.mag()/mom.z();
      //      smean   = s1; // checking...

      point   = krep->position(smean);
      mom     = krep->momentum(smean);
      trk_mom = mom.mag();

      if(_diagLevel > 2){
	cout<<"vane_id = "<<vane_id<<endl;
	cout<<"trk_mom = "<<trk_mom<<" [MeV]" << endl;
      }

      nx      = mom.x()/sqrt(mom.x()*mom.x()+mom.y()*mom.y());
      ny      = mom.y()/sqrt(mom.x()*mom.x()+mom.y()*mom.y());

					// transform position to Mu2e detector frame

      fromTrkToMu2eFrame(point, tmpPosVaneFrame);
	  
      if (_diagLevel > 2){
	cout<<"Mu2e general frame: tmpPosVaneFrame = "<< tmpPosVaneFrame <<  " [mm]" << endl;
      }
    
      tmpV    = cg->toSectionFrame(vane_id,tmpPosVaneFrame);
      trk_v   = tmpV.x();
      trk_w   = tmpV.y();
	  
      if(_diagLevel > 2){
	cout << "tmpPosVaneFrame = "<< tmpV << " [mm]" << endl;
      }
//-----------------------------------------------------------------------------
// loop over clusters
//-----------------------------------------------------------------------------
      for (int icl=0; icl<nclusters; icl++) {
	cl  = &(*caloClusters).at(icl);
	CaloClusterTools cluTool(*cl);

	if (cl->vaneId()    != vane_id           )                                    goto NEXT_CLUSTER;
	if (cl->energyDep() < _emcEnergyThreshold)                                    goto NEXT_CLUSTER;
	if (std::fabs(cluTool.timeFasterCrystal() - trk_time) > _deltaTimeTollerance) goto NEXT_CLUSTER;
//-----------------------------------------------------------------------------
// 2013-03-24 P.Murat: cluster center of gravity is determined in the GLOBAL 
//                     coordinate system, 
// clustering needs to be fixed to define cluster coordinates in the local frame system
// for now - transform back to the local coordinate system
//-----------------------------------------------------------------------------
	cogVaneFrame = cg->toSectionFrame(vane_id, cl->cog3Vector());
      
	cl_time      = cl->time();
	cl_v         = cogVaneFrame.x();
	cl_w         = cogVaneFrame.y();
	cl_energy    = cl->energyDep();
//-----------------------------------------------------------------------------
// V and W - coordinates in the disk system
// rotate them in the direction perpendicular to the track
//-----------------------------------------------------------------------------
	dv  = trk_v-cl_v;
	dw  = trk_w-cl_w;

	dvv = dv*nx+dw*ny;
	dww = dv*ny-dw*nx;
//-----------------------------------------------------------------------------
// ad-hoc corrections 
// at this point do not understand their origin
// the numbers, obviously, come from analysis distributions and correspond 
// to the BaF2 crystal length of 21cm
//-----------------------------------------------------------------------------
// 	dvv -= 70;
// 	dww += 18;
//-----------------------------------------------------------------------------
// 2014-05-14 P.Murat: use 10 MeV as the matching resolution, 
//            for a 10 MeV cluster the energy term in the chi2 would be (90/10)^2
//            also set coordinate resolution to 5cm , need to try using dv only
//-----------------------------------------------------------------------------
	sigmaE = 10.;			// sigma(E) = 10 MeV
	sigmaV = 40.;			// 40 mm
	sigmaW = 10.; 			// 10 mm
	sigmaT = 0.5; 			// 0.5 ns
	  
	if (_diagLevel > 2){
	  cout << " trk_v = "   << trk_v 
	       << ", cl_v = "   << cl_v 
	       << ", sigmaV = " << sigmaV << endl 
	       << ", trk_w = "  << trk_w 
	       << ", cl_w = "   << cl_w 
	       << ", sigmaW = " << sigmaW << endl 
	       << ", cl_time = "<< cl_time 
	       << ", sigmaT = " << sigmaT 
	       << endl;
	}
					// need to handle energy part properly, later! 
	xv = dvv/sigmaV;
	xw = dww/sigmaW;
	xt = (trk_time-cl_time)/sigmaT;
	xe = (trk_mom-cl_energy)/sigmaE;

	_posVChiSquare   = xv*xv;
	_posWChiSquare   = xw*xw;
	_timeChiSquare   = xt*xt;
	_energyChiSquare = xe*xe;

	chiQ = _posVChiSquare + _posWChiSquare + _timeChiSquare;
    
	if (_addEnergyToChiSquare == 1) {
	  chiQ += _energyChiSquare;
	}
    
	if (chiQ < chi2_best[ltrk][vane_id]) {
//-----------------------------------------------------------------------------
// new best match
//-----------------------------------------------------------------------------
	  chi2_best[ltrk][vane_id] = chiQ;
	  icl_best [ltrk][vane_id] = icl;
	  iex_best [ltrk][vane_id] = jex;
	}
      NEXT_CLUSTER:;
      }
    NEXT_INTERSECTION:;
    }

//-----------------------------------------------------------------------------
// form output list of matches
//-----------------------------------------------------------------------------
    nmatches = 0;
    //    TrackClusterMatch tcm;

    for (int it=0; it<ntracks; it++) {
      for (int iv=0; iv<4; iv++) {
	if (chi2_best[it][iv] < chi2_max) {
	  iex = iex_best [it][iv];
	  icl = icl_best [it][iv];

	  art::Ptr<TrkToCaloExtrapol> artPtrTex    (trjExtrapols,iex);
	  art::Ptr<CaloCluster>       artPtrCluster(caloClusters,icl);
	  TrackClusterMatch           tcm(artPtrTex,artPtrCluster,chi2_best[it][iv]);
	  tcmcoll->push_back(tcm);
	  nmatches += 1;
	}
      }
    }
    
    END:;

    evt.put(std::move(tcmcoll));

    if (evt.id().event() %100 == 0) {
      printf("Event %d TrackCaloMatching done...",evt.id().event() );
    }
  }

}

using mu2e::TrackCaloMatching;
DEFINE_ART_MODULE(TrackCaloMatching);
