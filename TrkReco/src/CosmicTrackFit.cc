// Calls fitting objects to perform track fit to ComboHits
// $Id: CosmicTrackFit code
// $Author: S Middleton
// $Date: Feb 2019
//
// Mu2e Cosmics:
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"
// art
#include "canvas/Persistency/Common/Ptr.h"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
//For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
//Fitting
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/BuildMatrixSums.hh"
#include "Mu2eUtilities/inc/LiklihoodFunctions.hh"
#include "Mu2eUtilities/inc/DriftFitUtils.hh"
//ROOT:
#include "TMatrixD.h"
#include "Math/VectorUtil.h"
#include "TH1F.h"
#include "TPolyMarker.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/DistFunc.h"
//Minit2:
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>


using namespace std;
using namespace boost::accumulators;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using namespace ROOT::Math;

std::vector<double> TrackerConstraints(1500, 1500);

struct ycomp : public std::binary_function<XYZVec,XYZVec,bool> {
    bool operator()(XYZVec const& p1, XYZVec p2) { return p1.y() > p2.y(); }
  };
  
namespace mu2e
{
  
  CosmicTrackFit::CosmicTrackFit(fhicl::ParameterSet const& pset) :
    _Npara(pset.get<unsigned>("Npara",4)),
    _diag(pset.get<int>("diagLevel",1)),
    _mcdiag(pset.get<int>("MCdiagLevel",1)),
    _debug(pset.get<int>("debugLevel",1)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnsh(pset.get<unsigned>("minNStrawHits",2)),
    _minCHHits(pset.get<unsigned>("minCHHits",4)),
    _n_outliers(pset.get<unsigned>("_n_outliers",5)),
    _maxniter(pset.get<unsigned>("maxniter",1000)),//10
    _maxpull(pset.get<float>("maxPull",3)),
    _maxd(pset.get<float>("maxd",300.0)),//max distance between hits at start of fit
    _maxDOCA(pset.get<float>("maxDOCA",1000)),//max distance of closest approach between a hit included in fit and one flag out-right as an outlier
    _maxchi2(pset.get<float>("maxchi2",10.0)) ,
    _max_chi2_change(pset.get<float>("max_chi2_change",0.001)),
    _max_position_deviation((pset.get<float>("max_position_deviation",200)))
    {}

    CosmicTrackFit::~CosmicTrackFit(){}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics) {
    
    bool is_ok(false);
    RunFitChi2(title, TrackData,diagnostics );
    is_ok = TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackOK);
    return is_ok;
  }

  /*-------------Line Direction-------------------------//
    Range of methods for finding track directions - some are redundent-check!
  //----------------------------------------------*/
  //ComboHits
  XYZVec CosmicTrackFit::InitLineDirection(const ComboHit *ch0, const ComboHit *chN) {
      float tx = chN->pos().x() - ch0->pos().x();
      float ty = chN->pos().y() - ch0->pos().y();
      float tz = chN->pos().z() - ch0->pos().z();       
      XYZVec track(tx,ty,tz);
    
      return track.Unit();
    } 
  //MCDigis: 
  XYZVec CosmicTrackFit::InitLineDirection( StrawDigiMC const& mc0,  StrawDigiMC const& mcN, XYZVec reco_dir, bool is_prime) {
       art::Ptr<StepPointMC> const& spmcp0 = mc0.stepPointMC(StrawEnd::cal);
       XYZVec pos0(spmcp0->position().x(), spmcp0->position().y(), spmcp0->position().z());//det->toDetector(spmcp0->position());
       art::Ptr<StepPointMC> const& spmcpN = mcN.stepPointMC(StrawEnd::cal);
       XYZVec posN(spmcpN->position().x(), spmcpN->position().y(), spmcpN->position().z());
       if(is_prime == true){
          std::vector<XYZVec> AxesList = ParametricFit::GetAxes(reco_dir); 
       	  posN.SetXYZ(posN.Dot(AxesList[0]), posN.Dot(AxesList[1]), posN.Dot(AxesList[2]));
       	  pos0.SetXYZ(pos0.Dot(AxesList[0]), pos0.Dot(AxesList[1]), pos0.Dot(AxesList[2]));
       }
       float tx = posN.x() -  pos0.x();
       float ty = posN.y() -  pos0.y();
       float tz = posN.z() -  pos0.z();
       XYZVec track(tx,ty,tz);
       return track.Unit();
    } 
  
    XYZVec CosmicTrackFit::LineDirection(double a1, double b1, const ComboHit *ch0, const ComboHit *chN, XYZVec ZPrime) {
      XYZVec track(a1,b1,1);
      return track.Unit();
    } 
  
//--------------Fit-----------------//
// This is the top level call to Fitting routines....
//-------------------------------------------// 
  void CosmicTrackFit::BeginFit(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics){ 
    //Clear Previous Flags:
    TrackData._tseed._status.clear(TrkFitFlag::StraightTrackOK); 
    // Initialize:
    bool init(false);
    if (!TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackInit)) {
      init = true;
      if (initCosmicTrack( title, TrackData, diagnostics))
	TrackData._tseed._status.merge(TrkFitFlag::StraightTrackInit);
      else
	return;
    } 
    //Start Chi2 Fitting:
    if (!init)RunFitChi2(title, TrackData, diagnostics);
  }
/*------------------------------Chi 2 Fit----------------------------//
//   Adds Chi-2 optimization to fitting routine //
//   Refits and adjusts track fit paramters by weights//
//------------------------------------------------------------------*/
void CosmicTrackFit::RunFitChi2(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics) {   
   CosmicTrack* track = &TrackData._tseed._track; 
   TrackData._tseed._status.merge(TrkFitFlag::StraightTrackOK);  
   TrackData._tseed._status.merge(TrkFitFlag::StraightTrackConverged);
   FitAll(title, TrackData, track, diagnostics); 
   
   TrackData._diag.CosmicTrackFitCounter += 1;//label as having a track
   diagnostics.nChFit = (TrackData._chHitsToProcess.size());
 
}


  
  std::vector<XYZVec> SortPoints(std::vector<XYZVec> pointY){
        std::vector<XYZVec> sortedPoints;
  	std::sort(pointY.begin(), pointY.end(),ycomp());
  	for (unsigned i=0; i<pointY.size(); ++i) { 
      		sortedPoints.push_back(pointY[i]);
      	}
      	return sortedPoints;
   }
   

 void CosmicTrackFit::FitAll(const char* title, CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, CosmicTrackFinderTypes::Data_t& diagnostics){
 
     ::BuildMatrixSums S;
     ComboHit   *hitP1(0), *hitP2(0); 
     size_t nHits (trackData._chHitsToProcess.size());   
     int DOF = (nHits);// - (_Npara);
     const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
     const ComboHit* chN = &trackData._chHitsToProcess[trackData._chHitsToProcess.size()-1]; 
    
     //Step 1: Get Initial Estimate of track direction
     XYZVec ZPrime = InitLineDirection(ch0, chN);  
     std::vector<XYZVec> AxesList = ParametricFit::GetAxes(ZPrime);
     TrackAxes InitAxes = ParametricFit::GetTrackAxes(ZPrime);
    
    //Step 2: Loop over hits and get track parameters based on above estimated track direction
    for (size_t f1=0; f1<nHits; ++f1){  
      hitP1 = &trackData._chHitsToProcess[f1];  
      if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;  
      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, InitAxes._XDoublePrime, InitAxes._YDoublePrime); 
      S.addPoint(point, InitAxes._XDoublePrime, InitAxes._YDoublePrime, InitAxes._ZPrime, ErrorsXY[0], ErrorsXY[1]); 
     }    
 
      //Step 3: Get the first estmiate of track parameters, get a updated track direction vector from these parameters
      double a0 = S.GetAlphaX()[0][0];
      double a1 = S.GetAlphaX()[1][0];
      double b0 = S.GetAlphaY()[0][0];
      double b1 = S.GetAlphaY()[1][0];
      XYZVec Direction(a1,b1,1);
      XYZVec UpdatedTrackDirection =Direction.Unit();
    
      //Step 4: Update axes and store them as initial axes for plotting
      ZPrime = UpdatedTrackDirection;
      TrackAxes Axes = ParametricFit::GetTrackAxes(ZPrime);
      TrackParams InitParams(a0,a1,b0,b1);
      
      cosmictrack->SetInitParams(InitParams);
      cosmictrack->SetTrackDirection(UpdatedTrackDirection);
      cosmictrack->SetInitTrackCoordSystem(Axes);
      cosmictrack->SetFitTrackCoOrdSystem(Axes);
      //Step 5: Loop for initial diagnostics
      if(_debug>0){
	cosmictrack->set_initchisq_dofY(S.GetChi2Y()/abs(DOF));
	cosmictrack->set_initchisq_dofX(S.GetChi2X()/abs(DOF));
	cosmictrack->set_initchisq_dof(S.GetTotalChi2()/abs(DOF));
	for (size_t f2=0; f2<nHits; ++f2){
         	if(isnan(cosmictrack->GetTrackDirection().Mag2()) == true) continue;     
         	hitP1 = &trackData._chHitsToProcess[f2];  
		if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;  
	        XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
		XYZVec point_prime(point.Dot(Axes._XDoublePrime), point.Dot(Axes._YDoublePrime), point.Dot(ZPrime));
		
		std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, Axes._XDoublePrime,Axes._YDoublePrime);   
		float Rx = ParametricFit::GetResidualX(a0,a1, point_prime);
		float Ry = ParametricFit::GetResidualY(b0,b1, point_prime);
		
		cosmictrack->set_init_fit_residualsX(Rx);
		cosmictrack->set_init_fit_residualsY(Ry); 
		cosmictrack->SetInitErrors(ErrorsXY[0],ErrorsXY[1]); 
		     
		} 
      }
      
      //Step 6: Begin iteration for finding the best track fit possible.
      unsigned niter(0);
      ::BuildMatrixSums S_niteration;
      bool converged = false;
      CosmicTrack* BestTrack = cosmictrack;	 
      double chi2_best_track = 10000000;//chosen arbitary high number
      //Step 7: find best prev chi2
      while(niter < _maxniter && converged==false){                 
      	niter +=1; 
 	double previous_chi2, changed_chi2;
      	if (niter == 1) {
 	     	 previous_chi2 = S.GetTotalChi2()/abs(DOF);
 	     	 chi2_best_track = previous_chi2;
       }
      	else {
      	 	previous_chi2 = S_niteration.GetTotalChi2()/abs(DOF);
      	}
     	//Step 8: Remove previous sums and use previously stored axes from updated last iteration:
     	S_niteration.clear(); 
     	Axes = ParametricFit::GetTrackAxes(cosmictrack->GetTrackDirection());
     	cosmictrack->set_niter(niter );
     	for (size_t f4=0; f4 < nHits; ++f4){ 
     	      if(isnan(cosmictrack->GetTrackDirection().Mag2()) == true) continue; 
     	      hitP2 = &trackData._chHitsToProcess[f4];
      	      if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;   
	      XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());	    	  
	      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP2, Axes._XDoublePrime, Axes._YDoublePrime);	 
	      S_niteration.addPoint(point, Axes._XDoublePrime, Axes._YDoublePrime, Axes._ZPrime, ErrorsXY[0],ErrorsXY[1]);
    	  }
          //Get new parameters
     	  a0 = S_niteration.GetAlphaX()[0][0];
          a1 = S_niteration.GetAlphaX()[1][0];
	  b0 = S_niteration.GetAlphaY()[0][0];
	  b1 = S_niteration.GetAlphaY()[1][0];
	 
    	  //Step 9: Get new direction:
	   XYZVec DirectionSecond(a1,b1,1);
	   XYZVec UpdatedTrackDirectionSecond = DirectionSecond.Unit();
	   
 	   //Step 10 - Get New Axes:
 	   //Axes = ParametricFit::GetTrackAxes(UpdatedTrackDirectionSecond);
	   cosmictrack->SetTrackDirection(UpdatedTrackDirectionSecond); 
	   cosmictrack->SetFitTrackCoOrdSystem(Axes);
	   //Step 11: Update Parameters:
	   TrackParams FitParams(a0,a1,b0,b1);
	   cosmictrack->SetFitParams(FitParams);
	
           //Step 12: Update Chi2:
           if(_debug > 0){ 
           	   cosmictrack->Diag.FinalChiTot = (S_niteration.GetTotalChi2()/abs(DOF)); 
		   cosmictrack->Diag.FinalChiY = S_niteration.GetChi2Y()/abs(DOF);
		   cosmictrack->Diag.FinalChiX = S_niteration.GetChi2X()/abs(DOF);
	   }
	
	   if(abs((a0) - (InitParams.A0)) > _max_position_deviation or abs((b0) - (InitParams.B0)) > _max_position_deviation){continue;}
	   
	   //Step13: Check if Chi2 is improved:, If so call this "BestTrack"
	   float updated_chi2 = S_niteration.GetTotalChi2()/abs(DOF);
	   changed_chi2 = chi2_best_track - updated_chi2;
	  
           if (updated_chi2 < chi2_best_track  ){
		
		BestTrack->clear_diag();
		TrackParams FitParams(a0,a1,b0,b1);
	   	BestTrack->SetFitParams(FitParams);
	        BestTrack->SetFitTrackCoOrdSystem(Axes);
	        BestTrack->SetTrackDirection(cosmictrack->GetTrackDirection());
	        
      		XYZVec EndTrackPosition(BestTrack->FitParams.A0,BestTrack->FitParams.B0,0);
      		TrackEquation EndTrack(EndTrackPosition, BestTrack->TrackCoordSystem._ZPrime);
      		BestTrack->SetTrackEquation(EndTrack);
      		
   		//Step 14: Save this Chi2:
                if(_debug > 0){
                	BestTrack->Diag.FinalChiTot = (S_niteration.GetTotalChi2()/abs(DOF)); 
		   	BestTrack->Diag.FinalChiY = S_niteration.GetChi2Y()/abs(DOF);
		        BestTrack->Diag.FinalChiX = S_niteration.GetChi2X()/abs(DOF);
                	
	        }
                chi2_best_track = updated_chi2;  
                cosmictrack=BestTrack;      	
	      }

         //If on any iteration the change in chi2 becomes small than we can assume we have converged on a track
         if( abs(changed_chi2) < _max_chi2_change ){
                  converged = true;
                  cosmictrack->converged = true;     
 	      }
 	      
 //If at end of the iteration process but the track is still not converging then remove it
               if(niter == _maxniter && converged ==false ){
 		    	 trackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
 		    	 trackData._tseed._status.clear(TrkFitFlag::StraightTrackConverged);
 		    	 continue;
		 }
		 
   }//end while 
   cosmictrack=BestTrack;
   if(cosmictrack->converged and  _diag>0){
	 for (size_t f5=0; f5<nHits; ++f5){
     		hitP2 = &trackData._chHitsToProcess[f5];
                if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;
		XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
		XYZVec point_prime(point.Dot(BestTrack->TrackCoordSystem._XDoublePrime), point.Dot(BestTrack->TrackCoordSystem._YDoublePrime), point.Dot(BestTrack->TrackCoordSystem._ZPrime));
		
		std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, BestTrack->TrackCoordSystem._XDoublePrime, BestTrack->TrackCoordSystem._YDoublePrime);  
		BestTrack->SetFinalErrors(ErrorsXY[0], ErrorsXY[1]);
	      	float newRx = ParametricFit::GetResidualX(BestTrack->FitParams.A0, BestTrack->FitParams.A1, point_prime);
		float newRy = ParametricFit::GetResidualY(BestTrack->FitParams.B0, BestTrack->FitParams.B1, point_prime);
		BestTrack->SetCovarience(S_niteration.GetCovX()[0][0],S_niteration.GetCovX()[1][1], S_niteration.GetCovY()[0][0], S_niteration.GetCovY()[1][1]);
		BestTrack->Diag.FinalResidualsX.push_back(newRx);
	      	BestTrack->Diag.FinalResidualsY.push_back(newRy);
	      	cosmictrack = BestTrack;      
	      	cout<<" size of diags  "<<BestTrack->Diag.FinalResidualsX.size()<<endl;		
      	}	                	
    }

     S.clear();
     if(_mcdiag > 0 and cosmictrack->converged == true){
          FitMC(trackData,cosmictrack, true, false);
     } 
     //cosmictrack->set_fit_phi(acos(cosmictrack->getZPrime().x()/sqrt(cosmictrack->getZPrime().Mag2())));
   }

//Applies fitting routine to MCDigis:
void CosmicTrackFit::FitMC(CosmicTrackFinderData& trackData, CosmicTrack* cosmictrack, bool XYZ, bool prime){	
        ::BuildMatrixSums S; 
    	XYZVec FitDirection = cosmictrack->GetTrackDirection();
    	size_t nHits (trackData._mcDigisToProcess.size());
        StrawDigiMC const& ch0 = trackData._mccol->at(0); 
        StrawDigiMC const& chN = trackData._mccol->at(nHits-1);
        XYZVec ZPrime = InitLineDirection(ch0, chN, FitDirection, prime);  
        StrawDigiMC *hitP1; 
        for (size_t f1=0; f1<nHits; ++f1){ 
            hitP1 = &trackData._mcDigisToProcess[f1]; 
            art::Ptr<StepPointMC> const& spmcp = hitP1->stepPointMC(StrawEnd::cal);
            XYZVec posN(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
            if(XYZ == false && prime ==false){ //gets the true X"Y"Z' and does in that frame
                XYZVec point(posN.x(), posN.y(), posN.z());
                TrackAxes FitAxes = ParametricFit::GetTrackAxes(ZPrime);
                S.addPoint(point, FitAxes._XDoublePrime,  FitAxes._YDoublePrime,  FitAxes._ZPrime, 1,1);
            }
            if(XYZ == true && prime ==false){ //gets the true XYZ and does in that frame USE THIS THEN TRANSFORM
                XYZVec point(posN.x(), posN.y(), posN.z());
                XYZVec X(1,0,0);//to evaluate for XYZ frame
                XYZVec Y(0,1,0);
                XYZVec Z(0,0,1);
                S.addPoint( point, X,Y,Z, 1,1);
            }
            if(XYZ == false && prime ==true){//true for calc in X"Y"Z' frame where these are from reco
                TrackAxes FitAxes = ParametricFit::GetTrackAxes(FitDirection);
                XYZVec point(posN.Dot(FitAxes._XDoublePrime), posN.Dot( FitAxes._YDoublePrime), posN.Dot( FitAxes._ZPrime));
                S.addPoint(point,  FitAxes._XDoublePrime,  FitAxes._YDoublePrime,  FitAxes._ZPrime, 1,1); 
            }
          
        }   
     //Transform into XYZ fitplane
     XYZVec Direction(S.GetAlphaX()[1][0],S.GetAlphaY()[1][0],1);
     XYZVec Position(S.GetAlphaX()[0][0],S.GetAlphaY()[0][0],0);

     XYZVec transformedDirection(Direction.Dot(cosmictrack->TrackCoordSystem._XDoublePrime)/Direction.Dot(cosmictrack->TrackCoordSystem._ZPrime), Direction.Dot(cosmictrack->TrackCoordSystem._YDoublePrime)/Direction.Dot(cosmictrack->TrackCoordSystem._ZPrime), 1);
     XYZVec unit_transformedDirection = transformedDirection.Unit();
     XYZVec transformedPosition(Position.Dot(cosmictrack->TrackCoordSystem._XDoublePrime), Position.Dot(cosmictrack->TrackCoordSystem._YDoublePrime), Position.Dot(cosmictrack->TrackCoordSystem._ZPrime));
     //Fill details in cosmictrack:
     double mcpos_mag = sqrt(Direction.Mag2());    
     TrackAxes TrueAxes = ParametricFit::GetTrackAxes(ZPrime);
     cosmictrack->SetTrueTrackCoordSystem(TrueAxes);
     TrackParams TrueParams(transformedPosition.X(),unit_transformedDirection.X(),transformedPosition.Y(),unit_transformedDirection.Y());
     cosmictrack->SetTrueParams(TrueParams);
     
     const double theta_start = acos(Direction.x()/mcpos_mag);
     const double phi_start = atan(Direction.y()/Direction.x());
     cosmictrack->set_true_finalchisq_dof(S.GetTotalChi2()/nHits);
     cosmictrack->set_true_phi(phi_start);
     cosmictrack->set_true_theta(theta_start);
     
     }


XYZVec CosmicTrackFit::MCInitHit(StrawDigiMC mcdigi){
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(StrawEnd::cal);
	XYZVec first(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
        return first;
}

XYZVec CosmicTrackFit::MCFinalHit(StrawDigiMC mcdigi){
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(StrawEnd::cal);
	XYZVec last(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
        return last;
}


void CosmicTrackFit::MCDirection(XYZVec first, XYZVec last, CosmicTrackFinderData& trackData){ 
      double tx = last.x() - first.x();
      double ty = last.y() - first.y();
      double tz = last.z() - first.z();       
      XYZVec track(tx,ty,tz);
      //trackData._tseed._track.set_true_track_direction(track.Unit());
      TrackAxes TrueAxes = ParametricFit::GetTrackAxes(track.Unit());
      trackData._tseed._track.SetTrueTrackCoordSystem(TrueAxes);
}
//Some Functions to store some analysis functions (Currently not doing a very good job)
  float CosmicTrackFit::PDF(float chisq, float ndf){
  	float pdf = ROOT::Math::chisquared_pdf(chisq, ndf);
  	return pdf;
  }
  
  float CosmicTrackFit::chi_sum(float chisq, float ndf){
  	float sum = TMath::ChisquareQuantile(chisq, ndf);
  	return sum;
  }
  	
  float CosmicTrackFit::CDF(float chisq, float ndf){	
  	float cdf;
  	if(chisq/ndf < 3){	
  	cdf = ROOT::Math::chisquared_cdf_c(chisq, ndf , 0);	
  	}
  	else{
  	cdf = ROOT::Math::chisquared_cdf(chisq, ndf , 0);
  	}
  	return cdf;
  }
bool CosmicTrackFit::goodTrack(CosmicTrack* track)
  { 
    if(track->Diag.FinalChiTot < _maxchi2) return true;
    else return false;
  }

/*--------------USE? ---------------------------//
//          Checks if flag                        //
//------------------------------------------------*/
bool CosmicTrackFit::use_hit(const ComboHit& thit) const 
  {
    return (!thit._flag.hasAnyProperty(_dontuseflag));
  }
bool CosmicTrackFit::use_track(double track_length) const 
  {
     return (track_length > _maxd) ? false : true ;
  }
  
void CosmicTrackFit::DriftFit(CosmicTrackFinderData& trackData){
	 cout<<"======== Drift Fit Diagnostics ============="<<endl;
         EndResult endresult = LiklihoodFunctions::DoFit(trackData._tseed,  _srep);
         //Store fits:
         trackData._tseed._track.MinuitFitParams.A0 =  endresult.bestfit[0];//a0
         trackData._tseed._track.MinuitFitParams.A1 =  endresult.bestfit[1];//a1
         trackData._tseed._track.MinuitFitParams.B0 =  endresult.bestfit[2];//b0
         trackData._tseed._track.MinuitFitParams.B1 =  endresult.bestfit[3];//b1
         trackData._tseed._track.MinuitFitParams.T0 =  endresult.bestfit[4];//t0
         //Get New Axesa:
	 XYZVec NewDirection(endresult.bestfit[1], endresult.bestfit[3], 1);
	 NewDirection = NewDirection.Unit();
	 TrackAxes Axes = ParametricFit::GetTrackAxes(NewDirection);
	 //Fill Diagnostics:
	 trackData._tseed._track.MinuitCoordSystem = Axes; 
	 trackData._tseed._track.DriftDiag.EndDOCAs = endresult.EndDOCAs;
	 trackData._tseed._track.DriftDiag.StartDOCAs = endresult.StartDOCAs;
	 trackData._tseed._track.DriftDiag.EndTimeResiduals = endresult.EndTimeResiduals;
	 trackData._tseed._track.DriftDiag.StartTimeResiduals = endresult.StartTimeResiduals;
	 if(_diag > 0){
		ComboHit *chit(0);
		for(size_t j=0; j<(trackData._chHitsToProcess.size()); j++){
		    //Get Hit:
		    chit = &trackData._chHitsToProcess[j];
		    //if(trackData._tseed._track.minuit_converged == false) continue;
		    if (((!use_hit(*chit) ) && (chit->nStrawHits() < _minnsh) )) continue;
		    //Get point:
		    XYZVec point(chit->pos().x(),chit->pos().y(),chit->pos().z());
		    //Convert Point
		    XYZVec chit_prime(point.Dot(trackData._tseed._track.MinuitCoordSystem._XDoublePrime), point.Dot(trackData._tseed._track.MinuitCoordSystem._YDoublePrime), point.Dot(trackData._tseed._track.MinuitCoordSystem._ZPrime));
		    
		    //Set Residuals:
		    trackData._tseed._track.DriftDiag.FinalResidualsX.push_back(ParametricFit::GetResidualX(trackData._tseed._track.MinuitFitParams.A0,  trackData._tseed._track.MinuitFitParams.A1, chit_prime ));
		    trackData._tseed._track.DriftDiag.FinalResidualsY.push_back(ParametricFit::GetResidualY( trackData._tseed._track.MinuitFitParams.B0,  trackData._tseed._track.MinuitFitParams.B1, chit_prime ));
		    
		   //Set Errors:
		    std::vector<double> ErrorsXY = ParametricFit::GetErrors(chit, trackData._tseed._track.MinuitCoordSystem._XDoublePrime, trackData._tseed._track.MinuitCoordSystem._YDoublePrime);
		     	
		    trackData._tseed._track.DriftDiag.FinalErrX.push_back(ErrorsXY[0]);
		    trackData._tseed._track.DriftDiag.FinalErrY.push_back(ErrorsXY[1]);
		}
      }
      /*
       //If minuit does not converge then remove "OK" flags: HERE LOSING ALL HITS
      if(trackData._tseed._track.minuit_converged ==false ){
	  trackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
	  trackData._tseed._status.clear(TrkFitFlag::StraightTrackConverged);
	  
      }
      */
}
/*
 TrkStrawHits  filler (not working)
 std::vector<TrkStrawHitSeed>const trkseeds = trackData._tseed.trkstrawhits();
     cout<<"trkseed size "<<trkseeds.size()<<endl;
     for(auto const& ths : trkseeds ){
      // create a TrkStrawHit from this seed.
      size_t index = ths.index();
      cout<<"index "<<index<<endl;
      const ComboHit& strawhit(trackData.chcol()->at(index));
      const Straw& straw = _tracker->getStraw(strawhit.strawId());
      StrawResponse::cptr_t strawResponse;
     
      TrkStrawHit* trkhit = new TrkStrawHit(strawResponse,strawhit,straw,ths.index(),ths.t0(),100., 5.,1.);
      cout<<" Phi "<<trkhit->driftPhi()<<" v drift "<<trkhit->driftVelocity()<<" time"<<trkhit->driftTime()<<endl;
        }
*/
}//end namespace
