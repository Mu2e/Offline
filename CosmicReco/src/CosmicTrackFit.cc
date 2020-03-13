// Author: S Middleton
// Date: March 2019
// Purpose: Holds functions for the fitting of Cosmic Tracks in tracker

// Mu2e Cosmics:
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"

// art
#include "canvas/Persistency/Common/Ptr.h"

//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

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
#include "Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"
#include "CosmicReco/inc/MinuitDriftFitter.hh"
#include "CosmicReco/inc/DriftFitUtils.hh"
#include "DataProducts/inc/XYZVec.hh"

//ROOT:
#include "TMatrixD.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/DistFunc.h"

// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

//String:
#include <string>

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
    CosmicTrackFit::CosmicTrackFit(const Config& conf) :
	_Npara (conf.Npara()),
	_diag (conf.diag()),
	_debug  (conf.debug()),
	_dontuseflag (conf.dontuseflag()),
	_minnsh   (conf.minnsh()),
	_minnch  (conf.minnch()),
	_n_outliers (conf.n_outliers()),
	_maxniter (conf.maxniter()),
	_max_seed_chi2 (conf.max_seed_chi2()),
	_max_chi2_change (conf.max_chi2_change()),
	_max_position_deviation (conf.max_position_deviation()),
	_maxHitDOCA (conf.maxHitDOCA()),
	_maxLogL (conf.maxLogL()),
	_gaussTres (conf.gaussTres()),
	_maxTres (conf.maxTres()),
	_maxd (conf.maxd()),
	_maxpull (conf.maxpull())
	{}
	 
    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(const char* title, CosmicTrackFinderData& TrackData) {
    
	bool is_ok(false);
	RunFitChi2(title, TrackData);
	is_ok = TrackData._tseed._status.hasAllProperties(TrkFitFlag::helixOK);
	return is_ok;
  }

  std::vector<XYZVec> SortPoints(std::vector<XYZVec> pointY){
        std::vector<XYZVec> sortedPoints;
  	std::sort(pointY.begin(), pointY.end(),ycomp());
  	for (unsigned i=0; i<pointY.size(); ++i) { 
      		sortedPoints.push_back(pointY[i]);
      	}
      	return sortedPoints;
   }
   
  /*-------------Line Direction-------------------------//
    Range of methods for finding track directions - some are redundent-check!
  //----------------------------------------------*/

  XYZVec CosmicTrackFit::InitLineDirection(const ComboHit *ch0, const ComboHit *chN) {
      float tx = chN->pos().x() - ch0->pos().x();
      float ty = chN->pos().y() - ch0->pos().y();
      float tz = chN->pos().z() - ch0->pos().z();       
      XYZVec track(tx,ty,tz);
    
      return track.Unit();
    } 

  
    XYZVec CosmicTrackFit::LineDirection(double a1, double b1, const ComboHit *ch0, const ComboHit *chN, XYZVec ZPrime) {
      XYZVec track(a1,b1,1);
      return track.Unit();
    } 
  
    XYZVec CosmicTrackFit::ConvertPointToDetFrame(XYZVec vec){
        Hep3Vector vec1(vec.x(),vec.y(),vec.z());
        GeomHandle<DetectorSystem> det;
        Hep3Vector vec2 = det->toDetector(vec1);
	XYZVec XYZ(vec2.x(), vec2.y(), vec2.z());
	return XYZ;

    }

  void CosmicTrackFit::BeginFit(const char* title, CosmicTrackFinderData& TrackData){ 
	TrackData._tseed._status.clear(TrkFitFlag::helixOK); 
	bool init(false);
	if (!TrackData._tseed._status.hasAllProperties(TrkFitFlag::circleInit)) {
		init = true;
		if (initCosmicTrack( title, TrackData)){
			TrackData._tseed._status.merge(TrkFitFlag::circleInit);
		}
		else { 
			return;
		}
	} 
	if (!init)RunFitChi2(title, TrackData);
    }

   void CosmicTrackFit::RunFitChi2(const char* title, CosmicTrackFinderData& TrackData) {   
	CosmicTrack* track = &TrackData._tseed._track; 
	TrackData._tseed._status.merge(TrkFitFlag::helixOK);  
	TrackData._tseed._status.merge(TrkFitFlag::helixConverged);
	FitAll(title, TrackData, track); 
   
    }



    void CosmicTrackFit::FitAll(const char* title, CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack){
 
     ::BuildLinearFitMatrixSums S;
     ComboHit   *hitP1(0), *hitP2(0); 
     size_t nHits (trackData._chHitsToProcess.size());   
     int DOF = (nHits);// - (_Npara);
     const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
     const ComboHit* chN = &trackData._chHitsToProcess[trackData._chHitsToProcess.size()-1]; 
     cosmictrack->SetFirstHitVec(ch0->pos().x(), ch0->pos().y(), ch0->pos().z());
     cosmictrack->SetLastHitVec(chN->pos().x(), chN->pos().y(), chN->pos().z());
     cosmictrack->Set_N(nHits);

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
      }
      
      //Step 6: Begin iteration for finding the best track fit possible.
      unsigned niter(0);
      ::BuildLinearFitMatrixSums S_niteration;
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
		
		//BestTrack->clear_diag();
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

         if( abs(changed_chi2) < _max_chi2_change ){
                  converged = true;
                  cosmictrack->converged = true;     
 	      }
 	  
           if(niter == _maxniter && converged ==false ){
 		    trackData._tseed._status.clear(TrkFitFlag::helixOK);
 		    trackData._tseed._status.clear(TrkFitFlag::helixConverged);
 		    continue;
		 }
		 
   }
   cosmictrack=BestTrack;
   if(cosmictrack->converged and  _diag > 0){
	 for (size_t f5=0; f5<nHits; ++f5){
     		hitP2 = &trackData._chHitsToProcess[f5];
                if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;
		XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
		XYZVec point_prime(point.Dot(BestTrack->TrackCoordSystem._XDoublePrime), point.Dot(BestTrack->TrackCoordSystem._YDoublePrime), point.Dot(BestTrack->TrackCoordSystem._ZPrime));
		
		std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, BestTrack->TrackCoordSystem._XDoublePrime, BestTrack->TrackCoordSystem._YDoublePrime);  
		
	      	float newRx = ParametricFit::GetResidualX(BestTrack->FitParams.A0, BestTrack->FitParams.A1, point_prime);
		float newRy = ParametricFit::GetResidualY(BestTrack->FitParams.B0, BestTrack->FitParams.B1, point_prime);
		
		if(abs(newRx)/ErrorsXY[0] > _maxpull or abs(newRy)/ErrorsXY[1] > _maxpull ){ cosmictrack->converged = false;}

		BestTrack->SetCovarience(S_niteration.GetCovX()[0][0],S_niteration.GetCovX()[0][1], S_niteration.GetCovX()[1][0], S_niteration.GetCovX()[1][1], S_niteration.GetCovY()[0][0], S_niteration.GetCovY()[0][1], S_niteration.GetCovY()[1][0], S_niteration.GetCovY()[1][1]);
		
	      	cosmictrack = BestTrack;  
                     		
      	}	                	
    }
     S.clear();
     
     if(cosmictrack->converged == true){
	  
          ConvertFitToDetectorFrame(trackData, cosmictrack->TrackCoordSystem, cosmictrack->GetTrackPosition(), cosmictrack->GetTrackDirection(), cosmictrack, true, false);
	 
     } 
    
   }

/*------------Translate fit back into XYZ and/or the detector frame------//
Using matrices to ctransform from local to global coordinates
//-----------------------------------------------------------------------*/
void CosmicTrackFit::ConvertFitToDetectorFrame(CosmicTrackFinderData& trackData, TrackAxes axes, XYZVec Position, XYZVec Direction, CosmicTrack* cosmictrack, bool seed, bool det){
	TMatrixD A(3,3);
	A[0][0] = axes._XDoublePrime.X();
	A[0][1] = axes._YDoublePrime.X();
	A[0][2] = axes._ZPrime.X();

	A[1][0] = axes._XDoublePrime.Y();
	A[1][1] = axes._YDoublePrime.Y();
	A[1][2] = axes._ZPrime.Y();

	A[2][0] = axes._XDoublePrime.Z();
	A[2][1] = axes._YDoublePrime.Z();
	A[2][2] = axes._ZPrime.Z();

	TMatrixD P(3,1);
	P[0][0] = Position.X();
	P[1][0] = Position.Y();
	P[2][0] = Position.Z();

	TMatrixD D(3,1);
	D[0][0] = Direction.X();
	D[1][0] = Direction.Y();
	D[2][0] = Direction.Z();

	TMatrixD PXYZ(A*P);
	TMatrixD DXYZ(A*D);

	DXYZ[0][0] = DXYZ[0][0]/ DXYZ[2][0];
	DXYZ[1][0] = DXYZ[1][0]/ DXYZ[2][0];
	DXYZ[2][0] = DXYZ[2][0]/ DXYZ[2][0];

	PXYZ[0][0] = PXYZ[0][0] - DXYZ[0][0]*PXYZ[2][0]/ DXYZ[2][0];
	PXYZ[1][0] = PXYZ[1][0]- DXYZ[1][0]*PXYZ[2][0]/ DXYZ[2][0];
	PXYZ[2][0] = PXYZ[2][0]- DXYZ[2][0]*PXYZ[2][0]/ DXYZ[2][0];
	/*
	//TMatrixD sigmaPos(3,1);
	TMatrixD Cov(3,3);
	TMatrixD At(A);
	At.T();
	for(unsigned i =0; i< 3; i++){
		for(unsigned j =0 ; j <3; j++){
		    
		    if(i==0 and j==0) {Cov[i][j] =  cosmictrack->FitParams.Covarience.sigA0; }
		    if(i==1 and j ==1) { Cov[i][j] = cosmictrack->FitParams.Covarience.sigA1;}
		    else {Cov[i][j] =0;}
		}
	}
	TMatrixD CovXYZ(A*Cov*At);
	for(unsigned i=0; i<3 ; i++){
		for(unsigned j=0; j<3; j++){
			cout<<"i= "<<i<<" j= "<<j<<" "<<Cov[i][j]<<"=======> "<<CovXYZ[i][j]<<endl;
			
		}
		
	}
	*/ 
	TMatrixD sigmaPos(3,1);
	sigmaPos[0][0] = cosmictrack->FitParams.Covarience.sigA0; 
	sigmaPos[1][0] = cosmictrack->FitParams.Covarience.sigB0;
	sigmaPos[2][0] = 0;
	TMatrixD sigmaPosXYZ(A*sigmaPos);

	TMatrixD sigmaDir(3,1);
	sigmaDir[0][0] = cosmictrack->FitParams.Covarience.sigA1; 
	sigmaDir[1][0] = cosmictrack->FitParams.Covarience.sigB1;
	sigmaDir[2][0] = 0;
	TMatrixD sigmaDirXYZ(A*sigmaDir);

	XYZVec Pos(PXYZ[0][0], PXYZ[1][0], PXYZ[2][0]);
	XYZVec Dir(DXYZ[0][0], DXYZ[1][0] , DXYZ[2][0]);

	if(seed == true){//is this the local coords at start?
		if (det == true){ // is this detector frame?
			XYZVec PosInDet = ConvertPointToDetFrame(Pos);
			TrackEquation XYZTrack(PosInDet, Dir);
			cosmictrack->SetTrackEquationXYZ(XYZTrack);
			cosmictrack->sigmaPos.SetXYZ(sigmaPosXYZ[0][0], sigmaPosXYZ[1][0], sigmaPosXYZ[2][0]);
			cosmictrack->sigmaDir.SetXYZ(sigmaDirXYZ[0][0], sigmaDirXYZ[1][0], sigmaDirXYZ[2][0]);
			
		
		} else {
			TrackEquation XYZTrack(Pos, Dir);
			cosmictrack->SetTrackEquationXYZ(XYZTrack);
			
		}}
	else{
		TrackEquation XYZTrack(Pos, Dir);
		cosmictrack->SetMinuitTrackEquation(XYZTrack);
		cosmictrack->set_fit_phi(acos(Dir.x()/Dir.Mag2()));
		cosmictrack->set_fit_theta(acos(Dir.y()/sqrt(Dir.Mag2())));
	}

    }


    bool CosmicTrackFit::goodTrack(CosmicTrack& track) 
    { 
	if(track.Diag.FinalChiTot < _max_seed_chi2) return true;
	else return false; 
    }

    bool CosmicTrackFit::use_hit(const ComboHit& thit) const 
    {
        return (!thit._flag.hasAnyProperty(_dontuseflag));
    }

    bool CosmicTrackFit::use_track(double track_length) const 
    {
	return (track_length > _maxd) ? false : true ;
    }

    void CosmicTrackFit::DriftFit(CosmicTrackFinderData& trackData, StrawResponse const& _srep ){
	 
	FitResult endresult = MinuitDriftFitter::DoFit(_diag, trackData, _srep, _tracker, _maxHitDOCA, _minnch, _maxLogL, _gaussTres, _maxTres);


	trackData._tseed._track.MinuitFitParams.A0 =  endresult.bestfit[0];//a0
	trackData._tseed._track.MinuitFitParams.A1 =  endresult.bestfit[1];//a1
	trackData._tseed._track.MinuitFitParams.B0 =  endresult.bestfit[2];//b0
	trackData._tseed._track.MinuitFitParams.B1 =  endresult.bestfit[3];//b1
	trackData._tseed._track.MinuitFitParams.T0 =  endresult.bestfit[4];//t0

	trackData._tseed._track.MinuitFitParams.deltaA0 =  endresult.bestfiterrors[0];//erra0
	trackData._tseed._track.MinuitFitParams.deltaA1 =  endresult.bestfiterrors[1];//erra1
	trackData._tseed._track.MinuitFitParams.deltaB0 =  endresult.bestfiterrors[2];//errb0
	trackData._tseed._track.MinuitFitParams.deltaB1 =  endresult.bestfiterrors[3];//errb1
	trackData._tseed._track.MinuitFitParams.deltaT0 =  endresult.bestfiterrors[4];//errt0

	if(endresult.bestfitcov.size() !=0 ){
		 TrackCov Cov(endresult.bestfitcov[0], 0., 0., endresult.bestfitcov[1], endresult.bestfitcov[2],0.,0., endresult.bestfitcov[3]);
		 trackData._tseed._track.MinuitFitParams.Covarience = Cov;
         }
         if(endresult.NLL !=0){ trackData._tseed._track.minuit_converged = true; }
	
	 XYZVec X(1,0,0);
	 XYZVec Y(0,1,0);
	 XYZVec Z(0,0,1);

	 TrackAxes XYZ(X,Y,Z);
	 trackData._tseed._track.MinuitCoordSystem = XYZ; 
         
	if(endresult.FullFitEndTimeResiduals.size() >0){
		for(unsigned i = 0; i< endresult.FullFitEndTimeResiduals.size()-1; i++){
			if( endresult.FullFitEndTimeResiduals[i] > _maxTres or isnan(endresult.FullFitEndTimeResiduals[i])==true){ 
				trackData._tseed._track.n_outliers +=1;
				trackData._tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier); 
				
		}
		}
		
	}
        if( trackData._tseed._track.n_outliers  > _n_outliers) {
		trackData._tseed._track.minuit_converged = false;
	  }
  
}

}//end namespace
