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

namespace mu2e
{
  
  CosmicTrackFit::CosmicTrackFit(fhicl::ParameterSet const& pset) :
    _Npara(pset.get<unsigned>("Npara",4)),
    _diag(pset.get<int>("diagLevel",0)),
    _mcdiag(pset.get<int>("MCdiagLevel",1)),
    _debug(pset.get<int>("debugLevel",1)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnsh(pset.get<unsigned>("minNStrawHits",2)),
    _minCHHits(pset.get<unsigned>("minCHHits",4)),
    _n_outliers(pset.get<unsigned>("_n_outliers",5)),
    _maxniter(pset.get<unsigned>("maxniter",100)),//10
    _maxpull(pset.get<float>("maxPull",3)),
    _maxd(pset.get<float>("maxd",300.0)),//max distance between hits at start of fit
    _maxDOCA(pset.get<float>("maxDOCA",1000)),//max distance of closest approach between a hit included in fit and one flag out-right as an outlier
    _maxchi2(pset.get<float>("maxchi2",10.0)) ,
    _max_chi2_change(pset.get<float>("max_chi2_change",0.5)) 
    {}

    CosmicTrackFit::~CosmicTrackFit(){}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics) {
    if(_debug>5){
    	std::cout<<"Initializing ST Fit ..."<<std::endl;
    }
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
      if(_debug>5){
      	std::cout<<" Beginning Fit ..." << std::endl;
      }
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
   DriftFit(TrackData, track); 
   TrackData._diag.CosmicTrackFitCounter += 1;//label as having a track
   diagnostics.npasses +=1;//total number of tracks
}

int n_lost = 0; //For Debug 
/*---------------Refine Fit ------ ----------------//
//    Refines the fit in and updates chi2 information   //
//-----------------------------------------------*/
void CosmicTrackFit::FitAll(const char* title, CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, CosmicTrackFinderTypes::Data_t& diagnostics){
    ::BuildMatrixSums S;
    ComboHit   *hitP1(0), *hitP2(0); 
    size_t nHits (trackData._chHitsToProcess.size());
    int DOF = (nHits);//-4;
    const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
    const ComboHit* chN = &trackData._chHitsToProcess[trackData._chHitsToProcess.size()-1]; 
    XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
    XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());
  
    //Step 1: Get Initial Estimate of track direction
    XYZVec ZPrime = InitLineDirection(ch0, chN);  
    std::vector<XYZVec> AxesList = ParametricFit::GetAxes(ZPrime);
   
    //Step 2: Loop over hits and get track parameters based on above estimated track direction
    for (size_t f1=0; f1<nHits; ++f1){  
      hitP1 = &trackData._chHitsToProcess[f1];  
      if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;  
      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1]); 
      S.addPoint(f1, point, AxesList[0], AxesList[1], AxesList[2], ErrorsXY[0], ErrorsXY[1]); //S will be  updated   
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
     AxesList = ParametricFit::GetAxes(ZPrime);
     cosmictrack->set_parameters(a0,a1,b0,b1);
     cosmictrack->set_track_direction(UpdatedTrackDirection);
 
     //Step 5: Loop for initial diagnostics
     if(_debug>0){
             cosmictrack->set_initial_parameters(a0,a1,b0,b1);
	     cosmictrack->setinitXPrime(AxesList[0]);
	     cosmictrack->setinitYPrime(AxesList[1]);
	     cosmictrack->setinitZPrime(AxesList[2]); 
	     cosmictrack->set_initchisq_dofY(S.GetChi2Y()/abs(DOF));
	     cosmictrack->set_initchisq_dofX(S.GetChi2X()/abs(DOF));
	     for (size_t f2=0; f2<nHits; ++f2){
	              if(isnan(cosmictrack->get_track_direction().Mag2()) == true) continue;     
	              hitP1 = &trackData._chHitsToProcess[f2];  
		      if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;  
		      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
		      XYZVec point_prime(point.Dot(AxesList[0]), point.Dot(AxesList[1]), point.Dot(ZPrime));
		      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1]);   
		      float Rx = ParametricFit::GetResidualX(a0,a1, point_prime);
		      float Ry = ParametricFit::GetResidualY(b0, b1, point_prime); 
		      cosmictrack->set_init_hit_errorsTotal(sqrt(pow(ErrorsXY[0],2)+pow(ErrorsXY[1],2)));
		      cosmictrack->set_init_fit_residual_errorsX(ErrorsXY[0]);
		      cosmictrack->set_init_fit_residual_errorsY(ErrorsXY[1]);
		      cosmictrack->set_init_pullsX(Rx/ErrorsXY[0]);
		      cosmictrack->set_init_pullsY(Ry/ErrorsXY[1]);
		      cosmictrack->set_init_fit_residualsX(Rx);
		      cosmictrack->set_init_fit_residualsY(Ry); 
		} 
      }
     //Step 6: Begin iteration for finding the best track fit possible.
     unsigned niter(0);
     ::BuildMatrixSums S_niteration;
     bool converged = false;
     //bool track_improved = false;
     CosmicTrack* BestTrack = cosmictrack;	 
     double chi2_best_track = 10000000;//chosen arbitary high number
     
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
     	S_niteration.clear(); 
     	
     	AxesList = ParametricFit::GetAxes(cosmictrack->get_track_direction());
        cosmictrack->set_niter(niter );
     	for (size_t f4=0; f4<nHits; ++f4){
     	      if(isnan(cosmictrack->get_track_direction().Mag2()) == true) continue;     
     	      hitP2 = &trackData._chHitsToProcess[f4];
      	      if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;   
	      XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());	    
	      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1]);  
      	      XYZVec point_prime(point.Dot(AxesList[0]), point.Dot(AxesList[1]), point.Dot(AxesList[2]));
	      S_niteration.addPoint(f4, point, AxesList[0], AxesList[1], AxesList[2], ErrorsXY[0],ErrorsXY[1]);
    	  }
           a0 = S_niteration.GetAlphaX()[0][0];
           a1 = S_niteration.GetAlphaX()[1][0];
	   b0 = S_niteration.GetAlphaY()[0][0];
	   b1 = S_niteration.GetAlphaY()[1][0];
	   XYZVec DirectionSecond(a1,b1,1);
	   XYZVec UpdatedTrackDirectionSecond =DirectionSecond.Unit();
           cosmictrack->set_track_direction(UpdatedTrackDirectionSecond); 
	   cosmictrack->setXPrime(AxesList[0]);
	   cosmictrack->setYPrime(AxesList[1]);
	   cosmictrack->setZPrime(AxesList[2]);
	   cosmictrack ->set_parameters(a0,a1,b0,b1);
           float updated_chi2 = S_niteration.GetTotalChi2()/abs(DOF);
           changed_chi2 = chi2_best_track - updated_chi2;
           if(_debug > 0){ //need to fill with something to avoid "0" peak (correct this later!)
		   cosmictrack->set_finalchisq_dof(S_niteration.GetTotalChi2()/abs(DOF)); 
		   cosmictrack->set_finalchisq_dofY(S_niteration.GetChi2Y()/abs(DOF));
		   cosmictrack->set_finalchisq_dofX(S_niteration.GetChi2X()/abs(DOF));
	   }
           if (updated_chi2 < chi2_best_track ){
           		//track_improved = true;
           		BestTrack->clear_all();
			BestTrack ->set_parameters(a0,a1,b0,b1);
			BestTrack->set_track_direction(cosmictrack->get_track_direction()); 
           		BestTrack->setXPrime(AxesList[0]);
	   		BestTrack->setYPrime(AxesList[1]);
	   		BestTrack->setZPrime(AxesList[2]);
           		cosmictrack = BestTrack;
                        if(_debug > 0){
                        	if(S.GetTotalChi2()/abs(DOF)<20){
                        	BestTrack->set_chi2_quant(PDF(S_niteration.GetTotalChi2()/abs(DOF), DOF));
                        	}
		                BestTrack->set_finalchisq_dof(S_niteration.GetTotalChi2()/abs(DOF)); 
		                BestTrack->set_finalchisq_dofY(S_niteration.GetChi2Y()/abs(DOF));
		                BestTrack->set_finalchisq_dofX(S_niteration.GetChi2X()/abs(DOF));   
		                
		        }
                        chi2_best_track = updated_chi2;  
                               	
        	 }
	     
             //If on any iteration the change in chi2 becomes small than we can assume we have converged on a track
             if( abs(changed_chi2) < _max_chi2_change ){
                 converged = true;
	      }
	      
	      //If at end of the iteration process but the track is still not converging then remove it
              if(niter == _maxniter && converged ==false ){
              		 n_lost+=1;
		    	 trackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
		    	 trackData._tseed._status.clear(TrkFitFlag::StraightTrackConverged);
		    	 continue;
    		 }
    		      	
	        }//end while 
	        cosmictrack=BestTrack;
		     //step 7 : final diagnostics stored
	             if(_debug>0){
	             
		     for (size_t f5=0; f5<nHits; ++f5){
		     		if(isnan(cosmictrack->get_track_direction().Mag2()) == true) continue;     
		                hitP2 = &trackData._chHitsToProcess[f5];
      	        		XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
      	        		XYZVec point_prime(point.Dot(BestTrack->getXPrime()), point.Dot(BestTrack->getYPrime()), point.Dot(BestTrack->getZPrime()));
      	        		
      	        		std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1]);  
				float newRx = ParametricFit::GetResidualX(BestTrack->get_parameter(0), BestTrack->get_parameter(1), point_prime);
				float newRy = ParametricFit::GetResidualY(BestTrack->get_parameter(2), BestTrack->get_parameter(3), point_prime);
				
				BestTrack->set_final_hit_errorsTotal(sqrt(pow(ErrorsXY[0],2)+ pow(ErrorsXY[1],2)));
				BestTrack->set_track_position(BestTrack->get_parameter(0),BestTrack->get_parameter(2),0); //in z'=0 plane
				BestTrack->set_final_fit_residual_errorsX(ErrorsXY[0]);
		      		BestTrack->set_final_fit_residual_errorsY(ErrorsXY[1]);
		      		BestTrack->set_final_pullsX(newRx/ErrorsXY[0]);
		      		BestTrack->set_final_pullsY(newRy/ErrorsXY[1]);
		      		BestTrack->set_final_fit_residualsX(newRx);
		      		BestTrack->set_final_fit_residualsY(newRy); 
		      		BestTrack->set_initchisq_dof(S.GetTotalChi2()/abs(DOF));
		      		BestTrack->set_cov(S.GetCovX()[0][0],S.GetCovX()[1][1], S.GetCovY()[0][0], S.GetCovY()[1][1]);  
		      		cosmictrack=BestTrack;      			
                   	  }               	
                   	}
    S.clear();
    if(_mcdiag > 0){
         FitMC(trackData,cosmictrack, true, false);
     }
    //cosmictrack->set_fit_phi(acos(cosmictrack->getZPrime().x()/sqrt(cosmictrack->getZPrime().Mag2())));
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
//Applies fitting routine to MCDigis:
void CosmicTrackFit::FitMC(CosmicTrackFinderData& trackData, CosmicTrack* cosmictrack, bool XYZ, bool prime){	
        ::BuildMatrixSums S;
         //GeomHandle<DetectorSystem> det;        
    	XYZVec FitDirection = cosmictrack->getZPrime();
    	
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
                std::vector<XYZVec> AxesList = ParametricFit::GetAxes(ZPrime);
                S.addPoint(f1, point, AxesList[0], AxesList[1], AxesList[2], 1,1);//AxesList[0], AxesList[1], AxesList[2], 1,1); 
            }
            if(XYZ == true && prime ==false){ //gets the true XYZ and does in that frame USE THIS THEN TRANSFORM
                XYZVec point(posN.x(), posN.y(), posN.z());
                std::vector<XYZVec> AxesList = ParametricFit::GetAxes(ZPrime);
                XYZVec X(1,0,0);//to evaluate for XYZ frame
                XYZVec Y(0,1,0);
                XYZVec Z(0,0,1);
                S.addPoint(f1, point, X,Y,Z, 1,1);//AxesList[0], AxesList[1], AxesList[2], 1,1); 
            }
            if(XYZ == false && prime ==true){//true for calc in X"Y"Z' frame where these are from reco
                std::vector<XYZVec> AxesList = ParametricFit::GetAxes(FitDirection); 
                XYZVec point(posN.Dot(AxesList[0]), posN.Dot(AxesList[1]), posN.Dot(AxesList[2]));
                S.addPoint(f1, point, AxesList[0], AxesList[1], AxesList[2], 1,1); 
            }
          
        }   
     //now transform into XYZ fitplane
     double a0 = S.GetAlphaX()[0][0];
     double a1 = S.GetAlphaX()[1][0];
     double b0 = S.GetAlphaY()[0][0];
     double b1 = S.GetAlphaY()[1][0];
     XYZVec Direction(a1,b1,1);
     XYZVec Position(a0,b0,0);
    // XYZVec UpdatedTrackDirection =Direction.Unit();
     XYZVec transformedDirection(Direction.Dot(cosmictrack->getXPrime())/Direction.Dot(cosmictrack->getZPrime()), Direction.Dot(cosmictrack->getYPrime())/Direction.Dot(cosmictrack->getZPrime()), 1);
     XYZVec unit_transformedDirection = transformedDirection.Unit();
     XYZVec transformedPosition(Position.Dot(cosmictrack->getXPrime()), Position.Dot(cosmictrack->getYPrime()), Position.Dot(cosmictrack->getZPrime()));
     
     double a0_trans = transformedPosition.X();
     double a1_trans = unit_transformedDirection.X();
     double b0_trans = transformedPosition.Y();
     double b1_trans = unit_transformedDirection.Y();
     
     double mcpos_mag = sqrt(Direction.Mag2());    
     cosmictrack->set_true_track_direction(ZPrime);
     //cosmictrack->set_true_phi(acos(ZPrime.x()/mcpos_mag)); 
     const double theta_start = acos(Direction.x()/mcpos_mag);
     const double phi_start = atan(Direction.y()/Direction.x());
     cosmictrack->set_true_phi(phi_start);
     cosmictrack->set_true_theta(theta_start);
     cosmictrack->set_true_parameters(a0_trans,a1_trans,b0_trans,b1_trans);
     cosmictrack->set_true_finalchisq_dof(S.GetTotalChi2()/nHits);
     
     }


XYZVec CosmicTrackFit::MCInitHit(StrawDigiMC mcdigi){
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(StrawEnd::cal);
	XYZVec first(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
        return first;
}
/*

//GeomHandle<DetectorSystem> det;
Hep3Vector pos = det->toDetector(spmcp->position());
double mcposx0 = spmcp->position().x();
double mcposy0 = spmcp->position().y();
double mcposz0 = spmcp->position().z();
*/
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
      trackData._tseed._track.set_true_track_direction(track.Unit());
   
}

bool CosmicTrackFit::goodTrack(CosmicTrack* track)
  { 
    if(track->get_finalchisq_dof() < _maxchi2) return true;
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
void CosmicTrackFit::DriftFit(CosmicTrackFinderData& trackData, CosmicTrack* cosmictrack){
         std::vector<XYZVec> positions;
         std::vector<double> times;
         std::vector<Straw> straws;
         std::vector<double> errorsX;
         std::vector<double> errorsY;
	 size_t nHits (trackData._chHitsToProcess.size());
         ComboHit   *hitP1(0);
         //double t0 = trackData._tseed._t0.t0() ; //Get t0 from tclust, from BTrk
         //double t0err = trackData._tseed._t0.t0Err() ; //error on t0 from BTrk
         cout<<"=================="<<endl;
         for (size_t f1=0; f1<nHits+1; ++f1){  
      		hitP1 = &trackData._chHitsToProcess[f1];
      		if (!use_hit(*hitP1) )    continue;
      		float hit_time = hitP1->time();
      		times.push_back(hit_time);
      		XYZVec pos(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      		positions.push_back(pos);
      		XYZVec track_direction = trackData._tseed.track().get_track_direction();
		std::vector<XYZVec> AxesList = ParametricFit::GetAxes(track_direction);
      		errorsX.push_back(ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1])[0]);
      		errorsY.push_back(ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1])[1]);   
      		//Straw const&  straw = _tracker->getStraw(hitP1->strawId());
      		//straws.push_back(straw);
      		/*
      	        Straw const&  straw = _tracker->getStraw(hitP1->strawId());
      	        XYZVec track_position(trackData._tseed.track().get_parameter(0), trackData._tseed.track().get_parameter(2), 0);
                XYZVec track_direction = trackData._tseed.track().get_track_direction();
		std::vector<XYZVec> TrackAxes = ParametricFit::GetAxes(track_direction);
		
      	        TrkPoca poca = DriftFitUtils::GetPOCA(straw, TrackAxes, track_position,  track_direction);
      	        double doca =  DriftFitUtils::GetDOCA(poca);
      	        cout<<"DOCA"<<doca<<endl;
      		
      		StrawId strawid = straw.id();
      		double phi =  0.;//LiklihoodFunctions::GetPhi(straw, track_position,  TrackAxes, doca);
      		//Use the StrawResponse::driftDistanceToTime() function to turn the transverse distance and t0 into an expected time that will let you calculate a time residual:
      	        double expectedtime= _srep.StrawResponse::driftDistanceToTime(strawid , fabs(doca), phi);
      	        //double time_residual_trans = expectedtime;
      	        double wirehalflen = straw.halfLength();              
      	        double propagation_time = wirehalflen/299.; 
      	        double time_residual_long = expectedtime+propagation_time;
                
                //Fill lists for Minuit:
                DOCAs.push_back(doca);
                times.push_back(hit_time);
                time_residuals.push_back(time_residual_long);
      	        */
      		}
      		
      		EndResult bestfit = LiklihoodFunctions::DoFit(trackData._tseed, times, positions, errorsX, errorsY, straws);
      		

}
 

}//end namespace
