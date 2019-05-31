// Object to perform track fit to combo hits
//
// $Id: CosmicTrackFit code
// $Author: S Middleton
// $Date: Feb 2019
//
// Mu2e Cosmics:
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"
//MC:
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
#include "RecoDataProducts/inc/DriftCircle.hh"

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
//Least Squares Fitter:

#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/BuildMatrixSums.hh"

//ROOT:
#include "TMatrixD.h"
#include "Math/VectorUtil.h"
#include "TH1F.h"
#include "TPolyMarker.h"

// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>




using namespace std;
using namespace boost::accumulators;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace mu2e
{
  CosmicTrackFit::CosmicTrackFit(fhicl::ParameterSet const& pset) :
    _Npara(pset.get<unsigned>("Npara",4)),
    _diag(pset.get<int>("diagLevel",0)),
    _mcdiag(pset.get<int>("MCdiagLevel",0)),
    _debug(pset.get<int>("debugLevel",1)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnsh(pset.get<unsigned>("minNStrawHits",2)),
    _minCHHits(pset.get<unsigned>("minCHHits",4)),
    _n_outliers(pset.get<unsigned>("_n_outliers",5)),
    _maxniter(pset.get<unsigned>("maxniter",100)),//10
    _maxpull(pset.get<double>("maxPull",5)),
    _maxd(pset.get<float>("maxd",300.0)),//max distance between hits at start of fit
    _maxDOCA(pset.get<float>("maxDOCA",1000)),//max distance of closest approach between a hit included in fit and one flag out-right as an outlier
    _maxchi2(pset.get<float>("maxchi2",1000.0)) ,
    _max_chi2_change(pset.get<float>("max_chi2_change",0.1)) 
    {}

//destructor
    CosmicTrackFit::~CosmicTrackFit(){}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics) {
    if(_debug>0){
    	std::cout<<"Initializing ST Fit ..."<<std::endl;
    }
    bool is_ok(false);
    RunFitChi2(title, TrackData,diagnostics );
    is_ok = TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackOK);
    return is_ok;
  }


  /*-------------Init Line-------------------------//
 Makes basic estimate of a line y=mx+c style with Cosmic line between first and last points on plane
  //----------------------------------------------*/
  XYZVec CosmicTrackFit::InitLineDirection(const ComboHit *ch0, const ComboHit *chN) {
     
      double tx = chN->pos().x() - ch0->pos().x();
      double ty = chN->pos().y() - ch0->pos().y();
      double tz = chN->pos().z() - ch0->pos().z();
       
      XYZVec track(tx,ty,tz);
   
      return track.Unit();
    } 
  
    XYZVec CosmicTrackFit::LineDirection(double a1, double b1, const ComboHit *ch0, const ComboHit *chN, XYZVec ZPrime) {
      
     // double tx = a1*(chN->pos().Dot(ZPrime) - ch0->pos().Dot(ZPrime));
      //double ty = b1*(chN->pos().Dot(ZPrime) - ch0->pos().Dot(ZPrime));
     // double tz = chN->pos().Dot(ZPrime) - ch0->pos().Dot(ZPrime);
      
      XYZVec track(a1,b1,1);//tx,ty,tz);
      return track.Unit();
    } 
    
    XYZVec CosmicTrackFit::GetTrackDirection(std::vector<XYZVec> hitXYZ, XYZVec XDoublePrime, XYZVec YDoublePrime, XYZVec ZPrime){
           std::vector<double> X, Y, Z;
    	   for (size_t j =0; j < hitXYZ.size(); j++){
    	  	 X.push_back( hitXYZ[j].Dot(XDoublePrime));
    	  	 Y.push_back(hitXYZ[j].Dot(YDoublePrime));
    	  	 Z.push_back(hitXYZ[j].Dot(ZPrime));
  	   }
           double minz = *std::min_element(Z.begin(), Z.end());
           double maxz = *std::max_element(Z.begin(),Z.end());
           unsigned first = 0;
           unsigned last = Z.size()-1;
           
           for(size_t i = 0 ; i < Z.size(); i++){
           	double zprime = Z.at(i);
           	if(zprime == maxz) first = i;
           	if(zprime == minz) last = i;
           	if (first >  last) {
           		unsigned temp;
           		temp = first;
           		first = last;
           		last = temp;
           	}
           	
           }
           
           double tx = (X.at(last) - Z.at(first));
      	   double ty = (Y.at(last) - Y.at(first));
           double tz = Z.at(last) - Z.at(first);
           
           XYZVec track(tx,ty,tz);
          
           return track.Unit();
      
    }
    XYZVec CosmicTrackFit::TransformTrueDirection(XYZVec TrueDirection, XYZVec ZPrime, XYZVec XDoublePrime, XYZVec YDoublePrime){
    return TrueDirection;
    
    }
//--------------Fit-----------------//
//Top call to Fitting routines....
//-------------------------------------------// 
  void CosmicTrackFit::BeginFit(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics){
      if(_debug>0){
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
   TrackData._diag.CosmicTrackFitCounter += 1;//label as having a track
   diagnostics.npasses +=1;//totsl number of tracks
   
}

/*---------------Refine Fit ------ ----------------//
//    Refines the fit in and updates chi2 information   //
//-----------------------------------------------*/
void CosmicTrackFit::FitAll(const char* title, CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, CosmicTrackFinderTypes::Data_t& diagnostics){
   
    ::BuildMatrixSums S, S_init, S_final;
    
    ComboHit   *hitP1(0), *hitP2(0); 
    size_t nHits (trackData._chHitsToProcess.size());
    int DOF = (nHits);
    const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
    const ComboHit* chN = &trackData._chHitsToProcess[trackData._chHitsToProcess.size()-1]; 
    XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
    XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());
  
    //Step 1: Get Initial Estimate of track direction
    XYZVec ZPrime = InitLineDirection(ch0, chN);  
    //cosmictrack->set_track_direction(ZPrime);
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
	     cosmictrack->set_initchisq_dof(S.GetTotalChi2()/abs(DOF)); 
	     cosmictrack->set_initchisq_dofY(S.GetChi2Y()/abs(DOF));
	     cosmictrack->set_initchisq_dofX(S.GetChi2X()/abs(DOF));
	     for (size_t f3=0; f3<nHits; ++f3){ 
		      std::cout<<"hit "<<f3<<" of "<<nHits<<std::endl;
		      if(isnan(cosmictrack->get_track_direction().Mag2()) == true) continue;     
		      hitP1 = &trackData._chHitsToProcess[f3];  
		      if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;  
		      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
		      XYZVec point_prime(point.Dot(AxesList[0]), point.Dot(AxesList[1]), point.Dot(ZPrime));
		      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1]);     
		     
		      double Rx = ParametricFit::GetResidualX(a0,a1, point_prime);
		      double Ry = ParametricFit::GetResidualY(b0, b1, point_prime); 
		      std::cout<<"Pull X'' "<<Rx/ErrorsXY[0]<<std::endl;
		      std::cout<<"Pull Y'' "<<Ry/ErrorsXY[1]<<std::endl;
		       
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
     //bool improved = false;
     CosmicTrack* BestTrack = cosmictrack;	 
     double chi2_best_track = 10000;
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
	
           double updated_chi2 = S_niteration.GetTotalChi2()/abs(DOF);
           changed_chi2 = chi2_best_track - updated_chi2;
          
           if (updated_chi2 < chi2_best_track ){
           		//improved =true;
           		BestTrack->clear_all();
           		std::cout<<"Improved....."<<std::endl;
			BestTrack ->set_parameters(a0,a1,b0,b1);
			BestTrack->set_track_direction(cosmictrack->get_track_direction()); 
			
           		BestTrack->setXPrime(AxesList[0]);
	   		BestTrack->setYPrime(AxesList[1]);
	   		BestTrack->setZPrime(AxesList[2]);
	   		 std::cout<<"Improved Z' "<<" "<<BestTrack->getZPrime()<<" X'' "<<BestTrack->getXPrime()<<" Y'' "<<BestTrack->getYPrime()<<std::endl;
           		cosmictrack = BestTrack;
                        if(_debug > 0){
                                std::cout<<"Filling chi2 with "<<S_niteration.GetTotalChi2()/abs(DOF)<<std::endl;
		                BestTrack->set_finalchisq_dof(S_niteration.GetTotalChi2()/abs(DOF));
		                BestTrack->set_finalchisq_dofY(S_niteration.GetChi2Y()/abs(DOF));
		                BestTrack->set_finalchisq_dofX(S_niteration.GetChi2X()/abs(DOF)); 
		                
		        }
                        
                        chi2_best_track = updated_chi2;                	
        	 }
	     
             //If on any iteration the change in chi2 becomes small than we can assume we have converged on a track
             if( abs(changed_chi2) < _max_chi2_change){
                 converged = true;
	      }
	      //If at end of the iteration process but the track is still not converging then remove it
              if(niter == _maxniter && converged ==false ){
		    	 trackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
		    	 trackData._tseed._status.clear(TrkFitFlag::StraightTrackConverged);
    		 }
    			
	        }//end while 
	        cosmictrack=BestTrack;
	        std::cout<<"Final axes Z' = "<<BestTrack->getZPrime()<<" Y'' = "<<BestTrack->getYPrime()<<" X'' = "<<BestTrack->getXPrime()<<std::endl;
	             if(_debug>0){
	             
		     for (size_t f5=0; f5<nHits; ++f5){
		     	          if(isnan(cosmictrack->get_track_direction().Mag2()) == true) continue;     
	     	      		hitP2 = &trackData._chHitsToProcess[f5];
      	        		XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
      	        		XYZVec point_prime(point.Dot(BestTrack->getXPrime()), point.Dot(BestTrack->getYPrime()), point.Dot(BestTrack->getZPrime()));
      	        		std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, AxesList[0], AxesList[1]);  
		        	
				double newRx = ParametricFit::GetResidualX(BestTrack->get_parameter(0), BestTrack->get_parameter(1), point_prime);
				double newRy = ParametricFit::GetResidualY(BestTrack->get_parameter(2), BestTrack->get_parameter(3), point_prime);
				  
				std::cout<<"Hit "<<f5<<" of "<<nHits<<std::endl;
				std::cout<<"Pull X'' "<<newRx/ErrorsXY[0]<<std::endl;
				std::cout<<"Pull Y'' "<<newRy/ErrorsXY[0]<<std::endl;
				BestTrack->set_final_hit_errorsTotal(sqrt(pow(ErrorsXY[0],2)+ pow(ErrorsXY[1],2)));
				BestTrack->set_final_fit_residual_errorsX(ErrorsXY[0]);
		      		BestTrack->set_final_fit_residual_errorsY(ErrorsXY[1]);
		      		BestTrack->set_final_pullsX(newRx/ErrorsXY[0]);
		      		BestTrack->set_final_pullsY(newRy/ErrorsXY[1]);
		      		BestTrack->set_final_fit_residualsX(newRx);
		      		BestTrack->set_final_fit_residualsY(newRy);                     	
                   	  }               	
                   	}
    S.clear();
    S_init.clear();

  }
  
  void CosmicTrackFit::Draw_Chi_Para( const char* title, TH2F *cha0, TH2F *cha1, TH2F *chb0, TH2F *chb1, double end_chi2, double a0, double a1, double b0, double b1){
        TCanvas canvas;
        canvas.Draw();
        canvas.Divide(2,2);
  	TPolyMarker a0Result( 1, &a0,&end_chi2);
	TPolyMarker a1Result( 1, &a1,&end_chi2);
	TPolyMarker b0Result( 1, &b0,&end_chi2);
	TPolyMarker b1Result( 1, &b1,&end_chi2);
	a0Result.SetMarkerColor(2);
	a0Result.SetMarkerStyle(5);
	a0Result.SetMarkerSize(3);
      	a1Result.SetMarkerColor(2);
      	a1Result.SetMarkerStyle(5);
      	a1Result.SetMarkerSize(3);
      	b0Result.SetMarkerColor(2);
      	b0Result.SetMarkerStyle(5);
      	b0Result.SetMarkerSize(3);
      	b1Result.SetMarkerColor(2);
      	b1Result.SetMarkerStyle(5);
      	b1Result.SetMarkerSize(3);
      	canvas.cd(1);
	cha0->SetMarkerStyle(5);
	cha0->SetMarkerSize(3);
	cha0->Draw();
	a0Result.Draw("same");
	canvas.cd(2);
	cha1->SetMarkerStyle(5);
	cha1->SetMarkerSize(3);
	cha1->Draw();
	a1Result.Draw("same");

	canvas.cd(3);
	chb0->SetMarkerStyle(5);
	chb0->SetMarkerSize(3);
	chb0->Draw();
	b0Result.Draw("same");
	canvas.cd(4);
	chb1->SetMarkerStyle(5);
	chb1->SetMarkerSize(3);
	chb1->Draw();
	b1Result.Draw("same");
	canvas.Modified();
        canvas.Update();
        canvas.SaveAs(title);
	
	
	}
//A Test of the Chi-2 minimization part of the code: (works)
  void CosmicTrackFit::FitXYZ(CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, CosmicTrackFinderTypes::Data_t& diagnostics){
    std::cout<<" -------------------------------"<<std::endl;
    ::BuildMatrixSums S, S_init;
    
    ComboHit   *hitP1(0);//, *hitP2(0); 
    size_t nHits (trackData._chHitsToProcess.size());
   
    const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
    const ComboHit* chN = &trackData._chHitsToProcess[trackData._chHitsToProcess.size()-1];
    
    XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
    XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());
    XYZVec ZPrime(0,0,1);
   
    cosmictrack->set_initial_track_direction(ZPrime);
    cosmictrack->set_track_direction(ZPrime);  
    
    XYZVec XDoublePrime(1,0,0);//ParametricFit::GetXDoublePrime(XPrime, YPrime, ZPrime);
    XYZVec YDoublePrime(0,1,0);//ParametricFit::GetYDoublePrime(XPrime, YPrime, ZPrime);
    //First loop to update track information and hit error estimate
    for (size_t f1=0; f1<nHits; ++f1){    
      hitP1 = &trackData._chHitsToProcess[f1];   
      if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;  
      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      XYZVec major_axis =  ParametricFit::MajorAxis(hitP1);
      XYZVec minor_axis =  ParametricFit::MinorAxis(hitP1);
      double errX =  ParametricFit::HitErrorX(hitP1, major_axis, minor_axis, XDoublePrime);
      double errY =  ParametricFit::HitErrorY(hitP1, major_axis, minor_axis, YDoublePrime);   
      S.addPoint(f1, point, XDoublePrime, YDoublePrime, ZPrime, errX, errY); //S will be  updated
      S_init.addPoint(f1, point, XDoublePrime, YDoublePrime, ZPrime, errX, errY);//S_init is store for reference     
     }
    
     double a0 = S_init.GetAlphaX()[0][0];
     double a1 = S_init.GetAlphaX()[1][0];
     double b0 = S_init.GetAlphaY()[0][0];
     double b1 = S_init.GetAlphaY()[1][0];
     cosmictrack->set_initial_parameters(a0,a1,b0,b1);
     cosmictrack->set_parameters(a0,a1,b0,b1);
     XYZVec updated_track_direction(a1,b1,1);  
     cosmictrack->set_track_direction(updated_track_direction.Unit());
     cosmictrack->setXPrime(XDoublePrime);
     cosmictrack->setYPrime(YDoublePrime);
     cosmictrack->setZPrime(ZPrime);
    
    S.clear();
    S_init.clear();
    
  }
  
void CosmicTrackFit::FitMC(CosmicTrackFinderData& trackData){	
        ::BuildMatrixSums SMC;
         GeomHandle<DetectorSystem> det;        
    	//size_t nHits (trackData._mccol->size());
    	size_t nHits (trackData._mcDigisToProcess.size());
        StrawDigiMC const& ch0 = trackData._mccol->at(0); 
        StrawDigiMC const& chN = trackData._mccol->at(nHits-1);
        //for(size_t i=0;i<nHits;i++) {
        //	StrawDigiMC const& chi = trackData._mccol->at(i);
        //	art::Ptr<StepPointMC> const& spmcpi = chi.stepPointMC(StrawEnd::cal);
        
        	
        //}
        art::Ptr<StepPointMC> const& spmcp0 = ch0.stepPointMC(StrawEnd::cal);
        
	Hep3Vector pos0 = det->toDetector(spmcp0->position());
	
	
        art::Ptr<StepPointMC> const& spmcpN = chN.stepPointMC(StrawEnd::cal);
        
        Hep3Vector posN = det->toDetector(spmcpN->position());
        
        double tx = posN.x() -  pos0.x();
        double ty = posN.y() -  pos0.y();
        double tz = posN.z() -  pos0.z();
        XYZVec track(tx,ty,tz);
        XYZVec ZPrime = track.Unit() ;
        //std::cout<<"Init MC Direction "<<ZPrime<<std::endl;
        XYZVec XPrime = ParametricFit::GetXPrime(ZPrime);
        XYZVec YPrime = ParametricFit::GetYPrime(XPrime, ZPrime);   
        XYZVec XDoublePrime =ParametricFit::GetXDoublePrime(XPrime, YPrime, ZPrime);
        XYZVec YDoublePrime = ParametricFit::GetYDoublePrime(XPrime, YPrime, ZPrime);
        //std::cout<<"MC Fitter "<<ZPrime<<" "<<nHits<<std::endl;
        for (size_t f1=0; f1<nHits-1; ++f1){    
          //std::cout<<"getting MC parameters "<<f1<<" of "<<nHits<<" length "<<trackData._mcDigisToProcess.size()<<std::endl;
	  StrawDigiMC* hitP1 = &trackData._mcDigisToProcess[f1];  
	  art::Ptr<StepPointMC> const& spmcp = hitP1->stepPointMC(StrawEnd::cal);
	  Hep3Vector pos = det->toDetector(spmcp->position());
	  XYZVec point(pos.x(),pos.y(),pos.z());
            
          SMC.addPoint(f1, point, XDoublePrime, YDoublePrime, ZPrime, 20,20); //S will be  updated
       }
       /*
     double a0 = SMC.GetAlphaX()[0][0];
     double a1 = SMC.GetAlphaX()[1][0];
     double b0 = SMC.GetAlphaY()[0][0];
     double b1 = SMC.GetAlphaY()[1][0];
     */
}


XYZVec CosmicTrackFit::MCInitHit(StrawDigiMC mcdigi){
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(StrawEnd::cal);
	//std::cout<<"first before "<<)<<" "<<<<" "<<spmcp->position().z()<<std::endl;
	//GeomHandle<DetectorSystem> det;
	//Hep3Vector pos = det->toDetector(spmcp->position());

	//double mcposx0 = spmcp->position().x();
	//double mcposy0 = spmcp->position().y();
	//double mcposz0 = spmcp->position().z();
	XYZVec first(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
	
        return first;
}

XYZVec CosmicTrackFit::MCFinalHit(StrawDigiMC mcdigi){
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(StrawEnd::cal);
	//GeomHandle<DetectorSystem> det;
	//Hep3Vector pos = det->toDetector(spmcp->position());
	//double mcposxN = spmcp->position().x();
	//double mcposyN = spmcp->position().y();
	//double mcposzN = spmcp->position().z();
	XYZVec last(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
	
        return last;
}


void CosmicTrackFit::MCDirection(XYZVec first, XYZVec last, CosmicTrackFinderData& trackData){ //int hitN, int N, CosmicTrackFinderData& trackData, StrawDigiMC mcdigi){
       
      //CosmicTrack* cosmictrack = &trackData._tseed._track; 
      double tx = last.x() - first.x();
      double ty = last.y() - first.y();
      double tz = last.z() - first.z();       
      XYZVec track(tx,ty,tz);
      
      trackData._tseed._track.setMCDirection(track.Unit());
      
        
}

void MulitpleTrackResolver(CosmicTrackFinderData& trackData,CosmicTrack* track){
	//if track has a significant amount of tracks with large chi2 and those large values are all similar...and resiudals within range of each other, this coule mean a second track....check for this although unlikely with cosmics....

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
  

  
 /*--------Drift Correction-------*/
void CosmicTrackFit::DriftCorrection(CosmicTrackFinderData& trackData){
	//Find SH in CH with smallest distance from wire centre, then need to find relative sign of that SH doca value. Then we know which side of wire the track approached. Then add on/minus off to z this distance
	
         size_t nHits (trackData._chHitsToProcess.size());
         signed ambig;
         ComboHit   *hitP1(0);
         const Straw*  straw;
         //TrkDef& mydef;
         //const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
         //const ComboHit* chN = &trackData._chHitsToProcess[nHits+1];
         //XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
         //XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());
         for (size_t f1=0; f1<nHits+1; ++f1){  
      		hitP1 = &trackData._chHitsToProcess[f1];
      		if (!use_hit(*hitP1) )    continue;
      		//XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      	        straw     = &_tracker->getStraw(hitP1->strawId());
      		const CLHEP::Hep3Vector& wpos = straw->getMidPoint();
      		const CLHEP::Hep3Vector& wdir = straw->getDirection();

      		HepPoint pointwire(wpos.x(),wpos.y(),wpos.z());
      		TrkLineTraj htraj(pointwire,wdir,-straw->halfLength(),straw->halfLength());
      		XYZVec tpos = trackData._tseed._track.get_track_position();
      		XYZVec tdir = trackData._tseed._track.get_track_direction();
      		double fltlen = (wpos.z()-tpos.z())/tdir.z();
	 	
	        if(fltlen < 0){
	        	 ambig =  -1;
	 
	        }	
	        else{
	        	 ambig =  1;
	        }
	        std::cout<<"AMBIG "<<ambig<<"Flt Len "<<fltlen<<std::endl;
	 }
	 
	}
        /*
	TrkT0 t0 = trackData._tseed._t0 ; //Get t0 from tclust
	//double            dz_max(1.e12) ; // closest_z(1.e12);

        HepPoint   tpos; //position (c?)
        //double     dt; //time difference
	double fltlen;
	int iambig = 0; //R = 1/L =-1 coeffient of r_drift from ambigity decision...
        double            dz_max(1.e12) ; // closest_z(1.e12);
        TrkStrawHit *closest(NULL);
        TrkStrawHit  *tsh;
	double doca =0;
        
	for (size_t index=0;index< trackData._chcol->size();++index) {

		ComboHit const& sh = trackData._chcol->at(index);
		//TrkStrawHit = _chcol.;
		//start with doca = 0 
		double 	    poca =0; 
       
		TrkHitVector hitvector = (sh.pos().x(), sh.pos().y(), sh.pos().z());
		//Drift Time :
		//dt        = trackData._chcol->at(index).time()-t0;
		//Get the Straw ID drom SH:
		Straw const&      straw = _tracker->getStraw(sh.strawId());
		//Get Straw "Mid point":
		CLHEP::Hep3Vector hpos  = straw.getMidPoint();

		//Get Straw Direction:
		CLHEP::Hep3Vector hdir  = straw.getDirection();

		//Find DOCA by finding TrkStrawHits with dz closest to the wire center:
		//z of hit = straw mispoint:
                double            zhit = hpos.z();
		//loop over TRKSTRAWHITS:
		TrkStrawHitVector tshv;
		convert(hitvector,tshv);
		vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
	        if(ifnd == tshv.end()){
			
		       Straw const&  trk_straw = _tracker->getStraw(tsh->comboHit().strawId());
                       double        ztrk      = trk_straw.getMidPoint().z();

	               double dz  = ztrk-zhit;
			//Find Point of Closest Approach in z:
		  	if (fabs(dz) < fabs(dz_max)) {
		    		closest   = tsh;
		    		dz_max    = dz;
	  		}
		} 
		poca   = dz_max; //wpoca.doca();
                //Get Track Direction:
		CLHEP::Hep3Vector tdir(trackData._tseed._track.get_track_direction().x(),trackData._tseed._track.get_track_direction().y(), trackData._tseed._track.get_track_direction().z());
		// Estimate flightlength of straw hit along track:
		fltlen = (hpos.z()-tpos.z())/tdir.z();
		doca = poca; //TODO find how the line goesthrough poca 
		if (doca > 0) iambig =  1; //correct by + r_drift onto distance
	        else    iambig = -1; //correct by -r_crift frm distance
        
	}//end chcol	

        //TrkPoca wpoca(trajectory,fltlen,wire,0.0);	 
        //TrkPoca     hitpoca(trajectory,fltlen,wire,0.0);
  	*/





}//end namespace
