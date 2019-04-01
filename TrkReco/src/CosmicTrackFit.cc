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
//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
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
    _dim(pset.get<int>("dimension",3)),
    _Npara(pset.get<int>("Npara",4)),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minresid(pset.get<unsigned>("minres",0)),
    _minnsh(pset.get<unsigned>("minNStrawHits",2)),
    _minCHHits(pset.get<unsigned>("minCHHits",4)),
    _n_outliers(pset.get<unsigned>("_n_outliers",5)),
    _D_error(pset.get<float>("Derror",0.0001)),
    _maxresid(pset.get<float>("maxresid",500)),
    _maxniter(pset.get<unsigned>("maxniter",1000)),
    _maxpull(pset.get<double>("maxPull",5)),
    //_minzsep(pset.get<float>("minzsep",stationwidth)),//at least 2 planes
    //_maxzsep(pset.get<float>("maxzsep",???)), //resolutiom limiter..(?)
    _maxd(pset.get<float>("maxd",300.0)),//max distance between hits at start of fit
    _maxDOCA(pset.get<float>("maxDOCA",1000)),//max distance of closest approach between a hit included in fit and one flag out-right as an outlier
    _maxchi2(pset.get<float>("maxchi2",1000.0))   
    {}

//destructor
    CosmicTrackFit::~CosmicTrackFit(){}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics) {
    if(_debug>0){
    	std::cout<<"Initializing ST Fit ..."<<std::endl;
    }
    bool is_ok(false);
    RunFitChi2(TrackData,diagnostics );
    is_ok = TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackOK);
    return is_ok;
  }


  /*-------------Init Line-------------------------//
 Makes basic estimate of a line y=mx+c style with Cosmic line between first and last points on plane
  //----------------------------------------------*/
  XYZVec CosmicTrackFit::InitLineDirection(const ComboHit *ch0, const ComboHit *chN,CosmicTrack* line) {
     
      double tx = chN->pos().x() - ch0->pos().x();
      double ty = chN->pos().y() - ch0->pos().y();
      double tz = chN->pos().z() - ch0->pos().z();
      
      XYZVec track(tx,ty,tz);
      return track.Unit();
    } 
    
    XYZVec CosmicTrackFit::LineDirection(double a1, double b1, const ComboHit *ch0, const ComboHit *chN) {
     
      double tx = a1*(chN->pos().z() - ch0->pos().z());
      double ty = b1*(chN->pos().z() - ch0->pos().z());
      double tz = chN->pos().z() - ch0->pos().z();
      
      XYZVec track(tx,ty,tz);
      return track.Unit();
    } 
//--------------Fit-----------------//
//Top call to Fitting routines....
//-------------------------------------------// 
  void CosmicTrackFit::BeginFit(CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics){
      if(_debug>0){
      	std::cout<<" Beginning Fit ..." << std::endl;
      }
      //Clear Previous Flags:
      TrackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
      
    // Initialize:
    bool init(false);
    if (!TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackInit)) {
      init = true;
      if (initCosmicTrack(TrackData, diagnostics))
	TrackData._tseed._status.merge(TrkFitFlag::StraightTrackInit);
      else
	return;
    } 
    //Start Chi2 Fitting:
    if (!init)RunFitChi2(TrackData, diagnostics);
  }
/*------------------------------Chi 2 Fit----------------------------//
//   Adds Chi-2 optimization to fitting routine //
//   Refits and adjusts track fit paramters by weights//
//------------------------------------------------------------------*/
void CosmicTrackFit::RunFitChi2(CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics) {   
   CosmicTrack* track = &TrackData._tseed._track; 
   //first perform the chi2 fit assuming all hits have same weight
   FitAll(TrackData, track, diagnostics);
    //if track is "good" add to list of good tracks:
    //if (goodTrack(all_hits_track)) 
   TrackData._tseed._status.merge(TrkFitFlag::StraightTrackOK);  
   //tag track as "converged" if it reaches this point in routine
   
   TrackData._diag.CosmicTrackFitCounter += 1;//label as having a track
   diagnostics.npasses +=1;//totsl number of tracks
   std::cout<<"passes all cuts "<<diagnostics.npasses<<std::endl;
}



/*---------------Refine Fit ------ ----------------//
//    Refines the fit in and updates chi2 information   //
//-----------------------------------------------*/
void CosmicTrackFit::FitAll(CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, CosmicTrackFinderTypes::Data_t& diagnostics){
    std::cout<<" -------------------------------"<<std::endl;
    ::BuildMatrixSums S, S_init;
    
    unsigned        n_outliers(0), nhits_passed(0);
    ComboHit   *hitP1(0), *hitP2(0); 
    size_t nHits (trackData._chHitsToProcess.size());
    //Get initial seed:
    const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
    const ComboHit* chN = &trackData._chHitsToProcess[nHits+1];
    XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
    XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());
    double track_length = sqrt(pow(chN->pos().x()-ch0->pos().x(),2)+pow(chN->pos().y()-ch0->pos().y(),2+pow(chN->pos().z()-ch0->pos().z(),2)));
    
    //Get Track Basis:
    XYZVec ZPrime = InitLineDirection(ch0, chN, cosmictrack); //Z'=track direction
    cosmictrack->set_initial_track_direction(ZPrime);
    cosmictrack->set_track_direction(ZPrime);
    XYZVec XPrime = ParametricFit::GetXPrime(ZPrime);
    XYZVec YPrime = ParametricFit::GetYPrime(XPrime, ZPrime);   
    XYZVec XDoublePrime = ParametricFit::GetXDoublePrime(XPrime, YPrime, ZPrime);
    XYZVec YDoublePrime = ParametricFit::GetYDoublePrime(XPrime, YPrime, ZPrime);
    
    //First loop to update track information and hit error estimate
    for (size_t f1=0; f1<nHits+1; ++f1){  
      
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
     std::cout<<"passed first cut "<<nhits_passed<<std::endl;
     //Get first estimate of track parameters
     int DOF = (nHits - _Npara);
  
     double a0 = S_init.GetAlphaX()[0][0];
     double a1 = S_init.GetAlphaX()[1][0];
     double b0 = S_init.GetAlphaY()[0][0];
     double b1 = S_init.GetAlphaY()[1][0];
     
     cosmictrack->set_initchisq_dof(S_init.GetTotalChi2()/abs(DOF)); 
     std::cout<<"Initial Parameters"<<a0<<" "<<a1<<" "<<b0<<" "<<b1<<std::endl;
     cosmictrack->set_initial_parameters(a0,a1,b0,b1);
     cosmictrack->set_parameters(a0,a1,b0,b1);
     XYZVec updated_track_direction(a1, b1, 1);  
     cosmictrack->set_track_direction(updated_track_direction.Unit());
     
     XYZVec init_XDoublePrime = XDoublePrime;
     XYZVec init_YDoublePrime = YDoublePrime;
     XYZVec init_ZPrime = ZPrime;
     //second loop to update errors and iterate:
     for (size_t f2=0; f2<nHits+1; ++f2){
     	      if(isnan(cosmictrack->get_track_direction().Mag2()) == true) continue;     
     	      hitP2 = &trackData._chHitsToProcess[f2];
      	      if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;
      	      // Store axes for later:
              
	      //std::cout<<"Initial ZPrime"<<ZPrime<<std::endl;
      	      //Get error axes again: 	   
	      XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
	      XYZVec major_axis =  ParametricFit::MajorAxis(hitP2);
      	      XYZVec minor_axis =  ParametricFit::MinorAxis(hitP2);
      	      double errX =  ParametricFit::HitErrorX(hitP2, major_axis, minor_axis, XDoublePrime);
      	      double errY =  ParametricFit::HitErrorY(hitP2, major_axis, minor_axis, YDoublePrime);
	      double hit_error =  ParametricFit::TotalHitError(hitP2, major_axis, minor_axis, XDoublePrime,YDoublePrime);
	      //Get initial residuals and store for analysis later:
	      
	      double Rx = ParametricFit::GetResidualX(a0,a1, init_XDoublePrime, point);
	      double Ry = ParametricFit::GetResidualY(b0, b1, init_YDoublePrime, point);
	       
	      cosmictrack->set_init_hit_errorsTotal(hit_error);
	      cosmictrack->set_init_fit_residual_errorsX(errX);
	      cosmictrack->set_init_fit_residual_errorsY(errY);
	      cosmictrack->set_init_fit_residualsX(Rx);
	      cosmictrack->set_init_fit_residualsY(Ry);
	      
	     //iterate errors:
	     unsigned niter(0);
	     bool errors_converged = false;
	     while(niter < _maxniter && errors_converged==false){
	        	        
		niter +=1;
		std::cout<<"++++++++++++++++++++++"<<std::endl;
		std::cout<<"Hit : "<<f2<<" Iteration: "<<niter<<std::endl;
		
		//To begin remove point which is being evaluated:
		if (niter == 1){
			S.removePoint(f2, point, init_XDoublePrime, init_YDoublePrime, init_ZPrime, errX, errY );
		}
		else{
			S.removePoint(f2, point, XDoublePrime, YDoublePrime, ZPrime, errX, errY );
		}
		//Define a place to store error/track information from previous iteration:
    		XYZVec previous_XDoublePrime = XDoublePrime;
	    	XYZVec previous_YDoublePrime = YDoublePrime;
	    	
		//Store errors from previous iteration:
		double previous_hit_errorX = errX;
    		double previous_hit_errorY = errY;
	    		
		//Get Track coordinate system informations: (first iteration uses initial guess of ZPrime)
		XPrime = ParametricFit::GetXPrime(ZPrime);
    		YPrime = ParametricFit::GetYPrime(XPrime, ZPrime);	
    		XDoublePrime = ParametricFit::GetXDoublePrime(XPrime, YPrime, ZPrime);
    		YDoublePrime = ParametricFit::GetYDoublePrime(XPrime, YPrime, ZPrime);
    		
    		//Get error projections on current axes:
    		errX = ParametricFit::HitErrorX(hitP2, major_axis, minor_axis, XDoublePrime);
                errY = ParametricFit::HitErrorY(hitP2, major_axis, minor_axis, YDoublePrime);
    		if(isnan(errX) == true || isnan(errY)==true) continue;	
		//Define error change
		double d_errorX = sqrt(pow(abs(errX) - abs(previous_hit_errorX),2));
		double d_errorY = sqrt(pow((abs(errY) - abs(previous_hit_errorY)),2));
		std::cout<<" errors X "<<errX<<" previous "<<previous_hit_errorX<<" errors Y "<<errY<<" previous "<<previous_hit_errorY<<std::endl;
		
		//If new error in X is worse than old X error but Y is better update fit as:		
		if (abs(errX) > abs(previous_hit_errorX) && abs(errY) <= abs(previous_hit_errorY)){
			//d_errorX = sqrt(pow(abs(errX) - abs(previous_hit_errorX),2));
			d_errorY = 0;
		        YDoublePrime = previous_YDoublePrime;
		        errY = previous_hit_errorY;
                	S.addPoint(f2, point, XDoublePrime, YDoublePrime, ZPrime, errX, errY);
                	std::cout<<" errX worse "<<std::endl;	
		}
		
		//If new error in Y is worse than old X error but X is better update fit as:
		if (abs(errY) > abs(previous_hit_errorY) && abs(errX) <= abs(previous_hit_errorX)){
			//d_errorY = sqrt(pow((abs(errY) - abs(previous_hit_errorY)),2));
			d_errorX = 0;
		        XDoublePrime = previous_XDoublePrime;
		        errX = previous_hit_errorX;
                	S.addPoint(f2, point, XDoublePrime, YDoublePrime,ZPrime, errX, errY );	
                	std::cout<<" errY worse "<<std::endl;
		}	
		
		//If both worse then update fit as:
		if(abs(errX) > abs(previous_hit_errorX) && abs(errY) > abs(previous_hit_errorY)){
		        //d_errorX = sqrt(pow(abs(errX) - abs(previous_hit_errorX),2));
			//d_errorY = sqrt(pow((abs(errY) - abs(previous_hit_errorY)),2));
			S.addPoint(f2, point, XDoublePrime, YDoublePrime, ZPrime, errX, errY );
			std::cout<<" both worse "<<std::endl;
		}
		//If neither is better add back in the old value and return to previous track information
		if(abs(errX) <= abs(previous_hit_errorX) && abs(errY) <= abs(previous_hit_errorY)){
			d_errorX = 0;
			d_errorY = 0;
			errX = previous_hit_errorX;
			errY = previous_hit_errorY;
			XDoublePrime = previous_XDoublePrime;
			YDoublePrime = previous_YDoublePrime;
			S.addPoint(f2, point, XDoublePrime, YDoublePrime, ZPrime, errX, errY);
			std::cout<<" errY and X same "<<std::endl;
		}
		
		//If the change is small stop iteration and claim track is converged:
        	if( (d_errorX) < _D_error && (d_errorY) < _D_error ){ 
        		std::cout<<"Converged!!!"<<std::endl;
           		errors_converged = true;
			cosmictrack->set_niter(niter);
                } 
                if( (errors_converged == false && niter == _maxniter)){
                        cosmictrack->add_outlier(hitP2);
                 }
                //if never converges just remove this point but not the track it could be an odd outlier:
                if( (errors_converged == false && niter == _maxniter)){
                       
                	nHits -=1;
                	hitP2->_flag.merge(StrawHitFlag::outlier);
                	n_outliers +=1;
                	continue;
                	
                }    	
                if(n_outliers > _n_outliers) continue;
                //update the track fit information before next iteration:
                XYZVec new_track_direction(S.GetAlphaX()[1][0],S.GetAlphaY()[1][0],1);//Relative to Z' 
                ZPrime = new_track_direction.Unit();                
	     	//ZPrime.SetXYZ(new_track_direction.X(),new_track_direction.Y(), new_track_direction.Z());
	     	cosmictrack->setXPrime(XDoublePrime);
    		cosmictrack->setYPrime(YDoublePrime);
    		cosmictrack->setZPrime(ZPrime);
	       }
	        if(nHits < _minCHHits) continue;
	        
	     	//Update stored track information
	     	double a0 = S.GetAlphaX()[0][0];
     		double a1 = S.GetAlphaX()[1][0];
     		double b0 = S.GetAlphaY()[0][0];
     		double b1 = S.GetAlphaY()[1][0];
     		cosmictrack->clear_parameters();
     		cosmictrack->set_parameters(a0,a1,b0,b1);
     		/*
     		std::cout<<"Final XDoublePrime"<<XDoublePrime<<std::endl;
     		std::cout<<"Final YDoublePrime"<<YDoublePrime<<std::endl;
     		std::cout<<"Final ZPrime"<<ZPrime<<std::endl;
     		std::cout<<"Final para"<<a0<<" "<<a1<<" "<<b0<<" "<<b1<<std::endl;
     		std::cout<<"NOut"<<n_outliers<<std::endl;
     		*/
	        double newRx = ParametricFit::GetResidualX(a0,a1, XDoublePrime,  point);
	        double newRy = ParametricFit::GetResidualY(b0, b1, YDoublePrime, point);
	        hit_error =  ParametricFit::TotalHitError(hitP2, major_axis, minor_axis, XDoublePrime, YDoublePrime);   
	        //cut on residuals of hit:
	        if(newRx > _maxresid || newRy > _maxresid){ 
		  hitP2->_flag.merge(StrawHitFlag::outlier);
		  cosmictrack->add_outlier(hitP2);
		  n_outliers +=1;
		  
		 } 
	        
		//Set Track Details:
		cosmictrack->set_track_length(track_length);
		cosmictrack->set_final_hit_errorsTotal(hit_error);
		cosmictrack->set_final_fit_residual_errorsX(errX);
	      	cosmictrack->set_final_fit_residual_errorsY(errY);
	      	cosmictrack->set_final_fit_residualsX(newRx);
	      	cosmictrack->set_final_fit_residualsY(newRy);
	      	cosmictrack->set_track_direction(ZPrime);
	      	cosmictrack->set_track_position(a0,b0,0); 	
    		cosmictrack->set_N(nHits);
    		DOF = (nHits - _Npara); 
   		
    		//New Diag Framwork (In progress)
    		diagnostics.nhits[diagnostics.nseeds-1] =nHits;
    		diagnostics.chi2d_track[diagnostics.nseeds-1] = S.GetTotalChi2()/abs(DOF);
    		diagnostics.hit_residualX[diagnostics.nseeds-1][f2]=(newRx);
    		diagnostics.hit_residualY[diagnostics.nseeds-1][f2] = (newRy);
    		diagnostics.hit_errorX[diagnostics.nseeds-1][f2]=(errX);
    		diagnostics.hit_errorY[diagnostics.nseeds-1][f2] = (errY);
    		diagnostics.nChFit[diagnostics.nseeds-1]= +1;
    		
    }//end second loop
    cosmictrack->set_finalchisq_dof(S.GetTotalChi2()/abs(DOF));
    cosmictrack->set_finalchisq_dofY(S.GetChi2Y()/abs(DOF));
    cosmictrack->set_finalchisq_dofX(S.GetChi2X()/abs(DOF));
    diagnostics.chi2numbers +=1;
    trackData._tseed._status.merge(TrkFitFlag::StraightTrackConverged);
    std::cout<<"chi2 number"<<diagnostics.chi2numbers<<std::endl;
    //DriftCorrection(trackData);   
    //S.clear();
    //return cosmictrack;
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
	 	//XYZVec trackPOCA =  ParametricFit::PointToLineCA(point, FirstPoint,  LastPoint);
	 	//XYZVec wirePOCA =  ParametricFit::PointToLineCA(pointwire, FirstPoint,  LastPoint);
	        //double trackDOCA = ParametricFit::PointToLineDCA(point, FirstPoint, LastPoint);
	        
	       
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
