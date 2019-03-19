// Object to perform track fit to combo hits
//
// $Id: CosmicTrackFit code
// $Author: S Middleton
// $Date: Feb 2019
//
// Mu2e
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/DriftCircle.hh"

//Least Squares Fitter:
#include "Mu2eUtilities/inc/LeastSquaresFitter.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/BuildMatrixSums.hh"
//For Drift:
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PanelStateIterator.hh"

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
    _minnsh(pset.get<unsigned>("minNStrawHits",4)),
    _minCHHits(pset.get<unsigned>("minCHHits",2)),
    _n_outliers(pset.get<float>("_n_outliers",0.5)),
    _maxresid(pset.get<float>("maxresid",500)),
    _maxniter(pset.get<unsigned>("maxniter",10)),
    //_minzsep(pset.get<float>("minzsep",stationwidth)),//at least 2 planes
    //_maxzsep(pset.get<float>("maxzsep",???)), //resolutiom limiter..(?)
    _maxd(pset.get<float>("maxd",1000.0)),//max distance between hits at start of fit
    _maxDOCA(pset.get<float>("maxDOCA",1000)),//max distance of closest approach between a hit included in fit and one flag out-right as an outlier
    _maxchi2(pset.get<float>("maxchi2",100.0))   
    {}

//destructor
    CosmicTrackFit::~CosmicTrackFit(){}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(CosmicTrackFinderData& TrackData) {
    if(_debug>0){
    	std::cout<<"Initializing ST Fit ..."<<std::endl;
    }
    bool is_ok(false);
    RunFitChi2(TrackData);
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

 //--------------Fit-----------------//
//Top call to Fitting routines....
//-------------------------------------------// 

  void CosmicTrackFit::BeginFit(CosmicTrackFinderData& TrackData){
      if(_debug>0){
      	std::cout<<" Beginning Fit ..." << std::endl;
      }
      //Clear Previous Flags:
      TrackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
      
    // Initialize:
    bool init(false);
    if (!TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackInit)) {
      init = true;
      if (initCosmicTrack(TrackData))
	TrackData._tseed._status.merge(TrkFitFlag::StraightTrackInit);
      else
	return;
    } 
    //Start Chi2 Fitting:
    if (!init)RunFitChi2(TrackData);
  }


/*------------------------------Chi 2 Fit----------------------------//
//   Adds Chi-2 optimization to fitting routine //
//   Refits and adjusts track fit paramters by weights//
//------------------------------------------------------------------*/
void CosmicTrackFit::RunFitChi2(CosmicTrackFinderData& TrackData) {   
   CosmicTrack* track = &TrackData._tseed._track; 
   //first perform the chi2 fit assuming all hits have same weight
   CosmicTrack* all_hits_track = FitAll(TrackData, track, 0);
   
    //if track is "good" add to list of good tracks:
   if (goodTrack(all_hits_track)) TrackData._tseed._status.merge(TrkFitFlag::StraightTrackOK);  

}

TH1F* hit_errors = new TH1F("Hit Total Errors","Hit Total Errors", 100, 0,100);
TH1F* updated_hit_errors = new TH1F("Hit Total Errors","Hit Total Errors", 100, 0,100);
TH1F* resX = new TH1F("resX","resX", 100,-300,300);
TH1F* resY =   new TH1F("resY","resY", 100,-300,300);
TH1F* updated_resX = new TH1F("newresX","newresX", 100,-1000,1000);
TH1F* updated_resY =   new TH1F("newresY","newresY", 100,-1000,1000);
TH1F* pullX = new TH1F("pullX","pullX", 100,-1,1);
TH1F* pullY =   new TH1F("pullY","pullY", 100,-1,1);
TH1F* chi2 = new TH1F("chi2", "chi2", 100, 0, 2000);
TH1F* updated_chi2 = new TH1F("chi2", "chi2", 100, 0, 2000);

/*---------------Refine Fit ------ ----------------//
//    Refines the fit in and updates chi2 information   //
//-----------------------------------------------*/
CosmicTrack* CosmicTrackFit::FitAll(CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, int WeightMode){
    ::BuildMatrixSums S;
   
    //double resid = 10;
    size_t nHits (trackData._chHitsToProcess.size());
    int DOF = nHits - _Npara;
    //CosmicTrack* track = &trackData._tseed._track; 
    
    const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
    const ComboHit* chN = &trackData._chHitsToProcess[nHits+1];
    XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
    XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());
    std::cout<<" New Track..."<<std::endl;
    //Get Track Basis:
    XYZVec ZPrime = InitLineDirection(ch0, chN, cosmictrack); //Z'=track direction
    std::cout<<"Initial Track Direction: "<<ZPrime.x()<<" "<<ZPrime.y()<<" "<<ZPrime.z()<<std::endl;
    XYZVec XPrime = ParametricFit::GetXPrime(ZPrime);
    XYZVec YPrime = ParametricFit::GetYPrime(XPrime, ZPrime);
    XYZVec XDoublePrime = ParametricFit::GetXDoublePrime(XPrime, YPrime, ZPrime);
    XYZVec YDoublePrime = ParametricFit::GetYDoublePrime(XPrime, YPrime, ZPrime);
    std::cout<<"X'': "<<XDoublePrime.x()<<" "<<XDoublePrime.y()<<" "<<XDoublePrime.z()<<std::endl;
    std::cout<<"Y'': "<<YDoublePrime.x()<<" "<<YDoublePrime.y()<<" "<<YDoublePrime.z()<<std::endl;
    
    ParametricFit::TestConditions( XDoublePrime, YDoublePrime, ZPrime);
    cosmictrack->setXPrime(XDoublePrime);
    cosmictrack->setYPrime(YDoublePrime);
    cosmictrack->setZPrime(ZPrime);
    
    cosmictrack->set_track_direction(ZPrime);
    int        nXYSh(0);
    ComboHit     *hitP1(0), *hitP2(0); 
    HitInfo_t  indexBestComboHit;
    //float      minResid(_minresid);
    //if (WeightMode == 0 ) minResid = _maxd;
    
    for (size_t f1=0; f1<nHits+1; ++f1){  
      hitP1 = &trackData._chHitsToProcess[f1];
      if (!use(*hitP1) )    continue;
      
      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      XYZVec major_axis =  ParametricFit::MajorAxis(hitP1);
      XYZVec minor_axis =  ParametricFit::MinorAxis(hitP1);
      double errX =  ParametricFit::HitErrorX(hitP1, major_axis, minor_axis, XDoublePrime);
      double errY =  ParametricFit::HitErrorY(hitP1, major_axis, minor_axis, YDoublePrime);
      nXYSh += hitP1->nStrawHits();
      
      S.addPoint(f1, point, XDoublePrime, YDoublePrime, errX, errY);
    }//end intial hit loop
     cosmictrack->set_chisq(S.GetTotalChi2()/DOF);
     chi2->Fill(S.GetTotalChi2()/DOF);
     
     double a0 = S.GetAlphaX()[0][0];
     double a1 = S.GetAlphaX()[1][0];
     double b0 = S.GetAlphaY()[0][0];
     double b1 = S.GetAlphaY()[1][0];
     
     cosmictrack->set_parameters(a0,a1,b0,b1);
     XYZVec updated_track_direction(a1, b1, 1); 
     cosmictrack->set_track_direction(updated_track_direction.Unit());
      std::cout<<"First Fit Track Direction: "<<updated_track_direction.Unit().x()<<" "<<updated_track_direction.Unit().y()<<" "<<updated_track_direction.Unit().z()<<std::endl;
     for (size_t f2=0; f2<nHits+1; ++f2){
     	     
     	      hitP2 = &trackData._chHitsToProcess[f2];
      	      if (!use(*hitP2) && nXYSh < _minnsh )    continue;     	   
	      XYZVec point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
	     
	      XYZVec major_axis =  ParametricFit::MajorAxis(hitP1);
      	      XYZVec minor_axis =  ParametricFit::MinorAxis(hitP1);
      	      double errX =  ParametricFit::HitErrorX(hitP1, major_axis, minor_axis, XDoublePrime);
      	      double errY =  ParametricFit::HitErrorY(hitP1, major_axis, minor_axis, YDoublePrime);
	      double hit_error =  ParametricFit::TotalHitError(hitP2, major_axis, minor_axis, XDoublePrime, YDoublePrime);
	      
	     
	      double Rx = ParametricFit::GetResidualX(a0,a1, XDoublePrime, point);
	      double Ry = ParametricFit::GetResidualY(b0, b1, YDoublePrime, point);
	      
	      //Fill diagnostics (best to use analyzer but for debugging leave in for now)
	      hit_errors->Fill(hit_error);
	      resX->Fill(Rx);
	      resY->Fill(Ry);
	      pullX->Fill(Rx/errX);
	      pullY->Fill(Ry/errY);
	      cosmictrack->set_hit_errorsTotal(hit_error);
	      cosmictrack->set_fit_residual_errorsX(errX);
	      cosmictrack->set_fit_residual_errorsY(errY);
	      cosmictrack->set_fit_residualsX(Rx);
	      cosmictrack->set_fit_residualsY(Ry);
//--------------------------------------------------------------//	
                 
	     //iterate errors:
	     unsigned ninter(0);
	     bool errors_converged = false;
	    
	     while(ninter < _maxniter && errors_converged==false){
		ninter +=1;
		//Define a place to store error/track information from previous iteration:
    		XYZVec previous_XDoublePrime = XDoublePrime;
	    	XYZVec previous_YDoublePrime = YDoublePrime;
	    	//XYZVec previous_ZPrime = ZPrime;
	    	
		//Store errors from previous iteration:
		double previous_hit_errorX = errX;
    		double previous_hit_errorY = errY;
	    		
		//Get Track coordinate system informations:
		ZPrime = cosmictrack->get_track_direction();
		XPrime = ParametricFit::GetXPrime(ZPrime);
    		YPrime = ParametricFit::GetYPrime(XPrime, ZPrime);
    		XDoublePrime = ParametricFit::GetXDoublePrime(XPrime, YPrime, ZPrime);
    		YDoublePrime = ParametricFit::GetYDoublePrime(XPrime, YPrime, ZPrime);
    		
    		//Get New Errors:
    		errX = ParametricFit::HitErrorX(hitP1, major_axis, minor_axis, XDoublePrime);
                errY = ParametricFit::HitErrorY(hitP1, major_axis, minor_axis, YDoublePrime);
    		
		//Define error change
		double d_errorX =0;
		double d_errorY =0;
		
		//If new error in X is worse than old X error but Y is better update fit:		
		if (abs(errX) > abs(previous_hit_errorX) && abs(errY) < abs(previous_hit_errorY)){
			d_errorX = sqrt(pow(abs(errX) - abs(previous_hit_errorX),2));
		        S.removePoint(f2, point, previous_XDoublePrime, previous_YDoublePrime, previous_hit_errorX, previous_hit_errorY);
                	S.addPoint(f2, point, XDoublePrime, YDoublePrime, errX, previous_hit_errorY);	
		}
		
		//If new error in Y is worse than old X error but X is better update fit:
		if (abs(errY) > abs(previous_hit_errorY) && abs(errX) < abs(previous_hit_errorX)){
			d_errorY = sqrt(pow((abs(errY) - abs(errY)),2));
		        S.removePoint(f2, point, previous_XDoublePrime, previous_YDoublePrime, previous_hit_errorX, previous_hit_errorY);
                	S.addPoint(f2, point, XDoublePrime, YDoublePrime, previous_hit_errorX, errY );	
		}
		//If both worse:
		if(abs(errX) > abs(previous_hit_errorX) && abs(errY) > abs(previous_hit_errorY)){
			d_errorX = sqrt(pow(abs(errX) - abs(previous_hit_errorX),2));
			d_errorY = sqrt(pow((abs(errY) - abs(previous_hit_errorY)),2));
			S.removePoint(f2, point, previous_XDoublePrime, previous_YDoublePrime, previous_hit_errorX, previous_hit_errorY);
			S.addPoint(f2, point, XDoublePrime, YDoublePrime, errX, errY );
		}
		//If the change is small stop iterations:
        	if( (d_errorX) < 0.1 && (d_errorY) < 0.1 ){ 
           		errors_converged = true;
           		
                } 
                
                //update the track fit information before next iteration:
                XYZVec new_track_direction(S.GetAlphaX()[1][0],S.GetAlphaY()[1][0],1);//Relative to Z'
	     	cosmictrack->set_track_direction(new_track_direction.Unit());
	     }//end while
	     
//--------------------------------------------------------------//	 
	     	//Get updated Parameters:
	     	double a0 = S.GetAlphaX()[0][0];
     		double a1 = S.GetAlphaX()[1][0];
     		double b0 = S.GetAlphaY()[0][0];
     		double b1 = S.GetAlphaY()[1][0];
     	        //Set final track details:
	     	XYZVec final_ZPrime = cosmictrack->get_track_direction();
		XYZVec final_XPrime = ParametricFit::GetXPrime(final_ZPrime);
    		XYZVec final_YPrime = ParametricFit::GetYPrime(final_XPrime, final_ZPrime);
    		XYZVec final_XDoublePrime = ParametricFit::GetXDoublePrime(final_XPrime, final_YPrime, final_ZPrime);
    		XYZVec final_YDoublePrime = ParametricFit::GetYDoublePrime(final_XPrime, final_YPrime, final_ZPrime);
	     	//set final parameters here
	        double newRx = ParametricFit::GetResidualX(a0,a1, final_XDoublePrime, point);
	        double newRy = ParametricFit::GetResidualY(b0, b1, final_YDoublePrime, point);
	        
	        //Cut on residuals/flag outliers:
	        if(newRx > _maxresid || newRy > _maxresid){ 
		  hitP2->_flag.merge(StrawHitFlag::outlier); 
		  S.removePoint(f2, point, final_XDoublePrime, final_YDoublePrime, errX, errY);    	   
		 }//end resid
		hit_error =  ParametricFit::TotalHitError(hitP2, major_axis, minor_axis, final_XDoublePrime, final_YDoublePrime);    	   
		updated_hit_errors->Fill(hit_error);
	        updated_resX->Fill(newRx);
	        updated_resY->Fill(newRy);
	        
    }//end second loop
    
    updated_chi2->Fill(S.GetTotalChi2()/DOF);
    updated_hit_errors->SaveAs("Final Errors");
    hit_errors->SaveAs("Combined Hit Error");
    resX->SaveAs("oldresX");
    resY->SaveAs("oldresY");
    updated_resX->SaveAs("newresX");
    updated_resY->SaveAs("newresY");
    pullX->SaveAs("oldpullX");
    pullY->SaveAs("oldpullY");
    chi2->SaveAs("oldchi2");
    updated_chi2->SaveAs("newchi2");
    S.clear();
    return cosmictrack;
  }




void MulitpleTrackResolver(CosmicTrackFinderData& trackData,CosmicTrack* track){
	//if track has a significant amount of tracks with large chi2 and those large values are all similar...and resiudals within range of each other, this coule mean a second track....check for this although unlikely with cosmics....

}


void CosmicTrackFit::UpdateFitErrors( std::vector<XYZVec> HitVec, std::vector<double> err, CosmicTrack* track, std::vector<XYZVec> major_axis,std::vector<XYZVec> minor_axis){

    for(int i=0; i< static_cast<int>(HitVec.size()); i++){
	bool errors_converged = false;
        
	while(errors_converged==false){

		//Get Current Track informations:
		XYZVec updated_track_dir = track->get_track_direction();

                //Getting hit errors:
                double hit_error = sqrt(major_axis[i].Dot(updated_track_dir)*major_axis[i].Dot(updated_track_dir)+minor_axis[i].Dot(updated_track_dir)*minor_axis[i].Dot(updated_track_dir));
		
		//Update Error:
		double d_error =0;
		if (hit_error > err[i]){
			d_error = sqrt((hit_error - err[i])*(hit_error - err[i]));
		        err[i] = hit_error; 	
		}
        	if( d_error < 0.1){ 
           		errors_converged = true;
                }         
       }//end while 
    }//end for
 
}//end update



bool CosmicTrackFit::goodTrack(CosmicTrack* track)
  { 
    
    
    if(track->get_chisq() < _maxchi2) return true;
    else return false;
    
  }

/*--------------USE HIT?---------------------------//
//          Checks if flag                        //
//------------------------------------------------*/
bool CosmicTrackFit::use(const ComboHit& thit) const 
  {
    return (!thit._flag.hasAnyProperty(_dontuseflag));
  }

float CosmicTrackFit::pointToLineDCA(CosmicTrack* track, StrawHit hit) {
        /*
	//here we consider only in 2D (YZ plane) i.e. this is the z direction we are correcting for not the XY plane
	float z_slope = -1*track->get_m_1();//zslope
	float x_slope = -1*track->get_m_0();
	float intercept = -1*track->get_c_0();
	float z_straw = 1.;//zStraw;
	float y_straw = 1.;//yStraw;
        float dca =1.0;
	//float dca = abs( z_slope * z_straw + x_slope * x_Straw + intercept ) / sqrt( a * a + b * b  ) ;
	
	return dca;
	*/
	return 1.0;
}

 /*--------Drift Correction-------*/
void CosmicTrackFit::DriftCorrection(CosmicTrackFinderData& trackData){
	//Find SH in CH with smallest distance from wire centre, then need to find relative sign of that SH doca value. Then we know which side of wire the track approached. Then add on/minus off to z this distance
        /*
	TrkT0 t0 = trackData._tseed._t0 ; //Get t0 from tclust
	//double            dz_max(1.e12) ; // closest_z(1.e12);

        HepPoint   tpos; //position (c?)
        //double     dt; //time difference

	int iambig = 0; //R = 1/L =-1 coeffient of r_drift from ambigity decision...
        double            dz_max(1.e12) ; // closest_z(1.e12);
        TrkStrawHit *closest(NULL);
        TrkStrawHit  *tsh;
	double doca =0;
        
	for (size_t index=0;index< trackData._chcol->size();++index) {

		ComboHit const& sh = trackData._chcol->at(index);
		TrkStrawHit = _chcol.;
		//start with doca = 0 
		double 	    poca =0; 
       
		//TrkHitVector hitvector = (sh.pos().x(), sh.pos().y(), sh.pos().z());
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
		convert(_chcol,tshv);
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
		CLHEP::Hep3Vector tdir(trackData._tseed._track.get_m_0(), 1, trackData._tseed._track.get_m_1());
		// Estimate flightlength of straw hit along track:
		//fltlen = (hpos.z()-tpos.z())/tdir.z();
		doca = poca; //TODO find how the line goesthrough poca
                
		if (doca > 0) iambig =  1; //correct by + r_drift onto distance
	        else    iambig = -1; //correct by -r_crift frm distance
               
                
	}//end chcol	



                //TrkPoca wpoca(trajectory,fltlen,wire,0.0);
                

		 
                //TrkPoca     hitpoca(trajectory,fltlen,wire,0.0);
	*/  	
}




}//end namespace
