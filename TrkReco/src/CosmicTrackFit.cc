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
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minresid(pset.get<unsigned>("minres",1)),
    _minnsh(pset.get<unsigned>("minNStrawHits",4)),
    _minCHHits(pset.get<unsigned>("minCHHits",2)),
    _n_outliers(pset.get<float>("_n_outliers",0.5)),
    _maxresid(pset.get<float>("maxresid",200)),
    _maxniter(pset.get<unsigned>("maxniter",100)),
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

     
      double track_length = sqrt((chN->pos().x()- ch0->pos().x())*(chN->pos().x()- ch0->pos().x()) + (chN->pos().y()- ch0->pos().y())*(chN->pos().y()- ch0->pos().y()) +(chN->pos().z()- ch0->pos().z())*(chN->pos().z()- ch0->pos().z()));

      XYZVec track_dir((chN->pos().x() - ch0->pos().x())/track_length,(chN->pos().y() - ch0->pos().y())/track_length, (chN->pos().z() - ch0->pos().z())/track_length);
      
      //XYVec Pos_Start = (XYVec(FirstP1->pos().x(),FirstP1->pos().y()));
      //XYVec Pos_End = (XYVec(LastP1->pos().x(),LastP1->pos().y()));
      //line->set_m_0(( Pos_End.x() - Pos_Start.x()) / (Pos_End.y() - Pos_Start.y() ));
      
      //line->set_c_0(Pos_Start.x() - ( Pos_Start.y() * line->get_m_0()) );
      
      return track_dir;
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
   //TODO Optimization Step here
    //if track is "good" add to list of good tracks:
   if (goodTrack(all_hits_track)) TrackData._tseed._status.merge(TrkFitFlag::StraightTrackOK);  

}


/*-----------------------Refine Fit ------ -----------------------//
//       Refines the fit in and updates chi2 information      //
//---------------------------------------------------------------*/
CosmicTrack* CosmicTrackFit::FitAll(CosmicTrackFinderData& trackData,  CosmicTrack* cosmictrack, int WeightMode){
    ::BuildMatrixSums S;
    //double resid;
    size_t nHits (trackData._chHitsToProcess.size());
    CosmicTrack* track = &trackData._tseed._track; 
    
    const ComboHit* ch0 = &trackData._chHitsToProcess[0]; 
    const ComboHit* chN = &trackData._chHitsToProcess[nHits];
    XYZVec FirstPoint(ch0->pos().x(),ch0->pos().y(),ch0->pos().z());
    XYZVec LastPoint(chN->pos().x(),chN->pos().y(),chN->pos().z());

    XYZVec track_dir = InitLineDirection(ch0, chN, track);
 
    int        nXYSh(0);
    ComboHit*     hitP1(0);
    
    //loop over hits
    for (size_t f1=1; f1<nHits-1; ++f1){
      hitP1 = &trackData._chHitsToProcess[f1];
      if (!use(*hitP1) )    continue;

      XYZVec point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      XYZVec major_axis =  ParametricFit::MajorAxis(hitP1, track_dir);
      XYZVec minor_axis =  ParametricFit::MinorAxis(hitP1, track_dir);
      double errX =  ParametricFit::HitErrorX(hitP1, major_axis, minor_axis, track_dir);
      double errY =  ParametricFit::HitErrorY(hitP1, major_axis, minor_axis, track_dir);
      S.addPoint(f1, point, track_dir, errX, errY);
      TMatrix P = S.GetAlphaX();
      std::cout<<P[0][0]<<"  "<<P[1][0]<<std::endl;
      // 1) Get X error (DONE)
      // 2) Fill X sums (DONE)
      // 3) Get Y error (DONE)
      // 4) Fill Y sums (DONE)
      // 5) Iterate Errors
      //....... Minimize chi2
      // 6) Set Cosmic Track Parameters as output
 
       //float      minResid(_minresid);
      //if (WeightMode == 0 ) minResid = _maxd;
      //if(resid > minResid){ 
       //	    hitP1->_flag.merge(StrawHitFlag::outlier);  
       //  }
	
	// hitP1->_flag.clear(StrawHitFlag::outlier);
	
	nXYSh += hitP1->nStrawHits();
      }
    return track;
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
