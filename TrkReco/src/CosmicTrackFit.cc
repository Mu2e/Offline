// Object to perform track fit to combo hits
//
// $Id: StriaghtTrackFit code
// $Author: S Middleton
// $Date: Nov 2018
//
// Mu2e
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"


//Least Squares Fitter:
#include "Mu2eUtilities/inc/LeastSquaresFitter.hh"

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
    _dim(pset.get<int>("dimension",2)),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minresid(pset.get<unsigned>("minres",0)),
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
  StraightTrack* CosmicTrackFit::InitLine(const ComboHit *FirstP1, const ComboHit *LastP1,StraightTrack* line) {
      
      XYVec Pos_Start = (XYVec(FirstP1->pos().x(),FirstP1->pos().y()));
      XYVec Pos_End = (XYVec(LastP1->pos().x(),LastP1->pos().y()));
      line->set_m_0(( Pos_End.x() - Pos_Start.x()) / (Pos_End.y() - Pos_Start.y() ));
      
      line->set_c_0(Pos_Start.x() - ( Pos_Start.y() * line->get_m_0()) );
      if(_debug>0){
	      
	      std::cout<<"initialized line .."<<line->get_m_0()<<"  "<<line->get_c_0()<<std::endl;
      }
      return line;
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
    //int  nHits(TrackData._chHitsToProcess.size());
    //ComboHit*     firstPt(0), *lastPt(0) ;
    
    CosmicTrack* track = &TrackData._tseed._track; 
    
    //loop over faces
    //firstPt = &TrackData._chHitsToProcess[0];
    //lastPt = &TrackData._chHitsToProcess[nHits];


    //track = InitLine(firstPt, lastPt,track);
    // might use later- set initial value high so it becomes best fit..
    
   //first perform the chi2 fit assuming all hits have same weight
   CosmicTrack* all_hits_track = FitAll(TrackData, track, 0);
   
    //if track is "good" add to list of good tracks:
   if (goodTrack(all_hits_track)) TrackData._tseed._status.merge(TrkFitFlag::StraightTrackOK);  


  
}


/*-----------------Refine Fit in XY----------------//
//Refines the fit in xy and updates chi2 information//
//-------------------------------------------------*/
CosmicTrack* CosmicTrackFit::FitAll(CosmicTrackFinderData& trackData,  CosmicTrack* track, int WeightMode){
 
    CosmicTrack* track = &trackData._tseed._track; 
    
    ComboHit ch0 = &trackData._chHitsToProcess[0]; 
    

    float      resid;
    int        minNReducedChi2Points(15);//TODO: Optimize
    int  nHits(trackData._chHitsToProcess.size());
    int        nXYSh(0);
    ComboHit*     hitP1(0);

    //For Fit:
    std::vector<XYZVec> ma, mi;
  
    //loop over hits
    for (int f1=1; f1<nHits-2; ++f1){
      hitP1 = &trackData._chHitsToProcess[f1];
      if (!use(*hitP1) )    continue;

      XYZ track_dir(hitP1.pos().x() - ch0.pos().x(),hitP1.pos().y() - ch0.pos().y(),hitP1.pos().z() - ch0.pos().z()); 

      double DOCA = ParametricFit::;

      //For Error Ellipses:
      XYZVec const& wdir = hitP1->wdir();//direction along wire
      XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire
      double werr_mag = hitP1->wireRes(); //hit major error axis 
      double terr_mag = hitP1->transRes(); //hit minor error axis
      XYZVec major_axis = werr_mag*wdir;
      XYZVec minor_axis = terr_mag*wtdir;
  
      HitInfo_t  indexBestComboHit;
      
      float      minResid(_minresid);
      if (WeightMode == 0 ) minResid = _maxd;
      
      //Get Initial Track Direction:
      XYZVec track_dir(1,1,1);// CosmicTrack.get_track_direction();
      
      //Get Hit Error:
      double hit_error = sqrt(major_axis.Cross(track_dir).mag2()+minor_axis.Cross(track_dir).mag2());
      
      //DCA
      resid =  1/hit_error;
       

	//if residual of new line is better than previous residual then added to hits:
	if (resid < minResid) { 
          
	  indexBestComboHit.face          = f1;
	  indexBestComboHit.panel         = hitP1->strawId().uniquePanel();
	  
	  //update minResid
	  minResid          = resid;
	  hitP1->_xyWeight  = hit_error;//or 1?
          if(_debug>0){
              
              std::cout<<"residual better F="<<indexBestComboHit.face <<" P = " <<indexBestComboHit.panel<<std::endl;
              
          }
	}
        //Reset fit:
        track->clear();

	if(resid > minResid){ 
            if(_debug>0){
                std::cout<<"Outlier ..."<<std::endl;
            }
	    hitP1->_flag.merge(StrawHitFlag::outlier);  
         }
      
	    
      //now add the best hit to least squares fit if found:
      if (indexBestComboHit.face >=0 ) {
	 
	 hitP1    = &trackData._chHitsToProcess[indexBestComboHit.face];
         
         ma.push_back(major_axis);
         mi.push_back(minor_axis);     
       
	 //remove the outlier flag
	 hitP1->_flag.clear(StrawHitFlag::outlier);

	 //increase the StawaHit counter
	 nXYSh += hitP1->nStrawHits();

      }
	    
    }//end loop over the faces

    //if we collected enough points update the results
    if (nXYSh >= _minnsh){

      //flag as StraightTRackOK
     } 
    //Update Errors
    return track;
  }



void MulitpleTrackResolver(CosmicTrackFinderData& trackData,CosmicTrack* track){
	//if track has a significant amount of tracks with large chi2 and those large values are all similar...and resiudals within range of each other, this coule mean a second track....check for this although unlikely with cosmics....

}


void CosmicTrackFit::UpdateFitErrors(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> err, StraightTrack* track,TMatrixD cov_x, std::vector<XYZVec> major,std::vector<XYZVec> minor){

    for(int i=0; i< static_cast<int>(x.size()); i++){
	bool errors_converged = false;
        //Get Hit Error:
	XYZVec major_axis = major[i];
	XYZVec minor_axis = minor[i];
	if(_diag >0){
        	std::cout<<"Starting error : "<<i<<" "<<err[i]<<std::endl;
        }
	while(errors_converged==false){
	        
		
		//Get Current Track informations:
		XYZVec updated_track_dir( track->get_m_0(), 1., track->get_m_1());

                //Getting hit errors:
                double hit_error = sqrt(major_axis.Dot(updated_track_dir)*major_axis.Dot(updated_track_dir)+minor_axis.Dot(updated_track_dir)*minor_axis.Dot(updated_track_dir));
		
		//Update Error:
		double d_error =0;
               
		if (hit_error > err[i]){
			d_error = sqrt((hit_error - err[i])*(hit_error - err[i]));
		        err[i] = hit_error; 
			if(_diag>0){
				std::cout<<" Error changed to : "<<err[i]<<std::endl;
			}
		}
		if(_diag>0){
			std::cout<<" Changed in errors of "<<d_error<<std::endl;
		}
        	if( d_error < 0.1){ 
			if(_diag>0){
           			std::cout<<"Considered converged ..."<<std::endl;
			}
           		errors_converged = true;
                }
    
                //Refit with Updated Error:
	       if (_dim==3){
			//Reset fit:
        		track->clear();
		      	LeastSquaresFitter::xyz_fit( _dim, x, y,z, err, track, cov_x);
		      }
	       if (_dim==2){
			//Reset fit:
        		track->clear();
			std::cout<<"clear check  "<<track->get_chisq()<<std::endl;
			LeastSquaresFitter::xy_fit( x, y, err, track, cov_x);
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
