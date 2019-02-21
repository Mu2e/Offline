// Object to perform track fit to straw hits
//
// $Id: StriaghtTrackFit code
// $Author: S Middleton
// $Date: Nov 2018
//
// Mu2e
#include "TrkReco/inc/StraightTrackFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"

//For Drift:
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PanelStateIterator.hh"
//Least Squares Fitter:
#include "Mu2eUtilities/inc/LeastSquaresFitter.hh"

//ROOT:
#include "TMatrixD.h"
#include "Math/VectorUtil.h"


using namespace std;
//using namespace boost::accumulators;
//using namespace ROOT::Math::VectorUtil;

namespace mu2e
{
  StraightTrackFit::StraightTrackFit(fhicl::ParameterSet const& pset) :
    _dim(pset.get<int>("dimension",3)),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnsh(pset.get<unsigned>("minNStrawHits",4)),
    _minxyresid(pset.get<float>("minXYResid",10.)),
    //_maxniter(pset.get<unsigned>("maxniter",100)),
    //_minzsep(pset.get<float>("minzsep",100.0)),
    //_maxzsep(pset.get<float>("maxzsep",500.0)),
    _maxdxy(pset.get<float>("maxdxy",1000.0)),
    _maxchi2xy(pset.get<float>("maxchi2xy",1000.0))   
    {}

//destructor
    StraightTrackFit::~StraightTrackFit(){}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool StraightTrackFit::initStraightTrack(StraightTrackFinderData& TrackData) {
    if(_debug>0){
    	std::cout<<"Initializing ST Fit ..."<<std::endl;
    }
    bool is_ok(false);
    FitChi2(TrackData);
    is_ok = TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackOK);
    return is_ok;
  }


  /*-------------Init Line-------------------------//
// Begin with straight line between first and last points on plane//
  //----------------------------------------------*/
  StraightTrack* StraightTrackFit::InitLine(const ComboHit *FirstP1, const ComboHit *LastP1,StraightTrack* line) {
      
      XYVec Pos_Start = (XYVec(FirstP1->pos().x(),FirstP1->pos().y()));
      XYVec Pos_End = (XYVec(LastP1->pos().x(),LastP1->pos().y()));
      line->set_m_0(( Pos_End.x() - Pos_Start.x()) / (Pos_End.y() - Pos_Start.y() ));
      //line->set_m_1(( Pos_End.z() - Pos_Start.z()) / (Pos_End.y() - Pos_Start.y() ));
      line->set_c_0(Pos_Start.x() - ( Pos_Start.y() * line->get_m_0()) );
      if(_debug>0){
	      
	      std::cout<<"initialized line .."<<line->get_m_0()<<"  "<<line->get_c_0()<<std::endl;
      }
      return line;
    } 

 //--------------Fit-----------------//
//Top call to Fitting routines....
//-------------------------------------------// 

  void StraightTrackFit::BeginFit(StraightTrackFinderData& TrackData){
      if(_debug>0){
      	std::cout<<" Beginning Fit ..." << std::endl;
      }
      //Clear Previous Flags:
      TrackData._tseed._status.clear(TrkFitFlag::StraightTrackOK);
      
    // Initialize:
    bool init(false);
    if (!TrackData._tseed._status.hasAllProperties(TrkFitFlag::StraightTrackInit)) {
      init = true;
      if (initStraightTrack(TrackData))
	TrackData._tseed._status.merge(TrkFitFlag::StraightTrackInit);
      else
	return;
    }
    
    //Start Chi2 Fitting:
    if (!init)FitChi2(TrackData);
  }


/*------------------------------Chi 2 Fit----------------------------//
//   Adds Chi-2 optimization to fitting routine //
//   Refits and adjusts track fit paramters by weights//
//------------------------------------------------------------------*/
void StraightTrackFit::FitChi2(StraightTrackFinderData& StraightTrackData) {   
    //int  nHits(StraightTrackData._chHitsToProcess.size());
    //ComboHit*     firstPt(0), *lastPt(0) ;
    
    StraightTrack* track = &StraightTrackData._tseed._track; 
    
    //loop over faces
    //firstPt = &StraightTrackData._chHitsToProcess[0];
    //lastPt = &StraightTrackData._chHitsToProcess[nHits];


    //track = InitLine(firstPt, lastPt,track);
    // might use later- set initial value high so it becomes best fit..
    
   //first perform the chi2 fit assuming all hits have same error
   StraightTrack* refined_track = refineFit(StraightTrackData, track, 0);
   if(_debug>0){
   	std::cout<<"Refined Final Track..."<<refined_track->get_m_0()<<" c "<<refined_track->get_c_0() << " chi2 "<<refined_track->get_chisq_dof()<<std::endl;
    }
   //now that gradient and offset are more accurate, repeat the fit using the orientation of the hit to estimate more accurately the expected uncertanty:
    
    //if track is "good" add to list of good tracks:
    
   if (goodTrack(refined_track)) StraightTrackData._tseed._status.merge(TrkFitFlag::StraightTrackOK);  


  
}

StraightTrack* StraightTrackFit::refineFit(StraightTrackFinderData& trackData,  StraightTrack* initial_fit, int WeightMode){
    //set up cosmic track:
    StraightTrack* track = &trackData._tseed._track; 
    
    int  nHits(trackData._chHitsToProcess.size());
    int        nXYSh(0);
    ComboHit*     hitP1(0);
   
    //For Fit:
    std::vector<double> x, y, err, z;
    // Store "good" hits:
    std::vector<ComboHit*> combo_hit_list;
    //For Error ellipses:
    std::vector<XYZVec> ma, mi;
    //covarience TODO:change for 3D
    TMatrixD cov_x(_dim, _dim); 
    
    //loop over hits
    for (int f1=1; f1<nHits; ++f1){
      hitP1 = &trackData._chHitsToProcess[f1];
      double tdrift = hitP1->driftTime();
      std::cout<<"Hit time "<<hitP1->time()<<" Drift Time of Hit "<<tdrift<<" TClust"<<trackData._tseed._t0._t0<<std::endl;
      if (!use(*hitP1) ) continue;
     
      //Fill points info:
      x.push_back(hitP1->pos().x());
      y.push_back(hitP1->pos().y());
      z.push_back(hitP1->pos().z()); 
      
      //caculate initial fit error estimate (ball-park....more accurately updated later):
      double initial_hit_error_guess = sqrt((hitP1->wireRes()*hitP1->wdir().x()*hitP1->wireRes()*hitP1->wdir().x()) +hitP1->wireRes()*hitP1->wdir().y()*hitP1->wireRes()*hitP1->wdir().y()+hitP1->wireRes()*hitP1->wdir().z()*hitP1->wireRes()*hitP1->wdir().z());
      
      err.push_back(initial_hit_error_guess);

      //Increase the StawaHit counter
      nXYSh += hitP1->nStrawHits();
      //Fill "good" hit list:
      combo_hit_list.push_back(hitP1);
      }
      //if we collected enough points update the results
    if (nXYSh >= _minnsh){

              //Initial Guess based on above:
	      if (_dim==2){
	      	LeastSquaresFitter::xy_fit( x, y, err, track, cov_x);
	      }
              if (_dim==3){
	      	LeastSquaresFitter::xyz_fit( x, y,z, err, track, cov_x);
	      }
              //Calculate Errors and Refine:
              for (int f1=0; f1<static_cast<int>(combo_hit_list.size()); ++f1){
		      ComboHit* hitP = combo_hit_list[f1];
		      //For Error Ellipses:
		      XYZVec const& wdir = hitP->wdir();//direction along wire
		      XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire
		      double werr_mag = hitP->wireRes(); //hit major error axis 
		      double terr_mag = hitP->transRes(); //hit minor error axis
		      XYZVec major_axis = werr_mag*wdir;
		      XYZVec minor_axis = terr_mag*wtdir;
		    
	      
		      ma.push_back(major_axis);
		      mi.push_back(minor_axis); 
              } 
  
	  
	     UpdateFitErrors(x,y,z, err, track, cov_x, ma, mi);
      
     } 

return track;
  }

/*

-----------------Refine Fit in XY----------------//
//Refines the fit in xy and updates chi2 information//
//-------------------------------------------------
StraightTrack* StraightTrackFit::refineFit(StraightTrackFinderData& trackData,  StraightTrack* initial_fit, int WeightMode){
    if(_debug>0){
        std::cout<<" Refining XY Fit..."<<std::endl;
        
    }
    
    StraightTrack* track = &trackData._tseed._track; 
    
    int  nHits(trackData._chHitsToProcess.size());
    float      wt(1.), resid;
    int        minNReducedChi2Points(15);//TODO: Optimize
    int        nXYSh(0);
    ComboHit*     hitP1(0);
    int nhits=0;

    //For Fit:
    std::vector<double> x, y, y_err, z;
  
    std::vector<XYZVec> ma, mi;
    TMatrixD cov_x(2,2); 
   
    if(_debug>0){
   	std::cout<<"Fitting in XY for Chi2..."<<std::endl;
   }
    //loop over hits
    for (int f1=1; f1<nHits-2; ++f1){
      hitP1 = &trackData._chHitsToProcess[f1];
      if (!use(*hitP1) )    continue;
      

      //For Error Ellipses:
      XYZVec const& wdir = hitP1->wdir();//direction along wire
      XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire
      double werr_mag = hitP1->wireRes(); //hit major error axis 
      double terr_mag = hitP1->transRes(); //hit minor error axis
      XYZVec major_axis = werr_mag*wdir;
      XYZVec minor_axis = terr_mag*wtdir;
  
      HitInfo_t  indexBestComboHit;
      
      float      minResid(_minxyresid);
      if (WeightMode == 0 ) minResid = _maxdxy;

      //Update Fit:
      double cx = track->get_c_0();
      double mx = track->get_m_0();
      //Get Initial Track Direction:
      XYZVec track_dir(1., mx, 0.);

      //Get Hit Error:
      double hit_error = sqrt(major_axis.Cross(track_dir).mag2()+minor_axis.Cross(track_dir).mag2());

      //get a perpendicular residual:
      resid = fabs(hitP1->pos().y()-(mx*hitP1->pos().x()+cx))/sqrt(1.0+mx*mx);
       if(_debug>0){ 
	    std::cout<<"Residual Is : "<<resid<<" Minimum Is : " <<minResid<<std::endl;
	}
	//if residual of new line is better than previous residual then added to hits:
	if (resid < minResid) { 
          
	  indexBestComboHit.face          = f1;
	  indexBestComboHit.panel         = hitP1->strawId().uniquePanel();
	  
	  //update minResid
	  minResid          = resid;
	  hitP1->_xyWeight  = wt;
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
         
         x.push_back(hitP1->pos().x());
         y.push_back(hitP1->pos().y());
         z.push_back(hitP1->pos().z()); 
         
         ma.push_back(major_axis);
         mi.push_back(minor_axis);     
       
         y_err.push_back(hit_error);
	
	 //remove the outlier flag
	 hitP1->_flag.clear(StrawHitFlag::outlier);

	 //increase the StawaHit counter
	 nXYSh += hitP1->nStrawHits();

	 if (nhits < minNReducedChi2Points) continue;
      }
	    
    }//end loop over the faces

    //if we collected enough points update the results
    if (nXYSh >= _minnsh){


      //LeastSquaresFitter::xy_fit( x, y, y_err, track, cov_x);
      LeastSquaresFitter::xyz_fit( x, y,z, y_err, track, cov_x);
      //UpdateFitErrors(x,y,y_err, track, cov_x, ma, mi);
      
     } 
    
    return track;
  }
*/
void StraightTrackFit::UpdateFitErrors(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> err, StraightTrack* track,TMatrixD cov_x, std::vector<XYZVec> major,std::vector<XYZVec> minor){

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
                double hit_error = sqrt(major_axis.Cross(updated_track_dir).mag2()+minor_axis.Cross(updated_track_dir).mag2());
		
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
		      	LeastSquaresFitter::xyz_fit( x, y,z, err, track, cov_x);
		      }
	       if (_dim==2){
			LeastSquaresFitter::xy_fit( x, y, err, track, cov_x);
		      }

       }//end while 
    }//end for
 
}//end update

void StraightTrackFit::add_drift(StraightTrackFinderData& trackData){
	//TrkT0 t0 = trackData._tseed._t0 ; //Get t0 from tclust
	//std::cout<<" Muon Time Cluster Time  = "<<t0<<std::endl;
}

bool StraightTrackFit::goodTrack(StraightTrack* track)
  { 
    
    
    if(track->get_chisq() < _maxchi2xy) return true;
    else return false;
    
  }

/*--------------USE HIT?---------------------------//
//          Checks if flag                        //
//------------------------------------------------*/
bool StraightTrackFit::use(const ComboHit& thit) const 
  {
    return (!thit._flag.hasAnyProperty(_dontuseflag));
  }



}//end namespace
