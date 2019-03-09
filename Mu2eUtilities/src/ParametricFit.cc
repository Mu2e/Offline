//S Middleton 2019
//Description: Stores all the functions associated with building and interpretting parametric line equations for the purpose of cosmic track based alignment

#include "Mu2eUtilities/inc/ParametricFit.hh"

#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/DriftCircle.hh"
//c++
#include<exception>
#include<bitset>
//ROOT
#include "TMath.h"
#include "TF1.h"

using namespace std;
using namespace mu2e;


namespace ParametricFit{

	double GettMin(XYZVec& point, XYZVec& starting_point, XYZVec& end_point){

            double tMin = -(starting_point-point).Dot(end_point-starting_point) /((end_point-starting_point).Mag2());
	    return tMin;

       }

      double LRAmbig(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point, double& dca, bool finiteLine){
        
	   XYZVec closestPointOnLine;
	   PointToLineCA( point, line_starting_point, line_end_point, closestPointOnLine, finiteLine);
	double LorR = dca > 0.0 ? 1.0 : -1.0;

	return LorR;
      }

	void PointToLineCA(XYZVec& point, XYZVec& starting_point, XYZVec& end_point, XYZVec& closestPointOnLine, bool finite){

           double tMin = GettMin(point, starting_point, end_point);
	 
	    if ( (tMin<0.) || (tMin>1.) ){
		std::cout<< "error- tMin not in range"<<std::endl;
	   }
         
	  double POCA_x = starting_point.x() + (end_point.x()-starting_point.x())*tMin;
	  double POCA_y = starting_point.y() + (end_point.x()-starting_point.y())*tMin;
	  double POCA_z = starting_point.z() + (end_point.z()-starting_point.z())*tMin;
	  //Put tMin back into equation to get point along line that is closest to p
	  closestPointOnLine.SetXYZ( POCA_x ,POCA_y ,POCA_z);

	}
	
	double PointToLineDCA(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point, double& dca, bool finiteLine){

	   XYZVec closestPointOnLine;
	   PointToLineCA( point, line_starting_point, line_end_point, closestPointOnLine, finiteLine);

	  //Return magnitude of vector between the two points as DCA
	  dca =  sqrt((closestPointOnLine-point).Mag2());
          return dca;
          }

        XYZVec ParellelVector(XYZVec lineStartPoint, XYZVec lineEndPoint){
		XYZVec parellel;
		double tx = lineStartPoint.x() -  lineEndPoint.x();
		double ty = lineStartPoint.y()-lineStartPoint.y();
		double tz = lineEndPoint.z() -  lineEndPoint.z();
		parellel.SetXYZ( tx,ty,tz);
		return parellel;

	} 
     
	
	XYZVec pointOnLineFromX(XYZVec lineStartPoint, XYZVec lineEndPoint, double x,XYZVec outputPoint,bool finiteLine)
	{

  	//Get t for given x coord (line eqn: r = r0 + mt)
 	XYZVec gradient = lineEndPoint - lineStartPoint;  
  	double t = ( x - lineStartPoint.x() ) / gradient.x();

  	
  	if(finiteLine) {
    	throw "t not in range";
  	}

  	//Use t to get all point coords
  	outputPoint.SetXYZ( x , lineStartPoint.y()+(gradient.y()*t) , lineStartPoint.z()+(gradient.z()*t));

  	return outputPoint;

	}



	bool LineToLineCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint, 
  XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, 
  XYZVec& closestPointOnFirstLine, XYZVec& closestPointOnSecondLine, bool finiteLine)
        {


	  XYZVec & p0 = firstLineStartPoint;
	  XYZVec u = ( firstLineEndPoint - firstLineStartPoint ).Unit();

	  XYZVec & q0 = secondLineStartPoint;
	  XYZVec v = ( secondLineEndPoint - secondLineStartPoint ).Unit();

	  XYZVec d = p0 - q0;

	  //Closest approach is when magnitude of d is minimized w.r.t. t and s
	  //These are given by:
	  //  t_min = [-d.u + (d.v)(u.v)] / [1 - (u.v)^2]
	  //  s_min = [ d.v - (d.u)(u.v)] / [1 - (u.v)^2]

	  //Pre-compute repeated terms for speed
	  double u_dot_v = u.Dot(v);
	  double d_dot_u = d.Dot(u);
	  double d_dot_v = d.Dot(v);
	  double denominator = 1 - (u_dot_v*u_dot_v);

	  //Calculate t and s values at closest approach 
	  double tMin = ( -d_dot_u + (d_dot_v*u_dot_v) ) / denominator;
	  double sMin = (  d_dot_v - (d_dot_u*u_dot_v) ) / denominator;
	
	  if(finiteLine) {
	    if( tMin > sqrt((firstLineEndPoint-firstLineStartPoint).Mag2()) ) return false;
	    if( tMin < 0. ) return false;
	    if( sMin > sqrt((secondLineEndPoint-secondLineStartPoint).Mag2()) ) return false;
	    if( sMin < 0. ) return false;
	   }

	  //Get the closest approach point on each line using the minimised parameteric scalars
	  closestPointOnFirstLine = firstLineStartPoint + ( tMin * u );
	  closestPointOnSecondLine = secondLineStartPoint + ( sMin * v );

          return true;
        }

	 double LineToLineDCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint,XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, double& dca, bool finiteLine)
{
	XYZVec closestPointOnFirstLine, closestPointOnSecondLine;
 	bool success = LineToLineCA(firstLineStartPoint, firstLineEndPoint, secondLineStartPoint, secondLineEndPoint, closestPointOnFirstLine, closestPointOnSecondLine,finiteLine);

  
         dca = success ? sqrt((closestPointOnSecondLine-closestPointOnFirstLine).Mag2()) : -1.;
         return success;

}  

/*-----------------------------Initial Error Estimate------------------------------//
//    Calculate major and minor error ellipse axes and find projected error for hit//
//---------------------------------------------------------------------------------*/
XYZVec MajorAxis(ComboHit* Hit, XYZVec track_dir){
      XYZVec const& wdir = Hit->wdir();//direction along wire
      double werr_mag = Hit->wireRes(); //hit major error axis 
      XYZVec major_axis = werr_mag*wdir;
      //double hit_error = sqrt(major_axis.Dot(track_dir)*major_axis.Dot(track_dir)+minor_axis.Dot(track_dir)*minor_axis.Dot(track_dir));
      return major_axis;
}

XYZVec MinorAxis(ComboHit* Hit, XYZVec track_dir){
      //For Error Ellipses:
      XYZVec const& wdir = Hit->wdir();//direction along wire
      XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire 
      double terr_mag = Hit->transRes(); //hit minor error axis
      XYZVec minor_axis = terr_mag*wtdir;
      return minor_axis;
}

/*--------------------------Hit Errors in X' and Y' ------------------------------//
The initial track is found and its direction. z' becomes the parametric variable and the fit factorizes into 2 2-D fits, with intercept and slopes (A0, A1 and B0, B1) along x' and y' respectively.  You have to project the hits and their error ellipses onto x' and y'. The following 2 equations do this:

sigma_x' = sqrt(sigma_Maj^2(x'.Maj)^2 + sigma_min^2(x'.Min)^2) i.e. sum or squares of error projections on to x' axis (same for y')
----------------------------------------------------------------------------------*/

double HitErrorX(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir){
	XYZVec x_track( track_dir.x(),0,0);
        double sigma_w_squared = major_axis.Mag2();
        double sigma_v_squared = minor_axis.Mag2();
	double sigma_x_track = sqrt(sigma_w_squared*pow((x_track.Dot(major_axis.Unit())),2)+sigma_v_squared*pow((x_track.Dot(minor_axis.Unit())),2));
 	return sigma_x_track;
}

double HitErrorY(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir){
	XYZVec y_track(0, track_dir.y(), 0);
        double sigma_w_squared = major_axis.Mag2();
        double sigma_v_squared = minor_axis.Mag2();
	double sigma_y_track = sqrt(sigma_w_squared*pow((y_track.Dot(major_axis.Unit())),2)+sigma_v_squared*pow((y_track.Dot(minor_axis.Unit())),2));
 	return sigma_y_track;
}




/*

double ResidualErrorX(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir){

	

}*/


/*---------------------------------Get 2D line ----------------------------------//
Each 2D line has 2 parameters A_0 and A_1 for x and B_0 and B_1 for y. 
A chi2 = sum of (A0+A1zi) - hi.zi / error_i^2  = 0
Can build 3 matrices for x' and 3 for y'. See documentation for details but alpha = gamma^-1 . beta where alpha are a given 2D 2 parameters.....
--------------------------------------------------------------------------------*/




/*
StraightTrack Calculate2DLineFits(std::vector<ComboHits> list_of_hits, std::vector<double> hit_errors, double pValCut)
{
    StraightTrack cal2DLines;
    // Number of hits
    int nHits = list_of_points.size();
    double S(0), Sz(0), Sxy(0), Szz(0), Sxyxy(0), Sxyz(0);
    for(int hit = 0; hit < nHits; hit++){

      double z = (list_of_points.at(hit)).z();
      double x = (list_of_points.at(hit)).x(); 
      double y = (list_of_points.at(hit)).x(); 
      double err2 = pow((list_of_errors.at(hit)), 2);
    
      S   += 1./err2;
      Sx  += x / err2;
      Sy  += y / err2;
      Sx2 += x*x / err2;
      Sxy += x*y / err2;
      Sy2 += y*y / err2;

    }
    double sigXX   =  Sx2 - (Sx*Sx); 
    double sigXY   = Sxy - Sx*Sy; 
    double sigYY   = Sy2 - Sy*Sy; 
    double chi2  = sigYY()*sigXX() - sigXY()*sigXY();
    // Number of LR combinations (2^N or 1 if using truth)
    int nLRCombos = pow(2,nHits);

}

*/

/*-------Convert the 2 fits into full track fit---------*/
//make a 3D comsic track from the 2 2D lines found from above
   void Construct3DTrack(StraightTrack* xyLineFit, StraightTrack* zrLineFit, CosmicTrack* track) {
	
       //return track;

  }




}//end namespace 



