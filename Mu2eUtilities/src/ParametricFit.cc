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
#include "TH1F.h"

using namespace std;
using namespace mu2e;


namespace ParametricFit{

	double GettMin(XYZVec& point, XYZVec& starting_point, XYZVec& end_point){

            double tMin = -(starting_point-point).Dot(end_point-starting_point) /((end_point-starting_point).Mag2());
	    return tMin;

       }

      double LRAmbig( double& dca){
        
	   //XYZVec closestPointOnLine=PointToLineCA( point, line_starting_point, line_end_point);
	double LorR = dca > 0.0 ? 1.0 : -1.0;

	return LorR;
      }

	XYZVec PointToLineCA(XYZVec& point, XYZVec& starting_point, XYZVec& end_point){

          double tMin = GettMin(point, starting_point, end_point);
	 
	  double POCA_x = starting_point.x() + (end_point.x()-starting_point.x())*tMin;
	  double POCA_y = starting_point.y() + (end_point.x()-starting_point.y())*tMin;
	  double POCA_z = starting_point.z() + (end_point.z()-starting_point.z())*tMin;
	  XYZVec closestPointOnLine;
	  closestPointOnLine.SetXYZ(POCA_x ,POCA_y ,POCA_z);
	  return closestPointOnLine;
	}
	
	double PointToLineDCA(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point){
	   double dca;
	   XYZVec closestPointOnLine=PointToLineCA( point, line_starting_point, line_end_point);

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
     
	
	XYZVec pointOnLineFromX(XYZVec lineStartPoint, XYZVec lineEndPoint, double x,XYZVec outputPoint)
	{

  	//Get t for given x coord (line eqn: r = r0 + mt)
 	XYZVec gradient = lineEndPoint - lineStartPoint;  
  	double t = ( x - lineStartPoint.x() ) / gradient.x();


  	//Use t to get all point coords
  	outputPoint.SetXYZ( x , lineStartPoint.y()+(gradient.y()*t) , lineStartPoint.z()+(gradient.z()*t));

  	return outputPoint;

	}



	bool LineToLineCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint, 
  XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, 
  XYZVec& closestPointOnFirstLine, XYZVec& closestPointOnSecondLine)
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
	
	  
	    if( tMin > sqrt((firstLineEndPoint-firstLineStartPoint).Mag2()) ) return false;
	    if( tMin < 0. ) return false;
	    if( sMin > sqrt((secondLineEndPoint-secondLineStartPoint).Mag2()) ) return false;
	    if( sMin < 0. ) return false;
	   

	  //Get the closest approach point on each line using the minimised parameteric scalars
	  closestPointOnFirstLine = firstLineStartPoint + ( tMin * u );
	  closestPointOnSecondLine = secondLineStartPoint + ( sMin * v );

          return true;
        }

	 double LineToLineDCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint,XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, double& dca)
{
	XYZVec closestPointOnFirstLine, closestPointOnSecondLine;
 	bool success = LineToLineCA(firstLineStartPoint, firstLineEndPoint, secondLineStartPoint, secondLineEndPoint, closestPointOnFirstLine, closestPointOnSecondLine);

  
         dca = success ? sqrt((closestPointOnSecondLine-closestPointOnFirstLine).Mag2()) : -1.;
         return success;

}  

/*----------------Initial Error Estimate------------------//
//    Calculate major and minor error ellipse axes and find projected error for hit//
//--------------------------- ----------------------------*/
XYZVec MajorAxis(ComboHit* Hit, XYZVec track_dir){
      XYZVec const& wdir = Hit->wdir();//direction along wire
      double werr_mag = Hit->wireRes(); //hit major error axis 
      XYZVec major_axis = werr_mag*wdir;
    
      return major_axis;
}

XYZVec MinorAxis(ComboHit* Hit, XYZVec track_dir){
    
      XYZVec const& wdir = Hit->wdir();//direction along wire
      XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire 
      double terr_mag = Hit->transRes(); //hit minor error axis
      XYZVec minor_axis = terr_mag*wtdir;
      return minor_axis;
}



/*-----------------Hit Errors in X' and Y' ------------//
The initial track is found and its direction. z' becomes the parametric variable and the fit factorizes into 2 2-D fits, with intercept and slopes (A0, A1 and B0, B1) along x' and y' respectively.  You have to project the hits and their error ellipses onto x' and y'. The following 2 equations do this:

sigma_x' = sqrt(sigma_Maj^2(x'.Maj)^2 + sigma_min^2(x'.Min)^2) i.e. sum or squares of error projections on to x' axis (same for y')
----------------------------------------------------------------------------------*/

double HitErrorX(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir){
	XYZVec track_x(track_dir.x(),0,0);
	XYZVec unit = track_x.Unit();
	//XYZVec x_track( track_dir.x(),0,0);
        double sigma_w_squared = major_axis.Mag2();
        double sigma_v_squared = minor_axis.Mag2();
	double sigma_x_track = sqrt(sigma_w_squared*pow(unit.Dot(major_axis.Unit()),2)+sigma_v_squared*pow(unit.Dot(minor_axis.Unit()),2));
	
	
 	return sigma_x_track;
}
 

double HitErrorY(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir){
	XYZVec y_track(0, track_dir.y(), 0);
        double sigma_w_squared = major_axis.Mag2();
        double sigma_v_squared = minor_axis.Mag2();
	double sigma_y_track = sqrt(sigma_w_squared*pow((y_track.Unit().Dot(major_axis.Unit())),2)+sigma_v_squared*pow((y_track.Unit().Dot(minor_axis.Unit())),2));
 	return sigma_y_track;
}

double TotalHitError(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir){
	double errX =  ParametricFit::HitErrorX(Hit, major_axis, minor_axis, track_dir);
      
      double errY =  ParametricFit::HitErrorY(Hit, major_axis, minor_axis, track_dir);
      
      return sqrt(pow(errX,2)+pow(errY,2));

}
/*-------------------------Get 2D line -------------------------//
Each 2D line has 2 parameters A_0 and A_1 for x and B_0 and B_1 for y. 
A chi2 = sum of (A0+A1zi) - hi.zi / error_i^2  = 0
Can build 3 matrices for x' and 3 for y'. See documentation for details but alpha = gamma^-1 . beta where alpha are a given 2D 2 parameters.....
----------------------------------------------------------------*/

int GetDOCASign(XYZVec track_dir, XYZVec point){
	
        int sign = (track_dir.y()-point.y() < 0) ?  -1 : 1 ;
	return sign;

}
double GetResidual(double A0, double A1, double B0, double B1, XYZVec track_dir, XYZVec point){

	XYZVec track_x(track_dir.x(),0,0);
	XYZVec track_y(0, track_dir.y(),0);
	XYZVec track_z(0,0, track_dir.z());
	double resid_x = point.Dot(track_x.Unit()) - A0 - A1*point.z();
	
	double resid_y = point.Dot(track_y.Unit()) - B0 - B1*point.z();
	std::cout<<" Resid "<<resid_x<<" "<<resid_y<<" "<<sqrt(pow(resid_x,2)+pow(resid_y,2))<<std::endl;
	return resid_x;//sqrt(pow(resid_x,2)+pow(resid_y,2));
	
}
double GetResidualError(XYZVec Major_Axis, XYZVec Minor_Axis, XYZVec track_direction, XYZVec& point, double DOCA){
        XYZVec track_x(track_direction.x(),0,0);
        double R = DOCA;
        double RA2 = Major_Axis.Mag2();
        double RB2 = Minor_Axis.Mag2();
        //Phi =  angle between x' and R
	double Phi = acos((point.Unit()).Dot(track_x.Unit()));
	
	//Theta =  angle between RA and x'i.e. around the vertial direction (polar)
	double Theta = acos((Major_Axis.Unit()).Dot(track_x.Unit()));
	
	double sigma_R_num =  pow(R,2)*((RA2/2)+(RB2/2)+((RA2+RB2)/2)*cos(2*Phi)*cos(2*Theta));
	//std::cout<<"sigma N"<<sigma_R_num<<std::endl;
	
	double sigma_R_den = (pow(RA2,2)*pow(cos(Theta),2)*pow(sin(Theta),2)+pow(RB2,2)*pow(sin(Theta),2)*pow(cos(Theta),2)+2*RA2*RB2*(pow(cos(Theta),2)+pow(sin(Theta),2)));
	//std::cout<<"Sigma D "<<sigma_R_den<<std::endl;
	double sigma_R = sqrt(sigma_R_num/sigma_R_den);
	std::cout<<" Resid Error "<<sigma_R<<std::endl;
	return sigma_R;
}


/*-------Convert the 2 fits into full track fit---------*/
//make a 3D comsic track from the 2 2D lines found from above
   void Construct3DTrack(StraightTrack* xyLineFit, StraightTrack* zrLineFit, CosmicTrack* track) {
	
       //return track;

  }




}//end namespace 



