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

int GetDOCASign(XYZVec track_dir, XYZVec point){
        int sign = (track_dir.y()-point.y() < 0) ?  -1 : 1 ;
	return sign;

}




/*-----------------Get X' and Y' and hit errors------------//
The initial track is found and its direction. z' becomes the parametric variable and the fit factorizes into 2 2-D fits, with intercept and slopes (A0, A1 and B0, B1) along x' and y' respectively.  You have to project the hits and their error ellipses onto x' and y'. The following 2 equations do this:

sigma_x' = sqrt(sigma_Maj^2(x'.Maj)^2 + sigma_min^2(x'.Min)^2) i.e. sum or squares of error projections on to x' axis (same for y')
----------------------------------------------------------------------------------*/
/*
We require:
1) x'.y = 0 (y component of x' = 0)
2) x'.z' = 0 (perp)
3) y'.x ' = 0 (perp)
4) y'.z' = 0 (perp)
*/

std::vector<XYZVec>GetAxes(XYZVec TrackDirection){
	std::vector<XYZVec> AxesList;
	
        XYZVec XPrime = ParametricFit::GetXPrime(TrackDirection);
        XYZVec YPrime = ParametricFit::GetYPrime(XPrime, TrackDirection); 
        XYZVec XDoublePrime = ParametricFit::GetXDoublePrime(XPrime, YPrime, TrackDirection);
        XYZVec YDoublePrime = ParametricFit::GetYDoublePrime(XPrime, YPrime, TrackDirection); 
	//ORDERED X'', Y'', Z' IMPORTANT!!!
	AxesList.push_back(XDoublePrime);
	AxesList.push_back(YDoublePrime);
	AxesList.push_back(TrackDirection);
	return AxesList;
}
XYZVec GetXPrime(XYZVec track_dir){
	double Ax= abs(track_dir.x());
	double Ay= abs(track_dir.y());
	double Az= abs(track_dir.z()); 
	if (Ax < Ay){
    		return Ax < Az ? XYZVec(0, -track_dir.z(), track_dir.y()).Unit() : XYZVec(-track_dir.y(), track_dir.x(), 0).Unit();
    	} else {
    		return Ay < Az ? XYZVec(track_dir.z(), 0, -track_dir.x()).Unit() : XYZVec(-track_dir.y(), track_dir.x(), 0).Unit();
    	} 	
} 

XYZVec GetYPrime(XYZVec OrthX, XYZVec track_dir){
	XYZVec OrthY = OrthX.Cross(-1*track_dir); //perp to x' and z'
        return OrthY.Unit();

}

XYZVec GetXDoublePrime(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime){	
	double theta = atan2(XPrime.y(),YPrime.y()); //enforces that X''Dot Y = 0
	XYZVec XDoublePrime = cos(theta)*(XPrime) -1*sin(theta)*(YPrime);
	return XDoublePrime.Unit();
}

XYZVec GetYDoublePrime(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime){
	double theta = atan2(XPrime.y(),YPrime.y());
	XYZVec YDoublePrime = cos(theta)*(YPrime) +sin(theta)*(XPrime);
	return YDoublePrime.Unit();
}

void TestConditions(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime){
	std::cout<<XPrime.y()<<"x'y"<<std::endl;
	std::cout<<XPrime.Dot(ZPrime)<<"x'z' "<<std::endl;;
	std::cout<<YPrime.Dot(XPrime) <<"y'x'"<<std::endl;
	std::cout<<(YPrime.Dot(ZPrime))<<"y'z'"<<std::endl;
	
}

/*----------------Initial Error Estimate------------------//
//    Calculate major and minor error ellipse axes and find projected error for hit//
//--------------------------- ----------------------------*/

//This function takes charge of the error calculation and calls the others, The list is such that it is consistant with axes so: "0" = X'', "1" = Y'' (and of course not used but "2" =Z'.)//
std::vector<double> GetErrors(ComboHit* Hit, XYZVec XAxis, XYZVec YAxis){
	XYZVec major_axis =  MajorAxis(Hit);
        XYZVec minor_axis =  MinorAxis(Hit);
        double errX =  HitErrorX(Hit, major_axis, minor_axis, XAxis);
        double errY =  HitErrorY(Hit, major_axis, minor_axis, YAxis);
        std::vector<double> ErrorsXY;
        ErrorsXY.push_back(errX);
        ErrorsXY.push_back(errY);
        return ErrorsXY;

}
XYZVec MajorAxis(ComboHit* Hit){
      XYZVec const& wdir = Hit->wdir();//direction along wire
      double werr_mag = Hit->wireRes(); //hit major error axis  
      XYZVec major_axis = werr_mag*wdir;
      return major_axis;
}

XYZVec MinorAxis(ComboHit* Hit){
      XYZVec const& wdir = Hit->wdir();//direction along wire
      XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire 
      double terr_mag = Hit->transRes(); //hit minor error axis
      XYZVec minor_axis = terr_mag*wtdir;
      return minor_axis;
}


double HitErrorX(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec XPrime){
        double sigma_w_squared = major_axis.Mag2();      
        double sigma_v_squared = minor_axis.Mag2();
	double sigma_x_track = sqrt(sigma_w_squared*pow(XPrime.Dot(major_axis.Unit()),2)+sigma_v_squared*pow(XPrime.Dot(minor_axis.Unit()),2));
 	return sigma_x_track;
}
 
double HitErrorY(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec YPrime){
        double sigma_w_squared = major_axis.Mag2();
        double sigma_v_squared = minor_axis.Mag2();
	double sigma_y_track = sqrt(sigma_w_squared*pow((YPrime.Dot(major_axis.Unit())),2)+sigma_v_squared*pow((YPrime.Dot(minor_axis.Unit())),2));
 	return sigma_y_track;
}

double TotalHitError(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec XPrime, XYZVec YPrime){
	double errX =  ParametricFit::HitErrorX(Hit, major_axis, minor_axis, XPrime);
        double errY =  ParametricFit::HitErrorY(Hit, major_axis, minor_axis, YPrime);
        return sqrt(pow(errX,2)+pow(errY,2));

}
/*-------------------------Get 2D line -------------------------//
Each 2D line has 2 parameters A_0 and A_1 for x and B_0 and B_1 for y. 
A chi2 = sum of (A0+A1zi) - hi.zi / error_i^2  = 0
Can build 3 matrices for x' and 3 for y'. See documentation for details but alpha = gamma^-1 . beta where alpha are a given 2D 2 parameters.....
----------------------------------------------------------------*/
double GetHitChi2(double A0, double A1, double errorX, XYZVec point, XYZVec XDoublePrime, XYZVec ZPrime){
	double chi_i = ((A0*A0))+((A1*A1))*pow(point.Dot(ZPrime),2)+pow(point.Dot(XDoublePrime),2)-2*A0*point.Dot(XDoublePrime)+2*A0*A1*point.Dot(ZPrime)-2*A1*point.Dot(ZPrime)*point.Dot(XDoublePrime);
	double chi2_i = (chi_i*chi_i)/(errorX*errorX);
	return chi2_i;

}

double GetGlobalChi2(double a0, double a1, double b0, double b1, XYZVec prime_point, double errX, double errY, int DOF){
	double pullX = (a0+(a1*prime_point.z())-prime_point.x())/errX;
	double pullY = (b0+(b1*prime_point.z())-prime_point.y())/errY;
	double chi2 = pullX*pullX + pullY*pullY;
	return chi2/DOF;
	}
	

double GetResidualX(double A0, double A1 ,XYZVec point_prime){
	double resid_x = A0 + A1*(point_prime.z())-point_prime.x();	
	return resid_x;//sqrt(pow(resid_x,2)+pow(resid_y,2));	
}

double GetResidualY( double B0, double B1,  XYZVec point_prime){
	double resid_y = B0 +(B1*point_prime.z())-point_prime.y();
	return resid_y;
}

double GetTotalResidual(double resid_x,double resid_y){
	return sqrt(pow(resid_x,2)+pow(resid_y,2));
}

//-----------------MC Diagnotics--------------------//

double GetMCResidualX(double A0, double A1, XYZVec MCTrackDirection, XYZVec MCX, XYZVec point){
	double resid_x = A0 + A1*point.Dot(MCTrackDirection)-point.Dot(MCX);	
	return resid_x;//sqrt(pow(resid_x,2)+pow(resid_y,2));

}

double GetMCResidualY(double B0, double B1, XYZVec MCTrackDirection, XYZVec MCY, XYZVec point){
	double resid_y = B0 + B1*point.Dot(MCTrackDirection)-point.Dot(MCY);		
	return resid_y;//sqrt(pow(resid_x,2)+pow(resid_y,2));

}


}//end namespace 



