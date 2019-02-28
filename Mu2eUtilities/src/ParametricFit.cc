#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
using namespace std;
//using namespace mu2e;


namespace ParametricFit{
 
	void PointToLineCA(XYZVec& point, XYZVec& starting_point, XYZVec& end_point, XYZVec& closestPointOnLine, bool finite){

	  double tMin = -(starting_point-point).Dot(end_point-starting_point) /((end_point-starting_point).Mag2());

	  //If line is finite, check that this closest point lies between the two points specified:

	  
	    if ( (tMin<0.) || (tMin>1.) ){
		throw "error- tMin not in range";
	   }
	  double POCA_x = starting_point.x() + (end_point.x()-starting_point.x())*tMin;
	  double POCA_y = starting_point.y() + (end_point.x()-starting_point.y())*tMin;
	  double POCA_z = starting_point.z() + (end_point.z()-starting_point.z())*tMin;
	  //Put tMin back into equation to get point along line that is closest to p
	  closestPointOnLine.SetXYZ( POCA_x ,POCA_y ,POCA_z);

	}
	
	double PointToLineDCA(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point, double dca, bool finiteLine){

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
        
	double GetResidual(){
 		return 0.;
	}

	/*
	CosmicTrack ConstructTrack(){
	
	}

	*/
	
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

}//end namespace 



