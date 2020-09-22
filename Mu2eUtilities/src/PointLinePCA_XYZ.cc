#include <iosfwd>

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

#include "Mu2eUtilities/inc/PointLinePCA_XYZ.hh"

using namespace std;
namespace mu2e{
	PointLinePCA_XYZ::PointLinePCA_XYZ( XYZVec const& point,
                          XYZVec const& start,
                          XYZVec const& end):
    _point(point),
    _start(start),
    _end(end),
    _dca(0.){

	  double tMin = -(start-point).Dot(end-start) /((end-start).Mag2());
	  double POCA_x = start.x() + (end.x()-start.x())*tMin;
	  double POCA_y = start.y() + (end.x()-start.y())*tMin;
	  double POCA_z = start.z() + (end.z()-start.z())*tMin;
	  XYZVec closestPointOnLine;
	  closestPointOnLine.SetXYZ(POCA_x ,POCA_y ,POCA_z);
	 
	 _pca = closestPointOnLine;
	 _dca = sqrt((closestPointOnLine-point).Mag2());
         }

	PointLinePCA_XYZ::~PointLinePCA_XYZ(){}
}
