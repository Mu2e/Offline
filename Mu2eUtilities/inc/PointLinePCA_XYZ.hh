#ifndef Mu2eUtilities_PointLinePCA_XYZ_hh
#define Mu2eUtilities_PointLinePCA_XYZ_hh

#include "DataProducts/inc/XYZVec.hh"

namespace mu2e {

  class PointLinePCA_XYZ{

  public:
	PointLinePCA_XYZ( XYZVec const& point,       
                XYZVec const& start ,
                XYZVec const& end);
        ~PointLinePCA_XYZ();
	double dca()   const { return _dca;};
        XYZVec const& pca() const { return _pca;}; 
    private:
	
	XYZVec _point;
	XYZVec _start;
	XYZVec _end;
	double _dca;
	XYZVec _pca;
   };
}
#endif
