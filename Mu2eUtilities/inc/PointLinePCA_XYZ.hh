#ifndef Mu2eUtilities_PointLinePCA_XYZ_hh
#define Mu2eUtilities_PointLinePCA_XYZ_hh

#include "Offline/DataProducts/inc/GenVector.hh"

namespace mu2e {

  class PointLinePCA_XYZ{

  public:
        PointLinePCA_XYZ( XYZVectorF const& point,
                XYZVectorF const& start ,
                XYZVectorF const& end);
        ~PointLinePCA_XYZ();
        double dca()   const { return _dca;};
        XYZVectorF const& pca() const { return _pca;};

    private:

        XYZVectorF _point;
        XYZVectorF _start;
        XYZVectorF _end;
        double _dca;
        XYZVectorF _pca;
   };
}
#endif
