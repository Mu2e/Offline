#ifndef RecoDataProducts_DriftCircle_hh
#define RecoDataProducts_DriftCircle_hh

using namespace std;

namespace mu2e {
  
  class DriftCircle{
  public:

//Cosmic Track can be constructed on drift circle - this struct defines the drift cricle parameters , x or y position, z position, radius and error. The radius should be DCA to the line in YZ plane..
            double t0;   //t0 from time cluster
	    double t_drift; //drift time =hit time - t0;
	    double x;   // hit cylinder length
	    double y;   // hit cylinder length 
	    double zc;    // Z co-ord of centre
	    double r;    // Drift Radius
	    double rErr; // Error or radius
            
	    DriftCircle(const double x, const double y, const double zc, const double r, const double rErr): x(x), y (y), zc(zc), r(r), rErr(rErr) {}

	    bool operator<(const DriftCircle& circle) {return zc < circle.zc;}


 };//end DriftCircle class

}//end namespace
#endif
