#ifndef RecoDataProducts_DriftCircle_hh
#define RecoDataProducts_DriftCircle_hh

using namespace std;

namespace mu2e {
  
  class DriftCircle{
  public:


//Cosmic Track can be constructed on drift circle - this struct defines the drift cricle parameters , x or y position, z position, radius and error. The radius should be DCA to the line in YZ plane..
            double t0;   //t0 from time cluster
	    double t_drift; //drift time =hit time - t0;
	    double xy;   // XY co-ord of centre
	    double z;    // Z co-ord of centre
	    double r;    // Drift Radius
	    double rErr; // Error or radius
            double panel; //panel ID

	    DriftCircle(const double xy, const double z, const double r, const double rErr, int panel): xy(xy), z(z), r(r), rErr(rErr), panel(panel) {}

	    bool operator<(const DriftCircle& circle) {return z < circle.z;}


 };//end DriftCircle class

}//end namespace
#endif
