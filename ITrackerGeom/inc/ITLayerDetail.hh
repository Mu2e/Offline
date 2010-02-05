#ifndef ITLAYERDETAIL_HH
#define ITLAYERDETAIL_HH

#include <string>

namespace mu2e {

class ITLayerDetail{

public:
  ITLayerDetail();

  ITLayerDetail( double center_radius_ringIn, double center_radius_ringOut, double epsilonIn,
     double epsilonOut, double halfLength, std::string materialNames
  );
  
  ~ITLayerDetail ();

  double      centerInnerRadiusRing()    const { return _center_radius_ringIn;}
  double      centerOuterRadiusRing()    const { return _center_radius_ringOut;}
  double      stereoAngleInnerRing()     const { return _epsilonIn;}
  double      stereoAngleOuterRing()     const { return _epsilonOut;}
  double      halfLength()               const { return _halfLength; }
  std::string materialName()             const { return _materialNames; }
  
private:

  double _center_radius_ringIn;   //Inner surface radius
  double _center_radius_ringOut;  //Outer surface radius
  double _epsilonIn;              //Inner surface stereo angle
  double _epsilonOut;             //Outer surface stereo angle
  double _halfLength;             //Half z length
  std::string _materialNames;     //G4 material name

};

}  //namespace mu2e

#endif /*ITLAYERDETAIL_HH*/
