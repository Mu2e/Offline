#ifndef Mu2eKinKal_KKFitUtilities_hh
#define Mu2eKinKal_KKFitUtilities_hh
//
//  untemplated utiltity classes and functions
//
#include "KinKal/Trajectory/Line.hh"
namespace mu2e {
  class ComboHit;
  class Straw;
  class StrawResponse;
  namespace Mu2eKinKal{
    // function to turn a StrawHit into a Line object
    KinKal::Line hitLine(ComboHit const& ch, Straw const& straw,StrawResponse const& strawresponse);
  }
}
#endif
