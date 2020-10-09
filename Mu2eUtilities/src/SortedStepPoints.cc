//
//
//
// Original author Rob Kutschke
//

#include <cmath>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"


#include "Mu2eUtilities/inc/SortedStepPoints.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"

using namespace std;

namespace mu2e {

  SortedStepPoints::SortedStepPoints ( SimParticleCollection::key_type trackId,
                                       StepPointMCCollection const&    allSteps ):
    trackId_(trackId),
    steps_(),
    midByZ_(0){

    // Z values of most upstream and most downdstream hits.
    double zmin( 999999.);
    double zmax(-999999.);

    // Find all StepPointMCs created by this track.
    for ( StepPointMCCollection::const_iterator i=allSteps.begin();
          i != allSteps.end(); ++i ){
      StepPointMC const& step = *i;
      double z = step.position().z();
      zmin = ( z < zmin ) ? z : zmin;
      zmax = ( z > zmax ) ? z : zmax;
//       cout << "Adding: " << z << endl;
      steps_.push_back( &step );
    }

    // Middle of the z range of the track.
    double zmid = (zmin + zmax)/2.;
//     cout << "Mid: "
//          << trackId_ << " "
//          << steps_.size() << " "
//          << allSteps.size() << " "
//          << zmin << " "
//          << zmax << " "
//          << zmid
//          << endl;

    for ( int i=0;i<5; ++i ){
//       cout << "Step check: "
//            << i << " "
//            << steps_[i]
//            << endl;
    }



    // Find the StepPointMC whose z is closest to zmid.
    double dmin(999999.);
    for ( vector<StepPointMC const*>::const_iterator i=steps_.begin();
          i != steps_.end(); ++i ){
      StepPointMC const* step = *i;
      double d = std::abs(step->position().z() - zmid);
      if ( d < dmin ){
        dmin = d;
        midByZ_ = step;
	//        cout << "             Mark" <<endl;
      }
//       cout << "Test: "
//            << step->position().z() << " "
//            << step << " "
//            << d << " "
//            << dmin
//            << endl;

    }
//     cout << "Answer: "
//          << midByZ_->position().z() << " "
//          << midByZ_->position().z()-dmin << " "
//          << endl;
  }

}
