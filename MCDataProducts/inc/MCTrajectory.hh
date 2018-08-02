#ifndef MCDataProducts_MCTrajectory_hh
#define MCDataProducts_MCTrajectory_hh
//
// A trajectory defined as a collection of 3D points + time + kinetic energy.
// The points are defined in the Mu2e coordinate system.
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/map_vector.h"
#include <vector>

namespace mu2e {

  class MCTrajectory{

  public:

    typedef cet::map_vector_key key_type;
    typedef MCTrajectoryPoint Point;

    // Default c'tor needed by ROOT persistency.
    MCTrajectory():sim_(),points_(){
    }

    // Usual c'tor.
    MCTrajectory( art::Ptr<SimParticle> const& sim,  std::vector<Point> const& points):
      sim_(sim),
      points_(points){
    }

    // Accept compiler written d'tor, copy c'tor and assignment operator.

    // Accessors
    art::Ptr<SimParticle> const& sim() const { return sim_; }
    art::Ptr<SimParticle>      & sim()       { return sim_; }

    int    simid() const { return sim_.key();     }
    size_t size()  const { return points_.size(); }

    std::vector<Point> const& points() const { return points_;}
    std::vector<Point>      & points()       { return points_;}

    // The following c'tor and addPoints method form a two-phase c'tor.
    // This is needed in addPointTrajectories to avoid copying the vector
    // of points twice.  Once we have move aware containers we can
    // get rid of these.
    MCTrajectory(  art::Ptr<SimParticle> const& sim ):
      sim_(sim),
      points_(){
    }

    // Second phase of the two phase constructor
    void addPoint(const Point& p){
      points_.push_back(p);
    }

  private:

    art::Ptr<SimParticle> sim_;
    std::vector<Point> points_;

  };

}

#endif /* MCDataProducts_MCTrajectory_hh */
