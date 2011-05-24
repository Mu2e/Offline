//
// Extract trajectories from the G4 internals and add them to the event.
// Skip trajectories with too few points.
//
// $Id: addPointTrajectories.cc,v 1.3 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) The data product is filled via two phase construction.
//    This avoids extra copies of the (often long) vector of points.
//    If and when we get move aware STL libraries we can go back
//    to single phase construction.
//

#include "G4AttDef.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4TrajectoryPoint.hh"
#include "Mu2eG4/inc/addPointTrajectories.hh"
#include <iostream>
#include <map>

using namespace std;

namespace mu2e{

  void addPointTrajectories ( const G4Event*             g4event,
                              PointTrajectoryCollection& pointTrajectories,
                              CLHEP::Hep3Vector const&   mu2eOriginInWorld ){

    typedef PointTrajectoryCollection::key_type    key_type;
    typedef PointTrajectoryCollection::mapped_type mapped_type;
    typedef std::map<key_type,mapped_type>         map_type;

    // Check that there is information to be copied.
    G4TrajectoryContainer const* trajectories  = g4event->GetTrajectoryContainer();
    if ( !trajectories ) return;
    TrajectoryVector const& vect = *trajectories->GetVector();
    if ( vect.empty() ) return;

    // Need to pre-sort before insertion into the collection; do it using this map.
    map_type tempMap;

    // Insert mostly empty PointTrajecotry objects into the temporary map.  These contain
    // only the ID and an empty std::vector so they are small.
    for ( size_t i=0; i<vect.size(); ++i){
      G4VTrajectory const& traj = *vect[i];

      int       id(traj.GetTrackID());
      key_type kid(id);

      // Cut if too few points.  Need to make this a variable.
      if ( traj.GetPointEntries() < 5 ) {
        continue;
      }

      tempMap[kid] = PointTrajectory(id);
    }

    // Phase 1 of construction of the data product.  See note 1.
    pointTrajectories.insert( tempMap.begin(), tempMap.end() );

    // Phase 2 of construction. Add the vector of points.
    for ( size_t i=0; i<vect.size(); ++i){
      G4VTrajectory const& traj = *vect[i];

      // Which trajectory are we looking for?
      int       id(traj.GetTrackID());
      key_type kid(id);

      // Locate this trajectory in the data product.
      // It is OK if we do not find it since we applied cuts above.
      PointTrajectoryCollection::iterator iter =  pointTrajectories.find(kid);
      if ( iter == pointTrajectories.end() ){
        continue;
      }
      PointTrajectory& ptraj = iter->second;

      // Add the points.
      for ( int j=0; j<traj.GetPointEntries(); ++j){
        G4VTrajectoryPoint const& pt = *traj.GetPoint(j);
        ptraj.addPoint( pt.GetPosition()-mu2eOriginInWorld );
      }

    } // end phase 2

  } // end addPointTrajectories.

} // end namespace mu2e
