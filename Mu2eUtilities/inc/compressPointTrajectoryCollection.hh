#ifndef Mu2eUtilities_compressPointTrajectories_hh
#define Mu2eUtilities_compressPointTrajectories_hh

//
// Compress a PointTrajectoryCollection.
//
// $Id: compressPointTrajectoryCollection.hh,v 1.1 2011/12/16 23:14:20 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/12/16 23:14:20 $
//
// Contact person Rob Kutschke
//
// Notes:
// 1) This function copies an existing PointTrajectoryCollection to a new one, but removes
//    PointTrajectories that are flagged as uninteresting.
//
// 2) The existing ( Dec 2011 ) PointTrajectory class contains on art::Ptr objects but the
//    next version will.  Therefore it may be necessary to reseat the art::Ptr objects
//    to point to new SimParticleCollections.  The first two arguments anticipate this feature.
//
// 3) The arguments are:
//    1 - the art ProductID of the output collection
//    2 - a productGetter object to get items from the output collection
//    3 - the input collection
//    4 - the object that knows whether to keep or delete each item - see note 7.
//    5 - the output collection.
//
// 4) The SELECTOR template argument can be any class that supports the following:
//       bool operator[]( cet::map_vector_key) const;
//
//    One example of a type that will work is:
//       cet::map_vector<bool>
//    For other examples see Filters/src/HitsInConversionTimeWindow_module.cc
//
//    Note that std::map<cet::map_vector_key,bool> will not work because it does not have
//    a const operator[].
//

#include "MCDataProducts/inc/PointTrajectoryCollection.hh"

#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Common/EDProductGetter.h"

namespace mu2e {

  // The first two arguments will be needed for the next revsion of the PointTrajectory class.
  template<typename SELECTOR>
  void compressPointTrajectoryCollection ( art::ProductID             const&   ,
                                           art::EDProductGetter       const*   ,
                                           PointTrajectoryCollection  const& in,
                                           SELECTOR                   const& keep,
                                           PointTrajectoryCollection&        out ){

    for ( PointTrajectoryCollection::const_iterator i=in.begin(), e=in.end(); i!=e; ++i ){
      if ( keep[i->first] ){
        PointTrajectory& sim = out[i->first];
        sim = i->second;
      }
    }
  }

}
#endif /* Mu2eUtilities_compressPointTrajectories_hh */
