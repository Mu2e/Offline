#ifndef RecoDataProducts_KalRepPayload_hh
#define RecoDataProducts_KalRepPayload_hh
//
// The persistent information from a KalRep
//
// $Id: KalRepPayload.hh,v 1.3 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
// Contact person Rob Kutschke
//
// Notes:
//
// 1) This class is made from the data within a transient KalRep.
//    On construction the member transientTrack_ points to the track from which
//    it was made.  On readback, transientTrack_ is default constructed: that is, it
//    does not point to anything.  When the transient tracks are reconstructed from
//    this data product, the transientTrack data member should be set to point to
//    reincarnated transient track.  When we write the class that reinstaniates the
//    transient track, we will make the mutator private and make that class a friend.
//
// 2) The mutators for trajectory_ and hots_ are used at construction time so that
//    the two data members can be built in place, thereby minimizing copies.
//    Come C++11 we can use move-aware arguments to initialize these data members.

// C++ includes
#include <iosfwd>
#include <vector>

#include "RecoDataProducts/inc/LocalHelixPayload.hh"
#include "RecoDataProducts/inc/HOTPayload.hh"

#include "canvas/Persistency/Common/Ptr.h"

  class KalRep;

namespace mu2e {

  class KalRepPayload{
  public:

    // Accessors
    LocalHelixPayload const& seed() const { return seed_; }

    bool success() const { return success_; }
    bool nHit()    const { return nHit_;    }
    bool nActive() const { return nActive_; }
    bool nDof()    const { return nDof_;    }
    bool cl()      const { return cl_;      }

    std::vector<LocalHelixPayload> const& trajectory() const { return trajectory_; }
    std::vector<HOTPayload> const& hots()              const { return hots_; }

    art::Ptr< const KalRep * const  > const& transientTrack () const { return transientTrack_; }

    // Constructors
    KalRepPayload():
      seed_(),
      success_(false),
      nHit_(0),
      nActive_(0),
      nDof_(0),
      cl_(0.),
      trajectory_(),
      hots_(),
      transientTrack_(){
    }
    // Accept compiler written d'tor, copy c'tor and copy assignment.

    KalRepPayload( LocalHelixPayload const& aseed,
                   bool asuccess,
                   bool anHit,
                   bool anActive,
                   bool anDof,
                   bool acl,
                   art::Ptr< const KalRep * const  > const& atransientTrack ):
      seed_(aseed),
      success_(asuccess),
      nHit_(anHit),
      nActive_(anActive),
      nDof_(anDof),
      cl_(acl),
      trajectory_(),
      hots_(),
      transientTrack_(atransientTrack) {

      // See note 2.
      trajectory_.reserve(nHit_);
      hots_.reserve(nHit_);

    }

    // Mutators

    // See note 2.
    std::vector<LocalHelixPayload>& trajectory() { return trajectory_; }
    std::vector<HOTPayload>& hots()       { return hots_; }

    // See note 1.
    art::Ptr< const KalRep * const>&
    transientTrack ( art::Ptr<const KalRep * const> const& p ) const {
      return transientTrack_=p;
    }

    // Print contents of the object.
    void print( std::ostream& ost, bool doEndl = true ) const;

  private:

    LocalHelixPayload seed_;

    bool success_;
    bool nHit_;
    bool nActive_;
    bool nDof_;
    bool cl_;

    std::vector<LocalHelixPayload> trajectory_;
    std::vector<HOTPayload>        hots_;

    // See note 1.
    mutable art::Ptr< const KalRep * const  > transientTrack_;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   KalRepPayload const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* RecoDataProducts_KalRepPayload_hh */
