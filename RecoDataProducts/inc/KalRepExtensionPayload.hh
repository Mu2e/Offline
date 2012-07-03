#ifndef RecoDataProducts_KalRepExtensionPayload_hh
#define RecoDataProducts_KalRepExtensionPayload_hh
//
// The extended persistent information from a KalRep - see note 1.
//
// $Id: KalRepExtensionPayload.hh,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Contact person Rob Kutschke
//
// Notes:
//
// 1) The use pattern for this class is as follows. The track fit module will write to the
//    event both a collection of transient tracks and a KalRepPayloadCollection; the former
//    is highly functional but not persistable while the latter is a data-only object that
//    can be used to reconstitute the transient tracks in a later job.  Subsequent modules
//    in the same job may ask the transient tracks to extend themselves; the extension
//    operation is not permitted to modify the information that was stored in the
//    KalRepPayloadCollection.  The purpose of this class is the hold the persistent information
//    associated with any trajectory information created by the extension process.
//
// 2) The mutators for fowardExtensionPayloads_ and reverseExtensionPayloads_ are used at construction
//    time so that the two data members can be built in place, thereby minimizing copies.
//    Come C++11 we can use move-aware arguments to initialize these data members.
//

// C++ includes
#include <iosfwd>
#include <vector>

#include "RecoDataProducts/inc/LocalHelixPayload.hh"

namespace mu2e {

  class KalRepExtensionPayload{

  public:

    // Accessors
    std::vector<LocalHelixPayload> const & forwardExtensions() const { return forwardExtensions_; }
    std::vector<LocalHelixPayload> const & reverseExtensions() const { return reverseExtensions_; }

    // Constructors - see note 2
    KalRepExtensionPayload():
      forwardExtensions_(),
      reverseExtensions_(){
    }
    // Accept compiler written d'tor, copy c'tor and copy assignment.


    // Mutators - See note 2.
    std::vector<LocalHelixPayload>& forwardExtensions(){ return forwardExtensions_; }
    std::vector<LocalHelixPayload>& reverseExtensions(){ return reverseExtensions_; }

    // Print contents of the object.
    void print( std::ostream& ost, bool doEndl = true ) const;

  private:

    std::vector<LocalHelixPayload> forwardExtensions_;
    std::vector<LocalHelixPayload> reverseExtensions_;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   KalRepExtensionPayload const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* RecoDataProducts_KalRepExtensionPayload_hh */
