//
// Convenience class to access provenance information about a KalRepCollection.
//
// Contact person, Rob Kutschke
//

#include "Mu2eUtilities/inc/KalRepCollectionInfo.hh"

namespace {

  // Helper function to access a handle; used to enable delegation of c'tor.
  art::Handle<mu2e::KalRepCollection> handleGetter(art::ProductID id, art::Event const& event ){
    art::Handle<mu2e::KalRepCollection> handle;
    event.get(id, handle);
    return handle;
  }

}

namespace mu2e {

  KalRepCollectionInfo::KalRepCollectionInfo( art::Handle<KalRepCollection> const& handle ):
    patRecType_(handle.provenance()->parameterSet().get<std::string>("module_type")),
    instance_(handle.provenance()->productInstanceName()){
  }

  KalRepCollectionInfo::KalRepCollectionInfo( art::ValidHandle<KalRepCollection> const& handle ):
    patRecType_(handle.provenance()->parameterSet().get<std::string>("module_type")),
    instance_(handle.provenance()->productInstanceName()){
  }

  // Other c'tors delegate.
  KalRepCollectionInfo::KalRepCollectionInfo( KalRepPtr const& ptr, art::Event const& event ):
    KalRepCollectionInfo( ptr.id(), event){
  }

  KalRepCollectionInfo::KalRepCollectionInfo( art::ProductID const& id, art::Event const& event ):
    KalRepCollectionInfo( handleGetter(id,event) ){
  }

}
