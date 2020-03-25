#ifndef TrkDiag_TrkQualCatalogConfig_hh
#define TrkDiag_TrkQualCatalogConfig_hh

//
// Configuration for the TrkQualCatalog ProditionsEntitiy
//

#include "fhiclcpp/types/Atom.h"

namespace mu2e {
  struct TrkQualCatalogConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    //    data initialization values ....
  };
}

#endif
