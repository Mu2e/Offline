#ifndef ToyDP_DPIndex_hh
#define ToyDP_DPIndex_hh

//
// A persistable index into another data product.
//
//
// $Id: DPIndex.hh,v 1.7 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//

// Framework includes.
#include "art/Persistency/Provenance/ProductID.h"

namespace mu2e {
  struct DPIndex{

    // The actual data for this struct.
    art::ProductID id;
    uint32_t       index;

    DPIndex():
      id(),
      index(0){}

    DPIndex( art::ProductID const& id_, int index_):
      id(id_),
      index(index_){
    } 

    // Compiler generated versions are OK for:
    // destructor, copy c'tor, assignment operator.

  };

  inline bool operator==(const DPIndex& lhs, 
                         const DPIndex& rhs){
      return ( lhs.id == rhs.id && lhs.index == rhs.index );
  }

  inline bool operator!=(const DPIndex& lhs, 
                         const DPIndex& rhs){
    return !(lhs==rhs);
  }

  // Sort first on ProductID and then on index.
  inline bool operator<(const DPIndex& lhs, 
                        const DPIndex& rhs){
    return ( lhs.id < rhs.id ) || 
      ( lhs.id == rhs.id && lhs.index < rhs.index );
  }

  // ProductID does not define operators >, <=, >= so we would need 
  // to fix that before defining those operators for this class.


  inline std::ostream& operator<<( std::ostream& ost,
                                   DPIndex const& dpi ){
    ost << "(" << dpi.id 
        << ","  << dpi.index 
        << ")";
    return ost;
  }


}

#endif /* ToyDP_DPIndex_hh */
