#ifndef TrkDiag_TrkQualCatalog_hh
#define TrkDiag_TrkQualCatalog_hh

//
// TrkQualCatalog is the ProditionsEntry for TrkQualDb
//

// Mu2e includes
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class TrkQualCatalog : virtual public ProditionsEntity {
  public:
    typedef std::shared_ptr<TrkQualCatalog> ptr_t;
    typedef std::shared_ptr<const TrkQualCatalog> cptr_t;

    TrkQualCatalog() : _name("TrkQualCatalog") {}

    // accessors
    void print(std::ostream& os) const;
    std::string const& name() const { return _name; }

  private:
    // data
    std::string _name;
  };
}

#endif
