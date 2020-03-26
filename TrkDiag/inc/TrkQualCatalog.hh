#ifndef TrkDiag_TrkQualCatalog_hh
#define TrkDiag_TrkQualCatalog_hh

//
// TrkQualCatalog is the ProditionsEntry for TrkQualDb
//

// Mu2e includes
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  struct TrkQualEntry {
    TrkQualEntry(std::string trainName, std::string xmlFileName) : _trainName(trainName), _mvaTool(new MVATools(xmlFileName)) {}
    
    std::string _trainName;
    MVATools* _mvaTool;
  };
  typedef std::vector<TrkQualEntry> TrkQualEntries;

  class TrkQualCatalog : virtual public ProditionsEntity {

  public:
    typedef std::shared_ptr<TrkQualCatalog> ptr_t;
    typedef std::shared_ptr<const TrkQualCatalog> cptr_t;

    TrkQualCatalog() : _name("TrkQualCatalog") {}
    TrkQualCatalog(TrkQualEntries entries) : _name("TrkQualCatalog"), _entries(entries) {}

    // accessors
    void print(std::ostream& os) const;
    std::string const& name() const { return _name; }
    TrkQualEntries& modifiableEntries() { return _entries;} 
    TrkQualEntries const& entries() const { return _entries;} 

  private:
    // data
    std::string _name;
    TrkQualEntries _entries;
  };
}

#endif
