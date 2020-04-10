#ifndef AnalysisConditions_MVACatalogMaker_hh
#define AnalysisConditions_MVACatalogMaker_hh

//
// Makes the MVACatalog ProditionsEntitiy
//

#include "AnalysisConditions/inc/MVACatalog.hh"
#include "AnalysisConfig/inc/MVACatalogConfig.hh"
#include "DbTables/inc/MVAToolDb.hh"

namespace mu2e {

  template <class T, class DB>
  class MVACatalogMaker {
  public:
    MVACatalogMaker(MVACatalogConfig const& config):_config(config) {}

    typename MVACatalog<T>::ptr_t fillEntries() {
      MVAEntries<T> mvaEntries;
      mvaEntries.reserve(_config.mvaConfigs().size()); // need this to avoid a segfault....
      for (const auto& i_entryConf : _config.mvaConfigs()) {
	mvaEntries.push_back(MVAEntry<T>(i_entryConf.trainName(), i_entryConf.xmlFileName(), i_entryConf.calibrated()));
      }

      auto ptr = std::make_shared<MVACatalog<T> >(mvaEntries);

      return ptr;
    }

    void initializeMVAs(typename MVACatalog<T>::ptr_t ptr) {
      // Now initialize the MVAs
      for (auto& i_mvaEntry : ptr->modifiableEntries()) {
	i_mvaEntry.initializeMVA();
      }
      if (_config.verbose() > 0) {
	ptr->print(std::cout);
      }
    }

    typename MVACatalog<T>::ptr_t fromFcl() {
      auto ptr = fillEntries();
      initializeMVAs(ptr);
      return ptr;
    }

    typename MVACatalog<T>::ptr_t fromDb(typename DB::cptr_t tqDb) {
      // fill the MVACatalog with initial values
      auto ptr = fillEntries();

      // Change or add entries to the configuration before we create all the MVATools pointers
      for (const auto& i_row : tqDb->rows()) {
	for (auto& i_mvaEntry : ptr->modifiableEntries()) {
	  if (i_row.mvaname() == i_mvaEntry._trainName) {
	    i_mvaEntry._xmlFileName = i_row.xmlfilename();
	    i_mvaEntry._calibrated = i_row.calibrated();
	    break; // break from inner for loop since this entry alread exists
	  }
	}
	// if we get here, then we need to add a new entry
	ptr->modifiableEntries().push_back(MVAEntry<T>(i_row.mvaname(), i_row.xmlfilename(), i_row.calibrated()));
      }
      initializeMVAs(ptr);
      return ptr;
    }

  private:
    // this object needs to be thread safe, 
    // _config should only be initialized once
    const MVACatalogConfig _config;
  };
}

#endif
