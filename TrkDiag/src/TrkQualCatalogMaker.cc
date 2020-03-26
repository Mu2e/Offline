#include "TrkDiag/inc/TrkQualCatalogMaker.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

namespace mu2e {

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fillEntries() {
    TrkQualEntries trkQualEntries;
    for (const auto& i_entryConf : _config.trkQualConfigs()) {
      trkQualEntries.push_back(TrkQualEntry(i_entryConf.trainName(), i_entryConf.xmlFileName(), i_entryConf.calibrated()));
    }

    auto ptr = std::make_shared<TrkQualCatalog>(trkQualEntries);//_config.deadTimeAnalog(), 

    return ptr;
  }

  void TrkQualCatalogMaker::initializeMVAs(TrkQualCatalog::ptr_t ptr) {
    // Now initialize the MVAs
    for (auto& i_trkQualEntry : ptr->modifiableEntries()) {
      i_trkQualEntry.initializeMVA();
    }

    if (_config.verbose() > 0) {
      ptr->print(std::cout);
    }
  }

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fromFcl() {

    auto ptr = fillEntries();
    initializeMVAs(ptr);
    return ptr;
  }

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fromDb(TrkQualDb::cptr_t tqDb) {

    // fill the TrkQualCatalog with initial values
    auto ptr = fillEntries();

    // Change or add entries to the configuration before we create all the MVATools pointers
    for (const auto& i_row : tqDb->rows()) {
      for (auto& i_trkQualEntry : ptr->modifiableEntries()) {
	if (i_row.mvaname() == i_trkQualEntry._trainName) {
	  i_trkQualEntry._xmlFileName = i_row.xmlfilename();
	  i_trkQualEntry._calibrated = i_row.calibrated();
	  break; // break from inner for loop since this entry alread exists
	}
      }
      // if we get here, then we need to add a new entry
      ptr->modifiableEntries().push_back(TrkQualEntry(i_row.mvaname(), i_row.xmlfilename(), i_row.calibrated()));
    }

    initializeMVAs(ptr);
    return ptr;
  }

};
