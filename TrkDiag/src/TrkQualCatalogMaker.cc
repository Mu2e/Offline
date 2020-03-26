#include "TrkDiag/inc/TrkQualCatalogMaker.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

namespace mu2e {

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fromFcl() {
    TrkQualEntries trkQualEntries;
    for (const auto& i_entryConf : _config.trkQualConfigs()) {
      trkQualEntries.push_back(TrkQualEntry(i_entryConf.trainName(), i_entryConf.xmlFileName()));
    }

    auto ptr = std::make_shared<TrkQualCatalog>(trkQualEntries);//_config.deadTimeAnalog(), 
    return ptr;
  }

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fromDb(TrkQualDb::cptr_t tqDb) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // now overwrite with db values
    for (auto& i_trkQualEntry : ptr->modifiableEntries()) {
      for (const auto& i_row : tqDb->rows()) {
	if (i_row.mvaname() == i_trkQualEntry._trainName) { // check the training name
	  i_trkQualEntry._mvaTool = new MVATools(i_row.xmlfilename());
	}
      }
    }

    return ptr;
  }

};
