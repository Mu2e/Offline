#include "Offline/CaloConditions/inc/CalSimCrystalMaker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "cetlib_except/exception.h"
#include <vector>
#include <fstream>


namespace mu2e {

  CalSimCrystal::ptr_t CalSimCrystalMaker::fromFcl() {

    if (_config.verbose()>0) {
      std::cout << "CalSimCrystalMaker::fromFcl making nominal CalCrystals\n";
    }

    ConfigFileLookupPolicy configFile;
    std::string fileName = configFile(_config.fileName());
    if (_config.verbose()>0) {
      std::cout << "CalSimCrystalMaker::fromFcl reading from " << fileName << "\n";
    }

    std::ifstream ordFile;
    ordFile.open(fileName);
    if (!ordFile.is_open()) {
      throw cet::exception("CalCrystals_OPEN_FAILED")
        << " failed to open file " << fileName << "\n";
    }

    CalSimCrystal::IdArray     crystalIds;
    CalSimCrystal::floatArray  LRUs;
    CalSimCrystal::floatVArray pePerMeVs;
    uint16_t sid,nRead(0);
    float lru,npe0,npe1;

    while (!ordFile.eof()){
      ordFile >> sid >> lru >> npe0 >> npe1;
      if (ordFile.eof()) break;

      if (sid != CaloConst::_invalid && sid >= CaloConst::_nCrystal) {
        throw cet::exception("CalSimCrystalMaker_RANGE") << "CalSimCrystalMaker read invalid offlineId " << sid << "\n";
      }

      crystalIds[sid]  = CrystalId(sid);
      LRUs[sid]        = lru;
      pePerMeVs[sid]   = std::vector<float>{npe0,npe1};
      ++nRead;
    }

    if (pePerMeVs.back().size() !=CaloConst::_nSiPMPerCrystal ) {
      throw cet::exception("CALSIMCRYSTALMAKER_RANGE")
      << "CalSimCrystalMaker found  too many nPE values\n";
    }

    if(nRead != CaloConst::_nCrystal) {
      throw cet::exception("CalSimCrystalMaker_COUNT")
        << "CalSimCrystalMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nCrystal << "\n";
    }

    auto ptr = make_shared<CalSimCrystal>(crystalIds,LRUs,pePerMeVs);
    return ptr;
  }


  //***************************************************
  CalSimCrystal::ptr_t CalSimCrystalMaker::fromDb(CalCrystals::cptr_t cch_p) {

    if (_config.verbose()>0) {
      std::cout << "CalSimCrystalMaker::fromDb making CalCrystals\n";
    }

    CalSimCrystal::IdArray     crystalIds;
    CalSimCrystal::floatArray  LRUs;
    CalSimCrystal::floatVArray pePerMeVs;
    uint16_t nRead(0);

    for (auto const& row : cch_p->rows()) {
      CrystalId  crid = row.crid();
      float      lru  = row.LRU();
      float      npe0 = row.pePerMeV_0();
      float      npe1 = row.pePerMeV_1();

      if (!(crid.isValid() || crid.id() == CaloConst::_invalid)) {
        throw cet::exception("CALSIMCRYSTALMAKER_RANGE")
        << "CalSimCrystalMaker found invalid offlineId " << crid << "\n";
      }

      crystalIds[crid.id()]  = crid;
      LRUs[crid.id()]        = lru;
      pePerMeVs[crid.id()]   = std::vector<float>{npe0,npe1};
      ++nRead;
    }

    if (pePerMeVs.back().size() !=CaloConst::_nSiPMPerCrystal ) {
      throw cet::exception("CALSIMCRYSTALMAKER_RANGE")
      << "CalSimCrystalMaker found  too many nPE values\n";
    }

    // check that all roid were filled
    if(nRead != CaloConst::_nCrystal) {
      throw cet::exception("CalSimCrystalMaker_COUNT")
        << "CalSimCrystalMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nCrystal << "\n";
    }

    auto ptr = make_shared<CalSimCrystal>(crystalIds,LRUs,pePerMeVs);
    return ptr;
  }
}
