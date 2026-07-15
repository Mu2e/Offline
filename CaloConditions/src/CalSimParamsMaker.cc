#include "Offline/CaloConditions/inc/CalSimParamsMaker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "cetlib_except/exception.h"
#include <vector>
#include <fstream>


namespace mu2e {

  CalSimParams::ptr_t CalSimParamsMaker::fromFcl() {

    CalSimParams::IdArray     crystalIds;
    CalSimParams::floatArray  LRUs;
    CalSimParams::floatVArray pePerMeVs;
    CalSimParams::floatVArray ADCPerMeVs;

    if (_config.verbose()>0) {
      std::cout << "CalSimParamsMaker::fromFcl making nominal CalSimCrystals\n";
    }

    ConfigFileLookupPolicy configFile;
    std::string fileName = configFile(_config.fileName());
    if (_config.verbose()>0) {
      std::cout << "CalSimParamsMaker::fromFcl reading from " << fileName << "\n";
    }

    std::ifstream ordFile(fileName);
    if (!ordFile.is_open()) {
      throw cet::exception("CALSIMPARAMSMAKER_OPEN_FAILED")
        << " failed to open file " << fileName << "\n";
    }

    uint16_t nRead(0);
    std::string line;
    while (std::getline(ordFile, line)) {
      uint16_t sid;
      float lru, npe0, npe1, adc0, adc1;

      std::istringstream iss(line);
      if (!(iss >> sid >> lru >> npe0 >> npe1 >> adc0 >> adc1)) {
        throw cet::exception("CALSIMPARAMSMAKER_RANGE")
        << "invalid format at line "<<nRead+1<<"\n";
      }

      // Check that there is nothing left on the line
      float extra;
      if (iss >> extra) {
        throw cet::exception("CALSIMPARAMSMAKER_RANGE")
        << "invalid format at line "<<nRead+1<<"\n";
      }

      //the file should contain crystals in ascending order
      if (sid != nRead){
        throw cet::exception("CALSIMPARAMSMAKER_RANGE")
        << "invalid sid values at line "<<nRead+1<<"\n";
      }

      //and the values should be positive
      if (lru<0 || npe0<0 || npe1<0){
        throw cet::exception("CALSIMPARAMSMAKER_RANGE")
        << "invalid lru, npe, or npe1 values at line "<<nRead+1<<"\n";
      }

      crystalIds[sid]  = CrystalId(sid);
      LRUs[sid]        = lru;
      pePerMeVs[sid]   = std::vector<float>{npe0,npe1};
      ADCPerMeVs[sid]  = std::vector<float>{adc0,adc1};
      ++nRead;
    }

    if (nRead != CaloConst::_nCrystal) {
      throw cet::exception("CALSIMPARAMSMAKER_COUNT")
        << "CalSimParamsMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nCrystal << "\n";
    }

    auto ptr = make_shared<CalSimParams>(crystalIds,LRUs,pePerMeVs,ADCPerMeVs);
    return ptr;
  }


  //***************************************************
  CalSimParams::ptr_t CalSimParamsMaker::fromDb(CalSimCrystals::cptr_t cch_p) {

    if (_config.verbose()>0) {
      std::cout << "CalSimParamsMaker::fromDb making CalSimCrystals\n";
    }

    CalSimParams::IdArray     crystalIds;
    CalSimParams::floatArray  LRUs;
    CalSimParams::floatVArray pePerMeVs;
    CalSimParams::floatVArray ADCPerMeVs;
    uint16_t nRead(0);

    for (auto const& row : cch_p->rows()) {
      CrystalId  crid = row.crid();
      float      lru  = row.LRU();
      float      npe0 = row.pePerMeV_0();
      float      npe1 = row.pePerMeV_1();
      float      adc0 = row.ADCPerMeV_0();
      float      adc1 = row.ADCPerMeV_1();

      //the file should contain crystals in ascending order
      if (crid.id() != nRead || !crid.isValid()){
        throw cet::exception("CALSIMPARAMSMAKER_RANGE")
        << "invalid crystalID values at line "<<nRead+1<<"\n";
      }

      crystalIds[crid.id()]  = crid;
      LRUs[crid.id()]        = lru;
      pePerMeVs[crid.id()]   = std::vector<float>{npe0,npe1};
      ADCPerMeVs[crid.id()]  = std::vector<float>{adc0,adc1};
      ++nRead;
    }


    // check that all roid were filled
    if (nRead != CaloConst::_nCrystal) {
      throw cet::exception("CALSIMPARAMSMAKER_COUNT")
        << "CalSimParamsMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nCrystal << "\n";
    }

    auto ptr = make_shared<CalSimParams>(crystalIds,LRUs,pePerMeVs,ADCPerMeVs);
    return ptr;
  }
}
