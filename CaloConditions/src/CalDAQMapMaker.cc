#include "Offline/CaloConditions/inc/CalDAQMapMaker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "cetlib_except/exception.h"
#include <vector>
#include <fstream>

using namespace std;

namespace mu2e {

  CalDAQMap::ptr_t CalDAQMapMaker::fromFcl() {

    CalDAQMap::RawArray raw2Offline;
    CalDAQMap::OfflineArray offline2Raw;

    if (_config.verbose()>0) {
      cout << "CalDAQMapMaker::fromFcl making nominal CalDAQMap\n";
    }

    ConfigFileLookupPolicy configFile;
    string fileSpec = configFile(_config.fileSpec());
    if (_config.verbose()>0) {
      cout << "CalDAQMapMaker::fromFcl reading from " << fileSpec << "\n";
    }

    ifstream ordFile(fileSpec);
    if (!ordFile.is_open()) {
      throw cet::exception("CALODAQMAP_OPEN_FAILED")
        << " failed to open file " << fileSpec << "\n";
    }

    uint16_t nRead(0);
    std::string line;
    while (std::getline(ordFile, line)) {
      uint16_t oid,rid;

      std::istringstream iss(line);
      if (!(iss >> rid >> oid)) {
        throw cet::exception("CALODAQMAPMAKER_RANGE")
        << "file format invalid on line "<<nRead+1<<"\n";
      }

      // Check that there is nothing left on the line
      float extra;
      if (iss >> extra) {
        throw cet::exception("CALODAQMAPMAKER_RANGE")
        << "file format invalid on line "<<nRead+1<<"\n";
      }

      if(rid >= CaloConst::_nRawChannel) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CalDAQMapMaker read invalid rawId" << rid << endl;
      }

      if(oid >= CaloConst::_nChannel) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CalDAQMapMaker read invalid offlineId " << oid << endl;
      }

      raw2Offline[rid] = CaloSiPMId(oid);
      offline2Raw[oid] = CaloRawSiPMId(rid);
      ++nRead;
    }

    if(nRead != CaloConst::_nRawChannel) {
      throw cet::exception("CALODAQMAPMAKER_COUNT")
        << "CalDAQMapMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nRawChannel << endl;
    }

    auto ptr = make_shared<CalDAQMap>(raw2Offline, offline2Raw);

    return ptr;

  } // end fromFcl

  //***************************************************

  CalDAQMap::ptr_t CalDAQMapMaker::fromDb(CalChannels::cptr_t cch_p) {

    if (_config.verbose()>0) {
      cout << "CalDAQMapMaker::fromDb making CalDAQMap\n";
    }

    CalDAQMap::RawArray raw2Offline;
    CalDAQMap::OfflineArray offline2Raw;

    for (auto const& row : cch_p->rows()) {
      CaloRawSiPMId rawid = row.rawid();
      CaloSiPMId roid = row.roid();
      if(!rawid.isValid()) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CalDAQMapMaker found invalid rawId" << rawid << endl;
      }
      if(!(roid.isValid() || roid.id() == CaloConst::_invalid)) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CalDAQMapMaker found invalid offlineId " << roid << endl;
      }

      raw2Offline[rawid.id()] = roid;

      if(roid.isValid()) {
        offline2Raw[roid.id()] = rawid;
      }
    } // end loop over table

    // check that all roid were filled since some table entries are set
    // to invalid and there is no other check that this was done correctly
    for (auto const& cc : offline2Raw) {
      if(!cc.isValid()) {
        throw cet::exception("CALODAQMAPMAKER_MISSING") << "CalDAQMapMaker found missing roid " << endl;
      }
    }

    auto ptr = make_shared<CalDAQMap>(raw2Offline, offline2Raw);

    return ptr;

  } // end fromDb

}
