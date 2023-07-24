#include "Offline/CaloConditions/inc/CaloDAQMapMaker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "cetlib_except/exception.h"
#include <vector>
#include <fstream>

using namespace std;

namespace mu2e {

  CaloDAQMap::ptr_t CaloDAQMapMaker::fromFcl() {

    if (_config.verbose()>0) {
      cout << "CaloDAQMapMaker::fromFcl making nominal CaloDAQMap\n";
    }

    ConfigFileLookupPolicy configFile;
    string fileSpec = configFile(_config.fileSpec());
    if (_config.verbose()>0) {
      cout << "CaloDAQMapMaker::fromFcl reading from " << fileSpec << "\n";
    }

    ifstream ordFile;
    ordFile.open(fileSpec);
    if (!ordFile.is_open()) {
      throw cet::exception("CALODAQMAP_OPEN_FAILED")
        << " failed to open file " << fileSpec << "\n";
    }

    string line;

    CaloDAQMap::RawArray raw2Offline;
    CaloDAQMap::OfflineArray offline2Raw;
    uint16_t oid,rid,nRead(0);

    ordFile >> oid >> rid;
    while (!ordFile.eof()) {

      if(rid >= CaloConst::_nRawChannel) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CaloDAQMapMaker read invalid rawId" << rid << endl;
      }
      if(oid != CaloConst::_invalid && oid >= CaloConst::_nChannel) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CaloDAQMapMaker read invalid offlineId " << oid << endl;
      }

      nRead++;

      raw2Offline[rid] = CaloSiPMId(oid);

      if(oid < CaloConst::_nChannel) {
        offline2Raw[oid] = CaloRawSiPMId(rid);
      }

      ordFile >> oid >> rid;

    }

    if(nRead != CaloConst::_nRawChannel) {
      throw cet::exception("CALODAQMAPMAKER_COUNT")
        << "CaloDAQMapMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nRawChannel << endl;
    }

    auto ptr = make_shared<CaloDAQMap>(raw2Offline, offline2Raw);

    return ptr;

  } // end fromFcl


  CaloDAQMap::ptr_t CaloDAQMapMaker::fromDb() {

    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // swaps table added here when needed

    return ptr;

  } // end fromDb

}
