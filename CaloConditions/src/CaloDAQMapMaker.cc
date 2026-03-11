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

    while (!ordFile.eof()){
      ordFile >> rid >> oid;
      if (ordFile.eof()) break;

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

    }

    if(nRead != CaloConst::_nRawChannel) {
      throw cet::exception("CALODAQMAPMAKER_COUNT")
        << "CaloDAQMapMaker read the wrong number of id's "
        << nRead << ", expected " << CaloConst::_nRawChannel << endl;
    }

    auto ptr = make_shared<CaloDAQMap>(raw2Offline, offline2Raw);

    return ptr;

  } // end fromFcl

  //***************************************************

  CaloDAQMap::ptr_t CaloDAQMapMaker::fromDb(CalChannels::cptr_t cch_p) {

    if (_config.verbose()>0) {
      cout << "CaloDAQMapMaker::fromDb making CaloDAQMap\n";
    }

    CaloDAQMap::RawArray raw2Offline;
    CaloDAQMap::OfflineArray offline2Raw;

    for (auto const& row : cch_p->rows()) {
      CaloRawSiPMId rawid = row.rawid();
      CaloSiPMId roid = row.roid();
      if(!rawid.isValid()) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CaloDAQMapMaker found invalid rawId" << rawid << endl;
      }
      if(!(roid.isValid() || roid.id() == CaloConst::_invalid)) {
        throw cet::exception("CALODAQMAPMAKER_RANGE") << "CaloDAQMapMaker found invalid offlineId " << roid << endl;
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
        throw cet::exception("CALODAQMAPMAKER_MISSING") << "CaloDAQMapMaker found missing roid " << endl;
      }
    }

    auto ptr = make_shared<CaloDAQMap>(raw2Offline, offline2Raw);

    return ptr;

  } // end fromDb

}
