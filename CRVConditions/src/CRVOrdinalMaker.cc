#include "Offline/CRVConditions/inc/CRVOrdinalMaker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

namespace mu2e {

//***************************************************

CRVOrdinal::ptr_t CRVOrdinalMaker::fromFcl() {
  if (_config.verbose()) {
    cout << "CRVOrdinalMaker::fromFcl making nominal CRVOrdinal\n";
  }

  size_t nChan =
      GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size() *
      CRVId::nChanPerBar;

  // both maps initialized to invalid
  // both are sparse, so this is used to catch invalid access
  CRVOrdinal::OfflineMap offMap;
  for (size_t i = 0; i <= CRVId::nROC; i++) {             //ROC numbers start at 1
    for (size_t j = 0; j <= CRVId::nFEBPerROC; j++) {     //FEB numbers start at 1
      for (size_t k = 0; k < CRVId::nChanPerFEB; k++) {
        offMap[i][j][k] = CRVId::nChannels;
      }
    }
  }
  CRVOrdinal::OnlineMap onMap(nChan, CRVROC(0, 0, CRVId::nChanPerFEB));

  // load the online-to-offline numbering
  // chose the file based on the CRV geometry name
  std::string fileStub = _config.filePath() + "/" +
                         GeomHandle<CosmicRayShield>()->getName() + ".txt";

  ConfigFileLookupPolicy configFile;
  std::string fileSpec = configFile(fileStub);
  if (_config.verbose()) {
    cout << "CRVOrdinalMaker::fromFcl reading from " << fileSpec << "\n";
  }

  std::ifstream ordFile;
  ordFile.open(fileSpec);
  if (!ordFile.is_open()) {
    throw cet::exception("CRVORDINAL_OPEN_FAILED")
        << " failed to open file " << fileSpec << "\n";
  }

  std::string line;

  // read the header line
  std::getline(ordFile, line);

  size_t nRead = 0, maxChan = 0;
  std::vector<std::string> words;
  while (std::getline(ordFile, line)) {
    boost::split(words, line, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    if (words.size() != 4) {
      throw cet::exception("CRVORDINAL_BAD_FILE")
          << " failed to read line " << line << "\n";
    }
    std::uint16_t channel = std::stoul(words[0]);
    std::uint16_t ROC = std::stoul(words[1]);
    std::uint16_t FEB = std::stoul(words[2]);
    std::uint16_t FEBchannel = std::stoul(words[3]);
    onMap.at(channel) = CRVROC(ROC, FEB, FEBchannel);
    offMap.at(ROC).at(FEB).at(FEBchannel) = channel;
    nRead++;
    if (channel > maxChan) maxChan = channel;
  }

  if (_config.verbose()) {
    cout << "CRVOrdinalMaker::fromFcl channels read: " << nRead
         << "  max: " << maxChan << "  geom: "  << nChan << "\n";
  }

  auto ptr = make_shared<CRVOrdinal>(onMap, offMap);
  return ptr;

}  // end fromFcl

//***************************************************

CRVOrdinal::ptr_t CRVOrdinalMaker::fromDb() {
  if (_config.verbose()) {
    cout << "CRVOrdinalMaker::fromDb making CRVOrdinal\n";
  }

  // no database dependence (i.e. cable swaps) yet, so just return nominal
  return fromFcl();

}  // end fromDb

}  // namespace mu2e
