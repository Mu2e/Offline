#include "Offline/CRVConditions/inc/CRVScintYieldMaker.hh"
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

CRVScintYield::ptr_t CRVScintYieldMaker::fromFcl() {
  if (_config.verbose()) {
    cout << "CRVScintYieldMaker::fromFcl making nominal CRVScintYield\n";
  }

  size_t nBarIndices =
      GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size();

  CRVScintYield::ScintYieldMap scintYieldMap(nBarIndices,0);  // initialized to 0 (no deviation from the nominal scintillation yield)

  // load the scintillation spread
  std::string fileStub = _config.fileName();

  ConfigFileLookupPolicy configFile;
  std::string fileSpec = configFile(fileStub);
  if (_config.verbose()) {
    cout << "CRVScintYieldMaker::fromFcl reading from " << fileSpec << "\n";
  }

  std::ifstream scintYieldFile;
  scintYieldFile.open(fileSpec);
  if (!scintYieldFile.is_open()) {
    throw cet::exception("CRVSCINTYIELD_OPEN_FAILED")
        << " failed to open file " << fileSpec << "\n";
  }

  std::string line;

  // read the header line
  std::getline(scintYieldFile, line);

  size_t nRead = 0, maxBars = 0;
  std::vector<std::string> words;
  while (std::getline(scintYieldFile, line)) {
    boost::split(words, line, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    if (words.size() != 2) {
      throw cet::exception("CRVSCINTYIELD_BAD_FILE")
          << " failed to read line " << line << "\n";
    }
    std::uint16_t barIndex = std::stoul(words[0]);
    if (barIndex >= scintYieldMap.size()) {
      throw cet::exception("CRVSCINTYIELD_BAD_SCINTILLATOR BAR INDEX")
          << "CRVScintYieldMaker::fromFcl read barIndex in file that doesn't exist in geometry: "
          << " barIndex=" << barIndex << "\n";
    }
    float scintYieldDeviation = std::stof(words[1]);
    scintYieldMap.at(barIndex) = scintYieldDeviation;
    nRead++;
    if (barIndex > maxBars) maxBars = barIndex;
  }

  if (_config.verbose()) {
    cout << "CRVScintYieldMaker::fromFcl bar indices read: " << nRead
         << "  max: " << maxBars << "  geom: "  << nBarIndices << "\n";
  }

  auto ptr = make_shared<CRVScintYield>(scintYieldMap);
  return ptr;

}  // end fromFcl

//***************************************************

CRVScintYield::ptr_t CRVScintYieldMaker::fromDb() {
  if (_config.verbose()) {
    cout << "CRVScintyieldMaker::fromDb making CRVScintYield\n";
  }

  // no database dependence yet, so just return nominal
  return fromFcl();

}  // end fromDb

}  // namespace mu2e
