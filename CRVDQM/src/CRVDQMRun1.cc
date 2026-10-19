// Frozen Run 1 CRV layout lookups.
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVDQMRun1.hh"

namespace mu2e {
namespace CRVDQMRun1 {

const char* configurationName(int configuration)
{
  if (configuration < 0 || configuration >= kNConfigurations) {
    return "unknown";
  }
  return kLayouts[configuration].name;
}

int configurationFor(const std::string& crsName)
{
  for (const GeometryAlias& g : kGeometries) {
    if (crsName == g.crsName) {
      return g.configuration;
    }
  }
  return -1;
}

int sectorIndex(int configuration, const std::string& sector)
{
  if (configuration < 0 || configuration >= kNConfigurations) {
    return -1;
  }
  //CRSScintillatorShield::getName() is "CRV_" + crs.sectorNames
  const std::string prefix = "CRV_";
  const std::string bare =
      sector.compare(0, prefix.size(), prefix) == 0 ? sector.substr(prefix.size()) : sector;
  const Layout& l = kLayouts[configuration];
  for (int i = 0; i < l.nSectors; ++i) {
    if (bare == l.sectors[i]) {
      return i;
    }
  }
  return -1;
}

std::string sectorTag(int configuration, int sector)
{
  return std::string("_") + kLayouts[configuration].name + "_" +
         kLayouts[configuration].sectors[sector];
}

} // namespace CRVDQMRun1
} // namespace mu2e
