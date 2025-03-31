#ifndef DbService_GrlList_hh
#define DbService_GrlList_hh

// holds a goodrun list as a vector of intervals of validity

#include "Offline/DbService/inc/GrlHeader.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

class GrlList {
 public:
  typedef std::vector<DbIoV> IoVVec;

  explicit GrlList(const GrlHeader& header);
  GrlList(const GrlHeader& header, const IoVVec& grl);
  // text file containing one IoV per line (see wiki for format)
  // lines may be blank or begin with "#' for comments
  GrlList(const GrlHeader& header, const std::string& filename);
  GrlList(const GrlList& other) = default;

  bool goodRun(uint32_t run);  // if part bad, this will return false
  bool goodSubRun(uint32_t run, uint32_t subrun);
  // int goodSubRun(const art::SubRun& subrun);

  void print(std::ostream& os = std::cout);

  const IoVVec& grl() const { return _grl; }
  const GrlHeader& header() const { return _header; }

  GrlHeader _header;

 private:
  IoVVec _grl;
};

}  // namespace mu2e

#endif
