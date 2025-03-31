#include "Offline/DbService/inc/GrlList.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace mu2e;

//**************************************************
GrlList::GrlList(const GrlHeader& header) : _header(header) {}

//**************************************************
GrlList::GrlList(const GrlHeader& header, const IoVVec& grl) :
    _header(header), _grl(grl) {}

//**************************************************
GrlList::GrlList(const GrlHeader& header, const std::string& filename) :
    _header(header) {
  if (filename.size() <= 0) {
    throw cet::exception("GRLLIST_NO_FILE_NAME")
        << "DbUtil::read called with no file name\n";
  }
  std::ifstream myfile;
  myfile.open(filename);
  if (!myfile.is_open()) {
    throw cet::exception("GRLLIST_OPEN_FAILED")
        << "GrlList failed to open " << filename << "\n";
  }

  std::string line;
  while (std::getline(myfile, line)) {
    boost::trim(line);  // remove whitespace
    if (line.size() <= 0 || line[0] == '#') continue;
    _grl.emplace_back(line);
  }

  return;
}

//**************************************************
bool GrlList::goodRun(uint32_t run) {
  // start with the whole run
  DbIoV riov(run, 0, run, DbIoV::maxSubrun());
  for (auto& iov : _grl) {
    // subtract good intervals
    riov.subtract(iov);
    // if nothing left, then all good
    if (riov.isNull()) return true;
  }
  return false;
}
//**************************************************
bool GrlList::goodSubRun(uint32_t run, uint32_t subrun) {
  for (auto& iov : _grl) {
    if (iov.inInterval(run, subrun)) return true;
  }
  return false;
}

//**************************************************
void GrlList::print(std::ostream& os) {
  for (auto& iov : _grl) {
    os << iov.to_string() << "\n";
  }
  return;
}
