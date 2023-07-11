#include "Offline/DbTables/inc/DbIoV.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <iostream>

void mu2e::DbIoV::subtract(DbIoV const& iov, uint32_t run, uint32_t subrun) {
  // check for no overlap
  if (iov.endRun() < startRun() ||
      (iov.endRun() == startRun() && iov.endSubrun() < startSubrun()))
    return;
  if (iov.startRun() > endRun() ||
      (iov.startRun() == endRun() && iov.startSubrun() > endSubrun()))
    return;

  // check if piece will be at begining or end
  bool sin =
      ((iov.startRun() > startRun()) ||
       (iov.startRun() == startRun() && iov.startSubrun() > startSubrun()));
  bool ein = ((iov.endRun() < endRun()) ||
              (iov.endRun() == endRun() && iov.endSubrun() < endSubrun()));

  int method = 0;
  if (!sin && !ein) {  // subtract all
    method = 0;
  } else if (sin && !ein) {  // subtract middle to end
    method = 1;
  } else if (ein && !sin) {  // subtract start to middle
    method = 2;
  } else {
    // subtract a piece in the middle, keep the early or
    // late fragment based on the run/subrun
    method = 1;
    if (run > iov.startRun() ||
        (run == iov.startRun() && subrun >= iov.startSubrun()))
      method = 2;
  }

  if (method == 0) {  // subtract the whole interval
    set(0, 0, 0, 0);
  } else if (method == 1) {  // subtract middle to end
    uint32_t er = iov.startRun();
    uint32_t es = iov.startSubrun();
    if (es == 0) {
      er--;
      es = maxSubrun();
    } else {
      es--;
    }
    set(startRun(), startSubrun(), er, es);
    return;
  } else if (method == 2) {  // subtract start to middle
    uint32_t sr = iov.endRun();
    uint32_t ss = iov.endSubrun();
    if (ss == maxSubrun()) {
      sr++;
      ss = 0;
    } else {
      ss++;
    }
    set(sr, ss, endRun(), endSubrun());
  }

  return;
}

void mu2e::DbIoV::overlap(DbIoV const& iov) {
  // check for no overlap
  bool null = false;
  if (iov.endRun() < startRun() ||
      (iov.endRun() == startRun() && iov.endSubrun() < startSubrun()))
    null = true;
  if (iov.startRun() > endRun() ||
      (iov.startRun() == endRun() && iov.startSubrun() > endSubrun()))
    null = true;
  if (null) {
    set(0, 0, 0, 0);
    return;
  }

  // check if iov start/end are inside this interval
  bool sin = inInterval(iov.startRun(), iov.startSubrun());
  bool ein = inInterval(iov.endRun(), iov.endSubrun());

  if (sin && ein) {  // overlap is iov
    set(iov.startRun(), iov.startSubrun(), iov.endRun(), iov.endSubrun());
  } else if (sin && !ein) {  // overlap middle to end
    set(iov.startRun(), iov.startSubrun(), endRun(), endSubrun());
  } else if (ein && !sin) {  // overlap start to middle
    set(startRun(), startSubrun(), iov.endRun(), iov.endSubrun());
  }
  // last case (!sin && !ein) is no change
  return;
}

int mu2e::DbIoV::isOverlapping(DbIoV const& iov) const {
  // these check that there is non-zero overlap
  bool sor = (iov.startRun() < endRun() ||
              (iov.startRun() == endRun() && iov.startSubrun() <= endSubrun()));
  bool eor = (iov.endRun() > startRun() ||
              (iov.endRun() == startRun() && iov.endSubrun() >= startSubrun()));

  if ((!eor) || (!sor)) return 0;  // no overlap

  // these check that the overlap is complete
  bool scp =
      (iov.startRun() < startRun() ||
       (iov.startRun() == startRun() && iov.startSubrun() <= startSubrun()));
  bool ecp = (iov.endRun() > endRun() ||
              (iov.endRun() == endRun() && iov.endSubrun() >= endSubrun()));

  if (scp && ecp) return 1;     // complete overlap
  if ((!scp) && ecp) return 2;  // non-ovelapping piece at begining
  if (scp && (!ecp)) return 3;  // non-overlapping piece at end

  return 4;  // non-overlapping piece at begining and end
}

std::string mu2e::DbIoV::simpleString() const {
  std::ostringstream ss;
  ss << startRun() << " " << startSubrun() << " " << endRun() << " "
     << endSubrun();
  return ss.str();
}

std::string mu2e::DbIoV::to_string(bool compress) const {
  std::ostringstream ss;

  if (compress) {
    ss << startRun();
    bool oner = (startRun() == endRun());
    bool onesr = (startSubrun() == endSubrun());
    if (startSubrun() != 0 || (oner && onesr)) {
      ss << ":" << startSubrun();
    }
    bool wholerun = (startSubrun() == 0 && endSubrun() == maxSubrun());
    bool drop = oner && (onesr || wholerun);
    if (!drop) {
      ss << "-" << endRun();
      if (endSubrun() != maxSubrun()) {
        ss << ":" << endSubrun();
      }
    }
  } else {
    ss << std::setw(6) << startRun() << ":";
    ss << std::setw(6) << startSubrun() << "-";
    ss << std::setw(6) << endRun() << ":";
    ss << std::setw(6) << endSubrun();
  }

  return ss.str();
}

void mu2e::DbIoV::setByString(std::string iovstr) {
  boost::trim(iovstr);  // remove leading, trailing whitespace
  boost::to_upper(iovstr);
  if (iovstr.empty() || iovstr == "MAX" || iovstr == "ALL") {
    setMax();
    return;
  }
  if (iovstr == "EMPTY") {
    set(0, 0, 0, 0);
    return;
  }
  std::vector<std::string> words;

  boost::split(words, iovstr, boost::is_any_of(" \t"),
               boost::token_compress_on);
  if (words.size() != 1) {
    throw cet::exception("DBIOV_MULTIWORD_INIT_STRING")
        << "DbIoV::setByString string has multiple words: " << iovstr << "\n";
  }

  std::string start, end;
  boost::split(words, iovstr, boost::is_any_of("-"), boost::token_compress_off);
  if (words.size() == 1) {
    start = iovstr;
    end = iovstr;
  } else if (words.size() == 2) {
    start = words[0];
    end = words[1];
  } else {
    throw cet::exception("DBIOV_DASH_INIT_STRING")
        << "DbIoV::setByString wrong number of dashed fields: " << iovstr
        << "\n";
  }

  // require a string for start and end
  if (start.empty() || end.empty()) {
    throw cet::exception("DBIOV_EMPTY_INIT_STRING")
        << "DbIoV::setByString found start or end point was blank: " << iovstr
        << "\n";
  }

  boost::split(words, start, boost::is_any_of(":"), boost::token_compress_on);
  uint32_t startr = 0, startsr = 0;

  if (words.size() >= 1) {
    if (words[0] != "MIN") {
      startr = std::stoi(words[0]);
    }
  }
  if (words.size() >= 2) {
    if (words[1] != "MIN") {
      startsr = std::stoi(words[1]);
    }
  }

  boost::split(words, end, boost::is_any_of(":"), boost::token_compress_on);
  uint32_t endr = maxRun(), endsr = maxSubrun();
  if (words.size() >= 1) {
    if (words[0] != "MAX") {
      endr = std::stoi(words[0]);
    }
  }
  if (words.size() >= 2) {
    if (words[1] != "MAX") {
      endsr = std::stoi(words[1]);
    }
  }

  if (endr < startr || (endr == startr && endsr < startsr)) {
    throw cet::exception("DBIOV_NEGATIVE_INIT_STRING")
        << "DbIoV::setByString found end point was before start point: "
        << iovstr << "\n";
  }

  set(startr, startsr, endr, endsr);
}
