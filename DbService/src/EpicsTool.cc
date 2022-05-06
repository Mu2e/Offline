#include "Offline/DbService/inc/EpicsTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

using namespace mu2e;
using namespace boost;

//**************************************************

EpicsTool::EpicsTool() {
  DbIdList idList;  // read file of db connection details
  _reader.setDbId(idList.getDbId("dcs_archive"));
  _reader.setVerbose(0);
}

//**************************************************

EpicsVar::EpicsVec EpicsTool::get(std::string const& name,
                                  std::string const& time) {
  EpicsVar::EpicsVec ev;

  // run query url, answer returned in csv
  int rc;
  std::string csv, select, table, order;
  StringVec where;

  select = "channel_id";
  table = "channel";
  where.push_back("name:eq:" + name);
  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0 || csv.empty()) {
    std::cout << "Error - EpicsTool could not look for channel id for name "
              << name << std::endl;
  }

  std::string id = csv;
  size_t inl = id.find('\n');
  if (inl != std::string::npos) {
    id.erase(inl, 1);
  }

  select.clear();
  table = "sample";
  where.clear();
  where.push_back("channel_id:eq:" + id);
  order = "smpl_time";
  csv.clear();
  rc = _reader.query(csv, select, table, where, order);

  if (rc != 0) {
    std::cout << "Error - EpicsTool could not lookup data for channel id " << id
              << std::endl;
  }

  escaped_list_separator<char> els('\\', '\n', '"');
  tokenizer<escaped_list_separator<char>> tok(csv, els);
  for (tokenizer<escaped_list_separator<char>>::iterator beg = tok.begin();
       beg != tok.end(); ++beg) {
    std::string ss(*beg);
    boost::trim(ss);
    if (!ss.empty()) ev.emplace_back(ss);
    EpicsVar& v = ev.back();
    v.setTime(parseTimeTZ(v.stime()));
  }

  if (ev.size() == 0) {
    return ev;
  }
  if (time.empty()) {
    return ev;
  }

  EpicsVar::EpicsVec evt;
  if (time.find('/') != std::string::npos) {
    // time interval
    size_t idiv = time.find('/');
    std::string st1 = time.substr(0, idiv);
    std::string st2 = time.substr(idiv + 1, std::string::npos);
    std::time_t t1 = parseTimeTZ(st1);
    std::time_t t2 = parseTimeTZ(st2);
    for (auto& e : ev) {
      if (e.ttime() >= t1 && e.ttime() <= t2) {
        evt.push_back(e);
      }
    }
  } else if (time.find('-') != std::string::npos) {
    // single time, find nearest entry before
    std::time_t tt = parseTimeTZ(time);
    std::cout << "parse t1 " << tt << std::endl;
    size_t i = ev.size() - 1;
    while (i > 0 && ev[i].ttime() > tt) {
      i--;
    }
    evt.push_back(ev[i]);
  } else {
    // integer or float back n days
    std::time_t now = std::time(nullptr);
    std::time_t delta = std::time_t(std::stof(time) * 86400.0);
    std::time_t tmin = now - delta;
    for (auto& e : ev) {
      if (e.ttime() > tmin) {
        evt.push_back(e);
      }
    }
  }

  return evt;
}

//**************************************************

int EpicsTool::names(StringVec& names) {
  // run query url, answer returned in csv
  int rc;
  std::string csv, select, table, order;
  StringVec where;

  names.clear();

  select = "name";
  table = "channel";
  order = "name";
  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0 || csv.empty()) {
    std::cout << "Error - EpicsTool could not list names" << std::endl;
    return 1;
  }

  escaped_list_separator<char> els('\\', '\n', '"');
  tokenizer<escaped_list_separator<char>> tok(csv, els);
  for (tokenizer<escaped_list_separator<char>>::iterator beg = tok.begin();
       beg != tok.end(); ++beg) {
    std::string ss(*beg);
    boost::trim(ss);
    if (!ss.empty()) names.emplace_back(ss);
  }

  return 0;
}

//**************************************************

std::time_t EpicsTool::parseTimeTZ(std::string const& stime) const {
  // allowed formats:
  //  2018-10-12
  //  2018-10-12T08:58:26-05:00
  //  2018-10-12 08:58:26-05:00
  //  2018-10-12T08:58:26.792518-05:00
  //  2018-10-12 08:58:26.792518-05:00

  // must at least be a date
  if (stime.size() < 10) {
    return 0;
  }

  std::string ftime(24, ' ');
  std::string zero = "2000-01-01 00:00:00";
  // copy data and time part
  for (size_t ipt = 0; ipt < 19; ipt++) {
    if (ipt != 10) {
      if (ipt < stime.size()) {
        ftime[ipt] = stime[ipt];
      } else {
        ftime[ipt] = zero[ipt];
      }
    }
  }

  // skip ns if present
  size_t ipt = 19;
  if (stime[19] == '.') {
    // advance past ns
    while (stime[ipt] != '-' && stime[ipt] != '+' && ipt < stime.size()) ipt++;
  }

  // sign and hour digits of time zone offset
  int tzi = 0;
  if (stime.size() >= ipt + 3) {
    std::string tz = stime.substr(ipt, 3);
    tzi = std::stoi(tz);
  }

  std::tm tt{};  // must initialize to zero
  std::istringstream ss(ftime);
  ss >> std::get_time(&tt, "%Y-%m-%d %H:%M:%S");
  // these don't usually get set, but make sure they aren't used
  tt.tm_isdst = 0;
  tt.tm_gmtoff = 0;
  // here is the main point - interpret date-time as UTC,
  // convert to epoch time, then correct for explicit tz offset
  std::time_t t = timegm(&tt) - tzi * 3600;

  return t;
}
