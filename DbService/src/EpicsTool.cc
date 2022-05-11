#include "Offline/DbService/inc/EpicsTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
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
                                  std::string const& time, float daysAgo) {
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
    throw cet::exception("EPICSVAR_BAD_CSV")
        << " EpicVar::ctor csv fields !=10 : " << csv << " \n";
    return ev;
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
    return ev;
  }

  escaped_list_separator<char> els('\\', '\n', '"');
  tokenizer<escaped_list_separator<char>> tok(csv, els);
  for (tokenizer<escaped_list_separator<char>>::iterator beg = tok.begin();
       beg != tok.end(); ++beg) {
    std::string ss(*beg);
    boost::trim(ss);
    if (!ss.empty()) ev.emplace_back(ss);
  }

  if (ev.empty()) {
    return ev;
  }

  EpicsVar::EpicsVec evt;
  if (time.find('/') != std::string::npos) {
    // time interval
    size_t idiv = time.find('/');
    std::string st1 = time.substr(0, idiv);
    std::string st2 = time.substr(idiv + 1, std::string::npos);
    std::time_t t1;
    if (TimeUtility::parseTimeTZ(st1, t1)) {
      throw cet::exception("EPICSTOOL_BAD_TIME")
          << " EpicsTool::get could not parse time " << st1 << " \n";
    }
    std::time_t t2;
    if (TimeUtility::parseTimeTZ(st2, t2)) {
      throw cet::exception("EPICSTOOL_BAD_TIME")
          << " EpicsTool::get could not parse time " << st2 << " \n";
    }
    for (auto& e : ev) {
      if (e.ttime() >= t1 && e.ttime() <= t2) {
        evt.push_back(e);
      }
    }
  } else if (time.find('-') != std::string::npos) {
    // single time, find nearest entry before
    std::time_t tt;
    if (TimeUtility::parseTimeTZ(time, tt)) {
      throw cet::exception("EPICSTOOL_BAD_TIME")
          << " EpicsTool::get could not parse time " << time << " \n";
    }
    size_t i = ev.size() - 1;
    while (i > 0 && ev[i].ttime() > tt) {
      i--;
    }
    evt.push_back(ev[i]);
  } else if (daysAgo > 0.0) {
    // integer or float back n days
    std::time_t now = std::time(nullptr);
    std::time_t delta = std::time_t(daysAgo * 86400.0);
    std::time_t tmin = now - delta;
    for (auto& e : ev) {
      if (e.ttime() > tmin) {
        evt.push_back(e);
      }
    }
  } else {
    throw cet::exception("EPICSTOOL_BAD_ARGS")
      << " EpicsTool::get could not parse arguments " << " \n";
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
