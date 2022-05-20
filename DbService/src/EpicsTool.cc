#include "Offline/DbService/inc/EpicsTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
#include "Offline/GeneralUtilities/inc/splitString.hh"
#include "cetlib_except/exception.h"

using namespace mu2e;

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
  std::string tcsv, select, table, order;
  StringVec where;

  select = "channel_id";
  table = "channel";
  where.push_back("name:eq:" + name);
  rc = _reader.query(tcsv, select, table, where, order);
  if (rc != 0 || tcsv.empty()) {
    throw cet::exception("EPICSTOOL_BAD_CHANNEL")
        << " EpicsTool::get failed to find channel number for name\n";
  }

  std::string id = tcsv;
  size_t inl = id.find('\n');
  if (inl != std::string::npos) {
    id.erase(inl, 1);
  }

  select.clear();
  table = "sample";
  where.clear();
  where.push_back("channel_id:eq:" + id);
  order = "smpl_time";
  tcsv.clear();
  rc = _reader.query(tcsv, select, table, where, order);

  if (rc != 0) {
    throw cet::exception("EPICSTOOL_BAD_SAMPLE")
        << "EpicsTool::get could not lookup data for channel id " << id << "\n";
  }

  // the table is delivered with rows separated by \n, split them
  StringVec rows = splitString(tcsv, "\n", "\"", "\\", true, false);

  for (auto const& csv : rows) {
    // Example csv:
    // channel_id,smpl_time,nanosecs,severity_id,status_id,num_val,float_val,str_val,datatype
    // 642,2022-05-05 14:10:02.089404-05:00,89404264,5,9,None,55.0,None, ,None

    // split on commas
    StringVec sv = splitString(csv);

    if (sv.size() != 10) {
      throw cet::exception("EPICSTOOL_BAD_CSV")
          << " EpicsTool::get number of csv fields !=10, csv: \n"
          << csv << " \n";
    }

    int64_t channel_id = std::stoll(sv[0]);
    std::time_t smpl_time_t;
    if (TimeUtility::parseTimeTZ(sv[1], smpl_time_t)) {
      throw cet::exception("EPICSTOOL_BAD_TIME")
          << " EpicsTool::get could not parse time " << sv[1] << " \n";
    }
    int64_t nanosecs = std::stoll(sv[2]);
    int64_t severity_id = std::stoll(sv[3]);
    int64_t status_id = std::stoll(sv[4]);
    EpicsVar::variVar value;
    if (sv[6] == "None" && sv[7] == "None") {
      value = std::stoi(sv[5]);
    } else if (sv[5] == "None" && sv[7] == "None") {
      value = std::stod(sv[6]);
    } else if (sv[5] == "None" && sv[6] == "None") {
      value = sv[7];
    } else {
      throw cet::exception("EPICSTOOL_BAD_CSV_ROW")
          << " EpicsTool::get could not parse csv row:\n"
          << csv << " \n";
    }

    ev.emplace_back(csv, channel_id, sv[1], smpl_time_t, nanosecs, severity_id,
                    status_id, value);

  }  // loop over db sample rows

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
        << " EpicsTool::get could not parse arguments "
        << " \n";
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

  names = splitString(csv, "\n", "\"", "\\", true, false);

  return 0;
}
