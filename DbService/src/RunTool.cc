#include "Offline/DbService/inc/RunTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
#include "Offline/GeneralUtilities/inc/splitString.hh"
#include "cetlib_except/exception.h"

using namespace mu2e;

//**************************************************

RunTool::RunTool() {
  DbIdList idList;  // read file of db connection details
  _reader.setDbId(idList.getDbId("run_info"));
  _reader.setVerbose(0);
}

//**************************************************
std::map<int, std::string> RunTool::flags(FlagType ftype) {
  std::map<int, std::string> flags;
  // run query url, answer returned in csv
  int rc;
  std::string tcsv, select, table, order;
  StringVec where;

  select = "";
  if (ftype == FlagType::run) {
    table = "production.run_type";
  } else if (ftype == FlagType::transition) {
    table = "production.transition_type";
  } else if (ftype == FlagType::cause) {
    table = "production.cause_type";
  }
  // where.push_back("name:eq:" + name);
  rc = _reader.query(tcsv, select, table, where, order);
  if (rc != 0 || tcsv.empty()) {
    throw cet::exception("RUNTOOL_BAD_FLAG_QUERY")
        << " RunTool::get failed to retrieve flags for " << table << "\n";
  }

  StringVec rows = splitString(tcsv, "\n", "\"", "\\", true, false);
  for (auto const& csv : rows) {
    StringVec sv = splitString(csv);
    flags.insert({std::stoi(sv[0]), sv[1]});
  }
  return flags;
}

//**************************************************
RunInfo::RunVec RunTool::listRuns(const RunSelect& runsel, bool conditions,
                                  bool transitions) {
  RunInfo::RunVec runv;
  std::string tcsv, select, table, order;
  StringVec where;
  std::string conditionss, transitionss;

  // fetch the list of all runs

  select = "";
  table = "production.run_configuration";
  order = "-run_number";
  int rc = _reader.query(tcsv, select, table, where, order);
  if (rc != 0 || tcsv.empty()) {
    throw cet::exception("RUNTOOL_BAD_RUN_QUERY")
        << " RunTool::listRuns failed to retrieve runs \n";
  }

  StringVec rows = splitString(tcsv, "\n", "\"", "\\", true, false);
  // if (rows.size() > 0) {
  //  std::cout << "HEADERS"<<rows[0] <<std::endl;
  //  // header = rows[0];  The first row is column names
  //  rows.erase(rows.begin());
  //}

  // apply the selection and fetch conditions and transitions if requested

  // setup run cuts
  std::string rr = runsel.run();
  int rmin = 0;
  int rmax = 999999;
  StringVec rlims = splitString(rr, "-");
  if (rlims.size() > 0 && rlims[0].size() > 0) {
    rmin = std::stoi(rlims[0]);
  }
  if (rlims.size() > 1 && rlims[1].size() > 0) {
    rmax = std::stoi(rlims[1]);
  }
  if (rlims.size() == 1) {
    rmax = rmin;
  }

  // setup time cuts
  std::string tt = runsel.time();
  std::time_t tmin = 0;
  std::time_t tmax = 10000000000;
  StringVec tlims = splitString(tt, "/");
  if (tlims.size() > 0 && tlims[0].size() > 0) {
    rc = TimeUtility::parseTimeTZ(tlims[0], tmin);
    if (rc != 0) {
      throw cet::exception("RUNTOOL_BAD_TIME")
          << " RunTool::listRuns can't parse time" << tlims[0] << " \n";
    }
  }
  if (tlims.size() > 1 && tlims[1].size() > 0) {
    rc = TimeUtility::parseTimeTZ(tlims[1], tmax);
    if (rc != 0) {
      throw cet::exception("RUNTOOL_BAD_TIME")
          << " RunTool::listRuns can't parse time" << tlims[1] << " \n";
    }
  }

  // the list of acceptable run types
  std::string ty = runsel.type();
  StringVec tys = splitString(ty, ",");
  std::vector<int> itypes;
  for (auto t : tys) itypes.emplace_back(std::stoi(t));

  std::time_t dtime = 0;
  if (runsel.days() > 0) {
    dtime = std::time(0) - runsel.days() * 24 * 60 * 60;
  }
  // loop over all runs, starting with highest run number
  for (auto const& csv : rows) {
    // apply cuts on this run
    bool pass = true;
    // use this temp RunInfo to parse the csv
    RunInfo ri(csv, std::string(), std::string());
    int run = ri.runNumber();
    // limit on the range
    if (run < rmin || run > rmax) pass = false;
    // limit on the number of runs
    if (runsel.last() > 0 && runv.size() >= size_t(runsel.last())) pass = false;
    std::string srtime = ri.commitTime();
    std::time_t rtime;
    rc = TimeUtility::parseTimeTZ(srtime, rtime);
    if (rc != 0) {
      throw cet::exception("RUNTOOL_BAD_TIME")
          << " RunTool::listRuns can't parse time" << srtime << " \n";
    }
    // limit on time range
    if (rtime < tmin || rtime > tmax) pass = false;
    // limt on daysage
    if (rtime < dtime) pass = false;
    // limits on run types
    if (itypes.size() > 0) {
      int type = ri.runType();
      if (std::find(itypes.begin(), itypes.end(), type) == itypes.end())
        pass = false;
    }

    if (pass) {
      if (conditions) {
        tcsv.clear();
        select = "condition";
        table = "production.run_condition";
        where.clear();
        std::string ww =
            std::string("condition_id:eq:") + std::to_string(ri.conditionId());
        where.emplace_back(ww);
        order.clear();
        int rc = _reader.query(tcsv, select, table, where, order);
        if (rc != 0) {
          throw cet::exception("RUNTOOL_BAD_CONDITIONS_QUERY")
              << " RunTool::listRuns failed to retrieve conditions \n";
        }
        conditionss = tcsv;
      }  // end conditions

      if (transitions) {
        StringVec sv = splitString(csv);
        tcsv.clear();
        select = "transition_type,cause_type,transition_time";
        table = "production.run_transition";
        where.clear();
        std::string ww =
            std::string("run_number:eq:") + std::to_string(ri.runNumber());
        where.emplace_back(ww);
        order = "transition_time";
        int rc = _reader.query(tcsv, select, table, where, order);
        if (rc != 0) {
          throw cet::exception("RUNTOOL_BAD_TRANSITION_QUERY")
              << " RunTool::listRuns failed to retrieve transitions for run "
              << sv[0] << "\n";
        }
        transitionss = tcsv;
      }  // end transitions

      runv.emplace_back(csv, conditionss, transitionss);
    }  // end if pass
  }

  return runv;
}

//**************************************************
void RunTool::printRun(const RunInfo& info, bool longp) {
  std::string committ = info.commitTime().substr(0, size_t(19));
  if (longp) {
    std::cout << "\n";
    std::cout << "runNumber             : " << info.runNumber() << "\n";
    std::cout << "runType               : " << info.runType() << "\n";
    std::cout << "commitTime            : " << committ << "\n";
    std::cout << "conditionId           : " << info.conditionId() << "\n";
    std::cout << "artdaqPartition       : " << info.artdaqPartition() << "\n";
    std::cout << "hostName              : " << info.hostName() << "\n";
    std::cout << "configurationName     : " << info.configurationName() << "\n";
    std::cout << "configurationVersion  : " << info.configurationVersion()
              << "\n";
    std::cout << "contextName           : " << info.contextName() << "\n";
    std::cout << "contextVersion        : " << info.contextVersion() << "\n";
    std::cout << "triggerTableName      : " << info.triggerTableName() << "\n";
    std::cout << "triggerTableVersion   : " << info.triggerTableVersion()
              << "\n";
    std::cout << "onlineSoftwareVersion : " << info.onlineSoftwareVersion()
              << "\n";
    std::cout << "shifterNote           : " << info.shifterNote() << "\n";
  } else {
    std::cout << std::setw(8) << info.runNumber();
    std::cout << std::setw(5) << info.runType();
    std::cout << std::setw(23) << committ;
    std::cout << "\n";
  }

  if (info.conditions().size() > 0) {
    std::cout << "conditions:\n" << info.conditions() << "\n";
  }
  if (info.transitions().size() > 0) {
    std::cout << "transitions:\n" << info.transitions();
  }
  // std::cout << " : " << info. << "\n";

  //  const std::string& conditions() const { return _conditions; }
  // const std::string& transitions() const { return _transitions; }
}
/*
EpicsVar::EpicsVec RunTool::get(std::string const& name,
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
        << " RunTool::get failed to find channel number for name\n";
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
        << "RunTool::get could not lookup data for channel id " << id << "\n";
  }

  // the table is delivered with rows separated by \n, split them
  StringVec rows = splitString(tcsv, "\n", "\"", "\\", true, false);

  for (auto const& csv : rows) {
    // Example csv:
    //
channel_id,smpl_time,nanosecs,severity_id,status_id,num_val,float_val,str_val,datatype
    // 642,2022-05-05 14:10:02.089404-05:00,89404264,5,9,None,55.0,None, ,None

    // split on commas
    StringVec sv = splitString(csv);

    if (sv.size() != 10) {
      throw cet::exception("EPICSTOOL_BAD_CSV")
          << " RunTool::get number of csv fields !=10, csv: \n"
          << csv << " \n";
    }

    int64_t channel_id = std::stoll(sv[0]);
    std::time_t smpl_time_t;
    if (TimeUtility::parseTimeTZ(sv[1], smpl_time_t)) {
      throw cet::exception("EPICSTOOL_BAD_TIME")
          << " RunTool::get could not parse time " << sv[1] << " \n";
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
          << " RunTool::get could not parse csv row:\n"
          << csv << " \n";
    }

    ev.emplace_back(csv, channel_id, sv[1], smpl_time_t, nanosecs, severity_id,
                    status_id, value);

  }  // loop over db sample rows

  if (ev.empty()) {
    return ev;
  }

  // if no constraints, return the whole list
  if ( time.empty() && daysAgo<=0.0 ) {
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
          << " RunTool::get could not parse time " << st1 << " \n";
    }
    std::time_t t2;
    if (TimeUtility::parseTimeTZ(st2, t2)) {
      throw cet::exception("EPICSTOOL_BAD_TIME")
          << " RunTool::get could not parse time " << st2 << " \n";
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
          << " RunTool::get could not parse time " << time << " \n";
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
        << " RunTool::get could not parse arguments "
        << " \n";
  }

  return evt;
}

// **************************************************

int RunTool::names(StringVec& names) {
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
    std::cout << "Error - RunTool could not list names" << std::endl;
    return 1;
  }

  names = splitString(csv, "\n", "\"", "\\", true, false);

  return 0;
}
*/
