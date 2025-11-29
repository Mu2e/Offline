#include "Offline/DbService/inc/GrlTool.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/DbService/inc/RunTool.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
#include "Offline/GeneralUtilities/inc/splitString.hh"
#include <boost/algorithm/string.hpp>
#include <forward_list>
#include <fstream>

using namespace mu2e;

//**************************************************

GrlTool::GrlTool() : _verbose(0) {
  DbIdList idList;  // read the connections info file
  _id = idList.getDbId("mu2e_dqm_prd");

  _reader.setDbId(_id);
  _reader.setUseCache(false);
  _reader.setVerbose(_verbose);
  _reader.setTimeVerbose(_verbose);

  _sql.setDbId(_id);
  _sql.setVerbose(_verbose);
}

//**************************************************
void GrlTool::setVerbose(int verbose) {
  _verbose = verbose;
  _reader.setVerbose(verbose);
  _sql.setVerbose(verbose);
}

//**************************************************
int GrlTool::defineWord(const std::string& name, const std::string& desc) {
  std::string command, result;
  int rc = 0;

  // validate the name
  StringVec words = splitString(name, "/");
  if (words.size() != 3) {
    std::string mess("Error - word name must be category/subcategory/name : " +
                     name);
    throw std::runtime_error(mess);
  }
  for (auto const& ww : words) {
    if (ww.size() == 0) {
      std::string mess(
          "Error - fields in word name category/subcategory/name must be "
          "non-zero: " +
          name);
      throw std::runtime_error(mess);
    }
  }

  rc = _sql.connect();
  if (rc) {
    std::cout << "ERROR - createList failed to connect " << std::endl;
    return 3;
  }

  command = "SET ROLE grlwrite;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  // insert table values, if duplicate, SQL will fail
  command =
      "INSERT INTO grl.words (name,description,create_time,create_user)  "
      "VALUES ('" +
      name + "','" + desc + "',CURRENT_TIMESTAMP,SESSION_USER);";

  rc = _sql.execute(command, result);
  if (_verbose > 1) {
    std::cout << command << std::endl;
    std::cout << result << std::endl;
  }
  if (_verbose > 0) {
    if (rc == 0) {
      std::cout << "failed to create " << name << std::endl;
    } else {
      std::cout << "created new word name " << name << std::endl;
    }
  }
  return rc;
}

//**************************************************
int GrlTool::defineBit(const std::string& name, const std::string& bitname,
                       const std::string& desc) {
  std::string command, result;
  int rc = 0;

  // validate the bitname
  if (bitname.size() == 0) {
    std::cout << "Error - bitname must be non-null" << std::endl;
    return 1;
  }

  // figure out the bit number by checking what already exists
  rc = bits();
  _result.clear();
  if (rc != 0) return rc;
  int bitnumber = 0;
  for (const auto& b : _bits) {
    if (b.name() == name) {
      if (b.bitnumber() >= bitnumber) bitnumber = b.bitnumber() + 1;
    }
  }

  if (bitnumber >= 32) {
    std::cout << "Error - word " << name << "has no room for more bits"
              << std::endl;
    return 1;
  }

  rc = _sql.connect();
  if (rc) {
    std::cout << "ERROR - createList failed to connect " << std::endl;
    return 3;
  }

  command = "SET ROLE grlwrite;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  // insert table values,SQL will fail
  // - if name is not in name tables
  // - if name,bitname combo is not unique
  command =
      "INSERT INTO grl.bits "
      "(name,bitname,bitnumber,description,create_time,create_user) "
      " VALUES ('" +
      name + "','" + bitname + "'," + std::to_string(bitnumber) + ",'" + desc +
      "',CURRENT_TIMESTAMP,SESSION_USER);";

  rc = _sql.execute(command, result);
  if (_verbose > 1) {
    std::cout << command << std::endl;
    std::cout << result << std::endl;
  }
  if (_verbose > 0) {
    if (rc == 0) {
      std::cout << "failed to create " << name << " : " << bitname << std::endl;
    } else {
      std::cout << "created new bit name " << name << " : " << bitname
                << std::endl;
    }
  }
  return rc;
}

//**************************************************
int GrlTool::setBits(const std::string& name, const DbIoV& iov,
                     const StringVec& bitnames) {
  std::string command, result;
  int rc = 0;

  rc = words();
  _result.clear();
  if (rc != 0) return rc;
  rc = bits(name);
  _result.clear();
  if (rc != 0) return rc;
  int value = 0;
  for (const auto& bn : bitnames) {
    int subvalue = 0;
    for (auto const& b : _bits) {
      if (b.name() == name && b.bitname() == bn) subvalue = 1 << b.bitnumber();
      if (_verbose > 2)
        std::cout << "confirmed match for " << name << " " << bn << std::endl;
    }
    if (subvalue == 0) {
      std::cout << "Error - word name " << name << " or bit name " << bn
                << " did not match declared names" << std::endl;
      return 1;
    }
    value |= subvalue;
    if (_verbose > 1) std::cout << "new word value " << value << std::endl;
  }

  // check if this new setting overlaps previous settings and if so
  // retire the old setting, reinstall any parts not overlapping,
  // then insert the new setting
  rc = entries();
  _result.clear();
  if (rc != 0) return rc;

  // collect retire and add commands
  StringVec commands;
  command = "SET ROLE grlwrite;";
  commands.emplace_back(command);

  for (const auto& ee : _entries) {
    if (ee.name() == name) {
      // 0 = no overlap
      // 1 = iov bigger than ee.iov,
      // 2 = piece of ee.iov left over at beginning
      // 2 = piece ee.iov left over at end
      // 4  = iov contained in ee.iov, piece at beginning and end
      int iover = ee.iov().isOverlapping(iov);
      if (iover != 0) {  // some degree of overlap
        // old entry must be retired and replaced
        command =
            "UPDATE grl.entries "
            "SET retired = 1 WHERE eid = " +
            std::to_string(ee.eid()) + ";";
        commands.emplace_back(command);

        if (iover == 2 || iover == 4) {
          // there is a early piece of the old iov that should continue
          DbIoV piov = ee.iov();
          // subtract iov from ee.iov, keeping the early part
          piov.subtract(iov, piov.startRun(), piov.startSubrun());
          command =
              "INSERT INTO grl.entries "
              "(name,iov,value,retired,create_time,create_user)  "
              "VALUES ('" +
              ee.name() + "','" + piov.to_string(true) + "'," +
              std::to_string(ee.value()) +
              ",0,CURRENT_TIMESTAMP,SESSION_USER);";
          commands.emplace_back(command);
        }
        if (iover == 3 || iover == 4) {
          // there is a late piece of the old iov that should continue
          DbIoV piov = ee.iov();
          // subtract iov from ee.iov, keeping the late part
          piov.subtract(iov, piov.endRun(), piov.endSubrun());
          command =
              "INSERT INTO grl.entries "
              "(name,iov,value,retired,create_time,create_user)  "
              "VALUES ('" +
              ee.name() + "','" + piov.to_string(true) + "'," +
              std::to_string(ee.value()) +
              ",0,CURRENT_TIMESTAMP,SESSION_USER);";
          commands.emplace_back(command);
        }
      }
    }
  }  // end loop over entries

  // always insert the new value
  command =
      "INSERT INTO grl.entries "
      "(name,iov,value,retired,create_time,create_user)  "
      "VALUES ('" +
      name + "','" + iov.to_string(true) + "'," + std::to_string(value) +
      ",0,CURRENT_TIMESTAMP,SESSION_USER);";
  commands.push_back(command);

  StringVec results;
  rc = _sql.transact(commands, results);
  if (rc) {
    std::cout << "ERROR - setBits failed to complete transaction " << std::endl;
    return 1;
  }

  return rc;
}

//**************************************************
int GrlTool::words() {
  std::string command, result;
  int rc = 0;

  std::string csv;
  std::string select;
  std::string table("grl.words");
  StringVec where;
  std::string order("name");

  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - words query failed " << std::endl;
    return 1;
  }
  StringVec lines = DbUtil::splitCsvLines(csv);
  for (auto& line : lines) {
    StringVec words = DbUtil::splitCsv(line);
    _words.emplace_back(words[0], words[1], words[2]);
  }
  _result = csv;
  return 0;
}

//**************************************************
int GrlTool::bits(const std::string& name) {
  _result.clear();
  std::string command, result;
  int rc = 0;

  std::string csv;
  std::string select;
  std::string table("grl.bits");
  StringVec where;
  std::string order("name,bitname");

  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - getBits query failed " << std::endl;
    return 1;
  }
  StringVec lines = DbUtil::splitCsvLines(csv);
  for (auto& line : lines) {
    StringVec words = DbUtil::splitCsv(line);

    //    GrlBit(std::string const& name, std::string const& bitname, int
    //    bitnumber, std::string const& description = "",
    //     std::string const& create_time = "",
    //     std::string const& create_user = "") :

    _bits.emplace_back(words[0], words[1], std::stoi(words[2]), words[3],
                       words[4], words[5]);
    if (name.empty() || name == words[0]) {
      _result.append(line + "\n");
    }
  }

  return 0;
}

//**************************************************
int GrlTool::entries(const std::string& name, const DbIoV& iov) {
  std::string command, result;
  _result.clear();
  int rc = 0;

  std::string csv;
  std::string select("name,iov,value,eid,retired,create_time,create_user");
  std::string table("grl.entries");
  StringVec where;
  where.emplace_back("retired:eq:0");
  std::string order("create_time");

  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - getEntries query failed " << std::endl;
    return 1;
  }

  StringVec lines = DbUtil::splitCsvLines(csv);

  for (auto& line : lines) {
    StringVec words = DbUtil::splitCsv(line);
    // name, iov, value,eid, retired,create_date,create_user
    _entries.emplace_back(words[0], words[1], std::stoi(words[2]),
                          std::stoi(words[3]), std::stoi(words[4]), words[5],
                          words[6]);
    if(name.empty() || words[0]==name) {
      if(iov.isNull() || iov.isOverlapping(DbIoV(words[1]))>0 ) {
        result.append(line+'\n');
      }
    }
  }
  _result = result;
  return 0;
}

//**************************************************
int GrlTool::genList(const std::string& comFile, const std::string& base) {
  int rc = 0;

  if (comFile.empty()) {
    std::cout << "ERROR - genList no command files " << std::endl;
    return 1;
  }

  std::ifstream myfile;
  myfile.open(comFile);
  if (!myfile.is_open()) {
    std::cout << "ERROR - genList could not open file " << comFile << std::endl;
    return 1;
  }

  std::string line;
  StringVec coms;
  while (std::getline(myfile, line)) {
    boost::trim(line);  // remove whitespace
    if (line.size() <= 0 || line[0] == '#') continue;
    coms.emplace_back(line);
  }
  if (_verbose > 0) {
    std::cout << "read " << coms.size() << " commands from " << comFile
              << std::endl;
  }

  GrlList::IoVVec iv;  // the list of iovs (runs) based on run type
  if (base.empty()) {
    DbIoV iov;
    iov.setMax();
    iv.emplace_back(iov);
  } else {
    RunSelect rs("0-999999", 0, base, "", 0);
    RunTool rtool;
    RunInfo::RunVec rv = rtool.listRuns(rs);
    if (_verbose > 0) {
      std::cout << "selection of " << base << " returned " << rv.size()
                << " runs" << std::endl;
    }
    for (const auto& rr : rv) {
      uint32_t run = rr.runNumber();
      DbIoV riov(run, 0, run, DbIoV::maxSubrun());
      iv.emplace_back(riov);
    }
  }
  words();
  _result.clear();
  if (rc != 0) return rc;
  bits();
  _result.clear();
  if (rc != 0) return rc;
  entries();
  _result.clear();
  if (rc != 0) return rc;

  GrlList::IoVVec ivout = iv;  // output iov list
  // apply each command
  for (auto const& com : coms) {
    // for each command start with previous result as baseline
    iv = ivout;
    ivout.clear();

    if (_verbose > 1) {
      std::cout << "processing " << com << std::endl;
    }
    StringVec words = splitString(com, " ", "\"", "\\", true, false);
    if (words.size() < 2) {
      std::cout << "Error parsing line: " << com << std::endl;
      return 1;
    }
    std::string action = words[0];
    std::string name = words[1];
    bool found = false;
    for (auto const& ww : _words) {
      if (name == ww.name()) found = true;
    }
    if (!found) {
      std::cout << "Error word " << name << " not found " << std::endl;
      return 1;
    }
    int mask = 0;
    // rest of the words are bits
    for (size_t i = 2; i < words.size(); i++) {
      found = false;
      for (auto const& bb : _bits) {
        if (bb.name() == name && bb.bitname() == words[i]) {
          mask |= (1 << bb.bitnumber());
        }
      }
    }  // loop over bits
    // mask can still be 0 if there were no bits
    if (_verbose > 1) {
      std::cout << "com " << words[0] << "  word " << name << "  mask 0x"
                << std::hex << mask << std::dec << std::endl;
    }
    if (action == "REQUIRE") {
      // this assume both the orginal runs list from the run type
      // and the bit list do not have overlaps among each list
      // so we just take anything that overlaps
      for (auto const& ee : _entries) {
        if (ee.name() == name && ee.retired() == 0) {
          for (const auto& ii : iv) {
            if ((ee.value() & mask) == mask) {
              DbIoV niov(ii);
              niov.overlap(ee.iov());
              if (_verbose > 2) {
                std::cout << "require compare " << ii.to_string(true) << " and "
                          << ee.iov().to_string(true) << " overlap "
                          << niov.to_string(true) << std::endl;
              }
              if (!niov.isNull()) ivout.emplace_back(niov);
            }
          }
        }
      }  // loop over entries

    } else if (action == "REJECT") {
      // this is trickier than require, since removals
      // and can remove and replace the iovs
      // start with the current list of good runs
      std::list<DbIoV> ivl(iv.begin(), iv.end());
      auto it = ivl.begin();
      while (it != ivl.end()) {
        if (_verbose > 5)
          std::cout << "ivl " << it->to_string(true) << std::endl;
        for (auto const& ee : _entries) {
          if (_verbose > 5)
            std::cout << "entry " << ee.iov().to_string(true) << std::endl;
          if (ee.name() == name && ee.retired() == 0 &&
              (ee.value() & mask) != 0) {
            int iover = it->isOverlapping(ee.iov());
            if (iover == 3 || iover == 4) {
              // there is a late piece of the grl iov that should continue
              DbIoV iiov = *it;
              // subtract bad runs, ee.iov, from grl, keeping the early part
              iiov.subtract(ee.iov(), iiov.endRun(), iiov.endSubrun());
              // insert the surviving piece
              if (_verbose > 5)
                std::cout << "insert1 " << iiov.to_string(true) << std::endl;
              ivl.insert(std::next(it), iiov);
            }
            if (iover == 2 || iover == 4) {
              // there is a early piece of the grl iov that should continue
              DbIoV iiov = *it;
              // subtract bad runs, ee.iov, from grl, keeping the early part
              iiov.subtract(ee.iov(), iiov.startRun(), iiov.startSubrun());
              // insert the surviving piece
              if (_verbose > 5)
                std::cout << "insert2 " << iiov.to_string(true) << std::endl;
              ivl.insert(std::next(it), iiov);
            }
            if (iover != 0) {
              // then the old grl iov was rejected or split, remove original
              if (_verbose > 5)
                std::cout << "erase " << it->to_string(true) << std::endl;
              it = ivl.erase(it);
            }
          }  // iff ee is relevant
        }    // loop over entries
        it++;
      }
      // after these mecanations, need to save the result for next action
      ivout.clear();
      for (const auto& x : ivl) {
        ivout.emplace_back(x);
      }
    } else {
      std::cout << "Error - unknown action " << words[0] << std::endl;
      return 1;
    }
  }  // loop over coms

  // save output
  for (const auto& ii : ivout) {
    _result.append(ii.to_string(true) + "\n");
  }

  /*
  std::string csv;
  std::string select("name,iov,value,eid,retired,create_time,create_user");
  std::string table("grl.entries");
  StringVec where;
  where.emplace_back("retired:eq:0");
  std::string order("create_time");

  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - getEntries query failed " << std::endl;
    return 1;
  }

  StringVec lines = DbUtil::splitCsvLines(csv);
  for (auto& line : lines) {
    StringVec words = DbUtil::splitCsv(line);
    // name, iov, value,eid, retired,create_date,create_user
    _entries.emplace_back(words[0], words[1], std::stoi(words[2]),
                          std::stoi(words[3]), std::stoi(words[4]), words[5],
                          words[6]);
  }
  _result = csv;
  */
  return 0;
}

//**************************************************
int GrlTool::lists() {
  if (_verbose > 1) std::cout << "starting lists" << std::endl;

  std::string command, result;
  int rc = 0;

  std::string csv;
  std::string select("name,locked,lid,create_time,create_user");
  std::string table("grl.lists");
  StringVec where;
  std::string order("name");

  //_reader.setVerbose(10);
  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - loadLists query failed " << std::endl;
    return 1;
  }
  if (_verbose > 1) std::cout << csv << std::endl;

  StringVec lines = DbUtil::splitCsvLines(csv);
  for (auto& line : lines) {
    StringVec words = DbUtil::splitCsv(line);
    _headers.emplace_back(words[0], std::stoi(words[1]), std::stoi(words[2]),
                          words[3], words[4]);
  }
  for (auto const& hh : _headers) {
    _result.append(hh.formatted() + "\n");
  }
  if (_verbose > 1) std::cout << "lists returning 0" << std::endl;
  return 0;
}

//**************************************************
GrlList GrlTool::list(const std::string& name) {
  if (_verbose > 1) std::cout << "starting list" << std::endl;
  int rc = 0;
  GrlHeader header("empty");

  // get definitions of available lists
  rc = lists();
  if (rc != 0) return GrlList(header, GrlList::IoVVec());

  _result.clear();

  for (auto& hh : _headers) {
    if (hh.name() == name) header = hh;
  }
  if (header.lid() < 0) {
    std::cout << "ERROR - GRL name not found in database: " << name
              << std::endl;
    return GrlList(header, GrlList::IoVVec());
  }
  if (_verbose > 1)
    std::cout << "list function found the list\n"
              << header.formatted() << std::endl;
  // read the database table of IoV's for this GRL

  std::string command, result;

  std::string csv;
  std::string select("iov,create_time,create_user");
  std::string table("grl.iovs");
  StringVec where;
  where.emplace_back("lid:eq:" + std::to_string(header.lid()));
  std::string order("iov");

  //_reader.setVerbose(10);
  rc = _reader.query(csv, select, table, where, order);
  if (rc != 0) {
    std::cout << "ERROR - getList query failed " << std::endl;
    return GrlList(header, GrlList::IoVVec());
  }
  if (_verbose > 1) {
    std::cout << "query returned " << std::endl;
    std::cout << csv;
  }
  StringVec lines = DbUtil::splitCsvLines(csv);
  GrlList::IoVVec iovs;
  for (auto& line : lines) {
    StringVec words = DbUtil::splitCsv(line);
    iovs.emplace_back(words[0]);
    _result.append(words[0] + std::string("\n"));
  }

  return GrlList(header, iovs);
}

//**************************************************
int GrlTool::createList(const GrlList& list) {
  if (_verbose > 1) std::cout << "running createList" << std::endl;

  int rc = 0;
  const GrlHeader& header = list.header();

  rc = _sql.connect();
  if (rc) {
    std::cout << "ERROR - createList failed to connect " << std::endl;
    return 3;
  }

  std::string command, result;
  command = "BEGIN";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  command = "SET ROLE grlwrite;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  // insert table values
  command =
      "INSERT INTO grl.lists ( "
      "name,locked,create_time,create_user) VALUES ('" +
      header.name() + "'," + std::to_string(header.locked()) +
      ",CURRENT_TIMESTAMP,SESSION_USER) RETURNING lid;";
  rc = _sql.execute(command, result);
  if (_verbose > 1) {
    std::cout << command << std::endl;
    std::cout << result << std::endl;
  }
  if (rc != 0) return rc;
  int lid = std::stoi(result);
  if (_verbose > 1) std::cout << "return lid = " << lid << std::endl;

  if (_verbose > 1) std::cout << "committing IoVs" << std::endl;

  for (const auto& iov : list.grl()) {
    command =
        "INSERT INTO grl.iovs ("
        "lid,iov,create_time,create_user) VALUES ( " +
        std::to_string(lid) + ",'" + iov.to_string(true) +
        "'"
        ",CURRENT_TIMESTAMP,SESSION_USER);";
    rc = _sql.execute(command, result);
    if (_verbose > 1) {
      std::cout << command << std::endl;
      std::cout << result << std::endl;
    }
    if (rc != 0) return rc;
  }

  std::stringstream ss;
  ss << "created " << header.name() << " with " << list.grl().size() << " iovs"
     << std::endl;
  _result = ss.str();

  //  if (_dryrun) {
  //    command = "ROLLBACK;";
  //  } else {
  command = "COMMIT;";
  //  }
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  rc = _sql.disconnect();
  if (rc != 0) return rc;

  if (_verbose > 1) std::cout << "returning 0" << std::endl;

  return rc;
}

//**************************************************
int GrlTool::extendList(const std::string& name, const std::string& siov) {
  if (_verbose > 1) std::cout << "running extendList" << std::endl;

  // establish the name is a declared list

  int rc = 0;
  GrlList nlist = list(name);
  const GrlHeader& header = nlist.header();

  if (header.lid() < 0) {
    std::cout << "ERROR - extendList did not find list " << name << std::endl;
    return 1;
  }

  // check the ne wiov doesn't overlap with exiting iovs

  DbIoV niov(siov);
  bool overlap = false;
  for (const auto& iov : nlist.grl()) {
    if (niov.isOverlapping(iov)) {
      std::cout << "New IoV " << niov.to_string(true) << " overlaps old IoV "
                << iov.to_string(true) << std::endl;
      overlap = true;
    }
  }

  if (overlap) return 1;

  // insert in the IoV table, attached to the right list via the lid

  rc = _sql.connect();
  if (rc) {
    std::cout << "ERROR - extendList failed to connect " << std::endl;
    return 3;
  }

  std::string command, result;
  command = "BEGIN";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  command = "SET ROLE grlwrite;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  if (_verbose > 1) std::cout << "committing extend IoV" << std::endl;

  command =
      "INSERT INTO grl.iovs ("
      "lid,iov,create_time,create_user) VALUES ( " +
      std::to_string(header.lid()) + ",'" + niov.to_string(true) +
      "'"
      ",CURRENT_TIMESTAMP,SESSION_USER);";
  rc = _sql.execute(command, result);
  if (_verbose > 1) {
    std::cout << command << std::endl;
    std::cout << result << std::endl;
  }
  if (rc != 0) return rc;

  //  if (_dryrun) {
  //    command = "ROLLBACK;";
  //  } else {
  command = "COMMIT;";
  //  }
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  if (_verbose > 0)
    std::cout << "extended " << header.name() << " with "
              << niov.to_string(true) << std::endl;
  _result.clear();

  rc = _sql.disconnect();
  if (rc != 0) return rc;

  return 0;
}

//**************************************************
int GrlTool::lockList(const std::string& name) {
  if (_verbose > 1) std::cout << "running lockList" << std::endl;
  // establish the name is a declared list

  int rc = 0;
  GrlList nlist = list(name);
  const GrlHeader& header = nlist.header();

  if (header.lid() < 0) {
    std::cout << "ERROR - lockList did not find list " << name << std::endl;
    return 1;
  }

  if (header.locked() > 0) {
    std::cout << "list " << name << " already locked " << std::endl;
    return 1;
  }

  rc = _sql.connect();
  if (rc) {
    std::cout << "ERROR - lockList failed to connect " << std::endl;
    return 3;
  }

  std::string command, result;
  command = "BEGIN";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  command = "SET ROLE grlwrite;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  if (_verbose > 1) std::cout << "setting lock flag" << std::endl;

  command =
      "UPDATE grl.lists "
      "SET locked = 1 WHERE lid = " +
      std::to_string(header.lid()) + ";";
  rc = _sql.execute(command, result);
  if (_verbose > 1) {
    std::cout << command << std::endl;
    std::cout << result << std::endl;
  }
  if (rc != 0) return rc;

  command = "COMMIT;";
  rc = _sql.execute(command, result);
  if (rc != 0) return rc;

  if (_verbose > 1) std::cout << "setting lock flag success" << std::endl;

  rc = _sql.disconnect();
  if (rc != 0) return rc;

  return 0;
}
