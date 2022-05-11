#include "Offline/DbService/inc/EpicsVar.hh"
#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <string>
#include <vector>

using namespace mu2e;
using namespace boost;

//**************************************************
// Example csv:
// 642,2022-05-05 14:10:02.089404-05:00,89404264,5,9,None,55.0,None, ,None
//
EpicsVar::EpicsVar(std::string const& csv) {
  _csv = csv;

  std::vector<std::string> sv;

  // split into words at commas
  escaped_list_separator<char> els('\\', ',', '"');
  tokenizer<escaped_list_separator<char>> tok(csv, els);
  for (tokenizer<escaped_list_separator<char>>::iterator beg = tok.begin();
       beg != tok.end(); ++beg) {
    std::string ss(*beg);
    boost::trim(ss);
    sv.emplace_back(ss);
  }

  if(sv.size()!=10) {
    throw cet::exception("EPICSVAR_BAD_CSV")
      << " EpicsVar::ctor number of csv fields !=10 : \n   " << csv << " \n";
  }

  _channel_id = std::stoll(sv[0]);
  _smpl_time_s = sv[1];
  if(TimeUtility::parseTimeTZ(sv[1],_smpl_time_t)) {
    throw cet::exception("EPICSVAR_BAD_TIME")
      << " EpicsVar::ctor could not parse time " << sv[1] << " \n";
  }
  _nanosecs = std::stoll(sv[2]);
  _severity_id = std::stoll(sv[3]);
  _status_id = std::stoll(sv[4]);

  if (sv[5] == "None") {
    _num_val = 0;
  } else {
    _num_val = std::stoi(sv[5]);
  }

  if (sv[6] == "None") {
    _float_val = 0.0;
  } else {
    _float_val = std::stod(sv[6]);
  }

  if (sv[7] == "None") {
    _str_val = "";
  } else {
    _str_val = sv[7];
  }

  if (!sv[8].empty()) {
    _datatype = sv[8][0];
  } else {
    _datatype = ' ';
  }

  return;
}
