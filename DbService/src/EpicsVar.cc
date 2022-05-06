#include "Offline/DbService/inc/EpicsVar.hh"
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

void EpicsVar::set(std::string const& csv) {
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

  _channel_id = std::stoll(sv[0]);
  _smpl_time_s = sv[1];
  _smpl_time_t = 0;
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

  std::string _str_val = sv[7];

  if (!sv[9].empty()) {
    _datatype = sv[9][0];
  } else {
    _datatype = ' ';
  }

  return;
}
