//
// main to drive a command line interface to EpicsTool,
// which reads data from the EPICS archive database
//

#include "Offline/DbService/inc/EpicsTool.hh"
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  std::vector<std::string> words;

  bool doHelp = false;
  bool doNames = false;
  std::string word, var, time;
  size_t i = 1;
  size_t lim(argc);
  while (i < lim) {
    // for (size_t i = 1; i < size_t(argc); ++i) {
    word = argv[i];
    if (word == "-h" || word == "--help") {
      doHelp = true;
    } else if (word == "--names") {
      doNames = true;
    } else if (word == "--print") {
      if (i + 1 == lim) {
        std::cout << "Error --print doesn't have an argument" << std::endl;
        doHelp = true;
      } else {
        var = argv[++i];
      }
    } else if (word == "--time") {
      if (i + 1 == lim) {
        std::cout << "Error --time doesn't have an argument" << std::endl;
        doHelp = true;
      } else {
        time = argv[++i];
      }
    } else {
      std::cout << "Error - unknown argument " << word << std::endl;
      doHelp = true;
    }
    i++;
  }

  if (doHelp) {
    std::cout <<
        R"(

  epicsTool [OPTIONS]

  A program to read EPICS variable data from the EPICS offline
  replica database and print it out.

  --names
       Print the list of variables, which can be grepped or
       examined for the one you want
  --print NAME
       Print the all the values for the variable with this
       name, can be restricted with --time
  --time TIME
       The time restriction in ISO8601 format. This can
       be one time in which case the variable point closes to the
       time will be printed, if it is a time range, then values in
       that range will be printed. If the time argument is a int or
       float number, thee results returned will be from that many
       days ago to now.

  epicsTool --names
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp \
            --time  2022-01-01T10:11:12
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp \
            --time  2022-01-01T10:11:12/2022-01-02T10:11:12
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp \
            --time 1.5

Variable prints will have the following format:
channel_id,smpl_time,nanosecs,severity_id,status_id,num_val,float_val,str_val,datatype

)" << std::endl;
    return 0;
  }

  mu2e::EpicsTool tool;
  int rc = 0;

  if (doNames) {
    mu2e::EpicsTool::StringVec names;
    rc = tool.names(names);
    if (rc != 0) return rc;
    for (auto const& n : names) {
      std::cout << n << std::endl;
    }
    return 0;
  } else {
    mu2e::EpicsVar::EpicsVec vec = tool.get(var, time);
    for (auto const& e : vec) {
      std::cout << e.csv() << std::endl;
    }
  }

  return 0;
}
