//
// main to drive a command line interface to EpicsTool,
// which reads data from the EPICS archive database
//

#include "Offline/DbService/inc/EpicsTool.hh"
#include <boost/program_options.hpp>
#include <iostream>

using namespace mu2e;
namespace po = boost::program_options;

int main(int argc, char** argv) {
  std::vector<std::string> words;

  // bool doHelp = false;
  bool doNames = false;
  std::string word, name, time;
  float daysAgo = 0.0;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "describe arguments")(
      "names,n", po::bool_switch(&doNames), "print names of all variables")(
      "print,p", po::value(&name), "print same for this variable name")(
      "time,t", po::value(&time), "Restrict times for print, see help")(
      "daysAgo,d", po::value(&daysAgo), "restrict time to dayaAgo to now");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
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
  --time TIME or TIME1/TIME2
       The time restriction in ISO8601 format. This can
       be one time in which case the variable point closes to the
       time will be printed, if it is a time range, then values in
       that range will be printed. If time zone is not specified,
       time is interpreted as UTC.
  --daysAgo NUMBER
       Number is an int or float, and the results returned will be
       from that many days ago to now.

  epicsTool --names
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp \
            --time  2022-01-01T10:11:12-0500
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp \
            --time  2022-01-01T10:11:12/2022-01-02T10:11:12
  epicsTool --print Mu2e:CompStatus:calo-01:DTC0:dtctemp \
            --daysAgo 1.5

Variable sample prints will have the following format:
channel_id,smpl_time,nanosecs,severity_id,status_id,num_val,float_val,str_val,datatype

)" << std::endl;

    return 1;
  }


  EpicsTool tool;
  int rc = 0;

  if (doNames) {
    StringVec names;
    rc = tool.names(names);
    if (rc != 0) return rc;
    for (auto const& name : names) {
      std::cout << name << std::endl;
    }
    return 0;
  } else if (!name.empty()) {
    EpicsVar::EpicsVec vec = tool.get(name, time, daysAgo);
    for (auto const& e : vec) {
      std::cout << e.csv() << std::endl;
    }
  }

  return 0;
}
