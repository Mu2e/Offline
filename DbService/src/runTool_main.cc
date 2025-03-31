//
// main to drive a command line interface to RunTool,
// which reads data from the run conditions database
//

#include "Offline/DbService/inc/RunTool.hh"
#include "Offline/DbService/inc/RunSelect.hh"
#include <boost/program_options.hpp>
#include <iostream>
#include <iomanip>

using namespace mu2e;
namespace po = boost::program_options;

int main(int argc, char** argv) {
  std::vector<std::string> words;

  bool flags = false, conditions = false, longp = false, transitions = false;
  std::string run, time, type;
  int last = 0, days = 0;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "describe arguments")(
      "flags,f", "print types of runs, ..")(
      "run,r", po::value(&run), "run range")(
      "last,n", po::value(&last), "print last N runs")(
      "type,y", po::value(&type), "comma-seperated list of run types")(
      "time,t", po::value(&time), "Restrict times for print, see help")(
      "days,d", po::value(&days), "restrict time to N days ago to now")(
      "long,l", "print longer report")(
      "conditions,c","also print conditions")(
      "transitions,a","also print transitions");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout <<
        R"(

  runTool [OPTIONS]

  A program to read run configuration data from the run_info
  replica database and print it out.

  -f, --flags
       print list of run, transition, and cause flags

  run selection:
  -r, --run=RUNRANGE
       Select a run range, like
       100000 (one run)
       100000-101000 (range)
  -n, --last=N
       Select N most recent runs
  -y, --type=A,B
       Select only runs of these types (see --flags)
  -t, --time=TIME1 or TIME1/TIME2
       The time restriction in ISO8601 format. This can be one time
       to select all runs since that time, or two times for a range
  -d, --days=N
       all runs since N days ago

    extend print:
  -l, --long
       Also print the details
  -c, --conditions
       Also print the conditions text block for selected runs
  -a, --transitions
       Also print the transitions

    Runs default print will have format:
   run   run_type   commit_time
 103025    4    2024-06-12 23:08:36
    Transitions will print in format
transition_flag, cause_flag, time
5,0,2024-06-06 17:27:42.666540-05:000

)" << std::endl;

    return 1;
  }

  if (vm.count("flags")) flags = true;

  RunTool tool;
  int rc = 0;

  if (flags) {
    std::map<int,std::string> flags;
    flags = tool.flags(RunTool::FlagType::run);
    std::cout << "run flags:"<<std::endl;
    for (auto pp : flags) {
      std::cout << std::setw(4)<<pp.first<< "  " <<pp.second <<std::endl;
    }

    flags = tool.flags(RunTool::FlagType::transition);
    std::cout << "transition flags:"<<std::endl;
    for (auto pp : flags) {
      std::cout << std::setw(4)<<pp.first<< "  " <<pp.second <<std::endl;
    }

    flags = tool.flags(RunTool::FlagType::cause);
    std::cout << "cause flags:"<<std::endl;
    for (auto pp : flags) {
      std::cout << std::setw(4)<<pp.first<< "  " <<pp.second <<std::endl;
    }

    return rc;
  }


  RunSelect runsel(run,last,type,time,days);

  if (vm.count("long")) longp = true;
  if (vm.count("conditions")) conditions = true;
  if (vm.count("transitions")) transitions = true;

  RunInfo::RunVec runvec = tool.listRuns(runsel,conditions,transitions);
  for(auto const& rr : runvec) {
    tool.printRun(rr,longp);
  }

  return rc;
}
