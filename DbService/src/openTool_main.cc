//
// main to drive a command line interface to OpenTool,
// which reads or writes open intervals for the ancillary iov system
//

#include "Offline/DbService/inc/OpenTool.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/GeneralUtilities/inc/splitString.hh"
#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>

using namespace mu2e;
namespace po = boost::program_options;

int main(int argc, char** argv) {
  int verbose = 0;
  std::string database{"mu2e_conditions_prd"}, filespec, runstr, iovstr, name,
      comment;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "describe arguments")(
      "verbose,v", po::value<int>(&verbose), "set verbose level (0-3)")(
      "database,d", po::value<std::string>(&database), "database name")(
      "run,r", po::value<std::string>(&runstr), "run or run:subrun")(
      "iov,i", po::value<std::string>(&iovstr), "iov in std format")(
      "file,f", po::value<std::string>(&filespec), "file spec")(
      "comment,c", po::value<std::string>(&comment), "comment")(
      "name,n", po::value<std::string>(&name), "table name")(
      "admin,a", "write through admin permissions");

  po::variables_map vm;
  po::parsed_options parsed = po::parse_command_line(argc, argv, desc);
  po::store(parsed, vm);
  po::notify(vm);

  std::vector<std::string> unregistered =
      po::collect_unrecognized(parsed.options, po::include_positional);

  if (vm.count("help") && unregistered.empty()) {
    std::cout <<
        R"(

  openTool [GLOBAL OPTIONS] command [OPTIONS]

  GLOBAL OPTIONS
  -h, --help   print help
  -d, --database sterr to database
  -v, --verbose INT increase verboseness

  commit-calibration
      Read a file containing a calibration, commmit, and mark with a given IoV.
      The file must be in the standard format, which include the table name.
      The Iov is in the stand text format
     -f, --file FILENAME  the calibration file name
     -i, --iov  IOV the interval in standard format
     -c, --comment "TEXT" ex. "a commit comment"
      Example:
        openTool commit-calibration --file gains.txt --iov 105000-MAX
  get-calibration
      Print the calibration table content for the given run
     -n,--name TABLENAME the table name
     -r, --run  RUN:SUBRUN  may be RUN or RUN:SBRUN
      Example:
        openTool get-calibration -n MyTable -r 105000
  print-iovs
      Print all IoVs for a given table
     -n, --name TABLENAME  (default is all tables)

)" << std::endl;

    return 1;
  }

  if (unregistered.empty()) {
    std::cout << "Please specify a command\n";
    return 1;
  }

  OpenTool tool(database);
  tool.setVerbose(verbose);

  if (vm.count("admin")) {
    tool.setAdmin(true);
  }

  std::string command = unregistered[0];

  int rc;
  if (command == "commit-calibration") {
    DbIoV iov(iovstr);
    rc = tool.commit(filespec, iov, comment);
    if (rc) return rc;
  } else if (command == "get-calibration") {
    uint32_t run{0}, subrun{0};
    StringVec words = splitString(runstr, ":", "\"", "\\", true, false);
    run = std::stoi(words[0]);
    if (words.size() > 1) subrun = std::stoi(words[1]);
    std::string csv;
    int cid;
    DbIoV iov;
    rc = tool.table(name, run, subrun, csv, cid, iov);
    if (rc) return rc;
    std::cout << "TABLE " << name << "\n";
    std::cout << csv;
  } else if (command == "print-iovs") {
    rc = tool.readIoVs(name);
    if (rc) return rc;
    const ValOpenIovs& iovs = tool.iovs();
    std::cout << iovs.csv();
  } else {
    std::cout << "unknown command: " << command << std::endl;
    return 1;
  }

  return 0;
}
