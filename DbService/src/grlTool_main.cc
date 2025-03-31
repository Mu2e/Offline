//
// main to drive a command line interface to GrlTool,
// which reads or writes data to the goodrun database
//

#include "Offline/DbService/inc/GrlTool.hh"
#include "Offline/DbService/inc/GrlList.hh"
#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>

using namespace mu2e;
namespace po = boost::program_options;

int main(int argc, char** argv) {
  // std::vector<std::string> words;

  // bool flags = false, conditions = false, longp = false, transitions = false;
  // std::string run, time, type;
  // int last = 0, days = 0;
  int verbose = 0;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "describe arguments")
    ("verbose,v", po::value<int>(&verbose), "set verbose level (0-3)");
    //  ("verbose,v", po::value<std::vector<int>(), "Verbosity level (can be repeated)");
  // ("verbose,v",  "more verbose");

  //(
  //    "run,r", po::value(&run), "run range")(
  //    "last,n", po::value(&last), "print last N runs")(
  //    "type,y", po::value(&type), "comma-seperated list of run types")(
  //    "time,t", po::value(&time), "Restrict times for print, see help")(
  //    "days,d", po::value(&days), "restrict time to N days ago to now")(
  //    "long,l", "print longer report")(
  //    "conditions,c","also print conditions")(
  //    "transitions,a","also print transitions");

  po::variables_map vm;
  po::parsed_options parsed = po::parse_command_line(argc, argv, desc);
  po::store(parsed, vm);
  po::notify(vm);

  std::vector<std::string> unregistered =
      po::collect_unrecognized(parsed.options, po::include_positional);

  if (vm.count("help") && unregistered.empty()) {
    std::cout <<
        R"(

  grTool [GLOBAL OPTIONS] command [OPTIONS]

  A program to read and write the goodrun database
    WORDNAMEs are category/subcategory/word and define an integer which
             may have different values for each run
    BITNAMEs are bits within a WORD, like "noise"=0x1

  -h, --help   print help

  **** bits ****

  define-word WORDNAME ["description"]
      ex.: define-word shift/dqm/trk
  print-words

  define-bit WORDNAME BITNAME ["description"]
      define a new bit associated to the given word
      ex.: define-bit shift/dqm/trk noise
  print-bits [WORDNAME]

  set-word WORDNAME RUNRANGE BITNAME [BITNAME]
      set a word with bits to indicate problems
      later entries can override earlier ones
      for run range text syntax see
         https://mu2ewiki.fnal.gov/wiki/ConditionsData#Intervals_of_validity
      # set this word to 0 (nothing bad), for one run, or run range
      ex.: set-bits shift/dqm/trk 10500
      ex.: set-bits shift/dqm/trk 10500-105010
      # set this word, overriding previous word, only for these subruns
      ex.: set-bits shift/dqm/trk 10500:50-10500:52 noise trips baddqm

  print-entries [WORDNAME] [RUNRANGE]

  gen-list COMMANDFILE [RUNTYPE]
      apply list of bits selections in command file to runtype runs
      and generate a goodrun list

  **** lists ****

  list-lists
      list all defined GRL lists
  create-list LISTNAME [FILENAME]
      define a list, if FILENAME is provided, then add those IoV's into the list
  extend-list LISTNAME IOV
      add an IoV to a list. IoV in standard text format
  lock-list LISTNAME
      prevent any new IoV's being added to the GRL
  print-list LISTNAME
      print the IoV's in the list

)" << std::endl;

    return 1;
  }

  if (unregistered.empty()) {
    std::cout << "Please specify a command\n";
    return 1;
  }

  GrlTool tool;
  tool.setVerbose(verbose);

  std::string command = unregistered[0];
  std::string arg1,arg2,arg3;
  if (unregistered.size()>1) arg1 = unregistered[1];
  if (unregistered.size()>2) arg2 = unregistered[2];
  if (unregistered.size()>3) arg3 = unregistered[3];

  if( command == "define-word" ) {
    tool.defineWord(arg1,arg2);
  } else if( command == "print-words" ) {
    tool.words();
  } else if( command == "define-bit" ) {
    tool.defineBit(arg1,arg2,arg3);
  } else if( command == "print-bits" ) {
    tool.bits(arg1);
  } else if( command == "set-word" ) {
    // any number of bit name, concatenate with commas
    StringVec bitnames;
    for (size_t i=3; i<unregistered.size(); i++) {
      bitnames.push_back(unregistered[i]);
    }
    tool.setBits(arg1,arg2,bitnames);
  } else if( command == "print-entries" ) {
    tool.entries();
  } else if( command == "gen-list" ) {
    tool.genList(arg1,arg2);
  } else if( command == "list-lists" ) {
    tool.lists();
  } else if( command == "create-list" ) {
    GrlHeader hh(arg1); // arg1=list name
    if (arg2.empty()) { // arg2 = file with IoVs
      GrlList newlist(hh);
      tool.createList(newlist);
    } else {
      GrlList newlist(hh,arg2);
      tool.createList(newlist);
    }
  } else if( command == "extend-list" ) {
    tool.extendList(arg1,arg2); // name, iov
  } else if( command == "lock-list" ) {
    tool.lockList(arg1); // name
  } else if( command == "print-list" ) {
    if (arg1.empty()) {
      std::cout << "print- list requires a list name as first orgument"<< std::endl;
    } else {
      GrlList ll = tool.list(arg1);
      //ll.print(std::cout);
    }
  } else {
    std::cout << "unknown command: " << command << std::endl;
    return 1;
  }
  std::cout << tool.result() <<std::flush;

  return 0;
}
