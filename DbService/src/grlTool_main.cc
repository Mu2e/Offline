//
// main to drive a command line interface to GrlTool,
// which reads or writes data to the goodrun database
//

#include "Offline/DbService/inc/GrlTool.hh"
#include "Offline/DbService/inc/GrlList.hh"
#include "Offline/GeneralUtilities/inc/ParseCLI.hh"
#include <iomanip>
#include <iostream>

using namespace mu2e;

int main(int argc, char** argv) {
  // Set up ParseCLI with help string
  std::string helpstr = R"(
A program to read and write the goodrun database
  WORDNAMEs are category/subcategory/word and define an integer which
           may have different values for each run
  BITNAMEs are bits within a WORD, like "noise"=0x1

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
)";

  ParseCLI cli(helpstr, true, false);

  // Register global options (empty subcommand "")
  cli.addSubcommand("", "Global options");
  cli.addSwitch("", "verbose", "v", "verbose", true, "set verbose level (0-3)", "0");

  // Register all command subcommands
  // Note: Each subcommand needs at least one switch registered for ParseCLI validation
  // We register a dummy switch with an unlikely name that won't conflict, just to satisfy validation
  cli.addSubcommand("define-word", "define a word");
  cli.addSwitch("define-word", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("print-words", "print all words");
  cli.addSwitch("print-words", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("define-bit", "define a bit");
  cli.addSwitch("define-bit", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("print-bits", "print bits");
  cli.addSwitch("print-bits", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("set-word", "set a word with bits");
  cli.addSwitch("set-word", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("print-entries", "print entries");
  cli.addSwitch("print-entries", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("gen-list", "generate a goodrun list");
  cli.addSwitch("gen-list", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("list-lists", "list all defined GRL lists");
  cli.addSwitch("list-lists", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("create-list", "create a list");
  cli.addSwitch("create-list", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("extend-list", "extend a list");
  cli.addSwitch("extend-list", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("lock-list", "lock a list");
  cli.addSwitch("lock-list", "_dummy", "", "_dummy", false, "", "", false, false);
  cli.addSubcommand("print-list", "print a list");
  cli.addSwitch("print-list", "_dummy", "", "_dummy", false, "", "", false, false);

  // Parse command line arguments
  int rc = cli.setArgs(argc, argv);
  if (rc != 0) {
    if (rc == 999) {
      // Help was printed by autohelp
      return 1;
    }
    return rc;
  }

  // Get verbose level
  int verbose = cli.getInt("", "verbose");

  // Get subcommand (command name)
  std::string command = cli.subcommand();
  if (command.empty()) {
    std::cout << "Please specify a command\n";
    return 1;
  }

  // Get positional arguments (command arguments)
  const StringVec& positionals = cli.positionals();
  std::string arg1, arg2, arg3;
  if (positionals.size() > 0) arg1 = positionals[0];
  if (positionals.size() > 1) arg2 = positionals[1];
  if (positionals.size() > 2) arg3 = positionals[2];

  GrlTool tool;
  tool.setVerbose(verbose);

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
    for (size_t i=2; i<positionals.size(); i++) {
      bitnames.push_back(positionals[i]);
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
      std::cout << "print-list requires a list name as first argument"<< std::endl;
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
