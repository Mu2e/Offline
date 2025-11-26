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

  // Register all command subcommands with proper help strings and switches
  cli.addSubcommand("define-word", "define a word\n    ex.: define-word shift/dqm/trk");
  cli.addSwitch("define-word", "wordname", "", "wordname", true, "WORDNAME (category/subcategory/word)", "", false, true);
  cli.addSwitch("define-word", "description", "", "description", true, "optional description", "", false, false);

  cli.addSubcommand("print-words", "print all words");

  cli.addSubcommand("define-bit", "define a new bit associated to the given word\n    ex.: define-bit shift/dqm/trk noise");
  cli.addSwitch("define-bit", "wordname", "", "wordname", true, "WORDNAME", "", false, true);
  cli.addSwitch("define-bit", "bitname", "", "bitname", true, "BITNAME", "", false, true);
  cli.addSwitch("define-bit", "description", "", "description", true, "optional description", "", false, false);

  cli.addSubcommand("print-bits", "print bits");
  cli.addSwitch("print-bits", "wordname", "", "wordname", true, "optional WORDNAME", "", false, false);

  cli.addSubcommand("set-word", "set a word with bits to indicate problems\n    later entries can override earlier ones\n    for run range text syntax see https://mu2ewiki.fnal.gov/wiki/ConditionsData#Intervals_of_validity\n    ex.: set-word --wordname shift/dqm/trk --runrange 10500 --bitname noise --bitname trips");
  cli.addSwitch("set-word", "wordname", "", "wordname", true, "WORDNAME", "", false, true);
  cli.addSwitch("set-word", "runrange", "", "runrange", true, "RUNRANGE", "", false, true);
  cli.addSwitch("set-word", "bitname", "", "bitname", true, "BITNAME (can be repeated)", "", true, true);

  cli.addSubcommand("print-entries", "print entries");
  cli.addSwitch("print-entries", "wordname", "", "wordname", true, "optional WORDNAME", "", false, false);
  cli.addSwitch("print-entries", "runrange", "", "runrange", true, "optional RUNRANGE", "", false, false);

  cli.addSubcommand("gen-list", "apply list of bits selections in command file to runtype runs\n    and generate a goodrun list");
  cli.addSwitch("gen-list", "commandfile", "", "commandfile", true, "COMMANDFILE", "", false, true);
  cli.addSwitch("gen-list", "runtype", "", "runtype", true, "optional RUNTYPE", "", false, false);

  cli.addSubcommand("list-lists", "list all defined GRL lists");

  cli.addSubcommand("create-list", "define a list, if FILENAME is provided, then add those IoV's into the list");
  cli.addSwitch("create-list", "listname", "", "listname", true, "LISTNAME", "", false, true);
  cli.addSwitch("create-list", "filename", "", "filename", true, "optional FILENAME with IoVs", "", false, false);

  cli.addSubcommand("extend-list", "add an IoV to a list. IoV in standard text format");
  cli.addSwitch("extend-list", "listname", "", "listname", true, "LISTNAME", "", false, true);
  cli.addSwitch("extend-list", "iov", "", "iov", true, "IOV", "", false, true);

  cli.addSubcommand("lock-list", "prevent any new IoV's being added to the GRL");
  cli.addSwitch("lock-list", "listname", "", "listname", true, "LISTNAME", "", false, true);

  cli.addSubcommand("print-list", "print the IoV's in the list");
  cli.addSwitch("print-list", "listname", "", "listname", true, "LISTNAME", "", false, true);

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

  GrlTool tool;
  tool.setVerbose(verbose);

  if( command == "define-word" ) {
    std::string wordname = cli.getString("define-word", "wordname");
    std::string description = cli.getString("define-word", "description");
    tool.defineWord(wordname, description);
  } else if( command == "print-words" ) {
    tool.words();
  } else if( command == "define-bit" ) {
    std::string wordname = cli.getString("define-bit", "wordname");
    std::string bitname = cli.getString("define-bit", "bitname");
    std::string description = cli.getString("define-bit", "description");
    tool.defineBit(wordname, bitname, description);
  } else if( command == "print-bits" ) {
    std::string wordname = cli.getString("print-bits", "wordname");
    tool.bits(wordname);
  } else if( command == "set-word" ) {
    std::string wordname = cli.getString("set-word", "wordname");
    std::string runrange = cli.getString("set-word", "runrange");
    StringVec bitnames = cli.getStrings("set-word", "bitname");
    tool.setBits(wordname, runrange, bitnames);
  } else if( command == "print-entries" ) {
    tool.entries();
  } else if( command == "gen-list" ) {
    std::string commandfile = cli.getString("gen-list", "commandfile");
    std::string runtype = cli.getString("gen-list", "runtype");
    tool.genList(commandfile, runtype);
  } else if( command == "list-lists" ) {
    tool.lists();
  } else if( command == "create-list" ) {
    std::string listname = cli.getString("create-list", "listname");
    std::string filename = cli.getString("create-list", "filename");
    GrlHeader hh(listname);
    if (filename.empty()) {
      GrlList newlist(hh);
      tool.createList(newlist);
    } else {
      GrlList newlist(hh, filename);
      tool.createList(newlist);
    }
  } else if( command == "extend-list" ) {
    std::string listname = cli.getString("extend-list", "listname");
    std::string iov = cli.getString("extend-list", "iov");
    tool.extendList(listname, iov);
  } else if( command == "lock-list" ) {
    std::string listname = cli.getString("lock-list", "listname");
    tool.lockList(listname);
  } else if( command == "print-list" ) {
    std::string listname = cli.getString("print-list", "listname");
    if (listname.empty()) {
      std::cout << "print-list requires a list name"<< std::endl;
      return 1;
    }
    GrlList ll = tool.list(listname);
    //ll.print(std::cout);
  } else {
    std::cout << "unknown command: " << command << std::endl;
    return 1;
  }
  std::cout << tool.result() <<std::flush;

  return 0;
}
