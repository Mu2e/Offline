//
// main to drive a command line interface to GrlTool,
// which reads or writes data to the goodrun database
//

#include "Offline/DbService/inc/GrlList.hh"
#include "Offline/DbService/inc/GrlTool.hh"
#include "Offline/GeneralUtilities/inc/ParseCLI.hh"
#include <iomanip>
#include <iostream>

using namespace mu2e;

int main(int argc, char** argv) {
  // Set up ParseCLI with brief help string
  // Detailed help is provided in each command's subcommand definition
  ParseCLI cli(R"(A program to read and write the goodrun database
   The first main interface is setting and reading bits within words.
   Bit settings may be superseded by making a new entry.
   WORDNAMEs are written as category/subcategory/word and define an integer
         which defines quality, and may have different values for each run
   BITNAMEs are bits within a WORD, like "noise" might be 0x1

   The second interface is to use the word and bit settings to create
   a list of runs. The lists can be extendable or locked. List can also be provided
   explicitly instead of dervied from bits)");

  // Register global options (empty subcommand "")
  cli.addSubcommand("", "Global options");
  cli.addSwitch("", "verbose", "v", "verbose", true, "set verbose level (0-3)",
                "0");

  // Register all command subcommands with proper help strings and switches
  cli.addSubcommand("define-word", "define a new word");
  cli.addSwitch("define-word", "wordname", "w", "wordname", true,
                "WORDNAME\n        word name (\"category/subcategory/word\")",
                "", false, true);
  cli.addSwitch("define-word", "description", "d", "description", true,
                "\"text\"\n        optional description", "", false, false);

  cli.addSubcommand("print-words", "print all defined words");

  cli.addSubcommand("define-bit",
                    "define a new bit associated to a given word\n   ex.: "
                    "define-bit -w shift/dqm/trk -b noise");
  cli.addSwitch("define-bit", "wordname", "w", "wordname", true,
                "WORDNAME\n        word name, (\"category/subcat/word\")", "",
                false, true);
  cli.addSwitch("define-bit", "bitname", "b", "bitname", true,
                "BITNAME\n        bit name within the word", "", false, true);
  cli.addSwitch("define-bit", "description", "d", "description", true,
                "\"text\"\n        optional description", "", false, false);

  cli.addSubcommand("print-bits", "print defined bits");
  cli.addSwitch("print-bits", "wordname", "", "wordname", true,
                "optional WORDNAME", "", false, false);

  cli.addSubcommand(
      "set-word",
      "set a word with bits to indicate data quality"
      "\n     Later entries can override earlier ones. For run range text "
      "syntax see "
      "\n     "
      "https://mu2ewiki.fnal.gov/wiki/ConditionsData#Intervals_of_validity"
      "\n     # set this word to 0 (nothing bad), for one run, or run range"
      "\n     ex.: set-bits -w shift/dqm/trk -r 10500"
      "\n     ex.: set-bits -w shift/dqm/trk -r 10500-105010"
      "\n     # set this word, overriding previous word, only for these subruns"
      "\n     ex.: set-bits -w shift/dqm/trk -r 10500:50-10500:52 -b noise -b "
      "trips");
  cli.addSwitch("set-word", "wordname", "w", "wordname", true,
                "WORDNAME\n        word name, (\"category/subcat/word\")", "",
                false, true);
  cli.addSwitch("set-word", "runrange", "r", "runrange", true,
                "RUNRANGE\n        run range, like 10500-105010", "", false,
                true);
  cli.addSwitch("set-word", "bitname", "b", "bitname", true,
                "BITNAME\n        bit name within the word", "", true, false);

  cli.addSubcommand("print-entries", "dump relevant database entries");
  cli.addSwitch("print-entries", "wordname", "w", "wordname", true,
                "WORDNAME\n        word name, (\"category/subcat/word\")", "",
                false, false);
  cli.addSwitch("print-entries", "runrange", "r", "runrange", true,
                "RUNRANGE\n        run range, like 10500-105010", "", false,
                false);

  cli.addSubcommand(
      "gen-list",
      "apply list of bits selections to generate a goodrun list"
      "\n   bit selections are written in a text file, see syntax at"
      "\n   https://mu2ewiki.fnal.gov/wiki/Goodruns"
      "\n   optionaally only allow runs of a given type");
  cli.addSwitch("gen-list", "commandfile", "c", "commandfile", true,
                "COMMANDFILE\n        filespec of text file of bit selections",
                "", false, true);
  cli.addSwitch("gen-list", "runtype", "r", "runtype", true,
                "RUNTYPE"
                "\n        integer run type from run_info database",
                "", false, false);

  cli.addSubcommand("list-lists", "list all defined GRL lists");

  cli.addSubcommand("create-list",
                    "create a new goodrun list"
                    "\n   optionally provide an initial set of runs to add to the list");
  cli.addSwitch("create-list", "listname", "l", "listname", true, "LISTNAME"
                "\n        a name", "",
                false, true);
  cli.addSwitch("create-list", "filename", "f", "filename", true,
                "FILENAME\n        file containing run ranges", "", false, false);

  cli.addSubcommand("extend-list",
                    "add a run range to a list");
  cli.addSwitch("extend-list", "listname", "l", "listname", true, "LISTNAME"
                "\n        the name of the list", "",
                false, true);
  cli.addSwitch("extend-list", "runrange", "r", "runrange", true, "RUNRANGE"
                "\n        run range, like 10500-105010", "", false, true);

  cli.addSubcommand("lock-list",
                    "prevent any new IoV's being added to the GRL");
  cli.addSwitch("lock-list", "listname", "l", "listname", true, "LISTNAME"
                "\n        the name of the list", "",
                false, true);

  cli.addSubcommand("print-list", "print the IoV's in the list");
  cli.addSwitch("print-list", "listname", "l", "listname", true, "LISTNAME"
                "\n        the name of the list", "",
                false, true);

  // Parse command line arguments
  int rc = cli.setArgs(argc, argv);
  if (rc != 0) {
    return (rc == 999) ? 1 : rc;  // 999 means help was printed
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

  if (command == "define-word") {
    std::string wordname = cli.getString("define-word", "wordname");
    std::string description = cli.getString("define-word", "description");
    tool.defineWord(wordname, description);
  } else if (command == "print-words") {
    tool.words();
  } else if (command == "define-bit") {
    std::string wordname = cli.getString("define-bit", "wordname");
    std::string bitname = cli.getString("define-bit", "bitname");
    std::string description = cli.getString("define-bit", "description");
    tool.defineBit(wordname, bitname, description);
  } else if (command == "print-bits") {
    std::string wordname = cli.getString("print-bits", "wordname");
    tool.bits(wordname);
  } else if (command == "set-word") {
    std::string wordname = cli.getString("set-word", "wordname");
    std::string runrange = cli.getString("set-word", "runrange");
    StringVec bitnames = cli.getStrings("set-word", "bitname");
    tool.setBits(wordname, runrange, bitnames);
  } else if (command == "print-entries") {
    std::string wordname = cli.getString("print-entries", "wordname");
    std::string bitname = cli.getString("print-entries", "runrange");
    tool.entries(wordname, bitname);
  } else if (command == "gen-list") {
    std::string commandfile = cli.getString("gen-list", "commandfile");
    std::string runtype = cli.getString("gen-list", "runtype");
    tool.genList(commandfile, runtype);
  } else if (command == "list-lists") {
    tool.lists();
  } else if (command == "create-list") {
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
  } else if (command == "extend-list") {
    std::string listname = cli.getString("extend-list", "listname");
    std::string runrange = cli.getString("extend-list", "runrange");
    tool.extendList(listname, runrange);
  } else if (command == "lock-list") {
    std::string listname = cli.getString("lock-list", "listname");
    tool.lockList(listname);
  } else if (command == "print-list") {
    std::string listname = cli.getString("print-list", "listname");
    GrlList ll = tool.list(listname);
    // ll.print(std::cout);
  } else {
    std::cout << "unknown command: " << command << std::endl;
    return 1;
  }
  std::cout << tool.result() << std::flush;

  return 0;
}
