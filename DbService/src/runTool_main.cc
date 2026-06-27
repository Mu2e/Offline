//
// main to drive a command line interface to RunTool,
// which reads data from the run conditions database
//

#include "Offline/DbService/inc/RunSelect.hh"
#include "Offline/DbService/inc/RunTool.hh"
#include "Offline/GeneralUtilities/inc/ParseCLI.hh"

#include <iomanip>
#include <iostream>

using namespace mu2e;

int main(int argc, char** argv) {
  ParseCLI cli("Run database query tool");

  // Add the global subcommand (no subcommands in this tool)
  cli.addSubcommand("", "");

  // Add switches
  cli.addSwitch("", "flags", "f", "flags", false,
                "print run and transition types");
  cli.addSwitch("", "run", "r", "run", true,
                "run range (e.g., 100000 or 100000-101000)");
  cli.addSwitch("", "last", "n", "last", true, "last N runs");
  cli.addSwitch("", "type", "y", "type", true,
                "comma-separated list of run types");
  cli.addSwitch("", "time", "t", "time", true,
                "time restriction (ISO8601 format)\n         since: "
                "2026-06-05T18:38:20-05:00 or range: 2026-06-01/2026-06-05");
  cli.addSwitch("", "days", "d", "days", true, "since N days ago to now");
  cli.addSwitch("", "configs", "c", "configs", false,
                "also print configuration records");
  cli.addSwitch("", "transitions", "a", "transitions", false,
                "also print transitions");
  cli.addSwitch("", "subruns", "s", "subruns", false, "also print subruns");
  cli.addSwitch("", "blob", "b", "blob", true,
                "print formatted JSON blob for specified subsystem");
  cli.addSwitch("", "dbtables", "q", "dbtables", false,
                "print the cat 3 DbService tables");

  // Parse command line
  int rc = cli.setArgs(argc, argv);
  if (rc != 0) {
    return rc;
  }

  // Check for help (autohelp handles this, but we return after)
  if (cli.getBool("", "help")) {
    return 0;
  }

  RunTool tool;

  // Handle flags request
  if (cli.getBool("", "flags")) {
    std::map<int, std::string> flags;

    flags = tool.flags(RunTool::FlagType::run);
    std::cout << "run flags:" << std::endl;
    for (auto pp : flags) {
      std::cout << std::setw(4) << pp.first << "  " << pp.second << std::endl;
    }

    flags = tool.flags(RunTool::FlagType::transition);
    std::cout << "transition flags:" << std::endl;
    for (auto pp : flags) {
      std::cout << std::setw(4) << pp.first << "  " << pp.second << std::endl;
    }

    return 0;
  }

  // Build RunSelect from command line arguments
  std::string run = cli.getString("", "run");
  std::string last_str = cli.getString("", "last");
  int last = last_str.empty() ? 0 : std::stoi(last_str);
  std::string type = cli.getString("", "type");
  std::string time = cli.getString("", "time");
  std::string days_str = cli.getString("", "days");
  int days = days_str.empty() ? 0 : std::stoi(days_str);

  RunSelect runsel(run, last, type, time, days);

  // Get print options
  bool configs = cli.getBool("", "configs");
  bool transitions = cli.getBool("", "transitions");
  bool subruns = cli.getBool("", "subruns");
  std::string blob_subsystem = cli.getString("", "blob");
  bool dbtables = cli.getBool("", "dbtables");

  // If blob subsystem or dbtables is specified, we need to fetch configs
  bool need_configs = configs || !blob_subsystem.empty() || dbtables;

  // Query and print runs
  RunInfo::RunVec runvec =
      tool.listRuns(runsel, need_configs, transitions, subruns);

  // Handle dbtables output: concatenated cat-3 DbService table values
  // for machine reading.  No labels or adornment - just the values,
  // looping over every config in every run.
  if (dbtables) {
    for (auto const& rr : runvec) {
      for (const auto& config : rr.configs()) {
        std::string tables = config.dbTables3(false);
        if (!tables.empty()) {
          std::cout << tables << "\n";
        }
      }
    }
    return 0;
  }

  // Handle JSON blob output for specific subsystem
  if (!blob_subsystem.empty()) {
    for (auto const& rr : runvec) {
      std::cout << "Run " << rr.runNumber() << ":\n";
      bool found = false;
      for (const auto& config : rr.configs()) {
        if (config.subsystem() == blob_subsystem) {
          found = true;
          // Get the settings string (already cleaned in RunTool.cc)
          std::string settings = config.settings();

          // Replace literal \n with actual newlines
          size_t pos = 0;
          while ((pos = settings.find("\\n", pos)) != std::string::npos) {
            settings.replace(pos, 2, "\n");
            pos += 1;
          }

          // Simple JSON pretty-printing: add indentation after { [ and newlines
          // after , ; :
          std::string formatted;
          int indent = 0;
          bool in_string = false;
          char prev_char = '\0';

          for (size_t i = 0; i < settings.length(); ++i) {
            char c = settings[i];

            // Track if we're inside a string
            if (c == '"' && prev_char != '\\') {
              in_string = !in_string;
            }

            if (!in_string) {
              // Add newline and indent after opening braces/brackets
              if (c == '{' || c == '[') {
                formatted += c;
                formatted += '\n';
                indent += 2;
                formatted += std::string(indent, ' ');
                prev_char = c;
                continue;
              }
              // Add newline and unindent before closing braces/brackets
              else if (c == '}' || c == ']') {
                formatted += '\n';
                indent -= 2;
                formatted += std::string(indent, ' ');
                formatted += c;
                prev_char = c;
                continue;
              }
              // Add newline and indent after commas
              else if (c == ',') {
                formatted += c;
                formatted += '\n';
                formatted += std::string(indent, ' ');
                prev_char = c;
                continue;
              }
              // Skip extra whitespace
              else if (c == ' ' && (prev_char == ',' || prev_char == '{' ||
                                    prev_char == '[')) {
                continue;
              }
            }

            formatted += c;
            prev_char = c;
          }

          std::cout << formatted << "\n";
          break;
        }
      }
      if (!found) {
        std::cout << "  Subsystem '" << blob_subsystem
                  << "' not found in configs\n";
      }
    }
  } else {
    // Normal print
    for (auto const& rr : runvec) {
      tool.printRun(rr);
    }
  }

  return 0;
}
