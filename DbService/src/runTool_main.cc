//
// main to drive a command line interface to RunTool,
// which reads data from the run conditions database
//

#include "Offline/DbService/inc/RunSelect.hh"
#include "Offline/DbService/inc/RunTool.hh"
#include "Offline/GeneralUtilities/inc/ParseCLI.hh"

#include "nlohmann/json.hpp"

#include <iomanip>
#include <iostream>

using namespace mu2e;

int main(int argc, char** argv) {
  // Single handler for the whole body: RunTool/RunInfo/RunConfig throw
  // cet::exception on a range of conditions (DB unreachable, a run's configs
  // disagreeing on a table key, malformed flag-table content, ...), and
  // std::stoi below throws std::invalid_argument on non-numeric --last/--days
  // input (e.g. a typo'd CLI arg). All derive from std::exception (verified
  // directly, including nlohmann::json::parse_error), so one catch here
  // covers every throw site in this file instead of the narrower per-branch
  // handling this used to have -- left uncaught, any of them aborts the
  // process with a core dump instead of a clean, parseable error.
  try {
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
    cli.addSwitch("", "cidtables", "e", "cidtables", false,
                  "print the cat 2 DbService CID tables (cid name)");
    cli.addSwitch("", "json", "j", "json", false,
                  "output content in JSON format");


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
      // JSON shaping lives on RunTool::flagsJson(), not here, so it's callable
      // from the Python/SWIG binding too, not just this CLI -- see its comment.
      if (cli.getBool("", "json")) {
        std::cout << tool.flagsJson() << "\n";
        return 0;
      }

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
    bool cidtables = cli.getBool("", "cidtables");
    bool json = cli.getBool("", "json");

    // If blob subsystem or dbtables is specified, we need to fetch configs
    bool need_configs =
        configs || !blob_subsystem.empty() || dbtables || cidtables;


    // Query and print runs
    RunInfo::RunVec runvec =
        tool.listRuns(runsel, need_configs, transitions, subruns);

    // Handle dbtables output: cat-3 DbService table values, collected across all
    // configs of each run by RunInfo::dbTables3.  With -j the merged content is
    // printed as formatted JSON; otherwise the values are printed one per line
    // for machine reading (no labels or adornment).
    if (dbtables) {
      // dbTables3() output carries no per-run label, and in JSON mode is one
      // dump() per run with no enclosing array -- concatenating more than
      // one is not valid JSON (or, in text mode, an unlabeled jumble no
      // reader could attribute back to a run). This is only meaningful for
      // a single run, and that's enforced here rather than just assumed --
      // a selection matching more than one run is rejected outright instead
      // of silently emitting output nothing downstream can parse.
      if (runvec.size() > 1) {
        std::cerr << "ERROR - -q/--dbtables requires a selection matching "
                     "exactly one run (matched "
                  << runvec.size()
                  << "); narrow with --run, --last 1, or a tighter "
                     "--time/--days window.\n";
        return 1;
      }

      // RunInfo::dbTables3()/RunConfig::dbTables3() throw
      // (RUNINFO_DUPLICATE_TABLE / RUNCONFIG_DUPLICATE_KEY) when a run's
      // configs disagree on a table key -- a rare user misconfiguration, but
      // still a fatal condition: merging past it would silently drop data.
      // Caught by the single handler around this whole function, same as
      // every other error path here.
      for (auto const& rr : runvec) {
        std::string tables = rr.dbTables3(json);
        if (!tables.empty()) {
          std::cout << tables << "\n";
        }
      }
      return 0;
    }

    // Handle cidtables output: cat-2 DbService CID table entries, collected by
    // RunInfo::dbTables2 (which restricts to the "TRG" config that holds them).
    // With -j the content is fetched in JSON format and printed as formatted
    // JSON; otherwise it is printed as "cid name" one pair per line.
    if (cidtables) {
      // Same single-run constraint as dbtables above and for the same
      // reason -- dbTables2() output carries no per-run label. Worth noting
      // this isn't only a JSON-mode concern: "[]" (no TRG config, or none
      // with CID entries) still counts as non-empty per the check below, so
      // even the JSON array case can emit multiple concatenated documents.
      if (runvec.size() > 1) {
        std::cerr << "ERROR - -e/--cidtables requires a selection matching "
                     "exactly one run (matched "
                  << runvec.size()
                  << "); narrow with --run, --last 1, or a tighter "
                     "--time/--days window.\n";
        return 1;
      }

      for (auto const& rr : runvec) {
        std::string tables = rr.dbTables2(json);
        if (!tables.empty()) {
          std::cout << tables << "\n";
        }
      }
      return 0;
    }


    // Handle JSON blob output for specific subsystem
    if (!blob_subsystem.empty()) {
      // Two separate output paths, not a shared one -- see below for why.
      nlohmann::json blobArr = nlohmann::json::array();  // used only when json

      for (auto const& rr : runvec) {
        if (!json) {
          std::cout << "Run " << rr.runNumber() << ":\n";
        }
        bool found = false;
        for (const auto& config : rr.configs()) {
          if (config.subsystem() == blob_subsystem) {
            found = true;
            // Get the settings string (already cleaned in RunTool.cc)
            std::string settings = config.settings();

            if (json) {
              // settings is already valid JSON text (a jsonb column) -- parse
              // and re-dump it through a real JSON library rather than the
              // human-readable formatter below. That formatter intentionally
              // unescapes \n into real newlines for readability, which is a
              // deliberate choice for the human-facing default, but it means
              // that output is NOT valid JSON (an unescaped control character
              // inside a JSON string isn't legal) -- hence a genuinely
              // separate path here rather than reusing/tweaking it.
              nlohmann::json entry = nlohmann::json::object();
              entry["run"] = rr.runNumber();
              entry["subsystem"] = blob_subsystem;
              entry["found"] = true;
              // Non-throwing overload, matching how RunConfig::dbTables2/3()
              // and RunInfo already parse this exact same settings column --
              // empty/malformed content degrades to a discarded (null) value
              // instead of an uncaught parse_error aborting the process (a
              // real, reachable case: DbUtil::simplfyQeString() passes empty/
              // single-character content through unvalidated).
              nlohmann::json parsed =
                  nlohmann::json::parse(settings, nullptr, false);
              entry["settings"] = parsed.is_discarded() ? nlohmann::json(nullptr) : parsed;
              blobArr.push_back(entry);
              break;
            }

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
          if (json) {
            nlohmann::json entry = nlohmann::json::object();
            entry["run"] = rr.runNumber();
            entry["subsystem"] = blob_subsystem;
            entry["found"] = false;
            entry["settings"] = nullptr;
            blobArr.push_back(entry);
          } else {
            std::cout << "  Subsystem '" << blob_subsystem
                      << "' not found in configs\n";
          }
        }
      }

      if (json) {
        std::cout << blobArr.dump(1, ' ') << "\n";
      }
    } else if (json) {
      // Normal run listing, JSON: one array of per-run objects.
      // RunTool::runJson() returns std::string (not nlohmann::json -- see its
      // own comment for why), so each run's JSON text is parsed back in here
      // to assemble a single valid array, rather than concatenating one
      // document per run the way -q/-e do (fine there since dbtables isn't
      // meant to be queried as a list; not fine here, since a run listing
      // routinely is).
      nlohmann::json arr = nlohmann::json::array();
      for (auto const& rr : runvec) {
        arr.push_back(nlohmann::json::parse(tool.runJson(rr)));
      }
      std::cout << arr.dump(1, ' ') << "\n";
    } else {
      // Normal print
      for (auto const& rr : runvec) {
        tool.printRun(rr);
      }
    }

    return 0;
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
