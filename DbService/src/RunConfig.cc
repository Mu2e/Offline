#include "Offline/DbService/inc/RunConfig.hh"

#include "cetlib_except/exception.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "nlohmann/json.hpp"

namespace {

using nlohmann::json;

// Convert a JSON value to a plain string: use the raw content for strings,
// and a compact dump() for anything else (number, bool, etc.).
std::string valueToString(const json& v) {
  return v.is_string() ? v.get<std::string>() : v.dump();
}

// Recursively walk a JSON node.  Whenever an object contains the key
// "DBServiceTables" whose value is an object, its key-value pairs are appended
// to `collected`.  Nested objects and arrays are always traversed so that
// "DBServiceTables" entries at any depth are found.
void walk(const json& node,
          std::vector<std::pair<std::string, std::string>>& collected) {
  if (node.is_object()) {
    for (auto it = node.begin(); it != node.end(); ++it) {
      if (it.key() == "DBServiceTables" && it.value().is_object()) {
        for (auto jt = it.value().begin(); jt != it.value().end(); ++jt) {
          collected.emplace_back(jt.key(), valueToString(jt.value()));
        }
      } else {
        walk(it.value(), collected);
      }
    }
  } else if (node.is_array()) {
    for (auto const& elem : node) {
      walk(elem, collected);
    }
  }
}

// Recursively walk a JSON node.  Whenever an object contains the key
// "DBServiceCIDTable" whose value is an array, each element (expected to be
// an object with "cid" and "name" fields) is appended to `collected` as a
// (cid, name) pair.  Nested objects and arrays are always traversed so that
// "DBServiceCIDTable" entries at any depth are found.
void walkCID(const json& node,
             std::vector<std::pair<std::string, std::string>>& collected) {
  if (node.is_object()) {
    for (auto it = node.begin(); it != node.end(); ++it) {
      if (it.key() == "DBServiceCIDTable" && it.value().is_array()) {
        for (auto const& elem : it.value()) {
          if (elem.is_object() && elem.contains("cid") &&
              elem.contains("name")) {
            collected.emplace_back(valueToString(elem["cid"]),
                                   valueToString(elem["name"]));
          }
        }
      } else {
        walkCID(it.value(), collected);
      }
    }
  } else if (node.is_array()) {
    for (auto const& elem : node) {
      walkCID(elem, collected);
    }
  }
}

// ---------------------------------------------------------------------------
// Trigger menu (TriggerConfigTable) support for RunConfig::triggerMenu
// ---------------------------------------------------------------------------

// One event-mode entry of one trigger path.
struct TrigMode {
  std::string mode;
  int prescale = 0;
  std::vector<std::string> streams;
};

// One trigger path: its bit, whether it is enabled, and its per-event-mode
// configuration.
struct TrigPath {
  std::string name;
  int bit = 0;
  bool enabled = false;
  std::vector<TrigMode> modes;
};

// Tolerant integer accessor.  The blob is not type-consistent: the
// TriggerConfigTable's own "version" is a JSON *string* ("19") while the
// menu's "version" one level down is a number, and "enabled" is an int 0/1
// rather than a bool.  A plain get<int>() throws nlohmann::json::type_error on
// the string form, so convert explicitly instead of trusting the type.
int toInt(const json& v, int def = 0) {
  if (v.is_number_integer()) return v.get<int>();
  if (v.is_number()) return static_cast<int>(v.get<double>());
  if (v.is_boolean()) return v.get<bool>() ? 1 : 0;
  if (v.is_string()) {
    try {
      return std::stoi(v.get<std::string>());
    } catch (std::exception const&) {
      return def;  // non-numeric string: treat as absent rather than throw
    }
  }
  return def;
}

// Recursively walk a JSON node for a "TriggerConfigTable" object carrying a
// "content" field, returning a pointer to it (nullptr if there is none).
// Recursive rather than a fixed lookup because the node's depth depends on
// where the blob came from -- it sits at custom.TriggerConfigTable in the
// config view's output and at config.custom.TriggerConfigTable in the raw
// config table's -- the same two shapes the Python readers of this column
// have to handle.  Follows walk()/walkCID() above in traversing objects and
// arrays alike so depth is never assumed.
const json* walkTrigger(const json& node) {
  if (node.is_object()) {
    for (auto it = node.begin(); it != node.end(); ++it) {
      if (it.key() == "TriggerConfigTable" && it.value().is_object() &&
          it.value().contains("content")) {
        return &it.value();
      }
      if (const json* found = walkTrigger(it.value())) return found;
    }
  } else if (node.is_array()) {
    for (auto const& elem : node) {
      if (const json* found = walkTrigger(elem)) return found;
    }
  }
  return nullptr;
}

// The menu payload under "content".  Normally an object, but it can arrive as
// a JSON string holding that object depending on how the blob was written, so
// a string is parsed a second time (non-throwing, as everywhere else here).
// Returns a discarded value when it is neither form.
json trigContent(const json& table) {
  json const& content = table.at("content");
  if (content.is_object()) return content;
  if (content.is_string()) {
    return json::parse(content.get<std::string>(), nullptr, false);
  }
  return json(json::value_t::discarded);
}

// Parse the "trigger_paths" object into TrigPath records, sorted by bit (then
// by name, so the order is total and stable even if two paths share a bit).
std::vector<TrigPath> trigPaths(const json& paths) {
  std::vector<TrigPath> out;
  if (!paths.is_object()) return out;

  for (auto it = paths.begin(); it != paths.end(); ++it) {
    json const& entry = it.value();
    if (!entry.is_object()) continue;

    TrigPath path;
    path.name = it.key();
    if (entry.contains("bit")) path.bit = toInt(entry.at("bit"));
    if (entry.contains("enabled")) {
      path.enabled = toInt(entry.at("enabled")) != 0;
    }
    if (entry.contains("eventModeConfig") &&
        entry.at("eventModeConfig").is_array()) {
      for (auto const& cfg : entry.at("eventModeConfig")) {
        if (!cfg.is_object()) continue;
        TrigMode mode;
        if (cfg.contains("eventMode") && cfg.at("eventMode").is_string()) {
          mode.mode = cfg.at("eventMode").get<std::string>();
        }
        if (cfg.contains("prescale")) mode.prescale = toInt(cfg.at("prescale"));
        if (cfg.contains("streams") && cfg.at("streams").is_array()) {
          for (auto const& stream : cfg.at("streams")) {
            if (stream.is_string()) {
              mode.streams.push_back(stream.get<std::string>());
            }
          }
        }
        path.modes.push_back(mode);
      }
    }
    out.push_back(path);
  }

  std::sort(out.begin(), out.end(),
            [](TrigPath const& a, TrigPath const& b) {
              if (a.bit != b.bit) return a.bit < b.bit;
              return a.name < b.name;
            });
  return out;
}

// The event modes present in the menu, in order of first appearance.
// Discovered from the data rather than hardcoded to OnSpill/OffSpill: the
// number of modes is a property of the menu, and a path is not required to
// configure every one of them.
std::vector<std::string> trigModes(std::vector<TrigPath> const& paths) {
  std::vector<std::string> modes;
  for (auto const& path : paths) {
    for (auto const& mode : path.modes) {
      if (std::find(modes.begin(), modes.end(), mode.mode) == modes.end()) {
        modes.push_back(mode.mode);
      }
    }
  }
  return modes;
}

// Parse a "*_streams" object -- {name: {"enabled": 0|1}} -- into (name,
// enabled) pairs.  A bare value instead of an object is tolerated as the
// enabled flag itself.
std::vector<std::pair<std::string, bool>> trigStreams(const json& node) {
  std::vector<std::pair<std::string, bool>> out;
  if (!node.is_object()) return out;
  for (auto it = node.begin(); it != node.end(); ++it) {
    bool enabled = false;
    if (it.value().is_object()) {
      if (it.value().contains("enabled")) {
        enabled = toInt(it.value().at("enabled")) != 0;
      }
    } else {
      enabled = toInt(it.value()) != 0;
    }
    out.emplace_back(it.key(), enabled);
  }
  return out;
}

// A prescale as the text format shows it.  -1 means the path is not taken in
// that event mode; the JSON form keeps the raw -1 instead (see triggerMenu).
std::string psText(int prescale) {
  return prescale == -1 ? std::string("off") : std::to_string(prescale);
}

// The fixed JSON schema, with every field at its empty default.  Both the
// found and not-found paths start here so the two can never disagree about
// which keys exist -- a consumer can rely on the shape without checking
// "found" first.
json trigSkeleton(int run_number, bool found) {
  json out = json::object();
  out["run"] = run_number;
  out["found"] = found;
  out["table_name"] = "";
  out["table_alias"] = "";
  out["table_version"] = "";
  out["menu_version"] = 0;
  out["event_modes"] = json::array();
  out["trigger_paths"] = json::array();
  out["agg_trigger_paths"] = json::array();
  out["dataLogger_streams"] = json::array();
  out["lumiLogger_streams"] = json::array();
  out["dataLogger_analyzers"] = json::object();
  out["lumiLogger_analyzers"] = json::object();
  out["crv_agg_process_name"] = "";
  out["trkcal_filter_process_name"] = "";
  return out;
}

// TrigPath records as a JSON array.  An array rather than an object keyed by
// path name, for the same reason flagsJson() uses one: this is a list of
// records, and the array preserves the bit ordering an object cannot.
json trigPathsJson(std::vector<TrigPath> const& paths) {
  json out = json::array();
  for (auto const& path : paths) {
    json entry = json::object();
    entry["path"] = path.name;
    entry["bit"] = path.bit;
    entry["enabled"] = path.enabled;
    json modes = json::array();
    for (auto const& mode : path.modes) {
      json m = json::object();
      m["mode"] = mode.mode;
      // Raw -1 rather than the text format's "off": unambiguous for a machine
      // reader, same philosophy as runJson() emitting Unix timestamps where
      // printRun() renders localized strings.
      m["prescale"] = mode.prescale;
      m["streams"] = mode.streams;
      modes.push_back(m);
    }
    entry["event_modes"] = modes;
    out.push_back(entry);
  }
  return out;
}

// (name, enabled) stream pairs as a JSON array of objects.
json trigStreamsJson(
    std::vector<std::pair<std::string, bool>> const& streams) {
  json out = json::array();
  for (auto const& stream : streams) {
    json entry = json::object();
    entry["stream"] = stream.first;
    entry["enabled"] = stream.second;
    out.push_back(entry);
  }
  return out;
}

// Append "label: a b c" to os, wrapping at kWrap columns with a hanging
// indent.  Streams are listed inline rather than one per line because in
// practice they are few and nearly always all enabled -- a table per stream
// would spend ten lines saying nothing.
void writeStreamLine(std::ostringstream& os, std::string const& label,
                     std::vector<std::pair<std::string, bool>> const& streams) {
  constexpr size_t kWrap = 78;
  if (streams.empty()) return;

  std::string const head =
      "  " + label + " (" + std::to_string(streams.size()) + "):";
  // Continuation lines hang under the first stream name, past the label.
  std::string const indent(head.size() + 1, ' ');

  std::string line = head;
  for (auto const& stream : streams) {
    // Only disabled streams are marked, so the common all-enabled case reads
    // as a plain list.
    std::string item = stream.first + (stream.second ? "" : "(off)");
    if (line.size() + 1 + item.size() > kWrap && line != head) {
      os << line << "\n";
      line = indent;
    }
    line += " " + item;
  }
  os << line << "\n";
}

}  // namespace

// ---------------------------------------------------------------------------
// mu2e::RunConfig::dbTables2
// ---------------------------------------------------------------------------
std::string mu2e::RunConfig::dbTables2(bool qjson) const {
  std::vector<std::pair<std::string, std::string>> pairs;

  // Parse the settings blob; tolerate malformed/empty content by returning
  // an empty result rather than throwing.
  nlohmann::json root = nlohmann::json::parse(_settings, nullptr, false);
  if (!root.is_discarded()) {
    walkCID(root, pairs);
  }

  std::string result;

  if (qjson) {
    // Build a JSON array of {"cid":..., "name":...} objects and serialize it,
    // reproducing the DBServiceCIDTable list.
    nlohmann::json out = nlohmann::json::array();
    for (auto const& p : pairs) {
      nlohmann::json entry = nlohmann::json::object();
      entry["cid"] = p.first;
      entry["name"] = p.second;
      out.push_back(entry);
    }
    result = out.dump(1, ' ');
  } else {
    // One (cid, name) pair per line: the CID followed by the name.
    for (auto const& p : pairs) {
      result += p.first;
      result += " ";
      result += p.second;
      result += "\n";
    }
    if (!result.empty() && result.back() == '\n') result.pop_back();
  }

  return result;
}

// ---------------------------------------------------------------------------
// mu2e::RunConfig::dbTables3
// ---------------------------------------------------------------------------

std::string mu2e::RunConfig::dbTables3(bool qjson) const {
  std::vector<std::pair<std::string, std::string>> pairs;

  // Parse the settings blob; tolerate malformed/empty content by returning
  // an empty result rather than throwing.
  nlohmann::json root = nlohmann::json::parse(_settings, nullptr, false);
  if (!root.is_discarded()) {
    walk(root, pairs);
  }

  std::string result;

  if (qjson) {
    // Build a JSON object and let nlohmann serialize it.  dump() correctly
    // escapes any control characters (e.g. newlines become \n in the text),
    // so the result is always valid JSON that json.loads can parse back into
    // values that retain their embedded newlines.
    nlohmann::json out = nlohmann::json::object();
    for (auto const& p : pairs) {
      // A duplicate key would silently overwrite the earlier value, losing
      // data, so refuse to continue if one is found.
      if (out.contains(p.first)) {
        throw cet::exception("RUNCONFIG_DUPLICATE_KEY")
            << " RunConfig::dbTables3 found duplicate DBServiceTables key "
            << p.first << " which would overwrite an existing value \n";
      }
      out[p.first] = p.second;
    }
    result = out.dump(1, ' ');
  } else {
    for (auto const& p : pairs) {
      result += p.second;
      result += "\n";
    }
    if (!result.empty() && result.back() == '\n') result.pop_back();
  }

  return result;
}

// ---------------------------------------------------------------------------
// mu2e::RunConfig::triggerMenuNotFound
// ---------------------------------------------------------------------------
std::string mu2e::RunConfig::triggerMenuNotFound(int run_number) {
  return trigSkeleton(run_number, false).dump(1, ' ');
}

// ---------------------------------------------------------------------------
// mu2e::RunConfig::triggerMenu
// ---------------------------------------------------------------------------
std::string mu2e::RunConfig::triggerMenu(bool qjson) const {
  // Parse the settings blob; tolerate malformed/empty content by returning
  // an empty result rather than throwing, as dbTables2/3 above do.
  nlohmann::json root = nlohmann::json::parse(_settings, nullptr, false);
  if (root.is_discarded()) return std::string();

  const nlohmann::json* table = walkTrigger(root);
  if (table == nullptr) return std::string();

  nlohmann::json content = trigContent(*table);
  if (!content.is_object()) return std::string();

  std::vector<TrigPath> paths =
      trigPaths(content.value("trigger_paths", nlohmann::json::object()));
  std::vector<std::string> modes = trigModes(paths);
  std::vector<std::pair<std::string, bool>> dataLoggers =
      trigStreams(content.value("dataLogger_streams", nlohmann::json::object()));
  std::vector<std::pair<std::string, bool>> lumiLoggers =
      trigStreams(content.value("lumiLogger_streams", nlohmann::json::object()));

  std::string tableName =
      table->contains("name") ? valueToString(table->at("name")) : std::string();
  std::string tableAlias = table->contains("alias")
                               ? valueToString(table->at("alias"))
                               : std::string();
  // Kept as text, not an int: this field is a JSON string ("19") in the blob
  // while the menu version below is a number, and coercing it would invent a
  // type the source does not have.
  std::string tableVersion = table->contains("version")
                                 ? valueToString(table->at("version"))
                                 : std::string();
  int menuVersion =
      content.contains("version") ? toInt(content.at("version")) : 0;
  std::string crvAgg = content.contains("crv_agg_process_name")
                           ? valueToString(content.at("crv_agg_process_name"))
                           : std::string();
  std::string trkcalFilter =
      content.contains("trkcal_filter_process_name")
          ? valueToString(content.at("trkcal_filter_process_name"))
          : std::string();

  if (qjson) {
    nlohmann::json out = trigSkeleton(_run_number, true);
    out["table_name"] = tableName;
    out["table_alias"] = tableAlias;
    out["table_version"] = tableVersion;
    out["menu_version"] = menuVersion;
    out["event_modes"] = modes;
    out["trigger_paths"] = trigPathsJson(paths);
    // agg_trigger_paths has been empty in every run inspected, so its element
    // shape is unconfirmed.  If it is an object it is assumed to match
    // trigger_paths and is parsed the same way; an array is passed through
    // untouched rather than reshaped into a form that might be wrong.
    if (content.contains("agg_trigger_paths")) {
      nlohmann::json const& agg = content.at("agg_trigger_paths");
      if (agg.is_object()) {
        out["agg_trigger_paths"] = trigPathsJson(trigPaths(agg));
      } else if (agg.is_array()) {
        out["agg_trigger_paths"] = agg;
      }
    }
    out["dataLogger_streams"] = trigStreamsJson(dataLoggers);
    out["lumiLogger_streams"] = trigStreamsJson(lumiLoggers);
    // Passed through as-is: these are empty objects in the runs inspected and
    // their contents are not interpreted here.
    if (content.contains("dataLogger_analyzers")) {
      out["dataLogger_analyzers"] = content.at("dataLogger_analyzers");
    }
    if (content.contains("lumiLogger_analyzers")) {
      out["lumiLogger_analyzers"] = content.at("lumiLogger_analyzers");
    }
    out["crv_agg_process_name"] = crvAgg;
    out["trkcal_filter_process_name"] = trkcalFilter;
    return out.dump(1, ' ');
  }

  // --- text format -------------------------------------------------------
  std::ostringstream os;

  os << "  table " << (tableName.empty() ? "?" : tableName);
  if (!tableAlias.empty()) os << "  alias " << tableAlias;
  if (!tableVersion.empty()) os << "  version " << tableVersion;
  os << "  menu version " << menuVersion << "\n";

  if (paths.empty()) {
    os << "\n  no trigger paths\n";
  } else {
    // Column widths from the content, so a long path name or a wide prescale
    // can't break the alignment.
    size_t pathw = std::string("path").size();
    size_t bitw = std::string("bit").size();
    for (auto const& path : paths) {
      pathw = std::max(pathw, path.name.size());
      bitw = std::max(bitw, std::to_string(path.bit).size());
    }
    size_t const onw = std::string("OFF").size();

    std::map<std::string, size_t> modew;
    for (auto const& mode : modes) {
      modew[mode] = std::max(mode.size(), psText(-1).size());
    }
    for (auto const& path : paths) {
      for (auto const& mode : path.modes) {
        modew[mode.mode] =
            std::max(modew[mode.mode], psText(mode.prescale).size());
      }
    }

    size_t total = 2 + pathw + 2 + bitw + 2 + onw;
    for (auto const& mode : modes) total += 2 + modew[mode];

    os << "\n";
    os << "  " << std::left << std::setw(pathw) << "path" << std::right << "  "
       << std::setw(bitw) << "bit" << "  " << std::setw(onw) << "on";
    for (auto const& mode : modes) {
      os << "  " << std::setw(modew[mode]) << mode;
    }
    os << "\n  " << std::string(total - 2, '-') << "\n";

    for (auto const& path : paths) {
      os << "  " << std::left << std::setw(pathw) << path.name << std::right
         << "  " << std::setw(bitw) << path.bit << "  " << std::setw(onw)
         << (path.enabled ? "ON" : "OFF");
      for (auto const& mode : modes) {
        auto it = std::find_if(path.modes.begin(), path.modes.end(),
                               [&mode](TrigMode const& m) {
                                 return m.mode == mode;
                               });
        // "-" where the path configures no entry for this event mode, which
        // is distinct from "off" (configured, prescale -1).
        os << "  " << std::setw(modew[mode])
           << (it == path.modes.end() ? std::string("-") : psText(it->prescale));
      }
      os << "\n";
    }
  }

  if (!dataLoggers.empty() || !lumiLoggers.empty() || !crvAgg.empty() ||
      !trkcalFilter.empty()) {
    os << "\n";
    writeStreamLine(os, "dataLogger streams", dataLoggers);
    writeStreamLine(os, "lumiLogger streams", lumiLoggers);
    if (!crvAgg.empty()) os << "  crv_agg_process_name: " << crvAgg << "\n";
    if (!trkcalFilter.empty()) {
      os << "  trkcal_filter_process_name: " << trkcalFilter << "\n";
    }
  }

  std::string result = os.str();
  if (!result.empty() && result.back() == '\n') result.pop_back();
  return result;
}
