#include "Offline/DbService/inc/RunConfig.hh"
#include "cetlib_except/exception.h"
#include "nlohmann/json.hpp"
#include <string>
#include <utility>
#include <vector>
#include <iostream>


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

}

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
      std::cout << "DEB " << p.first << "\n";
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
