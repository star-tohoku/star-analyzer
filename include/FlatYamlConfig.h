#ifndef FLAT_YAML_CONFIG_H
#define FLAT_YAML_CONFIG_H

#include "YamlParser.h"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Strict reader for flat "key: value" YAML files (the YamlParser format used by FemtoConfig and
// the other cut configs). Required keys are checked, not defaulted; a malformed value is an error;
// Unused() lists keys that were never read (typos). Lists are comma separated: "a, b, c".
// Used by StCandidateFemtoMaker and its parts, whose YAML is not registered in ConfigManager.
class FlatYamlConfig {
public:
  FlatYamlConfig() : mOk(true) {}

  bool Load(const std::string& path) {
    mPath = path;
    mValues.clear();
    mUsed.clear();
    mOk = YamlParser::ParseFile(path.c_str(), mValues);
    if (!mOk) Fail("cannot parse");
    return mOk;
  }

  // Resolve a file referenced by a mainconf key, like ConfigManager: <project>/config/<value>,
  // where <project> is the directory that contains config/mainconf/. Prints the resolved path.
  static bool ResolveMainConfigReference(const std::string& mainConfigPath, const std::string& key,
                                         std::string& path) {
    path.clear();
    std::map<std::string, std::string> values;
    if (mainConfigPath.empty() || !YamlParser::ParseFile(mainConfigPath.c_str(), values)) {
      std::cerr << "ERROR: cannot parse mainconf '" << mainConfigPath << "'" << std::endl;
      return false;
    }
    std::map<std::string, std::string>::const_iterator it = values.find(key);
    if (it == values.end() || it->second.empty()) {
      std::cerr << "ERROR: mainconf '" << mainConfigPath << "' has no key '" << key << "'" << std::endl;
      return false;
    }
    const std::string marker = "config/mainconf/";
    const std::string::size_type pos = mainConfigPath.rfind(marker);
    if (pos == std::string::npos) {
      std::cerr << "ERROR: mainconf path must be under config/mainconf/: " << mainConfigPath << std::endl;
      return false;
    }
    path = mainConfigPath.substr(0, pos) + "config/" + it->second;
    std::ifstream probe(path.c_str());
    if (!probe) {
      std::cerr << "ERROR: mainconf key '" << key << "' points to a missing file: " << path << std::endl;
      path.clear();
      return false;
    }
    std::cout << "[config] " << key << " = " << path << " (from " << mainConfigPath << ")" << std::endl;
    return true;
  }

  bool Has(const std::string& key) const { return mValues.find(key) != mValues.end(); }

  std::string RequireString(const std::string& key) {
    std::string v;
    if (!Find(key, v) || v.empty()) Fail("required key '" + key + "' missing");
    return v;
  }
  double RequireDouble(const std::string& key) {
    std::string v;
    if (!Find(key, v)) { Fail("required key '" + key + "' missing"); return 0.; }
    return ToDouble(key, v);
  }
  int RequireInt(const std::string& key) {
    std::string v;
    if (!Find(key, v)) { Fail("required key '" + key + "' missing"); return 0; }
    return ToInt(key, v);
  }
  bool RequireBool(const std::string& key) {
    std::string v;
    if (!Find(key, v)) { Fail("required key '" + key + "' missing"); return false; }
    return ToBool(key, v);
  }
  // Optional keys: the caller documents the value used when the key is absent.
  double OptionalDouble(const std::string& key, double absent) {
    std::string v;
    return Find(key, v) ? ToDouble(key, v) : absent;
  }
  bool OptionalBool(const std::string& key, bool absent) {
    std::string v;
    return Find(key, v) ? ToBool(key, v) : absent;
  }
  std::string OptionalString(const std::string& key, const std::string& absent) {
    std::string v;
    return Find(key, v) ? v : absent;
  }
  std::vector<std::string> RequireList(const std::string& key) {
    std::vector<std::string> out;
    std::string v;
    if (!Find(key, v)) { Fail("required key '" + key + "' missing"); return out; }
    out = SplitList(v);
    if (out.empty()) Fail("key '" + key + "' is an empty list");
    return out;
  }
  std::vector<std::string> OptionalList(const std::string& key) {
    std::string v;
    return Find(key, v) ? SplitList(v) : std::vector<std::string>();
  }

  // Keys present in the file but never read.
  std::vector<std::string> Unused() const {
    std::vector<std::string> out;
    for (std::map<std::string, std::string>::const_iterator it = mValues.begin(); it != mValues.end(); ++it)
      if (!mUsed.count(it->first)) out.push_back(it->first);
    return out;
  }
  // Fail on keys that were never read (call after reading everything).
  void RejectUnused() {
    const std::vector<std::string> unused = Unused();
    for (size_t i = 0; i < unused.size(); ++i) Fail("unknown key '" + unused[i] + "'");
  }

  void Fail(const std::string& message) {
    mOk = false;
    std::cerr << "ERROR: " << mPath << ": " << message << std::endl;
  }
  bool Ok() const { return mOk; }
  const std::string& Path() const { return mPath; }

  static std::vector<std::string> SplitList(const std::string& value) {
    std::vector<std::string> out;
    std::string item;
    std::istringstream in(value);
    while (std::getline(in, item, ',')) {
      const std::string t = Trim(item);
      if (!t.empty()) out.push_back(t);
    }
    return out;
  }
  static std::string Trim(const std::string& s) {
    const std::string::size_type b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const std::string::size_type e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
  }
  // Strict number parsing for list items etc.; ok = false if the whole string is not a number.
  static double ParseDouble(const std::string& s, bool& ok) {
    const std::string t = Trim(s);
    char* end = 0;
    errno = 0;
    const double v = std::strtod(t.c_str(), &end);
    ok = !t.empty() && end && *end == '\0' && errno == 0;
    return v;
  }

private:
  bool Find(const std::string& key, std::string& value) {
    std::map<std::string, std::string>::const_iterator it = mValues.find(key);
    if (it == mValues.end()) return false;
    mUsed.insert(key);
    value = it->second;
    return true;
  }
  double ToDouble(const std::string& key, const std::string& v) {
    bool ok = false;
    const double d = ParseDouble(v, ok);
    if (!ok) Fail("key '" + key + "' is not a number: '" + v + "'");
    return d;
  }
  int ToInt(const std::string& key, const std::string& v) {
    const std::string t = Trim(v);
    char* end = 0;
    errno = 0;
    const long l = std::strtol(t.c_str(), &end, 10);
    if (t.empty() || !end || *end != '\0' || errno != 0) Fail("key '" + key + "' is not an integer: '" + v + "'");
    return static_cast<int>(l);
  }
  bool ToBool(const std::string& key, const std::string& v) {
    const std::string t = Trim(v);
    if (t == "true" || t == "True" || t == "1") return true;
    if (t == "false" || t == "False" || t == "0") return false;
    Fail("key '" + key + "' is not a boolean: '" + v + "'");
    return false;
  }

  std::string mPath;
  std::map<std::string, std::string> mValues;
  std::set<std::string> mUsed;
  bool mOk;
};

#endif
