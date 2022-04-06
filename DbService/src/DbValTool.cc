
#include "Offline/DbService/inc/DbValTool.hh"
#include "cetlib_except/exception.h"

namespace mu2e {

//**************************************************

int DbValTool::findPid(std::string purpose) const {
  if (purpose.empty()) return -1;

  for (auto const& pp : _valcache.valPurposes().rows()) {
    if (pp.name() == purpose) return pp.pid();
  }
  return -1;
}

//**************************************************

void DbValTool::findPidVid(std::string purpose, std::string version, int& pid,
                           int& vid) const {
  pid = -1;
  vid = -1;

  if (purpose.empty() && !version.empty()) {
    throw cet::exception("DBVALTOOL_VERISON_WITHOUT_PURPOSE")
        << " DbValTool::findPidVid version but no purpose \n";
  }

  // both are empty
  if (purpose.empty()) return;

  pid = findPid(purpose);

  if (pid < 0) {
    throw cet::exception("DBVALTOOL_PURPOSE_LOOKUP_FAILED")
        << " DbValTool::findPidVid failed to find purpose " << purpose << "\n";
  }

  // purpose found, version empty
  if (version.empty()) return;

  DbVersion dbver(purpose, version);
  if (dbver.major() < 0 || dbver.minor() < 0) {
    throw cet::exception("DBVALTOOL_VERISON_MISSING_CONTENT")
        << " DbValTool::findPidVid version " << version
        << " does not have major and minor numbers\n";
  }

  for (auto const& vv : _valcache.valVersions().rows()) {
    if (vv.pid() == pid && vv.major() == dbver.major() &&
        vv.minor() == dbver.minor()) {
      vid = vv.vid();
    }
  }

  if (vid < 0) {
    throw cet::exception("DBVALTOOL_VERSION_LOOKUP_FAILED")
        << " DbValTool::findPidVid failed to find verison " << version
        << " for purpose " << purpose << "\n";
  }

  return;
}

//**************************************************

bool DbValTool::tidByName(std::string const& name, int& tid) const {
  for (auto const& r : _valcache.valTables().rows()) {
    if (r.name() == name) {
      tid = r.tid();
      return true;
    }
  }
  tid = -1;
  return false;
}

//**************************************************

bool DbValTool::nameByTid(int tid, std::string& name) const {
  for (auto const& r : _valcache.valTables().rows()) {
    if (r.tid() == tid) {
      name = r.name();
      return true;
    }
  }
  name.clear();
  return false;
}

//**************************************************

void DbValTool::printSet(DbSet const& dbset) const {
  auto const& emap = dbset.emap();
  std::cout << "            Table        N IoV" << std::endl;
  std::string name;
  for (auto const& p : emap) {
    nameByTid(p.first, name);
    std::cout << std::setw(20) << name << std::setw(8) << p.second.size()
              << std::endl;
  }
}

//**************************************************

void DbValTool::fillSetVer(DbVersion const& dver, DbSet& dbset) const {
  int pid = findPid(dver.purpose());

  if (pid < 0) {
    throw cet::exception("DBVALTOOL_BAD_PURPOSE")
        << " DbValTool::fillSet calibration purpose string not found in the "
           "DB: "
        << dver.purpose() << "\n";
  }

  // confirm version numbers and find version number (vid)
  // and table list number (lid)
  int vid = -1;
  int major = dver.major();
  int minor = dver.minor();
  int extension = dver.extension();

  // if the major verison is not set, find highest major verison
  // if it is set, make sure it is in the DB

  auto const& versions = _valcache.valVersions();

  bool qOK = false;
  if (major < 0) {
    for (auto const& r : versions.rows()) {
      if (r.pid() == pid) {
        if (r.major() > major) major = r.major();
        qOK = true;
      }
    }
  } else {
    for (auto const& r : versions.rows()) {
      if (r.pid() == pid) {
        if (r.major() == major) qOK = true;
      }
    }
  }

  if (major < 0 || !qOK) {
    throw cet::exception("DBVALTOOL_BAD_MAJOR")
        << " DbValTool::fillSet bad calibration major version number" << major
        << "\n";
  }

  // if the minor verison is not set, find highest minor version
  // if it is set, make sure it is in the DB
  qOK = false;
  if (minor < 0) {
    for (auto const& r : versions.rows()) {
      if (r.pid() == pid && r.major() == major) {
        if (r.minor() > minor) {
          minor = r.minor();
          vid = r.vid();
        }
        qOK = true;
      }
    }
  } else {
    for (auto const& r : versions.rows()) {
      if (r.pid() == pid && r.major() == major) {
        if (r.minor() == minor) {
          qOK = true;
          vid = r.vid();
        }
      }
    }
  }

  if (minor < 0 || !qOK) {
    throw cet::exception("DBVALTOOL_BAD_MINOR")
        << " DbValTool::fillSet bad calibration minor version number" << minor
        << "\n";
  }

  // loop over the extensions to this version,
  // to eventually collect groups consistent with
  // with the version number

  auto const& extensions = _valcache.valExtensions();

  // first collect the extension id's, which will go into
  // relational table extensionlists, to get gid's
  std::vector<int> eids;
  int max_extension = -1;
  for (auto const& r : extensions.rows()) {
    // keep if we are accepting all extensions,
    // or up to or equal the requested extension
    if (r.vid() == vid && (extension < 0 || r.extension() <= extension)) {
      eids.push_back(r.eid());
      if (r.extension() > max_extension) max_extension = r.extension();
    }
  }

  // now loop over extensionlists and collect groups (gid's) associated
  // with each extension (eid)
  std::vector<int> gids;
  auto const& extensionlists = _valcache.valExtensionLists();

  for (auto eid : eids) {
    for (auto const& r : extensionlists.rows()) {
      if (r.eid() == eid) gids.push_back(r.gid());
    }
  }

  if (gids.size() == 0) {
    throw cet::exception("DBVALTOOL_NO_EXTENSION")
        << " DbValTool::fillSet found no calibration groups for version "
        << dver.to_string() << "\n";
  }

  // make the list of tables in this purpose/version
  dbset.clear();
  fillSetGid(gids, dbset);

  return;
}

//**************************************************

void DbValTool::fillSetGid(std::vector<int> const& gids, DbSet& dbset) const {
  // take the list of groups and loop over the grouplists
  // which gives IOVs for a group
  // these should be sorted so this code could use that
  auto const& gls = _valcache.valGroupLists();
  auto const& iids = _valcache.valIovs();
  auto const& cids = _valcache.valCalibrations();
  for (auto g : gids) {
    for (auto const& r : gls.rows()) {
      if (r.gid() == g) {
        auto const& irow = iids.row(r.iid());
        auto const& crow = cids.row(irow.cid());
        dbset.add(crow.tid(), irow.cid(), irow.iov());
      }
    }
  }
}

}  // namespace mu2e
