#ifndef DbTables_DbValCache_hh
#define DbTables_DbValCache_hh

#include <iostream>
#include <iomanip>
#include <vector>
#include "DbTables/inc/ValTables.hh"
#include "DbTables/inc/ValCalibrations.hh"
#include "DbTables/inc/ValIovs.hh"
#include "DbTables/inc/ValGroups.hh"
#include "DbTables/inc/ValGroupLists.hh"
#include "DbTables/inc/ValPurposes.hh"
#include "DbTables/inc/ValLists.hh"
#include "DbTables/inc/ValTableLists.hh"
#include "DbTables/inc/ValVersions.hh"
#include "DbTables/inc/ValExtensions.hh"
#include "DbTables/inc/ValExtensionLists.hh"

namespace mu2e {

  class DbValCache {
  public:
    DbValCache():_verbose(0) { }

    void clear();
    size_t size();
    void print();

    ValTables const& valTables() const {return _tables;}
    ValCalibrations const& valCalibrations() const {return _calibrations;}
    ValIovs const& valIovs() const {return _iovs;}
    ValGroups const& valGroups() const {return _groups;}
    ValGroupLists const& valGroupLists() const {return _grouplists;}
    ValPurposes const& valPurposes() const {return _purposes;}
    ValLists const& valLists() const {return _lists;}
    ValTableLists const& valTableLists() const {return _tablelists;}
    ValVersions const& valVersions() const {return _versions;}
    ValExtensions const& valExtensions() const {return _extensions;}
    ValExtensionLists const& valExtensionLists() const {return _extensionlists;}

    DbTable const& asTable(std::string const& tableName) const;

    void setVerbose(int v) { _verbose = v; }
    void setValTables(ValTables const& tables) { _tables = tables; }
    void setValCalibrations(ValCalibrations const& calibrations) { 
      _calibrations = calibrations; }
    void setValIovs(ValIovs const& iovs) { _iovs = iovs; }
    void setValGroups(ValGroups const& groups) { _groups = groups; }
    void setValGroupLists(ValGroupLists const& grouplists) { 
                   _grouplists = grouplists; }
    void setValPurposes(ValPurposes const& purposes) { _purposes = purposes; }
    void setValLists(ValLists const& lists ) {_lists = lists;}
    void setValTableLists(ValTableLists const& tablelists ) {
      _tablelists = tablelists;}
    void setValVersions(ValVersions const& versions ) {_versions = versions;}
    void setValExtensions(ValExtensions const& extensions ) {
      _extensions = extensions;}
    void setValExtensionLists(ValExtensionLists const& extensionlists ) {
      _extensionlists = extensionlists;}

  private:
    int _verbose;

    ValTables _tables;
    ValCalibrations _calibrations;
    ValIovs _iovs;
    ValGroups _groups;
    ValGroupLists _grouplists;
    ValPurposes _purposes;
    ValLists _lists;
    ValTableLists _tablelists;
    ValVersions _versions;
    ValExtensions _extensions;
    ValExtensionLists _extensionlists;

  };

}

#endif
