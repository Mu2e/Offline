//
// tables for tracker element status.  An
// element is a straw, panel, or plane
//  These tables list just defects, so are variable length
//  Original author: Dave Brown (LBNL)

#ifndef DbTables_TrkElementStatus_hh
#define DbTables_TrkElementStatus_hh

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include <regex>
#include "cetlib_except/exception.h"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "DataProducts/inc/StrawStatus.hh"
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkElementStatus : public DbTable {
  public:

    // define a struct for each table row
    struct TrkElementStatusRow {
      StrawId _sid;
      StrawStatus _status;
      StrawId const& id() const { return _sid; }
      StrawStatus const& status() const { return _status; }
      TrkElementStatusRow(StrawId const& sid, StrawStatus const& status) : _sid(sid), _status(status) {}
    };

    typedef std::shared_ptr<TrkElementStatus> ptr_t;
    typedef std::shared_ptr<const TrkElementStatus> cptr_t;

    TrkElementStatus(const char* Name, const char* DbName, StrawIdMask sidmask, StrawStatus statusmask):
      DbTable(Name,DbName, "strawid,strawstatus"), _sidmask(sidmask), _statusmask(statusmask) {}
    const TrkElementStatusRow& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<TrkElementStatusRow> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    // this is a variable-size table, so don't overrido nrowFix()
    size_t size() const override { return baseSize() + nrow()*sizeof(TrkElementStatusRow); }
    // table-specific info
    StrawIdMask const& sidMask() const { return _sidmask; }
    StrawStatus const& statusMask() const { return _statusmask; }
    // build from text table format
    void addRow(const std::vector<std::string>& columns) override {
      auto sid = StrawId(columns[0]);
      StrawStatus status(columns[1]);
      // verify the status is allowed
      if(!_statusmask.hasAllProperties(status))
	throw cet::exception(name()) << "Illegal status specified " << status << std::endl;
      _rows.emplace_back(TrkElementStatusRow(sid,status));
    }
    // printout, used to fill db content (?)
    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      TrkElementStatusRow const& r = _rows.at(irow);
      sstream << r.id().plane() << "_" << r.id().panel() << "_" << r.id().straw() << ","; 
      sstream << r.status().hex(); // should this be the text string??
    }

    void clear() override { baseClear(); _rows.clear(); }

  private:
    StrawIdMask _sidmask; // defines matching for this element
    StrawStatus _statusmask; // status bits allowed for this element
    std::vector<TrkElementStatusRow> _rows;
  };

// unique classes for Db usage.  Each defines a mask and for StrawId and the StrawStatus
  class TrkPanelStatus : public TrkElementStatus {
    public:
      constexpr static const char* cxname = "TrkPanelStatus";
      TrkPanelStatus() : TrkElementStatus(cxname,"trk.panelstatus", StrawIdMask("uniquepanel"), StrawStatus("Absent:NoHV:NoGas:NoLV:LowGasGain")) {}
  };

  class TrkPlaneStatus : public TrkElementStatus {
    public:
      constexpr static const char* cxname = "TrkPlaneStatus";
      TrkPlaneStatus() : TrkElementStatus(cxname,"trk.planestatus", StrawIdMask("plane"), StrawStatus("Absent")) {}
  };
 // split individual straw status tables in 2: one for short-term, one for long-term 
  class TrkStrawStatusShort : public TrkElementStatus {
    public:
      constexpr static const char* cxname = "TrkStrawStatusShort";
      TrkStrawStatusShort() : TrkElementStatus(cxname,"trk.strawstatusshort", StrawIdMask("uniquestraw"), StrawStatus("Sparking:Suppress:Noise:Pickup")) {}
  };
  class TrkStrawStatusLong : public TrkElementStatus {
    public:
      constexpr static const char* cxname = "TrkStrawStatusLong";
      TrkStrawStatusLong() : TrkElementStatus(cxname,"trk.strawstatuslong", StrawIdMask("uniquestraw"), StrawStatus("Absent:NoWire:NoHV:NoPreamp:NoADC:NoTDC")) {}
  };
  
};
#endif
