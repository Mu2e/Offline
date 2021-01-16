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
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "DataProducts/inc/StrawStatus.hh"
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkElementStatus : public DbTable {
  public:

    // define a struct for each table row
    struct TrkElementStatusRow {
      int _index;//  I'm unclear what value this provides to a variable-length table, but the interface requires this
      StrawId _sid;
      StrawStatus _status;
      int index() const { return _index; }
      StrawId const& id() const { return _sid; }
      StrawStatus const& status() const { return _status; }
      TrkElementStatusRow(int index, StrawId const& sid, StrawStatus const& status) : _index(index), _sid(sid), _status(status) {}
    };

    typedef std::shared_ptr<TrkElementStatus> ptr_t;
    typedef std::shared_ptr<const TrkElementStatus> cptr_t;

    TrkElementStatus(const char* Name, const char* DbName, StrawIdMask sidmask, StrawStatus statusmask):
      DbTable(Name,DbName, "index,strawid,strawstatus"), _sidmask(sidmask), _statusmask(statusmask) {}
    const TrkElementStatusRow& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<TrkElementStatusRow> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    // this is a variable-size tabale, so no nrowFix()
    size_t size() const override { return baseSize() + nrow()*sizeof(TrkElementStatusRow); }
    // table-specific info
    StrawIdMask const& sidMask() const { return _sidmask; }
    StrawStatus const& statusMask() const { return _statusmask; }

    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(std::stoi(columns[0]),
			 StrawId(columns[1]),
			 StrawStatus(columns[2]));
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      TrkElementStatusRow const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << r.id().plane() << "_" << r.id().panel() << "_" << r.id().straw() << ",";
      sstream << r.status().hex();
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
      TrkPanelStatus() : TrkElementStatus(cxname,"trk.panelstatus", StrawIdMask("panel"), StrawStatus("Absent:NoHV:NoGas:NoLV:LowGasGain")) {}
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
      TrkStrawStatusShort() : TrkElementStatus(cxname,"trk.strawstatusshort", StrawIdMask("straw"), StrawStatus("Sparking:Suppress:Noise:Pickup")) {}
  };
  class TrkStrawStatusLong : public TrkElementStatus {
    public:
      constexpr static const char* cxname = "TrkStrawStatusLong";
      TrkStrawStatusLong() : TrkElementStatus(cxname,"trk.strawstatuslong", StrawIdMask("straw"), StrawStatus("Absent:NoWire:NoHV:NoPreamp:NoADC:NoTDC")) {}
  };
  
};
#endif
