// clang-format off

#ifndef DbTables_TrkPanelMap_hh
#define DbTables_TrkPanelMap_hh
#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TrkPanelMap : public DbTable {
 public:
  typedef std::shared_ptr<TrkPanelMap> ptr_t;
  typedef std::shared_ptr<const TrkPanelMap> cptr_t;

  class Row {
  public:
    Row(int MnID, int DtcID, int Link, int UniquePlane, int PpID, int Panel, int ZFace) {
      _mnid    = MnID;
      _dtc     = DtcID;
      _link    = Link ;
      _uniquePlane   = UniquePlane;
      _ppid    = PpID;
      _panel   = Panel;
      _zface   = ZFace;
    }

    int   mnid         () const { return _mnid;          }
    int   dtc          () const { return _dtc;           }
    int   link         () const { return _link;          }
    int   station      () const { return _uniquePlane/2; }
    int   uniquePlane  () const { return _uniquePlane;   }
    int   ppid         () const { return _ppid;          }
    int   panel        () const { return _panel;         }
    int   zface        () const { return _zface;         }
  private:
    int _mnid;                        // panel production ('Minnesota') ID
    int _dtc;                         // DTC ID, not just PCIE slot ID
    int _link;                        // ROC link of this panel
    int _uniquePlane;                 // 'geo' (or 'unique') plane ID (0-35)
    int _ppid;                        // production plane ID, as in 'plane 21'
    int _panel;                       // 'geo' (offline) panel ID
    int _zface;                       // Z index of a face in the installed station, 0-3
  };

  constexpr static const char* cxname = "TrkPanelMap";

  TrkPanelMap() : DbTable(cxname, "trk.panelmap",
                               "mnid,dtc,link,plane,ppid,panel,zface") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }

  std::vector<Row> const& rows   () const          { return _rows; }
  std::size_t             nrow   () const override { return _rows.size(); };
  // virtual std::size_t     nrowFix() const override { return 0; };  // the number of rows is variable
  size_t                  size   () const override { return baseSize() + nrow() * sizeof(Row); };
  const std::string       orderBy() const          { return std::string("mnid"); }

  void addRow(const std::vector<std::string>& columns) override {
    //    int index = std::stoi(columns[0]);
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]),
                       std::stoi(columns[2]), std::stoi(columns[3]),
                       std::stoi(columns[4]), std::stoi(columns[5]),
                       std::stoi(columns[6]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.mnid()    << ","
            << r.dtc()     << ","
            << r.link()    << ","
            << r.uniquePlane()   << ","
            << r.ppid()    << ","
            << r.panel()   << ","
            << r.zface();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<Row> _rows;
};

}  // namespace mu2e
#endif
