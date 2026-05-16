#ifndef DbTables_TrkHVTripStatus_hh
#define DbTables_TrkHVTripStatus_hh

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/StrawStatus.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <regex>
#include <sstream>
#include <string>

namespace mu2e {

class TrkHVTripStatus : public DbTable {
 public:

   class Row {
     public:
       Row(int index, StrawId sid, uint32_t startevent, uint32_t endevent) :
         _index(index), _sid(sid), _startevent(startevent), _endevent(endevent) {}
       int index() const { return _index; }
       StrawId const& id() const { return _sid; }
       uint32_t startevent() const { return _startevent; }
       uint32_t endevent() const { return _endevent; }

     private:
       int _index;
       StrawId _sid;
       uint32_t _startevent;
       uint32_t _endevent;
   };



  typedef std::shared_ptr<TrkHVTripStatus> ptr_t;
  typedef std::shared_ptr<const TrkHVTripStatus> cptr_t;

  constexpr static const char* cxname = "TrkHVTripStatus";

  TrkHVTripStatus() :
      DbTable(cxname, "trk.hvtripstatus",
              "index,strawid,startevent,endevent") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  // this is a variable-size table, so don't overrido nrowFix()
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
  const std::string orderBy() const { return std::string("index"); }

  void addRow(const std::vector<std::string>& columns) override {
    int index = std::stoi(columns[0]);
    auto sid = StrawId(columns[1]);
    uint32_t startevent = std::stoul(columns[2]);
    uint32_t endevent = std::stoul(columns[3]);
    // enforce a strict sequential order
    if (index != int(_rows.size())) {
      throw cet::exception("TRKHVTRIPSTATUS_BAD_INDEX")
          << "TrkHVTripStatus::addRow found index out of order: " << index
          << " != " << _rows.size() << "\n";
    }
    _rows.emplace_back(index, sid, startevent, endevent);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.index() << ",";
    sstream << r.id().plane() << "_" << r.id().panel() << "_" << r.id().straw()
            << ",";
    sstream << r.startevent() << ",";
    sstream << r.endevent();
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
