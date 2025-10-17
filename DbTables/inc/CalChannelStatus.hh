#ifndef DbTables_CalChannelStatus_hh
#define DbTables_CalChannelStatus_hh

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalChannelStatus : public DbTable {
  public:
  typedef std::shared_ptr<CalChannelStatus> ptr_t;
  typedef std::shared_ptr<const CalChannelStatus> cptr_t;

    class Row {
    public:
      Row(CaloSiPMId roid, std::string Status):_roid(roid),_status(Status) {}
      CaloSiPMId   roid() const { return _roid;}
      std::string Status() const { return _status; }

    private:
      CaloSiPMId   _roid;
      std::string _status; //"good", "dead", "hot", etc..
    };

    constexpr static const char* cxname = "CalChannelStatus";

    CalChannelStatus():DbTable(cxname,"cal.channelstatus","roid,status"){}

    const Row& row(CaloSiPMId id) const {
                return _rows[id.id()];
    }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannelDB; };
    const std::string orderBy() const { return std::string("roid"); }

    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index!=int(_rows.size())) {
        throw cet::exception("CALOCHANNELSTATUS_BAD_INDEX")<<"CalChannelStatus::addRow found index out of order:"<<index << " != " << int(_rows.size()) <<"\n";
      }
       _rows.emplace_back(CaloSiPMId(index),columns[1]);

    }


    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.roid()<<",";
      sstream << r.Status();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;

  };

}
#endif
