#ifndef DbTables_CalChannels_hh
#define DbTables_CalChannels_hh


#include <string>
#include <iomanip>
#include <sstream>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/DataProducts/inc/CaloRawSiPMId.hh"

namespace mu2e {

  class CalChannels : public DbTable {
  public:
    typedef std::shared_ptr<CalChannels> ptr_t;
    typedef std::shared_ptr<const CalChannels> cptr_t;

    class Row {
    public:
      Row(CaloRawSiPMId rawid, CaloSiPMId roid):
        _rawid(rawid),
        _roid(roid){}
      CaloRawSiPMId  rawid() const { return _rawid;}
      CaloSiPMId  roid() const { return _roid;}

    private:
      CaloRawSiPMId _rawid;
      CaloSiPMId _roid;
    };

    constexpr static const char* cxname = "CalChannels";

    CalChannels():DbTable(cxname,"cal.channels","rawid,roid"){}

    const Row& row(CaloRawSiPMId  rawid) const {
                return _rows.at(rawid.id()); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize()  + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nRawChannel; };
    const std::string orderBy() const { return std::string("rawid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
      // enforce order, so channels can be looked up by index
      if (index >= CaloConst::_nRawChannel  || index != _rows.size()) {
        throw cet::exception("CALOCHANNELS_BAD_INDEX")
        << "CalChannels::addRow found index out of order: "
        <<index<< " != " << _rows.size() <<"\n";

      }
   _rows.emplace_back(CaloRawSiPMId(index),
                      CaloSiPMId(std::stoul(columns[1])));

  }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.rawid() <<",";
      sstream << r.roid() ;
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
  };
}
#endif
