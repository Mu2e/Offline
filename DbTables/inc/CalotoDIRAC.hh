#ifndef DbTables_CalotoDIRAC_hh
#define DbTables_CalotoDIRAC_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"
#include "DataProducts/inc/CaloId.hh"

namespace mu2e {

  class CalotoDIRAC : public DbTable {
  public:

    typedef std::shared_ptr<CalotoDIRAC> ptr_t;
    typedef std::shared_ptr<const CalotoDIRAC> cptr_t;
    constexpr static const char* cxname = {"CalotoDIRAC"};
    
    class Row {
    public:
      Row(int offlineID, uint16_t diracID):_offlineID(offlineID),_diracID(diracID) {}
      int       offlineID() const { return _offlineID;}
      uint16_t  diracID()   const { return _diracID;}
    private:
      int      _offlineID;
      uint16_t _diracID;
    };
    
    
    CalotoDIRAC():DbTable(cxname,"calo.calotodirac",
			  "offlineID,diracID") {}
    const Row&              rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows()    const { return _rows;}
    std::size_t             nrow()    const { return _rows.size(); };
    virtual std::size_t     nrowFix() const { return CaloId::_nCrystalChannel; }; 
    size_t                  size()    const { return baseSize() + nrow()*sizeof(Row); };
    
    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stoi(columns[1]) );
    }
    
    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.offlineID()<<","<<r.diracID();
    }
    
    virtual void clear() override { baseClear(); _rows.clear(); }
    
  private:
    std::vector<Row> _rows;
  };
  
};
#endif
