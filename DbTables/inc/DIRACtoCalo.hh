#ifndef DbTables_DIRACtoCalo_hh
#define DbTables_DIRACtoCalo_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class DIRACtoCalo : public DbTable {
  public:

    typedef std::shared_ptr<DIRACtoCalo> ptr_t;
    typedef std::shared_ptr<const DIRACtoCalo> cptr_t;
    constexpr static const char* cxname = {"DIRACtoCalo"};
    
    class Row {
    public:
      Row(int index, uint16_t caloRoId):_index(index),_caloRoId(caloRoId) {}
      int  index() const { return _index;}
      uint16_t caloRoId() const {return _caloRoId;}
    private:
      int _index;
      uint16_t _caloRoId;
    };


    DIRACtoCalo():DbTable(cxname,"calo.diractocalo",
			  "index,caloRoId") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    virtual std::size_t nrowFix() const { return 216; }; 
    size_t size() const { return baseSize() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stoi(columns[1]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<","<<r.caloRoId();
    }

    virtual void clear() override { baseClear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
