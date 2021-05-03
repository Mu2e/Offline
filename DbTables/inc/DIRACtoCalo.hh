#ifndef DbTables_DIRACtoCalo_hh
#define DbTables_DIRACtoCalo_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"
#include "DataProducts/inc/CaloId.hh"

namespace mu2e {

  class DIRACtoCalo : public DbTable {
  public:

    typedef std::shared_ptr<DIRACtoCalo> ptr_t;
    typedef std::shared_ptr<const DIRACtoCalo> cptr_t;
    constexpr static const char* cxname = {"DIRACtoCalo"};
    
    class Row {
    public:
      Row(int diracID, uint16_t caloRoID):_diracID(diracID),_caloRoID(caloRoID) {}
      int        diracID()  const { return _diracID;}
      uint16_t   caloRoID() const { return _caloRoID;}
    private:
      int      _diracID;
      uint16_t _caloRoID;
    };


    DIRACtoCalo():DbTable(cxname,"calo.diractocalo",
			  "diracID,caloRoID") {}
    const Row&              rowAt(const std::size_t diracID) const { return _rows.at(diracID);}
    std::vector<Row> const& rows()    const { return _rows;}
    std::size_t             nrow()    const { return _rows.size(); };
    virtual std::size_t     nrowFix() const { return CaloId::_nTotChannel; }; 
    size_t                  size()    const { return baseSize() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      int channel = std::stoi(columns[0]);
      // enforce a strict sequential order - optional
      if(channel!=int(_rows.size())) {
	throw cet::exception("DIRACtoCalo_BAD_INDEX") 
	  << "DIRACtoCalo::addRow found index out of order: " 
	  << channel << " != " << _rows.back().caloRoID()+1 <<"\n";
      }
      _rows.emplace_back(channel,
			 std::stoi(columns[1]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.diracID()<<","<<r.caloRoID();
    }

    virtual void clear() override { baseClear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
