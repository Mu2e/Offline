#ifndef DbTables_TrkDRACtoStraw_hh
#define DbTables_TrkDRACtoStraw_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkDRACtoStraw : public DbTable {
  public:

    typedef std::shared_ptr<TrkDRACtoStraw> ptr_t;
    typedef std::shared_ptr<const TrkDRACtoStraw> cptr_t;

    class Row {
    public:
      Row(int index, uint16_t straw):_index(index),_straw(straw) {}
      int  index() const { return _index;}
      uint16_t straw() const {return _straw;}
    private:
      int _index;
      uint16_t _straw;
    };


    TrkDRACtoStraw():DbTable("TrkDRACtoStraw","trk.dractostraw",
	 "index,straw") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    virtual std::size_t nrowFix() const { return 216; }; 
    size_t size() const { return _csv.capacity() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stof(columns[1]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<","<<r.straw();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
