#ifndef DbTables_TrkROCtoPanel_hh
#define DbTables_TrkROCtoPanel_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkROCtoPanel : public DbTable {
  public:

    typedef std::shared_ptr<TrkROCtoPanel> ptr_t;
    typedef std::shared_ptr<const TrkROCtoPanel> cptr_t;

    class Row {
    public:
      Row(int index, uint16_t panel):_index(index),_panel(panel) {}
      int  index() const { return _index;}
      uint16_t panel() const {return _panel;}
    private:
      int _index;
      uint16_t _panel;
    };


    TrkROCtoPanel():DbTable("TrkROCtoPanel","trk.roctopanel",
	 "index,panel") {}
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
      sstream << r.index()<<","<<r.panel();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
