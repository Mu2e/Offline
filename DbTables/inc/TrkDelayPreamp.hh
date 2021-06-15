#ifndef DbTables_TrkDelayPreamp_hh
#define DbTables_TrkDelayPreamp_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkDelayPreamp : public DbTable {
  public:

    typedef std::shared_ptr<TrkDelayPreamp> ptr_t;
    typedef std::shared_ptr<const TrkDelayPreamp> cptr_t;

    class Row {
    public:
      Row(int index, float delay_hv, float delay_cal):_index(index),
        _delay_hv(delay_hv), _delay_cal(delay_cal) {}
      int  index() const { return _index;}
      float delayHv() const {return _delay_hv;}
      float delayCal() const {return _delay_cal;}
    private:
      int _index;
      float _delay_hv;
      float _delay_cal;
    };

    constexpr static const char* cxname = "TrkDelayPreamp";

    TrkDelayPreamp():DbTable(cxname,"trk.delaypreamp",
	 "index,delay_hv,delay_cal") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    virtual std::size_t nrowFix() const override { return 96; }; 
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stof(columns[1]),
                         std::stof(columns[2]));
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << std::fixed << std::setprecision(1);
      sstream << r.delayHv()<<",";
      sstream << r.delayCal();
    }

    virtual void clear() override { baseClear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
