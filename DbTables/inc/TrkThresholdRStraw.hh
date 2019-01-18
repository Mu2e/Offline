#ifndef DbTables_TrkThresholdRStraw_hh
#define DbTables_TrkThresholdRStraw_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkThresholdRStraw : public DbTable {
  public:

    typedef std::shared_ptr<TrkThresholdRStraw> ptr_t;
    typedef std::shared_ptr<const TrkThresholdRStraw> cptr_t;

    class Row {
    public:
      Row(int index, float threshold_hv, float threshold_cal):
	_index(index),
	_threshold_hv(threshold_hv),_threshold_cal(threshold_cal) {}
      int  index() const { return _index;}
      float thresholdHv() const {return _threshold_hv;}
      float thresholdCal() const {return _threshold_cal;}
    private:
      int _index;
      float _threshold_hv;
      float _threshold_cal;
    };


    TrkThresholdRStraw():DbTable("TrkThresholdRStraw","trk.thresholdrstraw",
	 "index,threshold_hv,threshold_cal") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    virtual std::size_t nrowFix() const { return 96; }; 
    size_t size() const { return _csv.capacity() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stof(columns[1]),
			 std::stof(columns[2]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << std::fixed << std::setprecision(1); 
      sstream << r.thresholdHv()<<",";
      sstream << r.thresholdCal();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
