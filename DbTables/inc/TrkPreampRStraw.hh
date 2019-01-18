#ifndef DbTables_TrkPreampRStraw_hh
#define DbTables_TrkPreampRStraw_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkPreampRStraw : public DbTable {
  public:

    typedef std::shared_ptr<TrkPreampRStraw> ptr_t;
    typedef std::shared_ptr<const TrkPreampRStraw> cptr_t;

    class Row {
    public:
      Row(int index, float delay_hv, float delay_cal,
	  float threshold_hv, float threshold_cal, float gain):
	_index(index),_delay_hv(delay_hv),_delay_cal(delay_cal),
	_threshold_hv(threshold_hv),_threshold_cal(threshold_cal),
	_gain(gain) {}
      int  index() const { return _index;}
      float delayHv() const {return _delay_hv;}
      float delayCal() const {return _delay_cal;}
      float thresholdHv() const {return _threshold_hv;}
      float thresholdCal() const {return _threshold_cal;}
      float gain() const {return _gain;}
    private:
      int _index;
      float _delay_hv;
      float _delay_cal;
      float _threshold_hv;
      float _threshold_cal;
      float _gain;
    };


    TrkPreampRStraw():DbTable("TrkPreampRStraw","trk.preamprstraw",
	 "index,delay_hv,delay_cal,threshold_hv,threshold_cal,gain") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    virtual std::size_t nrowFix() const { return 96; }; 
    size_t size() const { return _csv.capacity() + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),
			 std::stof(columns[1]),
			 std::stof(columns[2]),
			 std::stof(columns[3]),
			 std::stof(columns[4]),
			 std::stof(columns[5]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << std::fixed << std::setprecision(1); 
      sstream << r.delayHv()<<",";
      sstream << r.delayCal()<<",";
      sstream << r.thresholdHv()<<",";
      sstream << r.thresholdCal()<<",";
      sstream << r.gain();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
