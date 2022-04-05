#ifndef DbTables_TrkThresholdRStraw_hh
#define DbTables_TrkThresholdRStraw_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"

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

    constexpr static const char* cxname = "TrkThresholdRStraw";

    TrkThresholdRStraw():DbTable(cxname,"trk.thresholdrstraw",
         "index,threshold_hv,threshold_cal") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    virtual std::size_t nrowFix() const override { return 96; };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    const std::string orderBy() const {return std::string("index");}

    void addRow(const std::vector<std::string>& columns) override {
      int index = std::stoi(columns[0]);
      // enforce a strict sequential order
      if(index!=int(_rows.size())) {
        throw cet::exception("TRKTHRESHOLDRSTRAW_BAD_INDEX")
          << "TrkThresholdRStraw::addRow found index out of order: "
          <<index << " != " << _rows.size() <<"\n";
      }
      _rows.emplace_back(index,
                         std::stof(columns[1]),
                         std::stof(columns[2]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.index()<<",";
      sstream << std::fixed << std::setprecision(1);
      sstream << r.thresholdHv()<<",";
      sstream << r.thresholdCal();
    }

    virtual void clear() override { baseClear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };

};
#endif
