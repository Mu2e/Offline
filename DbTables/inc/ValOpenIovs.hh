#ifndef DbTables_ValOpenIovs_hh
#define DbTables_ValOpenIovs_hh

// intervals of validity for ancillary open IoV system

#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class ValOpenIovs : public DbTable {
 public:
  typedef std::shared_ptr<ValOpenIovs> ptr_t;
  typedef std::shared_ptr<const ValOpenIovs> cptr_t;

  class Row {
   public:
    Row(const std::string& name, int cid, int start_run, int start_subrun,
        int end_run, int end_subrun, const std::string& comment,
        const std::string& create_time, const std::string& create_user) :
        _name(name),
        _cid(cid), _start_run(start_run), _start_subrun(start_subrun),
        _end_run(end_run), _end_subrun(end_subrun), _comment(comment),
        _create_time(create_time), _create_user(create_user) {}
    const std::string& name() const { return _name; }
    int cid() const { return _cid; }
    int start_run() const { return _start_run; }
    int start_subrun() const { return _start_subrun; }
    int end_run() const { return _end_run; }
    int end_subrun() const { return _end_subrun; }
    const std::string& comment() const { return _comment; }
    DbIoV iov() const {
      return DbIoV(start_run(), start_subrun(), end_run(), end_subrun());
    }
    std::string const& create_time() const { return _create_time; }
    std::string const& create_user() const { return _create_user; }

   private:
    std::string _name;
    int _cid;
    int _start_run;
    int _start_subrun;
    int _end_run;
    int _end_subrun;
    std::string _comment;
    std::string _create_time;
    std::string _create_user;
  };

  constexpr static const char* cxname = "ValOpenIovs";

  ValOpenIovs() :
      DbTable(cxname, "val.openiovs",
              "name,cid,start_run,start_subrun,end_run,end_subrun,comment,"
              "create_time,create_user") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
  const std::string orderBy() const { return std::string("create_time"); }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(columns[0], std::stoi(columns[1]), std::stoi(columns[2]),
                       std::stoi(columns[3]), std::stoi(columns[4]),
                       std::stoi(columns[5]), columns[6], columns[7],
                       columns[8]);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.name() << ",";
    sstream << r.cid() << ",";
    sstream << r.start_run() << ",";
    sstream << r.start_subrun() << ",";
    sstream << r.end_run() << ",";
    sstream << r.end_subrun() << ",";
    sstream << r.comment() << ",";
    sstream << r.create_time() << ",";
    sstream << r.create_user();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<Row> _rows;
};

}  // namespace mu2e
#endif
