#ifndef DbTables_DbTable_hh
#define DbTables_DbTable_hh

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mu2e {

class DbTable {
 public:
  typedef std::shared_ptr<mu2e::DbTable> ptr_t;
  typedef std::shared_ptr<const mu2e::DbTable> cptr_t;
  enum tableType {Calibration=0,Validity=1,Adhoc=2};

  DbTable(const char* name = "DbTable", const char* dbname = "dne.dbtable",
          const char* query = "noquery") :
      _name(name),
      _dbname(dbname), _query(query) {}

  virtual ~DbTable() = default;

  // table name as it appears offline
  const std::string& name() const { return _name; }
  // table name as it appears in the db
  const std::string& dbname() const { return _dbname; }
  // the column names, written as in the db
  const std::string& query() const { return _query; }
  // the table data in string format
  const std::string& csv() const { return _csv; }
  // number of rows - overridden by derived class
  virtual std::size_t nrow() const = 0;
  // expected nrows - overridden by derived class
  virtual std::size_t nrowFix() const { return 0; }
  // approx size in bytes - overridden by derived class
  virtual std::size_t size() const = 0;
  std::size_t baseSize() const { return _csv.capacity(); }
  // if table needs to be sorted
  virtual const std::string orderBy() const { return std::string(); }
  virtual tableType type() const { return Calibration; }

  // take the cvs text from a query and build out the table contents
  int fill(const std::string& csv, bool saveCsv = true);
  // in case table was filled with binary values, convert to csv
  int toCsv();

  // part of building content, convert list of strings to binary row
  virtual void addRow(const std::vector<std::string>& columns) = 0;
  // convert a row in a binary format to a string
  virtual void rowToCsv(std::ostringstream& stream, size_t irow) const = 0;
  // remove all rows
  virtual void clear() = 0;
  void baseClear() { _csv.clear(); }

 private:
  std::string _name;
  std::string _dbname;
  std::string _query;
  std::string _csv;
};

}  // namespace mu2e
#endif
