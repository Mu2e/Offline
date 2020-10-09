#ifndef ConfigTools_src_SimpleConfigRecord_hh
#define ConfigTools_src_SimpleConfigRecord_hh
//
// A class to hold one record within the primitive
// SimpleConfig utility.
//
//
// Contact person Rob Kutschke
//

#include <iosfwd>
#include <string>
#include <vector>

namespace mu2e {

//
// This class holds a single record from the configuration file.
// It parses the record, checks for internal consistency and provides
// accesors for information in the record.
//
// Notes:
// 1) Supported types:
//    string, int, double, bool
//    vector<string>, vector<int>, vector<double>
//
// 2) Supports empty vectors and empty strings.
//
//
// Work list:
// 1) Do the type conversion in c'tor so that we fail immediately
//    and not later on when someone tries to read the value.
//

class SimpleConfigRecord {

 public:

  // Constructor.
  SimpleConfigRecord( const std::string& record_a );

  /**
   * Returns a copy of the input record as it was found in the input file
   * but with line breaks removed.
   *
   * @return A copy of the raw record.
   */
  std::string getRecord() const { return record; }

  /**
   * Returns the type field of this record.
   * @return The type.
   */
  std::string getType () const { return Type; }

  /**
   * Returns the variable name field of this record.
   * @return The variable name.
   */
  std::string getName () const { return Name; }

  /**
   * Returns the comment field, if any, of this record.
   * @return The comment, if any, that is part of this record.
   */
  std::string getComment() const { return comment; }

  /**
   * Returns true if the record contains nothing more than a comment or if the
   * the record is blank.
   *
   * @return true if the record is a pure comment or is blank; false otherwise.
   */
  bool isCommentOrBlank() const { return _isCommentOrBlank;}

  /**
   * Returns true if the record is one of the vector types.
   *
   */
  bool isVector() const { return _isVector; }

  /**
   * Return the number of values in the record; 1 for a scalar type
   * and 0 or more for a vector type.
   *
   */
  int size() const;

  // Accessors to return supported data types.

  /**
   * Return the value as string.  This will work for any data type.
   * This method is used inside the toString method and should not count when
   * called from there.
   * @return The value as a string.
   */
  std::string getString( bool count=true) const;

  /**
   * Return the value as an int.  Only works for recrods that are of type int.
   *
   * @return The value of an int record.
   */
  int getInt (bool count=true) const;

  /**
   * Return the value as an double.  Only works for records that are of type double.
   *
   * @return The value of a double record.
   */
  double getDouble (bool count=true) const;

  /**
   * Return the value as a bool.  Only works for records that are of type bool.
   *
   * @return The value of a bool record.
   */
  bool getBool(bool count=true) const;

  /**
   * Return the value as a vector of strings.  Works for all record types.
   * This method is used inside the toString method and should not count when
   * called from there.
   *
   * @return The value of a the record.
   */
  // Can return any type of vector as a vector of strings.
  void getVectorString( std::vector<std::string>& v, bool count=true) const;

  /**
   * Return the value as a std::vector<int>.
   * Only works for records that are of type std::vector<int>.
   *
   * @return The value of a std::vector<int> record.
   */
  void getVectorInt( std::vector<int>& v, bool count=true) const;

  /**
   * Return the value as a std::vector<Double>.
   * Only works for records that are of type std::vector<Double>.
   *
   * @return return the value of a std::vector<double> record.
   */
  void getVectorDouble( std::vector<double>& v, bool count=true ) const;

  /**
   *
   * Format the record as a string with standard spacing, ignoring the spacing
   * on the input line. If the data are strings, then enclose each string in
   * quotes, even if it has no embedded spaces.
   *
   * @return A formatted copy of the record.
   */
  std::string toString() const;

  /**
   *
   * Formatted printing of the record to the specified stream.
   *
   * @return A formatted copy of the record.
   */
  void print( std::ostream& ) const;

  /**
   *
   * Number of times this record has been accessed.
   *
   * @return the number of times this record has been accessed.
   */
  int accessCount() const { return _accessCount; }
  void clearAccessCount() { _accessCount = 0; }

  bool isSuperceded() const { return _superceded; }

  void setSuperceded() {
    _superceded=true;
  }

private:
  // Private instance data.


  // A copy of the record as it came in.
  // An external class does the concatenation of records that span lines into a single
  // string.
  std::string record;

  // Record with comments and enclosing white space stripped out.
  std::string barerecord;

  // Comment field, if any, including the // delimiter.
  std::string comment;

  // Data type: int, double, bool ...
  std::string Type;

  // Name of the datum.
  std::string Name;

  // The value field - not yet parsed into components.
  std::string Value;

  // The value field, parsed into components.
  std::vector<std::string> Values;

  // State data.
  mutable int _accessCount;
  bool _isCommentOrBlank;
  bool _isVector;
  bool _superceded;

  //  Private methods.

  /**
   *
   * Parse this record, starting from its input string.
   *
   */
  void Parse ();

  /**
   *
   * Split this record into 3 fields: Type Name = Value;
   *
   */
  void SplitLine();

  /**
   *
   * Parse the Value part of the record.
   *
   */
  void ParseValue();

  /**
   *
   * Check that the type of the current record matches the specified type.
   *
   */
  void CheckType( std::string s) const;

  /**
   *
   * Check that the type of the current record matches one of the known types.
   * Throw if it does not.
   *
   */
  void KnownType() const;


  /**
   *
   * Check that the type of the current record is one of the vector types.
   *
   */
  void AnyVector() const;

  /**
   *
   * Return value of string as an int, or throw.
   *
   */
  int toInt ( const std::string& s ) const ;

  /**
   *
   * Return value of string as a double, or throw.
   *
   */
  double toDouble ( const std::string& s ) const;

  /**
   *
   * Return value of string as a bool, or throw.
   *
   */
  bool toBool ( const std::string& s ) const;

};


inline std::ostream& operator<<(std::ostream& ost, const SimpleConfigRecord& r){
  r.print(ost);
  return ost;
}

}

#endif /* ConfigTools_src_SimpleConfigRecord_hh */
