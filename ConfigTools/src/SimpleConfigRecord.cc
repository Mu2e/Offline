//
// A class to hold one record within the primitive SimpleConfig utility.
//
// $Id: SimpleConfigRecord.cc,v 1.3 2012/08/06 19:09:38 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/08/06 19:09:38 $
//
// Contact person Rob Kutschke

#include "ConfigTools/src/SimpleConfigRecord.hh"
#include "GeneralUtilities/inc/trimInPlace.hh"

#include "cetlib_except/exception.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

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
  // 1) Warnings for:
  //    int i = 1.23;
  // 2)

  using namespace std;

  // Constructor.
  SimpleConfigRecord::SimpleConfigRecord( const string& record_a ):
    record(record_a),
    _accessCount(0),
    _isCommentOrBlank(false),
    _isVector(false),
    _superceded(false){
    Parse();
  }

  // Accessors to return supported data types.
  string SimpleConfigRecord::getString( bool count) const {
    if ( count ){
      ++_accessCount;
    }
    return Values.at(0);
  }

  int SimpleConfigRecord::getInt (bool count) const {
    if ( count ){
      ++_accessCount;
    }
    CheckType("int");
    return toInt( Values.at(0) );
  }

  double SimpleConfigRecord::getDouble (bool count) const {
    if ( count ){
      ++_accessCount;
    }
    CheckType("double");
    return toDouble( Values.at(0) );
  }

  bool SimpleConfigRecord::getBool(bool count) const {
    if ( count ){
      ++_accessCount;
    }
    CheckType("bool");
    return toBool( Values.at(0) );
  }

  void SimpleConfigRecord::getVectorString( vector<string>& v, bool count) const{
    if ( count ) {
      ++_accessCount;
    }
    AnyVector();
    for ( size_t i=0; i<Values.size(); ++i){
      v.push_back(Values[i]);
    }
  }

  void SimpleConfigRecord::getVectorInt( vector<int>& v, bool count) const {
    if ( count ){
      ++_accessCount;
    }
    CheckType("vector<int>");
    vector<string>::const_iterator b = Values.begin();
    vector<string>::const_iterator e = Values.end();
    for ( ; b!=e; ++b ){
      // ?? Why can't I write: V.push_back( toInt(*b) );
      // Probably an issue of temporaries?
      string s(*b);
      v.push_back(toInt(s));
    }
  }

  void SimpleConfigRecord::getVectorDouble( vector<double>& V, bool count) const {
    if ( count ){
      ++_accessCount;
    }
    CheckType("vector<double>");
    vector<string>::const_iterator b = Values.begin();
    vector<string>::const_iterator e = Values.end();
    for ( ;b!=e; ++b ){
      string s(*b);
      V.push_back(toDouble(s));
    }
  }


  /**
   *
   * Format the record as a string with standard spacing, ignoring the spacing
   * on the input line. If the data are strings, then enclose each string in
   * quotes, even if it has no embedded spaces.
   *
   *
   * @return A formatted copy of the record.
   */
  string SimpleConfigRecord::toString() const {

    if ( _isCommentOrBlank ){
      return comment;
    }

    ostringstream s;
    s << Type << " " << Name << " = ";

    //  try{

    if ( _isVector ) {
      s << "{ ";
      if ( Type == "vector<int>" ){
        vector<int> V;
        getVectorInt(V,false);
        vector<int>::const_iterator b = V.begin();
        vector<int>::const_iterator e = V.end();
        bool first = true;
        for ( ;b!=e; ++b){
          if( !first ){
            s << ", ";
          } else{
            first = false;
          }
          s << *b;
        }
      } else if ( Type == "vector<double>"){
        vector<double> V;
        getVectorDouble(V,false);
        vector<double>::const_iterator b = V.begin();
        vector<double>::const_iterator e = V.end();
        bool first = true;
        for ( ;b != e; ++b ){
          if( !first ){
            s << ", ";
          } else{
            first = false;
          }
          s << *b;
        }
      } else {

        vector<string> V;
        getVectorString(V,false);
        vector<string>::const_iterator b = V.begin();
        vector<string>::const_iterator e = V.end();

        bool first = true;
        for ( ;b != e; ++b ){
          if( !first ){
            s << ", \"";
          } else{
            s << "\"";
            first = false;
          }
          s << *b;
          s << "\"";
        }
      }
      s << " }";

    }else{

      // Scalar types.
      if ( Type == "int" ){
        s << getInt(false);
      } else if ( Type == "double" ){
        s << getDouble(false);
      } else if ( Type == "bool"){
        if ( getBool(false) ){
          s << "true";
        }else{
          s << "false";
        }
      } else {
        s << "\"";
        s << getString(false);
        s << "\"";
      }
    }
    // } catch (SimpleConfigException e) {
    //System.out.println (e.getMessage() );
    //s.append(" [Error formating this item], ");
    //}
    s << ";";
    return s.str();

  }

  void SimpleConfigRecord::print( ostream& ost) const{
    ost << toString();
  }

  //  Private methods.

  // Parse this record, starting from its input string.
  void SimpleConfigRecord::Parse () {

    // Find comment delimiter if present.
    string::size_type islash = record.find("//");

    // Separate the comment from the main body of the record.
    string tmp;
    if ( islash != string::npos ) {
      comment = record.substr(islash);
      tmp = record.substr(0,islash);
    }else{
      tmp = record;
    }
    trimInPlace(tmp);

    // Line is blank or contains only a comment.
    if ( tmp.size() == 0 ){
      _isCommentOrBlank = true;
      return;
    }

    // Check for syntax of a complete record.
    if ( tmp[tmp.size()-1] != ';' ){
      // Test: fail01.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: Not terminated by semicolon: "
        << record;
    }

    // Strip the trailing semicolon and the leading and trailing whitespace.
    barerecord = tmp.substr(0,tmp.size()-1);
    trimInPlace(barerecord);

    // Split the line into: type name = value;
    SplitLine();

    // Parse the value part of the record.
    ParseValue();

    // Check that the type is one of the known types.
    KnownType();

  }


  /**
   *
   * Split this record into 3 fields: Type Name = Value;
   *
   */
  void SimpleConfigRecord::SplitLine(){

    // Is there an equals sign with enough space before and after it?
    // The minimal line is:
    // t n=v
    // where,
    // t = type
    // n = name
    // v = value
    // the space between t and n is significant.
    // No embedded whitespace allowed within t or n.
    // So the first legal spot for the equals sign is:
    //   - must be at index 3 or greater
    //   - the first equals sign in the line must not be the last
    //      non-whitespace character in the line.
    // Remember that barerecord has leading and trailing spaces trimmed.
    size_t iequal = barerecord.find('=');
    if ( iequal < 3 || iequal+1 >= barerecord.size() ){

      // Test: fail04.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: Impossible position for equals sign in record: "
        << record;
    }

    // The value part of the field, to be parsed elsewhere.
    Value = barerecord.substr(iequal+1);
    trimInPlace(Value);

    // Temporary containing the type and name.
    string tmp = barerecord.substr(0,iequal);
    trimInPlace(tmp);

    // White space that separates the type and name.
    size_t white = tmp.find_first_of(" \t");
    if ( white == string::npos ||
         white < 1             ||
         white+2 > tmp.size() ){

      // Test: fail05.conf, fail06.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: Missing space between type and name: "
        << record;
    }

    // Data type.
    Type = tmp.substr(0,white);

    // Variable name.
    for ( string::size_type i = white+1;
          i < tmp.size(); ++i ){
      char c = tmp[i];
      if ( c != ' ' && c != '\t' ){
        Name = tmp.substr(i);
        break;
      }
    }

    // There should be no more white space;
    // Trailing white space has already been trimmed.
    string::size_type check = Name.find_first_of(" \t");
    if ( check != string::npos ){

      // Test: fail07.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: Too many fields before equals sign in record: "
        << record;
    }

  }

  /**
   *
   * Parse the Value part of the record.
   *
   */
  void SimpleConfigRecord::ParseValue(){

    // Check for a record that is a vector.
    _isVector = ( Type.find("vector<") != string::npos ) ? true: false;

    // Part of the string to parse.
    // Default is for non-vectors.
    string::size_type iopen  = 0;
    string::size_type iclose = Value.size();

    // If this is a vector, strip the enclosing {}.
    if ( _isVector ){
      iopen  = Value.find("{");
      iclose = Value.find_last_of("}");
      if ( ( iopen  == string::npos ) ||
           ( iclose == string::npos ) ||
           ( iclose < (iopen+1) ) ){

        // Test: fail08.conf, fail09.conf
        throw cet::exception("SimpleConfig")
          << "SimpleConfigRecord: "
          << "Missing or malformed {} in list of values for a vector: "
          << record;
      }
      iopen += 1;
    }

    // Remove {}, if present, and any leading and trailing whitespace.
    // What remains is the part of the line that holds the comma separated list of values.
    string listpart = Value.substr(iopen,iclose-iopen);
    trimInPlace(listpart);

    // Accumulate the next value in this variable.
    //stringBuffer next = new stringBuffer();
    string next;

    // Some predefined characters that have special meaning.
    char quote  = '\"';
    char bslash = '\\';
    char comma  = ',';

    // States:
    // 0: not within any value
    // 1: within a value that is not started by a quotation mark.
    // 2: within a value that is started by a quotation mark.
    // 3: into white space delimitation but have not yet found comma.
    // 10: last character was an escape, otherwise in state 0
    // 11: last character was an escape, otherwise in state 1
    // 12: last character was an escape, otherwise in state 2
    // 13: last character was an escape, otherwise in state 3

    int state(0);
    for ( string::size_type i=0; i<listpart.size(); ++i ){

      // Next character, both as a char and as a string.
      char c = listpart[i];
      string sc(1,c);

      // If this character was escaped, add it to the next field
      // and drop out of escape mode.
      if ( state > 9 ){
        next.append(sc);
        state = state - 10;
        continue;

      }
      // Not within an value in the list of values.
      else if ( state == 0 ){

        // Skip white space between tokens.
        if ( isspace(c) ) {
          continue;

        }

        // Starting a new item with a quote, strip the quote.
        else if ( c == quote ){
          state = 2;
          continue;

        }

        // Starting a new item with a escape.
        else if ( c == bslash  ){
          next.append(sc);
          state = 11;
          continue;

        }

        // Consecutive commas add an empty item to the vector.
        // State stays at 0.
        else if ( c == comma && next.size() == 0 ){
          Values.push_back("");
          continue;

        }
        // Starting a new item with neither escape nor quote.
        else {
          next.append(sc);
          state = 1;
          continue;
        }

      }

      // In the middle of a field not started by a quote.
      else if ( state == 1 ) {

        // End of item is marked by white space.
        if ( isspace(c) ){
          state = 3;
          continue;
        }

        // End of item marked by comma without preceeding whitespace.
        else if ( c == comma ){
          Values.push_back(next);
          next.clear();
          state = 0;
          continue;
        }

        // Do not allow an unescaped quote in mid word ...
        else if ( c == quote ){

          // Test: fail10.conf
          throw cet::exception("SimpleConfig")
            << "SimpleConfigRecord: "
            << "Unexpected \" character in record: "
            << record;
        }

        // Next character is to be escaped.
        else if ( c == bslash ) {
          next.append(sc);
          state = 11;
          continue;

        }

        // Not a special character, just add it to the string.
        else{
          next.append(sc);
          continue;
        }

      }

      // In the middle of a field started by a quote.
      else if ( state == 2 ){


        // Terminal quote marks the end of an item.
        if ( c == quote ){
          state = 3;

        }

        // Escape the next character.
        else if ( c == bslash ) {
          next.append(sc);
          state = 12;
          continue;

        }

        // Add character. Comma and white space are normal characters in this case.
        else {
          next.append(sc);
          continue;
        }

      }

      // Finished with a field but have not yet seen a comma or end of record.
      else if ( state == 3 ){

        // Skip white space between fields.
        if ( isspace(c) ) {
          continue;

        }

        // Add the previous item to the output vector.
        else if ( c == comma ) {
          Values.push_back(next);
          next.clear();
          state = 0;
          continue;
        }

        else{
          // Test: fail12.conf
          throw cet::exception("SimpleConfig")
            << "SimpleConfigRecord: "
            << "Unexpected white space inside an item in the value list: "
            << record;
        }
      }

      // There should be no way to reach this else.
      else{
        // No test for this error.
        throw cet::exception("SimpleConfig")
          << "SimpleConfigRecord: "
          << "Logic bug in the SimpleConfigRecord value parser: "
          << record;

      }// end main branch, starting with: "if ( state > 9 )"

    } // end loop over characters


    // Record ended with an unterminated quote
    if ( state == 2 || state == 12 ){

      // Test: fail11.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: "
        << "Unclosed quotes in this record: "
        << record;

      // Treat an end of record as the trailing delimiter for the last field.
    } else if ( state == 1 || state == 3 ){
      Values.push_back(next);

      // Last non-blank item was a non-quoted comma.
      // So add a blank item to the vector.
    } else if ( _isVector && state ==0 ){
      if ( Values.size() > 0 ){
        Values.push_back("");
      }
    }

    // On a legal record there is no way to reach this else.
    else{
      // No test for this error.
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: "
        << "Logic bug at final state in SimpleConfigRecord parser: "
        << record;
    }

    // For a scalar record, make sure that there was exactly 1 item.
    if ( !_isVector ){
      if( Values.size() != 1 ){
        // Test: fail13.conf, fail19.conf
        throw cet::exception("SimpleConfig")
          << "SimpleConfigRecord: "
          << "Scalar type record has more than one value: "
          << record;
      }
    }

  }

  /**
   *
   * Check that the type of the current record matches the specified type.
   *
   */
  void SimpleConfigRecord::CheckType( string s) const{
    if ( Type != s  ){
      // Tests: fail14.conf,
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: "
        << "Type mismatch: "
        << "Requested: " << s
        << "  Actual:  " << Type;
    }
  }

  /**
   *
   * Check that the type of the current record is one of the vector types.
   *
   */
  void SimpleConfigRecord::AnyVector() const{
    if ( Type.substr(0,7) != "vector<"  ){

      // Test: fail20.conf
      throw cet::exception("SimpleConfig")
        << "SimpleConfigRecord: "
        << "Requested a vector for a non-vector record: "
        << record;
    }
  }

  // Return value of string as int, or throw.
  // ?? Catch errors like: 1234ccc ?
  // ?? Should I call it an error if I truncate a float to integer ?
  int SimpleConfigRecord::toInt ( const string& s ) const {

    // Check for presence for a decimal point.
    // Allow numbers with no digits after the decimal point.
    // This also throws out numbers like: 1.000 which maybe I should allow.
    string::size_type dot = s.find(".");
    if ( dot != string::npos ){

      // Only OK if it is the last character in the string.
      if ( dot != s.size()-1 ){
        // Test: fail16.conf, pass01.conf
        throw cet::exception("SimpleConfig")
          << "SimpleConfigRecord: "
          << "Floating point value for int parameter: "
          << record;
      }
    }

    // Convert to integer
    istringstream is(s);
    int i;
    is >> i;
    if ( is ) return i;

    // Test: fail17.conf
    throw cet::exception("SimpleConfig")
      << "SimpleConfigRecord: "
      << "Died converting value to int: "
      << s
      << " In record: "
      << record;
  }

  // Return value of string as double, or throw.
  // ?? Catch errors like: 1.234ccc ?
  double SimpleConfigRecord::toDouble ( const string& s ) const {
    istringstream is(s);
    double d;
    is >> d;
    if ( is ) return d;
    // Test: fail15.conf
    throw cet::exception("SimpleConfig")
      << "SimpleConfigRecord: "
      << "Died converting value to double: "
      << s
      << "In record: "
      << record;
  }

  // Return value of string as bool, or throw.
  bool SimpleConfigRecord::toBool ( const string& s ) const {

    // Convert to lower case.
    string tmp(s);
    transform ( tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

    // Look for true/false, t/f, or 1/0.
    if ( tmp == "true" || tmp == "t" || tmp == "1"){
      return  true;
    } else if ( tmp == "false" || tmp == "f" || tmp == "0" ){
      return false;
    }

    // All others are errors.
    // Test: fail18.conf
    throw cet::exception("SimpleConfig")
      << "SimpleConfigRecord: "
      << "Died converting value to bool: "
      << s
      << "In record: "
      << record;
  }

  void SimpleConfigRecord::KnownType() const {

    // Define known types.
    static vector<string> types;
    if ( types.size() == 0 ){
      types.push_back("string");
      types.push_back("int");
      types.push_back("double");
      types.push_back("bool");
      types.push_back("vector<string>");
      types.push_back("vector<int>");
      types.push_back("vector<double>");
    }

    for ( vector<string>::size_type i=0;
          i<types.size(); ++i ){
      if ( Type == types[i] ) return;
    }

    // Test: fail02.conf
    throw cet::exception("SimpleConfig")
      << "SimpleConfigRecord: "
      << "Unrecognized data type in record: "
      << record;

  }

  int SimpleConfigRecord::size() const{
    if ( _isVector ) return Values.size();
    return 1;
  } // end size()

} // end namespace mu2e
