#ifndef GeneralUtilities_CsvReader_hh
#define GeneralUtilities_CsvReader_hh


//  g++ -O0 -g -std=c++17 -o csv csv.cc
// Read a csv file by rows and return a vector of strings for columns
//   StringVec row;
//   CsvReader cr("file.csv" <options>)
//   while(cr.getRow(row)) {
//     //proces row
//   }
// of strings containing the colums for a row
// and the separation between strings in the output vector
// This code implements standard
//    https://datatracker.ietf.org/doc/html/rfc4180
// with a few extensions:
// - whitespace lines are ignored (optional)
// - lines beginning with "#" are comments and ignored (optional)
// - escaped double quotes in a string can be \" as well as ""
// - the standard does not prescribe how to process embedded quotes
//   We leave them as escape sequences if the option is to return
//   the qoutes around the column, and to transform them into a single
//   quote if the choice is to remove the double quotes arounf a column.
//
// To demo the algorithm, this file
/*
# comment
 0, 1, 2
 col0, col 1, col2
,,

# comment " stuff'
"word0, ,word1", "abc\" " , " don't ""def"""

start,"break
here
and here", follow
*/
// with defaults, will return 5 rows (using "|" for a delimiter)
/*
|0|1|2|
|col0|col 1|col2|
||||
|word0, ,word1|abc"|don't "def"|
|start|break
here
and here|follow|
*/
//

#include "Offline/GeneralUtilities/inc/StringVec.hh"
#include <iostream>
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>

namespace mu2e {

class CsvReader {

public:
  CsvReader(const std::string& fileName,
            // allow blank and comment lines which have # as the first char
            bool allowComments = true,
            // check that all rows have the same number of columns
            bool checkColumnCount = true,
            // columns |"abc "," ""OK"""| return |abc | "OK"|
            bool stripOuterQuotes = true,
            // columns | "abc ", x  z  | return |"abc "|x  z|
            bool stripOuterWhitespace = true  );

  // caller provides StringVec to fill,
  // return is false when there is no row left to return
  bool getRow(StringVec& row);
  // the number of valid rows returned to current point
  size_t nRow() { return _nRow; }

private:

  bool _allowComments;
  bool _checkColumnCount;
  bool _stripOuterQuotes;
  bool _stripOuterWhitespace;

  std::ifstream _fstream;
  size_t _ncol;
  size_t _nRow;
};

}

#endif
