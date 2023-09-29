#include "Offline/GeneralUtilities/inc/CsvReader.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string.hpp>

using namespace mu2e;
using namespace std;


CsvReader::CsvReader(const string& fileName,
                     bool allowComments,
                     bool checkColumnCount,
                     bool stripOuterQuotes,
                     bool stripOuterWhitespace) :
  _allowComments(allowComments),
  _checkColumnCount(checkColumnCount),
  _stripOuterQuotes(stripOuterQuotes),
  _stripOuterWhitespace(stripOuterWhitespace),
  _ncol(0),
  _nRow(0)
{
  _fstream.open(fileName);
  if (!_fstream.is_open()) {
    throw cet::exception("CVSREADER_OPEN_FAILED")
      << "CvsReader failed to open " << fileName << "\n";
  }
  return;
}


bool CsvReader::getRow(StringVec& row) {

  row.clear();
  string line;

  bool more = true;
  while(more) {
    if(!std::getline(_fstream,line)) {
      _fstream.close();
      return false;
    }
    if(_allowComments) {
      string temp = line;
      boost::trim(temp);
      if( temp.size() > 0 && line[0] != '#') more = false;
    } else {
      more = false;
    }
  }

  _nRow++;

  std::size_t j;
  j = 0;  // current position in parse
  bool quote = false;     // deal with special char in quotes
  string value; // builds the current column in this row
  value.reserve(100);

  more = true;
  while (more) {  // while not the end of the row string
    if (j == line.size()) {  // j past end of line
      if (quote) {  // if in a quote, then a newline was part of quote
        // read the next line and append, so we can finish the row
        string temp;
        if(!std::getline(_fstream,temp)) {
          throw cet::exception("CVSREADER_OPEN_QUOTE")
            << "CvsReader could not complete multi-line quote\n";
        }
        line.append(1,'\n');
        line.append(temp);

      } else {  // non-quoted, non-comment newline - end of row
        row.push_back(value);
        value.clear();
        j++;
        more = false;
      }
    } else if (line[j] == '"') {
      if (!quote) {
        quote = true; // this char starts a quote
        value.append(1,line[j]);
        j++;
      } else {
        // in quote, check if embedded quote
        if ( j>0 && line[j - 1] == '\\') {  // has form \"
          // if stripping outer quotes, remove embedded quote format
          if(_stripOuterQuotes) value.pop_back(); // remove the backslash
          value.append(1,line[j]);  // add the "
          j++;       // just continue to next, slash already processed
        } else if (j < line.size() - 1 && line[j + 1] == '"') {  // has form ""
          value.append(1,line[j]);
          // if not stripping outer quotes, leave embedded quote format
          if(!_stripOuterQuotes) value.append(1,line[j]);
          j = j + 2;  // move past second quote
        } else {      // must be end quote
          quote = false;
          value.append(1,line[j]);
          j++;
        }
      }
    } else if (line[j] == ',') {  // commas separate the columns
      if (quote) {    // if the comma was in a quote, then add it to value
        value.append(1,line[j]);
        j++;
      } else {  // non-quoted comma, make new column
        row.push_back(value);
        value.clear();
        j++;
      } // end if in quote

    } else {  // just a char in a column
      value.append(1,line[j]);
      j++;
    }

  } // end while more loop

  if ( _checkColumnCount ) {
    if( _ncol == 0 ) _ncol = row.size();
    if(row.size() != _ncol) {
      throw cet::exception("CVSREADER_COLUMN_COUNT")
        << "CvsReader row had wrong column count\n";
    }
  }

  if (_stripOuterQuotes) {
    for(auto& col : row) {
      string temp = col;
      boost::trim(temp);
      if(temp.size()>1 && temp[0]=='"' && temp[temp.size()-1]=='"' ) {
        col = string( temp.c_str() + 1, temp.size()-2 );
      }
    }
  }

  if(_stripOuterWhitespace) {
    for(auto& col : row) {
      boost::trim(col);
    }
  }

  return true;
}
