//
// Parse and validate the argument list for the artProductSizes executable.
//
// See the implementation of usage() for the documentation.
//

#include "Arguments.hh"
#include <sstream>

using namespace std;


mu2e::Arguments::Arguments( int argc, char**argv ):
  minimumFractionDefault(0.05),
  minimumFraction(minimumFractionDefault),
  fileNames(){

  int i=1;

  while ( i< argc ){
    string a(argv[i]);

    // Properly formed minimumFraction option
    if ( a.find("-f=") == 0 ) {
      string sub = a.substr(3);
      minimumFraction=toDouble(sub.c_str(), 0., 1.);
    }

    // Malformed minimum fraction option.
    else if ( a.find("-f") == 0 ) {
      usage();
    }

    // Request for usage message.
    else if ( a.find("-?") == 0) {
      usage();
    }

    // Request for usage message.
    else if ( a.find("--help") == 0 ){
      usage();
    }

    // Interpret all other arguments as file names.
    else{
      fileNames.push_back( a );
    }
    ++i;
  }

  if ( fileNames.empty() ) usage();
}

double mu2e::Arguments::toDouble( char const* s, double low, double high){

  // Convert the string to a double.
  istringstream b(s);
  double f;
  b >> f;

  // Conversion to a double failed.
  if ( b.fail() ) usage();

  // Bounds check.
  if ( f < low || f > high ) usage();

  return f;
}

void  mu2e::Arguments::usage() const{
  cerr << "usage: artProductSizes [-f=<f>] filename1 [filename2] [filename3] ...\n\n"
       << " filenameN - are the names of art format ROOT event-data files to be processed.\n"
       << "             At least one filename must be specified.\n"
       << " <f>       - is a floating point number on the range [0,1]\n"
       << "             If a TTree occupies a fraction on disk of the total space in the file\n"
       << "             that is less than <f>, then a detailed analysis of its branches will not be done.\n"
       << "             The default value is: " << minimumFractionDefault
       << endl;
  exit(-1);
}
