//
// Parse and validate the argument list for the artProductSizes executable.
//
// See the implementation of usage() for the documentation.
//

#include "Print/artProductSizes/Arguments.hh"
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

    else if ( a.find("-?") == 0) {
      usage();
    }

    else if ( a.find("--help") == 0 ){
      usage();
    }

    // Reject anything else that starts with a -
    else if ( a.find("-") == 0 ){
      cerr << "Unrecognized command line option: " << a << endl;
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
       << "             that is more than <f>, then a detailed description of its branches will be printed.\n"
       << "             The default value is: " << minimumFractionDefault << "\n"
       << " -?        - print this message and exit\n"
       << " --help    - print this message and exit\n\n"
       << "For each named file, the output has two sections.\n"
       << "  - The total size on disk; the same as given by ls -l\n"
       << "  - Information about all of the top level objects in the file.\n"
       << "  - For each top level object that is a tree and is above the size threshold\n"
       << "    print information about all of its branches\n\n"
       << "Details of the information printed for top level objects\n"
       << "  - If a top level object is a TTree, the following information is printed\n"
       << "     - The size on disk of the TTree, in bytes\n"
       << "     - The size on disk per entry of the TTree, in bytes\n"
       << "     - The number of entries in the TTree\n"
       << "     - The size on disk of the TTree, as a fraction of the file size\n"
       << "     - The name of the TTree\n"
       << "  - If a top level object is a TKey, the information printed is the same as for\n"
       << "    a TTree with the execption that the size per entry and the number of entries\n"
       << "    are printed as --\n"
       << "  - If a top level object is of any other type, an error message is printed\n"
       << "  - The Total line is obtained by adding up the lines above it\n"
       << "  - The size on disk in the Total is less than the number of bytes in the disk file\n"
       << "    and the sum of fractions may be less than 100%.  This is because the ROOT header\n"
       << "    and some of the other ROOT bookkeeping is not counted in the sum of the top level objects.\n\n"
       << "Details of the information printed for TBranches in top level TTrees\n"
       << "  - The size on disk of the TBranch, in bytes\n"
       << "  - The size on disk per entry of the TBranch, in bytes\n"
       << "  - The size on disk of the TBranch, as a fraction of the size on disk of the TTree\n"
       << "  - The name of the TBranch\n"
       << endl;
  exit(-1);
}
