//
// Parse and validate the argument list for the artProductSizes executable.
//
// See the implementation of usage() for the documentation.
//

#include "Print/eventCount/Arguments.hh"
#include <iostream>

using namespace std;


mu2e::Arguments::Arguments( int argc, char**argv ):
  fileNames(){

  int i=1;

  while ( i< argc ){
    string a(argv[i]);

    if ( a.find("-?") == 0) {
      usage();
    }

    else if ( a.find("--help") == 0 ){
      usage();
    }

    else if ( a.find("--full") == 0 ){
      style = full;
    }

    else if ( a.find("--minimal") == 0 ){
      style = minimal;
    }

    else if ( a.find("--events") == 0 ){
      style = events;
    }

    else if ( a.find("--subruns") == 0 ){
      style = subruns;
    }

    else if ( a.find("--sam") == 0 ){
      style = sam;
    }

    // Reject any other argument that starts with a -
    else if ( a.find("-") == 0 ){
      cerr << "Unrecognized argument:: " << a << endl;
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

void  mu2e::Arguments::usage() const{

  cerr << "usage: eventCount [--full] [--minimal] filename1 [filename2] [filename3] ...\n\n"
       << " filenameN - are the names of art format ROOT event-data files to be processed.\n"
       << "             At least one filename must be specified.\n"
       << " -?        - print this message and exit\n"
       << " --help    - print this message and exit\n"
       << " --minimal - minimal printout [default]\n"
       << " --full    - full printout\n"
       << " --events  - print run/subrun/event for all events\n"
       << " --subruns - print run/subrun for all subruns, including those with no events \n"
       << " --sam     - print json format for sam record \n"
       << "\nIf both --minimal and --full are present, the last one wins.\n\n"
       << "Output format for option minimal: \n"
       << "   filename OK/BAD  nRuns nSubRuns nEvents\n\n"
       << "   BAD means that the file could not be opened or that one or more TTrees is missing.\n"
       << "   To understand which part is bad, look at the full output.\n\n"
       << "Output format for the full format is self describing.\n"
       << endl;
  exit(-1);
}
