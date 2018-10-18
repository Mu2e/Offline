//
// Main program for a utility that prints the number of events/subruns/runs
// that are present in an art format event-data file.
//

#include "Print/eventCount/Arguments.hh"
#include "Print/eventCount/FileInfo.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main( int argc, char** argv ){

  // Parse and validate command line arguments.
  mu2e::Arguments arg(argc, argv);

  // Separates the output for multiple files.
  string separator = "\n" + string(60,'=') + "\n";

  std::vector<mu2e::FileInfo> infos;
  infos.reserve( arg.fileNames.size() );

  int level = 1;
  if ( arg.style == mu2e::Arguments::minimal || 
       arg.style == mu2e::Arguments::full ) level = 0;

  for ( size_t i=0; i<arg.fileNames.size(); ++i ){
    auto const& filename = arg.fileNames.at(i);

    //if ( i!=0 && full ) cout << separator << endl;

    infos.emplace_back( filename, level );
    auto const& info(infos.back());

    if ( arg.style == mu2e::Arguments::minimal ){
      info.minimalPrint(cout);
    } else if ( arg.style == mu2e::Arguments::full ) {
      info.fullPrint(cout);
    } else if ( arg.style == mu2e::Arguments::events ) {
      info.eventPrint(cout);
    } else if ( arg.style == mu2e::Arguments::subruns ) {
      info.subrunPrint(cout);
    } else if ( arg.style == mu2e::Arguments::sam ) {
      info.samPrint(cout);
    }

  }

  return 0;
}
