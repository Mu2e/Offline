//
// Main program for a utility that looks at an art format event-data file
// and reports on on which objects use how much disk space.
//
// Fixme:
// 1) Add an option to print sizes in kB, MB, GB
// 2) Add an option to choose "sensible units" for each number.
// 3) Add options to print size of branches of named trees, regardless of size"
//      -r -s -e --run --subrun --event --tree=name
//
#include "Print/artProductSizes/Arguments.hh"
#include "Print/artProductSizes/RootSizeOnDisk.hh"

#include "TError.h"
#include "TFile.h"

#include <iostream>
#include <string>

using namespace std;

int main( int argc, char** argv ){

  // Parse and validate command line arguments.
  mu2e::Arguments arg(argc, argv);

  // Separates the output for multiple files.
  string separator = "\n" + string(60,'=') + "\n";

  int n(-1);
  for ( auto const& filename : arg.fileNames ) {
    ++n;

    if ( arg.fileNames.size() > 1 && n!=0 ) cout << separator << endl;

    // Suppress warnings messages about "no dictionary".
    // This is a little dangerous since it might suppress other warnings too ...
    int errorSave     = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    TFile* file = new TFile(filename.c_str());
    gErrorIgnoreLevel = errorSave;

    // Extract and print the information.
    mu2e::RootSizeOnDisk info( filename, file );
    file->Close();
    info.print(cout,arg.minimumFraction);

    // Mark end of output for this file.
    if ( arg.fileNames.size() > 1 ) cout << "Done: " << filename << endl;

  }

  return 0;
}
