//
// Find number of runs, subruns and events in an art format ROOT file.
//
// See the implementation of usage() for the documentation.
//

#include "ROOTtools/eventCount/FileInfo.hh"

#include "TError.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

namespace {

  // Helper function to fill in the info for one tree;
  void treeInfo ( std::string const& treeName,   // input argument
                  std::string const& filename,   // input argument
                  TFile *            file,       // input argument
                  bool&              exists,     // output argument
                  unsigned long&     nEntries    // output argument
                  ){

    TTree * tree;
    file->GetObject( treeName.c_str(), tree);
    if ( !tree ) return;

    exists   = true;
    nEntries = tree->GetEntries();
  }
}

mu2e::FileInfo::FileInfo( std::string const& name ):filename(name){

  // Suppress warning messages about "no dictionary" and error messages about "file does not exist"
  // This is a little dangerous since it might suppress other warnings too ...
  int errorSave     = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kBreak;
  TFile* file = new TFile(filename.c_str());
  gErrorIgnoreLevel = errorSave;

  if ( file->IsZombie() )return;
  openable = true;

  treeInfo ( "Events",  filename, file, hasEvents,  events);
  treeInfo ( "SubRuns", filename, file, hasSubRuns, subRuns  );
  treeInfo ( "Runs",    filename, file, hasRuns,    runs  );

}

void mu2e::FileInfo::fullPrint( ostream& ) const{
  if ( openable ){
    std::cout << "Filename:   " << filename << std::endl;
  }else{
    std::cout << "Filename:   " << filename << " could not be opened."
              << std::endl;
    std::cout << "             check that the file exists, is a root file and has read permission."
              << std::endl;
  }
  if ( hasRuns ){
    std::cout << "   Runs:    " << runs << std::endl;
  }else{
    std::cout << "   Runs tree is missing or not readable" << std::endl;
  }
  if ( hasSubRuns ){
    std::cout << "   SubRuns: " << subRuns<< std::endl;
  }else{
    std::cout << "   SubRuns tree is missing or not readable" << std::endl;
  }
  if ( hasEvents ){
    std::cout << "   Events:  " << events << std::endl;
  }else{
    std::cout << "   Events tree is missing or not readable" << std::endl;
  }

}

void mu2e::FileInfo::minimalPrint( ostream& ) const{
  std::cout << filename << " "
            << status() << " "
            << runs     << " "
            << subRuns  << " "
            << events
            << endl;
}
