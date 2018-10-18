//
// Find number of runs, subruns and events in an art format ROOT file.
//
// See the implementation of usage() for the documentation.
//

#include "Print/eventCount/FileInfo.hh"
#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "canvas/Persistency/Provenance/SubRunAuxiliary.h"

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

using namespace std;

mu2e::FileInfo::FileInfo( std::string const& name, int level ):filename(name){

  // Suppress warning messages about "no dictionary" and error messages about "file does not exist"
  // This is a little dangerous since it might suppress other warnings too ...
  int errorSave     = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kBreak;
  TFile* file = new TFile(filename.c_str());
  gErrorIgnoreLevel = errorSave;

  if ( file->IsZombie() )return;
  openable = true;

  treeInfo ( "Events",  file, hasEvents,  events );
  treeInfo ( "SubRuns", file, hasSubRuns, subRuns);
  treeInfo ( "Runs",    file, hasRuns,    runs   );
  if(level>0) makeVectors( file );
  file->Close();

}

void mu2e::FileInfo::fullPrint( ostream& os) const{
  if ( openable ){
    os << "Filename:   " << filename << std::endl;
  }else{
    os << "Filename:   " << filename << " could not be opened."
              << std::endl;
    os << "             check that the file exists, is a root file and has read permission."
              << std::endl;
  }
  if ( hasRuns ){
    os << "   Runs:    " << runs << std::endl;
  }else{
    os << "   Runs tree is missing or not readable" << std::endl;
  }
  if ( hasSubRuns ){
    os << "   SubRuns: " << subRuns<< std::endl;
  }else{
    os << "   SubRuns tree is missing or not readable" << std::endl;
  }
  if ( hasEvents ){
    os << "   Events:  " << events << std::endl;
  }else{
    os << "   Events tree is missing or not readable" << std::endl;
  }

}

void mu2e::FileInfo::minimalPrint( ostream& os) const{
  os << filename << " "
     << status() << " "
     << runs     << " "
     << subRuns  << " "
     << events
     << endl;
}

void mu2e::FileInfo::eventPrint( ostream& os) const{
  for ( auto e : eventList ) {
    os 
      << setw(9) << e.run() 
      << setw(9) << e.subRun() 
      << setw(9) << e.event()
      << endl;
  }
}

void mu2e::FileInfo::subrunPrint( ostream& os) const{
  for ( auto s : subrunList ) {
    os 
      << setw(9) << s.run()
      << setw(9) << s.subRun() 
      << endl;
  }
}

void mu2e::FileInfo::samPrint( ostream& os) const{

  // order the events
  //std::sort(eventList.begin(),eventList.end());
  art::EventID e_min = *(std::min_element(eventList.begin(),eventList.end()));
  art::EventID e_max = *(std::max_element(eventList.begin(),eventList.end()));

  //order the subruns
  //std::sort(subrunList.begin(),subrunList.end());
  art::SubRunID s_min = *(std::min_element(subrunList.begin(),subrunList.end()));
  art::SubRunID s_max = *(std::max_element(subrunList.begin(),subrunList.end()));

  os << "{" << std::endl;
  os << "  \"event_count\"      : " << events   << "," << std::endl;

  os << "  \"dh.first_run_subrun\" : " << s_min.run() << "," << std::endl;
  os << "  \"dh.first_subrun\"     : " << s_min.subRun() << "," << std::endl;

  os << "  \"dh.first_run_event\"  : " << e_min.run() << "," << std::endl;
  os << "  \"dh.first_subrun_event\" : " << e_min.subRun() << "," << std::endl;
  os << "  \"dh.first_event\"      : " << e_min.event()   << "," << std::endl;

  os << "  \"dh.last_run_subrun\"  : " << s_max.run() << "," << std::endl;
  os << "  \"dh.last_subrun\"      : " << s_max.subRun()   << "," << std::endl;

  os << "  \"dh.last_run_event\"   : " << e_max.run() << "," << std::endl;
  os << "  \"dh.last_subrun_event\"  : " << e_max.subRun() << "," << std::endl;
  os << "  \"dh.last_event\"       : " << e_max.event()   << "," << std::endl;

  os << "  \"runs\"             : [\n";
  for (auto sr : subrunList ){
    os << "\t[";
    os << sr.run() << "," << sr.subRun() << ",\"unknown\"]";
    if(sr!=subrunList.back()) os << ",\n" ;
  }
  os << "\n\t\t]" << std::endl;
  os << "}" << std::endl;

}


void mu2e::FileInfo::treeInfo ( std::string const& treeName,
				TFile* file,
				bool&              exists,
				unsigned long&     nEntries
				){
  exists = false;
  nEntries = 0;
  TTree * tree;
  file->GetObject( treeName.c_str(), tree);
  if ( !tree ) return;
  
  exists   = true;
  nEntries = tree->GetEntries();
}

void mu2e::FileInfo::makeVectors(TFile* file) {

  eventList.clear();
  subrunList.clear();

  TTree* tr = nullptr;
  file->GetObject( "Events", tr);
  art::EventAuxiliary* ea = new art::EventAuxiliary();
  tr->SetBranchAddress("EventAuxiliary", &ea);
  for (int i=0; i<tr->GetEntries(); ++i) {
    tr->GetEntry(i);
    eventList.emplace_back(ea->id());
  }

  tr = nullptr;
  file->GetObject( "SubRuns", tr);
  art::SubRunAuxiliary* sa = new art::SubRunAuxiliary();
  tr->SetBranchAddress("SubRunAuxiliary", &sa);
  for (int i=0; i<tr->GetEntries(); ++i) {
    tr->GetEntry(i);
    subrunList.emplace_back(sa->id());
  }
  
}
