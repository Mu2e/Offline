//
// An example of how run-dependent
// conditions sets can be accessed in a stand-alone bin
//
// To build this, copy it to an apropriate repo src dir
// You can see what you need in your SConscript by looking at
// the build for dbTool in Offline/DbService/src
//
// arguments are
// PURPOSE: the conditions set purpose
// VERSION: the conditions set version
// RUN: the run number to lookup
// TABLE: the table to dump
//
// example:
// dbBin TRACKER_VST v1_2 100000 TrkStrawStatusLong
//

#include <string>
#include <iostream>
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/DbService/inc/DbEngine.hh"

using namespace mu2e;
using namespace std;

int main(int argc, char**argv) {

  if(argc!=5) { // first arg is exe name
    cout << "Four required arguments: PURPOSE VERSION RUN TABLE" << endl;
    return 1;
  }

  string purpose(argv[1]);
  string version(argv[2]);
  uint32_t run = stoul(string(argv[3]));
  uint32_t subrun = 0;
  string table(argv[4]);

cout << "Using "<<purpose << " " << version << " "
       << run << " " << table<< endl;

DbEngine engine;

engine.setVerbose(2);

DbIdList idList; // read file of db connection details
engine.setDbId( idList.getDbId("mu2e_conditions_prd") );

engine.setVersion( DbVersion(purpose,version) );

// add overriding text file
//auto coll = DbUtil::readFile("myfile.txt");
//_engine.addOverride(coll);

int tid =  engine.tidByName(table);

DbLiveTable liveTable = engine.update(tid, run, subrun);

cout <<liveTable.table().csv() << endl;

  return 0;
}
