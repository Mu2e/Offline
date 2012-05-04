// ConsolidateMultipleTrackFiles.cpp

#include "FigureOfMerit/inc/ConsolidateMultipleTrackFiles.h"

// Root includes.
#include "TFile.h"
#include "TTree.h"

// C++ includes.
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>

// mu2e includes

namespace mu2e {

ConsolidateMultipleTrackFiles::ConsolidateMultipleTrackFiles()
  : OUTPUT_fileNames(true)
{}


void ConsolidateMultipleTrackFiles::consume() {
  rawCount = 0;
  processedCount = 0;
  std::vector<std::string> fnames = obtainFileNames();
  establishRootFile();
  for (unsigned int i = 0; i < fnames.size(); ++i) {
    processOneFile(fnames[i]);
  }
  completeRootFile();
  std::cout << "File contains " << processedCount 
            << " tracks, out of " << rawCount << " ("
            << 100.0*processedCount/rawCount << "%)\n";
} // consume

std::vector<std::string> ConsolidateMultipleTrackFiles::obtainFileNames() 
{
  std::vector<std::string> fileNames;
  std::cout << "1 - Explicit entry of several file names \n"
            << "2 - Name of a file containing a list of file names \n"
            << "3 - Directory starting string, number range, and file name \n"
            << "Enter file name input mode: "; 
  bool validInputStyleChoice = false;
  int inputStyle;
  std::string name;
  while (!validInputStyleChoice) {
   std::cin >> inputStyle;
    switch (inputStyle) {
      case 1:
        validInputStyleChoice = true;
        fileNames = explicitFileNames();
        break;
      case 2: 
        validInputStyleChoice = true;
        fileNames = listFileNames();
        break;
      case 3: 
        validInputStyleChoice = true;
        fileNames = directoriesFileNames();
        break;
      default:
        std::cout << "The input mode should be 1, 2, or 3\n";
    } // end of switch
  } // end of input style choice while-loop
  if (OUTPUT_fileNames) {
    std::cout << "List of files: \n";
    for (unsigned int j = 0; j < fileNames.size(); ++j) {
      std::cout << fileNames[j] << "\n";
    }
  }  
  return fileNames;
} // obtainFileNames

std::vector<std::string> 
ConsolidateMultipleTrackFiles::explicitFileNames() const
{
  std::vector<std::string> fileNames;
  std::cout << "Enter one file name per line. \n";
  std::cout << "Terminate the list with an entry reading end\n";
  std::string name;
  for (int i = 1; i < 100000; ++i) {
    std::cout << "File " << i << ": ";
    std::cin >> name;
    if (name == "end") break;
    { ifstream test(name.c_str());
      if (!test) {
        std::cout << "Cannot open file " << name.c_str() << "\n";
        continue;
      }
      fileNames.push_back(name);
    } 
  } // end of for loop getting names
  return fileNames;
} // explicitFileNames

std::vector<std::string> ConsolidateMultipleTrackFiles::listFileNames() const
{
  std::vector<std::string> fileNames;
  std::cout << "Enter name of a file. ";
  std::cout << "The file should contain one root file name per line\n";
  std::string listName;
  std::cin >> listName;
  ifstream list(listName.c_str());
  if (!list) {
      std::cout << "Cannot open file " << listName << "\n";
      std::exit(1);
  }
  std::string name;
  for (int i = 1; i < 100000; ++i) {
    std::cout << "File " << i << ": ";
    list >> name;
    if ( (!list) || (name == "end") ) {
      std::cout << "List of files has ended\n"; 
      break;
    }
    { ifstream test(name.c_str());
      if (!test) {
        std::cout << "Cannot open file " << name.c_str() << "\n";
        continue;
      }
      fileNames.push_back(name);
    } 
  } // end of for loop getting names
  return fileNames;
} //listFileNames


std::vector<std::string> 
ConsolidateMultipleTrackFiles::directoriesFileNames() const
{
  std::vector<std::string> fileNames;
  std::string startingString;
  std::cout << "Enter directory starting string \n"
            << "(include trailing underscore if applicable): ";
  std::cin >> startingString;
  int nLow  = 999999;
  int nHigh = 999999;
  std::cout << "Enter first number to be appended: ";
  std::cin >> nLow;
  std::cout << "Enter last number to be appended:  ";
  std::cin >> nHigh;
  if ( nLow > nHigh ) {
    std::cout << "Invalid range\n";
    std::exit(1);
  }
  std::string filename;
  std::cout << "Enter file name \n"
            << "(one file of this name should be in each directory): ";
  std::cin >> filename;
  for (int n = nLow; n <= nHigh; ++n) {
    std::ostringstream nameStream;
    nameStream << startingString << n << "/" << filename;
    std::cout << "---- will look for file " << nameStream.str() << "\n";
    { ifstream test(nameStream.str().c_str());
      if (!test) {
        std::cout << "Cannot open file " << nameStream.str() << "\n";
        continue;
      }
      fileNames.push_back(nameStream.str());
    } 
  } // end of loop over directory range n
  return fileNames;
} // directoriesFileNames

void ConsolidateMultipleTrackFiles::establishRootFile() 
{
  std::string outputFileName;
  std::cout << "Enter desired output file name: ";
  std::cin >> outputFileName;
  outputFile = new TFile (outputFileName.c_str(), "NEW");
  if (outputFile->IsZombie()) {
    std::cerr << "Failed to open the output file processedTracks.root\n";
    std::exit(1);
  }
  outputTree = new TTree ("trkdiag","properties of tracks with fitmom>0");
  outputTree->SetBranchStyle(0);

  outputTree->Branch("fitstatus", &fitstatus, "fitstatus/I");    
  outputTree->Branch("t0err"    , &t0err,     "t0err/F"    );    
  outputTree->Branch("t0"       , &t0,        "t0/F"       );      
  outputTree->Branch("nhits"    , &nhits,     "nhits/I"    );   
  outputTree->Branch("ndof"     , &ndof,      "ndof/I"     );   
  outputTree->Branch("nactive"  , &nactive,   "nactive/I"  ); 
  outputTree->Branch("chisq"    , &chisq,     "chisq/F"    );   
  outputTree->Branch("fitcon"   , &fitcon,    "fitcon/F"   );  
  outputTree->Branch("fitmom"   , &fitmom,    "fitmom/F"   );  
  outputTree->Branch("fitmomerr", &fitmomerr, "fitmomerr/F");
  outputTree->Branch("mct0"     , &mct0,      "mct0/F"     );   
  outputTree->Branch("mcmom"    , &mcmom,     "mcmom/F"    );   
  outputTree->Branch("diowt"    , &diowt,     "diowt/F"    );   

#ifdef FOLOWS_EXAMPLE_BUT_FAILS
  outputTree->Branch("fitstatus/I", &fitstatus);    
  outputTree->Branch("t0err/F"    , &t0err);    
  outputTree->Branch("t0/F"       , &t0);    
  outputTree->Branch("nhits/I"    , &nhits);    
  outputTree->Branch("ndof/I"     , &ndof);    
  outputTree->Branch("nactive/I"  , &nactive);
  outputTree->Branch("chisq/F"    , &chisq);    
  outputTree->Branch("fitcon/F"   , &fitcon);    
  outputTree->Branch("fitmom/F"   , &fitmom);    
  outputTree->Branch("fitmomerr/F", &fitmomerr);    
  outputTree->Branch("mct0/F"     , &mct0);    
  outputTree->Branch("mcmom/F"    , &mcmom);    
#endif

#ifdef WRONG_I_THINK
  outputTree->Branch("fitstatus:I", &fitstatus);    
  outputTree->Branch("t0err:F"    , &t0err);    
  outputTree->Branch("t0:F"       , &t0);    
  outputTree->Branch("nhits:I"    , &nhits);    
  outputTree->Branch("ndof:I"     , &ndof);    
  outputTree->Branch("nactive:I"  , &nactive);
  outputTree->Branch("chisq:F"    , &chisq);    
  outputTree->Branch("fitcon:F"   , &fitcon);    
  outputTree->Branch("fitmom:F"   , &fitmom);    
  outputTree->Branch("fitmomerr:F", &fitmomerr);    
  outputTree->Branch("mct0:F"     , &mct0);    
  outputTree->Branch("mcmom:F"    , &mcmom);    
#endif

  fitstatus_br = outputTree->GetBranch("fitstatus");
  t0err_br     = outputTree->GetBranch("t0err");
  t0_br        = outputTree->GetBranch("t0");
  nhits_br     = outputTree->GetBranch("nhits");
  ndof_br      = outputTree->GetBranch("ndof");
  nactive_br   = outputTree->GetBranch("nactive");
  chisq_br     = outputTree->GetBranch("chisq");
  fitcon_br    = outputTree->GetBranch("fitcon");
  fitmom_br    = outputTree->GetBranch("fitmom");
  fitmomerr_br = outputTree->GetBranch("fitmomerr");
  mct0_br      = outputTree->GetBranch("mct0");
  mcmom_br     = outputTree->GetBranch("mcmom");
  diowt_br     = outputTree->GetBranch("diowt");

}  // establishRootFile


void  ConsolidateMultipleTrackFiles::processOneFile(std::string const & fname) 
{
  if (OUTPUT_fileNames) {
    std::cout << "Processing file " << fname << "\n";
  }
  // Open the file
  TFile* inputFile = new TFile(fname.c_str(), "READ");
  if (inputFile->IsZombie()) {
    std::cerr << "Failed to open the input file " << fname << "\n";
    std::exit(1);
  }
  
  // get the data by hook or crook, into trkdiagDataRow
  TTree* inputTree = (TTree*)inputFile->Get("ReadKalFits/trkdiag");
  if (inputTree == 0) {
    std::cerr << "Failed to open tracks at ReadKalFits/trkdiag "
              << "for the input file " << fname << "\n";
    std::exit(1);
  }
  inputTree->SetBranchStyle(0);
  inputTree->SetBranchAddress("fitstatus", &fitstatus);    
  inputTree->SetBranchAddress("t0err"    , &t0err);    
  inputTree->SetBranchAddress("t0"       , &t0);    
  inputTree->SetBranchAddress("nhits"    , &nhits);    
  inputTree->SetBranchAddress("ndof"     , &ndof);    
  inputTree->SetBranchAddress("nactive"  , &nactive);
  inputTree->SetBranchAddress("chisq"    , &chisq);    
  inputTree->SetBranchAddress("fitcon"   , &fitcon);    
  inputTree->SetBranchAddress("fitmom"   , &fitmom);    
  inputTree->SetBranchAddress("fitmomerr", &fitmomerr);    
  inputTree->SetBranchAddress("mct0"     , &mct0);    
  inputTree->SetBranchAddress("mcmom"    , &mcmom);    
  inputTree->SetBranchAddress("diowt"    , &diowt);    

  long entriesCount = inputTree->GetEntries();
  
  for (long i=0; i<entriesCount; ++i) { 
    inputTree->GetEntry(i);
    if (fitstatus > 0) {
      outputTree->Fill();
      ++processedCount;
    }
  }
  rawCount += entriesCount;

} 

void ConsolidateMultipleTrackFiles::completeRootFile() 
{
  outputFile->Write();
  outputFile->Close();
  delete outputFile;
  outputFile = 0;
}  // completeRootFile


} // end namespace mu2e

int main()
{
  mu2e::ConsolidateMultipleTrackFiles cmtf;
  cmtf.consume();
  return 0;
}
