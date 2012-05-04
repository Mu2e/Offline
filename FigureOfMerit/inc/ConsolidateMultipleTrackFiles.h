#ifndef CONSOLIDATEMULTIPLETRACKFILES_H
#define CONSOLIDATEMULTIPLETRACKFILES_H

// ----------------------------------------------------------------------
//
// ConsolidateMultipleTrackFiles.h
//
// Tool to use ConsolidateMultipleTrackFiles class for reduction of DIO, 
// RPC or CE file collections preliminary to applicatoin of FMtool.
//
// ----------------------------------------------------------------------

#include <vector>
#include <string>

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTree.h"

namespace mu2e {

class ConsolidateMultipleTrackFiles
{
public:
  ConsolidateMultipleTrackFiles();
  void consume();
  std::vector<std::string> obtainFileNames();  
  void establishRootFile();
  void processOneFile(std::string const & filename);  
  void completeRootFile();
  
public:
  int rawCount;
  int processedCount;

private:

  std::vector<std::string> explicitFileNames()    const;
  std::vector<std::string> listFileNames()        const;
  std::vector<std::string> directoriesFileNames() const;

private:
  bool OUTPUT_fileNames;

private:
  TFile*   outputFile;
  TTree*   outputTree;

private:
    int      fitstatus;    
    float    t0err;     
    float    t0;     
    int      nhits;
    int      ndof;
    int      nactive;      
    float    chisq;
    float    fitcon;
    float    fitmom;
    float    fitmomerr;    
    float    mct0;
    float    mcmom;  
    float    diowt;  

    TBranch* fitstatus_br;
    TBranch* t0err_br;
    TBranch* t0_br;
    TBranch* nhits_br;
    TBranch* ndof_br;
    TBranch* nactive_br;
    TBranch* chisq_br;
    TBranch* fitcon_br;
    TBranch* fitmom_br;
    TBranch* fitmomerr_br;
    TBranch* mct0_br;
    TBranch* mcmom_br;
    TBranch* diowt_br;

  
}; // ConsolidateMultipleTrackFiles

} // end namespace mu2e

#endif // CONSOLIDATEMULTIPLETRACKFILES_H
