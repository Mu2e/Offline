//
#ifndef ROOT_TRint
#include "TRint.h"
#endif
// #include "RootUtils/Utils/TCdfRoot.hh"
#include "mod/TAnaRint.hh"

TAnaRint* TAnaRint::fgInstance = 0;
TRint*    TAnaRint::fgRint     = 0;

ClassImp(TAnaRint)
//______________________________________________________________________________
TAnaRint::TAnaRint() {
  //  TCdfRoot::Instance();
}

//______________________________________________________________________________
TAnaRint::~TAnaRint() {
}


//------------------------------------------------------------------------------
TAnaRint* TAnaRint::Instance(int argc, char** argv) {
  // initialize ROOT interpreter (Rint). Filter out of the parameter
  // list `-p' option and all the files with `.tcl' extension

  static TAnaRint::Cleaner cleaner;

				// for the sake of simplicity assume that the 
				// number of parameters is less than 100
  //  char* argvv[100];
  char* argvv[100];
  int   argcc;

  if  (!  fgInstance) {
    // need to make sure that ROOT is initialized

    //    TCdfRoot::Instance();
				// loop over the parameters and filter out
				// all the .tcl files and also -p flag
    argcc = 0;
    if (argc > 0) {
      argvv[0]  = argv[0];
      for (int i=1; i<argc; i++) {
	if (strcmp(argv[i],"-p") == 0) {
	}
	else if (strstr(argv[i],".tcl") != 0) {
	}
	else {
	  argvv[argcc++] = argv[i];
	}
      }
    }

    fgInstance = new TAnaRint();
    fgRint     = new TRint("Rint@MU2E",&argcc,argvv,NULL,0);
  }
  return fgInstance;
}


//------------------------------------------------------------------------------
void TAnaRint::Delete() {
  if (fgInstance) {
    // ***********************************    delete fgInstance;
    fgInstance = 0;
  }
}

//------------------------------------------------------------------------------
TAnaRint::Cleaner::Cleaner() {
}

//------------------------------------------------------------------------------
TAnaRint::Cleaner::~Cleaner() {
  TAnaRint::Delete();
}


