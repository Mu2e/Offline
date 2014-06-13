#include "base/TCdf2Filesets.hh"


ClassImp(TCdf2Filesets)

//_____________________________________________________________________________
TCdf2Filesets::TCdf2Filesets() {
}

//_____________________________________________________________________________
TCdf2Filesets::~TCdf2Filesets() {
}

//_____________________________________________________________________________
void TCdf2Filesets::Clear(Option_t* Opt) {
}

//_____________________________________________________________________________
void TCdf2Filesets::Print(Option_t* Opt) const {

  if (strstr(Opt,"banner") != 0) {
  }

  if ((strstr(Opt,"data") != 0) || (strcmp(Opt,"") == 0)) {
    printf(" %10s %10i %10s %10s %5i %5i \n",
	   fFILESET_NAME.Data(),
	   fCREATE_TIME,
	   fDS_NAME_ID.Data(),
	   fTAPE_LABEL.Data(),
	   fFILE_COUNT,
	   fTAPE_PARTITION);
  }
}

