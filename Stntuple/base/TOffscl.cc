#include "base/TOffscl.hh"


ClassImp(TOffscl)

//_____________________________________________________________________________
TOffscl::TOffscl() {
}

//_____________________________________________________________________________
TOffscl::~TOffscl() {
}

//_____________________________________________________________________________
void TOffscl::Print(const char* Opt) const {
  if ((strstr(Opt,"") == 0) || strstr(Opt,"banner") != 0) {
    printf("---------------------------------------------------------------\n");
    printf("   COMPNAME CID   COMPONENT  VALUE COMPNAME ");
    printf("---------------------------------------------------------------\n");
  }

  if (strstr(Opt,"data") != 0) {
    printf("%20s %8i %8i %12.5f\n",
	   fCOMPNAME.Data(),
	   fCID,
	   fCOMPONENT,
	   fCALIB);
  }

}
