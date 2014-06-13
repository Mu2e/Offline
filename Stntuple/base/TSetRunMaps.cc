#include "base/TSetRunMaps.hh"



ClassImp(TSetRunMaps)

//_____________________________________________________________________________
TSetRunMaps::TSetRunMaps() {
}

//_____________________________________________________________________________
TSetRunMaps::~TSetRunMaps() {
}

//_____________________________________________________________________________
void TSetRunMaps::Print(const char* Opt) const {
  if ((strstr(Opt,"") == 0) || strstr(Opt,"banner") != 0) {
    printf("-----------------------------------------------------------------");
    printf("--------------------------------\n");
    printf("   process_name      CID   JOBSET   create_date create_user TSID ");
    printf(" update_date update_user  create_host update_host   finished\n");
    printf("-----------------------------------------------------------------");
    printf("--------------------------------\n");
  }

  if (strstr(Opt,"data") != 0) {
    printf("%20s %8i %8i %8i %12s %8i %12s %12s %12s %2i\n",
	   fPROCESS_NAME.Data(),
	   fCID,
	   fJOBSET,
	   fCREATE_DATE,
	   fCREATE_USER.Data(),
	   fUPDATE_DATE,
	   fUPDATE_USER.Data(),
	   fCREATE_HOST.Data(),
	   fUPDATE_HOST.Data(),
	   fFINISHED);
  }

}
