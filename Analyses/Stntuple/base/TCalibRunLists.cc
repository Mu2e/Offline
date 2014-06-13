///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "base/TCalibRunLists.hh"


ClassImp(TCalibRunLists)

//_____________________________________________________________________________
TCalibRunLists::TCalibRunLists() {
}

//_____________________________________________________________________________
TCalibRunLists::~TCalibRunLists() {
}

//_____________________________________________________________________________
void TCalibRunLists::Print(const char* Opt) const {
  if ((strstr(Opt,"") == 0) || strstr(Opt,"banner") != 0) {
    printf("-----------------------------------------------------------------");
    printf("--------------------------------\n");
    printf("    CID   TABLE   STATUS   RUN   VERSION   DATE   reate_user ");
    printf(" TSID  update_date update_user  create_host update_host  DAQ_STATUS");
    printf(" ALGORITHM   PERIOD");
    printf("-----------------------------------------------------------------");
    printf("--------------------------------\n");
  }

  if (strstr(Opt,"data") != 0) {
    printf("%7i %20s %10s %8i %3i %8i %12s %8i %12s %12s %12s %12s %12s %3i\n",
	   fCID,
	   fCALIB_TABLE.Data(),
	   fDATA_STATUS.Data(),
	   fCALIB_RUN,
	   fCALIB_VERSION,
	   fCREATE_DATE,
	   fCREATE_USER.Data(),
	   fUPDATE_DATE,
	   fUPDATE_USER.Data(),
	   fCREATE_HOST.Data(),
	   fUPDATE_HOST.Data(),
	   fDAQ_STATUS.Data(),
	   fALGORITHM.Data(),
	   fPERIOD);
  }

}
