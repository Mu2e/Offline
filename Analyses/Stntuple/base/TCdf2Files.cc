#include "base/TCdf2Files.hh"


#include "TUrl.h"

ClassImp(TCdf2Files)

//_____________________________________________________________________________
TCdf2Files::TCdf2Files() {
  Clear();
}

//_____________________________________________________________________________
TCdf2Files::~TCdf2Files() {
}

//_____________________________________________________________________________
void TCdf2Files::Clear(Option_t* Opt) {
  fLUM_SUM_ONLINE = 0.0;
  fLUM_SUM_OFFLINE = 0.0;
}

//_____________________________________________________________________________
void TCdf2Files::Print(Option_t* Opt) const {
  // assume that file number is printed...

  float lum_online, lum_offline;

  if ((strstr(Opt,"") == 0) || strstr(Opt,"banner") != 0) {
    printf("-----------------------------------------------------------------");
    printf("--------------------------------------------\n");
    printf(" numb filename           fileset     size       nev ");
    printf(" lo_run    lo_evt  hi_run    hi_evt     L(onl)    L(offl)\n");
    printf("-----------------------------------------------------------------");
    printf("--------------------------------------------\n");
  }

  if (strstr(Opt,"data") != 0) {

    TUrl nm(fFILE_NAME.Data());

    if (strstr(Opt,"full") == 0) printf(" %17s",fFILE_NAME.Data());
    else                         printf(" %17s",nm.GetUrl ());


    lum_online = fLUM_SUM_ONLINE;
    if (fLUM_SUM_ONLINE > 1.e6) lum_online  = 0.;

    lum_offline = fLUM_SUM_OFFLINE;
    if (fLUM_SUM_OFFLINE > 1.e6) lum_offline = 0.;

    printf(" %8s %9.3f %8i %7i %9i %7i %9i %10.3f %10.3f %2d\n",
	   fFILESET_NAME.Data(),
	   fFILE_SIZE,
	   fEVENT_COUNT,
	   fLOW_RUN,
           fLOW_EVENT,
	   fHIGH_RUN,
	   fHIGH_EVENT,
	   lum_online,
	   lum_offline,
	   fSTATUS);
  }
}
