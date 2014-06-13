//-----------------------------------------------------------------------------
//  Mar 09 2004 P.Murat: 
//-----------------------------------------------------------------------------
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRegexp.h"

#include "Stntuple/base/TStnFileset.hh"

ClassImp(TStnFileset)
//_____________________________________________________________________________
TStnFileset::TStnFileset(const char* Name    , 
			 const char* Server  ,
			 Int_t       Location,
			 Int_t       NEvents ,
			 Int_t       MinRunNumber  , 
			 Int_t       MaxRunNumber  ): 
  TNamed(Name,Name)
{
  Set(Server,Location,NEvents,MinRunNumber,MaxRunNumber);
}

//_____________________________________________________________________________
TStnFileset::~TStnFileset() {
  // it is more efficient to delete objects in the order opposite to 
  // their creation

  delete fDataServer;
}

//-----------------------------------------------------------------------------
void TStnFileset::Set(const char* Server, Int_t Location, Int_t NEvents, 
		      Int_t MinRunNumber, Int_t MaxRunNumber) {

  fLocation     = Location;    

  if (Server) fDataServer = new TUrl(Server);
  else        fDataServer = 0;

  fNEvents      = 0;
  fNFiles       = 0;
  fMinRunNumber = MinRunNumber;
  fMaxRunNumber = MaxRunNumber;
}

//-----------------------------------------------------------------------------
void TStnFileset::Clear(Option_t* Opt) {
}

//_____________________________________________________________________________
void TStnFileset::Print(Option_t* Opt) const {
  Error("Print","Not implemented yet");
}
