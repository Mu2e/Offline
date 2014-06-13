#ifndef TTxtCatalogServer_hh
#define TTxtCatalogServer_hh

#include "TNamed.h"
#include "TUrl.h"
#include "TString.h"
class TChain;
class TStnDataset;

#include "TStnCatalogServer.hh"

class TTxtCatalogServer: public TStnCatalogServer {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTxtCatalogServer(const char* Server = "txt://nowhere.fnal.gov/cdf/nodir",
		    const char* Rsh    = "rsh -N",
		    int         Print  = 0);

  virtual ~TTxtCatalogServer();

  int InitChain(TChain*     Chain  ,
		const char* Book   ,
		const char* Dataset, 
		const char* Fileset,
		const char* File   ,
		Int_t       MinRun ,
		Int_t       MaxRun );

  Int_t          InitDataset(TStnDataset* Dataset            ,
			     const char*  Book    = ""       ,
			     const char*  Name    = ""       ,
			     const char*  Fileset = ""       ,
			     const char*  File    = ""       ,
			     Int_t        MinRun  = 1        ,
			     Int_t        MaxRun  = 100000000,
			     const char*  Type    = "STNTUPLE");
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  virtual int    FindDataset (const char* Book, const char* Dataset);

  virtual Int_t  GetNEvents  (const char* Book, 
			      const char* Dataset, 
			      const char* Fileset = 0,
			      const char* File    = 0);
  
  virtual Int_t  GetNFiles   (const char* Book       , 
			      const char* Dataset    ,
			      const char* Fileset = 0);

  virtual Int_t  GetNFilesets(const char* Book, const char* Dataset);

  virtual Int_t  GetRemoteServer(const char* Book,
				 const char* Dataset, 
				 const char* Fileset,
				 char*       Server,
				 char*       RemoteDir);
//-----------------------------------------------------------------------------
// commands
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  virtual Int_t  AddFiles    (TChain*     Chain  ,
			      const char* Book   ,
			      const char* Dataset, 
			      const char* Fileset,
			      Int_t       Run1   ,
			      Int_t       Run2   );
  
  ClassDef(TTxtCatalogServer,0)
};

#endif
