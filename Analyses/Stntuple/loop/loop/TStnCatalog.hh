#ifndef TStnCatalog_hh
#define TStnCatalog_hh

#include "TNamed.h"
#include "TString.h"
class TChain;
class TStnDataset;
class TStnCatalogServer;

class TStnCatalog: public TNamed {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TObjArray*  fListOfCatalogServers;
  TString*    fUser;		// user name (get it via klist)
  Int_t       fPrintLevel;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStnCatalog(const char* Name  = "StntupleCatalog");
  virtual ~TStnCatalog();

  Int_t Init(const char* Url, const char* Rsh = "rsh -N", int Print = 0);
				// dataset = "book:dataset:fileset:file"

  Int_t          InitDataset(TStnDataset* Dataset);

  Int_t          InitDataset(TStnDataset* Dataset            ,
			     const char*  Book               ,
			     const char*  Name               ,
			     const char*  Fileset = ""       ,
			     const char*  File    = ""       ,
			     Int_t        Run1    = 0        ,
			     Int_t        Run2    = 100000000,
			     const char*  Type    = "STNTUPLE");

  Int_t InitChain(TChain*     Chain              , 
		  const char* Dataset            ,
		  Int_t       Run1   = 0         , 
		  Int_t       Run2   = 100000000 );
  
  Int_t InitChain(TChain*     Chain              ,
		  const char* Book               ,
		  const char* Dataset            ,
		  const char* Fileset = 0        ,
		  const char* File    = 0        ,
		  Int_t       Run1    = 0        ,
		  Int_t       Run2    = 100000000);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t  GetPrintLevel() { return fPrintLevel; }

  TStnCatalogServer* GetCatalogServer(const char* Book, const char* Dataset);

				// conveniences

  Int_t  GetNEvents  (const char* Book       , 
		      const char* Dataset    , 
		      const char* Fileset = 0,
		      const char* File    = 0,
		      const char* Server  = 0) ;

  Int_t  GetNFiles   (const char* Book      ,
		      const char* Dataset   , 
		      const char* Fileset = 0,
		      const char* Server  = 0) ;

  Int_t  GetNFilesets(const char* Book       , 
		      const char* Dataset    , 
		      const char* Server  = 0) ;

  void   Print       (const char* Book               ,
		      const char* Dataset            ,
		      const char* Fileset = 0        ,
		      const char* File    = 0        ,
		      Int_t       Run1    = 0        , 
		      Int_t       Run2    = 100000000,
		      const char* Server  = 0        ) const ;
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void  SetPrintLevel(Int_t Level);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Clear(Option_t* Opt = "");
  void  Print(Option_t* Opt = "") const ;

  ClassDef(TStnCatalog,0)
};

#endif
