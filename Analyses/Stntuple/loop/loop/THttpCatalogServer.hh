#ifndef THttpCatalogServer_hh
#define THttpCatalogServer_hh

#include "TNamed.h"
#include "TUrl.h"
#include "TString.h"
#include "TObjArray.h"
class TChain;
class TStnDataset;

#include "TStnCatalogServer.hh"

class THttpCatalogServer: public TStnCatalogServer {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  THttpCatalogServer(const char* Server = "html://nowhere.fnal.gov/cdf/nodir",
		     const char* Rsh    = "rsh -N",
		     int         Print  = 0);

  virtual ~THttpCatalogServer();

  int            InitChain(TChain*     Chain  ,
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

  int            InitListOfFilesets(TStnDataset* Dataset, 
				    const char*  Fileset,
				    const char*  File,
				    int          RMin,
				    int          RMax,
				    TObjArray*   ListOfFilesets,
				    TObjArray*   ListOfFiles);

  int            InitListOfBadFiles(TStnDataset* Dataset);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  virtual int    FindDataset (const char* Book, const char* Dataset);

  // read the AAA_FILEs.html file, which replaces sam db calls
  int            LoadSamDBFromHtml(const char* Book, const char* Dataset);

  int            GetDatasetKey(const char* Book, 
			       const char* Dataset,
			       const char* Key,
			       char*       Buffer);

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

  // contents of AAA_FILES.html
  TObjArray fAAAFilesHtml;
  
  ClassDef(THttpCatalogServer,0)
};

#endif
