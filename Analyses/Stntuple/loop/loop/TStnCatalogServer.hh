#ifndef TStnCatalogServer_hh
#define TStnCatalogServer_hh

#include "TNamed.h"
#include "TUrl.h"
#include "TString.h"
class TChain;
class TStnDataset;

class TStnCatalogServer: public TUrl {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TString     fKerberosUser;            // whatever it is
  TString     fRshCommand;              // "rsh -N -l " default, allow for ssh
  Int_t       fType;                    // 0:text,1:oracle,2:root,3:mysql
  Int_t       fPrintLevel;
  TString     fCatalogServer;	        // name of the STNTUPLE catalog server
  TString     fGetIPCommand;
  TString     fGetDataServerNameCommand;
  TString     fGetNEventsCommand;
  TString     fGetListOfFilesCommand;
  TString     fGetListOfBadFilesCommand;
  TString     fGetListOfFilesetsCommand;
  TString     fGetListOfDatasetsCommand;
  TString     fGetKeyCommand;
  TString     fOracleServer;		// name of the oracle server
  TString     fCafName;		// name of the caf where this is running
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStnCatalogServer(const char* Server = "txt://nowhere.fnal.gov/cdf/nodir",
		    const char* Rsh    = "rsh -N",
		    int         Print  = 0);

  virtual ~TStnCatalogServer();

  //  Int_t Init(const char* Url);

  const char*  GetDCacheDoor (const char* HeadNode);

  Int_t GetDCacheFileName(const char* HeadNode,
				  const char* FilesetPath,
				  const char* Fileset       ,
				  const char* Filename      ,
				  char*       DCacheFileName);

  virtual int InitChain     (TChain*     Chain  ,
			     const char* Book   ,
			     const char* Dataset, 
			     const char* Fileset,
			     const char* File   ,
			     Int_t       MinRun ,
			     Int_t       MaxRun );

  virtual Int_t  InitDataset(TStnDataset* Dataset            ,
			     const char*  Book    = ""       ,
			     const char*  Name    = ""       ,
			     const char*  Fileset = ""       ,
			     const char*  File    = ""       ,
			     Int_t        MinRun  = 1        ,
			     Int_t        MaxRun  = 100000000,
			     const char*  Type    = "STNTUPLE");

  virtual Int_t  AddFiles    (TChain*     Chain  ,
			      const char* Book   ,
			      const char* Dataset, 
			      const char* Fileset,
			      Int_t       Run1   ,
			      Int_t       Run2   ) = 0;
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  virtual int    FindDataset (const char* Book, const char* Dataset) = 0;

  virtual Int_t  GetNEvents  (const char* Book, 
			      const char* Dataset, 
			      const char* Fileset = 0,
			      const char* File    = 0) = 0;
  
  virtual Int_t  GetNFiles   (const char* Book       , 
			      const char* Dataset    ,
			      const char* Fileset = 0) = 0;
  
  virtual Int_t  GetNFilesets(const char* Book, const char* Dataset) = 0;

  virtual Int_t  GetRemoteServer(const char* Book,
				 const char* Dataset, 
				 const char* Fileset,
				 char*       Server,
				 char*       RemoteDir) = 0; 
//-----------------------------------------------------------------------------
// commands
//-----------------------------------------------------------------------------
  virtual const char* GetDataServerNameCommand() {
    return fGetDataServerNameCommand.Data() ;  
  }

  virtual const char* GetListOfFilesCommand() { 
    return fGetListOfFilesCommand.Data(); 
  }

  virtual const char* GetListOfBadFilesCommand() { 
    return fGetListOfBadFilesCommand.Data(); 
  }

  virtual const char* GetListOfFilesetsCommand() { 
    return fGetListOfFilesetsCommand.Data(); 
  }

  virtual const char* GetListOfDatasetsCommand() { 
    return fGetListOfDatasetsCommand.Data(); 
  }

  virtual const char* GetIPCommand     () { return fGetIPCommand.Data     (); }
  virtual const char* GetNEventsCommand() { return fGetNEventsCommand.Data(); }
  virtual const char* GetKeyCommand    () { return fGetKeyCommand.Data    (); }
  virtual const char* GetOracleServer  () { return fOracleServer.Data     (); }
  virtual const char* GetCatalogServer () { return fCatalogServer.Data    (); }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void  SetPrintLevel(Int_t Level) { fPrintLevel = Level; }

  ClassDef(TStnCatalogServer,0)
};

#endif
