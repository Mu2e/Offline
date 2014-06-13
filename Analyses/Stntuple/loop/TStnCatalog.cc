//-----------------------------------------------------------------------------
//  Dec 28 2000 P.Murat: base class for STNTUPLE input module
//  TopDir syntax: host:directory , i.e. fcdfsgi2.fnal.gov:/cdf/data/cafdfc ,
//  network catalog not implemented yet
//-----------------------------------------------------------------------------
#include "TChain.h"
#include "TChainElement.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "Stntuple/loop/TTxtCatalogServer.hh"
#include "Stntuple/loop/THttpCatalogServer.hh"
#include "Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/base/TStnDataset.hh"

ClassImp(TStnCatalog)

//_____________________________________________________________________________
TStnCatalog::TStnCatalog(const char* name) :
  TNamed     (name,name),
  fPrintLevel(0)
{
  char        url [1000];
  const char *env;
  TString     cmd;
  TString     rsh("rsh -n -N"); // Default is set here
//-----------------------------------------------------------------------------
// modify the default Rsh command?
//-----------------------------------------------------------------------------
  cmd  = "fgrep -s 'Stntuple.Catalog.Rsh'  $HOME/.rootrc $PWD/.rootrc ";
  cmd += " | sed 's/#.*//' | awk '{if(NF>1) print $0}'"; //remove comments
  cmd += " | cut -d' ' -f2-300 | uniq";
  FILE* f = gSystem->OpenPipe(cmd.Data(),"r");
  while (fscanf(f,"%s",url) != EOF) {
    rsh = TString(url);
    rsh.ReplaceAll(TString(":"),TString(" ")); // TString sucks (hack it)
    //printf(" New rsh command: %s -> %s\n",url,rsh.Data());
  }
  gSystem->ClosePipe(f);

//-----------------------------------------------------------------------------
// browse through .rootrc files to find server
//-----------------------------------------------------------------------------
  cmd  = "fgrep -s 'Stntuple.Catalog ' $HOME/.rootrc  $PWD/.rootrc ";
  cmd += " | sed 's/#.*//' | awk '{if(NF>1) print $0}'"; //remove comments
  cmd += " | awk '{print $2}' | uniq";
  f    = gSystem->OpenPipe(cmd.Data(),"r");
  fListOfCatalogServers = new TObjArray();
  TStnCatalogServer* server;
  while (fscanf(f,"%s",url) != EOF) {
    server  = (TStnCatalogServer*) fListOfCatalogServers->FindObject(url);
    if (! server) {

      TUrl u(url);

      if (strcmp(u.GetProtocol(),"txt") == 0) {
	fListOfCatalogServers->Add(new TTxtCatalogServer(url,rsh.Data()));
      }
      else if (strcmp(u.GetProtocol(),"http") == 0) {
	fListOfCatalogServers->Add(new THttpCatalogServer(url,rsh.Data()));
      }
    }
  }
  gSystem->ClosePipe(f);
//-----------------------------------------------------------------------------
// in addition to that see if STNTUPLE_CATALOG variable is defined in the 
// environment - it may add one more server
//-----------------------------------------------------------------------------
  env  = gSystem->Getenv("STNTUPLE_CATALOG");
  if (env) {
    server  = (TStnCatalogServer*) fListOfCatalogServers->FindObject(env);
    if (! server) {

      TUrl u(env);

      if (strcmp(u.GetProtocol(),"txt") == 0) {
	fListOfCatalogServers->Add(new TTxtCatalogServer(env));
      }
      else if (strcmp(u.GetProtocol(),"http") == 0) {
	fListOfCatalogServers->Add(new THttpCatalogServer(env));
      }
    }
  }
}

//_____________________________________________________________________________
TStnCatalog::~TStnCatalog() {
  // destructor: module owns its histograms

  fListOfCatalogServers->Delete();
  delete fListOfCatalogServers;
}

//_____________________________________________________________________________
void TStnCatalog::SetPrintLevel(Int_t Level) {
  fPrintLevel = Level;
  TStnCatalogServer* server;
  TIter it(fListOfCatalogServers);
  while (server = (TStnCatalogServer*) it.Next()) {
    server->SetPrintLevel(Level);
  }
  return;
}


//_____________________________________________________________________________
int TStnCatalog::Init(const char* TopUrl, const char* Rsh, int Print) {
  // TopUrl - the topmost directory of the data catalog
  TStnCatalogServer* server;

  server = (TStnCatalogServer*) fListOfCatalogServers->FindObject(TopUrl);

  TUrl url(TopUrl);

  if (! server) {
    if (strcmp(url.GetProtocol(),"txt") == 0) {
      fListOfCatalogServers->Add(new TTxtCatalogServer(TopUrl,Rsh,Print));
    }
    else if (strcmp(url.GetProtocol(),"http") == 0) {
      fListOfCatalogServers->Add(new THttpCatalogServer(TopUrl,Rsh,Print));
    }
  }
  else
    Warning("init","Server already exists. Rsh command not updated.");
  
  return 0;
}


//_____________________________________________________________________________
int TStnCatalog::InitChain(TChain*     Chain, 
			   const char* Book,
			   const char* Dataset, 
			   const char* Fileset,
			   const char* File,
			   Int_t       Run1,
			   Int_t       Run2) {
//-----------------------------------------------------------------------------
// find server and pass the job to it
//-----------------------------------------------------------------------------
  char   fn[200];

  int                rc = -1;

  if (strcmp(Book,"") == 0) {
					// non-cataloged local file
    Chain->AddFile(File,TChain::kBigNumber);
    return 0;
  }
  else if (strcmp(Book,"root") == 0) {
					// non-cataloged remote file
    sprintf(fn,"root:%s",File);
    Chain->AddFile(fn,TChain::kBigNumber);
    return 0;
  }

  TStnCatalogServer* s  = GetCatalogServer(Book,Dataset);
  if (s) rc = s->InitChain(Chain,Book,Dataset,Fileset,File,Run1,Run2);
  return rc;
}


//_____________________________________________________________________________
int TStnCatalog::InitChain(TChain* Chain, const char* Dataset, 
			   Int_t Run1, Int_t Run2) 
{
  // Dataset = "book:dataset[:fileset[:file]]"
  // parse dataset name and pass the rest to another routine

  int  i, npos;
  char book[1000], dataset[1000], fileset[1000], file[1000];

  int  first = 0;

  TString s(Dataset);
  
  book   [0] = 0;
  dataset[0] = 0;
  fileset[0] = 0;
  file   [0] = 0;

  if (s[0] == '/') {
//-----------------------------------------------------------------------------
//  dataset name starts from '/' : local file
//-----------------------------------------------------------------------------
    strcpy(file,Dataset);
  }
  else if ((s.Index("root:") == 0)) {
//-----------------------------------------------------------------------------
//  there is "root:" in the name - either local or remote file
//-----------------------------------------------------------------------------
    if ((s[5] == '/') && (s[6] == '/')) {
				// "root://" - remote file
      strcpy(file,Dataset);
    }
    else {
      strcpy(file,Dataset+5);
    }
  }
  else {
//-----------------------------------------------------------------------------
//  dataset
//-----------------------------------------------------------------------------
    i = s.Index(":",first);
    if (i > 0) {
      npos = i;
      strncpy(book,s.Data(),npos);
      book[i] = 0;
				// figure dataset
      first = i+1;
      i = s.Index(":",first);
      if (i > 0) {
	npos = i-first;
	strncpy(dataset,s.Data()+first,npos);
	dataset[npos] = 0;
				// figure fileset
	first = i+1;
	i = s.Index(":",first);
	if (i > 0) {
	  npos = i-first;
	  strncpy(fileset,s.Data()+first,npos);
	  fileset[npos] = 0;
                                // figure file
	  first = i+1;
	  npos = s.Length()-first;
	  strncpy(file,s.Data()+first,npos);
	  file[npos] = 0;
	}
	else {
	  strcpy(fileset,s.Data()+first);
	}
      }
      else {
	strcpy(dataset,s.Data()+first);
      }
    }
  }

  return InitChain(Chain,book,dataset,fileset,file,Run1,Run2);
}

//_____________________________________________________________________________
int TStnCatalog::InitDataset(TStnDataset* Dataset) {

  TStnCatalogServer* s;

  s = GetCatalogServer(Dataset->GetBook(),Dataset->GetName());

  int rc = -1;
  if (s) {
    rc = s->InitDataset(Dataset);
  }
  return rc;
}

//_____________________________________________________________________________
int TStnCatalog::InitDataset(TStnDataset* Dataset            ,
			     const char*  Book               ,
			     const char*  Name               ,
			     const char*  Fileset ,
			     const char*  File    ,
			     Int_t        MinRun  ,
			     Int_t        MaxRun  ,
			     const char*  Type    )
{
  int                rc = -1;
  TStnCatalogServer* s  = GetCatalogServer(Book,Name);

  if (s) {
    rc = s->InitDataset(Dataset,Book,Name,Fileset,File,MinRun,MaxRun,Type);
  }
  else {
    Error("InitDataset",Form("Cant find book=%s dataset=%s\n",Book,Name));
  }
  return rc;
}

//_____________________________________________________________________________
TStnCatalogServer* TStnCatalog::GetCatalogServer(const char* Book   , 
						 const char* Dataset) {
  // for Book="file" return the 1st available

  TStnCatalogServer* server;
  TIter it(fListOfCatalogServers);

  while (server = (TStnCatalogServer*) it.Next()) {
    if ((strcmp(Book,"file") == 0) || (strcmp(Book,"dir") == 0))
      break;
    if (server->FindDataset(Book,Dataset) != 0)
      break;
  }

  return server;
}

//_____________________________________________________________________________
Int_t TStnCatalog::GetNEvents(const char* Book   ,
			      const char* Dataset,
			      const char* Fileset,
			      const char* File   ,
			      const char* ServerName) {
  int nev = -1;
  TStnCatalogServer* server;

  if (ServerName == 0)
    server = GetCatalogServer(Book,Dataset);
  else
    server = (TStnCatalogServer*) fListOfCatalogServers->FindObject(ServerName);

  if (server)
    nev = server->GetNEvents(Book,Dataset,Fileset,File);

  return nev;
}

//_____________________________________________________________________________
Int_t TStnCatalog::GetNFiles(const char* Book   ,
			     const char* Dataset,
			     const char* Fileset,
			     const char* ServerName)
{
  int nfiles = -1;
  TStnCatalogServer* server;

  if (ServerName == 0)
    server = GetCatalogServer(Book,Dataset);
  else
    server = (TStnCatalogServer*) fListOfCatalogServers->FindObject(ServerName);

  if (server)
    nfiles= server->GetNFiles(Book,Dataset,Fileset);

  return nfiles;
}


//_____________________________________________________________________________
Int_t TStnCatalog::GetNFilesets(const char* Book   ,
				const char* Dataset,
				const char* ServerName)
{
  int nfilesets = -1;
  TStnCatalogServer* server;

  if (ServerName == 0)
    server = GetCatalogServer(Book,Dataset);
  else
    server = (TStnCatalogServer*) fListOfCatalogServers->FindObject(ServerName);

  if (server)
    nfilesets= server->GetNFilesets(Book,Dataset);

  return nfilesets;
}

//_____________________________________________________________________________
void TStnCatalog::Clear(Option_t* Opt)
{
}

//_____________________________________________________________________________
void  TStnCatalog::Print(const char* Book   ,
			 const char* Dataset,
			 const char* Fileset,
			 const char* File   ,
			 Int_t       Run1   , 
			 Int_t       Run2   ,
			 const char* Server ) const
{
}

//_____________________________________________________________________________
void TStnCatalog::Print(Option_t* Opt) const
{
}
