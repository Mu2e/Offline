//-----------------------------------------------------------------------------
//  Dec 28 2000 P.Murat: base class for STNTUPLE input module
//  TopDir syntax: host:directory , i.e. fcdfsgi2.fnal.gov:/cdf/data/cafdfc ,
//  network catalog not implemented yet
//-----------------------------------------------------------------------------
#include "TEnv.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "Stntuple/base/TCdf2Files.hh"
#include "Stntuple/base/TStnFileset.hh"
#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/loop/TTxtCatalogServer.hh"

ClassImp(TTxtCatalogServer)
//_____________________________________________________________________________
TTxtCatalogServer::TTxtCatalogServer(const char* Url, const char* Rsh,
				     int Print): 
TStnCatalogServer(Url,Rsh,Print)
{
  TString  prefix; // initializes to ""
  TString  cmd; 
  char     hostname[200];
  //-----------------------------------------------------------------------------
//  now can use $USER on caf
//-----------------------------------------------------------------------------
  if (strcmp(GetHost(),hostname) == 0) {
    prefix = Form("/%s/scripts",GetFile());
  }
  else {
    prefix = Form("%s -l %s %s /%s/scripts",fRshCommand.Data(),
		  gSystem->Getenv("USER"),GetHost(),GetFile());
  }

  fGetDataServerNameCommand = Form("%s/get_server_name"         ,prefix.Data());
  fGetNEventsCommand        = Form("%s/get_nevents"             ,prefix.Data());
  fGetListOfFilesCommand    = Form("%s/get_list_of_files"       ,prefix.Data());
  fGetListOfBadFilesCommand = Form("%s/get_list_of_bad_files"   ,prefix.Data());
  fGetListOfFilesetsCommand = Form("%s/get_list_of_filesets"    ,prefix.Data());
  fGetListOfDatasetsCommand = Form("%s/get_list_of_datasets"    ,prefix.Data());
  fGetKeyCommand            = Form("%s/get_key"                 ,prefix.Data());

  // What are we going to use as retrieval commands
  if (fPrintLevel > 0)
    printf(" DEBUG: Retrieving data with: %s\n",prefix.Data());
}

//_____________________________________________________________________________
TTxtCatalogServer::~TTxtCatalogServer() {
  // destructor: module owns its histograms
}


//_____________________________________________________________________________
int TTxtCatalogServer::AddFiles(TChain*     Chain  , 
				const char* Book   ,
				const char* Dataset, 
				const char* Fileset,
				Int_t       Run1   ,
				Int_t       Run2   ) 
{
  char cmd[200], fs[200], fn[200];
  char size[50], date[50], time[50];
  int  nevents,  rlow, elow, rhigh, ehigh;
//-----------------------------------------------------------------------------
// now get list of files
//-----------------------------------------------------------------------------
  sprintf(cmd,"%s -b %s -d %s -s %s -r %i:%i",
	  GetListOfFilesCommand(),
	  Book,Dataset,Fileset,Run1,Run2);

  FILE* f = gSystem->OpenPipe(cmd,"r");
  while (fscanf(f,"%s %s %s %s %s %i %i %i %i %i",
		fs,fn,size,date,time,&nevents,&rlow,&elow,&rhigh,&ehigh) > 0) {
//-----------------------------------------------------------------------------
//  time to add file - make sure we're not adding the same file 2nd time
//-----------------------------------------------------------------------------
    TObjArray* list_of_files = Chain->GetListOfFiles();

    TObjArrayIter it(list_of_files);
    TChainElement* found;

    while ((found = (TChainElement*) it.Next())) {
      if (strcmp(fn,found->GetTitle()) == 0) break;
    }
    
    if (! found) Chain->AddFile(fn,nevents);
  }
  gSystem->ClosePipe(f);

  return 0;
}


//_____________________________________________________________________________
int TTxtCatalogServer::FindDataset(const char* Book, const char* Dataset) {

  int found = 0;
  char cmd[2000], dataset[2000];

  sprintf(cmd,"%s -b %s",GetListOfDatasetsCommand(),Book);
  if (fPrintLevel > 0) printf("GetListOfDatasets cmd= %s\n",cmd);

  FILE* pipe = gSystem->OpenPipe(cmd,"r");
  while ((! found) && (fscanf(pipe,"%s",dataset) != EOF)) {
    if (fPrintLevel > 0) printf("found dataset %s\n",dataset);
    if (strcmp(dataset,Dataset) == 0) found = 1;
  }
  gSystem->ClosePipe(pipe);

  return found;
}


//_____________________________________________________________________________
int TTxtCatalogServer::InitChain(TChain*     Chain, 
				 const char* Book,
				 const char* Dataset, 
				 const char* Fileset,
				 const char* File   ,
				 Int_t       Run1   ,
				 Int_t       Run2   ) 
{
//-----------------------------------------------------------------------------
// assume that all the files in the directory should be chained
// examples: 
// TStnCatalog::InitChain(chain,"stntuple/gqcd1g")
// TStnCatalog::InitChain(chain,"stntuple/gqcd1g","GI0783,GI0888")
// TStnCatalog::InitChain(chain,"stntuple/gqcd1g","GI0783",
//                              "stn_gqcd1g_GI0783.0.root")
//-----------------------------------------------------------------------------
  int        lorun, hirun, nev, loevt, hievt;
  float      size;
  char       buf[1000], fn[100], fs[100], date[100], time[100];
  char       remote_server[200], remote_dir[200], remote_file[200];
  char*      line;
  TObjArray  list_of_filesets;
  TString    cmd;
  TObjArray  files   (100);
  TObjString *ostr;
//-----------------------------------------------------------------------------
// first treat special cases
//-----------------------------------------------------------------------------

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

  if ((Fileset == 0) || (strcmp(Fileset,"") == 0)) {
//-----------------------------------------------------------------------------
//  the whole dataset
//-----------------------------------------------------------------------------
    cmd = Form("%s -b %s -d %s -r %i:%i",GetListOfFilesCommand(),
	       Book,Dataset,Run1,Run2);

    FILE* pipe = gSystem->OpenPipe(cmd,"r");
    while (fgets(buf,10000,pipe)) { files.Add(new TObjString(buf)); }
    gSystem->ClosePipe(pipe);

    TString prefix;
    TIter itt(&files);
    while ((ostr = (TObjString*) itt.Next())) {
      line = (char*) ostr->String().Data();
      sscanf(line,"%s %s %f %s %s %i %i %i %i %i", 
	     fs,fn,&size,date,time,&nev,&lorun,&loevt,&hirun,&hievt);
//-----------------------------------------------------------------------------
// file name is already formed
//-----------------------------------------------------------------------------
      TObjArray* list_of_files = Chain->GetListOfFiles();
      TObjArrayIter it(list_of_files);
      TChainElement* found;
      while ((found = (TChainElement*) it.Next())) {
	if (strcmp(fn,found->GetTitle()) == 0) break;
      }
      if (! found)
	Chain->AddFile(fn,nev);
    }
  }
  else if ((File == 0) || (strcmp(File,"") == 0)) {
//-----------------------------------------------------------------------------
// one or several filesets, for each fileset get its server, directory and 
// list of files. All files belonging to the fileset are supposed to reside 
// in the same directory
//-----------------------------------------------------------------------------
    char*    list_of_fs = new char[strlen(Fileset)+1];
    strcpy(list_of_fs,Fileset);

    char*    fs = strtok(list_of_fs,", ");
    do {
      AddFiles(Chain,Book,Dataset,fs,Run1,Run2);
    } while ((fs = strtok(0,", ")));

    delete [] list_of_fs;
  }
  else {
//-----------------------------------------------------------------------------
// one or several file names, separators: " ," (space and comma)
// assume that all the files belong to the same Book, Dataset and Fileset
// for each file define server and remote directory and add file to the chain
//-----------------------------------------------------------------------------
    GetRemoteServer(Book,Dataset,Fileset,remote_server,remote_dir);

    char*    list_of_files = new char[strlen(File)+1];
    strcpy(list_of_files,File);

    char*    file = strtok(list_of_files,", ");
    do {
      sprintf(remote_file,"root://%s/%s/%s",remote_server,remote_dir,file);
      Chain->AddFile(remote_file,TChain::kBigNumber);
    } while ((file = strtok(0,", ")));

    delete [] list_of_files;
  }
//-----------------------------------------------------------------------------
// in principle can foresee check for the dublicates... 
// in the end print list of files
//-----------------------------------------------------------------------------
  Chain->GetListOfFiles()->Print();
  return 0;
}

//_____________________________________________________________________________
int TTxtCatalogServer::InitDataset(TStnDataset*     Dataset,
				   const char*      Book   ,
				   const char*      Name   ,
				   const char*      Fileset,
				   const char*      File   ,
				   Int_t            MinRun ,
				   Int_t            MaxRun ,
				   const char*      Type   )
{
  FILE*        pipe;
  int          time, lorun, hirun, rmin, rmax, nev, loevt, hievt, status;
  int          n_filesets, mc_flag;
  float        size;
  char         buf[10000], date[1000], ctime[2000], fn[2000], fs[1000];
  char         full_name[200], directory[200], server[200], pnfs_path[1000];
  const char   *line, *book, *dset, *dir;
  TObjArray    *list_of_filesets;
  TObjString   *ostr;
  TStnFileset  *fileset, *fset;
  TString      cmd, s_fileset, s_file;
//-----------------------------------------------------------------------------
  if (strcmp(Book,"") != 0) {
    // inititalize dataset, otherwise assume it to be already initialized
    int rc = Dataset->Init(Book,Name,MinRun,MaxRun,Type);
    if (rc < 0) return rc;
  }

  book             = Dataset->GetBook();
  dset             = Dataset->GetName();
  s_fileset        = Fileset;
  s_file           = File;
//-----------------------------------------------------------------------------
//  handle case of a single file
//-----------------------------------------------------------------------------
  if (strcmp(book,"file") == 0) {
    Dataset->AddFile(dset);
    return 0;
  }
  else if (strcmp(book,"dir") == 0) {
//-----------------------------------------------------------------------------
// Fileset is the directory name , File is the filename patters
//-----------------------------------------------------------------------------
    dir = Fileset;
    cmd = Form("ls -al %s | grep %s | awk '{print $9}'",dir,File);
    pipe = gSystem->OpenPipe(cmd,"r");

    while (fgets(buf,10000,pipe)) { 
      sscanf(buf,"%s",fn);
      sprintf(fs,"%s/%s",dir,fn);
      Dataset->AddFile(fs);
    }
    return 0;
  }


  int found = FindDataset(book,dset);
  if (! found) return -1;
//-----------------------------------------------------------------------------
//  dataset is found , initalize metadata part of the dataset, so far all we 
//  need is a list of files with the corresponding "URL-type" prefixes
//  step 1: read all of the filesets at once - one transaction is better!
//-----------------------------------------------------------------------------
  list_of_filesets = Dataset->GetListOfFilesets();
  rmin             = Dataset->GetMinRunNumber  ();
  rmax             = Dataset->GetMaxRunNumber  ();

  TObjArray filesets(100);
  TObjArray files   (100);
//-----------------------------------------------------------------------------
// retrieve MC type
//-----------------------------------------------------------------------------
  buf[0] = 0;
  cmd = Form("%s -b %s -d %s -k mc_flag",GetKeyCommand(),book,dset);
  if (fPrintLevel > 0)
    printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());
  pipe = gSystem->OpenPipe(cmd,"r");
  fgets(buf,10000,pipe);
  if (strlen(buf) > 0) {
    sscanf(buf,"%i",&mc_flag);
    Dataset->SetMcFlag(mc_flag);
  }
  gSystem->ClosePipe(pipe);
//-----------------------------------------------------------------------------
// retrieve list of bad files - hopefully short - do it only once - so far can
// do it multiple times
// kludge: allow to do reading multiple times for empty list
//-----------------------------------------------------------------------------
  if (Dataset->DoneBadFiles() == 0) {
    cmd = Form("%s -b %s -d %s", GetListOfBadFilesCommand(),book,dset);
    if (fPrintLevel > 0)
      printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());
    
    pipe = gSystem->OpenPipe(cmd,"r");
    while (fgets(buf,10000,pipe)) { 
      sscanf(buf,"%s",fn);
      Dataset->GetListOfBadFiles()->Add(new TObjString(fn));
    }
    gSystem->ClosePipe(pipe);
    Dataset->SetDoneBadFiles();
  }
//-----------------------------------------------------------------------------
// retrieve list of filesets
//-----------------------------------------------------------------------------
  if (Dataset->GetListOfFilesets()->GetEntries() == 0) {
    cmd = Form("%s -b %s -d %s -r %i:%i", 
	       GetListOfFilesetsCommand(),
	       book,dset,rmin,rmax);
  }
  else {
				// OK, assume 1 fileset... redundant...

    fset = (TStnFileset*) Dataset->GetListOfFilesets()->At(0);
    cmd = Form("%s -b %s -d %s -s %s", 
	       GetListOfFilesetsCommand(),
	       book,dset,fset->GetName());
  }

  TString grep_filesets;

  if (fPrintLevel > 0)
    printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());
  pipe = gSystem->OpenPipe(cmd,"r");
  int first = 1;
  while (fgets(buf,10000,pipe)) { 
    sscanf(buf,"%s",fs);
//-----------------------------------------------------------------------------
// if s_fileset != "" only one fileset has been requested
//-----------------------------------------------------------------------------
    if (fPrintLevel > 0)
      printf(" DEBUG: fileset: %s -> %s\n",s_fileset.Data(),fs);
    if ((s_fileset == "") || (s_fileset == fs)) {
      filesets.Add(new TObjString(buf)); 
      if (first) { 
	grep_filesets += Form("%s",fs); 
	first = 0; 
      }
      else {
	if (fGetListOfFilesetsCommand.CompareTo(fRshCommand) >= 0) {
	  grep_filesets += Form("%c",0x5c);
	}
	grep_filesets += Form("|%s",fs);
      }
    }
  }
  gSystem->ClosePipe(pipe);

  n_filesets = filesets.GetEntries();
//-----------------------------------------------------------------------------
// and also the list of the non-DCache files, only for the requested filesets
//-----------------------------------------------------------------------------
  cmd = Form("%s -b %s -d %s -s \"%s\" -r %i:%i",
	     GetListOfFilesCommand(),
	     book,dset,grep_filesets.Data(),
	     rmin,rmax);

  //printf("GREPPING: %s\n",cmd.Data());
  if (fPrintLevel > 0)
    printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());
  pipe = gSystem->OpenPipe(cmd,"r");
  while (fgets(buf,10000,pipe)) { 
//-----------------------------------------------------------------------------
// skip comment lines
//-----------------------------------------------------------------------------
    if ( buf[0] != '#') {
      sscanf(buf,"%s %s ", fs,fn);
//-----------------------------------------------------------------------------
// handle case when a single file has been requested
//-----------------------------------------------------------------------------
      if ((s_file == "") || (strstr(fn,s_file.Data()) != 0)) {
	files.Add(new TObjString(buf)); 
      }
    }
  }
  gSystem->ClosePipe(pipe);

  //  int n_files = files.GetEntriesFast();
//-----------------------------------------------------------------------------
// loop again over the fileset definition lines and parse the information
// only the files corresponding to the requested filesets are stored. 
// Loop over the filesets - some may be Oracle
//-----------------------------------------------------------------------------
  for (int i=0; i<n_filesets; i++) {
    line = (char*) ((TObjString*) filesets.At(i))->String().Data();
    sscanf(line,"%s %s %s %i %i %i",fs,server,directory,&nev,&lorun,&hirun);
//-----------------------------------------------------------------------------
//  new fileset, make sure we are not adding it twice - loop over the filesets
//-----------------------------------------------------------------------------
    if (fPrintLevel > 0)
      printf(" DEBUG: line: %s at %s in %s\n",fs,server,directory);
    TString srv = server;
//-----------------------------------------------------------------------------
// deal with different TUrl handling by ROOT v4 and ROOT v5
//-----------------------------------------------------------------------------
    if ((gROOT->GetVersionInt() > 50000) && (directory[1] != '/')) srv += "/";
    srv   += directory;

    fileset = (TStnFileset*) list_of_filesets->FindObject(fs);
    if (fileset == 0) {
      fileset = new TStnFileset(fs,srv,0,nev,lorun,hirun);
      list_of_filesets->Add(fileset);
    }
    else {
      if (fileset->GetDataServer() == 0) {
	fileset->Set(srv,0,nev,lorun,hirun);
      }
      else {
	Error("InitDataset",Form("attempt to reinitialize fileset %s",fs));
	goto NEXT_FILESET;
      }
    }
//-----------------------------------------------------------------------------
// run range is OK by selection
//-----------------------------------------------------------------------------
    if (strcmp(server,"dfc") == 0) {
//-----------------------------------------------------------------------------
// query Oracle, directory="book:fileset"
//-----------------------------------------------------------------------------
      TString ss = directory;
      TString s1 = ss(0,ss.Index(":"));
      TString s2 = ss(ss.Index(":")+1,ss.Length());

      cmd  = Form("sqlplus cdf_reader/reader@%s <<EOF | ",GetOracleServer());
      cmd += Form("grep -s %s\n",s2.Data());
      cmd += Form("set linesize 1000;\n");
      cmd += Form("select fileset_name, file_name, file_size, ");
      cmd += Form("creation_time, ");
      cmd += Form("event_count, low_run, low_event, high_run, high_event ");
      cmd += Form("from %s.cdf2_files where ",s1.Data());
      cmd += Form("fileset_name='%s' ",s2.Data());
      cmd += Form("and low_run  <= %i and high_run >= %i",rmax,rmin);
      cmd += Form(";\nEOF\n");
      
      if (fPrintLevel > 0)
	printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());      
      pipe = gSystem->OpenPipe(cmd.Data(),"r");
      while (fscanf(pipe,"%s %s %f %i %i %i %i %i %i",
  		    fs,fn,&size,&time,&nev,&lorun,&loevt,&hirun,&hievt)!=EOF) {
	if (fPrintLevel > 20)
	  printf("Result of DFC query: %s %s %f %i %i %i %i %i %i\n",
		 fs,fn,size,time,nev,lorun,loevt,hirun,hievt);
	
//-----------------------------------------------------------------------------
// so far list of bad files is kept locally ...
// make sure we're not adding a bad file
//-----------------------------------------------------------------------------
	if (Dataset->GetListOfBadFiles()->FindObject(fn) != 0) status = -1;
	else                                                   status =  0;
//-----------------------------------------------------------------------------
// form DCache name for this file and add it to the chain, handle the case when 
// a file has been requested by name explicitly
// files cataloged in DFC are read through general DCACHE pool (cdfdca)
//-----------------------------------------------------------------------------
	if ((s_file == "") || (s_file == fn)) {
	  const char* path = 0;
	  GetDCacheFileName("cdfdca.fnal.gov",path,fs,fn,full_name);
	  Dataset->AddFile(full_name,fs,size,nev,loevt,lorun,hievt,hirun,status);
	}
      }
      gSystem->ClosePipe(pipe);
    }
    else if (strcmp(server,"sam") == 0) {
//-----------------------------------------------------------------------------
// query Oracle, directory = 
//-----------------------------------------------------------------------------
      TString ss = directory;
      TString s2 = ss+"."+fs;

      cmd  = Form("sqlplus cdf_reader/reader@%s <<EOF | ",GetOracleServer());
      cmd += Form("grep %s\n",s2.Data());
      cmd += Form("set linesize 500;\n");
      cmd += Form("set heading off;\n");
      cmd += Form("column location_id        format a40;\n");
      cmd += Form("column file_name          format a20;\n");
      cmd += Form("column proj_def_name      format a20;\n");
      cmd += Form("column full_path          format a60;\n");
      cmd += Form("column first_event_number format 999999999 ;\n");
      cmd += Form("column last_event_number  format 999999999 ;\n");
      cmd += Form("column proj_snap_version  format 99999 ;\n");

      cmd += Form("select smpjd.proj_def_name, smdf.file_name, smdf.file_size_in_bytes fsize, \n");
      cmd += Form(" 1 /* smdf.create_date */, \n");

      cmd += Form("smdf.event_count, \n");
      cmd += Form("(select min(smr.run_number) from data_files_runs smdfr, runs smr \n"); 
      cmd += Form(" where smdf.file_id = smdfr.file_id and smdfr.run_id = smr.run_id) low_run, \n");
      cmd += Form("smdf.first_event_number loevt, \n");

      cmd += Form("(select max(smr.run_number) from data_files_runs smdfr, runs smr \n"); 
      cmd += Form(" where smdf.file_id = smdfr.file_id and smdfr.run_id = smr.run_id) high_run, \n");
      cmd += Form("smdf.last_event_number hievt, sl.full_path \n");

      cmd += Form("from  project_definitions smpjd, project_snapshots smpjsnp, project_files smpjf, \n");
      cmd += Form("      data_files smdf, data_storage_locations sl, data_file_locations dfl \n");
      cmd += Form("where smpjd.proj_def_id = smpjsnp.proj_def_id and smpjsnp.proj_snap_id = smpjf.proj_snap_id and\n");
      cmd += Form("      smpjf.file_id = smdf.file_id and smpjd.proj_def_name = '%s' and \n",s2.Data());
      cmd += Form("      smdf.file_id = dfl.file_id and dfl.location_id = sl.location_id\n");
      cmd += Form("      and 'tape' in sl.location_type and substr(sl.full_path,1,20) = '/pnfs/cdfen/filesets'\n");
      //cmd += Form("and   low_run  <= %i and high_run >= %i",rmax,rmin);
      cmd += Form(";\nEOF\n");
      
      if (fPrintLevel > 0)
	printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());      
      pipe = gSystem->OpenPipe(cmd.Data(),"r");
      while (fscanf(pipe,"%s %s %f %i %i %i %i %i %i %s",
  		    fs,fn,&size,&time,&nev,&lorun,&loevt,&hirun,&hievt,pnfs_path)!=EOF) {
	if (fPrintLevel > 20)
	  printf("Result of SAM query: %s %s %f %i %i %i %i %i %i\n",
		 fs,fn,size,time,nev,lorun,loevt,hirun,hievt);
	

//-----------------------------------------------------------------------------
// so far list of bad files is kept locally ...
// make sure we're not adding a bad file
//-----------------------------------------------------------------------------
	if (Dataset->GetListOfBadFiles()->FindObject(fn) != 0) status = -1;
	else                                                   status =  0;
//-----------------------------------------------------------------------------
// form DCache name for this file and add it to the chain, handle the case when 
// a file has been requested by name explicitly
// files cataloged in DFC are read through general DCACHE pool (cdfdca)
//-----------------------------------------------------------------------------
	if ((s_file == "") || (s_file == fn)) {
	  //	  const char* path = 0;
	  GetDCacheFileName("cdfdca.fnal.gov",pnfs_path,0,fn,full_name);
	  if(lorun<=rmax && hirun>=rmin) {
	    Dataset->AddFile(full_name,fs,size,nev,loevt,lorun,hievt,hirun,status);
	  }
	}
      }
      gSystem->ClosePipe(pipe);
    }
    else {
//-----------------------------------------------------------------------------
// non-DFC case
// list of files for this fileset already retrieved and stored in 'files'
// filesets from this list have already passed all the name checks
//-----------------------------------------------------------------------------
      TString prefix;

      TIter itt(&files);
      while ((ostr = (TObjString*) itt.Next())) {
	line = (char*) ostr->String().Data();
	
	sscanf(line,"%s %s %f %s %s %i %i %i %i %i", 
	       fs,fn,&size,date,ctime,&nev,&lorun,&loevt,&hirun,&hievt);
//-----------------------------------------------------------------------------
// so far list of bad files is kept locally ...
// make sure we're not adding a bad file
//-----------------------------------------------------------------------------
	status =  0;
	for(int i=0; i<Dataset->GetListOfBadFiles()->GetEntries(); i++) {
	  TObjString* str = (TObjString*)Dataset->GetListOfBadFiles()->At(i);
	  TString fns(fn);
	  if(fns.Index(str->String().Data()) >=0) status = -1;
	}
	if (strcmp(fs,fileset->GetName()) == 0) {
//-----------------------------------------------------------------------------
//  file belongs to the fileset, form URL-like name
//-----------------------------------------------------------------------------
	  TUrl* srv = fileset->GetDataServer();
	  if (strcmp(srv->GetProtocol(),"cdf_dcap") == 0)    {
//-----------------------------------------------------------------------------
// file in the analysis DCACHE pool 
//-----------------------------------------------------------------------------
	    TString s    = srv->GetFile();
	    s.Prepend("/");
	    TString path = s(1,s.Length()-1);

	    GetDCacheFileName(srv->GetHost(),path.Data(),0,fn,full_name);
	  }
	  else if (strcmp(srv->GetProtocol(),"cdf") == 0) {
//-----------------------------------------------------------------------------
// file in the analysis DCACHE pool 
//-----------------------------------------------------------------------------
	    TString s    = srv->GetFile();
	    s.Prepend("/");
	    TString path = s(1,s.Length()-1);

	    GetDCacheFileName(srv->GetHost(),path.Data(),0,fn,full_name);
	  }
	  else {
//-----------------------------------------------------------------------------
// anything else, so far: ROOTD
//-----------------------------------------------------------------------------
	    sprintf(full_name,"%s",fn);
	  }

	  Dataset->AddFile(full_name,fs,size,nev,loevt,lorun,hievt,hirun,status);
	}
      }
    }
  NEXT_FILESET:;
  }

  filesets.Delete();
  files.Delete();

  return 0;
}


//_____________________________________________________________________________
Int_t TTxtCatalogServer::GetNEvents(const char* Book, 
				    const char* Dataset, 
				    const char* Fileset,
				    const char* File) 
{
  Int_t   nev;
  TString cmd;

  Error("GetNEvents","this works only in case of non-DFC file");

  cmd = Form("%s -b %s -d %s",GetNEventsCommand(),Book,Dataset);
  if (Fileset) cmd += Form(" -s %s",Fileset);
  if (File   ) cmd += Form(" -f %s",File   );

  FILE* pipe = gSystem->OpenPipe(cmd.Data(),"r");
  fscanf(pipe,"%i",&nev);
  gSystem->ClosePipe(pipe);

  return nev;
}

//_____________________________________________________________________________
Int_t TTxtCatalogServer::GetNFiles(const char* Book, 
				   const char* Dataset, 
				   const char* Fileset)
{
  Int_t   n = -1;

  Error("GetNFiles","Not implemented yet");
  return n;
}


//_____________________________________________________________________________
Int_t TTxtCatalogServer::GetNFilesets(const char* Book, const char* Dataset) {
  Int_t   n = -1;

  Error("GetNFilesets","Not implemented yet");
  return n;
}

//_____________________________________________________________________________
int TTxtCatalogServer::GetRemoteServer(const char* Book     ,
				       const char* Dataset  ,
				       const char* Fileset  ,
				       char*       Server   ,
				       char*       RemoteDir)
{
  // returns  smth like Server="ncdf131.fnal.gov" and 
  // RemoteDir="/cdf/scratch/data131" - address of the directory 
  // where the fileset is located

  char cmd[200];

  sprintf(cmd,"%s -b %s -d %s -s %s",
	  GetDataServerNameCommand(),
	  Book,Dataset,Fileset);

  FILE* f = gSystem->OpenPipe(cmd,"r");
  fscanf(f,"%s %s",Server,RemoteDir);
  gSystem->ClosePipe(f);

  return 0;
}

