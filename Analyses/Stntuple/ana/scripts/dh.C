///////////////////////////////////////////////////////////////////////////////
// example:
//
// root
// .L Stntuple/scripts/dh.C
//  dh::catalog_dataset("stntuple/dev_240","bpel08") ; > .bpel08.log
//
// cat results/bhel08-hpte-catalog.log | grep bhel08-hpte.0 \
//    | grep -v OBJ | grep -v Stntuple | grep -v mv | grep -v root://
///////////////////////////////////////////////////////////////////////////////
#include "string.h"


namespace dh {

  TStnAna*     x     = NULL;
  TChain*      chain;
  TStnDataset* dset;
  TDFCModule*  m_dfc;
  TStnCatalog* catalog;
  // Helping to organize filesets
  int          nSetStart    = 0;
  int          nFilesPerSet = 20;
  int          thisFileset  = 0;

//_____________________________________________________________________________
int  catalog_file(const char* InputFile,
                  Int_t       PrintLevel = 1,
                  const char* OutputDir  = ".",
                  const char* Dataset    = 0,
                  const char* Book       = "",
                  const char* DbID       = "production_file_catalog_write")
{
  // always process all the events
  // PrintLevel = 111 also would try to execute the commands

  catalog = new TStnCatalog();
  
  dset    = new TStnDataset();
  catalog->InitDataset(dset,"file",InputFile,"","");

  x = new TStnAna(dset);

  m_dfc = new TDFCModule;
  x->AddModule(m_dfc);

  m_dfc->SetOutputDir(OutputDir);

  if (Dataset)
    m_dfc->SetDataSet(Dataset,Book,DbID);

  m_dfc->SetPrintLevel(PrintLevel);
  m_dfc->SetPrintOpt("data");
  m_dfc->SetNFileset(thisFileset);
  if (PrintLevel > 100)
    m_dfc->SetExecCommands(PrintLevel/100);

  x->SetNEventsToReport(999999999);
  x->SetPrintLevel     (-999);

  if (PrintLevel > 9) {
    printf(" - run, InputFile=%s,OutputDir=%s,Dataset=%s,Book=%s,DbID=%s \n",
	   InputFile,OutputDir,Dataset,Book,DbID);
  }

  x->Run();

  int rc = m_dfc->ReturnCode();

  delete x;
  delete chain;
  delete catalog;

  return rc;
}


//_____________________________________________________________________________
void catalog_dataset(const char* Book,
                     const char* Dataset,
                     const char* Fileset = 0,
                     const char* RemoteDir = "ftp://ewk@fcdfdata030.fnal.gov:/cdf/scratch/ewk/gmbs08/stntuple/dev_240",
                     const char* DfcDsID = "sewk01",
                     const char* DfcBook = "cdfpewk") {

  const char* script = "/cdf/home/cdfopr/cafdfc/scripts/get_list_of_files";

  char  cmd[1000], fn[200];

  if ((Fileset != 0) && (Fileset[0] != 0)) {
    sprintf(cmd,"ssh fcdflnx2 %s -b %s -d %s -s %s  \| awk \'{print $2}\'",
            script,Book,Dataset,Fileset);
  }
  else {
    sprintf(cmd,"ssh fcdflnx2 %s -b %s -d %s \| awk \'{print $2}\'",
            script,Book,Dataset);
  }

  //  printf("executing: %s\n",cmd);

  FILE* f = gSystem->OpenPipe(cmd,"r");

  int i = 0;
  while (fscanf(f,"%s",fn) > 0) {
    printf("------------------- %s --------------------\n", fn);
    catalog_file(fn,11,RemoteDir,DfcDsID,DfcBook);
    i++;
  }

  fclose(f);
}

//_____________________________________________________________________________
void catalog_list_of_files(const char* RequestFile, Int_t Mode=11) {
  // catalog list of files returned by the shell command

  char  data_server[200], book     [100], dataset_id[100];
  char  output_dir [200], input_dir[200];

  FILE* pipe;

  char  cmd[1000];

  sprintf(cmd,"cat %s | awk '{ if ($2 == \"DATA_SERVER\") print $3}'",RequestFile);
  pipe = gSystem->OpenPipe(cmd,"r");
  fscanf(pipe,"%s",data_server);
  gSystem->ClosePipe(pipe);

  //  printf(">>> data_server = %s\n",data_server);

  sprintf(cmd,"cat %s | awk '{ if ($2 == \"BOOK\") print $3}'",RequestFile);
  pipe = gSystem->OpenPipe(cmd,"r");
  fscanf(pipe,"%s",book);
  gSystem->ClosePipe(pipe);

  sprintf(cmd,"cat %s | awk '{ if ($2 == \"DATASET\") print $3}'",RequestFile);
  pipe = gSystem->OpenPipe(cmd,"r");
  fscanf(pipe,"%s",dataset_id);
  gSystem->ClosePipe(pipe);

  sprintf(cmd,"cat %s | awk '{ if ($2 == \"OUTPUT_DIR\") print $3}'",RequestFile);
  pipe = gSystem->OpenPipe(cmd,"r");
  fscanf(pipe,"%s",output_dir);
  gSystem->ClosePipe(pipe);

  sprintf(cmd,"cat %s | awk '{ if ($2 == \"INPUT_DIR\") print $3}'",RequestFile);
  pipe = gSystem->OpenPipe(cmd,"r");
  fscanf(pipe,"%s",input_dir);
  gSystem->ClosePipe(pipe);

//   printf("output_dir    =  %s\n",output_dir   );
//   printf("input_dir     =  %s\n",input_dir    );
//   printf("book          =  %s\n",book         );
//   printf("dataset_id    =  %s\n",dataset_id   );
//
//     return 0;

  TChain*    chain = new TChain("STNTUPLE");
  TObjArray* list  = new TObjArray(100);

  char c[1000], file[1000], fn[1000];
  int done = 0;

  FILE* f = fopen(RequestFile,"r");
  if (! f) {
    return -1;
  }

  while ( ((c[0]=getc(f)) != EOF) && !done) {
                                        // check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
                                        // read next filename
      if (fscanf(f,"%s",file) > 0) {
//-----------------------------------------------------------------------------
// form remote file name 
//-----------------------------------------------------------------------------
        if (strcmp(data_server,".") != 0) {
          sprintf(fn,"%s/%s/%s",data_server,input_dir,file);
	}
	else {
          sprintf(fn,"%s/%s",input_dir,file);
	}

	//      printf("--- cataloging %s\n",fn);

	catalog_file(fn,Mode,output_dir,dataset_id,book);
      }
    }
    else {
      fgets(c,1000,f);
    }
  }
}


//_____________________________________________________________________________
void catalog_directory(const char* Dir,
                       const char* Pattern,
                       Int_t       PrintLevel = 11,
                       const char* OutputDir  = "",
                       const char* Dataset    = 0,
                       const char* Book       = "",
                       const char* DbID       ="production_file_catalog_write")
 {
  char cmd[200], fn[200];
  sprintf(cmd,"find %s/ -name \\*%s\\*",Dir,Pattern);
  printf("\nexecuting: %s\n",cmd);

  FILE* pipe = gSystem->OpenPipe(cmd,"r");
  int nFiles = 0;
  while (fscanf(pipe,"%s",fn) > 0) { nFiles++; }
  gSystem->ClosePipe(pipe);
  printf("NFiles found:   %d\n",nFiles);
  printf("NFiles per Set: %d\n",nFilesPerSet);
  int nSets = nFiles/nFilesPerSet;
  if (nFiles%nFilesPerSet) nSets++;
  printf("NSets:          %d  starting from: %d\n\n",nSets,nSetStart);

  const char* ds;
  if (Dataset)
    ds = Dataset;
  else
    ds = "UNKNOWN";
      
  printf("# -----------------------------------------------------------\n");
  printf("# Entry in - AAA_CATALOG.html\n");
  printf("# -----------------------------------------------------------\n");
  printf("# MetaName      Server       Directory ....\n");
  for (int i=nSetStart; i<nSetStart+nSets; i++)
    printf("%s.%04d local %s 1 1 999999\n",ds,i,Dir);
  printf("# -----------------------------------------------------------\n\n");

  printf("# -----------------------------------------------------------\n");
  printf("# Entries in - %s.files\n",ds);
  printf("# -----------------------------------------------------------\n");
  pipe = gSystem->OpenPipe(cmd,"r");
  nFiles = 0;
  thisFileset = nSetStart-1;
  while (fscanf(pipe,"%s",fn) > 0) {
    if (nFiles++%nFilesPerSet == 0)
      thisFileset++;
    if (PrintLevel > 9)
      printf("file: %s\n",fn);
    catalog_file(fn,PrintLevel,OutputDir,Dataset,Book,DbID);
  }
  gSystem->ClosePipe(pipe);
  printf("# -----------------------------------------------------------\n\n");
}


//_____________________________________________________________________________
// dir = usr@node:/dir for remote files, /dir for local
// output is AAA_CATALOG.html, AAA_FILES.txt
//
// if dataset is also given, then catalog will be based 
// on the renamed files and a set of commands to rename
// the files will be produced in  move_commands.csh
// If flags&1 then read existing AAA files and continue from that point
  void catalog_any_directory(const char* rdir, const char* pattern ="",
			   const char* dataset = NULL, int flag=0)
 {
   char cmd[200], fn[200], node[200], barenode[200],dir[200],tempstr[100];
   TObjArray doneFiles;
   bool local = (strstr(rdir,":")==0);
   bool rename = (dataset!=NULL);
   strncpy(node,"",0);
   strncpy(barenode,"",0);
   if(local) {
     strncpy(dir,rdir,200);
     sprintf(cmd,"find %s/ -name \\*%s\\*",dir,pattern);
   } else {
     strncpy(dir,strstr(rdir,":")+1,200);
     strncpy(node,rdir,200);
     *(strstr(node,":")) = '\0';
     if(strstr(node,"@")!=0) strncpy(barenode,strstr(node,"@")+1,200);
     printf("\nnode=%s    dir=%s\n",node,dir);
     sprintf(cmd,"ssh %s 'find %s/ -name \\*%s\\*'",node, dir,pattern);
   }
   printf("executing: %s\n",cmd);

   // if flag is set, then continue where previous catalog was left off
   int begfs = 0;
   doneFiles.Clear();
   if(flag&0x1) {
     char line[1000];
     FILE* cat = fopen("AAA_CATALOG.html","r");
     //while (fscanf(cat,"%s\n",line) > 0) {
     while (fgets(line,1000,cat) > 0) {
       if(line[0]!='#' && line[0]!='<') {
	 //printf("processing input line: %20s\n",line);
	 begfs++;
       }
     }
     fclose(cat);
     FILE* fil = fopen("AAA_FILES.txt","r");
     while (fgets(line,1000,fil) > 0) {
       if(line[0]!='#' && line[0]!='<') {
	 //printf("processing input file line: %20s\n",line);
	 sscanf(line,"%s %s",tempstr,fn);
	 doneFiles.Add(new TObjString(fn));
       }
     }
     fclose(fil);
     doneFiles.Compress();
     printf("Starting at: fs=%d   file=%d\n",begfs,doneFiles.GetEntries());
   }

   // preread the directory
   FILE* pipe = gSystem->OpenPipe(cmd,"r");
   int nFiles = 0;
   fscanf(pipe,"%s",fn); // first one is the directory
   while (fscanf(pipe,"%s",fn) > 0) { nFiles++; }
   gSystem->ClosePipe(pipe);

   printf("NFiles found:   %d\n",nFiles);
   printf("NFiles per Set: %d\n",nFilesPerSet);

   int nSets = nFiles/nFilesPerSet;
   if (nFiles%nFilesPerSet) nSets++;

   char mode[10]="w";
   if(flag&0x1) sprintf(mode,"%s","a");
   FILE* cat = fopen("AAA_CATALOG.html",mode);
   FILE* fil = fopen("AAA_FILES.txt",mode);
   FILE* mov = fopen("move_commands.csh",mode);

   if(flag&0x1==0) {
     fprintf(cat,"# -----------------------------------------------------------\n");
     fprintf(cat,"# catalog of dir = %s\n",rdir);
     if(dataset) {
       fprintf(cat,"# dataset = %s\n",dataset);
     }
     fprintf(cat,"# -----------------------------------------------------------\n");
     fprintf(cat,"# FileSet      Server       Directory   Nevt     LRun   LEvt  HRun  HEvt\n");
   }

   pipe = gSystem->OpenPipe(cmd,"r");
   int thisfs = begfs;
   int lrunfs = 999999;
   int levtfs = 1e9;
   int hrunfs = 0;
   int hevtfs = 0;
   int nevtfs = 0;
   int ifile = 0;
   int iskip = 0;

   //TDFCModule* dfc = new TDFCModule;
   TDFCModule* dfc;
   TStnAna* ana;
   char fullfn[500];

   fscanf(pipe,"%s",fn); // first one is the directory

   while (fscanf(pipe,"%s",fn) > 0) {
     //printf("x%sx%sx",((TObjString*)doneFiles[0])->String().Data(),fn);
     // extract the file name from the filespec
     TString temp(fn);
     TObjArray* oa = temp.Tokenize("/");
     temp = ((TObjString*)oa->Last())->String().Data();
     oa->Delete();
     // if this already is in the catalog, then skip it
     if(doneFiles.FindObject(temp.Data())) {
       printf("skipping %04d %s\n",iskip,fn);
       iskip++;
     } else {
       printf("cataloging %4d  %s\n",ifile,fn);
       if(local) {
	 sprintf(fullfn,"%s",fn);
       } else {
	 sprintf(fullfn,"root://%s/%s",barenode,fn);
       }
       //printf("fullfn %s\n",fullfn);
       dfc = new TDFCModule;
       ana = new TStnAna(fullfn);
       ana->AddModule(dfc);
       
       dfc->SetPrintLevel(0);
       
       if(rename) {
	 dfc->SetDataSet(dataset,"book","dbid");
       }
       ana->SetNEventsToReport(999999999);
       ana->SetPrintLevel(-999);
       ana->Run();
       
       int lrunf = dfc->fMinRunNumber;
       int levtf = dfc->fMinEventNumber;
       int hrunf = dfc->fMaxRunNumber;
       int hevtf = dfc->fMaxEventNumber;
       int nevtf = dfc->fNEvents;
       TString date = dfc->GetFileDate();
       TString fname;
       if(rename) {
	 fname = dfc->GetNewFileName();
	 fprintf(mov,"mv %s %s/%s",fn,dir,fname);
       } else {
	 fname = dfc->GetOldFileName();
       }
       
       float fsize = dfc->GetFileSize();
       
       fprintf(fil,"fs%04d     %s   %7.3f  %s   %7d %6d %8d %6d %8d\n",
	       thisfs,fname.Data(),fsize,date.Data(),
	       nevtf,lrunf,levtf,hrunf,hevtf);
       
       if(lrunf<lrunfs) {
	 lrunfs = lrunf;
	 levtfs = levtf;
       } elseif (lrunf==lrunfs && levtf<levtfs) {
	 levtfs = levtf;
       }
       if(hrunf>hrunfs) {
	 hrunfs = hrunf;
	 hevtfs = hevtf;
       } elseif (hrunf==hrunfs && hevtf>hevtfs) {
	 hevtfs = hevtf;
       }
       nevtfs += nevtf;
       
       ifile++;
       // update fileset
       if (ifile%nFilesPerSet == 0 || ifile==nFiles) {
	 fprintf(cat,"fs%04d     %s    %s   %7d %6d %8d %6d %8d\n",
		 thisfs,barenode,dir,nevtfs,lrunfs,levtfs,hrunfs,hevtfs);
	 lrunfs = 999999;
	 levtfs = 1e9;
	 hrunfs = 0;
	 hevtfs = 0;
	 nevtfs = 0;
	 fflush(cat);
	 fflush(fil);
	 thisfs++;
       }
       
       delete ana;
       //delete dfc;

     } // end if file already in doneFiles
   } // end big while loop over new files

   gSystem->ClosePipe(pipe);
   fclose(cat);
   fclose(fil);
   fclose(mov);

 }


//_____________________________________________________________________________
int concatenate(const char* ListOfFiles,
                const char* OutputFile ,
                const char* Book       ,
                const char* DatasetID  )
{
  // example: ListOfFiles = "/cdf/opr2/cdfopr/val/results/aa*"
  // note: this routine deletes concatenated files
  // do not concatenate if total size of the files pointed to
  // is less than 10 GB

  TStnAna* x;
  TChain*  chain;

  char cmd[200], fn[200];
  int  size;
  int  kMinTotSize = 10.;
//-----------------------------------------------------------------------------
//  loop over the files in results and move them to ./buffer in chunks of 1 GB
//-----------------------------------------------------------------------------
  int sum_size        = 0;
  int max_size        = 1024000000;
  float  tot_size;

  sprintf(cmd,"find %s -printf \"%s %p\n\"|awk '{n =n+$1}END{print n/1e9}'",
            ListOfFiles);
  FILE* f = gSystem->OpenPipe(cmd,"r");
  fscanf(f,"%f",&tot_size);

  printf("tot_size = %f\n",tot_size);

  if (tot_size < kMinTotSize) return -1; // finish

  sprintf(cmd,"find  %s -printf \"%s %p\n\"",ListOfFiles);
  FILE* f = gSystem->OpenPipe(cmd,"r");

  TObjArray* list_of_files = new TObjArray(100);

  while (fscanf(f,"%i %s",&size,fn) > 0) {
    if ((size < 1000000000) && (sum_size+size < 1500000000)) {
      sum_size += size;
      list_of_files->Add(new TObjString(fn));
      if (sum_size > max_size) break;
    }
    else {
                                // large file, deal with it separately
      if (sum_size < 500000000) {
        list_of_files->Add(new TObjString(fn));
      }
      break;
    }
  }

  fclose(f);
  int nfiles = list_of_files->GetEntriesFast();
  if (nfiles > 1) {
//-----------------------------------------------------------------------------
// concatenate

    chain = new TChain("STNTUPLE");

    for (int i=0; i<nfiles; i++) {
      TObjString* s = (TObjString*) list_of_files->At(i);
      chain->AddFile(s->String().Data(),TChain::kBigNumber);
    }

    x = new TStnAna(chain);
    x->GetInputModule()->SetPrintLevel(1);

    om = new TStnOutputModule(OutputFile);

    om->SetMaxFileSize(1800);
    x->SetOutputModule(om);

    //    om->DropDataBlock("L3SummaryBlock");

    x->Run();

    delete x;

//   chain->Delete();
//   delete chain;
  }
  else {
//-----------------------------------------------------------------------------
//  1 file do not need to concatenate, just move
//-----------------------------------------------------------------------------
    TObjString* s = (TObjString*) list_of_files->At(0);

    sprintf(cmd,"mv  %s  %s",s->String().Data(),OutputFile);
    gSystem->Exec(cmd);
  }
//-----------------------------------------------------------------------------
// supposedly done with the concatenation
// finally catalog the file and move it to sewk00
//-----------------------------------------------------------------------------
  int rc = catalog_file(OutputFile, 111, "",DatasetID,Book);
  //  printf(" --- dfc_ana: return code RC=%i\n",rc);

  if (rc == 0) {
//-----------------------------------------------------------------------------
// last step: if everything is OK, delete the input files
//-----------------------------------------------------------------------------
    for (int i=0; i<nfiles; i++) {
      TObjString* s = (TObjString*) list_of_files->At(i);
      sprintf(cmd,"rm  %s",s->String().Data());
      printf("executing: %s\n",cmd);
      gSystem->Exec(cmd);
    }
  }

  return 0;
}

//_____________________________________________________________________________
int concatenate() {

  char output_file[200], pid[100];

  const char* data_server   = gSystem->Getenv("DATA_SERVER"  );
  const char* list_of_files = gSystem->Getenv("LIST_OF_FILES");

  const char* book         = gSystem->Getenv("BOOK"          );
  const char* dataset_id   = gSystem->Getenv("DATASET_ID"    );
  const char* output_dir   = gSystem->Getenv("OUTPUT_DIR"    );

  if (data_server) output_dir = gSystem->Getenv("PWD"    );;

  sprintf(output_file,"%s/output.%i",output_dir,gSystem->GetPid());

  printf("list_of_files =  %s\n",list_of_files);
  printf("output_file   =  %s\n",output_file  );
  printf("output_dir    =  %s\n",output_dir   );
  printf("book          =  %s\n",book         );
  printf("dataset_id    =  %s\n",dataset_id   );

  //  return 0;

  TChain*    chain = new TChain("STNTUPLE");
  TObjArray* list  = new TObjArray(100);

  char  *cp, *file;
  char  delimitors[] = " \n";

  int n=strlen(list_of_files);

  cp = new char[n+1];

  strcpy(cp,list_of_files);           /* Make writable copy.  */
  file = strtok (cp, delimitors);      /* token => "words" */

  int nfiles=0;

  while (file) {
    printf(" new file : %s\n",file);

    if (data_server) {
                                // modify file name to read data remotely

      chain->AddFile(Form("root://%s/%s",data_server,file),
                     TChain::kBigNumber);
    }
    else {
      chain->AddFile(file,TChain::kBigNumber);
    }

    nfiles++;
    file = strtok (NULL,delimitors);
  }

  TStnAna* x = new TStnAna(chain);
  x->GetInputModule()->SetPrintLevel(1);

  TStnOutputModule* om = new TStnOutputModule(output_file);

  om->SetMaxFileSize(2000);
  x->SetOutputModule(om);

  x->Run();

  delete x;
//-----------------------------------------------------------------------------
// done with the concatenation, catalog the file and move it to sewk00
//-----------------------------------------------------------------------------
  int rc = dfc_ana(output_file, 111,output_dir,dataset_id,book);
  printf(" --- dfc_ana: return code RC=%i\n",rc);

return 0;
}

//_____________________________________________________________________________
int concatenate(const char* RequestFile, int JobNumber) {

  // request file contains all the information in it
  // output_dir is the remote dir to where the file has to be moved

  TString     cmd;
  TString     book;
  TString     output_dir;
  TString     dataset_id;

  char output_file[200], pid[100];

  FILE* pipe;

  const char  data_server[200], output_tcl_dir[200];
  char        c[1000], file[1000], fn[1000];
  int         done = 0;

  gInterpreter->LoadMacro(RequestFile);

  TChain*    chain = new TChain("STNTUPLE");
  TObjArray* list  = new TObjArray(100);

  init_chain(chain,JobNumber,output_dir,book,dataset_id);

  TUrl od(output_dir);

  if (strcmp(od.GetHost(),"") == 0) {
//------------------------------------------------------------------------------
// local output directory
//------------------------------------------------------------------------------
    sprintf(output_file,"%s/output.%i",output_dir.Data(),gSystem->GetPid());
  }
  else {
//------------------------------------------------------------------------------
// remote output directory
//------------------------------------------------------------------------------
    sprintf(output_file,"%s/output.%i",gSystem->Getenv("PWD"),gSystem->GetPid());
  }

  printf("output_file    =  %s\n",output_file       );
  printf("output_dir     =  %s\n",output_dir.Data() );
  printf("book           =  %s\n",book.Data()       );
  printf("dataset_id     =  %s\n",dataset_id.Data() );

  //  return 0;
  //  chain->GetListOfFiles()->Print();
  //  return 0;

  TStnAna* x = new TStnAna(chain);
  x->GetInputModule()->SetPrintLevel(1);

  TStnOutputModule* om = new TStnOutputModule(output_file);
//-----------------------------------------------------------------------------
//  work around I/O issue with L3 summary block
//-----------------------------------------------------------------------------
  om->SetMaxFileSize(1800);
  x->SetOutputModule(om);

  x->Run();

  delete x;
//-----------------------------------------------------------------------------
// done with the concatenation, file resides in a local directory, rename it
//-----------------------------------------------------------------------------
  int rc = catalog_file(output_file,111,output_dir,
			dataset_id.Data(),book.Data());

  printf(" --- dh::concatenate: return code RC=%i\n",rc);

  return 0;
}


//_____________________________________________________________________________
int concatenate_dataset(const char* Dataset) {

  TString ds = Dataset;

  if ((ds == "sewk00") || (ds == "express")) {

    while (concatenate("/cdf/opr2/cdfopr/val/results/aa*",
                       "/cdf/opr2/cdfopr/val/sewk00/aaaaaa.output",
                       "cdfpewk",
                       "sewk00") == 0) {
    }
  }
  else if (ds == "jewk00") {
    while (concatenate("/cdf/opr2/cdfopr/jpmm08/JPsi_*",
                       "/cdf/opr2/cdfopr/jewk00/aaaaaa.output",
                       "cdfpewk",
                       "jewk00") == 0) {};
  }

 }

//_____________________________________________________________________________
  int strip(const char* Book, const char* Dataset, int Run, int Event) {
    // prerequisites: - cdfopr package checked out
    //                - Stntuple shared libraries available
    //                - libStntuple_oracle.so available, cdfofprd defined

    TString    _cmd;
    TCdf2Files file;

    if (Book == "") Book = "filecatalog";

    file.fFILE_NAME = "";
    cdfofprd.FindFile(Book,Dataset,Run,Event,&file);

    if (file.fFILE_NAME == "") {
      return -1;
    }

    _cmd = Form("./cdfopr/scripts/strip_event -b %s -D %s -f %s -e %i:%i &",
               Book,Dataset,file.fFILE_NAME.Data(),Run,Event);
    printf(" executing %s\n",_cmd.Data());
    gSystem->Exec(_cmd.Data());
    return 0;
  }



//_____________________________________________________________________________
  int strip(const char* Book, const char* Dataset, int* EventList) {
    // usage - call from the higher level script, event list  looks like
    // [run,event,  run,event, run,event, ...., -1]
    // prerequisites: - cdfopr package checked out
    //                - Stntuple shared libraries available
    //                - libStntuple_oracle.so available, cdfofprd defined

    TString    _cmd;
    TCdf2Files file;

    int run, event, first(1);

    FILE*   f = fopen("strip_tmp.tcl","w");

    fprintf(f,"talk DHInput\n");
    fprintf(f,"  cache set DCACHE\n");

    if (Book == "") Book = "filecatalog";

    file.fFILE_NAME = "";

    int i=0;
    while (EventList[i] > 0) {
      run   = EventList[i ];
      event = EventList[i+1];
      cdfofprd.FindFile(Book,Dataset,run,event,&file);

      if (file.fFILE_NAME != "") {
	fprintf(f,"  include file %s book=%s\n",file.fFILE_NAME.Data(),Book);
	if (first) {
	  fprintf(f,"  selectEvents set      run=%i event=%i\n",run,event);
	  first = 0;
	}
	else       fprintf(f,"  selectEvents include  run=%i event=%i\n",run,event);
      }

      i+=2;
    }

    fprintf(f,"  show\n");
    fprintf(f,"exit\n");

    fprintf(f,"path enable AllPath\n");
    fprintf(f,"module disable PuffModule\n");
    fprintf(f,"module enable FileOutput\n");
    fprintf(f,"module output FileOutput\n");
    fprintf(f,"talk FileOutput\n");
    fprintf(f,"  output create AA ./%s_strip.root\n",Dataset);
    fprintf(f,"  output paths  AA AllPath\n");
    fprintf(f,"show\n");
    fprintf(f,"exit\n");
    fprintf(f,"begin\n");
    fprintf(f,"show all\n");
    fprintf(f,"exit\n");

    fclose(f);

    _cmd = Form("./cdfopr/scripts/strip_event -i %s &","strip_tmp.tcl");
    printf(" executing %s\n",_cmd.Data());
    //    gSystem->Exec(_cmd.Data());
    return 0;
  }



//_____________________________________________________________________________
  int print_dataset(const char* Book, const char* Dataset,
                    const char* Fileset = "", const char* File = "") {
///////////////////////////////////////////////////////////////////////////////
//  examples:
//                 print_dataset("stntuple/dev_240","gmbs08")
//                 print_dataset("stntuple","sewk00")
//                 print_dataset("stntuple/dev_241","gqcd0c")
//
///////////////////////////////////////////////////////////////////////////////
    TStnCatalog* catalog = new TStnCatalog();
    TStnDataset* dataset = new TStnDataset(Book,Dataset,Fileset,File);

    catalog->InitDataset(dataset);
    dataset->Print();

    printf(" --- dataset has %i events\n",dataset->GetNEvents());

//      TStnAna* x = new TStnAna(dataset);
//      x->PrintStat(100);
  }

//_____________________________________________________________________________
  int print_dataset(const char* Book, const char* Dataset,
                    int MinRun = 0, int MaxRun = 999999) {
///////////////////////////////////////////////////////////////////////////////
//  examples:
//                 print_dataset("stntuple/dev_240","gmbs08")
//                 print_dataset("stntuple","sewk00")
//                 print_dataset("stntuple/dev_241","gqcd0c")
//
///////////////////////////////////////////////////////////////////////////////
    TStnCatalog* catalog = new TStnCatalog();
    TStnDataset* dataset = new TStnDataset(Book,Dataset,MinRun,MaxRun);

    catalog->InitDataset(dataset);
    dataset->Print();

    printf(" --- dataset has %i events\n",dataset->GetNEvents());

//      TStnAna* x = new TStnAna(dataset);
//      x->PrintStat(100);
  }

}
