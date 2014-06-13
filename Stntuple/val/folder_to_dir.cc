///////////////////////////////////////////////////////////////////////////////
// read file with the folders in the , convert each folder into a
// subdirectory and write file out
///////////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TKey.h"

#include "Stntuple/val/stntuple_val_functions.hh"
//_____________________________________________________________________________
int make_new_tree(TDirectory* Dir, TFolder* Folder) {
  Dir->cd();

  TDirectory* dir = new TDirectory(Folder->GetName(),Folder->GetName(),"");

  TIter    it(Folder->GetListOfFolders());
  TObject* o;

//   printf(" ------------------- Dir: %s, new dir: %s\n",
// 	 Dir->GetName(),dir->GetName());

  dir->cd();

  while ((o = it.Next())) {
//     printf(" o->GetName, o->ClassName : %-20s %-20s\n",
// 	   o->GetName(),
// 	   o->ClassName());

    if (strcmp(o->ClassName(),"TFolder") == 0) {
      make_new_tree(dir,(TFolder*) o);
      //      dir->cd();
    }
    else if (! o->InheritsFrom("TStnModule")) {
      //      printf("gDirectory->GetPath = %s\n",gDirectory->GetPath());
      o->Write();
      //      gDirectory->GetListOfKeys()->Print();
    }
  }

  Dir->cd();
  return 0;
}


//_____________________________________________________________________________
int folder_to_dir(const char* InputFile, const char* OutputFile) {

  TDirectory  *dir1, *dir2;
  TFolder     *folder;

  dir1 = new TFile(InputFile);

  TKey* k = (TKey*) dir1->GetListOfKeys()->FindObject("Ana");
  if (strcmp(k->GetClassName(),"TFolder") != 0) return 0;
  folder = (TFolder*) dir1->Get("Ana");

  dir2 = new TFile(OutputFile,"recreate");
//-----------------------------------------------------------------------------
// found top folder called Ana, start conversion
// step 1: read the data in
//-----------------------------------------------------------------------------
  TH1::AddDirectory(0);
  make_new_tree(dir2,folder);

  //  dir2->Write();

  dir2->Close();

  delete dir2;

  return 0;
}
