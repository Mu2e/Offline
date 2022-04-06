#include "Offline/Validation/inc/TValCompare.hh"
#include "TApplication.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TPDF.h"
#include "TROOT.h"
#include <fstream>
#include <iomanip>
#include <iostream>

// ClassImp(TValCompare)

//_____________________________________________________________________________
void TValCompare::Delete(Option_t* Opt) {
  if (fFile1) {
    if (fFile1->IsOpen()) fFile1->Close();
    delete fFile1;
    fFile1 = nullptr;
  }
  if (fFile2) {
    if (fFile2->IsOpen()) fFile2->Close();
    delete fFile2;
    fFile2 = nullptr;
  }
  // fPar.Clear();
  fList.Delete();
  fDirs.Delete();
}

//_____________________________________________________________________________
Int_t TValCompare::GetDirs() {
  fFile1 = TFile::Open(fFileN1.Data(), "READ");
  if (!fFile1) {
    printf("Error opening one file 1\n");
    return 1;
  }
  if (!fFile1->IsOpen()) {
    printf("Error opening file 1\n");
    return 1;
  }
  if (fFileN2.Length() > 0) {
    fFile2 = TFile::Open(fFileN2.Data(), "READ");
    if (!fFile2) {
      printf("Error opening file 2\n");
      return 1;
    }
    if (!fFile2->IsOpen()) {
      printf("Error opening file 2\n");
      return 1;
    }
  }

  // don't create objects in the file
  gROOT->cd();

  // scan the first file, make a list of directories

  TObjArray todo;    // temp list of directotires to check for subdirectories
  todo.Add(fFile1);  // prime with the top of file1

  int itodo = 0;
  TObject* oo;
  TKey* kk;

  while (itodo < todo.GetEntries()) {
    TDirectory* dd = (TDirectory*)todo[itodo];
    if (fVerbose > 9) printf("scanning %s\n", dd->GetName());
    fDirs.AddLast((TObject*)dd);

    // look in this directory for subdirectories
    TIter it(dd->GetListOfKeys());
    while ((kk = (TKey*)it.Next())) {
      oo = dd->Get(kk->GetName());
      if (oo->ClassName() == TString("TDirectoryFile") ||
          oo->ClassName() == TString("TDirectory")) {
        todo.AddLast(oo);
      }
    }
    itodo++;
  }

  // List the directories, if requested
  TIter itd = TIter(&fDirs);
  if (fVerbose > 5) {
    printf("List of %d directories in file 1:\n", fDirs.GetEntries());
    itd = TIter(&fDirs);
    TDirectory* di;
    while ((di = (TDirectoryFile*)itd.Next())) {
      printf("%s\n", di->GetPath());
    }
  }

  return 0;
}

//_____________________________________________________________________________
Int_t TValCompare::OneFile(Option_t* Opt) {
  Delete();  // remove any previous analysis

  // traverse the first file and get subdirectories
  int rc = GetDirs();
  if (rc != 0) return rc;

  // now process all objects in the list of directories
  TDirectory* di;
  TObject* hh;
  TKey* kk;
  TString path, combo;

  TIter itd = TIter(&fDirs);
  while ((di = (TDirectory*)itd.Next())) {
    path = di->GetPath();
    int ind = path.Index(":");  // strip leading file name
    path = path(ind + 2, path.Length());
    if (fVerbose > 1) printf("Processing %s\n", path.Data());

    // look in this directory for histograms
    TIter ith(di->GetListOfKeys());
    while ((kk = (TKey*)ith.Next())) {
      hh = di->Get(kk->GetName());

      if (hh->ClassName() == TString("TH1F") ||
          hh->ClassName() == TString("TH1D") ||
          hh->ClassName() == TString("TProfile") ||
          hh->ClassName() == TString("TEfficiency") ||
          hh->ClassName() == TString("TH2F") ||
          hh->ClassName() == TString("TH2D")) {
        // have to pack path and name together
        combo = path + "/" + TString(hh->GetName());
        ((TNamed*)hh)->SetName(combo);
        fList.Add(hh);
      }
    }

  }  // end loop over list of directories in file 1

  return 0;
}

//_____________________________________________________________________________
Int_t TValCompare::Analyze(Option_t* Opt) {
  Delete();  // remove any previous analysis

  // traverse the first file and get subdirectories
  int rc = GetDirs();
  if (rc != 0) return rc;

  // now process all objects in the list of directories
  TDirectory *di, *dj;
  TObject *o1, *o2;
  TKey* kk;

  TIter itd = TIter(&fDirs);
  while ((di = (TDirectory*)itd.Next())) {
    TString path = di->GetPath();
    int ind = path.Index(":");  // strip leading file name
    path = path(ind + 2, path.Length());
    if (fVerbose > 1) printf("Processing %s\n", path.Data());

    if (path.Length() == 0) {
      dj = fFile2;
    } else {
      dj = (TDirectory*)fFile2->GetObjectUnchecked(path.Data());
    }

    if (!dj) {
      if (fVerbose > 0)
        printf("Warning: did not find directory %s in file 2\n", path.Data());
      continue;
    }

    // look in this directory for histograms
    TIter ith(di->GetListOfKeys());
    while ((kk = (TKey*)ith.Next())) {
      o1 = di->Get(kk->GetName());
      o2 = dj->Get(kk->GetName());
      if (o2) {
        bool ok1, ok2;
        int htype = 0;
        TValHist* hh = nullptr;
        ok1 = o1->ClassName() == TString("TH1F") ||
              o1->ClassName() == TString("TH1D");
        ok2 = o2->ClassName() == TString("TH1F") ||
              o2->ClassName() == TString("TH1D");
        if (ok1 && ok2) {
          hh = new TValHistH((TH1*)o1, (TH1*)o2);
          htype = 1;
        }
        ok1 = o1->ClassName() == TString("TProfile");
        ok2 = o2->ClassName() == TString("TProfile");
        if (ok1 && ok2) {
          hh = new TValHistP((TProfile*)o1, (TProfile*)o2);
          htype = 2;
        }
        ok1 = o1->ClassName() == TString("TEfficiency");
        ok2 = o2->ClassName() == TString("TEfficiency");
        if (ok1 && ok2) {
          hh = new TValHistE((TEfficiency*)o1, (TEfficiency*)o2);
          htype = 3;
        }
        ok1 = o1->ClassName() == TString("TH2F") ||
              o1->ClassName() == TString("TH2D");
        ok2 = o2->ClassName() == TString("TH2F") ||
              o2->ClassName() == TString("TH2D");
        if (ok1 && ok2) {
          hh = new TValHist2((TH2*)o1, (TH2*)o2);
          htype = 4;
        }

        if (htype > 0) {
          hh->SetPar(fPar);
          hh->SetTag(path);
          hh->Analyze();
          fList.Add(hh);
        }
      }
    }

  }  // end loop over list of directories in file 1

  return 0;
}

//_____________________________________________________________________________
void TValCompare::Report(Option_t* Opt) {
  TIter it(&fList);
  TValHist* hh;

  printf("   KS     Frac   I     Sum1       Sum2        Name        Title\n");

  while ((hh = (TValHist*)it.Next())) {
    bool qStat = (hh->GetStatus() >= fMinStat && hh->GetStatus() <= fMaxStat);
    if (qStat) hh->Summary();
  }
}

//_____________________________________________________________________________
void TValCompare::Summary(Option_t* Opt) {
  int nPerfect = 0, nSkip = 0, nEmpty = 0, nTight = 0, nLoose = 0, nFail = 0,
      nCantCompare = 0, nUnknown = 0;
  TIter it(&fList);
  TValHist* hh;

  while ((hh = (TValHist*)it.Next())) {
    // if the title contains "[info]" then it is for info only,
    // not in comparison summary, for example, CPU time is expected to change
    TString title(hh->GetTitle());
    bool useInSummary = true;
    if (title.Index("[info]") >= 0) useInSummary = false;

    if (useInSummary) {
      if (hh->GetStatus() == TValHist::fPerfect)
        nPerfect++;
      else if (hh->GetStatus() == TValHist::fTight)
        nTight++;
      else if (hh->GetStatus() == TValHist::fLoose)
        nLoose++;
      else if (hh->GetStatus() == TValHist::fFail)
        nFail++;
      else if (hh->GetStatus() == TValHist::fCantCompare)
        nCantCompare++;
      else
        nUnknown++;
      if (hh->GetEmpty()) nEmpty++;
    } else {
      nSkip++;
    }
  }

  printf("TValCompare Status Summary:\n");
  printf("%5d Compared\n", fList.GetEntries());
  printf("%5d marked to skip\n", nSkip);
  printf("%5d had unknown status\n", nUnknown);
  printf("%5d could not be compared\n", nCantCompare);
  printf("%5d had at least one histogram empty\n", nEmpty);
  printf("%5d failed loose comparison\n", nFail + nCantCompare + nUnknown);
  printf("%5d passed loose comparison, failed tight\n", nLoose);
  printf("%5d passed tight comparison, not perfect match\n", nTight);
  printf("%5d had perfect match\n", nPerfect);
  printf("%5d passed loose or better\n", nPerfect + nTight + nLoose);
  printf("%5d passed tight or better\n", nPerfect + nTight);
}

//_____________________________________________________________________________
void TValCompare::Display(Option_t* Opt) {
  TString opt(Opt);
  opt.ToUpper();

  bool q12 = opt.Contains("1X2");
  bool q22 = opt.Contains("2X2");
  int cx = 700;  // canvas size
  int cy = 500;
  float sf = 1.0;  // scaling for font size

  int clim = 1;  // plots per page
  if (q12) clim = 2;
  if (q22) clim = 4;

  if (q12) {
    cy = (11.0 / 8.5) * cx;
    sf = 1.15;
  } else if (q22) {
    cx = 800;
    cy = cx;
    sf = 0.82;
  }

  TValHist* hh;
  TIter it(&fList);
  while ((hh = (TValHist*)it.Next())) hh->SetFontScale(sf);

  TCanvas* ccc = nullptr;

  // when running from a linked executable, graphics is turned off
  // and this turns it back on
  Bool_t qBatch = gROOT->IsBatch();
  gROOT->SetBatch(kFALSE);
  if (!gApplication) gApplication = new TApplication("App", 0, NULL);

  ccc = new TCanvas("ccc", "TValCompare", cx, cy);

  int ind = -1;
  int iccc = 0;
  std::string str;
  while (1) {
    printf("<CR>,b,q : ");
    std::getline(std::cin, str);
    if (str == "q" || str == "Q") {
      break;
    } else if (str == "b" || str == "B") {
      ind = (ind <= 0 ? 0 : ind - 1);
      hh = (TValHist*)fList[ind];
      while (ind > 0 &&
             (hh->GetStatus() < fMinStat || hh->GetStatus() > fMaxStat)) {
        ind--;
        hh = (TValHist*)fList[ind];
      }
    } else {
      ind++;
      hh = (TValHist*)fList[ind];
      while (ind < fList.GetEntries() &&
             (hh->GetStatus() < fMinStat || hh->GetStatus() > fMaxStat)) {
        ind++;
        hh = (TValHist*)fList[ind];
      }
      if (ind >= fList.GetEntries()) break;
    }
    // if(!ccc) {
    //        ccc = new TCanvas("ccc","TValCompare",cx,cy);
    //}
    if (iccc % clim == 0) {
      ccc->Clear();
      if (q22) ccc->Divide(2, 2);
      if (q12) ccc->Divide(1, 2);
    }
    ccc->cd((iccc++) % clim + 1);
    fList[ind]->Draw(Opt);
    ccc->Modified();
    ccc->Update();
  }
  delete ccc;
  ccc = nullptr;
  gROOT->SetBatch(qBatch);
}

//_____________________________________________________________________________
void TValCompare::SaveAs(const char* filename, Option_t* option) const {
  TString file(filename);
  TString opt(option);
  opt.ToUpper();

  bool q12 = opt.Contains("1X2");
  bool q22 = opt.Contains("2X2");
  bool qPDF = file.EndsWith("pdf");
  bool qWeb = file.EndsWith("html");
  if (!(qPDF || qWeb)) {
    printf("ERRR - SaveAs file does not end with pdf or html");
    return;
  }

  int cx = 700;  // canvas size
  int cy = 500;
  float sf = 1.0;  // scaling for font size

  int clim = 1;  // plots per page
  if (qPDF && !q22) q12 = true;
  if (q12) clim = 2;
  if (q22) clim = 4;

  if (q12) {
    cy = (11.0 / 8.5) * cx;
    sf = 1.15;
  } else if (q22) {
    cx = 800;
    cy = cx;
    sf = 0.82;
  }

  TValHist* hh;
  TIter it(&fList);
  while ((hh = (TValHist*)it.Next())) hh->SetFontScale(sf);

  TCanvas* ccc = nullptr;

  Bool_t qBatch = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  if (qPDF) {
    ccc = new TCanvas("ccc-p", "TValCompare-pdf", cx, cy);
    TPDF* pdf = new TPDF(file);
    TIter it(&fList);
    int ind = 0;
    while ((hh = (TValHist*)it.Next())) {
      if (hh->GetStatus() >= fMinStat && hh->GetStatus() <= fMaxStat) {
        if ((ind % clim) == 0) {
          ccc->Clear();
          if (q12) ccc->Divide(1, 2);
          if (q22) ccc->Divide(2, 2);
        }
        ccc->cd(ind % clim + 1);
        ccc->Update();
        hh->Draw();
        ind++;
      }
    }
    // ccc->Destructor();
    delete ccc;
    ccc = nullptr;
    pdf->Close();
  }

  if (qWeb) {
    std::ofstream inf;
    inf.open(file);
    if (!inf.is_open()) {
      printf("ERROR - could not open web file %s\n", file.Data());
      return;
    }
    TString dir("./");
    int c = int(file.Last('/'));
    if (c >= 0) dir = file(0, c + 1);

    inf << "<html>\n<body>\n ";
    inf << "<title>valCompare</title>\n";
    inf << "<BR><BR>\n";
    inf << "<h2>" << fFileN1 << " (hist)<BR>&nbsp&nbsp vs &nbsp&nbsp<BR>";
    inf << fFileN2 << " (dots)</h2>\n";
    inf << "<BR><BR>\n";
    inf << "<TABLE>\n";
    inf << "<TR><TD width=120 align=left>KS</TD>\n";
    inf << "<TD width=120 align=left>Fraction</TD>\n";
    inf << "<TD width=50 align=left>Status</TD>\n";
    inf << "<TD align=center>Name</TD>\n";
    inf << "<TD align=center>Title</TD>\n";
    inf << "</TR>\n";

    ccc = new TCanvas("ccc-w", "TValCompare", 700, 500);
    TIter it(&fList);
    TString color;
    TString gifFile, gifName;
    TString gifFileLog, gifNameLog;
    while ((hh = (TValHist*)it.Next())) {
      // TEfficiency does not handle log scale well
      bool qDoLog = (hh->ClassName() != TString("TValHistE"));
      if (hh->GetStatus() >= fMinStat && hh->GetStatus() <= fMaxStat) {
        hh->Draw();
        gifName = hh->GetTag() + "/" + hh->GetName();
        gifName.ReplaceAll("/", "_");
        gifNameLog = gifName;
        gifName.Append(".gif");
        gifNameLog.Append("_log.gif");
        gifFile = dir + gifName;
        gifFileLog = dir + gifNameLog;
        ccc->SaveAs(gifFile);
        if (qDoLog) {
          hh->Draw("log");
          ccc->SaveAs(gifFileLog);
        }

        inf << "<TR><TD>";
        for (int io = 0; io < 2; io++) {
          color = "Black";
          float res = (io == 0 ? hh->GetKsProb() : hh->GetFrProb());
          // if Fraction test and samples are independent
          bool frblack = (io == 1 && fPar.GetIndependent() != 0);
          if (!frblack) {
            color = "Red";
            if (res > fPar.GetLoose()) color = "Orange";
            if (res > fPar.GetTight()) color = "Green";
            if (hh->GetStatus() == TValHist::fPerfect) color = "DarkGreen";
          }
          inf << "&nbsp&nbsp<font color=" << color << "> " << std::setw(8)
              << std::setprecision(6) << res << "</TD><TD>";
        }
        inf << hh->GetStatus() << "</TD><TD>";
        inf << "<a href=\"" << gifName << "\">" << hh->GetTag() << "/"
            << hh->GetName() << "</a> ";
        if (qDoLog) {
          inf << " &nbsp <a href=\"" << gifNameLog << "\">log</a> ";
        }
        inf << "</TD><TD>";
        inf << hh->GetTitle();
        inf << "</TD></TR>\n";
      }  // if within range
    }    // hist loop
    delete ccc;
    ccc = nullptr;
  }

  gROOT->SetBatch(qBatch);
}

//_____________________________________________________________________________
void TValCompare::SaveAs1(const char* filename, Option_t* option) const {
  TString file(filename);
  TString opt(option);
  opt.ToUpper();

  bool qPDF = file.EndsWith("pdf");
  bool qWeb = file.EndsWith("html");
  if (!(qPDF || qWeb)) {
    printf("ERRR - SaveAs1 file does not end with pdf or html");
    return;
  }

  /*  only in ValHist
  TH1* hh;
  TIter it(&fList);
  while ( (hh = (TH1*) it.Next()) ) hh->SetFontScale(sf);
  */

  TCanvas* ccc = nullptr;

  Bool_t qBatch = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  if (qPDF) {
    bool q12 = opt.Contains("1X2");
    bool q22 = opt.Contains("2X2");

    int cx = 700;  // canvas size
    int cy = 500;

    int clim = 1;  // plots per page
    if (qPDF && !q22) q12 = true;
    if (q12) clim = 2;
    if (q22) clim = 4;

    if (q12) {
      cy = (11.0 / 8.5) * cx;
    } else if (q22) {
      cx = 800;
      cy = cx;
    }

    TH1* hh;

    ccc = new TCanvas("ccc-p", "TValCompare-pdf", cx, cy);
    TPDF* pdf = new TPDF(file);
    TIter it(&fList);
    int ind = 0;
    while ((hh = (TH1*)it.Next())) {
      if ((ind % clim) == 0) {
        ccc->Clear();
        if (q12) ccc->Divide(1, 2);
        if (q22) ccc->Divide(2, 2);
      }
      ccc->cd(ind % clim + 1);
      ccc->Update();
      hh->Draw();
      ind++;
    }
    // ccc->Destructor();
    delete ccc;
    ccc = nullptr;
    pdf->Close();
  }

  if (qWeb) {
    std::ofstream inf;
    inf.open(file);
    if (!inf.is_open()) {
      printf("ERROR - could not open web file %s\n", file.Data());
      return;
    }
    TString dir("./");
    int c = int(file.Last('/'));
    if (c >= 0) dir = file(0, c + 1);

    inf << "<html>\n<body>\n ";
    inf << "<title>valCompare one file</title>\n";
    inf << "<BR><BR>\n";
    inf << "<h2>" << fFileN1 << " &nbsp&nbsp </h2> <BR>";

    inf << "<TABLE>\n";
    inf << "<TD align=center>Name</TD>\n";
    inf << "<TD align=center>Title</TD>\n";
    inf << "</TR>\n";

    ccc = new TCanvas("ccc-w", "TValCompare", 700, 500);
    TIter it(&fList);
    TString color;
    TString gifFile, gifName;
    TString gifFileLog, gifNameLog;
    TH1* hh;
    while ((hh = (TH1*)it.Next())) {
      // TEfficiency does not handle log scale well
      bool qDoLog = (hh->ClassName() != TString("TEfficiency"));
      hh->Draw();
      // gifName = hh->GetTag()+"/"+hh->GetName();
      gifName = hh->GetName();
      gifName.ReplaceAll("/", "_");
      gifNameLog = gifName;
      gifName.Append(".gif");
      gifNameLog.Append("_log.gif");
      gifFile = dir + gifName;
      gifFileLog = dir + gifNameLog;
      ccc->SaveAs(gifFile);
      if (qDoLog) {
        gPad->SetLogy(1);
        ccc->SaveAs(gifFileLog);
        gPad->SetLogy(0);
      }

      inf << "<TR><TD>";

      inf << "<a href=\"" << gifName << "\">" << hh->GetName() << "</a> ";
      if (qDoLog) {
        inf << " &nbsp <a href=\"" << gifNameLog << "\">log</a> ";
      }
      inf << "</TD><TD>";
      inf << hh->GetTitle();
      inf << "</TD></TR>\n";

    }  // hist loop
    delete ccc;
    ccc = nullptr;
  }

  gROOT->SetBatch(qBatch);
}

//_____________________________________________________________________________
TValHist* TValCompare::GetHist(TString str) {
  TIter it(&fList);
  TValHist *hh, *ret = nullptr;

  bool match = false;
  TString temp;
  while ((hh = (TValHist*)it.Next())) {
    temp = hh->GetTag() + "/" + hh->GetName() +
           " "
           "" +
           hh->GetTitle() +
           ""
           "";
    if (temp.Contains(str, TString::kIgnoreCase)) {
      if (!match) {
        match = true;
        printf("matches:\n");
      }
      printf("%s\n", temp.Data());
      ret = hh;
    }
  }
  if (!match) {
    printf("no matches\n");
  }
  return ret;
}

/*
//_____________________________________________________________________________
Int_t TValCompare::Write(const char *name, Int_t option, Int_t bufsize) {


}

*/
