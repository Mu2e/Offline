#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cctype>

#include "TGraph.h"
#include "TH2Poly.h"
#include "TList.h"
#include "TMultiGraph.h"
#include "TFormula.h"

#include "Offline/DataProducts/inc/CaloConst.hh"

#include "Offline/CaloVisualizer/inc/THMu2eCaloDisk.hh"

ClassImp(mu2e::THMu2eCaloDisk);

mu2e::THMu2eCaloDiskBin::THMu2eCaloDiskBin() {
  fCrystalId = -1;
  fOfflineIdL = -1;
  fOfflineIdR = -1;
  fBoardL = -1;
  fBoardR = -1;
  fChannelL = -1;
  fChannelR = -1;
  fType = -1;
  fRow = -1;
  fColumn = -1;
  fCenterX = 0.0;
  fCenterY = 0.0;
  fWidthX = 0.0;
  fWidthY = 0.0;
  fContentL = 0.0;
  fContentR = 0.0;
  fCombineMode = kSum;
  ClearStats();
}

mu2e::THMu2eCaloDiskBin::THMu2eCaloDiskBin(TObject* poly, Int_t bin_number) :
    TH2PolyBin(poly, bin_number) {
  fCrystalId = -1;
  fOfflineIdL = -1;
  fOfflineIdR = -1;
  fBoardL = -1;
  fBoardR = -1;
  fChannelL = -1;
  fChannelR = -1;
  fType = -1;
  fRow = -1;
  fColumn = -1;
  fCenterX = 0.0;
  fCenterY = 0.0;
  fWidthX = 0.0;
  fWidthY = 0.0;
  fContentL = 0.0;
  fContentR = 0.0;
  fCombineMode = kSum;
  ClearStats();
}

void mu2e::THMu2eCaloDiskBin::Update() {
  fSum = fContentL + fContentR;
  fAverage = 0.5 * (fContentL + fContentR);
  fDifference = fContentL - fContentR;
  fAsymmetry = (fSum != 0 ? fDifference / fSum : 0.0);
  if (_formula.IsValid()){
    fFormula = _formula.Eval(fContentL,fContentR);
  } else {
    fFormula = 0;
  }

  switch (fCombineMode) {
  case kLeft:
    fContent = fContentL;
    break;
  case kRight:
    fContent = fContentR;
    break;
  case kSum:
    fContent = fSum;
    break;
  case kAverage:
    fContent = fAverage;
    break;
  case kDifference:
    fContent = fDifference;
    break;
  case kAsymmetry:
    fContent = fAsymmetry;
    break;
  case kFormula:
    fContent = fFormula;
    break;
  default:
    fContent = fSum;
  }

  SetChanged(true);
}

void mu2e::THMu2eCaloDiskBin::SetContentL(Double_t content) {
  fContentL = content;
  Update();
  SetChanged(true);
}

void mu2e::THMu2eCaloDiskBin::SetContentR(Double_t content) {
  fContentR = content;
  Update();
  SetChanged(true);
}

void mu2e::THMu2eCaloDiskBin::SetFormula(const char* formula) {
  //Replace L and R with x and y
  std::string tformula(formula);
  for (size_t i=0; i<tformula.length(); i++) {
    if (tformula[i] == 'L' || tformula[i] == 'R') {
      //Check it's not surrounded by alphanumeric symbols
      bool left_alnum = (i > 0) && std::isalnum(static_cast<unsigned char>(tformula[i-1]));
      bool right_alnum = (i < tformula.length()-1) && std::isalnum(static_cast<unsigned char>(tformula[i+1]));
      if (!left_alnum && !right_alnum) {
        if (tformula[i] == 'L'){
          tformula[i] = 'x';
        } else if (tformula[i] == 'R'){
          tformula[i] = 'y';
        }
      }
    }
  }
  _formula.Compile(tformula.c_str());
  Update();
  SetChanged(true);
}

void mu2e::THMu2eCaloDiskBin::Merge(const mu2e::THMu2eCaloDiskBin* toMerge) {
  this->fContentL += toMerge->fContentL;
  this->fContentR += toMerge->fContentR;
  Update();
  SetChanged(true);
}

void mu2e::THMu2eCaloDiskBin::ClearStats() {
  fContentL = 0;
  fContentR = 0;
  Update();
  SetChanged(true);
}

mu2e::THMu2eCaloDiskBin* mu2e::THMu2eCaloDisk::CreateBin(TObject* poly) {
  if (!poly)
    return nullptr;

  if (fBins == nullptr) {
    fBins = new TList();
    fBins->SetOwner();
  }

  fNcells++;
  Int_t ibin = fNcells - kNOverflow;
  return new THMu2eCaloDiskBin(poly, ibin);
}

void mu2e::THMu2eCaloDiskBin::FillL(Double_t w) {
  fContentL += w;
  SetChanged(true);
  this->Update();
}

void mu2e::THMu2eCaloDiskBin::FillR(Double_t w) {
  fContentR += w;
  SetChanged(true);
  this->Update();
}

mu2e::THMu2eCaloDisk::THMu2eCaloDisk(const char* name, const char* title, Int_t disk) :
    TH2Poly(name, title, xmin0, xmax0, ymin0, ymax0) {
  fDisk = disk;

  std::vector<int> col_ind = {kSpring, kOrange, kBlue, kRed, kBlack};

  std::map<int, std::map<int, mu2e::channelInfo>> diskmap;

  if (!LoadMapFile(diskmap)) {
    std::cout << "Failed constructing the THMu2eCaloDisk\n";
    return;
  }

  std::map<int, int> cryId_map;

  for (auto boardmap : diskmap) {
    for (auto channelmap : boardmap.second) {
      mu2e::channelInfo thisChannel = channelmap.second;

      if (thisChannel.type == "EMPTY")
        continue;
      if (thisChannel.board < 80 && disk == 1)
        continue;
      if (thisChannel.board >= 80 && disk == 0)
        continue;

      if (cryId_map.find(thisChannel.cryId) != cryId_map.end()) { // We created this crystal aready
        THMu2eCaloDiskBin* existingBin =
            (THMu2eCaloDiskBin*)fBins->At(cryId_map[thisChannel.cryId] - 1);
        if (thisChannel.roid == 0) {
          existingBin->fOfflineIdL = thisChannel.cryId * 2 + thisChannel.roid;
          existingBin->fBoardL = thisChannel.board;
          existingBin->fChannelL = thisChannel.chan;
        } else {
          existingBin->fOfflineIdR = thisChannel.cryId * 2 + thisChannel.roid;
          existingBin->fBoardR = thisChannel.board;
          existingBin->fChannelR = thisChannel.chan;
        }
        continue; // No need to go further
      }

      double y1 = thisChannel.chy - 0.5 * wcry;
      double y2 = thisChannel.chy + 0.5 * wcry;
      double x1 = thisChannel.chx - 0.5 * wcry;
      double x2 = thisChannel.chx + 0.5 * wcry;
      double x[5] = {x1, x2, x2, x1, x1};
      double y[5] = {y2, y2, y1, y1, y2};

      TGraph* g_bin = new TGraph(5, x, y);
      g_bin->SetLineWidth(1);
      g_bin->SetName(Form("#frac{%03d}{%02d}", thisChannel.board, thisChannel.chan));
      Int_t ibin = AddBin(g_bin);
      cryId_map[thisChannel.cryId] = ibin;

      THMu2eCaloDiskBin* thisBin = (THMu2eCaloDiskBin*)fBins->At(ibin - 1);
      thisBin->fCrystalId = thisChannel.cryId;
      if (thisChannel.roid == 0) {
        thisBin->fOfflineIdL = thisChannel.cryId * 2 + thisChannel.roid;
        thisBin->fBoardL = thisChannel.board;
        thisBin->fChannelL = thisChannel.chan;
      } else {
        thisBin->fOfflineIdR = thisChannel.cryId * 2 + thisChannel.roid;
        thisBin->fBoardR = thisChannel.board;
        thisBin->fChannelR = thisChannel.chan;
      }
      thisBin->fType = thisChannel.getTypeInt();
      thisBin->fRow = thisChannel.row;
      thisBin->fColumn = thisChannel.column;
      thisBin->fCenterX = thisChannel.chx;
      thisBin->fCenterY = thisChannel.chy;
      thisBin->fWidthX = wcry;
      thisBin->fWidthY = wcry;
      thisBin->fColor = col_ind[(thisChannel.board / 2) % col_ind.size()];
      thisBin->fCombineMode = combineMode;
    }
  }

  GetXaxis()->SetTitle("X [mm]");
  GetYaxis()->SetTitle("Y [mm]");
  SetDrawOption("COLZL");
}

bool mu2e::THMu2eCaloDisk::LoadMapFile(std::map<int, std::map<int, mu2e::channelInfo>>& output,
                                       const char* filename) {
  // Read channel map
  std::ifstream fmap;
  fmap.open(filename);
  if (!fmap.is_open()) {
    std::cout << "Couldn't open file " << filename << "\n";
    return false;
  }

  fmap.ignore(256, '\n');
  while (!fmap.eof()) {
    double temp, x, y;
    int row, column;
    int board, channel, cryid, roid;
    std::string type;
    fmap >> row >> column >> temp >> temp >> temp >> temp >> temp;
    fmap >> board >> temp >> temp >> channel;
    fmap >> temp >> x >> y >> cryid >> roid >> type;
    if (fmap.eof())
      break;

    output[board][channel] =
        mu2e::channelInfo(board, channel, x, y, row, column, cryid, roid, type);
  }
  fmap.close();

  return true;
}

bool mu2e::THMu2eCaloDisk::LoadMapDB(std::map<int, std::map<int, mu2e::channelInfo>>& output) {
  // TODO
  return false;
}

void mu2e::THMu2eCaloDisk::SetBinCombineMode(Int_t bin, Int_t mode) {
  THMu2eCaloDiskBin* thisBin = (THMu2eCaloDiskBin*)fBins->At(bin - 1);
  thisBin->fCombineMode = mode;
  thisBin->Update();
}

void mu2e::THMu2eCaloDisk::SetCombineMode(Int_t mode) {
  TIter next(fBins);
  TObject* obj;
  THMu2eCaloDiskBin* bin;
  while ((obj = next())) {
    bin = (THMu2eCaloDiskBin*)obj;
    bin->fCombineMode = mode;
    bin->Update();
  }
}

void mu2e::THMu2eCaloDisk::SetCombineMode(const char* formula) {
  TIter next(fBins);
  TObject* obj;
  THMu2eCaloDiskBin* bin;
  while ((obj = next())) {
    bin = (THMu2eCaloDiskBin*)obj;
    bin->fCombineMode = kFormula;
    bin->SetFormula(formula);
    bin->Update();
  }
}

Int_t mu2e::THMu2eCaloDisk::FillOffline(int SiPMId, Double_t w) {
  TIter next(fBins);
  TObject* obj;
  THMu2eCaloDiskBin* bin;

  while ((obj = next())) {
    bin = (THMu2eCaloDiskBin*)obj;
    if (SiPMId == bin->GetOfflineIdL()) {
      bin->FillL(w);
      fEntries++;
      SetBinContentChanged(kTRUE);
      return bin->GetBinNumber();
    } else if (SiPMId == bin->GetOfflineIdR()) {
      bin->FillR(w);
      fEntries++;
      SetBinContentChanged(kTRUE);
      return bin->GetBinNumber();
    }
  }

  return 0;
}

Int_t mu2e::THMu2eCaloDisk::FillRaw(int board, int channel, Double_t w) {
  TIter next(fBins);
  TObject* obj;
  THMu2eCaloDiskBin* bin;

  while ((obj = next())) {
    bin = (THMu2eCaloDiskBin*)obj;
    if (board == bin->GetBoardL() && channel == bin->GetChannelL()) {
      bin->FillL(w);
      fEntries++;
      SetBinContentChanged(kTRUE);
      return bin->GetBinNumber();
    } else if (board == bin->GetBoardR() && channel == bin->GetChannelR()) {
      bin->FillR(w);
      fEntries++;
      SetBinContentChanged(kTRUE);
      return bin->GetBinNumber();
    }
  }

  return 0;
}

void mu2e::THMu2eCaloDisk::Scale(Double_t c1, Option_t* option) {
  for (int i = 0; i < this->GetNumberOfBins(); i++) {
    this->SetBinContentL(i + 1, c1 * this->GetBinContentL(i + 1));
    this->SetBinContentR(i + 1, c1 * this->GetBinContentR(i + 1));
  }
  for (int i = 0; i < kNOverflow; i++) {
    this->SetBinContentL(-i - 1, c1 * this->GetBinContentL(-i - 1));
    this->SetBinContentR(-i - 1, c1 * this->GetBinContentR(-i - 1));
  }
}

Double_t mu2e::THMu2eCaloDisk::GetBinContentL(Int_t bin) const {
  if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow)
    return 0;
  if (bin < 0)
    return fOverflow[-bin - 1];
  return ((THMu2eCaloDiskBin*)fBins->At(bin - 1))->GetContentL();
}

Double_t mu2e::THMu2eCaloDisk::GetBinContentR(Int_t bin) const {
  if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow)
    return 0;
  if (bin < 0)
    return fOverflow[-bin - 1];
  return ((THMu2eCaloDiskBin*)fBins->At(bin - 1))->GetContentR();
}

void mu2e::THMu2eCaloDisk::SetBinContentL(Int_t bin, Double_t content) {
  if (bin > GetNumberOfBins() || bin == 0 || bin < -9)
    return;
  if (bin > 0) {
    ((THMu2eCaloDiskBin*)fBins->At(bin - 1))->SetContentL(content);
  } else {
    fOverflow[-bin - 1] = content;
  }

  SetBinContentChanged(kTRUE);
  fEntries++;
}

void mu2e::THMu2eCaloDisk::SetBinContentR(Int_t bin, Double_t content) {
  if (bin > GetNumberOfBins() || bin == 0 || bin < -9)
    return;
  if (bin > 0) {
    ((THMu2eCaloDiskBin*)fBins->At(bin - 1))->SetContentR(content);
  } else {
    fOverflow[-bin - 1] = content;
  }

  SetBinContentChanged(kTRUE);
  fEntries++;
}

void mu2e::THMu2eCaloDisk::Reset(Option_t* opt) {
  TIter next(fBins);
  TObject* obj;
  THMu2eCaloDiskBin* bin;

  // Clears the bin contents
  while ((obj = next())) {
    bin = (THMu2eCaloDiskBin*)obj;
    bin->ClearContent();
  }

  TH2::Reset(opt);
}
