#ifndef ROOT_THMu2eCaloDisk_H
#define ROOT_THMu2eCaloDisk_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "TGraph.h"
#include "TH2Poly.h"

namespace mu2e {

enum ECombineMode { kLeft = 0, kRight = 1, kSum = 2, kAverage = 3, kDifference = 4, kAsymmetry = 5 };

struct channelInfo {
	int         board;
	int         chan;
	double      chx;
	double      chy;
	int         row;
	int         column;
	int         cryId;
	int         roid;
	std::string type;

	channelInfo() {}
	channelInfo(int b, int c, double x, double y, int row, int col, int cry, int ro, std::string type) : board(b), chan(c), chx(x), chy(y), row(row), column(col), cryId(cry), roid(ro), type(type) {}

	int getTypeInt() {
		if(type == "CAL" || type == "CALO" || type == "CsI") {
			return 0;
		} else if(type == "PIN" || type == "PIN-DIODE") {
			return 1;
		} else if(type == "LYSO" || type == "CAPHRI") {
			return 2;
		} else {
			return -1;
		}
	}
};

class THMu2eCaloDiskBin : public TH2PolyBin {
  public:
	friend class THMu2eCaloDisk;

	THMu2eCaloDiskBin();
	THMu2eCaloDiskBin(TObject* poly, Int_t bin_number);
	~THMu2eCaloDiskBin() override {}

	void Update();
	void Merge(const THMu2eCaloDiskBin* toMerge);
	void ClearStats();

	Int_t    GetCrystalId() const { return fCrystalId; }
	Int_t    GetOfflineIdL() const { return fOfflineIdL; }
	Int_t    GetOfflineIdR() const { return fOfflineIdR; }
	Int_t    GetBoardL() const { return fBoardL; }
	Int_t    GetBoardR() const { return fBoardR; }
	Int_t    GetChannelL() const { return fChannelL; }
	Int_t    GetChannelR() const { return fChannelR; }
	Int_t    GetType() const { return fType; }
	Int_t    GetRow() const { return fRow; }
	Int_t    GetColumn() const { return fColumn; }
	Double_t GetCenterX() const { return fCenterX; }
	Double_t GetCenterY() const { return fCenterY; }
	Double_t GetWidthX() const { return fWidthX; }
	Double_t GetWidthY() const { return fWidthY; }
	Int_t    GetColor() const { return fColor; }

	Double_t GetContentL() const { return fContentL; }
	Double_t GetContentR() const { return fContentR; }
	Double_t GetSum() const { return fSum; }
	Double_t GetAverage() const { return fAverage; }
	Double_t GetDifference() const { return fDifference; }
	Double_t GetAsymmetry() const { return fAsymmetry; }

	void SetDrawColor(Bool_t doColor = true) {
		if(doColor) {
			((TGraph*)fPoly)->SetLineColor(fColor);
		} else {
			((TGraph*)fPoly)->SetLineColor(kBlack);
		}
	}

  private:
	Int_t    fCrystalId;
	Int_t    fOfflineIdL;
	Int_t    fOfflineIdR;
	Int_t    fBoardL;
	Int_t    fBoardR;
	Int_t    fChannelL;
	Int_t    fChannelR;
	Int_t    fType;
	Int_t    fRow;
	Int_t    fColumn;
	Double_t fCenterX;
	Double_t fCenterY;
	Double_t fWidthX;
	Double_t fWidthY;
	Int_t    fColor;

	Double_t fContentL;
	Double_t fContentR;
	Double_t fSum;
	Double_t fAverage;
	Double_t fDifference;
	Double_t fAsymmetry;

  protected:
	int  fCombineMode;
	void FillL(Double_t w);
	void FillR(Double_t w);

	ClassDefOverride(THMu2eCaloDiskBin, 1)
};

class THMu2eCaloDisk : public TH2Poly {
  public:
	THMu2eCaloDisk();
	THMu2eCaloDisk(const char* name, const char* title, Int_t disk);
	~THMu2eCaloDisk() {}

	THMu2eCaloDiskBin* CreateBin(TObject* poly) override;

	void SetBinCombineMode(Int_t bin, Int_t mode);
	void SetCombineMode(Int_t mode);
	bool LoadMapFile(std::map<int, std::map<int, mu2e::channelInfo>>& output, const char* filename = "/home/mu2edaq/dev_calo/pgirotti_testing/srcs/Offline/CaloConditions/data/caloDMAP_latest.dat");
	bool LoadMapDB(std::map<int, std::map<int, mu2e::channelInfo>>& output);

	Int_t FillOffline(int SiPMId, Double_t w);
	Int_t FillRaw(int board, int channel, Double_t w);

  protected:
	constexpr static float wcry  = 34.3;
	constexpr static float xmin0 = -755;
	constexpr static float xmax0 = 755;
	constexpr static float ymin0 = -755;
	constexpr static float ymax0 = 755;
	int                    fDisk;
	int                    combineMode;

	ClassDefOverride(mu2e::THMu2eCaloDisk, 1)
};
}  // namespace mu2e

#endif
