# THMu2eCaloDisk documentation

`THMu2eCaloDisk` is a ROOT-like library inheriting from `TH2Poly`.
It's a simple and portable way to have calorimeter disk displays of any kind.
Each bin is an expanded `THMu2eCaloDiskBin` class and knows everything about the crystal and hardware properties.
Each bin holds both channels (left and right) and can calculate and display any quantity based on the two channels.
Each class object instance is a single calo disk, and can therefore be either Disk 0 or Disk 1.
Laser pin-diodes and LYSO crystals are represented too.

Code example:
```
mu2e::THMu2eCaloDisk* h_disk0 = new mu2e::THMu2eCaloDisk("h_disk0","Disk 0",0);
h_disk0->FillOffline(sipmId,3.14);       //Fill a channel based on offline SiPM ID
h_disk0->SetCombineMode(kAsymmetry);     //Set the quantity to be displayed for each crystal
h_disk0->Draw("colzl");                  //Draw the disk with borders and colors

mu2e::THMu2eCaloDisk* h_disk1 = new mu2e::THMu2eCaloDisk("h_disk1","Disk 1",1);
for (int board=80; board<160; board++){
  for (int channel=0; channel<20; channel++){
    h_disk1->FillRaw(board,channel,value);           //Fill a channel based on online board and channel numbers
                                                     //Non-valid board-channel combinations will be ignored
  }
}
h_disk1->SetCombineMode("(R<0.1)?0:L/R");   //Set the quantity to be displayed for each crystal with a custom formula
h_disk1->Draw("colzl");                              //Draw the disk with borders and colors
```

**Important**: The `SetCombineMode()` function changes the legacy bin content but not the content of the left and right channels.
The combination mode can be changed arbitrarily multiple times without affecting the histogram stored information.

## Constructor

```
THMu2eCaloDisk(const char* name, const char* title, Int_t disk)
```

The constructor behaves like a `TH2Poly`, requiring the object name, the title, but also requires the disk index (either 0 or 1).
The disk geometry and online/offline information is currently read from the `Offline/CaloConditions/data/caloDMAP_latest.dat` file.
Database and more file options will be added in the future.

## Filling modes

The bins can be filled with the following functions:

```
  void SetBinContentL(Int_t bin, Double_t content);
  void SetBinContentR(Int_t bin, Double_t content);
```
Set the bin content (either left or right) via the bin number.

```
  Int_t FillOffline(int SiPMId, Double_t w);
  Int_t FillRaw(int board, int channel, Double_t w);
```
Fill the channel via offline or online numbers.
It behaves like a histogram `Fill` function, i.e. increasing the content by 1 or by `w`.
The `Int_t` value returned is the bin number.
The function does nothing and returns `0` if the specified offlineID or board/channel combination do not exist for the disk.

## Display modes

The filling functions alter the internal `fContentL` and `fContentR` variables.
However, the quantity displayed in the plots is the legacy `fContent`.
The value of `fContent` is computed from `fContentL` and `fContentR` according to user's choice, with the following functions:

```
void SetBinCombineMode(Int_t bin, Int_t mode);
void SetCombineMode(Int_t mode);
```
Set the L/R combination mode according to one of the predefined modes:
```
ECombineMode::kLeft         // Content of the left channel
ECombineMode::kRight        // Content of the right channel
ECombineMode::kSum          // L+R, sum of the two channels
ECombineMode::kAverage      // 0.5*(L+R), average of the two channels
ECombineMode::kDifference   // L-R, difference of the two channels
ECombineMode::kAsymmetry    // (L-R)/(L+R), channel asymmetry. Zero if L+R==0.
```
The first function allows to set the combination mode to a single bin.

```
void SetCombineMode(const char* formula);
```
Set the combination according to a custom formula following the rules of `ROOT::TFormula`.
The two channels are represented by the symbols `L` and `R` respectively.
Use of `x` and `y`, or `x[0]` and `x[1]` is also fine but discouraged.
Invalid formulas will throw a compilation error.
Examples:
```
h_disk0->SetCombineMode("L*R");
h_disk0->SetCombineMode("sqrt(L^2+R^2)");
h_disk0->SetCombineMode("max(L,R)");
h_disk0->SetCombineMode("TMath::Landau(0.5*(L+R))");
h_disk0->SetCombineMode("(L>R)?L/R:R/L");
```

