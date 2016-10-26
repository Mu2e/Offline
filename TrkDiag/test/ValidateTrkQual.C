#include <string>
void ValidateTrkQual(std::string filename="TrkQual.root"){
  TMVA::TMVAGui(filename.c_str());
}
