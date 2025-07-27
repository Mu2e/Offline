cet_make_library(
    SOURCE
      src/Angles.cc
      src/Binning.cc
      src/CombineTwoDPoints.cc
      src/CsvReader.cc
      src/DigitalFiltering.cc
      src/HepTransform.cc
      src/LinePointPCA.cc
      src/LineSegmentPCA.cc
      src/MinMax.cc
      src/NUBinning.cc
      src/OrientationResolver.cc
      src/ParameterSetFromFile.cc
      src/ParameterSetHelpers.cc
      src/ParseCLI.cc
      src/PathnameWithNextVersion.cc
      src/RMS.cc
      src/RootNameTitleHelper.cc
      src/SplineInterpolation.cc
      src/splitLine.cc
      src/splitString.cc
      src/sqrtOrThrow.cc
      src/TimeUtility.cc
      src/toHex.cc
      src/to_string.cc
      src/trimInPlace.cc
      src/TwoBodyKinematics.cc
      src/TwoDPoint.cc
      src/TwoDWeight.cc
      src/VMInfo.cc
    LIBRARIES PUBLIC
      ROOT::MathMore
      ROOT::GenVector
      Boost::headers
      Boost::filesystem
      cetlib_except::cetlib_except
      CLHEP::CLHEP
      fhiclcpp::fhiclcpp
      canvas::canvas      
)

cet_make_exec(NAME TwoDPointTest
    SOURCE src/TwoDPointTest_main.cc
    LIBRARIES
      Offline::GeneralUtilities
      
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
