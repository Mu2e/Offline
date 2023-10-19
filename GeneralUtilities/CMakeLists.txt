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
src/TwoDPointTest_main.cc
src/TwoDWeight.cc
src/VMInfo.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
