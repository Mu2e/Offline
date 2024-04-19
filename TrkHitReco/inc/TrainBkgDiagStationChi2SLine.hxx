//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiagStationChi2SLine.h5] at [Thu Apr 11 20:23:54 2024] 

#ifndef ROOT_TMVA_SOFIE_TRAINBKGDIAGSTATIONCHI2SLINE
#define ROOT_TMVA_SOFIE_TRAINBKGDIAGSTATIONCHI2SLINE

#include <algorithm>
#include <cmath>
#include <vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_TrainBkgDiagStationChi2SLine{
namespace BLAS{
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_dense3bias0 = std::vector<float>(1);
float * tensor_dense3bias0 = fTensor_dense3bias0.data();
std::vector<float> fTensor_dense3kernel0 = std::vector<float>(24);
float * tensor_dense3kernel0 = fTensor_dense3kernel0.data();
std::vector<float> fTensor_dense2bias0 = std::vector<float>(24);
float * tensor_dense2bias0 = fTensor_dense2bias0.data();
std::vector<float> fTensor_dense2kernel0 = std::vector<float>(576);
float * tensor_dense2kernel0 = fTensor_dense2kernel0.data();
std::vector<float> fTensor_dense1bias0 = std::vector<float>(24);
float * tensor_dense1bias0 = fTensor_dense1bias0.data();
std::vector<float> fTensor_dense1kernel0 = std::vector<float>(576);
float * tensor_dense1kernel0 = fTensor_dense1kernel0.data();
std::vector<float> fTensor_densebias0 = std::vector<float>(24);
float * tensor_densebias0 = fTensor_densebias0.data();
std::vector<float> fTensor_densekernel0 = std::vector<float>(288);
float * tensor_densekernel0 = fTensor_densekernel0.data();
std::vector<float> fTensor_dense3Dense = std::vector<float>(1);
float * tensor_dense3Dense = fTensor_dense3Dense.data();
std::vector<float> fTensor_dense2Dense = std::vector<float>(24);
float * tensor_dense2Dense = fTensor_dense2Dense.data();
std::vector<float> fTensor_dense3bias0bcast = std::vector<float>(1);
float * tensor_dense3bias0bcast = fTensor_dense3bias0bcast.data();
std::vector<float> fTensor_dense1Relu0 = std::vector<float>(24);
float * tensor_dense1Relu0 = fTensor_dense1Relu0.data();
std::vector<float> fTensor_dense3Sigmoid0 = std::vector<float>(1);
float * tensor_dense3Sigmoid0 = fTensor_dense3Sigmoid0.data();
std::vector<float> fTensor_dense2Relu0 = std::vector<float>(24);
float * tensor_dense2Relu0 = fTensor_dense2Relu0.data();
std::vector<float> fTensor_dense1Dense = std::vector<float>(24);
float * tensor_dense1Dense = fTensor_dense1Dense.data();
std::vector<float> fTensor_dense1bias0bcast = std::vector<float>(24);
float * tensor_dense1bias0bcast = fTensor_dense1bias0bcast.data();
std::vector<float> fTensor_denseRelu0 = std::vector<float>(24);
float * tensor_denseRelu0 = fTensor_denseRelu0.data();
std::vector<float> fTensor_denseDense = std::vector<float>(24);
float * tensor_denseDense = fTensor_denseDense.data();
std::vector<float> fTensor_dense2bias0bcast = std::vector<float>(24);
float * tensor_dense2bias0bcast = fTensor_dense2bias0bcast.data();
std::vector<float> fTensor_densebias0bcast = std::vector<float>(24);
float * tensor_densebias0bcast = fTensor_densebias0bcast.data();


Session(std::string filename ="") {
   if (filename.empty()) filename = "TrainBkgDiagStationChi2SLine.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()) {
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   size_t length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense3bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense3bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense3bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense3kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense3kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense3kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense2bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense2bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense2bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense2kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense2kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 576) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 576 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense2kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense1bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense1bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense1bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense1kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense1kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 576) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 576 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense1kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_densebias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_densebias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_densebias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_densekernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_densekernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 288) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 288 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_densekernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_densebias0,{ 24 }, { 1 , 24 });
      std::copy(data, data + 24, tensor_densebias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense1bias0,{ 24 }, { 1 , 24 });
      std::copy(data, data + 24, tensor_dense1bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense2bias0,{ 24 }, { 1 , 24 });
      std::copy(data, data + 24, tensor_dense2bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense3bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense3bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input179){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 24;
   int op_0_k = 12;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 12;
   int op_0_ldb = 24;
   std::copy(tensor_densebias0bcast, tensor_densebias0bcast + 24, tensor_denseDense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_densekernel0, &op_0_ldb, tensor_input179, &op_0_lda, &op_0_beta, tensor_denseDense, &op_0_n);

//------ RELU
   for (int id = 0; id < 24 ; id++){
      tensor_denseRelu0[id] = ((tensor_denseDense[id] > 0 )? tensor_denseDense[id] : 0);
   }

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 1;
   int op_2_n = 24;
   int op_2_k = 24;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 24;
   int op_2_ldb = 24;
   std::copy(tensor_dense1bias0bcast, tensor_dense1bias0bcast + 24, tensor_dense1Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense1kernel0, &op_2_ldb, tensor_denseRelu0, &op_2_lda, &op_2_beta, tensor_dense1Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 24 ; id++){
      tensor_dense1Relu0[id] = ((tensor_dense1Dense[id] > 0 )? tensor_dense1Dense[id] : 0);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 1;
   int op_4_n = 24;
   int op_4_k = 24;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 24;
   int op_4_ldb = 24;
   std::copy(tensor_dense2bias0bcast, tensor_dense2bias0bcast + 24, tensor_dense2Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense2kernel0, &op_4_ldb, tensor_dense1Relu0, &op_4_lda, &op_4_beta, tensor_dense2Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 24 ; id++){
      tensor_dense2Relu0[id] = ((tensor_dense2Dense[id] > 0 )? tensor_dense2Dense[id] : 0);
   }

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 1;
   int op_6_n = 1;
   int op_6_k = 24;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 24;
   int op_6_ldb = 1;
   std::copy(tensor_dense3bias0bcast, tensor_dense3bias0bcast + 1, tensor_dense3Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense3kernel0, &op_6_ldb, tensor_dense2Relu0, &op_6_lda, &op_6_beta, tensor_dense3Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense3Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense3Dense[id]));
	}
   std::vector<float> ret (tensor_dense3Sigmoid0, tensor_dense3Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiagStationChi2SLine

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAGSTATIONCHI2SLINE
