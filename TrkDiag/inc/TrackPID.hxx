//Code generated automatically by TMVA for Inference of Model file [TrackPID.onnx] at [Tue Sep 23 21:14:35 2025] 

#ifndef ROOT_TMVA_SOFIE_TRACKPID
#define ROOT_TMVA_SOFIE_TRACKPID

#include <algorithm>
#include <cmath>
#include <vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_TrackPID{
namespace BLAS{
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_sequential1dense1BiasAddReadVariableOp0 = std::vector<float>(5);
float * tensor_sequential1dense1BiasAddReadVariableOp0 = fTensor_sequential1dense1BiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense12CastReadVariableOp0 = std::vector<float>(50);
float * tensor_sequential1dense12CastReadVariableOp0 = fTensor_sequential1dense12CastReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense12BiasAddReadVariableOp0 = std::vector<float>(10);
float * tensor_sequential1dense12BiasAddReadVariableOp0 = fTensor_sequential1dense12BiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense21BiasAddReadVariableOp0 = std::vector<float>(5);
float * tensor_sequential1dense21BiasAddReadVariableOp0 = fTensor_sequential1dense21BiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense1CastReadVariableOp0 = std::vector<float>(20);
float * tensor_sequential1dense1CastReadVariableOp0 = fTensor_sequential1dense1CastReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense21CastReadVariableOp0 = std::vector<float>(50);
float * tensor_sequential1dense21CastReadVariableOp0 = fTensor_sequential1dense21CastReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense31AddReadVariableOp0 = std::vector<float>(1);
float * tensor_sequential1dense31AddReadVariableOp0 = fTensor_sequential1dense31AddReadVariableOp0.data();
std::vector<float> fTensor_sequential1dense31CastReadVariableOp0 = std::vector<float>(5);
float * tensor_sequential1dense31CastReadVariableOp0 = fTensor_sequential1dense31CastReadVariableOp0.data();

//--- declare and allocate the intermediate tensors
std::vector<float> fTensor_sequential1dense31AddReadVariableOp0bcast = std::vector<float>(32);
float * tensor_sequential1dense31AddReadVariableOp0bcast = fTensor_sequential1dense31AddReadVariableOp0bcast.data();
std::vector<float> fTensor_sequential1dense21Relu0 = std::vector<float>(160);
float * tensor_sequential1dense21Relu0 = fTensor_sequential1dense21Relu0.data();
std::vector<float> fTensor_sequential1dense21BiasAddReadVariableOp0bcast = std::vector<float>(160);
float * tensor_sequential1dense21BiasAddReadVariableOp0bcast = fTensor_sequential1dense21BiasAddReadVariableOp0bcast.data();
std::vector<float> fTensor_sequential1dense12Relu0 = std::vector<float>(320);
float * tensor_sequential1dense12Relu0 = fTensor_sequential1dense12Relu0.data();
std::vector<float> fTensor_sequential1dense12MatMulGemm80 = std::vector<float>(320);
float * tensor_sequential1dense12MatMulGemm80 = fTensor_sequential1dense12MatMulGemm80.data();
std::vector<float> fTensor_sequential1dense12BiasAddReadVariableOp0bcast = std::vector<float>(320);
float * tensor_sequential1dense12BiasAddReadVariableOp0bcast = fTensor_sequential1dense12BiasAddReadVariableOp0bcast.data();
std::vector<float> fTensor_sequential1dense21MatMulGemm90 = std::vector<float>(160);
float * tensor_sequential1dense21MatMulGemm90 = fTensor_sequential1dense21MatMulGemm90.data();
std::vector<float> fTensor_sequential1dense1Relu0 = std::vector<float>(160);
float * tensor_sequential1dense1Relu0 = fTensor_sequential1dense1Relu0.data();
std::vector<float> fTensor_output = std::vector<float>(32);
float * tensor_output = fTensor_output.data();
std::vector<float> fTensor_sequential1dense31MatMulGemm60 = std::vector<float>(32);
float * tensor_sequential1dense31MatMulGemm60 = fTensor_sequential1dense31MatMulGemm60.data();
std::vector<float> fTensor_sequential1dense1MatMulGemm70 = std::vector<float>(160);
float * tensor_sequential1dense1MatMulGemm70 = fTensor_sequential1dense1MatMulGemm70.data();
std::vector<float> fTensor_sequential1dense1BiasAddReadVariableOp0bcast = std::vector<float>(160);
float * tensor_sequential1dense1BiasAddReadVariableOp0bcast = fTensor_sequential1dense1BiasAddReadVariableOp0bcast.data();


Session(std::string filename ="TrackPID.dat") {

//--- reading weights from file
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()) {
      throw std::runtime_error("tmva-sofie failed to open file " + filename + " for input weights");
   }
   std::string tensor_name;
   size_t length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense1BiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense1BiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 5) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 5 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense1BiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense12CastReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense12CastReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 50) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 50 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense12CastReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense12BiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense12BiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 10) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 10 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense12BiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense21BiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense21BiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 5) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 5 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense21BiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense1CastReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense1CastReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 20) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 20 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense1CastReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense21CastReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense21CastReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 50) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 50 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense21CastReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense31AddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense31AddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense31AddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequential1dense31CastReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequential1dense31CastReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 5) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 5 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequential1dense31CastReadVariableOp0[i];
   f.close();

//---- allocate the intermediate dynamic tensors
//--- broadcast bias tensor sequential1dense1BiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequential1dense1BiasAddReadVariableOp0,{ 5 }, { 32 , 5 });
      std::copy(data, data + 160, tensor_sequential1dense1BiasAddReadVariableOp0bcast);
      delete [] data;
   }
//--- broadcast bias tensor sequential1dense12BiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequential1dense12BiasAddReadVariableOp0,{ 10 }, { 32 , 10 });
      std::copy(data, data + 320, tensor_sequential1dense12BiasAddReadVariableOp0bcast);
      delete [] data;
   }
//--- broadcast bias tensor sequential1dense21BiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequential1dense21BiasAddReadVariableOp0,{ 5 }, { 32 , 5 });
      std::copy(data, data + 160, tensor_sequential1dense21BiasAddReadVariableOp0bcast);
      delete [] data;
   }
//--- broadcast bias tensor sequential1dense31AddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequential1dense31AddReadVariableOp0,{ 1 }, { 32 , 1 });
      std::copy(data, data + 32, tensor_sequential1dense31AddReadVariableOp0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_inputlayer){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 32;
   int op_0_n = 5;
   int op_0_k = 4;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 4;
   int op_0_ldb = 5;
   std::copy(tensor_sequential1dense1BiasAddReadVariableOp0bcast, tensor_sequential1dense1BiasAddReadVariableOp0bcast + 160, tensor_sequential1dense1MatMulGemm70);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_sequential1dense1CastReadVariableOp0, &op_0_ldb, tensor_inputlayer, &op_0_lda, &op_0_beta, tensor_sequential1dense1MatMulGemm70, &op_0_n);

//------ RELU
   for (int id = 0; id < 160 ; id++){
      tensor_sequential1dense1Relu0[id] = ((tensor_sequential1dense1MatMulGemm70[id] > 0 )? tensor_sequential1dense1MatMulGemm70[id] : 0);
   }

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 32;
   int op_2_n = 10;
   int op_2_k = 5;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 5;
   int op_2_ldb = 10;
   std::copy(tensor_sequential1dense12BiasAddReadVariableOp0bcast, tensor_sequential1dense12BiasAddReadVariableOp0bcast + 320, tensor_sequential1dense12MatMulGemm80);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_sequential1dense12CastReadVariableOp0, &op_2_ldb, tensor_sequential1dense1Relu0, &op_2_lda, &op_2_beta, tensor_sequential1dense12MatMulGemm80, &op_2_n);

//------ RELU
   for (int id = 0; id < 320 ; id++){
      tensor_sequential1dense12Relu0[id] = ((tensor_sequential1dense12MatMulGemm80[id] > 0 )? tensor_sequential1dense12MatMulGemm80[id] : 0);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 32;
   int op_4_n = 5;
   int op_4_k = 10;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 10;
   int op_4_ldb = 5;
   std::copy(tensor_sequential1dense21BiasAddReadVariableOp0bcast, tensor_sequential1dense21BiasAddReadVariableOp0bcast + 160, tensor_sequential1dense21MatMulGemm90);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_sequential1dense21CastReadVariableOp0, &op_4_ldb, tensor_sequential1dense12Relu0, &op_4_lda, &op_4_beta, tensor_sequential1dense21MatMulGemm90, &op_4_n);

//------ RELU
   for (int id = 0; id < 160 ; id++){
      tensor_sequential1dense21Relu0[id] = ((tensor_sequential1dense21MatMulGemm90[id] > 0 )? tensor_sequential1dense21MatMulGemm90[id] : 0);
   }

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 32;
   int op_6_n = 1;
   int op_6_k = 5;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 5;
   int op_6_ldb = 1;
   std::copy(tensor_sequential1dense31AddReadVariableOp0bcast, tensor_sequential1dense31AddReadVariableOp0bcast + 32, tensor_sequential1dense31MatMulGemm60);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_sequential1dense31CastReadVariableOp0, &op_6_ldb, tensor_sequential1dense21Relu0, &op_6_lda, &op_6_beta, tensor_sequential1dense31MatMulGemm60, &op_6_n);
	for (int id = 0; id < 32 ; id++){
		tensor_output[id] = 1 / (1 + std::exp( - tensor_sequential1dense31MatMulGemm60[id]));
	}
   return fTensor_output;
}
};
} //TMVA_SOFIE_TrackPID

#endif  // ROOT_TMVA_SOFIE_TRACKPID
