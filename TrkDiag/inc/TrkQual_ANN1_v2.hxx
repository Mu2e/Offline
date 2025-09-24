//Code generated automatically by TMVA for Inference of Model file [TrkQual_ANN1_v2.onnx] at [Wed Sep 24 16:57:16 2025] 

#ifndef ROOT_TMVA_SOFIE_TRKQUAL_ANN1_V2
#define ROOT_TMVA_SOFIE_TRKQUAL_ANN1_V2

#include <algorithm>
#include <cmath>
#include <vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_TrkQual_ANN1_v2{
namespace BLAS{
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_sequentialdenseBiasAddReadVariableOp0 = std::vector<float>(7);
float * tensor_sequentialdenseBiasAddReadVariableOp0 = fTensor_sequentialdenseBiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequentialdense1MatMulReadVariableOp0 = std::vector<float>(49);
float * tensor_sequentialdense1MatMulReadVariableOp0 = fTensor_sequentialdense1MatMulReadVariableOp0.data();
std::vector<float> fTensor_sequentialdense1BiasAddReadVariableOp0 = std::vector<float>(7);
float * tensor_sequentialdense1BiasAddReadVariableOp0 = fTensor_sequentialdense1BiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequentialdense2BiasAddReadVariableOp0 = std::vector<float>(6);
float * tensor_sequentialdense2BiasAddReadVariableOp0 = fTensor_sequentialdense2BiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequentialdense2MatMulReadVariableOp0 = std::vector<float>(42);
float * tensor_sequentialdense2MatMulReadVariableOp0 = fTensor_sequentialdense2MatMulReadVariableOp0.data();
std::vector<float> fTensor_sequentialdense3BiasAddReadVariableOp0 = std::vector<float>(1);
float * tensor_sequentialdense3BiasAddReadVariableOp0 = fTensor_sequentialdense3BiasAddReadVariableOp0.data();
std::vector<float> fTensor_sequentialdenseMatMulReadVariableOp0 = std::vector<float>(49);
float * tensor_sequentialdenseMatMulReadVariableOp0 = fTensor_sequentialdenseMatMulReadVariableOp0.data();
std::vector<float> fTensor_sequentialdense3MatMulReadVariableOp0 = std::vector<float>(6);
float * tensor_sequentialdense3MatMulReadVariableOp0 = fTensor_sequentialdense3MatMulReadVariableOp0.data();

//--- declare and allocate the intermediate tensors
std::vector<float> fTensor_sequentialdense3BiasAddReadVariableOp0bcast = std::vector<float>(1);
float * tensor_sequentialdense3BiasAddReadVariableOp0bcast = fTensor_sequentialdense3BiasAddReadVariableOp0bcast.data();
std::vector<float> fTensor_sequentialdense2Sigmoid0 = std::vector<float>(6);
float * tensor_sequentialdense2Sigmoid0 = fTensor_sequentialdense2Sigmoid0.data();
std::vector<float> fTensor_sequentialdense2MatMulGemm80 = std::vector<float>(6);
float * tensor_sequentialdense2MatMulGemm80 = fTensor_sequentialdense2MatMulGemm80.data();
std::vector<float> fTensor_sequentialdense2BiasAddReadVariableOp0bcast = std::vector<float>(6);
float * tensor_sequentialdense2BiasAddReadVariableOp0bcast = fTensor_sequentialdense2BiasAddReadVariableOp0bcast.data();
std::vector<float> fTensor_dense3 = std::vector<float>(1);
float * tensor_dense3 = fTensor_dense3.data();
std::vector<float> fTensor_sequentialdense1MatMulGemm70 = std::vector<float>(7);
float * tensor_sequentialdense1MatMulGemm70 = fTensor_sequentialdense1MatMulGemm70.data();
std::vector<float> fTensor_sequentialdense1Sigmoid0 = std::vector<float>(7);
float * tensor_sequentialdense1Sigmoid0 = fTensor_sequentialdense1Sigmoid0.data();
std::vector<float> fTensor_sequentialdense1BiasAddReadVariableOp0bcast = std::vector<float>(7);
float * tensor_sequentialdense1BiasAddReadVariableOp0bcast = fTensor_sequentialdense1BiasAddReadVariableOp0bcast.data();
std::vector<float> fTensor_sequentialdenseSigmoid0 = std::vector<float>(7);
float * tensor_sequentialdenseSigmoid0 = fTensor_sequentialdenseSigmoid0.data();
std::vector<float> fTensor_sequentialdenseMatMulGemm60 = std::vector<float>(7);
float * tensor_sequentialdenseMatMulGemm60 = fTensor_sequentialdenseMatMulGemm60.data();
std::vector<float> fTensor_sequentialdense3MatMulGemm90 = std::vector<float>(1);
float * tensor_sequentialdense3MatMulGemm90 = fTensor_sequentialdense3MatMulGemm90.data();
std::vector<float> fTensor_sequentialdenseBiasAddReadVariableOp0bcast = std::vector<float>(7);
float * tensor_sequentialdenseBiasAddReadVariableOp0bcast = fTensor_sequentialdenseBiasAddReadVariableOp0bcast.data();


Session(std::string filename ="TrkQual_ANN1_v2.dat") {

//--- reading weights from file
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()) {
      throw std::runtime_error("tmva-sofie failed to open file " + filename + " for input weights");
   }
   std::string tensor_name;
   size_t length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdenseBiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdenseBiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 7) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 7 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdenseBiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdense1MatMulReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdense1MatMulReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 49) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 49 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdense1MatMulReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdense1BiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdense1BiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 7) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 7 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdense1BiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdense2BiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdense2BiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 6) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 6 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdense2BiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdense2MatMulReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdense2MatMulReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 42) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 42 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdense2MatMulReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdense3BiasAddReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdense3BiasAddReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdense3BiasAddReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdenseMatMulReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdenseMatMulReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 49) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 49 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdenseMatMulReadVariableOp0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_sequentialdense3MatMulReadVariableOp0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_sequentialdense3MatMulReadVariableOp0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 6) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 6 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_sequentialdense3MatMulReadVariableOp0[i];
   f.close();

//---- allocate the intermediate dynamic tensors
//--- broadcast bias tensor sequentialdenseBiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequentialdenseBiasAddReadVariableOp0,{ 7 }, { 1 , 7 });
      std::copy(data, data + 7, tensor_sequentialdenseBiasAddReadVariableOp0bcast);
      delete [] data;
   }
//--- broadcast bias tensor sequentialdense1BiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequentialdense1BiasAddReadVariableOp0,{ 7 }, { 1 , 7 });
      std::copy(data, data + 7, tensor_sequentialdense1BiasAddReadVariableOp0bcast);
      delete [] data;
   }
//--- broadcast bias tensor sequentialdense2BiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequentialdense2BiasAddReadVariableOp0,{ 6 }, { 1 , 6 });
      std::copy(data, data + 6, tensor_sequentialdense2BiasAddReadVariableOp0bcast);
      delete [] data;
   }
//--- broadcast bias tensor sequentialdense3BiasAddReadVariableOp0for Gemm op
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_sequentialdense3BiasAddReadVariableOp0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_sequentialdense3BiasAddReadVariableOp0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 7;
   int op_0_k = 7;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 7;
   int op_0_ldb = 7;
   std::copy(tensor_sequentialdenseBiasAddReadVariableOp0bcast, tensor_sequentialdenseBiasAddReadVariableOp0bcast + 7, tensor_sequentialdenseMatMulGemm60);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_sequentialdenseMatMulReadVariableOp0, &op_0_ldb, tensor_input, &op_0_lda, &op_0_beta, tensor_sequentialdenseMatMulGemm60, &op_0_n);
	for (int id = 0; id < 7 ; id++){
		tensor_sequentialdenseSigmoid0[id] = 1 / (1 + std::exp( - tensor_sequentialdenseMatMulGemm60[id]));
	}

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 1;
   int op_2_n = 7;
   int op_2_k = 7;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 7;
   int op_2_ldb = 7;
   std::copy(tensor_sequentialdense1BiasAddReadVariableOp0bcast, tensor_sequentialdense1BiasAddReadVariableOp0bcast + 7, tensor_sequentialdense1MatMulGemm70);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_sequentialdense1MatMulReadVariableOp0, &op_2_ldb, tensor_sequentialdenseSigmoid0, &op_2_lda, &op_2_beta, tensor_sequentialdense1MatMulGemm70, &op_2_n);
	for (int id = 0; id < 7 ; id++){
		tensor_sequentialdense1Sigmoid0[id] = 1 / (1 + std::exp( - tensor_sequentialdense1MatMulGemm70[id]));
	}

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 1;
   int op_4_n = 6;
   int op_4_k = 7;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 7;
   int op_4_ldb = 6;
   std::copy(tensor_sequentialdense2BiasAddReadVariableOp0bcast, tensor_sequentialdense2BiasAddReadVariableOp0bcast + 6, tensor_sequentialdense2MatMulGemm80);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_sequentialdense2MatMulReadVariableOp0, &op_4_ldb, tensor_sequentialdense1Sigmoid0, &op_4_lda, &op_4_beta, tensor_sequentialdense2MatMulGemm80, &op_4_n);
	for (int id = 0; id < 6 ; id++){
		tensor_sequentialdense2Sigmoid0[id] = 1 / (1 + std::exp( - tensor_sequentialdense2MatMulGemm80[id]));
	}

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 1;
   int op_6_n = 1;
   int op_6_k = 6;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 6;
   int op_6_ldb = 1;
   std::copy(tensor_sequentialdense3BiasAddReadVariableOp0bcast, tensor_sequentialdense3BiasAddReadVariableOp0bcast + 1, tensor_sequentialdense3MatMulGemm90);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_sequentialdense3MatMulReadVariableOp0, &op_6_ldb, tensor_sequentialdense2Sigmoid0, &op_6_lda, &op_6_beta, tensor_sequentialdense3MatMulGemm90, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense3[id] = 1 / (1 + std::exp( - tensor_sequentialdense3MatMulGemm90[id]));
	}
   return fTensor_dense3;
}
};
} //TMVA_SOFIE_TrkQual_ANN1_v2

#endif  // ROOT_TMVA_SOFIE_TRKQUAL_ANN1_V2
