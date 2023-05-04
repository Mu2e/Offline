//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Thu Apr 27 17:54:50 2023] 

#ifndef TMVA_SOFIE_TRAINBKGDIAG
#define TMVA_SOFIE_TRAINBKGDIAG

#include<algorithm>
#include<vector>
#include<cmath>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_TrainBkgDiag{
namespace BLAS{
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_dense7bias0 = std::vector<float>(1);
float * tensor_dense7bias0 = fTensor_dense7bias0.data();
std::vector<float> fTensor_dense7kernel0 = std::vector<float>(22);
float * tensor_dense7kernel0 = fTensor_dense7kernel0.data();
std::vector<float> fTensor_dense6bias0 = std::vector<float>(22);
float * tensor_dense6bias0 = fTensor_dense6bias0.data();
std::vector<float> fTensor_dense4kernel0 = std::vector<float>(242);
float * tensor_dense4kernel0 = fTensor_dense4kernel0.data();
std::vector<float> fTensor_dense4bias0 = std::vector<float>(22);
float * tensor_dense4bias0 = fTensor_dense4bias0.data();
std::vector<float> fTensor_dense5kernel0 = std::vector<float>(484);
float * tensor_dense5kernel0 = fTensor_dense5kernel0.data();
std::vector<float> fTensor_dense5bias0 = std::vector<float>(22);
float * tensor_dense5bias0 = fTensor_dense5bias0.data();
std::vector<float> fTensor_dense6kernel0 = std::vector<float>(484);
float * tensor_dense6kernel0 = fTensor_dense6kernel0.data();
std::vector<float> fTensor_dense7Sigmoid0 = std::vector<float>(32);
float * tensor_dense7Sigmoid0 = fTensor_dense7Sigmoid0.data();
std::vector<float> fTensor_dense4bias0bcast = std::vector<float>(704);
float * tensor_dense4bias0bcast = fTensor_dense4bias0bcast.data();
std::vector<float> fTensor_dense4Relu0 = std::vector<float>(704);
float * tensor_dense4Relu0 = fTensor_dense4Relu0.data();
std::vector<float> fTensor_dense5Dense = std::vector<float>(704);
float * tensor_dense5Dense = fTensor_dense5Dense.data();
std::vector<float> fTensor_dense6Dense = std::vector<float>(704);
float * tensor_dense6Dense = fTensor_dense6Dense.data();
std::vector<float> fTensor_dense4Dense = std::vector<float>(704);
float * tensor_dense4Dense = fTensor_dense4Dense.data();
std::vector<float> fTensor_dense5bias0bcast = std::vector<float>(704);
float * tensor_dense5bias0bcast = fTensor_dense5bias0bcast.data();
std::vector<float> fTensor_dense5Relu0 = std::vector<float>(704);
float * tensor_dense5Relu0 = fTensor_dense5Relu0.data();
std::vector<float> fTensor_dense6bias0bcast = std::vector<float>(704);
float * tensor_dense6bias0bcast = fTensor_dense6bias0bcast.data();
std::vector<float> fTensor_dense7Dense = std::vector<float>(32);
float * tensor_dense7Dense = fTensor_dense7Dense.data();
std::vector<float> fTensor_dense6Relu0 = std::vector<float>(704);
float * tensor_dense6Relu0 = fTensor_dense6Relu0.data();
std::vector<float> fTensor_dense7bias0bcast = std::vector<float>(32);
float * tensor_dense7bias0bcast = fTensor_dense7bias0bcast.data();


Session(std::string filename ="") {
   if (filename.empty()) filename = "TrainBkgDiag.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense7bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense7bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense7bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense7kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense7kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 22) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 22 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense7kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense6bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense6bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 22) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 22 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense6bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense4kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense4kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 242) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 242 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense4kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense4bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense4bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 22) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 22 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense4bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense5kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense5kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 484) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 484 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense5kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense5bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense5bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 22) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 22 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense5bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense6kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense6kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 484) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 484 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense6kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense4bias0,{ 22 }, { 32 , 22 });
      std::copy(data, data + 704, tensor_dense4bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense5bias0,{ 22 }, { 32 , 22 });
      std::copy(data, data + 704, tensor_dense5bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense6bias0,{ 22 }, { 32 , 22 });
      std::copy(data, data + 704, tensor_dense6bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense7bias0,{ 1 }, { 32 , 1 });
      std::copy(data, data + 32, tensor_dense7bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input2){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 32;
   int op_0_n = 22;
   int op_0_k = 11;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 11;
   int op_0_ldb = 22;
   std::copy(tensor_dense4bias0bcast, tensor_dense4bias0bcast + 704, tensor_dense4Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense4kernel0, &op_0_ldb, tensor_input2, &op_0_lda, &op_0_beta, tensor_dense4Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 704 ; id++){
      tensor_dense4Relu0[id] = ((tensor_dense4Dense[id] > 0 )? tensor_dense4Dense[id] : 0);
   }

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 32;
   int op_2_n = 22;
   int op_2_k = 22;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 22;
   int op_2_ldb = 22;
   std::copy(tensor_dense5bias0bcast, tensor_dense5bias0bcast + 704, tensor_dense5Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense5kernel0, &op_2_ldb, tensor_dense4Relu0, &op_2_lda, &op_2_beta, tensor_dense5Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 704 ; id++){
      tensor_dense5Relu0[id] = ((tensor_dense5Dense[id] > 0 )? tensor_dense5Dense[id] : 0);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 32;
   int op_4_n = 22;
   int op_4_k = 22;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 22;
   int op_4_ldb = 22;
   std::copy(tensor_dense6bias0bcast, tensor_dense6bias0bcast + 704, tensor_dense6Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense6kernel0, &op_4_ldb, tensor_dense5Relu0, &op_4_lda, &op_4_beta, tensor_dense6Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 704 ; id++){
      tensor_dense6Relu0[id] = ((tensor_dense6Dense[id] > 0 )? tensor_dense6Dense[id] : 0);
   }

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 32;
   int op_6_n = 1;
   int op_6_k = 22;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 22;
   int op_6_ldb = 1;
   std::copy(tensor_dense7bias0bcast, tensor_dense7bias0bcast + 32, tensor_dense7Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense7kernel0, &op_6_ldb, tensor_dense6Relu0, &op_6_lda, &op_6_beta, tensor_dense7Dense, &op_6_n);
	for (int id = 0; id < 32 ; id++){
		tensor_dense7Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense7Dense[id]));
	}
   std::vector<float> ret (tensor_dense7Sigmoid0, tensor_dense7Sigmoid0 + 32);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // TMVA_SOFIE_TRAINBKGDIAG
