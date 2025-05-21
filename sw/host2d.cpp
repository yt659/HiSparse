//--------------------------------------------------------------------------------------------------
// Part1 Header File Inclusions and Macro Definitions
//--------------------------------------------------------------------------------------------------
#include "common.h"
#include "data_loader.h"
#include "data_formatter.h"
#include "matrix_preprocessor.h"

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <thread>
#include <chrono>
#include <string>

#include "xcl2.hpp"

// device memory channels
#define MAX_HBM_CHANNEL_COUNT 32
#define CHANNEL_NAME(n) n | XCL_MEM_TOPOLOGY
const int HBM[MAX_HBM_CHANNEL_COUNT] = {
    CHANNEL_NAME(0),  CHANNEL_NAME(1),  CHANNEL_NAME(2),  CHANNEL_NAME(3),  CHANNEL_NAME(4),
    CHANNEL_NAME(5),  CHANNEL_NAME(6),  CHANNEL_NAME(7),  CHANNEL_NAME(8),  CHANNEL_NAME(9),
    CHANNEL_NAME(10), CHANNEL_NAME(11), CHANNEL_NAME(12), CHANNEL_NAME(13), CHANNEL_NAME(14),
    CHANNEL_NAME(15), CHANNEL_NAME(16), CHANNEL_NAME(17), CHANNEL_NAME(18), CHANNEL_NAME(19),
    CHANNEL_NAME(20), CHANNEL_NAME(21), CHANNEL_NAME(22), CHANNEL_NAME(23), CHANNEL_NAME(24),
    CHANNEL_NAME(25), CHANNEL_NAME(26), CHANNEL_NAME(27), CHANNEL_NAME(28), CHANNEL_NAME(29),
    CHANNEL_NAME(30), CHANNEL_NAME(31)};

const int DDR[2] = { CHANNEL_NAME(32), CHANNEL_NAME(33) };

template<typename T>
using aligned_vector = std::vector<T, aligned_allocator<T> >;

//--------------------------------------------------------------------------------------------------
// Part2 Reference and Verify Utils
//--------------------------------------------------------------------------------------------------
void compute_ref(spmv::io::CSRMatrix<float> &mat,
                 std::vector<float> &vector,
                 std::vector<float> &ref_result) {
    ref_result.resize(mat.num_rows);
    std::fill(ref_result.begin(), ref_result.end(), 0);
    for (size_t row_idx = 0; row_idx < mat.num_rows; row_idx++) {
        IDX_T start = mat.adj_indptr[row_idx];
        IDX_T end = mat.adj_indptr[row_idx + 1];
        for (size_t i = start; i < end; i++) {
            IDX_T idx = mat.adj_indices[i];
            ref_result[row_idx] += mat.adj_data[i] * vector[idx];
        }
    }
}

bool verify(std::vector<float> reference_results,
            std::vector<VAL_T> kernel_results) {
    float epsilon = 0.0001;
    if (reference_results.size() != kernel_results.size()) {
        std::cout << "Error: Size mismatch" << std::endl;
        std::cout << "  Reference result size: " << reference_results.size()
                  << "  Kernel result size: " << kernel_results.size() << std::endl;
        return false;
    }
    for (size_t i = 0; i < reference_results.size(); i++) {
        bool match = (std::abs(float(kernel_results[i]) - reference_results[i]) < epsilon);
        if (!match) {
            std::cout << "Error: Result mismatch" << std::endl;
            std::cout << "  i = " << i
                      << "  Reference result = " << reference_results[i]
                      << "  Kernel result = " << kernel_results[i] << std::endl;
            return false;
        }
    }
    return true;
}

void unpack_vector(aligned_vector<PACKED_VAL_T> &pdv, std::vector<VAL_T> &dv) {
    dv.resize(pdv.size() * PACK_SIZE);
    for (size_t i = 0; i < pdv.size(); i++) {
        for (size_t k = 0; k < PACK_SIZE; k++) {
            dv[i * PACK_SIZE + k] = pdv[i].data[k];
        }
    }
}

//---------------------------------------------------------------
// Part3 Test Harness Utils
//---------------------------------------------------------------
#define CL_CREATE_EXT_PTR(name, data, channel)  \
    cl_mem_ext_ptr_t name;                      \
    name.obj = data;                            \
    name.param = 0;                             \
    name.flags = channel;

#define CL_BUFFER_RDONLY(context, size, ext, err)               \
    cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR, size, &ext, &err);

#define CL_BUFFER_WRONLY(context, size, ext, err)               \
    cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR, size, &ext, &err);

#define CL_BUFFER(context, size, ext, err)                      \
    cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR, size, &ext, &err);

#define CHECK_ERR(err)                                          \
    if (err != CL_SUCCESS) {                                    \
        printf("OCL Error at %s:%d, error code is: %d\n",       \
               __FILE__, __LINE__, err);                        \
        exit(EXIT_FAILURE);                                     \
    }

struct cl_runtime {
    cl::Context context;
    cl::CommandQueue command_queue;
    cl::Kernel spmv_sk0;
    cl::Kernel spmv_sk1;
    cl::Kernel spmv_sk2;
    cl::Kernel vector_loader;
    cl::Kernel result_drain;
};

//---------------------------------------------------------------
// Part4 Dual FPGA Test Harness
//---------------------------------------------------------------
bool spmv_test_harness_dual(cl_runtime &runtime0,
                            cl_runtime &runtime1,
                            spmv::io::CSRMatrix<float> &ext_matrix,
                            bool skip_empty_rows) {
    using namespace spmv::io;
    std::cout << "INFO : Test started (Dual FPGA)" << std::endl;
    
    // Use matrix preprocessor to preprocess matrix with CSR-level row cutting
    MatrixPreprocessResult prep_result = preprocess_sparse_matrix(ext_matrix, skip_empty_rows);
    
    // Use preprocessing results
    spmv::io::CSRMatrix<VAL_T> &mat = prep_result.mat;
    size_t num_row_partitions = prep_result.num_row_partitions;
    size_t num_col_partitions = prep_result.num_col_partitions;
    size_t num_partitions = prep_result.num_partitions;
    
    // Convert channel data packages to aligned_vector format
    std::vector<aligned_vector<SPMV_MAT_PKT_T>> channel_packets_fpga0(NUM_HBM_CHANNELS);
    std::vector<aligned_vector<SPMV_MAT_PKT_T>> channel_packets_fpga1(NUM_HBM_CHANNELS);
    
    for (size_t c = 0; c < NUM_HBM_CHANNELS; c++) {
        // Convert FPGA-0 packets
        channel_packets_fpga0[c].resize(prep_result.channel_packets_fpga0[c].size());
        for (size_t i = 0; i < prep_result.channel_packets_fpga0[c].size(); i++) {
            channel_packets_fpga0[c][i] = prep_result.channel_packets_fpga0[c][i];
        }
        
        // Convert FPGA-1 packets
        channel_packets_fpga1[c].resize(prep_result.channel_packets_fpga1[c].size());
        for (size_t i = 0; i < prep_result.channel_packets_fpga1[c].size(); i++) {
            channel_packets_fpga1[c][i] = prep_result.channel_packets_fpga1[c][i];
        }
    }

    //-------------------------------------------------------------------------
    // Generate input vector and create two copies for two FPGAs
    //-------------------------------------------------------------------------
    std::vector<float> vector_f(ext_matrix.num_cols);
    std::generate(vector_f.begin(), vector_f.end(), [](){ return float(rand() % 2); });
    aligned_vector<PACKED_VAL_T> vector(mat.num_cols / PACK_SIZE);
    for (size_t i = 0; i < vector.size(); i++) {
        for (size_t k = 0; k < PACK_SIZE; k++) {
            vector[i].data[k] = VAL_T(vector_f[i * PACK_SIZE + k]);
        }
    }
    aligned_vector<PACKED_VAL_T> vector_fpga0 = vector;
    aligned_vector<PACKED_VAL_T> vector_fpga1 = vector;

    //-------------------------------------------------------------------------
    // allocate space for results (full result buffer for each device)
    //-------------------------------------------------------------------------
    aligned_vector<PACKED_VAL_T> result_full0(mat.num_rows / PACK_SIZE);
    aligned_vector<PACKED_VAL_T> result_full1(mat.num_rows / PACK_SIZE);
    for (size_t i = 0; i < result_full0.size(); i++) {
        for (size_t k = 0; k < PACK_SIZE; k++) {
            result_full0[i].data[k] = 0;
            result_full1[i].data[k] = 0;
        }
    }
    std::cout << "INFO : Input/result initialization complete!" << std::endl;

    //-------------------------------------------------------------------------
    // allocate memory on FPGA for both devices and move data
    //-------------------------------------------------------------------------
    cl_int err;
    std::vector<cl::Buffer> channel_packets_buf0(NUM_HBM_CHANNELS);
    std::vector<cl::Buffer> channel_packets_buf1(NUM_HBM_CHANNELS);
    cl_mem_ext_ptr_t channel_packets_ext0[NUM_HBM_CHANNELS];
    cl_mem_ext_ptr_t channel_packets_ext1[NUM_HBM_CHANNELS];
    for (size_t c = 0; c < NUM_HBM_CHANNELS; c++) {
        // Setup memory for FPGA-0
        channel_packets_ext0[c].obj = channel_packets_fpga0[c].data();
        channel_packets_ext0[c].param = 0;
        channel_packets_ext0[c].flags = HBM[c];
        size_t channel_packets_size0 = sizeof(SPMV_MAT_PKT_T) * channel_packets_fpga0[c].size();
        if (channel_packets_size0 >= 256 * 1000 * 1000) {
            std::cout << "Error: Trying to allocate " << channel_packets_size0/1000/1000
                      << " MB on HBM channel " << c << " for FPGA-0" << std::endl
                      << ", but the capacity of one HBM channel is 256 MB." << std::endl;
            exit(EXIT_FAILURE);
        }
        channel_packets_buf0[c] =
            CL_BUFFER_RDONLY(runtime0.context, channel_packets_size0, channel_packets_ext0[c], err);
        CHECK_ERR(err);
        
        // Setup memory for FPGA-1
        channel_packets_ext1[c].obj = channel_packets_fpga1[c].data();
        channel_packets_ext1[c].param = 0;
        channel_packets_ext1[c].flags = HBM[c];
        size_t channel_packets_size1 = sizeof(SPMV_MAT_PKT_T) * channel_packets_fpga1[c].size();
        if (channel_packets_size1 >= 256 * 1000 * 1000) {
            std::cout << "Error: Trying to allocate " << channel_packets_size1/1000/1000
                      << " MB on HBM channel " << c << " for FPGA-1" << std::endl
                      << ", but the capacity of one HBM channel is 256 MB." << std::endl;
            exit(EXIT_FAILURE);
        }
        channel_packets_buf1[c] =
            CL_BUFFER_RDONLY(runtime1.context, channel_packets_size1, channel_packets_ext1[c], err);
        CHECK_ERR(err);
    }

    // According to hardware requirements, both devices use HBM[20] to store vectors and HBM[21] to store results
    CL_CREATE_EXT_PTR(vector_ext0, vector_fpga0.data(), HBM[20]);
    CL_CREATE_EXT_PTR(vector_ext1, vector_fpga1.data(), HBM[20]);
    CL_CREATE_EXT_PTR(result_ext0, result_full0.data(), HBM[21]);
    CL_CREATE_EXT_PTR(result_ext1, result_full1.data(), HBM[21]);
    size_t vector_size = sizeof(VAL_T) * mat.num_cols;
    size_t result_size = sizeof(VAL_T) * mat.num_rows;
    cl::Buffer vector_buf0 = CL_BUFFER_RDONLY(runtime0.context, vector_size, vector_ext0, err);
    CHECK_ERR(err);
    cl::Buffer vector_buf1 = CL_BUFFER_RDONLY(runtime1.context, vector_size, vector_ext1, err);
    CHECK_ERR(err);
    cl::Buffer result_buf0 = CL_BUFFER_WRONLY(runtime0.context, result_size, result_ext0, err);
    CHECK_ERR(err);
    cl::Buffer result_buf1 = CL_BUFFER_WRONLY(runtime1.context, result_size, result_ext1, err);
    CHECK_ERR(err);

    for (size_t c = 0; c < NUM_HBM_CHANNELS; c++) {
        OCL_CHECK(err, err = runtime0.command_queue.enqueueMigrateMemObjects({channel_packets_buf0[c]}, 0));
        OCL_CHECK(err, err = runtime1.command_queue.enqueueMigrateMemObjects({channel_packets_buf1[c]}, 0));
    }
    OCL_CHECK(err, err = runtime0.command_queue.enqueueMigrateMemObjects({vector_buf0}, 0));
    OCL_CHECK(err, err = runtime1.command_queue.enqueueMigrateMemObjects({vector_buf1}, 0));
    runtime0.command_queue.finish();
    runtime1.command_queue.finish();
    std::cout << "INFO : Host -> Device data transfer complete!" << std::endl;

    //-------------------------------------------------------------------------
    // Set fixed kernel arguments for both FPGAs
    //-------------------------------------------------------------------------
    for (size_t c = 0; c < SK0_CLUSTER; c++) {
        OCL_CHECK(err, err = runtime0.spmv_sk0.setArg(c, channel_packets_buf0[c]));
    }
    for (size_t c = 0; c < SK1_CLUSTER; c++) {
        OCL_CHECK(err, err = runtime0.spmv_sk1.setArg(c, channel_packets_buf0[c + SK0_CLUSTER]));
    }
    for (size_t c = 0; c < SK2_CLUSTER; c++) {
        OCL_CHECK(err, err = runtime0.spmv_sk2.setArg(c, channel_packets_buf0[c + SK0_CLUSTER + SK1_CLUSTER]));
    }
    OCL_CHECK(err, err = runtime0.spmv_sk0.setArg(SK0_CLUSTER + 4, (unsigned)num_col_partitions));
    OCL_CHECK(err, err = runtime0.spmv_sk0.setArg(SK0_CLUSTER + 5, (unsigned)(prep_result.fpga0_row_partitions * num_col_partitions)));
    OCL_CHECK(err, err = runtime0.spmv_sk1.setArg(SK1_CLUSTER + 4, (unsigned)num_col_partitions));
    OCL_CHECK(err, err = runtime0.spmv_sk1.setArg(SK1_CLUSTER + 5, (unsigned)(prep_result.fpga0_row_partitions * num_col_partitions)));
    OCL_CHECK(err, err = runtime0.spmv_sk2.setArg(SK2_CLUSTER + 4, (unsigned)num_col_partitions));
    OCL_CHECK(err, err = runtime0.spmv_sk2.setArg(SK2_CLUSTER + 5, (unsigned)(prep_result.fpga0_row_partitions * num_col_partitions)));
    OCL_CHECK(err, err = runtime0.vector_loader.setArg(0, vector_buf0));
    OCL_CHECK(err, err = runtime0.vector_loader.setArg(1, (unsigned)mat.num_cols));
    OCL_CHECK(err, err = runtime0.result_drain.setArg(0, result_buf0));

    for (size_t c = 0; c < SK0_CLUSTER; c++) {
        OCL_CHECK(err, err = runtime1.spmv_sk0.setArg(c, channel_packets_buf1[c]));
    }
    for (size_t c = 0; c < SK1_CLUSTER; c++) {
        OCL_CHECK(err, err = runtime1.spmv_sk1.setArg(c, channel_packets_buf1[c + SK0_CLUSTER]));
    }
    for (size_t c = 0; c < SK2_CLUSTER; c++) {
        OCL_CHECK(err, err = runtime1.spmv_sk2.setArg(c, channel_packets_buf1[c + SK0_CLUSTER + SK1_CLUSTER]));
    }
    OCL_CHECK(err, err = runtime1.spmv_sk0.setArg(SK0_CLUSTER + 4, (unsigned)num_col_partitions));
    OCL_CHECK(err, err = runtime1.spmv_sk0.setArg(SK0_CLUSTER + 5, (unsigned)(prep_result.fpga1_row_partitions * num_col_partitions)));
    OCL_CHECK(err, err = runtime1.spmv_sk1.setArg(SK1_CLUSTER + 4, (unsigned)num_col_partitions));
    OCL_CHECK(err, err = runtime1.spmv_sk1.setArg(SK1_CLUSTER + 5, (unsigned)(prep_result.fpga1_row_partitions * num_col_partitions)));
    OCL_CHECK(err, err = runtime1.spmv_sk2.setArg(SK2_CLUSTER + 4, (unsigned)num_col_partitions));
    OCL_CHECK(err, err = runtime1.spmv_sk2.setArg(SK2_CLUSTER + 5, (unsigned)(prep_result.fpga1_row_partitions * num_col_partitions)));
    OCL_CHECK(err, err = runtime1.vector_loader.setArg(0, vector_buf1));
    OCL_CHECK(err, err = runtime1.vector_loader.setArg(1, (unsigned)mat.num_cols));
    OCL_CHECK(err, err = runtime1.result_drain.setArg(0, result_buf1));

    //-------------------------------------------------------------------------
    // Submit tasks for both FPGAs simultaneously (MODIFIED PART)
    //-------------------------------------------------------------------------
    // First, queue all tasks for FPGA-0 without waiting for completion
    for (size_t row_part_id = 0; row_part_id < prep_result.fpga0_row_partitions; row_part_id++) {
        unsigned part_len = LOGICAL_OB_SIZE / NUM_HBM_CHANNELS;
        if (row_part_id == prep_result.fpga0_row_partitions - 1) {
            part_len = prep_result.rows_per_ch_in_last_row_part_fpga0;
        }
        
        std::cout << "DEBUG: Queuing FPGA-0 row partition " << row_part_id 
                  << ", part_len = " << part_len << std::endl;
                  
        OCL_CHECK(err, err = runtime0.spmv_sk0.setArg(SK0_CLUSTER + 2, (unsigned)row_part_id));
        OCL_CHECK(err, err = runtime0.spmv_sk0.setArg(SK0_CLUSTER + 3, (unsigned)part_len));
        OCL_CHECK(err, err = runtime0.spmv_sk1.setArg(SK1_CLUSTER + 2, (unsigned)row_part_id));
        OCL_CHECK(err, err = runtime0.spmv_sk1.setArg(SK1_CLUSTER + 3, (unsigned)part_len));
        OCL_CHECK(err, err = runtime0.spmv_sk2.setArg(SK2_CLUSTER + 2, (unsigned)row_part_id));
        OCL_CHECK(err, err = runtime0.spmv_sk2.setArg(SK2_CLUSTER + 3, (unsigned)part_len));
        OCL_CHECK(err, err = runtime0.result_drain.setArg(1, (unsigned)row_part_id));
        
        OCL_CHECK(err, err = runtime0.command_queue.enqueueTask(runtime0.vector_loader));
        OCL_CHECK(err, err = runtime0.command_queue.enqueueTask(runtime0.spmv_sk0));
        OCL_CHECK(err, err = runtime0.command_queue.enqueueTask(runtime0.spmv_sk1));
        OCL_CHECK(err, err = runtime0.command_queue.enqueueTask(runtime0.spmv_sk2));
        OCL_CHECK(err, err = runtime0.command_queue.enqueueTask(runtime0.result_drain));
        
        // No finish() here - we continue to queue tasks without waiting
    }
    
    // Then, queue all tasks for FPGA-1 without waiting for completion
    for (size_t row_part_id = 0; row_part_id < prep_result.fpga1_row_partitions; row_part_id++) {
        unsigned part_len = LOGICAL_OB_SIZE / NUM_HBM_CHANNELS;
        if (row_part_id == prep_result.fpga1_row_partitions - 1) {
            part_len = prep_result.rows_per_ch_in_last_row_part_fpga1;
        }
        
        std::cout << "DEBUG: Queuing FPGA-1 row partition " << row_part_id 
                  << ", part_len = " << part_len << std::endl;
                  
        OCL_CHECK(err, err = runtime1.spmv_sk0.setArg(SK0_CLUSTER + 2, (unsigned)row_part_id));
        OCL_CHECK(err, err = runtime1.spmv_sk0.setArg(SK0_CLUSTER + 3, (unsigned)part_len));
        OCL_CHECK(err, err = runtime1.spmv_sk1.setArg(SK1_CLUSTER + 2, (unsigned)row_part_id));
        OCL_CHECK(err, err = runtime1.spmv_sk1.setArg(SK1_CLUSTER + 3, (unsigned)part_len));
        OCL_CHECK(err, err = runtime1.spmv_sk2.setArg(SK2_CLUSTER + 2, (unsigned)row_part_id));
        OCL_CHECK(err, err = runtime1.spmv_sk2.setArg(SK2_CLUSTER + 3, (unsigned)part_len));
        OCL_CHECK(err, err = runtime1.result_drain.setArg(1, (unsigned)row_part_id));
        
        OCL_CHECK(err, err = runtime1.command_queue.enqueueTask(runtime1.vector_loader));
        OCL_CHECK(err, err = runtime1.command_queue.enqueueTask(runtime1.spmv_sk0));
        OCL_CHECK(err, err = runtime1.command_queue.enqueueTask(runtime1.spmv_sk1));
        OCL_CHECK(err, err = runtime1.command_queue.enqueueTask(runtime1.spmv_sk2));
        OCL_CHECK(err, err = runtime1.command_queue.enqueueTask(runtime1.result_drain));
        
        // No finish() here - continue queuing tasks
    }
    
    // Now wait for all tasks on both FPGAs to complete
    std::cout << "DEBUG: Waiting for both FPGAs to complete all tasks..." << std::endl;
    runtime0.command_queue.finish();
    runtime1.command_queue.finish();
    
    std::cout << "INFO : SpMV kernel complete on both devices!" << std::endl;

    //-------------------------------------------------------------------------
    // Compute reference results
    //-------------------------------------------------------------------------
    std::vector<float> ref_result;
    compute_ref(ext_matrix, vector_f, ref_result);
    std::cout << "INFO : Compute reference complete!" << std::endl;

    //-------------------------------------------------------------------------
    // Transfer device results back to Host
    //-------------------------------------------------------------------------
    OCL_CHECK(err, err = runtime0.command_queue.enqueueMigrateMemObjects({result_buf0}, CL_MIGRATE_MEM_OBJECT_HOST));
    runtime0.command_queue.finish();
    
    OCL_CHECK(err, err = runtime1.command_queue.enqueueMigrateMemObjects({result_buf1}, CL_MIGRATE_MEM_OBJECT_HOST));
    runtime1.command_queue.finish();
    
    std::cout << "INFO : Device -> Host data transfer complete!" << std::endl;

    //-------------------------------------------------------------------------
    // Merge results - for row partitioning, we need to keep the first half of results from FPGA-0, second half from FPGA-1
    //-------------------------------------------------------------------------
    aligned_vector<PACKED_VAL_T> result_merged(mat.num_rows / PACK_SIZE);
    
    // Row partition point (aligned to PACK_SIZE)
    size_t row_split_packed = (prep_result.row_partition_point + PACK_SIZE - 1) / PACK_SIZE;
    
    // Copy FPGA-0 results (first part of rows)
    for (size_t i = 0; i < row_split_packed; i++) {
        for (size_t k = 0; k < PACK_SIZE; k++) {
            result_merged[i].data[k] = result_full0[i].data[k];
        }
    }
    
    // Copy FPGA-1 results (second part of rows)
    for (size_t i = row_split_packed; i < result_merged.size(); i++) {
        for (size_t k = 0; k < PACK_SIZE; k++) {
            // Note: FPGA-1 results start at index 0, need to adjust
            size_t fpga1_idx = i - row_split_packed;
            if (fpga1_idx < result_full1.size()) {
                result_merged[i].data[k] = result_full1[fpga1_idx].data[k];
            }
        }
    }
    
    std::vector<VAL_T> upk_result;
    unpack_vector(result_merged, upk_result);
    return verify(ref_result, upk_result);
}
//---------------------------------------------------------------
// Part5 test case utils
//---------------------------------------------------------------
spmv::io::CSRMatrix<float> create_dense_CSR(unsigned num_rows, unsigned num_cols) {
    spmv::io::CSRMatrix<float> mat_f;
    mat_f.num_rows = num_rows;
    mat_f.num_cols = num_cols;
    mat_f.adj_data.resize(num_rows * num_cols);
    mat_f.adj_indices.resize(num_rows * num_cols);
    mat_f.adj_indptr.resize(num_rows + 1);
    for (auto &x : mat_f.adj_data) { x = 1; }
    for (size_t i = 0; i < num_rows; i++) {
        for (size_t j = 0; j < num_cols; j++) {
            mat_f.adj_indices[i * num_cols + j] = j;
        }
    }
    for (size_t i = 0; i < num_rows + 1; i++) {
        mat_f.adj_indptr[i] = num_cols * i;
    }
    return mat_f;
}

spmv::io::CSRMatrix<float> create_uniform_sparse_CSR(unsigned num_rows, unsigned num_cols, unsigned nnz_per_row) {
    spmv::io::CSRMatrix<float> mat_f;
    mat_f.num_rows = num_rows;
    mat_f.num_cols = num_cols;
    mat_f.adj_data.resize(num_rows * nnz_per_row);
    mat_f.adj_indices.resize(num_rows * nnz_per_row);
    mat_f.adj_indptr.resize(num_rows + 1);
    for (auto &x : mat_f.adj_data) { x = 1; }
    unsigned indice_step = num_cols / nnz_per_row;
    for (size_t i = 0; i < num_rows; i++) {
        for (size_t j = 0; j < nnz_per_row; j++) {
            mat_f.adj_indices[i * nnz_per_row + j] = (indice_step * j + i) % num_cols;
        }
    }
    for (size_t i = 0; i < num_rows + 1; i++) {
        mat_f.adj_indptr[i] = nnz_per_row * i;
    }
    return mat_f;
}

//---------------------------------------------------------------
// test cases (Dual FPGA)
//---------------------------------------------------------------
std::string GRAPH_DATASET_DIR = "../datasets/graph/";
std::string NN_DATASET_DIR = "../datasets/pruned_nn/";

bool test_basic(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on basic dense matrix (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = create_dense_CSR(128, 128);
    for (auto &x : mat_f.adj_data) {x = 1;}
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, false)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

bool test_basic_sparse(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on basic sparse matrix (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = create_uniform_sparse_CSR(1000, 1024, 10);    
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, false)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

bool test_medium_sparse(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on uniform 10K 10 (10K, 1M) (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = create_uniform_sparse_CSR(10000, 10000, 10);
    for (auto &x : mat_f.adj_data) { x = 1; }
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, false)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

bool test_gplus(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on google_plus (108K, 13M) (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = spmv::io::load_csr_matrix_from_float_npz(GRAPH_DATASET_DIR + "gplus_108K_13M_csr_float32.npz");
    for (auto &x : mat_f.adj_data) { x = 1 / mat_f.num_cols; }
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, false)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

bool test_ogbl_ppa(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on ogbl_ppa (576K, 42M) (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = spmv::io::load_csr_matrix_from_float_npz(GRAPH_DATASET_DIR + "ogbl_ppa_576K_42M_csr_float32.npz");
    for (auto &x : mat_f.adj_data) { x = 1 / mat_f.num_cols; }
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, false)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

bool test_transformer_50_t(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on transformer-50-t (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = spmv::io::load_csr_matrix_from_float_npz(NN_DATASET_DIR + "transformer_50_512_33288_csr_float32.npz");
    for (auto &x : mat_f.adj_data) { x = 1 / mat_f.num_cols; }
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, true)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

bool test_transformer_95_t(cl_runtime &runtime0, cl_runtime &runtime1) {
    std::cout << "------ Running test: on transformer-95-t (dual FPGA)" << std::endl;
    spmv::io::CSRMatrix<float> mat_f = spmv::io::load_csr_matrix_from_float_npz(NN_DATASET_DIR + "transformer_95_512_33288_csr_float32.npz");
    for (auto &x : mat_f.adj_data) { x = 1 / mat_f.num_cols; }
    if (spmv_test_harness_dual(runtime0, runtime1, mat_f, true)) {
        std::cout << "INFO : Testcase passed." << std::endl;
        return true;
    } else {
        std::cout << "INFO : Testcase failed." << std::endl;
        return false;
    }
}

//---------------------------------------------------------------
// Part6 main
//---------------------------------------------------------------
int main (int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <hw_emu/hw> <xclbin>" << std::endl;
        return 0;
    }
    std::string target = argv[1];
    std::string xclbin = argv[2];
    if (target != "hw_emu" && target != "hw") {
        std::cout << "This host program only supports hw_emu and hw!" << std::endl;
        return 1;
    }
    cl_runtime runtime0, runtime1;
    cl_int err;
    if (target == "sw_emu" || target == "hw_emu") {
        setenv("XCL_EMULATION_MODE", target.c_str(), true);
    }
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    std::vector<cl::Device> target_devices;
    for (size_t i = 0; i < devices.size(); i++) {
        if (devices[i].getInfo<CL_DEVICE_NAME>() == "xilinx_u280_gen3x16_xdma_base_1") {
            target_devices.push_back(devices[i]);
            if (target_devices.size() == 2)
                break;
        }
    }
    if (target_devices.size() < 2) {
        std::cout << "ERROR : Failed to find two devices: xilinx_u280_gen3x16_xdma_base_1" << std::endl;
        exit(EXIT_FAILURE);
    }
    runtime0.context = cl::Context(target_devices[1], NULL, NULL, NULL, &err);
    CHECK_ERR(err);
    runtime1.context = cl::Context(target_devices[0], NULL, NULL, NULL, &err);
    CHECK_ERR(err);
    auto file_buf = xcl::read_binary_file(xclbin);
    cl::Program::Binaries binaries{{file_buf.data(), file_buf.size()}};
    cl::Program program0(runtime0.context, {target_devices[1]}, binaries, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cout << "ERROR : Failed to program device0 with xclbin file" << std::endl;
        return 1;
    } else {
        std::cout << "INFO : Successfully programmed device0 with xclbin file" << std::endl;
    }
    cl::Program program1(runtime1.context, {target_devices[0]}, binaries, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cout << "ERROR : Failed to program device1 with xclbin file" << std::endl;
        return 1;
    } else {
        std::cout << "INFO : Successfully programmed device1 with xclbin file" << std::endl;
    }
    OCL_CHECK(err, runtime0.spmv_sk0 = cl::Kernel(program0, "spmv_sk0", &err));
    OCL_CHECK(err, runtime0.spmv_sk1 = cl::Kernel(program0, "spmv_sk1", &err));
    OCL_CHECK(err, runtime0.spmv_sk2 = cl::Kernel(program0, "spmv_sk2", &err));
    OCL_CHECK(err, runtime0.vector_loader = cl::Kernel(program0, "spmv_vector_loader", &err));
    OCL_CHECK(err, runtime0.result_drain = cl::Kernel(program0, "spmv_result_drain", &err));
    OCL_CHECK(err, runtime1.spmv_sk0 = cl::Kernel(program1, "spmv_sk0", &err));
    OCL_CHECK(err, runtime1.spmv_sk1 = cl::Kernel(program1, "spmv_sk1", &err));
    OCL_CHECK(err, runtime1.spmv_sk2 = cl::Kernel(program1, "spmv_sk2", &err));
    OCL_CHECK(err, runtime1.vector_loader = cl::Kernel(program1, "spmv_vector_loader", &err));
    OCL_CHECK(err, runtime1.result_drain = cl::Kernel(program1, "spmv_result_drain", &err));

    // Create a command queue that supports out-of-order (two FPGAs execute in parallel at the same time)
    OCL_CHECK(err, runtime0.command_queue = cl::CommandQueue(
        runtime0.context, target_devices[1],
        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err));
    OCL_CHECK(err, runtime1.command_queue = cl::CommandQueue(
        runtime1.context, target_devices[0],
        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err));

    bool passed = true;
    passed = passed && test_basic(runtime0, runtime1);
    passed = passed && test_basic_sparse(runtime0, runtime1);
    passed = passed && test_medium_sparse(runtime0, runtime1);
    if (target != "hw_emu") {
        passed = passed && test_gplus(runtime0, runtime1);
        passed = passed && test_ogbl_ppa(runtime0, runtime1);
        passed = passed && test_transformer_50_t(runtime0, runtime1);
    }
    passed = passed && test_transformer_95_t(runtime0, runtime1);
    std::cout << (passed ? "===== All Test Passed! =====" : "===== Test FAILED! =====") << std::endl;
    return passed ? 0 : 1;
}