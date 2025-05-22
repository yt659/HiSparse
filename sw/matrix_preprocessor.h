#ifndef MATRIX_PREPROCESSOR_H
#define MATRIX_PREPROCESSOR_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>

// Define data structure types
using partition_indptr_t = struct { IDX_T start; PACKED_IDX_T nnz; };
using ch_partition_indptr_t = std::vector<partition_indptr_t>;
using ch_packed_idx_t = std::vector<PACKED_IDX_T>;
using ch_packed_val_t = std::vector<PACKED_VAL_T>;

// Extract a submatrix from CSR matrix with specified row range, ensuring row count is multiple of 128
inline spmv::io::CSRMatrix<VAL_T> extract_csr_rows(
    spmv::io::CSRMatrix<VAL_T> &input, 
    size_t row_start, 
    size_t row_end
) {
    spmv::io::CSRMatrix<VAL_T> result;
    
    // Set basic properties
    size_t actual_rows = row_end - row_start;
    
    // Calculate row count padded to multiple of 128
    const size_t align = 128;  // Row count must be multiple of 128
    size_t padded_rows = ((actual_rows + align - 1) / align) * align;
    
    result.num_rows = padded_rows;
    result.num_cols = input.num_cols;
    
    // Calculate non-zero elements count for submatrix
    size_t nnz_start = input.adj_indptr[row_start];
    size_t nnz_end = input.adj_indptr[row_end];
    size_t nnz_count = nnz_end - nnz_start;
    
    // Allocate memory
    result.adj_indptr.resize(result.num_rows + 1);
    result.adj_indices.resize(nnz_count);
    result.adj_data.resize(nnz_count);
    
    // Copy and adjust indptr for actual rows
    for (size_t i = 0; i <= actual_rows; i++) {
        result.adj_indptr[i] = input.adj_indptr[row_start + i] - nnz_start;
    }
    
    // For padded rows, set their indptr equal to the last actual row's end index
    for (size_t i = actual_rows + 1; i <= padded_rows; i++) {
        result.adj_indptr[i] = result.adj_indptr[actual_rows];
    }
    
    // Copy indices and data
    for (size_t i = 0; i < nnz_count; i++) {
        result.adj_indices[i] = input.adj_indices[nnz_start + i];
        result.adj_data[i] = input.adj_data[nnz_start + i];
    }
    
    return result;
}

// Matrix preprocessing result structure with partition information
struct MatrixPreprocessResult {
    spmv::io::CSRMatrix<VAL_T> mat;                                  // reference matrix
    size_t num_row_partitions;                                       // Total number of row partitions
    size_t num_col_partitions;                                       // Number of column partitions
    size_t num_partitions;                                           // Total number of partitions
    std::vector<std::vector<SPMV_MAT_PKT_T>> channel_packets_fpga0;  // Channel data packets for FPGA-0
    std::vector<std::vector<SPMV_MAT_PKT_T>> channel_packets_fpga1;  // Channel data packets for FPGA-1
    size_t rows_per_ch_in_last_row_part_fpga0;                       // Rows per channel in last partition of FPGA-0
    size_t rows_per_ch_in_last_row_part_fpga1;                       // Rows per channel in last partition of FPGA-1
    size_t row_partition_point;                                      // Row partition point
    size_t fpga0_row_partitions;                                     // Number of row partitions for FPGA-0
    size_t fpga1_row_partitions;                                     // Number of row partitions for FPGA-1
    size_t original_total_rows;                                      // Original total rows (before padding)
    size_t fpga0_actual_rows;                                        // Actual rows processed by FPGA-0 (before padding)
    size_t fpga1_actual_rows;                                        // Actual rows processed by FPGA-1 (before padding)
    size_t fpga0_nnz;                                                // Number of non-zeros processed by FPGA-0
    size_t fpga1_nnz;                                                // Number of non-zeros processed by FPGA-1
};

/**
 * Preprocess sparse matrix for SpMV computation using CSR-level row partitioning based on non-zero element balance
 */
inline MatrixPreprocessResult preprocess_sparse_matrix(
    spmv::io::CSRMatrix<float> &ext_matrix,
    bool skip_empty_rows
) {
    using namespace spmv::io;
    MatrixPreprocessResult result;
    
    // Record original matrix row count
    result.original_total_rows = ext_matrix.num_rows;
    
    // Step 1: Round matrix row count up to multiple of PACK_SIZE * NUM_HBM_CHANNELS
    util_round_csr_matrix_dim<float>(ext_matrix, PACK_SIZE * NUM_HBM_CHANNELS, PACK_SIZE);
    
    // Step 2: Convert float matrix to specified numeric type matrix
    result.mat = csr_matrix_convert_from_float<VAL_T>(ext_matrix);
    
    // Step 3: Determine row partition point based on non-zero element distribution
    size_t total_nnz = result.mat.adj_data.size();
    size_t target_nnz_per_fpga = total_nnz / 2;
    size_t current_nnz = 0;
    size_t row_split = 0;
    
    // Scan matrix rows to find row with cumulative non-zero elements closest to target value
    for (size_t r = 0; r < result.mat.num_rows; r++) {
        size_t row_nnz = result.mat.adj_indptr[r + 1] - result.mat.adj_indptr[r];
        current_nnz += row_nnz;
        
        if (current_nnz >= target_nnz_per_fpga) {
            // Check which of current or previous row is closer to target value
            if (r > 0) {
                size_t prev_nnz = current_nnz - row_nnz;
                if ((target_nnz_per_fpga - prev_nnz) < (current_nnz - target_nnz_per_fpga)) {
                    row_split = r;  // Previous row is closer to target
                } else {
                    row_split = r + 1;  // Current row is closer to target
                }
            } else {
                row_split = r + 1;  // Target reached in first row
            }
            break;
        }
    }
    
    // If no suitable partition point found (rare case), use middle row
    if (row_split == 0 || row_split >= result.mat.num_rows) {
        row_split = result.mat.num_rows / 2;  // Fallback to middle row partition
    }
    
    // Ensure partition point is multiple of PACK_SIZE for data alignment
    row_split = (row_split + PACK_SIZE - 1) / PACK_SIZE * PACK_SIZE;
    if (row_split >= result.mat.num_rows) {
        row_split = result.mat.num_rows - PACK_SIZE;
    }
    
    result.row_partition_point = row_split;
    
    // Calculate non-zero elements for each FPGA
    result.fpga0_nnz = result.mat.adj_indptr[row_split];
    result.fpga1_nnz = total_nnz - result.fpga0_nnz;
    
    // Record actual rows processed by each FPGA (before padding)
    result.fpga0_actual_rows = row_split;
    result.fpga1_actual_rows = result.mat.num_rows - row_split;
    
    // Step 4: Split CSR matrix into two parts
    CSRMatrix<VAL_T> csr_fpga0 = extract_csr_rows(result.mat, 0, row_split);
    CSRMatrix<VAL_T> csr_fpga1 = extract_csr_rows(result.mat, row_split, result.mat.num_rows);
    
    std::cout << "DEBUG: Original matrix: " << result.mat.num_rows << " rows, " 
              << total_nnz << " non-zeros" << std::endl;
    std::cout << "DEBUG: Optimal split at row " << row_split 
              << " (target NNZ per FPGA: " << target_nnz_per_fpga << ")" << std::endl;
    std::cout << "DEBUG: FPGA0 matrix: " << csr_fpga0.num_rows << " rows (padded), " 
              << result.fpga0_actual_rows << " rows (actual), " 
              << result.fpga0_nnz << " non-zeros (" 
              << (result.fpga0_nnz * 100.0 / total_nnz) << "%), starting from row 0" << std::endl;
    std::cout << "DEBUG: FPGA1 matrix: " << csr_fpga1.num_rows << " rows (padded), " 
              << result.fpga1_actual_rows << " rows (actual), " 
              << result.fpga1_nnz << " non-zeros (" 
              << (result.fpga1_nnz * 100.0 / total_nnz) << "%), starting from row " 
              << row_split << std::endl;
    
    // Step 5: Convert each CSR part to CPSR format
    // FPGA-0 part
    CPSRMatrix<PACKED_VAL_T, PACKED_IDX_T, PACK_SIZE> cpsr_fpga0 =
        csr2cpsr<PACKED_VAL_T, PACKED_IDX_T, VAL_T, IDX_T, PACK_SIZE>(
            csr_fpga0, IDX_MARKER, LOGICAL_OB_SIZE, LOGICAL_VB_SIZE, NUM_HBM_CHANNELS, skip_empty_rows);
            
    // FPGA-1 part
    CPSRMatrix<PACKED_VAL_T, PACKED_IDX_T, PACK_SIZE> cpsr_fpga1 =
        csr2cpsr<PACKED_VAL_T, PACKED_IDX_T, VAL_T, IDX_T, PACK_SIZE>(
            csr_fpga1, IDX_MARKER, LOGICAL_OB_SIZE, LOGICAL_VB_SIZE, NUM_HBM_CHANNELS, skip_empty_rows);
    
    // Step 6: Calculate partition information for each FPGA
    result.fpga0_row_partitions = (csr_fpga0.num_rows + LOGICAL_OB_SIZE - 1) / LOGICAL_OB_SIZE;
    result.fpga1_row_partitions = (csr_fpga1.num_rows + LOGICAL_OB_SIZE - 1) / LOGICAL_OB_SIZE;
    
    size_t num_col_partitions = (result.mat.num_cols + LOGICAL_VB_SIZE - 1) / LOGICAL_VB_SIZE;
    
    result.num_row_partitions = result.fpga0_row_partitions + result.fpga1_row_partitions;
    result.num_col_partitions = num_col_partitions;
    result.num_partitions = result.num_row_partitions * num_col_partitions;
    
    std::cout << "DEBUG: FPGA0 row partitions: " << result.fpga0_row_partitions << std::endl;
    std::cout << "DEBUG: FPGA1 row partitions: " << result.fpga1_row_partitions << std::endl;
    
    // Step 7: Prepare channel data packets for FPGA-0
    std::vector<ch_partition_indptr_t> channel_partition_indptr_fpga0(NUM_HBM_CHANNELS);
    std::vector<ch_packed_idx_t> channel_indices_fpga0(NUM_HBM_CHANNELS);
    std::vector<ch_packed_val_t> channel_vals_fpga0(NUM_HBM_CHANNELS);
    result.channel_packets_fpga0.resize(NUM_HBM_CHANNELS);
    
    // Initialize channel partition pointers for FPGA-0
    for (size_t c = 0; c < NUM_HBM_CHANNELS; c++) {
        channel_partition_indptr_fpga0[c].resize(result.fpga0_row_partitions * num_col_partitions);
        channel_partition_indptr_fpga0[c][0].start = 0;
    }
    
    // Process channel data for FPGA-0
    for (size_t pc = 0; pc < NUM_HBM_CHANNELS; pc++) {
        for (size_t j = 0; j < result.fpga0_row_partitions; j++) {
            for (size_t i = 0; i < num_col_partitions; i++) {
                auto indptr_partition = cpsr_fpga0.get_packed_indptr(j, i, pc);
                uint32_t num_packets = *std::max_element(indptr_partition.back().data,
                                                      indptr_partition.back().data + PACK_SIZE);
                
                auto indices_partition = cpsr_fpga0.get_packed_indices(j, i, pc);
                auto vals_partition = cpsr_fpga0.get_packed_data(j, i, pc);
                
                channel_indices_fpga0[pc].insert(channel_indices_fpga0[pc].end(), 
                                             indices_partition.begin(), indices_partition.end());
                channel_vals_fpga0[pc].insert(channel_vals_fpga0[pc].end(), 
                                           vals_partition.begin(), vals_partition.end());
                
                channel_indices_fpga0[pc].resize(
                    channel_partition_indptr_fpga0[pc][j * num_col_partitions + i].start + num_packets);
                channel_vals_fpga0[pc].resize(
                    channel_partition_indptr_fpga0[pc][j * num_col_partitions + i].start + num_packets);
                
                assert(channel_indices_fpga0[pc].size() == channel_vals_fpga0[pc].size());
                
                channel_partition_indptr_fpga0[pc][j * num_col_partitions + i].nnz = indptr_partition.back();
                
                if (!((j == (result.fpga0_row_partitions - 1)) && (i == (num_col_partitions - 1)))) {
                    channel_partition_indptr_fpga0[pc][j * num_col_partitions + i + 1].start =
                        channel_partition_indptr_fpga0[pc][j * num_col_partitions + i].start + num_packets;
                }
            }
        }
        
        // Create data packets for FPGA-0
        size_t fpga0_partitions = result.fpga0_row_partitions * num_col_partitions;
        result.channel_packets_fpga0[pc].resize(fpga0_partitions * 2 + channel_indices_fpga0[pc].size());
        
        // Fill partition pointer information for FPGA-0
        for (size_t ij = 0; ij < fpga0_partitions; ij++) {
            result.channel_packets_fpga0[pc][ij * 2].indices.data[0] = 
                channel_partition_indptr_fpga0[pc][ij].start;
            result.channel_packets_fpga0[pc][ij * 2 + 1].indices = 
                channel_partition_indptr_fpga0[pc][ij].nnz;
        }
        
        // Fill actual matrix data for FPGA-0
        uint32_t offset = fpga0_partitions * 2;
        for (size_t i = 0; i < channel_indices_fpga0[pc].size(); i++) {
            result.channel_packets_fpga0[pc][offset + i].indices = channel_indices_fpga0[pc][i];
            result.channel_packets_fpga0[pc][offset + i].vals = channel_vals_fpga0[pc][i];
        }
    }
    
    // Step 8: Prepare channel data packets for FPGA-1
    std::vector<ch_partition_indptr_t> channel_partition_indptr_fpga1(NUM_HBM_CHANNELS);
    std::vector<ch_packed_idx_t> channel_indices_fpga1(NUM_HBM_CHANNELS);
    std::vector<ch_packed_val_t> channel_vals_fpga1(NUM_HBM_CHANNELS);
    result.channel_packets_fpga1.resize(NUM_HBM_CHANNELS);
    
    // Initialize channel partition pointers for FPGA-1
    for (size_t c = 0; c < NUM_HBM_CHANNELS; c++) {
        channel_partition_indptr_fpga1[c].resize(result.fpga1_row_partitions * num_col_partitions);
        channel_partition_indptr_fpga1[c][0].start = 0;
    }
    
    // Process channel data for FPGA-1
    for (size_t pc = 0; pc < NUM_HBM_CHANNELS; pc++) {
        for (size_t j = 0; j < result.fpga1_row_partitions; j++) {
            for (size_t i = 0; i < num_col_partitions; i++) {
                auto indptr_partition = cpsr_fpga1.get_packed_indptr(j, i, pc);
                uint32_t num_packets = *std::max_element(indptr_partition.back().data,
                                                      indptr_partition.back().data + PACK_SIZE);
                
                auto indices_partition = cpsr_fpga1.get_packed_indices(j, i, pc);
                auto vals_partition = cpsr_fpga1.get_packed_data(j, i, pc);
                
                channel_indices_fpga1[pc].insert(channel_indices_fpga1[pc].end(), 
                                             indices_partition.begin(), indices_partition.end());
                channel_vals_fpga1[pc].insert(channel_vals_fpga1[pc].end(), 
                                           vals_partition.begin(), vals_partition.end());
                
                channel_indices_fpga1[pc].resize(
                    channel_partition_indptr_fpga1[pc][j * num_col_partitions + i].start + num_packets);
                channel_vals_fpga1[pc].resize(
                    channel_partition_indptr_fpga1[pc][j * num_col_partitions + i].start + num_packets);
                
                assert(channel_indices_fpga1[pc].size() == channel_vals_fpga1[pc].size());
                
                channel_partition_indptr_fpga1[pc][j * num_col_partitions + i].nnz = indptr_partition.back();
                
                if (!((j == (result.fpga1_row_partitions - 1)) && (i == (num_col_partitions - 1)))) {
                    channel_partition_indptr_fpga1[pc][j * num_col_partitions + i + 1].start =
                        channel_partition_indptr_fpga1[pc][j * num_col_partitions + i].start + num_packets;
                }
            }
        }
        
        // Create data packets for FPGA-1
        size_t fpga1_partitions = result.fpga1_row_partitions * num_col_partitions;
        result.channel_packets_fpga1[pc].resize(fpga1_partitions * 2 + channel_indices_fpga1[pc].size());
        
        // Fill partition pointer information for FPGA-1
        for (size_t ij = 0; ij < fpga1_partitions; ij++) {
            result.channel_packets_fpga1[pc][ij * 2].indices.data[0] = 
                channel_partition_indptr_fpga1[pc][ij].start;
            result.channel_packets_fpga1[pc][ij * 2 + 1].indices = 
                channel_partition_indptr_fpga1[pc][ij].nnz;
        }
        
        // Fill actual matrix data for FPGA-1
        uint32_t offset = fpga1_partitions * 2;
        for (size_t i = 0; i < channel_indices_fpga1[pc].size(); i++) {
            result.channel_packets_fpga1[pc][offset + i].indices = channel_indices_fpga1[pc][i];
            result.channel_packets_fpga1[pc][offset + i].vals = channel_vals_fpga1[pc][i];
        }
    }
    
    // Step 9: Calculate rows per channel in last partition for each FPGA
    if (csr_fpga0.num_rows % LOGICAL_OB_SIZE == 0) {
        result.rows_per_ch_in_last_row_part_fpga0 = LOGICAL_OB_SIZE / NUM_HBM_CHANNELS;
    } else {
        result.rows_per_ch_in_last_row_part_fpga0 = (csr_fpga0.num_rows % LOGICAL_OB_SIZE) / NUM_HBM_CHANNELS;
    }
    
    if (csr_fpga1.num_rows % LOGICAL_OB_SIZE == 0) {
        result.rows_per_ch_in_last_row_part_fpga1 = LOGICAL_OB_SIZE / NUM_HBM_CHANNELS;
    } else {
        result.rows_per_ch_in_last_row_part_fpga1 = (csr_fpga1.num_rows % LOGICAL_OB_SIZE) / NUM_HBM_CHANNELS;
    }
    
    std::cout << "INFO : Matrix loading/preprocessing complete!" << std::endl;
    
    return result;
}

#endif // MATRIX_PREPROCESSOR_H