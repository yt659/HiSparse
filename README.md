# HiSparse

High-Performance Sparse Linear Algebra on HBM-Equipped FPGAs Using HLS: A Case Study on SpMV

**Note**: This is an updated version that supports newer Vitis toolchains (2023.2) and includes compatibility fixes for modern GCC compilers.

## Overview

HiSparse is a high-performance accelerator for sparse-matrix vector multiplication (SpMV) implemented on multi-die HBM-equipped FPGA devices. The design achieves 237MHz and delivers promising speedup with increased bandwidth efficiency compared to prior implementations on CPUs, GPUs, and FPGAs.

For detailed technical information, please refer to our [FPGA 2022 paper](https://github.com/cornell-zhang/HiSparse/blob/master/fpgafp193a-du.pdf).

## Citation

```bibtex
@article{du2022hisparse,
title={{High-Performance Sparse Linear Algebra on HBM-Equipped FPGAs Using HLS: A Case Study on SpMV}},
author={Du, Yixiao and Hu, Yuwei and Zhou, Zhongchun and Zhang, Zhiru},
journal={{Int'l Symp. on Field-Programmable Gate Arrays (FPGA)}},
year={2022}
}
```

## Platform Requirements

- **FPGA Platform**: Xilinx Alveo U280
- **Toolchain**: Xilinx Vitis 2023.2 (updated from original 2020.2)
- **OS**: Linux (tested on CentOS/RHEL)
- **Compiler**: GCC 11+ (for compatibility with modern C++ standards)

## Quick Start

### 1. Clone the Repository

```bash
git clone -b 2021+ https://github.com/yt659/HiSparse.git
cd HiSparse
```

**Important**: Use the `2021+` branch which includes compatibility updates for newer Vitis versions.

### 2. Compile Local CNPY Library

Due to GCC version compatibility issues, you need to compile a local version of the CNPY library:

```bash
./compile_cnpy.sh
```

This script will:
- Create a `local_cnpy` directory
- Copy CNPY source files
- Compile CNPY with your current GCC version

### 3. Setup Environment

```bash
source setup.sh
```

This will configure:
- Vitis 2023.2 environment variables
- HLS include paths
- Local CNPY library paths
- Required library paths

### 4. Download Datasets

```bash
cd datasets
source download.sh
cd ..
```

This downloads two dataset categories:
- `graph`: Graph-based sparse matrices
- `pruned_nn`: Pruned neural network matrices

### 5. Build and Run Demo

```bash
cd sw
make demo
```

The demo will run with a pre-compiled bitstream and output performance metrics in the format:
```
{Preprocessing: 0.64566 s | SpMV: 0.77102 ms | 49.4087 GBPS | 12.9698 GOPS }
```

Where:
- **Preprocessing**: Data preparation time
- **SpMV**: Sparse matrix-vector multiplication runtime  
- **GBPS**: Data throughput (Gigabytes per second)
- **GOPS**: Operation throughput (Giga operations per second)

Note: Data throughput = Operation throughput / 2 * 8

## Advanced Usage

### Building Host Application

```bash
cd sw
make host
```

### Hardware Emulation

```bash
cd sw
make hw_emu
```

### Hardware Implementation

```bash
cd sw
make hw
```

### Running Benchmarks

```bash
cd sw
make benchmark IMPL=<fixed/float_pob/float_stall>
```

Implementation options:
- `fixed`: Fixed-point design (default)
- `float_pob`: Floating-point design using partial output buffers
- `float_stall`: Floating-point design using stall + row interleaving

## Troubleshooting

### Environment Setup Issues

If you encounter Vitis setup errors:

1. **Check Vitis Installation**:
   ```bash
   ls -la /opt/xilinx/Vitis/
   ```

2. **Verify HLS Include Path**:
   ```bash
   echo $HLS_INCLUDE
   ls $HLS_INCLUDE/ap_fixed.h
   ```

3. **Check CNPY Library**:
   ```bash
   echo $CNPY_LIB
   ls $CNPY_LIB/libcnpy.so
   ```

### Compilation Issues

If you get linking errors related to CNPY:
1. Ensure you've run `./compile_cnpy.sh`
2. Re-source the setup script: `source setup.sh`
3. Clean and rebuild: `make cleanall && make host`

### Runtime Library Issues

If you encounter `GLIBCXX` version errors:
- The setup script automatically configures library paths to avoid conflicts
- If issues persist, check that `LD_LIBRARY_PATH` prioritizes system libraries

## Project Structure

```
HiSparse/
├── spmv/                   # Fixed-point SpMV implementation
├── spmv-fp/               # Floating-point SpMV implementation  
├── sw/                    # Host software and drivers
├── datasets/              # Benchmark datasets
├── local_cnpy/           # Local CNPY library (generated)
├── xrt/                  # XRT utilities
├── demo_spmv.xclbin      # Pre-compiled demo bitstream
├── setup.sh              # Environment setup script
└── compile_cnpy.sh       # CNPY compilation script
```

## Key Improvements in This Version

1. **Updated Toolchain Support**: Compatible with Vitis 2023.2
2. **GCC Compatibility**: Resolves C++ ABI issues with modern compilers
3. **Automated Setup**: Simplified environment configuration
4. **Local CNPY Build**: Eliminates version conflicts with system libraries
5. **Improved Documentation**: Clear setup and troubleshooting instructions

## Original Repository

This is a fork of the original [cornell-zhang/HiSparse](https://github.com/cornell-zhang/HiSparse) repository with compatibility updates for modern development environments.

## Support

For technical issues:
1. Check the troubleshooting section above
2. Verify your environment meets the platform requirements  
3. Open an issue with detailed error messages and system information

## License

Please refer to the original repository for licensing information.