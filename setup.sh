echo setting up Vitis-2023.2 environment ...
unset LM_LICENSE_FILE
export XILINXD_LICENSE_FILE=2100@flex.ece.cornell.edu
source scl_source enable devtoolset-8
source /opt/xilinx/Vitis/2023.2/settings64.sh > /dev/null
source /opt/xilinx/xrt/setup.sh > /dev/null
export HLS_INCLUDE=/opt/xilinx/Vitis_HLS/2023.2/include
export LD_LIBRARY_PATH=/opt/xilinx/Vitis/2023.2/lib/lnx64.o:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/xilinx/Vitis/2023.2/lib/lnx64.o/Default:$LD_LIBRARY_PATH
echo Vitis-2023.2 setup finished

echo setting up svpp ...
export PATH=$PWD/svpp:$PATH
echo svpp setup finished

echo setting up cnpy ...
# Always use local cnpy compiled with current GCC version
# Get the HiSparse root directory (parent of current sw directory)
HISPARSE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CNPY_INCLUDE=$HISPARSE_ROOT/local_cnpy
export CNPY_LIB=$HISPARSE_ROOT/local_cnpy
# Set simplified LD_LIBRARY_PATH to avoid libstdc++ version conflicts
export LD_LIBRARY_PATH=$HISPARSE_ROOT/local_cnpy:/opt/xilinx/xrt/lib
echo "Using local cnpy library: $CNPY_LIB/libcnpy.so"
echo cnpy setup finished

echo setup finished:
echo " v++ is: $(which v++)"
echo " xbutil is: $(which xbutil)"
echo " svpp is: $(which svpp)"
echo " libcnpy.so is: $CNPY_LIB/libcnpy.so"