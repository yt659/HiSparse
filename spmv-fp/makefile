IMPL := float_pob
TARGET := hw_emu

# see if we can find svpp executable
VPP := $(shell which svpp)
ifeq ($(VPP),)
$(error svpp not found, please setup svpp \(https://github.com/yxd97/svpp\) before running make)
endif

HW_HLS_VERSION = 2020.2
HW_HLS_DEVICE = /work/shared/common/CAD_tool/xilinx/platforms/xilinx_u280_xdma_201920_3/xilinx_u280_xdma_201920_3.xpfm
HW_BACKEND_VERSION = 2022.1
HW_BACKEND_DEVICE = xilinx_u280_gen3x16_xdma_1_202211_1.xpfm

HLS_DIR := ./_x.$(TARGET)
IMP_DIR := ./build_dir.$(TARGET)

VPP += --hls-version $(HW_HLS_VERSION) --backend-version $(HW_BACKEND_VERSION)

HLS_CLFLAGS += -t $(TARGET) --platform $(HW_HLS_DEVICE) --save-temps
HLS_CLFLAGS += -I./libfpga

BACKEND_CLFLAGS += -t $(TARGET) --platform $(HW_BACKEND_DEVICE) --save-temps

ifneq ($(TARGET), hw)
	HLS_CLFLAGS += -g
else
	BACKEND_CLFLAGS += --optimize 3
endif

ifeq ($(IMPL), float_pob)
	HLS_CLFLAGS += -DFP_POB
else
	HLS_CLFLAGS += -DFP_STALL
endif

SPMV_XCLBIN = $(IMP_DIR)/spmv.xclbin

emconfig.json:
	emconfigutil --platform $(DEVICE)

build: $(SPMV_XCLBIN) emconfig.json

BACKEND_CLFLAGS += --config spmv.ini
KERNEL_OBJS += $(HLS_DIR)/spmv_sk0.xo

$(HLS_DIR)/spmv_sk0.xo: spmv_sk0.cpp
	mkdir -p $(HLS_DIR)
	$(VPP) $(HLS_CLFLAGS) --temp_dir $(HLS_DIR) -c -k spmv_sk0 -I $(<D) -o $@ $^

KERNEL_OBJS += $(HLS_DIR)/spmv_sk1.xo

$(HLS_DIR)/spmv_sk1.xo: spmv_sk1.cpp
	mkdir -p $(HLS_DIR)
	$(VPP) $(HLS_CLFLAGS) --temp_dir $(HLS_DIR) -c -k spmv_sk1 -I $(<D) -o $@ $^

KERNEL_OBJS += $(HLS_DIR)/spmv_sk2.xo

$(HLS_DIR)/spmv_sk2.xo: spmv_sk2.cpp
	mkdir -p $(HLS_DIR)
	$(VPP) $(HLS_CLFLAGS) --temp_dir $(HLS_DIR) -c -k spmv_sk2 -I $(<D) -o $@ $^

KERNEL_OBJS += $(HLS_DIR)/spmv_vector_loader.xo

$(HLS_DIR)/spmv_vector_loader.xo: spmv_vector_loader.cpp
	mkdir -p $(HLS_DIR)
	$(VPP) $(HLS_CLFLAGS) --temp_dir $(HLS_DIR) -c -k spmv_vector_loader -I $(<D) -o $@ $^

KERNEL_OBJS += $(HLS_DIR)/spmv_result_drain.xo

$(HLS_DIR)/spmv_result_drain.xo: spmv_result_drain.cpp
	mkdir -p $(HLS_DIR)
	$(VPP) $(HLS_CLFLAGS) --temp_dir $(HLS_DIR) -c -k spmv_result_drain -I $(<D) -o $@ $^

KERNEL_OBJS += $(HLS_DIR)/k2k_relay.xo

$(HLS_DIR)/k2k_relay.xo: k2k_relay.cpp
	mkdir -p $(HLS_DIR)
	$(VPP) $(HLS_CLFLAGS) --temp_dir $(HLS_DIR) -c -k k2k_relay -I $(<D) -o $@ $^

$(SPMV_XCLBIN): $(KERNEL_OBJS)
	mkdir -p $(IMP_DIR)
	$(VPP) $(BACKEND_CLFLAGS) --temp_dir $(IMP_DIR) -l -o $@ $^


.PHONY: clean
clean:
	rm -rf *.log
	rm -rf .run

.PHONY: cleanall
cleanall: clean
	rm -rf emconfig.json
	rm -rf ./_x.sw_emu ./build_dir.sw_emu
	rm -rf ./_x.hw_emu ./build_dir.hw_emu
	rm -rf ./_x.hw ./build_dir.hw
