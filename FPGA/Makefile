CPP=g++
XCL_LIBS_DIR=../../../libs

COMMON_COMPILE_OPTIONS=-g -Wall -O0
CPP_COMPILE_OPTIONS=$(COMMON_COMPILE_OPTIONS) -Wall -O3 -std=c++11 -I$(XCL_LIBS_DIR)/xcl2 -DNOT_AWS_FPGA
C_COMPILE_OPTIONS=$(COMMON_COMPILE_OPTIONS) -std=c11
LINK_OPTIONS=-lOpenCL

.PHONY: clean

all:	crackalack_fpga_gen

%.o: %.cpp
	$(CPP) $(CPP_COMPILE_OPTIONS) -o $@ -c $<

crackalack_fpga_gen: crackalack_fpga_gen.o xcl.o xcl2.o
	$(CPP) -o crackalack_fpga_gen crackalack_fpga_gen.o $(LINK_OPTIONS)
xcl.o:
	$(CC) $(C_COMPILE_OPTIONS) -o $@ -c $(XCL_LIBS_DIR)/xcl/xcl.c
xcl2.o:
	$(CPP) $(CPP_COMPILE_OPTIONS) -o $@ -c $(XCL_LIBS_DIR)/xcl2/xcl2.cpp
clean:
	rm -f *~ *.o crackalack_fpga_gen
