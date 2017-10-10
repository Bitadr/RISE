############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################
open_project OMPQR_final
set_top omp
add_files OMPQR_final/omp.h
add_files OMPQR_final/omp.cpp
add_files -tb OMPQR_final/main.cpp
open_solution "solution1"
set_part {xc7vx690tffg1761-2}
create_clock -period 10 -name default
config_schedule -effort medium -verbose
#source "./OMPQR_final/solution1/directives.tcl"
csim_design
csynth_design
cosim_design -rtl systemc
export_design -format sysgen_ise
