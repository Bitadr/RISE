############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
############################################################
open_project OMPQR_baseline
set_top omp
add_files omp.h
add_files omp.cpp
add_files -tb OMPQR_baseline/main.cpp
open_solution "solution1"
set_part {xc7vx485tffg1761-2} -tool vivado
create_clock -period 10 -name default
#source "./OMPQR_baseline/solution1/directives.tcl"
csim_design
csynth_design
cosim_design
export_design -format ip_catalog
