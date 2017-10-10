# RISE

This repository provides the source code for RISE framework (http://www.aceslab.org/sites/default/files/RISE.pdf). RISE can be used as a generic framework Real-Time Intelligent Video Surveillance on FPGA.

### Program Required
* Vivado HLS

## File Hierarchy of Repository

### OMPQR_baseline:
Implementation of Orthogonal Matching Pursuit (OMP) algorithm with no specific hardware optimization

### OMPQR_optimized:
Optimized implementation of OMP algorithm with specific memory management, pipelining, etc. to improve the throughput of the underlying system.

![ScreenShot](https://raw.github.com/Bitadr/RISE/master/opt.PNG)
  
## Code Description
**RISE_CodeDescription.pdf** describes the source code and different optimizations in detail.
