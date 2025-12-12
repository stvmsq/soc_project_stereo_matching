# System on Project Project - Stereo Matching
This project aims to develop a reproducible and modular stereo matching pipeline on a Zedboard as a foundation for future research in partial reconfiguration (PR) and hardware-accelerated image processing.
The system implements the Semi-Global Matching (SGM) algorithm for depth estimation from stereo image pairs, both as a pure software solution on the ARM processor and as a hardware-accelerated variant with a custom cost computation module in the programmable logic (PL).
A dedicated test platform will evaluate different configurations using standard datasets such as KITTI 2012/2015 and Middlebury, providing a reliable baseline for performance, accuracy, and resource efficiency in reconfigurable stereo vision systems.

## Project structure
- **HostScript_Server**  
  Contains the Python-based test platform executed on the host PC.  
  A `client.py` script is included, which can simulate a client.
- **SemiGlobalMatching**  
  Contains a working C implementation of the Semi-Global Matching (SGM) algorithm.
- **ZedBoard/Vivado**  
  Contains the hardware design of the ZedBoard project.
- **ZedBoard/Vitis**  
  Contains the software for the ZedBoard, including  
  - the communication interface to the host PC, and  
  - the actual SGM implementation.
- **Documentation**  
  Contains the project report.

## Objective
The aim of this project is to develop a reproducible and modular stereo matching pipeline on a Zedboard with the Zynq SoC, which will serve as a foundation for future work in the field of partial reconfiguration (PR) and hardware-accelerated image processing.
The focus is not on developing new algorithms, but on establishing a solid reference platform that allows different approaches, parameters, and hardware/software combinations to be objectively compared.

## Application
The application performs depth estimation from stereo image pairs based on the well-known Semi-Global Matching (SGM) algorithm.
SGM determines the disparity i.e., the horizontal shift between the left and right images—for each pixel by aggregating local matching costs across multiple directions, thereby achieving a balance between accuracy and computational effort.
It is considered a highly efficient compromise between simple block-matching methods and complex global optimization techniques such as belief propagation.

## Task
### Stereo Matching
The SGM algorithm will be fully implemented on the ARM processor of the Zynq SoC.
The system receives a stereo image pair as input and outputs a corresponding disparity or depth map.
In addition to the pure software implementation, a hardware-accelerated variant will also be developed.
A custom instruction or dedicated hardware kernel in the programmable logic (PL) will perform the computationally intensive cost calculation.

### Test Platform
A modular test platform will be developed to automatically evaluate each implementation variant.
It will handle input data management, execute the processing, and analyze the resulting output against the ground truth data.

The following metrics will be recorded by the platform:
1. Bad Pixel Rate (BPR)
2. Root Mean Square Error (RMSE)
3. Frame Rate (FPS)
4. Latency per Frame
5. Resource Utilization (LUTs, FFs, DSPs, BRAM)
Optionally, dynamic and static power consumption can be estimated or, if possible, measured for each approach.

## Discussion
The results obtained from the test platform for the different implementation approaches will be compared and analyzed.
Based on these results, the advantages and disadvantages of the approaches will be evaluated.
In addition, suggestions will be made on how partial reconfiguration can be integrated into the system.

## Evaluation
The different implementations will be evaluated using the developed test platform and the following datasets:
 - KITTI 2012 / 2015
 - Middlebury
Optionally, the DrivingStereo datasets can also be used for additional testing.
