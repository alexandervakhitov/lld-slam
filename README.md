# LLD-SLAM
**Authors:** 
[Alexander Vakhitov](https://alexandervakhitov.github.io/), based on the [ORB-SLAM2](https://github.com/raulmur/ORB_SLAM2.git) system by
[Raul Mur-Artal](http://webdiis.unizar.es/~raulmur/), [Juan D. Tardos](http://webdiis.unizar.es/~jdtardos/), [J. M. M. Montiel](http://webdiis.unizar.es/~josemari/) and [Dorian Galvez-Lopez](http://doriangalvez.com/) ([DBoW2](https://github.com/dorian3d/DBoW2))


**03.06.2019** Release of the LLD-SLAM source code

LLD-SLAM is a real-time point+line **Stereo SLAM** library based on a well known ORB-SLAM2 system.
It computes the camera trajectory and a sparse 3D reconstruction consisting of points and lines. 

[![LLD-SLAM](https://img.youtube.com/vi/ntFFiwXIhoA/0.jpg)](https://www.youtube.com/watch?v=ntFFiwXIhoA)

### Related Publications:

[Point+Line LLD SLAM] A. Vakhitov and V. Lempitsky. **Learnable Line Segment Descriptor for Visual SLAM**. *IEEE Access, 2019.*

# 1. License

LLD-SLAM is released under a [GPLv3 license](https://github.com/raulmur/ORB_SLAM2/blob/master/License-gpl.txt). For a list of all code/library dependencies (and associated licenses), please see [Dependencies.md](https://github.com/raulmur/ORB_SLAM2/blob/master/Dependencies.md).

If you use LLD-SLAM in an academic work, please cite:

    @article{murTRO2015,
      title={Learnable Line Segment Descriptor for Visual SLAM},
      author={Vakhitov, Alexander and Lempitsky, Victor},
      journal={IEEE Access},
      year={2019}
     }

# 2. Prerequisites
We have tested the library in **Ubuntu 12.04**, **14.04** and **16.04**, but it should be easy to compile in other platforms. A powerful computer (e.g. i7) will ensure real-time performance and provide more stable and accurate results.

## C++11 or C++0x Compiler
We use the new thread and chrono functionalities of C++11.

## LBDMOD library
We use [LBDMOD](git@github.com:alexandervakhitov/lbdmod.git) for line detection and description. Please download, make, install and insert path to the CMakelists.txt.

## Pangolin
We use [Pangolin](https://github.com/stevenlovegrove/Pangolin) for visualization and user interface. Dowload and install instructions can be found at: https://github.com/stevenlovegrove/Pangolin.

## OpenCV
We use [OpenCV](http://opencv.org) to manipulate images and features. Dowload and install instructions can be found at: http://opencv.org. **Required at leat 2.4.3. Tested with OpenCV 2.4.11 and OpenCV 3.2**.

## Eigen3
Required by g2o (see below). Download and install instructions can be found at: http://eigen.tuxfamily.org. **Required at least 3.1.0**.

## DBoW2 and g2o (Included in Thirdparty folder)
We use modified versions of the [DBoW2](https://github.com/dorian3d/DBoW2) library to perform place recognition and [g2o](https://github.com/RainerKuemmerle/g2o) library to perform non-linear optimizations. Both modified libraries (which are BSD) are included in the *Thirdparty* folder.

## ROS (optional)
We provide some examples to process the live input of a monocular, stereo or RGB-D camera using [ROS](ros.org). Building these examples is optional. In case you want to use ROS, a version Hydro or newer is needed.

# 3. Building LLD-SLAM library and examples

Clone the repository:
```
git clone git@github.com:alexandervakhitov/lld-slam.git LLD_SLAM
```

We provide a script `build.sh` to build the *Thirdparty* libraries and *LLD-SLAM*. Please make sure you have installed all required dependencies (see section 2). 
Please modify the CMakelists.txt: insert a correct path to the LBDMOD library instead of <LBDMOD LIBRARY DIR>. Then execute:
```
cd LLD_SLAM
chmod +x build.sh
./build.sh
```

This will create **libLLD_SLAM.so**  at *lib* folder and the executables **stereo_kitti** and **stereo_euroc** in *Examples* folder.

# 4. Stereo Examples

## KITTI+EuRoC combined Dataset

1. Download the KITTI dataset (grayscale images) from http://www.cvlibs.net/datasets/kitti/eval_odometry.php 

2. Download a sequence (ASL format) from http://projects.asl.ethz.ch/datasets/doku.php?id=kmavvisualinertialdatasets

3. Download the precomputed line detections and LBD descriptors from https://yadi.sk/d/D5QEuced7y5I1w. Modify the `KITTIX_Y.yaml`, 'EuRoC_Y.yaml' files to include the paths to the downloaded and unpacked dataset.

4. Execute one of the following commands depending on the dataset you are going to use. Change `KITTIX_Y.yaml` to 
KITTI00-02_Y.yaml, KITTI03_Y.yaml or KITTI04-12_Y.yaml for Y = {LBD,Empty} for sequence 0 to 2, 3, and 4 to 12 respectively. Change `PATH_TO_DATASET_FOLDER` to the uncompressed dataset folder. Change `SEQUENCE_NUMBER` to 00, 01, 02,.., 11. Substitute $SSS with SEQUENCE_NUMBER in the *.yaml file.
```
./Examples/Stereo/stereo_kitti Vocabulary/ORBvoc.txt Examples/Stereo/KITTIX.yaml PATH_TO_DATASET_FOLDER/dataset/sequences/SEQUENCE_NUMBER
```
```
./Examples/Stereo/stereo_euroc Vocabulary/ORBvoc.txt Examples/Stereo/EuRoC_Y.yaml PATH_TO_SEQUENCE/mav0/cam0/data PATH_TO_SEQUENCE/mav0/cam1/data Examples/Stereo/EuRoC_TimeStamps/SEQUENCE.txt
```
```
./Examples/Stereo/stereo_euroc Vocabulary/ORBvoc.txt Examples/Stereo/EuRoC_Y.yaml PATH_TO_SEQUENCE/cam0/data PATH_TO_SEQUENCE/cam1/data Examples/Stereo/EuRoC_TimeStamps/SEQUENCE.txt
```

## KITTI+EuRoC dataset composition
|No. in combined | No. in KITTI | No. in EuRoC |
|----------------|--------------|--------------|
| 0              | 0            |              |
| 1              | 1            |              |
| 2              | 2            |              |
| 3              | 3            |              |
| 4              | 4            |              |
| 5              |              | 0            |
| 6              |              | 1            |
| 7              |              | 3            |
| 8              |              | 5            |
| 9              |              | 7            |
| 10             | 7            |              |
| 11             |              | 9            |
| 12             | 6            |              |
| 13             | 5            |              |
| 14             | 8            |              |
| 15             | 9            |              |
| 16             | 10           |              |
| 17             |              | 2            |
| 18             |              | 4            |
| 19             |              | 6            |
| 20             |              | 8            |
| 21             |              | 10           |
  
  
# 5. Processing your own sequences
You will need to create a settings file with the calibration of your camera. 
See the settings file provided for the KITTI and EuRoC datasets for stereo cameras. We use the calibration model of OpenCV. See the examples to learn how to create a program that makes use of the LLD-SLAM library and how to pass images to the SLAM system. 
Stereo input must be synchronized and rectified. Use the LBD descriptor and detector available in the LBDMOD library, as we show in the [lld-public](git@github.com:alexandervakhitov/lld-public.git). 

# 6. SLAM and Localization Modes
You can change between the *SLAM* and *Localization mode* using the GUI of the map viewer.

### SLAM Mode
This is the default mode. The system runs in parallal three threads: Tracking, Local Mapping and Loop Closing. The system localizes the camera, builds new map and tries to close loops.

### Localization Mode
This mode can be used when you have a good map of your working area. In this mode the Local Mapping and Loop Closing are deactivated. The system localizes the camera in the map (which is no longer updated), using relocalization if needed. 

