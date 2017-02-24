# TFT_vs_Fund
Pose estimation of 3 views based on the trifocal tensor and the fundamental matrix.

This MATLAB(R) directory contains the code implemented for the ICIP 2017 submission "A Critical Review of the Trifocal Tensor Estimation" by Laura F. Julia and Pascal Monasse.

## CONTENT
This directory should contain the following files and folders:

 Files                 | Description
 :-------------------- | :---------------------------------------------------------------
 TFT_methods/          | Pose estimation methods based on the TFT and auxiliary functions
 F_methods/            | Pose estimation methods based on the F and auxiliary functions
 Optimization/         | Optimization functions: Gauss_Helmert.m and BundleAdjustment.m
 auxiliar_functions/   | Auxiliary functuions for triangulation, error computation, etc
 README.md             | This file
 LICENSE               | License file
 experiments.m         | experiments script
 example.m             | example script

## SETUP
1. You must have MATLAB software installed on your computer.
2. Copy/move the 'TFT_vs_Fund' folder to the MATLAB 'work' directory.
3. Open MATLAB and add the 'TFT_vs_Fund' directory to the Path.
4. Run one of the example scripts detailed in usage section below.

## USAGE
Two scripts are provided to easily use the code: 
example.m     - Gives an example on how to use the code for pose estimation of three views.
experiments.m - Script to recreate the experiments with synthetic data
