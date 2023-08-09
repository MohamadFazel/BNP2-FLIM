# BNP2-FLIM
Bayesian nonparametric analysis of fluorescence lifetime imaging (BNP-FLIM) is a wide-filed FLIM data analysis software package developed based on the framework described in “Building Fluorescence Lifetime Maps Photon-by-photon by Leveraging Spatial Correlations” by Mohamadreza Fazel, Sina Jazani, Lorenzo Scipioni, Alexander Vallmitjana, Songning Zhu, Enrico Gratton, Michelle A. Digman and Steve Presse; BioRxiv. 

BNP-FLIM is capable of deducing the number of lifetime components present within an input FLIM data and learn the corresponding lifetimes maps with resolutions better than the pixel size leveraging spatial correlations, over a wide range of lifetimes and lifetimes with small differences. Moreover, BNP-FLIM takes full advantages of all the available data using direct photon arrival times and empty pulses. 

Running the example script:
1) It requires MATLAB 2020a or higher versions.
2) Download the software package and open the "runFLIM_Example.m" script. All the parameters are described in the scripts.
3) Add path to the "Functions" directory
4) Run the script.
5) The scripts will load "Data_20Ph_per_pixel.mat" containing lifetimes of 1ns 1.8ns and 4.5ns simulated over an area of 20x5 pixels.
6) It takes ~3 hours to run this example on a desktop with AMD Ryzen 9 3900X 3.8 GHz CPU.

To analyze any new data set, it has to be first saved in the format compatible with BNP2-FLIM. That is a structure array where the number of elements are the same as the pixels. Each element of the array containt information of a corresponding pixel, including: pixel center (mu); photon arrival times collected from the pixel (ns); a binary vector with the same length as the number of laser pulses used to illuminate the pixel where each element represent a pulse, zero and one, respectively, stand for an empty or occupied pulse.  

