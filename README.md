# AMATH 482: Computational Methods for Data Processing 

This class covered exploratory and objective data analysis methods applied to the physical, engineering, and biological sciences. Brief review of statistical methods and their computational implementation for studying time series analysis, spectral analysis, filtering methods, principal component analysis, orthogonal mode decomposition, and image processing and compression.

## Repository Description: 
This repository contains the homework projects for AMATH 482. Each folder contains: the homework specification, the final report, and code for the project. 

#### HW 1: Fourier Analysis
Given a noisy Ultrasound dataset, the goal of this problem is to determine the frequency
signature of a marble, its trajectory, and its location at any point in time.
The original noisy data is in the spatial domain. Using Fourier transforms to convert the data
into the frequency domain, the data is then averaged to find the frequency signature (center
frequency) of the marble. Next, a Gaussian filter is used to remove noise. Finally, by using an
inverse Fourier Transform to convert the data back into the spatial domain, we have a clean
data set to observe the marble

#### HW 2: Gabor Transforms
This report has two parts. Both are focused on time-frequency analysis using Gabor Transforms
to build spectrograms. Part 1, we will analyze a music sample and explore how changing
different properties of the Gabor Transform effects the spectrogram. Additionally, we will see
how a Mexican hat wavelet filter effects the spectrogram.
In Part 2, we will analyze music samples of two different instruments playing the same song.
The goal is to determine what notes are being played, their order, and the relative duration of
each note, and filter out overtones for each music sample.

#### HW 3: Principal Component Analysis
In this report is about using PCA on different datasets. The goal is to illustrate the usefulness
and effects of noise on the PCA algorithm. A mass is attached to a spring and is bouncing up and
down. There are four different tests of the mass on the spring. Given video taken from multiple
different perspectives, for each test we will extract the position of the mass and then use PCA
to analyze the dynamics of this system. By extracting the energies from each dataset, we should
be able to tell what dimension the motion of the mass is in.

#### HW 4: Machine Learning
The objective of this homework is to explore the basics of machine learning by writing code that
can classify a song by what genre it is using a five second sample clip. To program the classifier,
we will use spectrograms, Principal Component Analysis, (PCA) and Linear Discriminant Analysis
(LDA). 

#### HW 5: Neural Networks
This homework is about exploring neural networks. We will use two different types of neural
networks to build a classifier for the MINST Fashion dataset. This dataset contains images of 10
different classes of fashion items. In part one we will build a classifier using a fully-connected
neural network, and in part two we will use convolutional neural network. Then we will use the
results from each to compare their accuracy.

