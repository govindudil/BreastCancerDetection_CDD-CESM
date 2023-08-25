//
// Created by govindu on 5/26/23.
//

// asm - Angular Second Moment
// contrast - Contrast
// corr - Correlation
// var - Sum of Squares: Variance
// idm - Inverse Difference Moment
// savg - Sum Average
// svar - Sum Variance
// sentropy - Sum Entropy
// entropy - Entropy
// dvar - Difference Variance
// dentropy - Difference Entropy
// icorr1 - Information Measures of Correlation
// icorr2 - Information Measures of Correlation
// maxcorr - Maximal Correlation Coefficient

#ifndef CANCERDETECTION_GLCMMAIN_H
#define CANCERDETECTION_GLCMMAIN_H

#include <cstdio>
#include <zconf.h>
#include "stdio.h"
#include "opencv2/opencv.hpp"

using namespace cv;
using namespace std;

//Initialize functions that have been used for measuring co-occurrence matrices for 0,45,90,135 degree angle
double** CoOcMat_Angle_0   (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_45  (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_90  (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_135 (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);


//Initializing functions used for quantifying co-occurance matrixes
double f1_asm (double **P, int Ng);
double f2_contrast (double **P, int Ng);
double f3_corr (double **P, int Ng);
double f4_var (double **P, int Ng);
double f5_idm (double **P, int Ng);
double f6_savg (double **P, int Ng);
double f7_svar (double **P, int Ng, double S);
double f8_sentropy (double **P, int Ng);
double f9_entropy (double **P, int Ng);
double f10_dvar (double **P, int Ng);
double f11_dentropy (double **P, int Ng);
double f12_icorr (double **P, int Ng);
double f13_icorr (double **P, int Ng);
double f14_maxcorr (double **P, int Ng);

//Supporting matrix for above calculation
double *allocate_vector (int nl, int nh);
double **allocate_matrix (int nrl, int nrh, int ncl, int nch);
void free_matrix(double **matrix,int nrh);

tuple<__decay_and_strip<double **&>::__type, int> grayco_matrcis(Mat img, int angle, int distance);

double  *calculator(Mat img,char* param[],int length,int angle,int distance);

#endif //CANCERDETECTION_GLCMMAIN_H
