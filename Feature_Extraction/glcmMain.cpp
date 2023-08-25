//
// Created by govindu on 5/26/23.
//

#include "glcmMain.h"
#include "opencv2/opencv.hpp"
#include <cstdlib>
#include <iostream>
#include "opencv2/core/core.hpp"
#include <opencv2/highgui.hpp>
#include <iomanip>
#include <typeinfo>
#include <time.h>

#define PGM_MAXMAXVAL 255
#define EPSILON 0.000000001
#define RADIX 2.0
#define SIGN(x,y) ((y)<0 ? -fabs(x) : fabs(x))
#define SWAP(a,b) {y=(a);(a)=(b);(b)=y;}

using namespace std;
using namespace cv;
typedef unsigned char u_int8_t;



/* Allocates a double matrix with range [nrl..nrh][ncl..nch] */
double **allocate_matrix (int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (double *));
    if (!m) fprintf (stderr, "memory allocation failure (allocate_matrix 1) "), exit (1);
    m -= ncl;

    /* allocate rows and set pointers to them */
    for (i = nrl; i <= nrh; i++) {
        m[i] = (double *) malloc ((unsigned) (nch - ncl + 1) * sizeof (double));
        if (!m[i]) fprintf (stderr, "memory allocation failure (allocate_matrix 2) "), exit (2);
        m[i] -= ncl;
    }

    /* return pointer to array of pointers to rows */
    return m;
}


double *allocate_vector (int nl, int nh) {
    double *v;

    v = (double *) calloc (1, (unsigned) (nh - nl + 1) * sizeof (double));
    if (!v) fprintf (stderr, "memory allocation failure (allocate_vector) "),  exit (1);

    return v - nl;
}


/* support functions to compute f14_maxcorr */
void mkbalanced (double **a, int n)
{
    int last, j, i;
    double s, r, g, f, c, sqrdx;

    sqrdx = RADIX * RADIX;
    last = 0;
    while (last == 0)
    {
        last = 1;
        for (i = 1; i <= n; i++)
        {
            r = c = 0.0;
            for (j = 1; j <= n; j++)
                if (j != i)
                {
                    c += fabs (a[j][i]);
                    r += fabs (a[i][j]);
                }
            if (c && r)
            {
                g = r / RADIX;
                f = 1.0;
                s = c + r;
                while (c < g)
                {
                    f *= RADIX;
                    c *= sqrdx;
                }
                g = r * RADIX;
                while (c > g)
                {
                    f /= RADIX;
                    c /= sqrdx;
                }
                if ((c + r) / f < 0.95 * s)
                {
                    last = 0;
                    g = 1.0 / f;
                    for (j = 1; j <= n; j++)
                        a[i][j] *= g;
                    for (j = 1; j <= n; j++)
                        a[j][i] *= f;
                }
            }
        }
    }
}


void reduction (double **a, int n)
{
    int m, j, i;
    double y, x;

    for (m = 2; m < n; m++)
    {
        x = 0.0;
        i = m;
        for (j = m; j <= n; j++)
        {
            if (fabs (a[j][m - 1]) > fabs (x))
            {
                x = a[j][m - 1];
                i = j;
            }
        }
        if (i != m)
        {
            for (j = m - 1; j <= n; j++)
            SWAP (a[i][j], a[m][j])
            for (j = 1; j <= n; j++)
            SWAP (a[j][i], a[j][m])
            a[j][i] = a[j][i];
        }
        if (x)
        {
            for (i = m + 1; i <= n; i++)
            {
                if ( (y = a[i][m - 1]) )
                {
                    y /= x;
                    a[i][m - 1] = y;
                    for (j = m; j <= n; j++)
                        a[i][j] -= y * a[m][j];
                    for (j = 1; j <= n; j++)
                        a[j][m] += y * a[j][i];
                }
            }
        }
    }
}


int hessenberg (double **a, int n, double wr[], double wi[])
{
    int nn, m, l, k, j, its, i, mmin;
    double z, y, x, w, v, u, t, s, r=0.0, q=0.0, p=0.0, anorm;

    anorm = fabs (a[1][1]);
    for (i = 2; i <= n; i++)
        for (j = (i - 1); j <= n; j++)
            anorm += fabs (a[i][j]);
    nn = n;
    t = 0.0;
    while (nn >= 1)
    {
        its = 0;
        do
        {
            for (l = nn; l >= 2; l--)
            {
                s = fabs (a[l - 1][l - 1]) + fabs (a[l][l]);
                if (s == 0.0)
                    s = anorm;
                if ((double) (fabs (a[l][l - 1]) + s) == s)
                    break;
            }
            x = a[nn][nn];
            if (l == nn)
            {
                wr[nn] = x + t;
                wi[nn--] = 0.0;
            }
            else
            {
                y = a[nn - 1][nn - 1];
                w = a[nn][nn - 1] * a[nn - 1][nn];
                if (l == (nn - 1))
                {
                    p = 0.5 * (y - x);
                    q = p * p + w;
                    z = sqrt (fabs (q));
                    x += t;
                    if (q >= 0.0)
                    {
                        z = p + SIGN (z, p);
                        wr[nn - 1] = wr[nn] = x + z;
                        if (z)
                            wr[nn] = x - w / z;
                        wi[nn - 1] = wi[nn] = 0.0;
                    }
                    else
                    {
                        wr[nn - 1] = wr[nn] = x + p;
                        wi[nn - 1] = -(wi[nn] = z);
                    }
                    nn -= 2;
                }
                else
                {
                    if (its == 30)
                    {
                        return 0;
                    }
                    if (its == 10 || its == 20)
                    {
                        t += x;
                        for (i = 1; i <= nn; i++)
                            a[i][i] -= x;
                        s = fabs (a[nn][nn - 1]) + fabs (a[nn - 1][nn - 2]);
                        y = x = 0.75 * s;
                        w = -0.4375 * s * s;
                    }
                    ++its;
                    for (m = (nn - 2); m >= l; m--)
                    {
                        z = a[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
                        q = a[m + 1][m + 1] - z - r - s;
                        r = a[m + 2][m + 1];
                        s = fabs (p) + fabs (q) + fabs (r);
                        p /= s;
                        q /= s;
                        r /= s;
                        if (m == l)
                            break;
                        u = fabs (a[m][m - 1]) * (fabs (q) + fabs (r));
                        v = fabs (p) * (fabs (a[m - 1][m - 1]) +
                                        fabs (z) + fabs (a[m + 1][m + 1]));
                        if ((double) (u + v) == v)
                            break;
                    }
                    for (i = m + 2; i <= nn; i++)
                    {
                        a[i][i - 2] = 0.0;
                        if (i != (m + 2))
                            a[i][i - 3] = 0.0;
                    }
                    for (k = m; k <= nn - 1; k++)
                    {
                        if (k != m)
                        {
                            p = a[k][k - 1];
                            q = a[k + 1][k - 1];
                            r = 0.0;
                            if (k != (nn - 1))
                                r = a[k + 2][k - 1];
                            if ( (x = fabs (p) + fabs (q) + fabs (r)) )
                            {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        if ( (s = SIGN (sqrt (p * p + q * q + r * r), p)) )
                        {
                            if (k == m)
                            {
                                if (l != m)
                                    a[k][k - 1] = -a[k][k - 1];
                            }
                            else
                                a[k][k - 1] = -s * x;
                            p += s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q /= p;
                            r /= p;
                            for (j = k; j <= nn; j++)
                            {
                                p = a[k][j] + q * a[k + 1][j];
                                if (k != (nn - 1))
                                {
                                    p += r * a[k + 2][j];
                                    a[k + 2][j] -= p * z;
                                }
                                a[k + 1][j] -= p * y;
                                a[k][j] -= p * x;
                            }
                            mmin = nn < k + 3 ? nn : k + 3;
                            for (i = l; i <= mmin; i++)
                            {
                                p = x * a[i][k] + y * a[i][k + 1];
                                if (k != (nn - 1))
                                {
                                    p += z * a[i][k + 2];
                                    a[i][k + 2] -= p * r;
                                }
                                a[i][k + 1] -= p * q;
                                a[i][k] -= p;
                            }
                        }
                    }
                }
            }
        } while (l < nn - 1);
    }
    return 1;
}

/* 
   This function frees the memory allocated for a double matrix.
   It takes in a double matrix and the number of rows as input parameters.
   It iterates through each column of the matrix and frees the memory allocated for it.
   Finally, it frees the memory allocated for the matrix itself.
*/
void free_matrix(double **matrix,int nrh)
{  
    int col_index;
    for (col_index=0;col_index<=nrh;col_index++)
        free(matrix[col_index]); // free memory allocated for each column
    free(matrix); // free memory allocated for the matrix itself
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// This function calculates the co-occurrence matrix at 0 degree angle for a given distance
// It takes in the distance, grayscale image, number of rows and columns, tone lookup table, and tone count as input parameters
double** CoOcMat_Angle_0 (int distance, u_int8_t **grays,
                          int rows, int cols, int* tone_LUT, int tone_count)
{
    int d = distance;
    int x, y;
    int row, col, itone, jtone;
    double count=0.0; /* normalizing factor */

    // Allocate memory for the co-occurrence matrix
    double** matrix = allocate_matrix (0, tone_count, 0, tone_count);

    // Initialize the co-occurrence matrix with zeros
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            matrix[itone][jtone] = 0.0;

    // Loop through the grayscale image to calculate the co-occurrence matrix
    for (row = 0; row < rows; ++row)
        for (col = 0; col < cols; ++col) {

            // Check if the current column + distance is within the image boundary
            if (col + d < cols) {		// previously stated condition(col + d < cols && grays[row][col + d])
                // Get the tone values for the current pixel and the pixel at a distance of d
                x = tone_LUT[grays[row][col]];
                y = tone_LUT[grays[row][col + d]];
                // Increment the co-occurrence matrix values for the corresponding tone values
                matrix[x][y]++;
                matrix[y][x]++;
                // Increment the normalizing factor
                count += 2.0 ;
            }
        }

    // Normalize the co-occurrence matrix by dividing each element by the normalizing factor
    for (itone = 0; itone < tone_count; ++itone){
        for (jtone = 0; jtone < tone_count; ++jtone){
            if (count==0.0)   /* protect from error */
                matrix[itone][jtone]=0.0;
            else matrix[itone][jtone] /= count;
        }
    }

    // Return the co-occurrence matrix
    return matrix;
}


// This function calculates the co-occurrence matrix at 90 degree angle for a given distance
// It takes in the distance, grayscale image, number of rows and columns, tone lookup table, and tone count as input parameters
double** CoOcMat_Angle_90 (int distance, u_int8_t **grays,
                           int rows, int cols, int* tone_LUT, int tone_count)
{
    int d = distance;
    int x, y;
    int row, col, itone, jtone;
    int count=0; /* normalizing factor */

    // Allocate memory for the co-occurrence matrix
    double** matrix = allocate_matrix (0, tone_count, 0, tone_count);

    // Initialize the co-occurrence matrix with zeros
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            matrix[itone][jtone] = 0;

    // Loop through the grayscale image to calculate the co-occurrence matrix
    for (row = 0; row < rows; ++row)
        for (col = 0; col < cols; ++col) {
            // Check if the current row + distance is within the image boundary
            if (row + d < rows) {
                // Get the tone values for the current pixel and the pixel at a distance of d
                x = tone_LUT [grays[row][col]];
                y = tone_LUT [grays[row + d][col]];
                // Increment the co-occurrence matrix values for the corresponding tone values
                matrix[x][y]++;
                matrix[y][x]++;
                // Increment the normalizing factor
                count += 2 ;
            }
        }

    // Normalize the co-occurrence matrix by dividing each element by the normalizing factor
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            if (count==0) matrix[itone][jtone]=0;
            else matrix[itone][jtone] /= count;

    // Return the co-occurrence matrix
    return matrix;
}


// This function calculates the co-occurrence matrix at 45 degree angle for a given distance
// It takes in the distance, grayscale image, number of rows and columns, tone lookup table, and tone count as input parameters
double** CoOcMat_Angle_45 (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count) {
    int d = distance;
    int x, y;
    int row, col, itone, jtone;
    int count=0; /* normalizing factor */

    // Allocate memory for the co-occurrence matrix
    double** matrix = allocate_matrix (0, tone_count, 0, tone_count);

    // Initialize the co-occurrence matrix with zeros
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            matrix[itone][jtone] = 0;

    // Loop through the grayscale image to calculate the co-occurrence matrix
    for (row = 0; row < rows; ++row)
        for (col = 0; col < cols; ++col) {

            // Check if the current row + distance and column + distance are within the image boundary
            if (row + d < rows && col + d < cols) {
                // Get the tone values for the current pixel and the pixel at a distance of d and 45 degree angle
                x = tone_LUT [grays[row][col]];
                y = tone_LUT [grays[row + d][col + d]];
                // Increment the co-occurrence matrix values for the corresponding tone values
                matrix[x][y]++;
                matrix[y][x]++;
                // Increment the normalizing factor
                count += 2 ;
            }
        }

    // Normalize the co-occurrence matrix by dividing each element by the normalizing factor
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            if (count==0) matrix[itone][jtone]=0;   /* protect from error */
            else matrix[itone][jtone] /= count;

    // Return the co-occurrence matrix
    return matrix;
}


// This function calculates the co-occurrence matrix at 135 degree angle for a given distance
// It takes in the distance, grayscale image, number of rows and columns, tone lookup table, and tone count as input parameters
double** CoOcMat_Angle_135 (int distance, u_int8_t **grays,
                            int rows, int cols, int* tone_LUT, int tone_count)
{
    // Initialize variables
    int d = distance;
    int x, y;
    int row, col, itone, jtone;
    int count=0; /* normalizing factor */

    // Allocate memory for the co-occurrence matrix
    double** matrix = allocate_matrix (0, tone_count, 0, tone_count);

    // Initialize the co-occurrence matrix with zeros
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            matrix[itone][jtone] = 0;

    // Loop through the grayscale image to calculate the co-occurrence matrix
    for (row = 0; row < rows; ++row)
        for (col = 0; col < cols; ++col) {
            // Check if the current row + distance and column - distance are within the image boundary
            if (row + d < rows && col - d >= 0) {
                // Get the tone values for the current pixel and the pixel at a distance of d and 135 degree angle
                x = tone_LUT [grays[row][col]];
                y = tone_LUT [grays[row + d][col - d]];
                // Increment the co-occurrence matrix values for the corresponding tone values
                matrix[x][y]++;
                matrix[y][x]++;
                // Increment the normalizing factor
                count += 2 ;
            }
        }

    // Normalize the co-occurrence matrix by dividing each element by the normalizing factor
    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            if (count==0) matrix[itone][jtone]=0;       /* protect from error */
            else matrix[itone][jtone] /= count;

    // Return the co-occurrence matrix
    return matrix;
}






// Defining all quantifying functions derived from co-occurance matrix

/* Angular Second Moment
* The angular second-moment feature (ASM) f1 is a measure of homogeneity
* of the image. In a homogeneous image, there are very few dominant
* gray-tone transitions. Hence the P matrix for such an image will have
* fewer entries of large magnitude.
*/
double f1_asm (double **P, int Ng) {
    int i, j;
    double sum = 0;

    // Loop through the co-occurrence matrix to calculate the ASM
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            sum += P[i][j] * P[i][j]; // Calculate the sum of the squared values of the co-occurrence matrix

    // Return the ASM
    return sum;
}



/* Contrast
* The contrast feature is a difference moment of the P matrix and is a
* measure of the contrast or the amount of local variations present in an
* image.
*/
double f2_contrast (double **P, int Ng) {
    int i, j, n;
    double sum = 0, bigsum = 0;

    for (n = 0; n < Ng; ++n) {
        for (i = 0; i < Ng; ++i)
            for (j = 0; j < Ng; ++j) {
                if ((i - j) == n || (j - i) == n)
                    sum += P[i][j];
            }
        bigsum += n * n * sum;
        sum = 0;
    }

    return bigsum;
}


/* Correlation
*
* This correlation feature is a measure of gray-tone linear-dependencies
* in the image.
*/
double f3_corr (double **P, int Ng) {
    int i, j;
    double sum_sqrx = 0, sum_sqry = 0, tmp, *px;
    double meanx =0 , meany = 0 , stddevx, stddevy;

    // Allocate memory for the marginal probability matrix
    px = allocate_vector (0, Ng);
    for (i = 0; i < Ng; ++i)
        px[i] = 0;

    // Calculate the marginal probability matrix by summing the rows of the co-occurrence matrix
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            px[i] += P[i][j];

    // Calculate the mean and standard deviation of px and py
    for (i = 0; i < Ng; ++i) {
        meanx += px[i]*i;
        sum_sqrx += px[i]*i*i;
    }

    // Calculate the mean and standard deviation of py
    meany = meanx;
    sum_sqry = sum_sqrx;
    stddevx = sqrt (sum_sqrx - (meanx * meanx));
    stddevy = stddevx;

    // Calculate the correlation feature
    for (tmp = 0, i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            tmp += i*j*P[i][j];

    // Free the memory allocated for the marginal probability matrix
    free(px);

    // Protect from error by checking if the denominator is zero
    if (stddevx * stddevy==0) return(1);
    else return (tmp - meanx * meany) / (stddevx * stddevy);
}

/* Sum of Squares: Variance */
double f4_var (double **P, int Ng) {
    int i, j;
    double mean = 0, var = 0;


    // Calculate the mean intensity level of the image
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            mean += i * P[i][j];

    // Calculate the variance of the co-occurrence matrix
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            var += (i - mean) * (i - mean) * P[i][j];

    // Return the variance value
    return var;
}

/* Inverse Difference Moment */
double f5_idm (double **P, int Ng) {
    int i, j;
    double idm = 0;

    // loop through the co-occurrence matrix and calculate the idm value
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            idm += P[i][j] / (1 + (i - j) * (i - j));

    // return the idm value
    return idm;
}
/* Sum Average */
double f6_savg (double **P, int Ng) {
    int i, j;
    double savg = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    // Initialize Pxpy vector to 0
    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    // Calculate Pxpy vector
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[i + j] += P[i][j];

    // Calculate sum average value
    for (i = 0; i <= (2 * Ng - 2); ++i)
        savg += i * Pxpy[i];

    // Free memory allocated for Pxpy vector
    free (Pxpy);

    // Return the sum average value
    return savg;
}

/* Sum Variance */
double f7_svar (double **P, int Ng, double S) {
    int i, j;
    double var = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    // Initialize Pxpy vector to 0
    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    // Calculate Pxpy vector
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[i + j] += P[i][j];

    // Calculate sum variance value
    for (i = 0; i <= (2 * Ng - 2); ++i)
        var += (i - S) * (i - S) * Pxpy[i];

    // Free memory allocated for Pxpy vector
    free (Pxpy);

    // Return the sum variance value
    return var;
}


/* Sum Entropy */
double f8_sentropy (double **P, int Ng) {
    int i, j;
    double sentropy = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    // Initialize Pxpy vector to 0
    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    // Calculate Pxpy vector
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[i + j + 2] += P[i][j];

    // Calculate sum entropy value
    for (i = 2; i <= 2 * Ng; ++i)
        sentropy -= Pxpy[i] * log10 (Pxpy[i] + EPSILON)/log10(2.0) ;

    // Free memory allocated for Pxpy vector
    free (Pxpy);

    // Return the sum entropy value
    return sentropy;
}

/* Entropy */
double f9_entropy (double **P, int Ng) {
    int i, j;
    double entropy = 0;

    // Calculate entropy value
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            entropy += P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0) ;

    // Return the entropy value
    return -entropy;
}

/* Difference Variance */
double f10_dvar (double **P, int Ng) {
    int i, j;
    double sum = 0, sum_sqr = 0, var = 0;
    double *Pxpy = allocate_vector (0, 2*Ng); // Allocate memory for Pxpy vector

    // Initialize Pxpy vector to 0
    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    // Calculate Pxpy vector
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[abs (i - j)] += P[i][j];

    // Calculate the variance of Pxpy (Px-y)
    for (i = 0; i < Ng; ++i) {
        sum += i * Pxpy[i] ;
        sum_sqr += i * i * Pxpy[i] ;
    }

    var = sum_sqr - sum*sum ; // Calculate the variance

    free (Pxpy); // Free memory allocated for Pxpy vector
    return var; // Return the variance value
}


/* Difference Entropy */
double f11_dentropy (double **P, int Ng) {
    int i, j;
    double sum = 0;
    double *Pxpy = allocate_vector (0, 2*Ng); // Allocate memory for Pxpy vector

    // Initialize Pxpy vector to 0
    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    // Calculate Pxpy vector
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[abs (i - j)] += P[i][j];

    // Calculate the difference entropy value
    for (i = 0; i < Ng; ++i)
        sum += Pxpy[i] * log10 (Pxpy[i] + EPSILON)/log10(2.0) ;

    free (Pxpy); // Free memory allocated for Pxpy vector
    return -sum; // Return the difference entropy value
}


/* Information Measures of Correlation */
double f12_icorr (double **P, int Ng) {
    int i, j;
    double *px, *py;
    double hx = 0, hy = 0, hxy = 0, hxy1 = 0, hxy2 = 0;

    px = allocate_vector (0, Ng);
    py = allocate_vector (0, Ng);

    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            px[i] += P[i][j];
            py[j] += P[i][j];
        }
    }

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j) {
            hxy1 -= P[i][j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy2 -= px[i] * py[j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy -= P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0);
        }

    /* Calculate entropies of px and py */
    for (i = 0; i < Ng; ++i) {
        hx -= px[i] * log10 (px[i] + EPSILON)/log10(2.0);
        hy -= py[i] * log10 (py[i] + EPSILON)/log10(2.0);
    }

    free(px);
    free(py);
    if ((hx > hy ? hx : hy)==0) return(1);
    else
        return ((hxy - hxy1) / (hx > hy ? hx : hy));
}


/* Information Measures of Correlation */
double f13_icorr (double **P, int Ng) {
    int i, j;
    double *px, *py;
    double hx = 0, hy = 0, hxy = 0, hxy1 = 0, hxy2 = 0;

    px = allocate_vector (0, Ng);
    py = allocate_vector (0, Ng);

    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            px[i] += P[i][j];
            py[j] += P[i][j];
        }
    }

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j) {
            hxy1 -= P[i][j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy2 -= px[i] * py[j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy -= P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0);
        }

    /* Calculate entropies of px and py */
    for (i = 0; i < Ng; ++i) {
        hx -= px[i] * log10 (px[i] + EPSILON)/log10(2.0);
        hy -= py[i] * log10 (py[i] + EPSILON)/log10(2.0);
    }

    free(px);
    free(py);
    return (sqrt (fabs (1 - exp (-2.0 * (hxy2 - hxy)))));
}


/* Returns the Maximal Correlation Coefficient */
double f14_maxcorr (double **P, int Ng) {
    int i, j, k;
    double *px, *py, **Q;
    double *x, *iy, tmp;
    double f=0.0;

    px = allocate_vector (0, Ng);
    py = allocate_vector (0, Ng);
    Q = allocate_matrix (1, Ng + 1, 1, Ng + 1);
    x = allocate_vector (1, Ng);
    iy = allocate_vector (1, Ng);

    /*
    * px[i] is the (i-1)th entry in the marginal probability matrix obtained
    * by summing the rows of p[i][j]
    */
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            px[i] += P[i][j];
            py[j] += P[i][j];
        }
    }

    /* Find the Q matrix */
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            Q[i + 1][j + 1] = 0;
            for (k = 0; k < Ng; ++k)
                if (px[i] && py[k])  /* make sure to protect division by zero */
                    Q[i + 1][j + 1] += P[i][k] * P[j][k] / px[i] / py[k];
        }
    }

    /* Balance the matrix */
    mkbalanced (Q, Ng);
    /* Reduction to Hessenberg Form */
    reduction (Q, Ng);
    /* Finding eigenvalue for nonsymetric matrix using QR algorithm */
    if (!hessenberg (Q, Ng, x, iy)) {
        /* Memmory cleanup */
        for (i=1; i<=Ng+1; i++) free(Q[i]+1);
        free(Q+1);
        free((char *)px);
        free((char *)py);
        free((x+1));
        free((iy+1));

        /* computation failed ! */
        return 0.0;
    }

    /* simplesrt(Ng,x); */
    /* Returns the sqrt of the second largest eigenvalue of Q */
    for (i = 2, tmp = x[1]; i <= Ng; ++i)
        tmp = (tmp > x[i]) ? tmp : x[i];

    if (x[Ng - 1]>=0)
        f = sqrt(x[Ng - 1]);

    for (i=1; i<=Ng+1; i++) free(Q[i]+1);
    free(Q+1);
    free((char *)px);
    free((char *)py);
    free((x+1));
    free((iy+1));

    return f;
}
// This function takes in an image, angle, and distance as input parameters and returns a tuple containing the grayco matrix and the number of different tones in the image
tuple<__decay_and_strip<double **&>::__type, int> grayco_matrcis(Mat img, int angle, int distance){

    // Initialize variables
    int row, col, rows, cols;
    cols = img.cols;
    rows = img.rows;

    clock_t start_gray_pre,end_gray_pre,start_co,end_co;
    double cpu_time_used_gray_pre,cpu_time_used_co;

    // Start timer for grayco matrix preprocessing
    start_gray_pre = clock();

    // Create a 2D array **pGray and copy all the image values in **pGray
    unsigned char **pGray;
    pGray = new unsigned char *[rows];

    for(int i = 0; i < rows; i++){
        pGray[i] = new unsigned char[cols];
    }

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < cols; x++)
        {
            // Copy image values to **pGray
            pGray[y][x] = img.data[img.step * y + x * 1];
        }
    }

    int toneLUT[PGM_MAXMAXVAL + 1];		// toneLUT is an array that can hold 256 values
    int toneCount = 0;
    int iTone;

    // Fill toneLUT with -1
    for(row = PGM_MAXMAXVAL; row >= 0; --row)
        toneLUT[row] = -1;

    // Fill toneLUT with those 8 bit values which are present in the image
    // Example: if the image has values 0,1,2,3,4, then toneLUT will have values only from 0 - 4.
    for(row = rows - 1; row >= 0; --row){
        for(col = 0; col < cols; ++col){
            toneLUT[(u_int8_t)img.data[img.step * row + col * 1]] = (u_int8_t)img.data[img.step * row + col * 1];
        }
    }

    // ToneCount contains the number of 8-bit value variations in the image
    for (row = PGM_MAXMAXVAL, toneCount = 0; row >= 0; --row){
        if (toneLUT[row] != -1)
            toneCount++;
        else
            ;
    }

    /* Use the number of different tones to build LUT */
    for (row = 0, iTone = 0; row <= PGM_MAXMAXVAL; row++){
        if (toneLUT[row] != -1)
            toneLUT[row] = iTone++;
    }

    // End timer for grayco matrix preprocessing
    end_gray_pre = clock();
    cpu_time_used_gray_pre= ((double) (end_gray_pre - start_gray_pre)) / CLOCKS_PER_SEC;

    double **pMatrix;

    // Start timer for calculating grayco matrices
    start_co = clock();

    // Calculate grayco matrices based on the angle
    if (angle == 0) {
        pMatrix = CoOcMat_Angle_0(distance, pGray, rows, cols, toneLUT, toneCount);
    }else if(angle == 45){
        pMatrix = CoOcMat_Angle_45(distance, pGray, rows, cols, toneLUT, toneCount);
    }else if(angle == 90){
        pMatrix = CoOcMat_Angle_90(distance, pGray, rows, cols, toneLUT, toneCount);
    }else if(angle == 135){
        pMatrix = CoOcMat_Angle_135(distance, pGray, rows, cols, toneLUT, toneCount);
    }

    // End timer for calculating grayco matrices
    end_co = clock();
    cpu_time_used_co= ((double) (end_co - start_co)) / CLOCKS_PER_SEC;

    // Free memory allocated for pGray
    delete [] pGray;

    // Return a tuple containing the grayco matrix and the number of different tones in the image
    return make_tuple(pMatrix, toneCount);
}
// This function takes in an image, an array of feature parameters, the length of the array, an angle, and a distance as input parameters.
// It calculates the grayco matrix and then calculates the specified features from the matrix.
// It returns an array of the calculated feature values.
double  *calculator(Mat img,char* param[],int length,int angle,int distance){

    int toneCount;
    double **pMatrix;

    // Start the clock to measure the time taken to generate the grayco matrix
    clock_t start_GR,end_GR,start_PR,end_PR;
    double cpu_time_used_GR, cpu_time_used_PR;
    start_GR = clock();

    // Calculate the grayco matrix and the number of different tones in the image
    tie (pMatrix,toneCount) = grayco_matrcis(img,angle,distance);

    // Set first row and column to zero
    for (int i = 0; i < toneCount; ++i) {
        pMatrix[0][i] = 0;
        pMatrix[i][0] = 0;
    }

    // Calculate the sum of all elements in the matrix
    double sum = 0;
    for (int i = 0; i < toneCount; ++i) {
        for (int j = 0; j < toneCount; ++j) {
            sum += pMatrix[i][j];
        }
    }

    // Stop the clock and calculate the time taken to generate the grayco matrix
    end_GR = clock();
    cpu_time_used_GR = ((double) (end_GR - start_GR)) / CLOCKS_PER_SEC;

    // Start the clock to measure the time taken to calculate the specified features
    start_PR = clock();

    // Create an array to store the calculated feature values
    double *result = new double[length];
    for (int i = 0;i<length;i++) {
        if (strcmp(param[i],"asm")==0) {
            result[i] = round(f1_asm(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "contrast") == 0) {
            result[i] = round(f2_contrast(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "corr") == 0) {
            result[i] = round(f3_corr(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "var") == 0) {
            result[i] = round(f4_var(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "idm") == 0) {
            result[i] = round(f5_idm(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "savg") == 0) {
            result[i] = round(f6_savg(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "svar") == 0) {
            result[i] = round(f7_svar(pMatrix, toneCount, f8_sentropy(pMatrix, toneCount)) * 10000) / 10000;
        } else if (strcmp(param[i], "sentropy") == 0) {
            result[i] = round(f8_sentropy(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "entropy") == 0) {
            result[i] = round(f9_entropy(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "dvar") == 0) {
            result[i] = round(f10_dvar(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "dentropy") == 0) {
            result[i] = round(f11_dentropy(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "icorr1") == 0) {
            result[i] = round(f12_icorr(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "icorr2") == 0) {
            result[i] = round(f13_icorr(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "maxcorr") == 0) {
            result[i] = round(f14_maxcorr(pMatrix, toneCount) * 10000) / 10000;
        } else if (strcmp(param[i], "sum") == 0) {
            result[i] = sum;
        }
    }

    // Stop the clock and calculate the time taken to calculate the specified features
    end_PR = clock();
    cpu_time_used_PR = ((double) (end_PR - start_PR)) / CLOCKS_PER_SEC;

    // Free the memory allocated for the grayco matrix
    free_matrix(pMatrix, toneCount);

    // Return the array of calculated feature values
    return result;
}