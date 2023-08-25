#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "GEcancer.h"


// This function calculates the F1 score for a given set of evolved solutions, dynamic fitness cases, and classes.
// The F1 score is a measure of a test's accuracy, calculated as the harmonic mean of the precision and recall.
// The function takes in three integer arrays: Evolved, Class, and DinamicFitnessCases.
// Evolved is an array of integers representing the evolved solutions.
// Class is an array of integers representing the classes.
// DinamicFitnessCases is an integer representing the number of dynamic fitness cases.
// The function returns a double representing the F1 score.
double F1(int * Evolved, int DinamicFitnessCases, int * Class) {
  int i, TP, FP, TN, FN;
  double F1, F1_divisor;
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  F1 = 0;

  // Loop through each dynamic fitness case.
  for(i = 0; i < DinamicFitnessCases; i ++)
  {
    // If the evolved solution is 1 and the class is 1, increment the true positive count.
    if (Evolved[i] == 1)
    {
      if (Class[i] == 1)
      {
        ++ TP;
      }
      // If the evolved solution is 1 and the class is 0, increment the false positive count.
      else
      {
        ++ FP;
      }
    }
    // If the evolved solution is 0 and the class is 0, increment the true negative count.
    else
    {
      if (Class[i] == 0)
      {
        ++ TN;
      }
      // If the evolved solution is 0 and the class is 1, increment the false negative count.
      else
      {
        ++ FN;
      }
    }    
  }

  // Calculate the F1 score.
  F1_divisor = 2.0 * TP + FP + FN;
  if (F1_divisor > 0)
  {
    F1 = (2.0 * TP) / F1_divisor;  
  }
  return F1;
}

// This function calculates the Balanced Accuracy (BA) for a given set of evolved solutions, dynamic fitness cases, and classes.
// The Balanced Accuracy is a measure of a test's accuracy, calculated as the arithmetic mean of the true positive rate (TPR) and true negative rate (TNR).
// The function takes in three integer arrays: Evolved, Class, and DinamicFitnessCases.
// Evolved is an array of integers representing the evolved solutions.
// Class is an array of integers representing the classes.
// DinamicFitnessCases is an integer representing the number of dynamic fitness cases.
// The function returns a double representing the Balanced Accuracy.
double BA(int * Evolved, int DinamicFitnessCases, int * Class) {
  int i, TP, FP, TN, FN;
  double BA, TPR, TNR;
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  BA = 0;
  TPR = 0;
  TNR = 0;

  // Loop through each dynamic fitness case.
  for(i = 0; i < DinamicFitnessCases; i ++)
  {
    // If the evolved solution is 1 and the class is 1, increment the true positive count.
    if (Evolved[i] == 1)
    {
      if (Class[i] == 1)
      {
        ++ TP;
      }
      // If the evolved solution is 1 and the class is 0, increment the false positive count.
      else
      {
        ++ FP;
      }
    }
    // If the evolved solution is 0 and the class is 0, increment the true negative count.
    else
    {
      if (Class[i] == 0)
      {
        ++ TN;
      }
      // If the evolved solution is 0 and the class is 1, increment the false negative count.
      else
      {
        ++ FN;
      }
    }    
  }

  // Calculate the true positive rate (TPR) and true negative rate (TNR).
  if (TP > 0)
  {
    TPR = (double)TP / ((double)TP + (double)FN);
  }
  if (TN > 0)
  {
    TNR = (double)TN / ((double)TN + (double)FP);
  }

  // Calculate the Balanced Accuracy (BA).
  BA = (TPR + TNR) / 2.0;
  return BA;
}

// This function calculates the True Positive Rate (TPR) for a given set of evolved solutions, dynamic fitness cases, and classes.
// The True Positive Rate is a measure of a test's ability to correctly identify positive results, calculated as the ratio of true positives to the sum of true positives and false negatives.
// The function takes in three integer arrays: Evolved, Class, and DinamicFitnessCases.
// Evolved is an array of integers representing the evolved solutions.
// Class is an array of integers representing the classes.
// DinamicFitnessCases is an integer representing the number of dynamic fitness cases.
// The function returns a double representing the True Positive Rate.
double TPR(int * Evolved, int DinamicFitnessCases, int * Class)
{
  int i, TP, FN;
  double TPR;
  TP = 0;
  FN = 0;
  TPR = 0;

  // Loop through each dynamic fitness case.
  for(i = 0; i < DinamicFitnessCases; i ++)
  {
    // If the evolved solution is 1 and the class is 1, increment the true positive count.
    if (Evolved[i] == 1)
    {
      if (Class[i] == 1)
      {
        ++ TP;
      }
      // If the evolved solution is 1 and the class is 0, increment the false negative count.
      else
      {
        ++ FN;
      }
    }
    // If the evolved solution is 0 and the class is 1, increment the false negative count.
    else
    {
      if (Class[i] == 1)
      {
        ++ FN;
      }
    }    
  }

  // Calculate the True Positive Rate (TPR).
  if (TP > 0)
  {
    TPR = (double)TP / ((double)TP + (double)FN);
  }
  return TPR;
}

// This function calculates the False Positive Rate (FPR) for a given set of evolved solutions, dynamic fitness cases, and classes.
// The False Positive Rate is a measure of a test's ability to correctly identify negative results, calculated as the ratio of false positives to the sum of false positives and true negatives.
// The function takes in three integer arrays: Evolved, Class, and DinamicFitnessCases.
// Evolved is an array of integers representing the evolved solutions.
// Class is an array of integers representing the classes.
// DinamicFitnessCases is an integer representing the number of dynamic fitness cases.
// The function returns a double representing the False Positive Rate.
double FPR(int * Evolved, int DinamicFitnessCases, int * Class)
{
  int i, FP, TN;
  double FPR;
  FP = 0;
  TN = 0;
  FPR = 0;

  // Loop through each dynamic fitness case.
  for(i = 0; i < DinamicFitnessCases; i ++)
  {
    // If the evolved solution is 1 and the class is 0, increment the false positive count.
    if (Evolved[i] == 1)
    {
      if (Class[i] == 0)
      {
        ++ FP;
      }
    }
    // If the evolved solution is 0 and the class is 0, increment the true negative count.
    else
    {
      if (Class[i] == 0)
      {
        ++ TN;
      }
    }    
  }

  // Calculate the False Positive Rate (FPR).
  if (FP > 0)
  {
    FPR = (double)FP / ((double)FP + (double)TN);
  }
  return FPR;
}


// This function returns the value of exp2 if exp1 is less than or equal to pivot, otherwise it returns the value of exp3.
// The function takes in four doubles: exp1, pivot, exp2, and exp3.
// exp1 is the first expression to compare with pivot.
// pivot is the value to compare exp1 with.
// exp2 is the value to return if exp1 is less than or equal to pivot.
// exp3 is the value to return if exp1 is greater than pivot.
double ifLessThan(double exp1, double pivot, double exp2, double exp3) {
  if (exp1 <= pivot)
    return exp2;
  else
    return exp3;
}

// This function performs a protected division of a1 by a2.
// If a2 is equal to 0, the function returns 0 to avoid division by zero errors.
// Otherwise, the function returns the result of dividing a1 by a2.
// The function takes in two doubles: a1 and a2.
double protected_div(double a1, double a2) {
  if(a2 == 0)
    return 0;
  else
    return a1/a2;
}

// This function performs a protected logarithm of a1.
// If a1 is less than or equal to 0, the function returns 0 to avoid logarithm errors.
// Otherwise, the function returns the natural logarithm of a1.
// The function takes in one double: a1.
double protected_log(double a1) {
  if(a1 <= 0)
    return 0;
  else
    return log(a1);
}

// This function performs a protected absolute value of a1.
// The function returns the absolute value of a1.
// The function takes in one double: a1.
double protected_abs(double a1) {
  return (a1 > 0.0) ? a1 : ((-1.0) * a1);
}

// This function performs a protected square root of a1.
// If the square root of a1 is not a finite number, the function returns 0 to avoid errors.
// Otherwise, the function returns the square root of a1.
// The function takes in one double: a1.
double protected_sqrt(double a1) {
  if(isfinite(sqrt(a1)))
    return sqrt(a1);
  else
    return 0;
}

// This function performs a protected exponential of a1.
// If the exponential of a1 is not a finite number, the function returns 0 to avoid errors.
// Otherwise, the function returns the exponential of a1.
// The function takes in one double: a1.
double protected_exp(double a1) {
  if(isfinite(exp(a1)))
    return exp(a1);
  else
    return 0;
}
