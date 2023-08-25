#include "../src/GEcancer.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>

// Define an array of expressions using lambda functions
std::array<std::function<double(double **Dim, int i)>, 25> expressions = {    [](double** Dim, int i) { return  protected_log(protected_log(protected_log(((((Dim[i][64]*Dim[i][20])-((Dim[i][22]-(sin((Dim[i][64]))*protected_div(Dim[i][63], Dim[i][61])))))-Dim[i][71])-Dim[i][47])))); },
    [](double** Dim, int i) { return  protected_log((((protected_div(Dim[i][18], Dim[i][40])*Dim[i][15])*Dim[i][80])-(Dim[i][55]+protected_sqrt((Dim[i][55]+protected_sqrt((ifLessThan(Dim[i][54], Dim[i][84], protected_div(Dim[i][16], ifLessThan(protected_sqrt(ifLessThan(Dim[i][18], Dim[i][54], Dim[i][84], (ifLessThan((Dim[i][54]), Dim[i][81], Dim[i][55], (Dim[i][9]))*Dim[i][97]))), protected_sqrt(protected_log(((Dim[i][93]-(Dim[i][9]*Dim[i][35]))+Dim[i][8]))), Dim[i][23], tanh((ifLessThan((Dim[i][45]), Dim[i][57], (Dim[i][90]*Dim[i][33]), protected_sqrt((Dim[i][80]))))))), Dim[i][70])+protected_sqrt((Dim[i][80]))))))))); },
    [](double** Dim, int i) { return  protected_log((Dim[i][19]+protected_log(((Dim[i][50]*(Dim[i][58]*Dim[i][60]))+protected_log((Dim[i][50]*(Dim[i][41]*Dim[i][11]))))))); },
    [](double** Dim, int i) { return  protected_log(((((Dim[i][110]*protected_div((Dim[i][30]*Dim[i][78]), Dim[i][24]))))+ifLessThan((Dim[i][80]-Dim[i][16]), protected_log((Dim[i][30]*Dim[i][78])), Dim[i][90], sin(protected_log((Dim[i][93]-Dim[i][20])))))); },
    [](double** Dim, int i) { return  protected_log((Dim[i][47]+(sin(protected_log(protected_div(Dim[i][55], (Dim[i][80]-Dim[i][9]))))*ifLessThan(protected_abs(protected_div(Dim[i][43], (Dim[i][47]*protected_div(Dim[i][55], (Dim[i][5]-Dim[i][97]))))), ifLessThan(Dim[i][87], Dim[i][1], (Dim[i][43]-protected_abs(ifLessThan(Dim[i][70], protected_log(protected_div(Dim[i][7], protected_div(Dim[i][35], (protected_log(protected_sqrt(Dim[i][47]))-sin(protected_log(protected_div(Dim[i][55], (Dim[i][80]-ifLessThan(Dim[i][25], protected_abs(sin(protected_log(Dim[i][113]))), sin(Dim[i][55]), (protected_log(Dim[i][108])*Dim[i][49])))))))))), protected_log((Dim[i][119]*ifLessThan(sin((Dim[i][88]-Dim[i][1])), Dim[i][38], Dim[i][108], Dim[i][113]))), protected_abs(Dim[i][115])))), Dim[i][83]), protected_abs((Dim[i][93]*Dim[i][64])), ifLessThan(Dim[i][2], Dim[i][51], Dim[i][55], Dim[i][7]))))); },
    [](double** Dim, int i) { return  protected_log(((((Dim[i][107]+protected_log(((Dim[i][3]*tanh(ifLessThan(Dim[i][30], Dim[i][3], Dim[i][105], Dim[i][101])))-Dim[i][56]))))*ifLessThan(protected_log(((Dim[i][3]*tanh(ifLessThan(Dim[i][30], Dim[i][3], Dim[i][25], Dim[i][14])))-Dim[i][4])), Dim[i][4], Dim[i][104], Dim[i][67])))); },
    [](double** Dim, int i) { return  protected_log(protected_log(protected_log((((Dim[i][36]*Dim[i][105])-Dim[i][38])-(Dim[i][101]*((Dim[i][50]-Dim[i][105])-(Dim[i][38]-((protected_log(((Dim[i][3])))-Dim[i][61]))))))))); },
    [](double** Dim, int i) { return  protected_log(protected_div(protected_div((Dim[i][33]-protected_sqrt(protected_div(Dim[i][46], Dim[i][105]))), (Dim[i][27])), protected_div(Dim[i][83], protected_div(Dim[i][34], protected_div((((Dim[i][47])+Dim[i][4])+(Dim[i][17]+Dim[i][109])), protected_log(Dim[i][97])))))); },
    [](double** Dim, int i) { return  ifLessThan(Dim[i][113], Dim[i][92], Dim[i][105], protected_div(tanh(ifLessThan(Dim[i][9], Dim[i][94], Dim[i][85], protected_log((Dim[i][2]*(protected_div(Dim[i][92], Dim[i][41])+protected_log((((Dim[i][104]*Dim[i][5]))))))))), Dim[i][38])); },
    [](double** Dim, int i) { return  sin((Dim[i][71]*((Dim[i][38]-protected_log(protected_log(protected_log(protected_div(Dim[i][82], sin(Dim[i][40]))))))))); },
    [](double** Dim, int i) { return  (sin(protected_log(Dim[i][33]))+ifLessThan(Dim[i][8], Dim[i][69], ifLessThan(Dim[i][8], Dim[i][69], protected_log(((protected_sqrt(protected_div((Dim[i][65]), Dim[i][115])))*Dim[i][64])), (Dim[i][112]-Dim[i][67])), ((sin(protected_log(Dim[i][33]))+ifLessThan(Dim[i][8], Dim[i][69], (sin(Dim[i][4])), Dim[i][44]))-protected_abs(((protected_log((sin(tanh(protected_abs(((Dim[i][23]*protected_div(Dim[i][8], Dim[i][69]))))))-Dim[i][56]))*ifLessThan(Dim[i][108], protected_sqrt(Dim[i][74]), Dim[i][112], protected_div(Dim[i][0], (Dim[i][112]-ifLessThan(Dim[i][4], sin(protected_sqrt(protected_div(Dim[i][65], Dim[i][115]))), Dim[i][72], Dim[i][106])))))))))); },
    [](double** Dim, int i) { return  (protected_log(protected_sqrt(protected_log(protected_log(protected_div(Dim[i][75], protected_div(protected_div(Dim[i][37], Dim[i][33]), protected_div(protected_div(Dim[i][75], protected_div(Dim[i][9], (Dim[i][66]))), protected_div(protected_div(Dim[i][37], Dim[i][33]), (Dim[i][38]+protected_log(protected_sqrt(ifLessThan((Dim[i][109]), Dim[i][107], Dim[i][18], (Dim[i][80])))))))))))))); },
    [](double** Dim, int i) { return  protected_log((Dim[i][30]+(Dim[i][30]+protected_abs(protected_log((Dim[i][93]*protected_sqrt(protected_div(Dim[i][15], (protected_div(protected_div(Dim[i][4], Dim[i][100]), Dim[i][59])-protected_div(Dim[i][56], sin((Dim[i][71]*Dim[i][23])))))))))))); },
    [](double** Dim, int i) { return  ((Dim[i][60]*Dim[i][78])-(protected_div((sin(protected_div(Dim[i][113], protected_div(Dim[i][50], Dim[i][40])))+Dim[i][60]), Dim[i][62]))); },
    [](double** Dim, int i) { return  protected_log(((protected_abs((Dim[i][40]*(Dim[i][41]+protected_div(Dim[i][52], protected_div(Dim[i][80], Dim[i][107])))))+ifLessThan(Dim[i][107], Dim[i][1], ((Dim[i][109]+ifLessThan(Dim[i][107], protected_log(protected_abs((Dim[i][94]))), protected_abs((Dim[i][40]*(Dim[i][41]+protected_div(Dim[i][75], Dim[i][97])))), (protected_abs(((Dim[i][41]*Dim[i][107])*ifLessThan(Dim[i][107], Dim[i][10], sin(Dim[i][32]), Dim[i][40])))*protected_abs(((Dim[i][41]*Dim[i][107])*ifLessThan(Dim[i][107], Dim[i][1], (protected_abs((Dim[i][22]*Dim[i][56]))), Dim[i][104]))))))), Dim[i][41])))); },
    [](double** Dim, int i) { return  protected_log((((Dim[i][71])+((protected_div(Dim[i][60], Dim[i][85])*Dim[i][33])*protected_abs(protected_div(Dim[i][20], Dim[i][97])))))); },
    [](double** Dim, int i) { return  protected_log(((protected_log((protected_sqrt(Dim[i][110])*sin(tanh(protected_div((protected_div(Dim[i][10], sin(((Dim[i][41]+Dim[i][75])*protected_sqrt(Dim[i][110]))))), Dim[i][55])))))))); },
    [](double** Dim, int i) { return  protected_log((protected_div((Dim[i][32]*(Dim[i][15]*protected_div(Dim[i][48], protected_div(Dim[i][76], Dim[i][63])))), protected_log((protected_log((protected_log((Dim[i][15]*protected_div(Dim[i][48], protected_div(Dim[i][76], (Dim[i][32]*Dim[i][63])))))*Dim[i][70]))))))); },
    [](double** Dim, int i) { return  sin((Dim[i][56]*ifLessThan(Dim[i][9], Dim[i][99], (Dim[i][37]+Dim[i][52]), (sin((Dim[i][56]*Dim[i][116]))*Dim[i][110])))); },
    [](double** Dim, int i) { return  (protected_log((protected_log(((Dim[i][30]+(Dim[i][4]*(protected_abs((Dim[i][5]*protected_abs((Dim[i][86]))))-Dim[i][98])))))))); },
    [](double** Dim, int i) { return  protected_log((Dim[i][75]+(((Dim[i][5]*protected_sqrt(protected_log((Dim[i][30]*Dim[i][44]))))+(Dim[i][5]*protected_div((Dim[i][30]*Dim[i][48]), Dim[i][25])))-Dim[i][112]))); },
    [](double** Dim, int i) { return  protected_div(protected_log((protected_log((protected_div((Dim[i][90]*Dim[i][78]), Dim[i][55])*(protected_div(ifLessThan((Dim[i][5]*Dim[i][94]), Dim[i][103], Dim[i][77], Dim[i][5]), Dim[i][25]))))+(Dim[i][90]*Dim[i][78]))), Dim[i][76]); },
    [](double** Dim, int i) { return  sin((Dim[i][23]*(ifLessThan(Dim[i][122], Dim[i][41], ifLessThan(Dim[i][67], Dim[i][23], protected_div(sin(Dim[i][30]), protected_div(sin((Dim[i][23]*ifLessThan(sin(Dim[i][90]), Dim[i][15], Dim[i][71], Dim[i][41]))), protected_div(protected_abs(Dim[i][85]), Dim[i][114]))), Dim[i][102]), ifLessThan(protected_div(protected_log(Dim[i][85]), Dim[i][3]), sin(Dim[i][90]), Dim[i][15], Dim[i][71]))))); },
    [](double** Dim, int i) { return  protected_log(((ifLessThan(Dim[i][19], Dim[i][89], Dim[i][108], Dim[i][75])*Dim[i][90])+protected_log(protected_log(((ifLessThan(Dim[i][19], Dim[i][113], Dim[i][36], Dim[i][91])*Dim[i][90])+(sin(protected_log(Dim[i][48]))*protected_sqrt(protected_div(Dim[i][108], Dim[i][91])))))))); },
    [](double** Dim, int i) { return  protected_log((Dim[i][90]+protected_log(protected_div((Dim[i][3]*ifLessThan(Dim[i][108], (Dim[i][9]*Dim[i][71]), Dim[i][10], Dim[i][30])), protected_div(Dim[i][53], (Dim[i][65]-Dim[i][53])))))); }
    
};

/// @tparam T1 Type of array elements, should be a numerical type
/// @tparam T2 Type of array elements, should be a numerical type
/// @param label Array of ground truth labels, 0 is negative, 1 is positive
/// @param score Array of predicted scores, can be any real finite number
/// @param n Number of elements in the array, I assume it's correct
/// @return AUROC/ROC-AUC score, range [0.0, 1.0]
template <class T1, class T2>
double AUROC(const T1 label[], const T2 score[], int n) {
  for (int i = 0; i < n; i++)
    if (!std::isfinite(score[i]) || label[i] != 0 && label[i] != 1)
      return std::numeric_limits<double>::signaling_NaN();

  const auto order = new int[n];
  std::iota(order, order + n, 0);
  std::sort(order, order + n,
            [&](int a, int b) { return score[a] > score[b]; });
  const auto y = new double[n];
  const auto z = new double[n];
  for (int i = 0; i < n; i++) {
    y[i] = label[order[i]];
    z[i] = score[order[i]];
  }

  const auto tp = y; // Reuse
  std::partial_sum(y, y + n, tp);

  int top = 0; // # diff
  for (int i = 0; i < n - 1; i++)
    if (z[i] != z[i + 1])
      order[top++] = i;
  order[top++] = n - 1;
  n = top; // Size of y/z -> sizeof tps/fps

  const auto fp = z; // Reuse
  for (int i = 0; i < n; i++) {
    tp[i] = tp[order[i]];         // order is mono. inc.
    fp[i] = 1 + order[i] - tp[i]; // Type conversion prevents vectorization
  }
  delete[] order;

  const auto tpn = tp[n - 1], fpn = fp[n - 1];
  for (int i = 0; i < n; i++) { // Vectorization
    tp[i] /= tpn;
    fp[i] /= fpn;
  }

  auto area = tp[0] * fp[0] / 2; // The first triangle from origin;
  double partial = 0;            // For Kahan summation
  for (int i = 1; i < n; i++) {
    const auto x = (fp[i] - fp[i - 1]) * (tp[i] + tp[i - 1]) / 2 - partial;
    const auto sum = area + x;
    partial = (sum - area) - x;
    area = sum;
  }

  delete[] tp;
  delete[] fp;

  return area;
}

double F1(int *Evolved, int DinamicValidationFitnessCases, int *Class) {
  int i;
  double TP, FP, TN, FN;
  double F1, F1_divisor;
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  F1 = 0;
  for (i = 0; i < DinamicValidationFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 1) {
        ++TP;
      } else {
        ++FP;
      }
    } else {
      if (Class[i] == 0) {
        ++TN;
      } else {
        ++FN;
      }
    }
  }
  F1_divisor = 2.0 * TP + FP + FN;
  if (F1_divisor > 0) {
    F1 = (2.0 * TP) / F1_divisor;
  }
  return F1;
}

double BA(int *Evolved, int DinamicValidationFitnessCases, int *Class) {
  int i;
  double TP, FP, TN, FN;
  double BA, TPR, TNR;
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  BA = 0;
  TPR = 0;
  TNR = 0;
  for (i = 0; i < DinamicValidationFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 1) {
        ++TP;
      } else {
        ++FP;
      }
    } else {
      if (Class[i] == 0) {
        ++TN;
      } else {
        ++FN;
      }
    }
  }
  if (TP > 0) {
    TPR = TP / (TP + FN);
  }
  if (TN > 0) {
    TNR = TN / (TN + FP);
  }
  BA = (TPR + TNR) / 2.0;
  return BA;
}

double MCC(int *Evolved, int DinamicFitnessCases, int *Class) {
  int i;
  double TP, FP, TN, FN;
  double MCC, num, den;
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  MCC = 0;
  for (i = 0; i < DinamicFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 1) {
        ++TP;
      } else {
        ++FP;
      }
    } else {
      if (Class[i] == 0) {
        ++TN;
      } else {
        ++FN;
      }
    }
  }
  num = TP * TN - FP * FN;
  double temp1 = (TP + FP) * (TP + FN);
  double temp2 = (TN + FP) * (TN + FN);
  den = sqrt(temp1 * temp2);
  if (den == 0) {
    MCC = 0;
  } else {
    MCC = num / den;
  }
  return MCC;
}

double Accuracy(int *Evolved, int DinamicFitnessCases, int *Class) {
  int i;
  double TP, FP, TN, FN;
  double ACC, num, den;
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  ACC = 0;
  for (i = 0; i < DinamicFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 1) {
        ++TP;
      } else {
        ++FP;
      }
    } else {
      if (Class[i] == 0) {
        ++TN;
      } else {
        ++FN;
      }
    }
  }
  num = TP + TN;
  den = TP + FP + TN + FN;
  if (den == 0) {
    ACC = 0;
  } else {
    ACC = num / den;
  }
  return ACC;
}

double TPR(int *Evolved, int DinamicValidationFitnessCases, int *Class) {
  int i;
  double TP, FN;
  double TPR;
  TP = 0;
  FN = 0;
  TPR = 0;
  for (i = 0; i < DinamicValidationFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 1) {
        ++TP;
      }
    } else {
      if (Class[i] == 1) {
        ++FN;
      }
    }
  }
  if (TP > 0) {
    TPR = (double)TP / ((double)TP + (double)FN);
  }
  return TPR;
}

double FPR(int *Evolved, int DinamicValidationFitnessCases, int *Class) {
  int i;
  double FP, TN;
  double FPR;
  FP = 0;
  TN = 0;
  FPR = 0;
  for (i = 0; i < DinamicValidationFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 0) {
        ++FP;
      }
    } else {
      if (Class[i] == 0) {
        ++TN;
      }
    }
  }
  if (FP > 0) {
    FPR = (double)FP / ((double)FP + (double)TN);
  }
  return FPR;
}

double *consusion_matrix(int *Evolved, int DinamicValidationFitnessCases,
                         int *Class) {
  int i;
  double TP, FP, TN, FN;
  double *matrix = new double[4];
  TP = 0;
  FP = 0;
  TN = 0;
  FN = 0;
  for (i = 0; i < DinamicValidationFitnessCases; i++) {
    if (Evolved[i] == 1) {
      if (Class[i] == 1) {
        ++TP;
      } else {
        ++FP;
      }
    } else {
      if (Class[i] == 0) {
        ++TN;
      } else {
        ++FN;
      }
    }
  }
  matrix[0] = TP;
  matrix[1] = FP;
  matrix[2] = TN;
  matrix[3] = FN;
  return matrix;
}

double ifLessThan(double exp1, double pivot, double exp2, double exp3) {
  if (exp1 <= pivot)
    return exp2;
  else
    return exp3;
}

double protected_div(double a1, double a2) {
  if (a2 == 0)
    return 0;
  else
    return a1 / a2;
}

double protected_log(double a1) {
  if (a1 <= 0)
    return 0;
  else
    return log(a1);
}

double protected_abs(double a1) { return (a1 > 0.0) ? a1 : ((-1.0) * a1); }

double protected_sqrt(double a1) {
  if (isfinite(sqrt(a1)))
    return sqrt(a1);
  else
    return 0;
}

double protected_exp(double a1) {
  if (isfinite(exp(a1)))
    return exp(a1);
  else
    return 0;
}

double AUC(int *Evolved, int DinamicValidationFitnessCases, int *Class) {
  return AUROC(Class, Evolved, DinamicValidationFitnessCases);
}

int Evolved[ValidationFitnessCases];
int Class[ValidationFitnessCases];

int main(int argc, char *argv[]) {
  std::string filename = "../dataset/cc_dm/TestData.txt";

  FILE *file;
  // double Dim[ValidationFitnessCases][Variables + 1];
  double **Dim = new double *[FitnessCases];
  for (int i = 0; i < FitnessCases; i++) {
    Dim[i] = new double[Variables + 1];
  }
  if (!(file = fopen(filename.c_str(), "r"))) {
    printf("Error opening file");
    exit(0);
  }

  for (int i = 0; i < ValidationFitnessCases; i++) {
    for (int j = 0; j < Variables + 1; j++) {
      fscanf(file, "%lf", &Dim[i][j]);
    }
    Class[i] = Dim[i][Variables];
  }

  fclose(file);

  for (int i = 0; i < ValidationFitnessCases; i++) {
    int sum = 0;
    for (const auto &expr : expressions) {
      double temp = expr(Dim, i);
      if (temp > 0) {
        sum += 1;
      }
    }
    if (sum > 12) {
      Evolved[i] = 1;
    } else {
      Evolved[i] = 0;
    }
  }
  printf("FPR: %f\n", FPR(Evolved, ValidationFitnessCases, Class));
  printf("1.0 - TPR: %f\n", 1.0 - TPR(Evolved, ValidationFitnessCases, Class));
  printf("1.0 - AUC: %f\n", 1.0 - AUC(Evolved, ValidationFitnessCases, Class));

  printf("F1: %f\n", F1(Evolved, ValidationFitnessCases, Class));
  printf("BA: %f\n", BA(Evolved, ValidationFitnessCases, Class));
  printf("MCC: %f\n", MCC(Evolved, ValidationFitnessCases, Class));
  printf("Accuracy: %f\n", Accuracy(Evolved, ValidationFitnessCases, Class));

  double *matrix;
  matrix = consusion_matrix(Evolved, ValidationFitnessCases, Class);
  printf("TP: %f\nFP: %f\nTN: %f\nFN: %f\n", matrix[0], matrix[1], matrix[2],
         matrix[3]);

  //   for (const auto &expr : expressions) {
  //     int i, countC1, countC2;
  //     countC1 = countC2 = 0;
  //     double contrast, variance, entropy;
  //     double boundaryC1, boundaryC2, boundary;
  //     boundaryC1 = boundaryC2 = boundary = 0.0;
  //     double temp[ValidationFitnessCases];
  //     for (i = 0; i < ValidationFitnessCases; i++) {
  //       contrast = Dim[i][0];
  //       variance = Dim[i][1];
  //       entropy = Dim[i][2];
  //       temp[i] = expr(Dim, i);
  //       if (Class[i] == 0) {
  //         boundaryC1 += temp[i];
  //         countC1++;
  //       } else {
  //         boundaryC2 += temp[i];
  //         countC2++;
  //       }
  //     }
  //     boundaryC1 = boundaryC1 / countC1;
  //     boundaryC2 = boundaryC2 / countC2;
  //     boundary = 0; //(boundaryC1 + boundaryC2) / 2;
  //     for (i = 0; i < ValidationFitnessCases; i++) {
  //       if (temp[i] > boundary)
  //         Evolved[i] = 1;
  //       else
  //         Evolved[i] = 0;
  //     }
  //     printf("%f\n", FPR(Evolved, ValidationFitnessCases, Class));
  //     printf("%f\n", 1.0 - TPR(Evolved, ValidationFitnessCases, Class));
  //     printf("%f\n", 1.0 - AUC(Evolved, ValidationFitnessCases, Class));
  //     // printf("FPR: %f\n", FPR(Evolved, ValidationFitnessCases, Class));
  //     // printf("1.0 - TPR: %f\n", 1.0 - TPR(Evolved, ValidationFitnessCases,
  //     // Class)); printf("1.0 - AUC: %f\n", 1.0 - AUC(Evolved,
  //     // ValidationFitnessCases, Class));

  //     // printf("F1: %f\n", F1(Evolved, ValidationFitnessCases, Class));
  //     // printf("BA: %f\n", BA(Evolved, ValidationFitnessCases, Class));
  //     // printf("MCC: %f\n", MCC(Evolved, ValidationFitnessCases, Class));
  //     // printf("Accuracy: %f\n", Accuracy(Evolved, ValidationFitnessCases,
  //     // Class));

  //     // double *matrix;
  //     // matrix = consusion_matrix(Evolved, ValidationFitnessCases, Class);
  //     // printf("TP: %f\nFP: %f\nTN: %f\nFN: %f\n", matrix[0], matrix[1],
  //     // matrix[2],
  //     //       matrix[3]);
  //   }

  for (int i = 0; i < FitnessCases; i++) {
    delete[] Dim[i];
  }
  delete[] Dim;
}
