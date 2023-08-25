// cc_cm
// #define FitnessCases 215840
// #define ValidationFitnessCases 215840 // 53960
// #define FitnessCases 67456
// #define ValidationFitnessCases 67456

// cc_dm
// #define FitnessCases 221739
// #define ValidationFitnessCases 55434
// #define FitnessCases 69297
// #define ValidationFitnessCases 69297

// mlo_cm
// #define FitnessCases 
// #define ValidationFitnessCases 
// #define FitnessCases 76606
// #define ValidationFitnessCases 76606

// mlo_dm
// #define FitnessCases 235452
// #define ValidationFitnessCases 58863
// #define FitnessCases 73584
// #define ValidationFitnessCases 73584



#define FitnessCases 1281
#define ValidationFitnessCases 1281

#define Variables 123

double TPR(int * Evolved, int DinamicFitnessCases, int * Class);
double FPR(int * Evolved, int DinamicFitnessCases, int * Class);
double F1(int * Evolved, int DinamicFitnessCases, int * Class);
double BA(int * Evolved, int DinamicFitnessCases, int * Class);
double AUC(int * Evolved, int DinamicFitnessCases, int * Class);
double ifLessThan(double exp1, double pivot, double exp2, double exp3);
double protected_div(double a1, double a2);
double protected_log(double a1);
double protected_abs(double a1);
double protected_sqrt(double a1);
double protected_exp(double a1);
