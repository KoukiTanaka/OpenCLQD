#ifndef INLINE_H
#define INLINE_H
#define HIGH 0
#define LOW 1
#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996
#define NULL 0
#define QD_RAND_MAX 0x07ffffff

double QuickTwoSum(double a, double b, double *err);
double QuickTwoDiff(double a, double b, double *err);
double TwoSum(double a, double b, double *err);
double TwoDiff(double a, double b, double *err);
void split(double a, double *hi, double *lo) ;
double TwoProd(double a, double b, double *err) ;
double TwoSqr(double a, double *err) ;
double nint(double d);
double aint(double d);
void sincosh(double t, double *sinh_t, double *cosh_t);
int OpenCLRand(int *seed);

double QuickTwoSum(double a, double b, double *err){
    double s = a + b;
    *err = b - ( s - a);
    
    return s;
}
double QuickTwoDiff(double a, double b, double *err){
    double s = a - b;
    *err = (a - s) - b;
    
    return s;
}
double TwoSum(double a, double b, double *err){
    double s = a + b;
    double bb = s - a;
    *err = (a - (s - bb)) + (b-bb);
    
    return s;
}

double TwoDiff(double a, double b, double *err){
    double s = a - b;
    double bb = s - a;
    *err = (a - (s - bb)) - (b + bb);
    return s;
}

#ifndef QD_FMS
void split(double a, double *hi, double *lo) {
    double temp;
    if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
        a *= 3.7252902984619140625e-09;  // 2^-28
        temp = _QD_SPLITTER * a;
        *hi = temp - (temp - a);
        *lo = a - *hi;
        *hi *= 268435456.0;          // 2^28
        *lo *= 268435456.0;          // 2^28
    } else {
        temp = _QD_SPLITTER * a;
        *hi = temp - (temp - a);
        *lo = a - *hi;
    }
}
#endif

double TwoProd(double a, double b, double *err) {
#ifdef QD_FMS
    double p = a * b;
    *err = QD_FMS(a, b, p);
    return p;
#else
    double a_hi, a_lo, b_hi, b_lo;
    double p = a * b;
    split(a, &a_hi, &a_lo);
    split(b, &b_hi, &b_lo);
    *err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    return p;
#endif
}

double TwoSqr(double a, double *err) {
#ifdef QD_FMS
    double p = a * a;
    *err = QD_FMS(a, a, p);
    return p;
#else
    double hi, lo;
    double q = a * a;
    split(a, &hi, &lo);
    *err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
    return q;
#endif
}
double nint(double d){
    if( d == floor(d)){
        return d;
    }
    return floor(d+0.5);
}
double aint(double d){
    if( d >= 0.0){
        return floor(d);
    }else{
        return ceil(d);
    }
}
void sincosh(double t, double *sinh_t, double *cosh_t){
    *sinh_t = sinh(t);
    *cosh_t = cosh(t);
}

#define QD_DEF_SEED 2463534242

// OpenCL環境下では乱数関数randの実装はベンダー次第のため仮としての乱数関数
int OpenCLRand(int *seed){
    if( *seed == 0.0){
        *seed = QD_DEF_SEED;
    }
    int y = *seed;
    y ^= (y << 13);
    y ^= (y >> 17);
    y ^= (y << 15);
    *seed = y;
    return y & 0x7fffffff;
}

#endif
