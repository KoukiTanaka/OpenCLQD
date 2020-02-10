#ifndef CLDD_MAIN
#define CLDD_MAIN
#include "inline.h"
typedef struct{
    double x[2];
}cl_dd;

bool DDEqual(cl_dd a, cl_dd b);
bool DDisOne(cl_dd a);
cl_dd DDabs(cl_dd a);
cl_dd DtoDD(double x_);
double DDtoD(cl_dd a);
cl_dd DD(double x0,double x1);
int DDisZero(cl_dd a);
cl_dd DDldexp(cl_dd a, int exp) ;
cl_dd DDmul_pwr2(cl_dd a, double b) ;
cl_dd DDminus(cl_dd a);
cl_dd DDieeeAdd(cl_dd a, cl_dd b);
cl_dd DDSloopyAdd(cl_dd a, cl_dd b);
cl_dd DaddD(double a,double b);
cl_dd DDadd(cl_dd a, cl_dd b);
cl_dd DDaddD(cl_dd a, double b);
cl_dd DaddDD(double a, cl_dd b);
cl_dd DsubD(double a,double b);
cl_dd DDsub(cl_dd a, cl_dd b);
cl_dd DsubDD(double a,cl_dd b);
cl_dd DDsubD(cl_dd a, double b);
cl_dd DmulD(double a,double b);
cl_dd DDmul(cl_dd a, cl_dd b);
cl_dd DDmulD(cl_dd a, double b);
cl_dd DdivD(double a, double b);
cl_dd DDdivD(cl_dd a, double b);
cl_dd DDSloopyDiv(cl_dd a,cl_dd b);
cl_dd DDAccurateDiv(cl_dd a, cl_dd b);
cl_dd DDdiv(cl_dd a, cl_dd b);
cl_dd DdivDD(double a, cl_dd b);
cl_dd DDinv(cl_dd a);
cl_dd DDnint(cl_dd a) ;
cl_dd DDdrem(cl_dd a, cl_dd b) ;
cl_dd DDdivrem(cl_dd a, cl_dd b, cl_dd *r) ;
cl_dd DDfloor(cl_dd a);
cl_dd DDceil(cl_dd a);
cl_dd DDaint(cl_dd a);
cl_dd DDsqr(cl_dd a) ;
cl_dd DDsqrD(double a) ;

bool DDEqual(cl_dd a, cl_dd b){
    if( a.x[0] == b.x[0] && a.x[1] == b.x[1]){
        return true;
    }else{
        return false;
    }
}

int DDisZero(cl_dd a){
    if( a.x[0] == 0.0){
        return 1;
    }
    return 0;
}

bool DDisOne(cl_dd a){
    if(a.x[0] == 1.0 && a.x[1] == 0){
        return true;
    }else{
        return false;
    }
   
}

cl_dd DDabs(cl_dd a){
    if( a.x[0] < 0){
        a.x[0] *= -1.0f;
        a.x[1] *= -1.0f;
    }
    
    return a;
}

cl_dd DDminus(cl_dd a){
    int i;
    a.x[0] *= -1;
    a.x[1] *= -1;
    
    return a;
}

cl_dd DD(double x0,double x1){
    cl_dd a;
    a.x[0] = x0;
    a.x[1] = x1;
    
    return a;
}

cl_dd DtoDD(double x_){
    cl_dd x;
    x.x[0] = x_;
    x.x[1] = 0.0;
    return x;
}

double DDtoD(cl_dd a){
    return a.x[0];
}


cl_dd DDldexp(cl_dd a, int exp) {
    return DD(ldexp(a.x[0], exp), ldexp(a.x[1], exp));
}

cl_dd DDmul_pwr2(cl_dd a, double b) {
    return DD(a.x[0] * b, a.x[1] * b);
}

cl_dd DaddD(double a,double b){
    double s,e;
    s = TwoSum(a,b,&e);
    return DD(s,e);
}

cl_dd DDaddD(cl_dd a, double b){
    double s1, s2;
    s1 = TwoSum(a.x[0], b, &s2);
    s2 += a.x[1];
    s1 = QuickTwoSum(s1, s2, &s2);
    return DD(s1, s2);
}

cl_dd DaddDD(double a, cl_dd b){
    return DDaddD(b,a);
}


cl_dd DDieeeAdd(cl_dd a, cl_dd b){
    double s1,s2,t1,t2;
    
    s1 = TwoSum(a.x[0], b.x[0], &s2);
    t1 = TwoSum(a.x[1], b.x[1], &t2);
    s2 += t1;
    s1 = QuickTwoSum(s1, s2, &s2);
    s2 += t2;
    s1 = QuickTwoSum(s1, s2, &s2);
    return DD(s1, s2);
}

cl_dd DDSloopyAdd(cl_dd a, cl_dd b){
    double s, e;
    
    s = TwoSum(a.x[0], b.x[0], &e);
    e += (a.x[1] + b.x[1]);
    s = QuickTwoSum(s, e, &e);
    return DD(s, e);
}

cl_dd DDadd(cl_dd a, cl_dd b){
#ifndef QD_IEEE_ADD
    return DDSloopyAdd(a,b);
#else
    return DDieeeAdd(a,b);
#endif
}

cl_dd DsubD(double a,double b){
    double s, e;
    s = TwoDiff(a, b, &e);
    return DD(s, e);
}

cl_dd DsubDD(double a,cl_dd b){
    double s1, s2;
    s1 = TwoDiff(a, b.x[0], &s2);
    s2 -= b.x[1];
    s1 = QuickTwoSum(s1, s2, &s2);
    return DD(s1, s2);
}
cl_dd DDsubD(cl_dd a, double b){
    return DDadd(a,DtoDD(-b));
}

cl_dd DDsub(cl_dd a, cl_dd b){
#ifndef QD_IEEE_ADD
    double s, e;
    s = TwoDiff(a.x[0], b.x[0], &e);
    e += a.x[1];
    e -= b.x[1];
    s = QuickTwoSum(s, e, &e);
    return DD(s, e);
#else
    double s1,s2,t1,t2;
    s1 = TwoDiff(a.x[0],b.x[0],&s2);
    t1 = TwoDiff(a.x[1],b.x[1],&t2);
    s2 += t1;
    s1 = QuickTwoSum(s1,s2,&s2);
    s2 += t2;
    s1 = QuickTwoSum(s1,s2,&s2);
    return DD(s1,s2);
#endif
}


cl_dd DmulD(double a,double b){
    double p,e;
    p = TwoProd(a,b,&e);
    return DD(p,e);
}

cl_dd DDmulD(cl_dd a, double b){
    double p1, p2;
    
    p1 = TwoProd(a.x[0], b, &p2);
    p2 += (a.x[1] * b);
    p1 = QuickTwoSum(p1, p2, &p2);
    return DD(p1, p2);
}

cl_dd DmulDD(double a, cl_dd b){
    return DDmulD(b,a);
}

cl_dd DDmul(cl_dd a, cl_dd b){
    double p1, p2;
    
    p1 = TwoProd(a.x[0], b.x[0], &p2);
    p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
    p1 = QuickTwoSum(p1, p2, &p2);
    return DD(p1, p2);
}


cl_dd DdivD(double a, double b){
    double q1, q2;
    double p1, p2;
    double s, e;
    
    q1 = a / b;
    
    /* Compute  a - q1 * b */
    p1 = TwoProd(q1, b, &p2);
    s = TwoDiff(a, p1, &e);
    e -= p2;
    
    /* get next approximation */
    q2 = (s + e) / b;
    
    s = QuickTwoSum(q1, q2, &e);
    
    return DD(s, e);
}

cl_dd DDdivD(cl_dd a, double b){
    double q1, q2;
    double p1, p2;
    double s, e;
    cl_dd r;
    
    q1 = a.x[0] / b;   /* approximate quotient. */
    
    /* Compute  this - q1 * d */
    p1 = TwoProd(q1, b, &p2);
    s = TwoDiff(a.x[0], p1, &e);
    e += a.x[1];
    e -= p2;
    
    /* get next approximation. */
    q2 = (s + e) / b;
    
    /* renormalize */
    r.x[0] = QuickTwoSum(q1, q2, &r.x[1]);
    
    return r;
}

cl_dd DdivDD(double a, cl_dd b){
    return DDdiv(DtoDD(a),b);
}

cl_dd DDSloopyDiv(cl_dd a,cl_dd b){
    double s1, s2;
    double q1, q2;
    cl_dd r;
    
    q1 = a.x[0] / b.x[0];  /* approximate quotient */
    
    /* compute  this - q1 * dd */
    r = DDmulD(b,q1);
    s1 = TwoDiff(a.x[0], r.x[0], &s2);
    s2 -= r.x[1];
    s2 += a.x[1];
    
    /* get next approximation */
    q2 = (s1 + s2) / b.x[0];
    
    /* renormalize */
    r.x[0] = QuickTwoSum(q1, q2, &r.x[1]);
    return r;
}

cl_dd DDAccurateDiv(cl_dd a, cl_dd b){
    double q1, q2, q3;
    cl_dd r;
    
    q1 = a.x[0] / b.x[0];  /* approximate quotient */
    
    r = DDsub(a,DDmulD(b,q1));
    
    q2 = r.x[0] / b.x[0];
    r = DDsub(r,DDmulD(b,q2));
    
    q3 = r.x[0] / b.x[0];
    
    q1 = QuickTwoSum(q1, q2, &q2);
    r = DDaddD(DD(q1, q2),q3);
    return r;
}

cl_dd DDdiv(cl_dd a, cl_dd b){
#ifdef QD_SLOPPY_DIV
    return DDSloopyDiv(a,b);
#else
    return DDAccurateDiv(a,b);
#endif
}

cl_dd DDinv(cl_dd a){
    return DdivDD(1.0,a);
}

cl_dd DDnint(cl_dd a) {
  double hi = nint(a.x[0]);
  double lo;

  if (hi == a.x[0]) {
    /* High word is an integer already.  Round the low word.*/
    lo = nint(a.x[1]);
    
    /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
    hi = QuickTwoSum(hi, lo,&lo);
    hi = QuickTwoSum(hi,lo,&lo);
  } else {
    /* High word is not an integer. */
    lo = 0.0;
    if (fabs(hi-a.x[0]) == 0.5 && a.x[1] < 0.0) {
      /* There is a tie in the high word, consult the low word 
         to break the tie. */
      hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
    }
  }

  return DD(hi, lo);
}

cl_dd DDdrem(cl_dd a, cl_dd b) {
  cl_dd n = DDnint(DDdiv(a,b));
  return DDsub(a,DDmul(n,b));
}

cl_dd DDdivrem(cl_dd a, cl_dd b, cl_dd *r) {
  cl_dd n = DDnint(DDdiv(a,b));
  *r = DDsub(a,DDmul(n,b));
  return n;
}

cl_dd DDfloor(cl_dd a){
    double hi = floor(a.x[0]);
    double lo = 0.0;

    if( hi  == a.x[0] ){
        lo = floor(a.x[1]);
        hi = QuickTwoSum(hi,lo,&lo);
    }
    return DD(hi,lo);
}

cl_dd DDceil(cl_dd a){
    double hi = ceil(a.x[0]);
    double lo = 0.0;

    if( hi  == a.x[0] ){
        lo = ceil(a.x[1]);
        hi = QuickTwoSum(hi,lo,&lo);
    }
    return DD(hi,lo);
}

cl_dd DDaint(cl_dd a){
    if( a.x[0] >= 0.0){
        return DDfloor(a);
    }else{
        return DDceil(a);
    }
}

cl_dd DDsqr(cl_dd a) {
  double p1, p2;
  double s1, s2;
  p1 = TwoSqr(a.x[0], &p2);
  p2 += 2.0 * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  s1 = QuickTwoSum(p1, p2, &s2);
  return DD(s1, s2);
}

cl_dd DDsqrD(double a) {
  double p1, p2;
  p1 = TwoSqr(a, &p2);
  return DD(p1, p2);
}

#endif
