#include "dd_inline.h"
#ifndef CLQD_INLINE
#define CLQD_INLINE


typedef struct{
    double x[4];
}cl_qd;

bool QDEqual(cl_qd a,cl_qd b);
cl_qd QDabs(cl_qd a);
cl_qd QDminus(cl_qd a);
int QDisZero(cl_qd a);
bool QDisOne(cl_qd a);
cl_qd DtoQD(double x_);
cl_qd DDtoQD(cl_dd a);
double QDtoD(cl_qd a);
cl_dd QDtoDD(cl_qd a);
cl_qd QD(double x0,double x1, double x2,double x3);
cl_qd QDmul_pwr2( cl_qd a, double b) ;
cl_qd QDldexp(cl_qd a, int n);
void QuickRenorm(double *c0,double *c1,double *c2,double *c3,double *c4);
void renormD4(double *c0,double *c1,double *c2,double *c3);
void renormD5(double *c0,double *c1,double *c2,double *c3,double *c4);
void renormQD(cl_qd *a);
void renormQDD(cl_qd *a,double *e);
void ThreeSum(double *a,double *b, double *c);
void ThreeSum2(double *a,double *b,double *c);
cl_qd QDaddD(cl_qd a,double b);
cl_qd QDaddDD(cl_qd a, cl_dd b);
cl_qd DaddQD(double a, cl_qd b);
cl_qd DDaddQD(cl_dd a, cl_qd b);
double QuickThreeAccum(double *a,double *b,double c);
cl_qd QDieeeAdd(cl_qd a, cl_qd b);
cl_qd QDSloopyAdd(cl_qd a, cl_qd b);
cl_qd QDadd(cl_qd a, cl_qd b);
cl_qd QDsubD(cl_qd a, double b);
cl_qd QDsubDD(cl_qd a, cl_dd b);
cl_qd QDsub(cl_qd a, cl_qd b);
cl_qd DsubQD(double a,cl_qd b);
cl_qd DDsubQD(cl_dd a, cl_qd b);
cl_qd QDmulD(cl_qd a, double b);
cl_qd DmulQD(double a, cl_qd b);
cl_qd QDmulDD(cl_qd a, cl_dd b);
cl_qd DDmulDD(cl_dd a, cl_qd b);
cl_qd QDSloppyMul(cl_qd a, cl_qd b);
cl_qd QDAccurateMul(cl_qd a,cl_qd b);
cl_qd QDmul(cl_qd a, cl_qd b);
cl_qd QDsqr(cl_qd a);
cl_qd QDdivD(cl_qd a,double b);
cl_qd QDSloopyDivDD(cl_qd a, cl_dd b);
cl_qd QDAccurateDivDD(cl_qd a, cl_dd b);
cl_qd QDdivDD(cl_qd a,cl_dd b);
cl_qd QDSloppyDiv(cl_qd a, cl_qd b);
cl_qd QDAccurateDiv(cl_qd a, cl_qd b);
cl_qd QDdiv(cl_qd a,cl_qd b);
cl_qd DdivQD(double a, cl_qd b);
cl_qd DDdidQD(cl_dd a,cl_qd b);
cl_qd QDinv(cl_qd a);
cl_qd QDquick_nint(cl_qd a);
cl_qd QDfloor(cl_qd a) ;
cl_qd QDceil(cl_qd a) ;
cl_qd QDaint(cl_qd a) ;

bool QDEqual(cl_qd a,cl_qd b){
    if( a.x[0] == b.x[0] && a.x[1] == b.x[1] && a.x[2] == b.x[2] && a.x[3] == b.x[3]){
        return true;
    }
    return false;
}

int QDisZero(cl_qd a){
    if( a.x[0] == 0 ){
        return 1;
    }
    
    return 0;
}
bool QDisOne(cl_qd a){
    if(a.x[0] == 1.0f && a.x[1] == 0.0 && a.x[2] == 0.0 && a.x[3] == 0.0){
        return true;
    }

    return false;
}

cl_qd QDabs(cl_qd a){
    if( a.x[0] < 0){
        a.x[0] *= -1.0f;
        a.x[1] *= -1.0f;
        a.x[2] *= -1.0f;
        a.x[3] *= -1.0f;
    }
    
    return a;
}
cl_qd QDminus(cl_qd a){
    a.x[0] *= -1;
    a.x[1] *= -1;
    a.x[2] *= -1;
    a.x[3] *= -1;
    
    return a;
}

cl_qd DtoQD(double x_){
    cl_qd x;
    x.x[0] = x_;
    x.x[1] = x.x[2] = x.x[3] = 0.0;
    return x;
}

cl_qd DDtoQD(cl_dd a){
    cl_qd x;
    x.x[0] = a.x[0];
    x.x[1] = a.x[1];
    x.x[2] = x.x[3] = 0.0;
    
    return x;
}

double QDtoD(cl_qd a){
    return a.x[0];
}

cl_dd QDtoDD(cl_qd a){
    return DD(a.x[0],a.x[1]);
}

cl_qd QD(double x0,double x1, double x2,double x3){
    cl_qd x;
    x.x[0] = x0;
    x.x[1] = x1;
    x.x[2] = x2;
    x.x[3] = x3;
    
    return x;
}

cl_qd QDldexp(cl_qd a, int n){
    return QD(ldexp(a.x[0],n),ldexp(a.x[1],n),ldexp(a.x[2],n),ldexp(a.x[3],n));
}

cl_qd QDmul_pwr2( cl_qd a, double b) {
    return QD(a.x[0] * b, a.x[1] * b, a.x[2] * b, a.x[3] * b);
}

void QuickRenorm(double *c0,double *c1,double *c2,double *c3,double *c4){
    double t0,t1,t2,t3;
    double s;
    
    s  = QuickTwoSum(*c3, *c4, &t3);
    s  = QuickTwoSum(*c2, s , &t2);
    s  = QuickTwoSum(*c1, s , &t1);
    *c0 = QuickTwoSum(*c0, s , &t0);
    
    s  = QuickTwoSum(t2, t3, &t2);
    s  = QuickTwoSum(t1, s , &t1);
    *c1 = QuickTwoSum(t0, s , &t0);
    
    s  = QuickTwoSum(t1, t2, &t1);
    *c2 = QuickTwoSum(t0, s , &t0);
    
    *c3 = t0 + t1;
}

void renormD4(double *c0,double *c1,double *c2,double *c3){
    double s0,s1,s2=0.0,s3=0.0;
    
    s0 = QuickTwoSum(*c2, *c3, c3);
    s0 = QuickTwoSum(*c1, s0, c2);
    *c0 = QuickTwoSum(*c0, s0, c1);
    
    s0 = *c0;
    s1 = *c1;
    if (s1 != 0.0) {
        s1 = QuickTwoSum(s1, *c2, &s2);
        if (s2 != 0.0)
            s2 = QuickTwoSum(s2, *c3, &s3);
        else
            s1 = QuickTwoSum(s1, *c3, &s2);
    } else {
        s0 = QuickTwoSum(s0, *c2, &s1);
        if (s1 != 0.0)
            s1 = QuickTwoSum(s1, *c3, &s2);
        else
            s0 = QuickTwoSum(s0, *c3, &s1);
    }
    
    *c0 = s0;
    *c1 = s1;
    *c2 = s2;
    *c3 = s3;
}

void renormD5(double *c0,double *c1,double *c2,double *c3,double *c4){
    double s0, s1, s2 = 0.0, s3 = 0.0;
    
    
    s0 = QuickTwoSum(*c3, *c4, c4);
    s0 = QuickTwoSum(*c2, s0, c3);
    s0 = QuickTwoSum(*c1, s0, c2);
    *c0 = QuickTwoSum(*c0, s0, c1);
    
    s0 = *c0;
    s1 = *c1;
    
    if (s1 != 0.0) {
        s1 = QuickTwoSum(s1, *c2, &s2);
        if (s2 != 0.0) {
            s2 = QuickTwoSum(s2, *c3,&s3);
            if (s3 != 0.0)
                s3 += *c4;
            else
                s2 = QuickTwoSum(s2, *c4, &s3);
        } else {
            s1 = QuickTwoSum(s1, *c3, &s2);
            if (s2 != 0.0)
                s2 = QuickTwoSum(s2, *c4, &s3);
            else
                s1 = QuickTwoSum(s1, *c4, &s2);
        }
    } else {
        s0 = QuickTwoSum(s0, *c2, &s1);
        if (s1 != 0.0) {
            s1 = QuickTwoSum(s1, *c3, &s2);
            if (s2 != 0.0)
                s2 = QuickTwoSum(s2, *c4, &s3);
            else
                s1 = QuickTwoSum(s1, *c4, &s2);
        } else {
            s0 = QuickTwoSum(s0, *c3, &s1);
            if (s1 != 0.0)
                s1 = QuickTwoSum(s1, *c4, &s2);
            else
                s0 = QuickTwoSum(s0, *c4, &s1);
        }
    }
    
    *c0 = s0;
    *c1 = s1;
    *c2 = s2;
    *c3 = s3;
}

void renormQD(cl_qd *a){
    renormD4(&(a->x[0]),&(a->x[1]),&(a->x[2]),&(a->x[3]));
}

void renormQDD(cl_qd *a,double *e){
    renormD5(&(a->x[0]),&(a->x[1]),&(a->x[2]),&(a->x[3]),e);
}

void ThreeSum(double *a,double *b, double *c){
    double t1,t2,t3;
    t1 = TwoSum(*a,*b,&t2);
    *a = TwoSum(*c,t1,&t3);
    *b = TwoSum(t2,t3,c);
}
void ThreeSum2(double *a,double *b,double *c){
    double t1,t2,t3;
    t1 = TwoSum(*a,*b,&t2);
    *a = TwoSum(*c,t1,&t3);
    *b = t2 + t3;
}
cl_qd QDaddD(cl_qd a,double b){
    double c0,c1,c2,c3;
    double e;
    
    c0 = TwoSum(a.x[0],b,&e);
    c1 = TwoSum(a.x[1],e,&e);
    c2 = TwoSum(a.x[2],e,&e);
    c3 = TwoSum(a.x[3],e,&e);
    
    renormD5(&c0,&c1,&c2,&c3,&e);
    return QD(c0,c1,c2,c3);
}

cl_qd DaddQD(double a, cl_qd b){
    return QDaddD(b,a);
}

cl_qd QDaddDD(cl_qd a, cl_dd b){
    double s0,s1,s2,s3;
    double t0,t1;
    
    s0 = TwoSum(a.x[0],b.x[0],&t0);
    s1 = TwoSum(a.x[1],b.x[1],&t1);
    
    s1 = TwoSum(s1,t0,&t0);
    s2 = a.x[2];
    ThreeSum(&s2,&t0,&t1);
    s3 = TwoSum(t0,a.x[3],&t0);
    t0 += t1;
    
    renormD5(&s0,&s1,&s2,&s3,&t0);
    return QD(s0,s1,s2,s3);
}

cl_qd DDaddQD(cl_dd a, cl_qd b){
    return QDaddDD(b,a);
}

double QuickThreeAccum(double *a,double *b,double c){
    double s;
    int za,zb;
    
    s = TwoSum(*b,c,b);
    s = TwoSum(*a,s,a);
    
    za = *a != 0.0;
    zb = *b != 0.0;
    
    if( za && zb )
        return s;
    
    if(!zb){
        *b = *a;
        *a = s;
    }else{
        *a = s;
    }
    
    return 0.0;
}

cl_qd QDieeeAdd(cl_qd a, cl_qd b){
    int i,j,k;
    double s,t;
    double u,v;
    double x[4] = {0.0,0.0,0.0,0.0};
    
    
    i = j = k = 0;
    if( fabs(a.x[i]) > fabs(b.x[j]))
        u = a.x[i++];
    else
        u = b.x[j++];
    if( fabs(a.x[i]) > fabs(b.x[j]))
        v = a.x[i++];
    else
        v = b.x[j++];
    
    u = QuickTwoSum(u,v,&v);
    
    while(k < 4){
        if( i >= 4 && j >= 4){
            x[k] = u;
            if( k < 3)
                x[++k] = v;
            break;
        }
        
        if( i >= 4){
            t = b.x[j++];
        }else if(j >= 4){
            t = a.x[i++];
        }else if( fabs(a.x[i]) > fabs(b.x[j]) ){
            t = a.x[i++];
        }else{
            t = b.x[j++];
        }
        
        s = QuickThreeAccum(&u,&v,t);
        if( s != 0.0){
            x[k++] = s;
        }
        
    }
    
    for(k = i; k < 4;k++)
        x[3] += a.x[k];
    for(k = j; k < 4;k++)
        x[3] += b.x[k];
    
    renormD4(&x[0],&x[1],&x[2],&x[3]);
    return QD(x[0],x[1],x[2],x[3]);

}

cl_qd QDSloopyAdd(cl_qd a, cl_qd b){
    double s0, s1, s2, s3;
    double t0, t1, t2, t3;
    
    double v0, v1, v2, v3;
    double u0, u1, u2, u3;
    double w0, w1, w2, w3;
    
    s0 = a.x[0] + b.x[0];
    s1 = a.x[1] + b.x[1];
    s2 = a.x[2] + b.x[2];
    s3 = a.x[3] + b.x[3];
    
    v0 = s0 - a.x[0];
    v1 = s1 - a.x[1];
    v2 = s2 - a.x[2];
    v3 = s3 - a.x[3];
    
    u0 = s0 - v0;
    u1 = s1 - v1;
    u2 = s2 - v2;
    u3 = s3 - v3;
    
    w0 = a.x[0] - u0;
    w1 = a.x[1] - u1;
    w2 = a.x[2] - u2;
    w3 = a.x[3] - u3;
    
    u0 = b.x[0] - v0;
    u1 = b.x[1] - v1;
    u2 = b.x[2] - v2;
    u3 = b.x[3] - v3;
    
    t0 = w0 + u0;
    t1 = w1 + u1;
    t2 = w2 + u2;
    t3 = w3 + u3;
    
    s1 = TwoSum(s1, t0, &t0);
    ThreeSum(&s2, &t0, &t1);
    ThreeSum2(&s3, &t0, &t2);
    t0 = t0 + t1 + t3;
    
    // renormalize
    renormD5(&s0, &s1, &s2, &s3, &t0);
    return QD(s0, s1, s2, s3);
}

// QD + QD
cl_qd QDadd(cl_qd a, cl_qd b){
#ifndef QD_IEEE_ADD
    return QDSloopyAdd(a,b);
#elif
    return QDieeeAdd(a,b);
#endif
}

// QD - D
cl_qd QDsubD(cl_qd a, double b){
    return QDaddD(a,-b);
}

cl_qd DsubQD(double a,cl_qd b){
    return DaddQD(a,QDminus(b));
}

cl_qd QDsubDD(cl_qd a, cl_dd b){
    return QDaddDD(a,DDminus(b));
}

cl_qd DDsubQD(cl_dd a, cl_qd b){
    return DDaddQD(a,QDminus(b));
}

cl_qd QDsub(cl_qd a, cl_qd b){
    return QDadd(a,QDminus(b));
}

cl_qd QDmulD(cl_qd a, double b){
    double p0, p1, p2, p3;
    double q0, q1, q2;
    double s0, s1, s2, s3, s4;
    
    p0 = TwoProd(a.x[0], b, &q0);
    p1 = TwoProd(a.x[1], b, &q1);
    p2 = TwoProd(a.x[2], b, &q2);
    p3 = a.x[3] * b;
    
    s0 = p0;
    
    s1 = TwoSum(q0, p1, &s2);
    
    ThreeSum(&s2, &q1, &p2);
    ThreeSum2(&q1, &q2, &p3);
    s3 = q1;
    
    s4 = q2 + p2;
    
    renormD5(&s0, &s1, &s2, &s3, &s4);
    return QD(s0, s1, s2, s3);
}

cl_qd DmulQD(double a, cl_qd b){
    return QDmulD(b,a);
}

cl_qd QDmulDD(cl_qd a, cl_dd b){
    double p0, p1, p2, p3, p4;
    double q0, q1, q2, q3, q4;
    double s0, s1, s2;
    double t0, t1;
    
    p0 = TwoProd(a.x[0], b.x[0], &q0);
    p1 = TwoProd(a.x[0], b.x[1], &q1);
    p2 = TwoProd(a.x[1], b.x[0], &q2);
    p3 = TwoProd(a.x[1], b.x[1], &q3);
    p4 = TwoProd(a.x[2], b.x[0], &q4);
    
    ThreeSum(&p1, &p2, &q0);
    
    // Five-Three-Sum
    ThreeSum(&p2, &p3, &p4);
    q1 = TwoSum(q1, q2, &q2);
    s0 = TwoSum(p2, q1, &t0);
    s1 = TwoSum(p3, q2, &t1);
    s1 = TwoSum(s1, t0, &t0);
    s2 = t0 + t1 + p4;
    p2 = s0;
    
    p3 = a.x[2] * b.x[0] + a.x[3] * b.x[1] + q3 + q4;
    ThreeSum2(&p3, &q0, &s1);
    p4 = q0 + s2;
    
    renormD5(&p0, &p1, &p2, &p3, &p4);
    return QD(p0, p1, p2, p3);
}

cl_qd DDmulQD(cl_dd a, cl_qd b){
    return QDmulDD(b,a);
}

cl_qd QDSloppyMul(cl_qd a, cl_qd b){
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double t0, t1;
    double s0, s1, s2;
    
    p0 = TwoProd(a.x[0], b.x[0], &q0);
    
    p1 = TwoProd(a.x[0], b.x[1], &q1);
    p2 = TwoProd(a.x[1], b.x[0], &q2);
    
    p3 = TwoProd(a.x[0], b.x[2], &q3);
    p4 = TwoProd(a.x[1], b.x[1], &q4);
    p5 = TwoProd(a.x[2], b.x[0], &q5);
    
    // Start Accumulation
    ThreeSum(&p1, &p2, &q0);
    
    // Six-Three Sum  of p2, q1, q2, p3, p4, p5.
    ThreeSum(&p2, &q1, &q2);
    ThreeSum(&p3, &p4, &p5);
    // compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5).
    s0 = TwoSum(p2, p3, &t0);
    s1 = TwoSum(q1, p4, &t1);
    s2 = q2 + p5;
    s1 = TwoSum(s1, t0, &t0);
    s2 += (t0 + t1);
    
    // O(eps^3) order terms
    s1 += a.x[0]*b.x[3] + a.x[1]*b.x[2] + a.x[2]*b.x[1] + a.x[3]*b.x[0] + q0 + q3 + q4 + q5;
    renormD5(&p0, &p1, &s0, &s1, &s2);
    return QD(p0, p1, s0, s1);
}

cl_qd QDAccurateMul(cl_qd a,cl_qd b){
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double p6, p7, p8, p9;
    double q6, q7, q8, q9;
    double r0, r1;
    double t0, t1;
    double s0, s1, s2;
    
    p0 = TwoProd(a.x[0], b.x[0], &q0);
    
    p1 = TwoProd(a.x[0], b.x[1], &q1);
    p2 = TwoProd(a.x[1], b.x[0], &q2);
    
    p3 = TwoProd(a.x[0], b.x[2], &q3);
    p4 = TwoProd(a.x[1], b.x[1], &q4);
    p5 = TwoProd(a.x[2], b.x[0], &q5);
    
    // Start Accumulation
    ThreeSum(&p1, &p2, &q0);
    
    // Six-Three Sum  of p2, q1, q2, p3, p4, p5.
    ThreeSum(&p2, &q1, &q2);
    ThreeSum(&p3, &p4, &p5);
    // compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5).
    s0 = TwoSum(p2, p3, &t0);
    s1 = TwoSum(q1, p4, &t1);
    s2 = q2 + p5;
    s1 = TwoSum(s1, t0, &t0);
    s2 += (t0 + t1);
    
    // O(eps^3) order terms
    p6 = TwoProd(a.x[0], b.x[3], &q6);
    p7 = TwoProd(a.x[1], b.x[2], &q7);
    p8 = TwoProd(a.x[2], b.x[1], &q8);
    p9 = TwoProd(a.x[3], b.x[0], &q9);
    
    // Nine-Two-Sum of q0, s1, q3, q4, q5, p6, p7, p8, p9.
    q0 = TwoSum(q0, q3, &q3);
    q4 = TwoSum(q4, q5, &q5);
    p6 = TwoSum(p6, p7, &p7);
    p8 = TwoSum(p8, p9, &p9);
    // Compute (t0, t1) = (q0, q3) + (q4, q5).
    t0 = TwoSum(q0, q4, &t1);
    t1 += (q3 + q5);
    // Compute (r0, r1) = (p6, p7) + (p8, p9).
    r0 = TwoSum(p6, p8, &r1);
    r1 += (p7 + p9);
    // Compute (q3, q4) = (t0, t1) + (r0, r1).
    q3 = TwoSum(t0, r0, &q4);
    q4 += (t1 + r1);
    // Compute (t0, t1) = (q3, q4) + s1.
    t0 = TwoSum(q3, s1, &t1);
    t1 += q4;
    
    // O(eps^4) terms -- Nine-One-Sum
    t1 += a.x[1] * b.x[3] + a.x[2] * b.x[2] + a.x[3] * b.x[1] + q6 + q7 + q8 + q9 + s2;
    
    renormD5(&p0, &p1, &s0, &t0, &t1);
    return QD(p0, p1, s0, t0);
}

cl_qd QDmul(cl_qd a, cl_qd b){
#ifdef QD_SLOPPY_MUL
    return QDSloppyMul(a,b);
#else
    return QDAccurateMul(a,b);
#endif
}

cl_qd QDsqr(cl_qd a){
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3;
    double s0, s1;
    double t0, t1;
    
    p0 = TwoSqr(a.x[0], &q0);
    p1 = TwoProd(2.0 * a.x[0], a.x[1], &q1);
    p2 = TwoProd(2.0 * a.x[0], a.x[2], &q2);
    p3 = TwoSqr(a.x[1], &q3);
    
    p1 = TwoSum(q0, p1, &q0);
    
    q0 = TwoSum(q0, q1, &q1);
    p2 = TwoSum(p2, p3, &p3);
    
    s0 = TwoSum(q0, p2, &t0);
    s1 = TwoSum(q1, p3, &t1);
    
    s1 = TwoSum(s1, t0, &t0);
    t0 += t1;
    
    s1 = QuickTwoSum(s1, t0, &t0);
    p2 = QuickTwoSum(s0, s1, &t1);
    p3 = QuickTwoSum(t1, t0, &q0);
    
    p4 = 2.0 * a.x[0] * a.x[3];
    p5 = 2.0 * a.x[1] * a.x[2];
    
    p4 = TwoSum(p4, p5, &p5);
    q2 = TwoSum(q2, q3, &q3);
    
    t0 = TwoSum(p4, q2, &t1);
    t1 = t1 + p5 + q3;
    
    p3 = TwoSum(p3, t0, &p4);
    p4 = p4 + q0 + t1;
    
    renormD5(&p0, &p1, &p2, &p3, &p4);
    return QD(p0, p1, p2, p3);
}

cl_qd QDdivD(cl_qd a,double b){
    double t0, t1;
    double q0, q1, q2, q3;
    cl_qd r;
    
    q0 = a.x[0] / b;  /* approximate quotient */
    
    /* Compute the remainder  a - q0 * b */
    t0 = TwoProd(q0, b, &t1);
    r = QDsubDD(a ,DD(t0, t1));
    
    /* Compute the first correction */
    q1 = r.x[0] / b;
    t0 = TwoProd(q1, b, &t1);
    r  = QDsubDD(r,DD(t0, t1));
    
    /* Second correction to the quotient. */
    q2 = r.x[0] / b;
    t0 = TwoProd(q2, b, &t1);
    r  = QDsubDD(r,DD(t0, t1));
    
    /* Final correction to the quotient. */
    q3 = r.x[0] / b;
    
    renormD4(&q0, &q1, &q2, &q3);
    return QD(q0, q1, q2, q3);
}

cl_qd DdivQD(double a, cl_qd b){
    return QDdiv(DtoQD(a),b);
}

cl_qd QDSloopyDivDD(cl_qd a, cl_dd b){
    double q0, q1, q2, q3;
    cl_qd r;
    cl_qd qd_b = DDtoQD(b);
    
    q0 = a.x[0] / b.x[0];
    r = QDsub(a,DmulQD(q0,qd_b));
    
    q1 = r.x[0] / b.x[0];
    r  = QDsub(r,DmulQD(q1,qd_b));
    
    q2 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q2,qd_b));
    
    q3 = r.x[0] / b.x[0];
    
    renormD4(&q0, &q1, &q2, &q3);
    return QD(q0, q1, q2, q3);
}
cl_qd QDAccurateDivDD(cl_qd a, cl_dd b){
    double q0, q1, q2, q3, q4;
    cl_qd r;
    cl_qd qd_b = DDtoQD(b);
    
    q0 = a.x[0] / b.x[0];
    r = QDsub(a,DmulQD(q0,qd_b));
    
    q1 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q1,qd_b));
    
    q2 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q2,qd_b));
    
    q3 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q3,qd_b));
    
    q4 = r.x[0] / b.x[0];
    
    renormD5(&q0, &q1, &q2, &q3, &q4);
    return QD(q0, q1, q2, q3);
}

cl_qd QDdivDD(cl_qd a,cl_dd b){
#ifdef QD_SLOPPY_DIV
    return QDSloopyDivDD(a,b);
#else
    return QDAccurateDivDD(a,b);
#endif
}

cl_qd DDdidQD(cl_dd a,cl_qd b){
    return QDdiv(DDtoQD(a),b);
}

cl_qd QDSloppyDiv(cl_qd a, cl_qd b){

    double q0, q1, q2, q3;
    
    cl_qd r;
    
    q0 = a.x[0] / b.x[0];
    r = QDsub(a,DmulQD(q0,b));
    
    q1 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q1,b));
    
    q2 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q2,b));
    
    q3 = r.x[0] / b.x[0];
    
    renormD4(&q0, &q1, &q2, &q3);
    
    return QD(q0, q1, q2, q3);
}
cl_qd QDAccurateDiv(cl_qd a, cl_qd b){
    double q0, q1, q2, q3,q4;
    
    cl_qd r;
    
    q0 = a.x[0] / b.x[0];
    r = QDsub(a,DmulQD(q0,b));
    
    q1 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q1,b));
    
    q2 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q2,b));
    
    q3 = r.x[0] / b.x[0];
    r = QDsub(r,DmulQD(q3,b));
    
    q4 = r.x[0] / b.x[0];
    
    renormD5(&q0, &q1, &q2, &q3,&q4);
    
    return QD(q0, q1, q2, q3);
}
cl_qd QDdiv(cl_qd a,cl_qd b){
#ifdef QD_SLOPPY_DIV
    return QDSloppyDiv(a,b);
#else
    return QDAccurateDiv(a,b);
#endif
}

cl_qd QDinv(cl_qd a){
    return DdivQD(1.0,a);
}

cl_qd QDquick_nint(cl_qd a){
    cl_qd r = QD(nint(a.x[0]),nint(a.x[1]),nint(a.x[2]),nint(a.x[3]));
    renormQD(&r);
    return r;
}

cl_qd QDfloor(cl_qd a) {
  double x0, x1, x2, x3;
  x1 = x2 = x3 = 0.0;
  x0 = floor(a.x[0]);

  if (x0 == a.x[0]) {
    x1 = floor(a.x[1]);
    
    if (x1 == a.x[1]) {
      x2 = floor(a.x[2]);

      if (x2 == a.x[2]) {
        x3 = floor(a.x[3]);
      }
    }

    renormD4(&x0, &x1, &x2, &x3);
    return QD(x0, x1, x2, x3);
  }

  return QD(x0, x1, x2, x3);
}

cl_qd QDceil(cl_qd a) {
  double x0, x1, x2, x3;
  x1 = x2 = x3 = 0.0;
  x0 = ceil(a.x[0]);

  if (x0 == a.x[0]) {
    x1 = ceil(a.x[1]);
    
    if (x1 == a.x[1]) {
      x2 = ceil(a.x[2]);

      if (x2 == a.x[2]) {
        x3 = ceil(a.x[3]);
      }
    }

    renormD4(&x0, &x1, &x2, &x3);
    return QD(x0, x1, x2, x3);
  }

  return QD(x0, x1, x2, x3);
}

cl_qd QDaint(cl_qd a) {
    if( a.x[0] >= 0){
        return QDfloor(a);
    }else{
        return QDceil(a);
    }
}
#endif
