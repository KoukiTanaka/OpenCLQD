#include "dd_inline.h"
#include "dd_const.h"
#ifndef CLDD_H
#define CLDD_H

cl_dd DDPowI(cl_dd a, int n);
cl_dd DDnpwr(cl_dd a,int n);
cl_dd DDsqrt(cl_dd a);
cl_dd DDnroot(cl_dd a, int n);
double DDinv_fact(int i,int j);
cl_dd DDexp(cl_dd a);
cl_dd DDlog(cl_dd a);
cl_dd DDpow(cl_dd a, cl_dd b) ;
cl_dd DDlog10(cl_dd a);
double DDsin_table(int i, int j);
double DDcos_table(int i, int j);
cl_dd DDsin_taylor(cl_dd a);
cl_dd DDcos_taylor(cl_dd a);
void DDsincos_taylor(cl_dd a, cl_dd *sin_a, cl_dd *cos_a);
cl_dd DDsin(cl_dd a);
cl_dd DDcos(cl_dd a);
void DDsincos(cl_dd a,cl_dd *sin_a,cl_dd *cos_a);
cl_dd DDtan(cl_dd a);
cl_dd DDatan(cl_dd a);
cl_dd DDatan2(cl_dd y, cl_dd x);
cl_dd DDasin(cl_dd a);
cl_dd DDacos(cl_dd a);
cl_dd DDsinh(cl_dd a);
cl_dd DDcosh(cl_dd a);
cl_dd DDtanh(cl_dd a);
void DDsincosh(cl_dd a, cl_dd *s, cl_dd *c);
cl_dd DDasinh(cl_dd a);
cl_dd DDacosh(cl_dd a);
cl_dd DDatanh(cl_dd a);
cl_dd DDfmod(cl_dd a, cl_dd b);
cl_dd DDrand(int *seed);
#ifdef QD_MALLOC
cl_dd DDpolyeval( cl_dd *c, int n, cl_dd x);
cl_dd DDpolyroot(cl_dd *c, int n, cl_dd x0, int max_iter, double thresh);
#else
cl_dd DDpolyeval( __global cl_dd *c, int n, cl_dd x);
cl_dd DDpolyevalDifferential(__global cl_dd *c, int n, cl_dd x);
cl_dd DDpolyroot(__global cl_dd *c, int n, cl_dd x0, int max_iter, double thresh);
#endif

cl_dd DDPowI(cl_dd a, int n){
    int i;
    cl_dd x;
    x = DtoDD(1.0);
    for(i=0;i<n;i++){
        x = DDmul(x,a);
    }
    
    return x;
    
}

cl_dd DDnpwr(cl_dd a,int n){
    cl_dd r = a;
    cl_dd s = DtoDD(1.0);
    int N = abs(n);
    
    if( N > 1){
        while(N > 0){
            if( N % 2 == 1){
                s = DDmul(s,r);
            }
            N /= 2;
            if( N > 0){
                r = DDsqr(r);
            }
        }
    }else{
        s = r;
    }
    if( n < 0){
        return DDdiv(DtoDD(1.0),s);
    }
    return s;
}

cl_dd DDsqrt(cl_dd a){
    double x = 1.0 / sqrt(a.x[0]);
    double ax = a.x[0] * x;
    return DDadd(DtoDD(DDsub(a,DDsqrD(ax)).x[0] * (x * 0.5)),DtoDD(ax));
}

cl_dd DDnroot(cl_dd a, int n){
    cl_dd r = DDabs(a);
    cl_dd x = DtoDD(exp(-log(r.x[0])/n));
    
    x = DDadd(x, DDdiv( DDmul(x , DDsub( DtoDD(1.0), DDmul( r, DDnpwr(x , n ) ) ) ) , DtoDD(n)));
    
    //x = DDabs(x);
    
    return DDdiv(DtoDD(1.0),x);
}

/*
OpenCL環境でグローバル変数を利用すると致命的なバグが発生したため、
関数を呼び出す形式に変更している。
*/
double DDinv_fact(int i,int j){
    double inv_fact[15][2] = {
        { 1.66666666666666657e-01,  9.25185853854297066e-18},
        { 4.16666666666666644e-02,  2.31296463463574266e-18},
        { 8.33333333333333322e-03,  1.15648231731787138e-19},
        { 1.38888888888888894e-03, -5.30054395437357706e-20},
        { 1.98412698412698413e-04,  1.72095582934207053e-22},
        { 2.48015873015873016e-05,  2.15119478667758816e-23},
        { 2.75573192239858925e-06, -1.85839327404647208e-22},
        { 2.75573192239858883e-07,  2.37677146222502973e-23},
        { 2.50521083854417202e-08, -1.44881407093591197e-24},
        { 2.08767569878681002e-09, -1.20734505911325997e-25},
        { 1.60590438368216133e-10,  1.25852945887520981e-26},
        { 1.14707455977297245e-11,  2.06555127528307454e-28},
        { 7.64716373181981641e-13,  7.03872877733453001e-30},
        { 4.77947733238738525e-14,  4.39920548583408126e-31},
        { 2.81145725434552060e-15,  1.65088427308614326e-31}
    };
    
    return inv_fact[i][j];
}

cl_dd DDexp(cl_dd a){
    const double k = 512.0;
    const double inv_k = 1.0/k;
    
    
    if (a.x[0] <= -709.0)
        return DtoDD(0.0);

    if (a.x[0] >=  709.0)
        return DD_inf;

    if (DDisZero(a))
        return DtoDD(1.0);

    if (DDisOne(a))
        return DD_e;
    
    
    double m = floor(a.x[0] / DD_log2.x[0] + 0.5);
    cl_dd r = DDmul_pwr2(DDsub(a,DDmulD(DD_log2,m)),inv_k);
    cl_dd s,t,p;
    
    p = DDsqr(r);
    s = DDadd(r,DDmul_pwr2(p,0.5));
    p = DDmul(p,r);
    t = DDmul(p,DD(DDinv_fact(0,0),DDinv_fact(0,1)));
    
    int i=0;
    do{
        s = DDadd(s,t);
        p = DDmul(p,r);
        ++i;
        t = DDmul(p,DD(DDinv_fact(i,0),DDinv_fact(i,1)));
    }while(fabs(t.x[0]) > inv_k * DD_eps && i < 5);
    s = DDadd(s,t);
    
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    s = DDadd(DDmul_pwr2(s,2.0),DDsqr(s));
    
    s = DDaddD(s,1.0);
    
    return DDldexp(s,m);
}

cl_dd DDlog(cl_dd a){
    cl_dd x = DtoDD(log(a.x[0]));
    x = DDsubD(DDadd(x,DDmul(a,DDexp(DDminus(x)))),1.0);
    return x;
}

cl_dd DDpow(cl_dd a, cl_dd b) {
    return DDexp(DDmul(b,DDlog(a)));
}

cl_dd DDlog10(cl_dd a){
    return DDdiv(DDlog(a),DD_log10);
}

/*
OpenCL環境でグローバル変数を利用すると致命的なバグが発生したため、
関数を呼び出す形式に変更している。
*/
double DDsin_table(int i, int j){
    double sin_table [4][2] = {
        {1.950903220161282758e-01, -7.991079068461731263e-18},
        {3.826834323650897818e-01, -1.005077269646158761e-17},
        {5.555702330196021776e-01,  4.709410940561676821e-17},
        {7.071067811865475727e-01, -4.833646656726456726e-17}
    };
    return sin_table[i][j];
}
double DDcos_table(int i, int j){
    double cos_table [4][2] = {
        {9.807852804032304306e-01, 1.854693999782500573e-17},
        {9.238795325112867385e-01, 1.764504708433667706e-17},
        {8.314696123025452357e-01, 1.407385698472802389e-18},
        {7.071067811865475727e-01, -4.833646656726456726e-17}
    };
    return cos_table[i][j];
}

cl_dd DDsin_taylor(cl_dd a){
    double thresh = 0.5 * fabs(a.x[0]) * DD_eps;
    
    cl_dd r,s,t,x;
    
    if( DDisZero(a) ){
        return DtoDD(0.0);
    }
    int i = 0;
    x = DDminus(DDsqr(a));
    s = a;
    r = a;
    
    
    do{
        r = DDmul(r,x);
        t = DDmul(r,DD(DDinv_fact(i,0),DDinv_fact(i,1)));
        s = DDadd(s,t);
        i+=2;
    }while(i < 15 && fabs(t.x[0]) > thresh);
    
    return s;
}

cl_dd DDcos_taylor(cl_dd a){
    const double thresh = 0.5 * DD_eps;
    cl_dd r,s,t,x;
    
    if( DDisZero(a)){
        return DtoDD(1.0);
    }
    
    x = DDminus(DDsqr(a));
    r = x;
    s = DDadd(DtoDD(1.0),DDmul_pwr2(r,0.5));
    int i=1;
    do{
        r = DDmul(r,x);
        t = DDmul(r,DD(DDinv_fact(i,0),DDinv_fact(i,1)));
        s = DDadd(s,t);
        i+=2;
    }while(i < 15 && fabs(t.x[0]) > thresh);
    
    return s;
}
void DDsincos_taylor(cl_dd a, cl_dd *sin_a, cl_dd *cos_a){
    if( DDisZero(a)){
        *sin_a = DtoDD(0.0);
        *cos_a = DtoDD(1.0);
        return;
    }
    double q;
    *sin_a = DDsin_taylor(a);
    *cos_a = DDsqrt(DDsub(DtoDD(1.0),DDsqr(*sin_a)));
}

cl_dd DDsin(cl_dd a){
    if( DDisZero(a)){
       return DtoDD(0.0);
    }
    cl_dd z = DtoDD(nint(DDtoD(DDdiv(a,DD_2pi))));
    cl_dd r = DDsub(a,DDmul(DD_2pi,z));
    
    
    cl_dd t;
    double q = floor(r.x[0] / DD_pi2.x[0]+0.5);
    t = DDsub(r,DDmulD(DD_pi2,q));
    int j = q;
    q = floor(t.x[0] / DD_pi16.x[0] + 0.5);
    t = DDsub(t,DDmulD(DD_pi16,q));
    int k = q;
    int abs_k = abs(k);
    
    if( k == 0 ){
        switch(j){
            case 0:
                return DDsin_taylor(t);
            case 1:
                return DDcos_taylor(t);
            case -1:
                return DDminus(DDcos_taylor(t));
            default:
                return DDminus(DDsin_taylor(t));
        }
    }
    cl_dd u,v;
    u = DD(DDcos_table(abs_k-1,0),DDcos_table(abs_k-1,1));
    v = DD(DDsin_table(abs_k-1,0),DDsin_table(abs_k-1,1));
    
    cl_dd sin_t,cos_t;
    DDsincos_taylor(t,&sin_t,&cos_t);
    
    if( j == 0){
        if( k > 0 ){
            r = DDadd(DDmul(u,sin_t),DDmul(v,cos_t));
        }else{
            r = DDsub(DDmul(u,sin_t),DDmul(v,cos_t));
        }
    }else if( j == 1){
        if( k > 0 ){
            r = DDsub(DDmul(u,cos_t),DDmul(v,sin_t));
        }else{
            r = DDadd(DDmul(u,cos_t),DDmul(v,sin_t));
        }
    }else if(j == - 1){
        if( k > 0 ){
            r = DDsub(DDmul(v,sin_t),DDmul(u,cos_t));
        }else{
            r = DDsub(DDminus(DDmul(u,cos_t)),DDmul(v,sin_t));
        }
    }else{
        if( k > 0 ){
            r = DDsub(DDminus(DDmul(u,sin_t)),DDmul(v,cos_t));
        }else{
            r = DDsub(DDmul(v,cos_t),DDmul(u,sin_t));
        }
    }
    return r;
}


cl_dd DDcos(cl_dd a){
    if( DDisZero(a)){
        return DtoDD(1.0);
    }
    
    cl_dd z = DtoDD(nint(DDtoD(DDdiv(a,DD_2pi))));
    cl_dd r = DDsub(a,DDmul(z,DD_2pi));
    
    cl_dd t;
    double q = floor(r.x[0] / DD_pi2.x[0] + 0.5);
    t = DDsub(r,DDmulD(DD_pi2,q));
    int j = q;
    q = floor(t.x[0]/DD_pi16.x[0] + 0.5);
    t = DDsub(t,DDmulD(DD_pi16,q));
    int k = q;
    int abs_k = abs(k);
    if (k == 0) {
        switch (j) {
            case 0:
                return DDcos_taylor(t);
            case 1:
                return DDminus(DDsin_taylor(t));
            case -1:
                return DDsin_taylor(t);
            default:
                return DDminus(DDcos_taylor(t));
        }
    }
    
    cl_dd sin_t,cos_t;
    DDsincos_taylor(t,&sin_t,&cos_t);
    cl_dd u,v;
    u = DD(DDcos_table(abs_k-1,0),DDcos_table(abs_k-1,1));
    v = DD(DDsin_table(abs_k-1,0),DDsin_table(abs_k-1,1));
    if( j== 0){
        if( k > 0 ){
            r = DDsub(DDmul(u,cos_t),DDmul(v,sin_t));
        }else{
            r = DDadd(DDmul(u,cos_t),DDmul(v,sin_t));
        }
    }else if( j == 1){
        if( k > 0 ){
            r = DDsub(DDminus(DDmul(u,sin_t)),DDmul(v,cos_t));
        }else{
            r = DDsub(DDmul(v,cos_t),DDmul(u,sin_t));
        }
    }else if(j == -1){
        if( k > 0 ){
            r = DDadd(DDmul(u,sin_t),DDmul(v,cos_t));
        }else{
            r = DDsub(DDmul(u,sin_t),DDmul(v,cos_t));
        }
    }else{
        if( k > 0 ){
            r = DDsub(DDmul(v,sin_t),DDmul(u,cos_t));
        }else{
            r = DDsub(DDminus(DDmul(u,cos_t)),DDmul(v,sin_t));
        }
    }
    
    return r;
}


void DDsincos(cl_dd a,cl_dd *sin_a,cl_dd *cos_a){
    if( DDisZero(a)){
        *sin_a = DtoDD(0.0);
        *cos_a = DtoDD(1.0);
    }
    cl_dd z = DtoDD(nint(DDtoD(DDdiv(a,DD_2pi))));
    cl_dd r = DDsub(a,DDmul(DD_2pi,z));
    
    cl_dd t;
    double q = floor(r.x[0]/DD_pi2.x[0] + 0.5);
    t = DDsub(r,DDmulD(DD_pi2,q));
    int j = q;
    int abs_j = abs(j);
    q = floor(t.x[0] / DD_pi16.x[0] + 0.5);
    t = DDsub(t,DDmulD(DD_pi16,q));
    int k = q;
    int abs_k = abs(k);
    
    cl_dd sin_t,cos_t;
    cl_dd s,c;
    
    DDsincos_taylor(t,&sin_t,&cos_t);
    if( abs_k == 0 ){
        s = sin_t;
        c = cos_t;
    }else{
        cl_dd u,v;
        u = DD(DDcos_table(abs_k-1,0),DDcos_table(abs_k-1,1));
        v = DD(DDsin_table(abs_k-1,0),DDsin_table(abs_k-1,1));
        
        if( k > 0 ){
            s = DDadd(DDmul(u,sin_t),DDmul(v,cos_t));
            c = DDsub(DDmul(u,cos_t),DDmul(v,sin_t));
        }else{
            s = DDsub(DDmul(u,sin_t),DDmul(v,cos_t));
            c = DDadd(DDmul(u,cos_t),DDmul(v,sin_t));
        }
    }
    
    if(abs_j == 0){
        *sin_a = s;
        *cos_a = c;
    }else if(j == 1){
        *sin_a = c;
        *cos_a = DDminus(s);
    }else if(j == -1){
        *sin_a = DDminus(c);
        *cos_a = s;
    }else{
        *sin_a = DDminus(s);
        *cos_a = DDminus(c);
    }
}

cl_dd DDtan(cl_dd a){
    cl_dd s,c;
    DDsincos(a,&s,&c);
    return DDdiv(s,c);
}

cl_dd DDatan(cl_dd a){
    return DDatan2(a,DtoDD(1.0));
}


cl_dd DDatan2(cl_dd y, cl_dd x){
    
    if(DDisZero(x)){
        if(DDisZero(y)){
            return DD_nan;
        }
        if( y.x[0] > 0 ){
            return DD_pi2;
        }else{
            return DDminus(DD_pi2);
        }
    }else if(DDisZero(y)){
        if( x.x[0] > 0 ){
            return DtoDD(0.0);
        }else{
            return DD_pi;
        }
    }

    if( DDEqual(x,y) == true){
        if( y.x[0] > 0 ){
            return DD_pi4;
        }else{
            return DDminus(DD_3pi4);
        }
    }
    if( DDEqual(x,DDminus(y)) == true){
        if(y.x[0] > 0){
            return DD_3pi4;
        }else{
            return DDminus(DD_pi4);
        }
        
    }
    cl_dd r = DDsqrt(DDadd(DDsqr(x),DDsqr(y)));
    cl_dd xx = DDdiv(x,r);
    cl_dd yy = DDdiv(y,r);
      
    cl_dd z = DtoDD(atan2(DDtoD(y),DDtoD(x)));
    cl_dd sin_z,cos_z;

    if( fabs(xx.x[0]) > fabs(yy.x[0])){
        DDsincos(z,&sin_z,&cos_z);
        z = DDadd(z,DDdiv(DDsub(yy,sin_z),cos_z));
    }else{
        DDsincos(z,&sin_z,&cos_z);
        z = DDsub(z,DDdiv(DDsub(xx,cos_z),sin_z));
    }
    return z;
}

cl_dd DDasin(cl_dd a){
    cl_dd abs_a = DDabs(a);

    
    if( abs_a.x[0] > 1.0){
        return DD_nan;
    }

    if( DDisOne(abs_a) == true){
        if( a.x[0] > 0){
            return DD_pi2;
        }else{
            return DDminus(DD_pi);
        }
    }
    return DDatan2(a,DDsqrt(DsubDD(1.0,DDsqr(a))));
}


cl_dd DDacos(cl_dd a){
    cl_dd abs_a = DDabs(a);
    
    if( abs_a.x[0] > 1.0){
        return DD_nan;
    }
    if( DDisOne(abs_a) == true){
        if( a.x[0] > 0){
            return DtoDD(0.0f);
        }else{
            return DD_pi;
        }
    }
    return DDatan2(DDsqrt(DsubDD(1.0,DDsqr(a))),a);
}

cl_dd DDsinh(cl_dd a){
    if( DDisZero(a)){
        return DtoDD(0.0);
    }
    if(fabs(a.x[0]) > 0.05){
        cl_dd ea = DDexp(a);
        return DDmul_pwr2(DDsub(ea,DDinv(ea)),0.5);
    }

    cl_dd s = a;
    cl_dd t = a;
    cl_dd r = DDsqr(t);
    double m = 1.0;
    double thresh = fabs(DDtoD(a) * DD_eps);
    
    do {
        m += 2.0;
        t = DDmul(t,r);
        t = DDdivD(t, (m-1) * m);
        s = DDadd(s,t);
    } while (fabs(t.x[0]) > thresh);

    return s;
}

cl_dd DDcosh(cl_dd a){
    if( DDisZero(a)){
        return DtoDD(1.0);
    }
    cl_dd ea = DDexp(a);
    return DDmul_pwr2(DDadd(ea,DDinv(ea)),0.5);
}


cl_dd DDtanh(cl_dd a){
    if( DDisZero(a)){
        return DtoDD(0.0);
    }
    if( fabs(DDtoD(a)) > 0.05 ){
        cl_dd ea = DDexp(a);
        cl_dd inv_ea = DDinv(ea);
        return DDdiv(DDsub(ea,inv_ea),DDadd(ea,inv_ea));
    }else{
        cl_dd s,c;
        s = DDsinh(a);
        c = DDsqrt(DaddDD(1.0,DDsqr(s)));
        return DDdiv(s,c);
    }
}

void DDsincosh(cl_dd a, cl_dd *s, cl_dd *c){
    if( fabs(DDtoD(a)) <= 0.05){
        *s = DDsinh(a);
        *c = DDsqrt(DaddDD(1.0,DDsqr(*s)));
    }else{
        cl_dd ea = DDexp(a);
        cl_dd inv_ea = DDinv(ea);
        *s = DDmul_pwr2(DDsub(ea,inv_ea),0.5);
        *c = DDmul_pwr2(DDadd(ea,inv_ea),0.5);
    }
}

void DDsincoshD(double a, cl_dd *s, cl_dd *c){
    DDsincosh(DtoDD(a),s,c);
}

cl_dd DDasinh(cl_dd a){
    return DDlog(DDadd(a,DDsqrt(DDaddD(DDsqr(a),1.0))));
}

cl_dd DDacosh(cl_dd a){
    if( a.x[0] < 1.0){
        return DD_nan;
    }
    return DDlog(DDadd(a,DDsqrt(DDsubD(DDsqr(a),1.0))));
}

cl_dd DDatanh(cl_dd a){
    if( a.x[0] >= 1.0){
        return DD_nan;
    }
    return DDmul_pwr2(DDlog(DDdiv(DaddDD(1.0,a),DsubDD(1.0,a))),0.5);
}
cl_dd DDfmod(cl_dd a, cl_dd b){
    cl_dd n = DDaint(DDdiv(a,b));
    return DDsub(a,DDmul(b,n));
}
//乱数生成
cl_dd DDrand(int *seed){
    double m_const = 4.6566128730773926e-10;  // = 2^{-31}
    double m = m_const;
    cl_dd r = DtoDD(0.0);
    double d = 1.0;
    int i;
    for(i=0;i<4;i++,m*=m_const){
        d = OpenCLRand(seed) * m;
        r = DDaddD(r,d);
    }

    return r;
}

#ifdef QD_MALLOC 
//malloc　freeができるデバイスの場合
//nは配列サイズじゃなくてn次式のnで c[0]が定数項
cl_dd DDpolyeval( cl_dd *c, int n, cl_dd x){
    cl_dd r = c[n];
    int i;
    
    for(i = n-1;i >= 0; i--){
        r = DDmul( r,x);
        r = DDadd(r,c[i]);
    }
    return r;
}
cl_dd DDpolyroot(cl_dd *c, int n, cl_dd x0, int max_iter, double thresh){
    cl_dd x = x0;
    cl_dd f;
    cl_dd *d;
    d = (cl_dd*)malloc(sizeof(cl_dd)*n);
    bool conv = false;
    int i;
    double max_c = fabs(DDtoD(c[0]));
    double v;

    if( thresh == 0.0){
        thresh = DD_eps;
    }
    for(i=1;i<=n;i++){
        v = fabs(DDtoD(c[i]));
        if( v > max_c){
            max_c = v;
        }
        d[i-1] = DDmulD(c[i],i);
    }
    thresh *= max_c;

    for(i=0;i<max_iter; i++){
        f = DDpolyeval(c,n,x);
        if( fabs(DDtoD(f)) < thresh){
            conv = true;
            break;
        }
        x = DDsub(x, DDdiv(f,DDpolyeval(d,n-1,x)));
    }
    free(d);
    if(!conv){
        return DD_nan;
    }
    return x;
}
#else
//malloc free が利用できないデバイスで実行する場合
cl_dd DDpolyeval( __global cl_dd *c, int n, cl_dd x){
    //cl_dd r = c[n];
    cl_dd r = DtoDD(0.0);
    int i;
    
    for(i = n;i >= 0; i--){
        r = DDmul(r,x);
        r = DDadd(r,c[i]);
    }
    return r;
}
//微分した値を計算 nはf(x)の次数
cl_dd DDpolyevalDifferential(__global cl_dd *c, int n, cl_dd x){
    cl_dd r = DtoDD(0.0);
    int i;
    
    for(i = n-1;i >= 0; i--){
        r = DDmul(r,x);
        r = DDadd(r,DDmulD(c[i+1],i+1));
    }
    return r;
}
cl_dd DDpolyroot(__global cl_dd *c, int n, cl_dd x0, int max_iter, double thresh){
    cl_dd x = x0;
    cl_dd f;
    cl_dd temp;
    bool conv = false;
    int i,j;
    double max_c = fabs(DDtoD(c[0]));
    double v;

    if( thresh == 0.0){
        thresh = DD_eps;
    }
    
    for(i=1;i<=n;i++){
        v = fabs(DDtoD(c[i]));
        if( v > max_c){
            max_c = v;
        }
    }
    thresh *= max_c;

    for(i=0;i<max_iter; i++){
        //f(x)を計算
        f = DDpolyeval(c,n,x);
        temp = DDabs(f);
        if( temp.x[0] < thresh ){
            conv = true;
            break;
        }
        //x = x - f(x) / f'(x)を計算
        x = DDsub(x, DDdiv(f,DDpolyevalDifferential(c,n,x)));
    }
    if(!conv){
        return DD_nan;
    }
    
    return x;
}
#endif

#endif
