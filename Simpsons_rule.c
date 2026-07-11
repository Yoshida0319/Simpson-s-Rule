#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>

#define N 16 //n次式近似
#define NN pow(10,6) //分割数
#define NNN 332 //mpfrのbit数(332だと100桁精度)
#define NNNN 100 //表示桁数
#define start 0 //積分範囲_下端
#define end 1 //積分範囲_上端

void function(mpfr_t function_value, mpfr_t x){//関数を入力//f(x)=sin(10x)
    mpfr_mul_si(x, x, 10, MPFR_RNDN);
    mpfr_sin(function_value, x, MPFR_RNDN);
}

void factorial(mpz_t big, int n){
    if(n != 0){
        factorial(big, n-1);
        mpz_mul_ui(big, big, n);
    }
}

int main(void){
    if(N == 0){
        return 0;
    }
    int i;
    mpfr_set_default_prec(NNN);//332
    mpfr_t xx, function_value, hs, integ, start_, end_, a[N+1];
    mpfr_inits(xx, function_value, hs, integ, start_, end_, NULL);
    mpfr_set_d(start_, start, MPFR_RNDN);
    mpfr_set_d(end_, end, MPFR_RNDN);
    mpfr_sub(hs, end_, start_, MPFR_RNDN);
    mpfr_clears(start_, end_, NULL);
    mpz_t s[N+1];
    for(i=0; i<=N; i++){
        mpfr_init(a[i]);
        mpz_init(s[i]);
    }
    mpz_set_ui(s[0], 1);
    int no, zi, mod, Nuu;
    int j, k, l, m;
    for(no=0; no<=N; no++){
        if(no != 0){
            for(m=N+1; m>0; m--){
                if(m == 1){
                    mpz_mul_ui(s[m-1], s[m-1], no);
                }else{
                    mpz_mul_ui(s[m-1], s[m-1], no);
                    mpz_add(s[m-1], s[m-1], s[m-2]);
                }
            }
        }
        for(i=0; i<=N; i++){
            mpfr_set_d(a[i], 0, MPFR_RNDN);
        }
        for (j=0; j<=no; j++) {
            for(k=1; k<=no+1; k++){
                for(l=k; l<=no+1; l++){
                    mpz_t big, ss;
                    mpz_inits(big, ss, NULL);
                    mpz_ui_pow_ui(ss, j, l-k);
                    mpz_ui_pow_ui(big, no, k-1);
                    mpz_mul(ss, ss, big);
                    mpz_mul(ss, ss, s[l-1]);
                    mpz_mul_si(ss, ss, (j+l+1 & 1) ? -1 : 1);
                    mpfr_t tmp;
                    mpfr_init(tmp);
                    mpfr_set_z(tmp, ss, MPFR_RNDN);
                    mpfr_div_si(tmp, tmp, k, MPFR_RNDN);
                    mpfr_add(a[j], a[j], tmp, MPFR_RNDN);
                    mpz_clears(big, ss, NULL);
                    mpfr_clear(tmp);
                }
            }
            mpz_t big;
            mpz_init(big);
            mpz_set_ui(big, 1);
            factorial(big, no-j);
            factorial(big, j);
            mpfr_t tmp;
            mpfr_init(tmp);
            mpfr_set_z(tmp, big, MPFR_RNDN);
            mpfr_div(a[j], a[j], tmp, MPFR_RNDN);
            mpz_clear(big);
            mpfr_clear(tmp);
        }
    }

    no = N;
    int Nu = NN/no;
    mpfr_t h, sum;
    mpfr_inits(h, sum, NULL);
    mpfr_div_si(h, hs, no*Nu, MPFR_RNDN);
    mpfr_set_d(integ, 0, MPFR_RNDN);
    for (mod=0; mod<=no; mod++){
        mpfr_set_d(sum, 0, MPFR_RNDN);
        for (Nuu=0; Nuu<Nu; Nuu++){
            mpfr_set_d(xx, start, MPFR_RNDN);
            mpfr_div(xx, xx, h, MPFR_RNDN);
            mpfr_add_si(xx, xx, no*Nuu + mod, MPFR_RNDN);
            mpfr_mul(xx, xx, h, MPFR_RNDN);
            function(function_value, xx);
            mpfr_add(sum, sum, function_value, MPFR_RNDN);
        }
        mpfr_mul(sum, sum, a[mod], MPFR_RNDN);
        mpfr_add(integ, integ, sum, MPFR_RNDN);
    }
    mpfr_mul_si(h, h, no, MPFR_RNDN);
    mpfr_mul(integ, integ, h, MPFR_RNDN);
    mpfr_clears(h, sum, NULL);
    printf("%d ", no*Nu);
    mpfr_printf("%.*Rf\n", NNNN, integ);
    mpfr_clears(xx, function_value, hs, integ, NULL);
    for(i=0; i<=N; i++){
        mpfr_clear(a[i]);
        mpz_clear(s[i]);
    }
    return 0;
}
