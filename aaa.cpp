//
//  aaa.cpp
//  aa
//
//  Created by j394chen on 2017-05-13.
//  Copyright © 2017 j394chen. All rights reserved.
//

#include <cstdio>
#define LL int
#define N 10005 // change this before submit

int n;
LL a[N], b[N], tmp1[2*N], tmp2[2*N], tmp3[2*N], tmp4[2*N], tmp5[2*N];

void poly_add (int size1, int size2, LL a[], LL b[], LL* result) {
    if (size1 > size2) {
        result = new LL[size1];
        for (int i  = 0; i < size2; ++i) result[i] = a[i] + b[i];
        for (int j = size2 + 1; j < size1; ++j) result[j] = a[j];
    } else {
        result = new LL[size2];
        for (int i  = 0; i < size1; ++i) result[i] = a[i] + b[i];
        for (int j = size1 + 1; j < size2; ++j) result[j] = b[j];
    }
    delete[] a;
    delete[] b;
}

void poly_minus (int size1, int size2, LL a[], LL b[], LL* result) {
    if (size1 > size2) {
        result = new LL[size1];
        for (int i  = 0; i < size2; ++i) result[i] = a[i] - b[i];
        for (int j = size2 + 1; j < size1; ++j) result[j] = - a[j];
    } else {
        result = new LL[size2];
        for (int i  = 0; i < size1; ++i) result[i] = a[i] - b[i];
        for (int j = size1 + 1; j < size2; ++j) result[j] = - b[j];
    }
    delete[] a;
    delete[] b;
}

void poly_add_help(int l, int r, int m, LL p[], LL* result) {
    LL* subarr1 = new LL[m-l];
    LL* subarr2 = new LL[r - m];
    for ( int i = l;  i < m; ++i) subarr1[i] = p[i];
    for ( int i = m + 1; i < r; ++i) subarr2[i] = p[i];
    poly_add(m-l, r-m, subarr1, subarr2, result);
}

void poly_mult(int l, int r, LL p[], LL q[], LL tmp[]) {
    if (r -l == 1) {
        tmp[0] = p[0] * q[0];
    } else {
        int m = r + l / 2;
        poly_mult(m, r, p, q, tmp1);
        poly_mult(l, m, p, q, tmp2);
        LL * x1_plus_x0, * y1_plus_y0;
        poly_add_help(l, r, m, p, x1_plus_x0);
        poly_add_help(l, r, m, q, y1_plus_y0);
        int size = r - m > m ? r - m : m;
        poly_mult(0, size, x1_plus_x0, y1_plus_y0, tmp3);
        
        LL* x0_times_y0 = new LL[2 * (m-l) - 1];
        for (int i = 0; i < 2*(m-l) -1; ++i) x0_times_y0[i] = tmp2[i];
        
        LL* x1_times_y1 = new LL(2 * (r - m) - 1);
        for (int i = 0; i < 2*(r-m) -1; ++i) x0_times_y0[i] = tmp1[i];
        
        LL* minus_result1;
        int size1 = size > 2*(m-l) - 1;
        poly_minus(size, 2*(m-l) - 1, tmp3, x0_times_y0, minus_result1);
        
        LL* minus_result2;
        poly_minus(size1, 2*(r-m)-1, minus_result1, x1_times_y1, minus_result2);
        
        
        for (int i = 0; i < m; ++i) tmp3[i] = tmp1[i];
        for (int i = m + 1; i <)
    }
}

int main() {
    scanf("%d", &n);
    for (int i = 0; i < n+1; ++i) scanf("%d", &a[i]);
    for (int i = 0; i < n+1; ++i) scanf("%d", &b[i]);
    poly_mult(0, n, a, b, tmp3);
    for (int i = 0; i < 2* n + 1; ++i) printf("%d", tmp3[i]);
    return 0;
}
