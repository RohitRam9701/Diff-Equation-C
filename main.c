#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef double (*func_yd) (double, double);

double g(double x) {
    return (-x*x*x*x + 6*x*x*x - 12*x*x + 14*x - 9) / ((1 + x)*(1 + x));
}

//f
double yd(double x, double y) {
    return y*y + g(x);
}

double* euler(double* x, double y0, int n, func_yd yd) {
    double* y = (double*)malloc(n*sizeof(double));
    y[0] = y0;
    for (int i = 1; i < n; i++)
        y[i] = y[i - 1] + (x[i] - x[i - 1])*yd(x[i - 1], y[i - 1]);
    return y;
}

double* solution(double* x, int n) {
    double* y = (double*)malloc(n*sizeof(double));
    for (int i = 0; i < n; i++)
        y[i] = (1 - x[i])*(2 - x[i])/(1 + x[i]);
    return y;
}

//a
double* one_step_a(double* x, double y0, double h, int n, func_yd yd) {
    double* y = (double*)malloc(n*sizeof(double));
    y[0] = y0;
    // h/2*(y[i+1])^2 - y[i+1] + c(i) = 0
    //      c(i) = y[i] + h/2*(f[i]+g[i])
    double a = h/2;
    double b = -1;
    double c;
    for (int i = 0; i < n - 1; i++) {
        c = y[i] + h/2*(yd(x[i], y[i]) + g(x[i + 1]));
        //printf("\n%f\n", b*b - 4*a*c);
        double y1 = -b/(2*a) + sqrt(b*b - 4*a*c)/(2*a);
        double y2 = -b/(2*a) - sqrt(b*b - 4*a*c)/(2*a);
        //printf("y1: %f y2: %f\n", y1, y2);
        if (fabs(y[i] - y1) < fabs(y[i] - y2))
            y[i + 1] = y1;
        else
            y[i + 1] = y2;
    }
    return y;
}

//b (newtone)
double* one_step_b(double* x, double y0, double h, int n, func_yd yd, int N /*iterations for eq*/) {
    double* y = (double*)malloc(n*sizeof(double));
    y[0] = y0;
    // h/2*(y[i+1])^2 - y[i+1] + c(i) = 0
    //      c(i) = y[i] + h/2*(f[i]+g[i])
    // derivative:
    //      h*y[i+1] - 1
    double a = h/2;
    double b = -1;
    double c;

    for (int i = 0; i < n - 1; i++) {
        c = y[i] + h/2*(yd(x[i], y[i]) + g(x[i + 1]));
        double _y = y[i];
        for (int j = 0; j < N; j++) {
            _y = _y - (a*_y*_y + b*_y + c)/(2*a*_y + b);
        }
        y[i + 1] = _y;
    }
    return y;
}

//v (newtone + euler)
double* one_step_v(double* x, double y0, double h, int n, func_yd yd, int N /*iterations for eq*/) {
    double* y = (double*)malloc(n*sizeof(double));
    y[0] = y0;
    double* euler_y = euler(x, y0, n, yd);
    // h/2*(y[i+1])^2 - y[i+1] + c(i) = 0
    //      c(i) = y[i] + h/2*(f[i]+g[i])
    // derivative:
    //      h*y[i+1] - 1
    double a = h/2;
    double b = -1;
    double c;

    for (int i = 0; i < n - 1; i++) {
        c = y[i] + h/2*(yd(x[i], y[i]) + g(x[i + 1]));
        double _y = euler_y[i + 1];
        for (int j = 0; j < N; j++) {
            _y = _y - (a*_y*_y + b*_y + c)/(2*a*_y + b);
        }
        y[i + 1] = _y;
    }
    free(euler_y);
    return y;
}


int main() {
    int n = 10;
    double a = 0;
    double b = 1;
    double y0 = 2;
    double h = (b - a) / n;
    double* x = (double*)malloc(n*sizeof(double));
    for (int i = 0; i < n; i++)
        x[i] = a + i*h;

    printf("\n");
    printf("exact_solution:\n");
    double* exact_solution = solution(x, n);
    for (int i = 0; i < n; i++)
        printf("%f ", exact_solution[i]);
    printf("\n\n");

    printf("euler:\n");
    double* y_euler = euler(x, y0, n, yd);
    for (int i = 0; i < n; i++)
        printf("%f ", y_euler[i]);
    printf("\n\n");

    printf("one_step(a):\n");
    double* y_one_step_a = one_step_a(x, y0, h, n, yd);
    for (int i = 0; i < n; i++)
        printf("%f ", y_one_step_a[i]);
    printf("\n\n");

    printf("one_step(b):\n");
    double* y_one_step_b = one_step_b(x, y0, h, n, yd, n);
    for (int i = 0; i < n; i++)
        printf("%f ", y_one_step_b[i]);
    printf("\n\n");

    printf("one_step(v):\n");
    double* y_one_step_v = one_step_v(x, y0, h, n, yd, n);
    for (int i = 0; i < n; i++)
        printf("%f ", y_one_step_v[i]);
    printf("\n\n");


    FILE *es_file, *eu_file, *osa_file, *osb_file, *osv_file;

    es_file = fopen("exact_solution.dat", "w");
    if (es_file != NULL) {
        for (int i = 0; i < n; i++)
            fprintf(es_file, "%f\t%f\n", x[i], exact_solution[i]);
        fclose(es_file);
    }
    eu_file = fopen("euler.dat", "w");
    if (eu_file != NULL) {
        for (int i = 0; i < n; i++)
            fprintf(eu_file, "%f\t%f\n", x[i], y_euler[i]);
        fclose(eu_file);
    }
    osa_file = fopen("one_step_a.dat", "w");
    if (osa_file != NULL) {
        for (int i = 0; i < n; i++)
            fprintf(osa_file, "%f\t%f\n", x[i], y_one_step_a[i]);
        fclose(osa_file);
    }
    osb_file = fopen("one_step_b.dat", "w");
    if (osb_file != NULL) {
        for (int i = 0; i < n; i++)
            fprintf(osb_file, "%f\t%f\n", x[i], y_one_step_b[i]);
        fclose(osb_file);
    }
    osv_file = fopen("one_step_v.dat", "w");
    if (osv_file != NULL) {
        for (int i = 0; i < n; i++)
            fprintf(osv_file, "%f\t%f\n", x[i], y_one_step_v[i]);
        fclose(osv_file);
    }

    return 0;
}