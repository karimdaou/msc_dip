#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "dip1out.txt"

#define a 10.
#define L (0.)

#define NUM 100
#define NUM2 (2*NUM)

#define C0_val 1.

#define z0 0. //на сколько поднимаем плоскость (относительно её начального положения z = -a)
#define M ((int)(z0/hz))

//#define hrho (as/(NUM-2))
#define hz (a/(NUM-1))

#define G 0.995

#if NUM == 100
#define dt 2.e-2
#define as  (3.0*10.0)
#elif NUM == 150
#define dt 2.e-1
#define as  (3.0*1.5)
#elif NUM == 200
#define dt 2.e-1
#define as  (3.0*2.0)
#elif NUM == 300
#define dt 2.e-1
#define as  (3.0*3.0)
#endif

double rho[NUM], z[NUM2], drho[NUM], dz[NUM2], psi_1[NUM][NUM2], psi_0[NUM][NUM2], vel[NUM][NUM2], V[NUM][NUM2];

double E, num, denum, maxgE, C0, E0, h, velC0;
double gE[NUM][NUM2], gC0;

void prepfuncs() {
    int i, j;

    C0 = C0_val;
    velC0 = 0;
    E0 = -0.5;

    for (i = 0; i < NUM; i++) {
        for (j = M; j < NUM2; j++) {
            if ((i == 0) && (j == NUM2 / 2)) {
                V[i][j] = 0;
                psi_1[i][j] = 0;         // дальше по psi[0][NUM2 / 2] будет происходить оптимизация C0
            } else {
                V[i][j] = -1. / (sqrt(rho[i] * rho[i] + z[j] * z[j]));
                psi_1[i][j] = sqrt(rho[i] * rho[i] + z[j] * z[j]) *
                              exp(-(rho[i] * rho[i] + z[j] * z[j]));        // psi_1 = r*exp(-r^2)
            }
            psi_0[i][j] = exp(-sqrt(rho[i] * rho[i] + z[j] * z[j]));
            vel[i][j] = 0;
        }
    }

    for (i = 0; i < NUM - 1; i++) {
        drho[i] = (rho[i + 1] - rho[i]);
    }
    for (j = M; j < NUM2 - 1; j++) {
        dz[j] = (z[j + 1] - z[j]);
    }
}


void calc_E() {
    double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16, sum17, sum18, sum19, sum20;

    int i, j;

    sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = sum11 = sum12 = sum13 = sum14 = sum15 = sum16 = sum17 = sum18 = 0;

    for (i = 0; i < NUM - 1; i++) {
        // СЛАГАЕМЫЕ от плоскости с множителем lambda
        sum11 += (2 * C0 * psi_0[i][M] * psi_1[i][M] + psi_1[i][M] * psi_1[i][M]) * rho[i] * drho[i];
        sum12 +=
                (2 * C0 * psi_0[i + 1][M] * psi_1[i + 1][M] + psi_1[i + 1][M] * psi_1[i + 1][M]) * rho[i + 1] * drho[i];

        // поверхностное слагаемое от формулы О-Г   ( dpsi_0[на границе] = (+ z / r) * psi_0 )
        sum13 += z[M] * (-V[i][M]) * psi_0[i][M] * psi_1[i][M] * rho[i] * drho[i];
        sum14 += z[M] * (-V[i + 1][M]) * psi_0[i + 1][M] * psi_1[i + 1][M] * rho[i + 1] * drho[i + 1];

        for (j = M; j < NUM2 - 1; j++) {

            // СЛАГАЕМЫЕ ВХОДЯЩИЕ В ЧИСЛИТЕЛЬ
            // (dpsi_1)^2
            sum1 += (psi_1[i + 1][j] - psi_1[i][j]) * (psi_1[i + 1][j] - psi_1[i][j]) * (rho[i + 1] + rho[i]) * dz[j] /
                    drho[i];
            sum2 += (psi_1[i + 1][j + 1] - psi_1[i][j + 1]) * (psi_1[i + 1][j + 1] - psi_1[i][j + 1]) *
                    (rho[i + 1] + rho[i]) *
                    dz[j] / drho[i];

            sum3 += (psi_1[i][j + 1] - psi_1[i][j]) * (psi_1[i][j + 1] - psi_1[i][j]) * rho[i] * drho[i] / dz[j];
            sum4 += (psi_1[i + 1][j + 1] - psi_1[i + 1][j]) * (psi_1[i + 1][j + 1] - psi_1[i + 1][j]) * rho[i + 1] *
                    drho[i] / dz[j];

            // psi_1^2 * V
            sum5 += psi_1[i][j] * psi_1[i][j] * rho[i] * V[i][j] * drho[i] * dz[j];
            sum6 += psi_1[i][j + 1] * psi_1[i][j + 1] * rho[i] * V[i][j + 1] * drho[i] * dz[j];
            sum7 += psi_1[i + 1][j] * psi_1[i + 1][j] * rho[i + 1] * V[i + 1][j] * drho[i] * dz[j];
            sum8 += psi_1[i + 1][j + 1] * psi_1[i + 1][j + 1] * rho[i + 1] * V[i + 1][j + 1] * drho[i] * dz[j];

            // psi_1 * E0 * psi_0
            sum9 += (psi_1[i][j] * psi_0[i][j] + psi_1[i][j + 1] * psi_0[i][j + 1]) * rho[i] * drho[i] * dz[j];
            sum10 += (psi_1[i + 1][j] * psi_0[i + 1][j] + psi_1[i + 1][j + 1] * psi_0[i + 1][j + 1]) * rho[i + 1] *
                     drho[i] *
                     dz[j];



            // СЛАГАЕМЫЕ ВХОДЯЩИЕ В ЗНАМЕНАТЕЛЬ
            // psi_1^2
            sum15 += (psi_1[i][j] * psi_1[i][j] + psi_1[i][j + 1] * psi_1[i][j + 1]) * rho[i] * drho[i] * dz[j];
            sum16 += (psi_1[i + 1][j] * psi_1[i + 1][j] + psi_1[i + 1][j + 1] * psi_1[i + 1][j + 1]) * rho[i + 1] *
                     drho[i] *
                     dz[j];

            // psi_0 * psi_1
            sum17 += (psi_0[i][j] * psi_1[i][j] + psi_0[i][j + 1] * psi_1[i][j + 1]) * rho[i] * drho[i] * dz[j];
            sum18 += (psi_0[i + 1][j] * psi_1[i + 1][j] + psi_0[i + 1][j + 1] * psi_1[i + 1][j + 1]) * rho[i + 1] *
                     drho[i] *
                     dz[j];

        }
    }


    sum1 = sum1 / 8.;
    sum2 = sum2 / 8.;

    sum3 = sum3 / 4.;
    sum4 = sum4 / 4.;

    sum5 = sum5 / 4.;
    sum6 = sum6 / 4.;
    sum7 = sum7 / 4.;
    sum8 = sum8 / 4.;

    sum9 = C0 * E0 * sum9 / 2.;
    sum10 = C0 * E0 * sum10 / 2.;

    sum11 = L * (sum11 / 2. + C0 * C0 / 4.) / 2.; // добавил аналитический интеграл psi_0^2 по поверхности
    sum12 = L * sum12 / 4.;

    sum13 = C0 * sum13 / 2.;
    sum14 = C0 * sum14 / 2.;

    sum15 = sum15 / 4.;
    sum16 = sum16 / 4.;

    sum17 = C0 * sum17 / 2.;
    sum18 = C0 * sum18 / 2.;

    h = (a - z0);

    //ИНТЕГРАЛЫ, посчитанные аналитически в числителе
    sum19 = C0 * C0 * (-2. + (1. - h) * exp(-2. * h)) / 8.;

    // ==//== в знаменателе
    sum20 = C0 * C0 * (2. - (1. + h) * exp(-2. * h)) / 4.;


    num = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 + sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum19;
    denum = sum15 + sum16 + sum17 + sum18 + sum20;

    E = num / denum;
}

void clear_gE() {
    int i, j;

    gC0 = 0;
    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM2; j++) {
            gE[i][j] = 0;
        }
    }
}

void calc_gE() {
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19, tmp20;
    double dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8, dv9, dv10, dv11, dv12;
    double dw1, dw2, dw3, dw4, dw5, dw6, dw7, dw8;
    int i, j;

    clear_gE();

    //ДАЛЬШЕ НАДО ПОМНИТЬ, ЧТО МЫ ЗАПИСЫВАЕМ (-1)*ГРАДИЕНТ
    //gE = (dw*num - dv*denum) / (denum*denum);

    for (i = 0; i < NUM - 1; i++) {

        gE[i][M] -= L * (C0 * psi_0[i][M] + psi_1[i][M]) * rho[i] * drho[i] / (2. * denum);
        gE[i + 1][M] -= L * (C0 * psi_0[i + 1][M] + psi_1[i + 1][M]) * rho[i + 1] * drho[i] / (2. * denum);
        gE[i][M] -= C0 * z[M] * (-V[i][M]) * psi_0[i][M] * rho[i] * drho[i] / (2. * denum);
        gE[i + 1][M] -= C0 * z[M] * (-V[i + 1][M]) * psi_0[i + 1][M] * rho[i + 1] * drho[i + 1] / (2. * denum);


        //gC0
        gC0 -= L * (psi_0[i][M] * psi_1[i][M]) * rho[i] * drho[i] / (2. * denum);
        gC0 -= L * (psi_0[i + 1][M] * psi_1[i + 1][M]) * rho[i + 1] * drho[i] / (2. * denum);

        gC0 -= z[M] * (-V[i][M]) * psi_0[i][M] * psi_1[i][M] * rho[i] * drho[i] / (2. * denum);
        gC0 -= z[M] * (-V[i + 1][M]) * psi_0[i + 1][M] * psi_1[i + 1][M] * rho[i + 1] * drho[i + 1] / (2. * denum);

        for (j = M; j < NUM2 - 1; j++) {

            //ПРОИЗВОДНАЯ ЧИСЛИТЕЛЯ

            tmp1 = (psi_1[i + 1][j] - psi_1[i][j]) * (rho[i + 1] + rho[i]) * dz[j] / (4. * drho[i]);
            // - (i,j) //+ (i+1,j)
            dv1 = -tmp1 / denum;
            gE[i][j] -= dv1;
            gE[i + 1][j] += dv1;
            ///////////////////////////////////////////////

            tmp2 = (psi_1[i + 1][j + 1] - psi_1[i][j + 1]) * (rho[i + 1] + rho[i]) * dz[j] / (4. * drho[i]);
            // - (i,j+1) //+(i+1,j+1)
            dv2 = -tmp2 / denum;
            gE[i][j + 1] -= dv2;
            gE[i + 1][j + 1] += dv2;
            ///////////////////////////////////////////////

            tmp3 = (psi_1[i][j + 1] - psi_1[i][j]) * rho[i] * drho[i] / (2. * dz[j]);
            //-(i,j) //+(i,j+1)
            dv3 = -tmp3 / denum;
            gE[i][j] -= dv3;
            gE[i][j + 1] += dv3;

            /////////////////////////////////////////////

            tmp4 = (psi_1[i + 1][j + 1] - psi_1[i + 1][j]) * rho[i + 1] * drho[i] / (2. * dz[j]);
            //-(i+1,j) //+(i+1,j+1)
            dv4 = -tmp4 / denum;
            gE[i + 1][j] -= dv4;
            gE[i + 1][j + 1] += dv4;

            tmp5 = psi_1[i][j] * rho[i] * V[i][j] * drho[i] * dz[j] / 2.;
            dv5 = -tmp5 / denum;
            gE[i][j] += dv5;

            tmp6 = psi_1[i][j + 1] * rho[i] * V[i][j + 1] * drho[i] * dz[j] / 2.;
            dv6 = -tmp6 / denum;
            gE[i][j + 1] += dv6;

            tmp7 = psi_1[i + 1][j] * rho[i + 1] * V[i + 1][j] * drho[i] * dz[j] / 2.;
            dv7 = -tmp7 / denum;
            gE[i + 1][j] += dv7;

            tmp8 = psi_1[i + 1][j + 1] * rho[i + 1] * V[i + 1][j + 1] * drho[i] * dz[j] / 2.;
            dv8 = -tmp8 / denum;
            gE[i + 1][j + 1] += dv8;

            /////////////////////////////////////////////

            tmp9 = C0 * E0 * psi_0[i][j] * rho[i] * drho[i] * dz[j] / 2.;
            dv9 = -tmp9 / denum;
            gE[i][j] += dv9;

            tmp10 = C0 * E0 * psi_0[i][j + 1] * rho[i] * drho[i] * dz[j] / 2.;
            dv10 = -tmp10 / denum;
            gE[i][j + 1] += dv10;

            tmp11 = C0 * E0 * psi_0[i + 1][j] * rho[i + 1] * drho[i] * dz[j] / 2.;
            dv11 = -tmp11 / denum;
            gE[i + 1][j] += dv11;

            tmp12 = C0 * E0 * psi_0[i + 1][j + 1] * rho[i + 1] * drho[i] * dz[j] / 2.;
            dv12 = -tmp12 / denum;
            gE[i + 1][j + 1] += dv12;


            //ПРОИЗВОДНАЯ ЗНАМЕНАТЕЛЯ
            tmp13 = psi_1[i][j] * rho[i] * drho[i] * dz[j] / 2.;
            dw1 = tmp13 * num / (denum * denum);
            gE[i][j] += dw1;

            tmp14 = psi_1[i][j + 1] * rho[i] * drho[i] * dz[j] / 2.;
            dw2 = tmp14 * num / (denum * denum);
            gE[i][j + 1] += dw2;

            tmp15 = psi_1[i + 1][j] * rho[i + 1] * drho[i] * dz[j] / 2.;
            dw3 = tmp15 * num / (denum * denum);
            gE[i + 1][j] += dw3;

            tmp16 = psi_1[i + 1][j + 1] * rho[i + 1] * drho[i] * dz[j] / 2.;
            dw4 = tmp16 * num / (denum * denum);
            gE[i + 1][j + 1] += dw4;

            tmp17 = C0 * psi_0[i][j] * rho[i] * drho[i] * dz[j] / 2.;
            dw5 = tmp17 * num / (denum * denum);
            gE[i][j] += dw5;

            tmp18 = C0 * psi_0[i][j + 1] * rho[i] * drho[i] * dz[j] / 2.;
            dw6 = tmp18 * num / (denum * denum);
            gE[i][j + 1] += dw6;

            tmp19 = C0 * psi_0[i + 1][j] * rho[i + 1] * drho[i] * dz[j] / 2.;
            dw7 = tmp19 * num / (denum * denum);
            gE[i + 1][j] += dw7;

            tmp20 = C0 * psi_0[i + 1][j + 1] * rho[i + 1] * drho[i] * dz[j] / 2.;
            dw8 = tmp20 * num / (denum * denum);
            gE[i + 1][j + 1] += dw8;


            //gC0
            gC0 -= E0 * (psi_1[i][j] * psi_0[i][j] + psi_1[i][j + 1] * psi_0[i][j + 1]) * rho[i] * drho[i] * dz[j] /
                   (2. * denum);
            gC0 -= E0 * (psi_1[i + 1][j] * psi_0[i + 1][j] + psi_1[i + 1][j + 1] * psi_0[i + 1][j + 1]) * rho[i + 1] *
                   drho[i] *
                   dz[j] / (2. * denum);

            gC0 -= (psi_0[i][j] * psi_1[i][j] + psi_0[i][j + 1] * psi_1[i][j + 1]) * rho[i] * drho[i] * dz[j] * num /
                   (2. * denum * denum);
            gC0 -= (psi_0[i + 1][j] * psi_1[i + 1][j] + psi_0[i + 1][j + 1] * psi_1[i + 1][j + 1]) * rho[i + 1] *
                   drho[i] *
                   dz[j] * num / (2. * denum * denum);

        }
    }
    // от интегралов, посчитанных аналитически в числителе
    gC0 -= C0 * (-2. + (1. - h) * exp(-2. * h)) / (4. * denum);

    gC0 -= L * C0 / (4. * denum);

    // ==//== в знаменателе
    gC0 -= C0 * (2. - (1. + h) * exp(-2. * h)) * num / (2. * denum * denum);
}

//нормировка ВФ на 1
void normto1() {
    int i, j;
    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM2; j++) {
            psi_1[i][j] = psi_1[i][j] / sqrt(denum);
        }
    }
    C0 = C0 / sqrt(denum);
    calc_E();
}

int main() {
    int i, j, k, l, step, sec;
    double F, v, n1, n2;
    double h1, hrho;

    time_t tstart = time(NULL);

    hrho = (as / (NUM - 2));
    h1 = a / ((NUM - 1) + ((NUM - 2) * (NUM - 1) * hrho) / 2.); // расстояние от нуля до первой точки

    rho[0] = 0;
    for (i = 1; i < NUM; i++) {
        rho[i] = rho[i - 1] + h1 * (1 + (i - 1) * hrho);
    }

    if (M != NUM - 1) {
        hrho = (as / (NUM - M - 1));
        h1 = ((a - z0) / ((NUM - M) + ((NUM - M - 1) * (NUM - M) * hrho) / 2.)); // расстояние от нуля до первой точки

        for (i = NUM2 / 2 + 1; i < NUM2 - M + 1; i++)
            z[i] = z[i - 1] + h1 * (1 + (i - NUM2 / 2 - 1) * hrho);

        for (i = NUM2 / 2 - 1; i > M - 1; i--)
            z[i] = z[i + 1] - h1 * (1 + (NUM2 / 2 - 1 - i) * hrho);

        for (i = NUM2 - M + 1; i < NUM2; i++)
            z[i] = z[i - 1] + z0 / (M - 1);
    } else {
        hrho = (as / (NUM - 2));
        h1 = (a / ((NUM - 1) + ((NUM - 2) * (NUM - 1) * hrho) / 2.));

        for (i = NUM2 / 2 + 1; i < NUM2; i++)
            z[i] = z[i - 1] + h1 * (1 + (i - NUM2 / 2 - 1) * hrho);
    }


    printf("M = %d	hz = %.14le	h1 = %.14le	dh = %.14le\n", M, hz, h1, hrho);
    printf("%.14le	%.14le\n", rho[0], rho[NUM - 1]);
    printf("%.14le	%.14le	%.14le	%.14le	%.14le\n\n", z[M],
           z[M + 1], z[NUM2 - M - 1], z[NUM2 - M], z[NUM2 - 1]);


////////////////////////////////////////////////////////////////////////////////////////////////////////////
    prepfuncs();

    n1 = n2 = 0;
    step = 0;

    calc_E();
    calc_gE();


    maxgE = gE[0][M];
    for (i = 0; i < NUM; i++) {
        for (j = M; j < NUM2; j++) //пров
        {
            if (maxgE < fabs(gE[i][j]))
                maxgE = fabs(gE[i][j]);
            n1 = i;
            n2 = j;
        }
    }

    //normto1();
    printf("%.14lf	%.14lf\n", E, maxgE);


    while (maxgE > 1e-10) {

        gE[0][NUM2 / 2] = 0;
        vel[0][NUM2 / 2] = 0;


        for (k = 0; k < NUM; k++) {
            for (l = M; l < NUM2; l++) {
                vel[k][l] = G * (vel[k][l] + gE[k][l] * dt);
                psi_1[k][l] = psi_1[k][l] + vel[k][l] * dt;
            }
        }
        velC0 = G * (velC0 + gC0 * dt);
        C0 = C0 + velC0 * dt;


        calc_E();
        calc_gE();
        //normto1();

        maxgE = gE[0][M];
        n1 = n2 = 0;
        for (i = 0; i < NUM; i++) {
            for (j = M; j < NUM2; j++) {
                if (maxgE < fabs(gE[i][j])) {
                    maxgE = fabs(gE[i][j]);
                    n1 = i;
                    n2 = j;
                }
            }
        }

        if (step == 5000) {
            //normto1();
            printf("%.14lf	%.14lf	%.14lf	%.14lf	%f	%f	%.14lf  \n", E, maxgE, num, denum, n1, n2, C0);
            step = 0;
        }

        step++;
    }

    printf("%.14lf	%.14lf	%.14lf	%.14le\n", E, maxgE, num, denum);
    normto1();
    printf("%.14lf	%.14lf	%.14lf	%.14le\n", E, maxgE, num, denum);

    time_t tstop = time(NULL);
    sec = tstop - tstart;
    printf("spent_time = %ld\n\a", sec);
/*	FILE *f;
	f = fopen(OUTFILE, "w");
	for (i = 0; i < NUM; i++)
		fprintf(f, "%.14le	%.14le\n", rho[i], z[i]);
*/
    getchar();
    return 0;
}