#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "output.txt"

#define H 10.   //характерный размер квантовой системы
//(пространственная бесконечность)
#define a 3.    //расстояние, на которое развели ядра
#define L (0.)  // lambda parameter

#define NUM 100
#define NUM2 (2*NUM)

#define step1 (H / (NUM + ((NUM-2) * (NUM - 1) * step_x) / 2.))  // расстояние от нуля до первой точки
//(использована формула суммы арифметической прогрессии) (сумма слагаемых дает H)
#define c 0.699 //коэффициент определяющий размер до ближайшей к атому точку вдоль оси z

#define z0 5.   ///на сколько поднимаем плоскость (относительно её начального положения z = -a)
//физическое расстояние между атомами и плоскостью = (a-z0)
#define M ((int)(z0 / step_z))    //индекс, начиная с которого заполняется массив точек z

#define H_var   4.  //чем больше параметр(переменного шага),
// тем меньше будут первые шаги вдоль соответстующей оси
// (график будет сильнее прогибаться)
#define step_x (H_var/(NUM - 2))
#define step_z (H/(NUM-1+c))    ////////CHECK DENUMENATOR m.b. NUM-2?????\

#define G 0.995
#define dt 2.e-1

//##### GLOBAL VARIABLES    #####
double x[NUM], y[NUM],
/*NUM, т.к. задача симметрична относительно начала координат относительно x и y*/
        z[NUM2], dx[NUM], dy[NUM], dz[NUM2], psi[NUM][NUM][NUM2], vel[NUM][NUM][NUM2], V[NUM][NUM][NUM2];

double E, numerator, denumerator, max_grad_E, E_kinetic, E_potential, surface_term;
double grad_E[NUM][NUM][NUM2];

void prepare_functions() {
    int i, j, k;
    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM; j++) {
            for (k = M; k < NUM2; k++) {
                V[i][j][k] = 1. / (sqrt((x[i] + a / 2.) * (x[i] + a / 2.) + y[j] * y[j] + z[k] * z[k])) + \
                                1. / (sqrt((x[i] - a / 2.) * (x[i] - a / 2.) + y[j] * y[j] + z[k] * z[k]));

                psi[i][j][k] = exp(-sqrt(x[i] * x[i] + y[j] * y[j] + z[k] * z[k]));
                vel[i][j][k] = 0;
            }

        }

    }

    for (i = 0; i < NUM - 1; i++) {
        dx[i] = x[i + 1] - x[i];
        dy[i] = y[i + 1] - y[i];
    }

    for (j = M; j < NUM2 - 1; j++) {
        dz[j] = z[j + 1] - z[j];
    }
}

void calculate_E() {
    double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10,
            sum11, sum12, sum13, sum14, sum15, sum16, sum17, sum18, sum19, sum20,
            sum21, sum22, sum23, sum24, sum25, sum26, sum27, sum28, sum29, sum30,
            sum31, sum32;

    int i, j, k;

    sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = \
    sum11 = sum12 = sum13 = sum14 = sum15 = sum16 = sum17 = sum18 = sum19 = sum20 = \
    sum21 = sum22 = sum23 = sum24 = sum25 = sum26 = sum27 = sum28 = sum29 = sum30 = \
    sum31 = sum32 = 0;

    for (i = 0; i < NUM - 1; i++) {
        for (j = 0; j < NUM - 1; j++) {
            //поверхностные члены (входят в числитель функционала энергии)
            // (из слагаемого с множителем L, отвечающего за ГУ)
            sum21 += psi[i][j][M] * psi[i][j][M] * dx[i] * dy[j];
            sum22 += psi[i + 1][j][M] * psi[i + 1][j][M] * dx[i] * dy[j];
            sum23 += psi[i][j + 1][M] * psi[i][j + 1][M] * dx[i] * dy[j];
            sum24 += psi[i + 1][j + 1][M] * psi[i + 1][j + 1][M] * dx[i] * dy[j];

            for (k = M; k < NUM2 - 1; k++) {

                //слагаемые, входящие в числитель функционала энергии
                sum1 += (psi[i + 1][j][k] - psi[i][j][k]) * (psi[i + 1][j][k] - psi[i][j][k]) * dy[j] * dz[k] / dx[i];

                sum2 += (psi[i + 1][j + 1][k] - psi[i][j + 1][k]) * \
                            (psi[i + 1][j + 1][k] - psi[i][j + 1][k]) * dy[j] * dz[k] / dx[i];

                sum3 += (psi[i + 1][j][k + 1] - psi[i][j][k + 1]) * \
                            (psi[i + 1][j][k + 1] - psi[i][j][k + 1]) * dy[j] * dz[k] / dx[i];

                sum4 += (psi[i + 1][j + 1][k + 1] - psi[i][j + 1][k + 1]) * \
                            (psi[i + 1][j + 1][k + 1] - psi[i][j + 1][k + 1]) * dy[j] * dz[k] / dx[i];

                sum5 += (psi[i][j + 1][k] - psi[i][j][k]) * (psi[i][j + 1][k] - psi[i][j][k]) * dx[i] * dz[k] / dy[j];
                sum6 += (psi[i + 1][j + 1][k] - psi[i + 1][j][k]) * \
                            (psi[i + 1][j + 1][k] - psi[i + 1][j][k]) * dx[i] * dz[k] / dy[j];

                sum7 += (psi[i][j + 1][k + 1] - psi[i][j][k + 1]) * \
                            (psi[i][j + 1][k + 1] - psi[i][j][k + 1]) * dx[i] * dz[k] / dy[j];

                sum8 += (psi[i + 1][j + 1][k + 1] - psi[i + 1][j][k + 1]) * \
                            (psi[i + 1][j + 1][k + 1] - psi[i + 1][j][k + 1]) * dx[i] * dz[k] / dy[j];

                sum9 += (psi[i][j][k + 1] - psi[i][j][k]) * (psi[i][j][k + 1] - psi[i][j][k]) * dx[i] * dy[j] / dz[k];
                sum10 += (psi[i + 1][j][k + 1] - psi[i + 1][j][k]) * \
                            (psi[i + 1][j][k + 1] - psi[i + 1][j][k]) * dx[i] * dy[j] / dz[k];

                sum11 += (psi[i][j + 1][k + 1] - psi[i][j + 1][k]) * \
                            (psi[i][j + 1][k + 1] - psi[i][j + 1][k]) * dx[i] * dy[j] / dz[k];

                sum12 += (psi[i + 1][j + 1][k + 1] - psi[i + 1][j + 1][k]) * \
                            (psi[i + 1][j + 1][k + 1] - psi[i + 1][j + 1][k]) * dx[i] * dy[j] / dz[k];

                sum13 += V[i][j][k] * psi[i][j][k] * psi[i][j][k] * dx[i] * dy[j] * dz[k];
                sum14 += V[i + 1][j][k] * psi[i + 1][j][k] * psi[i + 1][j][k] * dx[i] * dy[j] * dz[k];
                sum15 += V[i][j + 1][k] * psi[i][j + 1][k] * psi[i][j + 1][k] * dx[i] * dy[j] * dz[k];
                sum16 += V[i][j][k + 1] * psi[i][j][k + 1] * psi[i][j][k + 1] * dx[i] * dy[j] * dz[k];
                sum17 += V[i + 1][j + 1][k] * psi[i + 1][j + 1][k] * psi[i + 1][j + 1][k] * dx[i] * dy[j] * dz[k];
                sum18 += V[i + 1][j][k + 1] * psi[i + 1][j][k + 1] * psi[i + 1][j][k + 1] * dx[i] * dy[j] * dz[k];
                sum19 += V[i][j + 1][k + 1] * psi[i][j + 1][k + 1] * psi[i][j + 1][k + 1] * dx[i] * dy[j] * dz[k];
                sum20 += V[i + 1][j + 1][k + 1] * psi[i + 1][j + 1][k + 1] * \
                            psi[i + 1][j + 1][k + 1] * dx[i] * dy[j] * dz[k];

                //слагаемые, входящие в знаменатель функционала энергии
                sum25 += psi[i][j][k] * psi[i][j][k] * dx[i] * dy[j] * dz[k];
                sum26 += psi[i + 1][j][k] * psi[i + 1][j][k] * dx[i] * dy[j] * dz[k];
                sum27 += psi[i][j + 1][k] * psi[i][j + 1][k] * dx[i] * dy[j] * dz[k];
                sum28 += psi[i][j][k + 1] * psi[i][j][k + 1] * dx[i] * dy[j] * dz[k];
                sum29 += psi[i + 1][j + 1][k] * psi[i + 1][j + 1][k] * dx[i] * dy[j] * dz[k];
                sum30 += psi[i + 1][j][k + 1] * psi[i + 1][j][k + 1] * dx[i] * dy[j] * dz[k];
                sum31 += psi[i][j + 1][k + 1] * psi[i][j + 1][k + 1] * dx[i] * dy[j] * dz[k];
                sum32 += psi[i + 1][j + 1][k + 1] * psi[i + 1][j + 1][k + 1] * dx[i] * dy[j] * dz[k];
            }
        }
    }

    sum1 = sum1 / 8.;
    sum2 = sum2 / 8.;
    sum3 = sum3 / 8.;
    sum4 = sum4 / 8.;

    sum5 = sum5 / 8.;
    sum6 = sum6 / 8.;
    sum7 = sum7 / 8.;
    sum8 = sum8 / 8.;

    sum9 = sum9 / 8.;
    sum10 = sum10 / 8.;
    sum11 = sum11 / 8.;
    sum12 = sum12 / 8.;

    sum13 = sum13 / 8.;
    sum14 = sum14 / 8.;
    sum15 = sum15 / 8.;
    sum16 = sum16 / 8.;
    sum17 = sum17 / 8.;
    sum18 = sum18 / 8.;
    sum19 = sum19 / 8.;
    sum20 = sum20 / 8.;

    sum21 = L * sum21 / 8.;
    sum22 = L * sum22 / 8.;
    sum23 = L * sum23 / 8.;
    sum24 = L * sum24 / 8.;

    sum25 = sum25 / 8.;
    sum26 = sum26 / 8.;
    sum27 = sum27 / 8.;
    sum28 = sum28 / 8.;
    sum29 = sum29 / 8.;
    sum30 = sum30 / 8.;
    sum31 = sum31 / 8.;
    sum32 = sum32 / 8.;

    E_kinetic = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 \
 + sum9 + sum10 + sum11 + sum12;
    E_potential = sum13 + sum14 + sum15 + sum16 + sum17 + sum18 \
 + sum19 + sum20;
    surface_term = sum21 + sum22 + sum23 + sum24;
    numerator = E_kinetic + E_potential + surface_term;
    denumerator = sum25 + sum26 + sum27 + sum28 + sum29 + sum30 + sum31 + sum32;

    E = numerator / denumerator;

    //printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t\n",\
     sum1, sum5, sum9, sum13, sum21, sum25, numerator, denumerator, E);
}

void clear_grad_E() {
    int i, j, k;
    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM; j++) {
            for (k = M; k < NUM2; k++) {
                grad_E[i][j][k] = 0;
            }
        }
    }
}

void calculate_grad_E() {
    double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10,
            temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20,
            temp25, temp26, temp27, temp28, temp29, temp30, temp31, temp32;
    //производные слагаемых, входящих в числитель
    double d_v1, d_v2, d_v3, d_v4, d_v5, d_v6, d_v7, d_v8, d_v9, d_v10,
            d_v11, d_v12, d_v13, d_v14, d_v15, d_v16, d_v17, d_v18, d_v19, d_v20;
    //производные слагаемых, входящих в знаменатель
    double d_w1, d_w2, d_w3, d_w4, d_w5, d_w6, d_w7, d_w8;

    int i, j, k;

    clear_grad_E();

    //ДАЛЬШЕ НАДО ПОМНИТЬ, ЧТО МЫ ЗАПИСЫВАЕМ (-1)*ГРАДИЕНТ
    //grad_E = (dw*num - dv*denum) / (denum*denum);
    // где dv - производная числителя, а dw - знаменателя

    for (i = 0; i < NUM - 1; i++) {
        for (j = 0; j < NUM - 1; j++) {
            grad_E[i][j][M] -= L * psi[i][j][M] * dx[i] * dy[j] / (4. * denumerator);
            grad_E[i + 1][j][M] -= L * psi[i + 1][j][M] * dx[i] * dy[j] / (4. * denumerator);
            grad_E[i][j + 1][M] -= L * psi[i][j + 1][M] * dx[i] * dy[j] / (4. * denumerator);
            grad_E[i + 1][j + 1][M] -= L * psi[i + 1][j + 1][M] * dx[i] * dy[j] / (4. * denumerator);

            for (k = M; k < NUM2 - 1; k++) {

                //слагаемые, входящие в производную числителя функционала энергии

                //КИНЕТИЧЕСКИЕ СЛАГАЕМЫЕ
                temp1 = (psi[i + 1][j][k] - psi[i][j][k]) * dy[j] * dz[k] / dx[i] / 4.;
                // += (i+1, j, k) // -= (i, j, k)
                d_v1 = -temp1 / denumerator;    //помним, что записыываем (-1)*(градиент энергии)
                grad_E[i + 1][j][k] += d_v1;
                grad_E[i][j][k] -= d_v1;

                temp2 = (psi[i + 1][j + 1][k] - psi[i][j + 1][k]) * dy[j] * dz[k] / dx[i] / 4.;
                // += (i+1, j+1, k) // -= (i, j+1, k)
                d_v2 = -temp2 / denumerator;
                grad_E[i + 1][j + 1][k] += d_v2;
                grad_E[i][j + 1][k] -= d_v2;

                temp3 = (psi[i + 1][j][k + 1] - psi[i][j][k + 1]) * dy[j] * dz[k] / dx[i] / 4.;
                // += (i+1, j, k+1) // -= (i, j, k+1)
                d_v3 = -temp3 / denumerator;
                grad_E[i + 1][j][k + 1] += d_v3;
                grad_E[i][j][k + 1] -= d_v3;

                temp4 = (psi[i + 1][j + 1][k + 1] - psi[i][j + 1][k + 1]) * dy[j] * dz[k] / dx[i] / 4.;
                // += (i+1, j+1, k+1) // -= (i, j+1, k+1)
                d_v4 = -temp4 / denumerator;
                grad_E[i + 1][j + 1][k + 1] += d_v4;
                grad_E[i][j + 1][k + 1] -= d_v4;
////////////////////////////////////////////////////////////////
                temp5 = (psi[i][j + 1][k] - psi[i][j][k]) * dx[i] * dz[k] / dy[j] / 4.;
                // += (i, j+1, k) // -= (i, j, k)
                d_v5 = -temp5 / denumerator;
                grad_E[i][j + 1][k] += d_v5;
                grad_E[i][j][k] -= d_v5;

                temp6 = (psi[i + 1][j + 1][k] - psi[i + 1][j][k]) * dx[i] * dz[k] / dy[j] / 4.;
                // += (i+1, j+1, k) // -= (i+1, j, k)
                d_v6 = -temp6 / denumerator;
                grad_E[i + 1][j + 1][k] += d_v6;
                grad_E[i + 1][j][k] -= d_v6;

                temp7 = (psi[i][j + 1][k + 1] - psi[i][j][k + 1]) * dx[i] * dz[k] / dy[j] / 4.;
                // += (i, j+1, k+1) // -= (i, j, k+1)
                d_v7 = -temp7 / denumerator;
                grad_E[i][j + 1][k + 1] += d_v7;
                grad_E[i][j][k + 1] -= d_v7;

                temp8 = (psi[i + 1][j + 1][k + 1] - psi[i + 1][j][k + 1]) * dx[i] * dz[k] / dy[j] / 4.;
                // += (i+1, j+1, k+1) // -= (i+1, j, k+1)
                d_v8 = -temp8 / denumerator;
                grad_E[i + 1][j + 1][k + 1] += d_v8;
                grad_E[i + 1][j][k + 1] -= d_v8;
/////////////////////////////////////////////////////////////////
                temp9 = (psi[i][j][k + 1] - psi[i][j][k]) * dx[i] * dy[j] / dz[k] / 4.;
                // += (i, j, k+1) // -= (i, j, k)
                d_v9 = -temp9 / denumerator;
                grad_E[i][j][k + 1] += d_v9;
                grad_E[i][j][k] -= d_v9;

                temp10 = (psi[i + 1][j][k + 1] - psi[i + 1][j][k]) * dx[i] * dy[j] / dz[k] / 4.;
                // += (i+1, j, k+1) // -= (i+1, j, k)
                d_v10 = -temp10 / denumerator;
                grad_E[i + 1][j][k + 1] += d_v10;
                grad_E[i + 1][j][k] -= d_v10;

                temp11 = (psi[i][j + 1][k + 1] - psi[i][j + 1][k]) * dx[i] * dy[j] / dz[k] / 4.;
                // += (i, j+1, k+1) // -= (i, j+1, k)
                d_v11 = -temp11 / denumerator;
                grad_E[i][j + 1][k + 1] += d_v11;
                grad_E[i][j + 1][k] -= d_v11;

                temp12 = (psi[i + 1][j + 1][k + 1] - psi[i + 1][j + 1][k]) * dx[i] * dy[j] / dz[k] / 4.;
                // += (i+1, j+1, k+1) // -= (i+1, j+1, k)
                d_v12 = -temp12 / denumerator;
                grad_E[i + 1][j + 1][k + 1] += d_v12;
                grad_E[i + 1][j + 1][k] -= d_v12;
////////////////////////////////////////////////////////////////
                //ПОТЕНЦИАЛЬНЫЕ СЛАГАЕМЫЕ
                temp13 = V[i][j][k] * psi[i][j][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j, k)
                d_v13 = -temp13 / denumerator;
                grad_E[i][j][k] += d_v13;

                temp14 = V[i + 1][j][k] * psi[i + 1][j][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j, k)
                d_v14 = -temp14 / denumerator;
                grad_E[i + 1][j][k] += d_v14;

                temp15 = V[i][j + 1][k] * psi[i][j + 1][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j+1, k)
                d_v15 = -temp15 / denumerator;
                grad_E[i][j + 1][k] += d_v15;

                temp16 = V[i][j][k + 1] * psi[i][j][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j, k+1)
                d_v16 = -temp16 / denumerator;
                grad_E[i][j][k + 1] += d_v16;

                temp17 = V[i + 1][j + 1][k] * psi[i + 1][j + 1][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j+1, k)
                d_v17 = -temp17 / denumerator;
                grad_E[i + 1][j + 1][k] += d_v17;

                temp18 = V[i + 1][j][k + 1] * psi[i + 1][j][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j, k+1)
                d_v18 = -temp18 / denumerator;
                grad_E[i + 1][j][k + 1] += d_v18;

                temp19 = V[i][j + 1][k + 1] * psi[i][j + 1][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j+1, k+1)
                d_v19 = -temp19 / denumerator;
                grad_E[i][j + 1][k + 1] += d_v19;

                temp20 = V[i + 1][j + 1][k + 1] * psi[i + 1][j + 1][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j+1, k+1)
                d_v20 = -temp20 / denumerator;
                grad_E[i + 1][j + 1][k + 1] += d_v20;
////////////////////////////////////////////////////////////////
                //слагаемые, входящие в производную знаменателя функционала энергии
                temp25 = psi[i][j][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j, k)
                d_w1 = temp25 * numerator / denumerator;
                grad_E[i][j][k] += d_w1;

                temp26 = psi[i + 1][j][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j, k)
                d_w2 = temp26 * numerator / denumerator;
                grad_E[i + 1][j][k] += d_w2;

                temp27 = psi[i][j + 1][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j+1, k)
                d_w3 = temp27 * numerator / denumerator;
                grad_E[i][j + 1][k] += d_w3;

                temp28 = psi[i][j][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j, k+1)
                d_w4 = temp28 * numerator / denumerator;
                grad_E[i][j][k + 1] += d_w4;

                temp29 = psi[i + 1][j + 1][k] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j+1, k)
                d_w5 = temp29 * numerator / denumerator;
                grad_E[i + 1][j + 1][k] += d_w5;

                temp30 = psi[i + 1][j][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j, k+1)
                d_w6 = temp30 * numerator / denumerator;
                grad_E[i + 1][j][k + 1] += d_w6;

                temp31 = psi[i][j + 1][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i, j+1, k+1)
                d_w7 = temp31 * numerator / denumerator;
                grad_E[i][j + 1][k + 1] += d_w7;

                temp32 = psi[i + 1][j + 1][k + 1] * dx[i] * dy[j] * dz[k] / 4.;
                // += (i+1, j+1, k+1)
                d_w8 = temp32 * numerator / denumerator;
                grad_E[i + 1][j + 1][k + 1] += d_w8;

            }
        }
    }
}

//нормировка волновой функции на 1
void normalize_to_1() {
    int i, j, k;
    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM; j++) {
            for (k = M; k < NUM2; k++) {
                psi[i][j][k] = psi[i][j][k] / sqrt(denumerator);
            }
        }
    }
    calculate_E();
}

int main() {
    int i, j, k, count_steps;
    time_t sec;
    double F, v;

    time_t tstart = time(NULL);

    //переменный шаг по x и y
    x[0] = y[0] = step1;
    for (i = 1; i < NUM; i++) {
        x[i] = y[i] = x[i - 1] + step1 * (1 + (i - 1) * step_x);
    }

    //постоянный шаг по z
    z[NUM2 / 2] = c * step_z;
    for (j = NUM2 / 2 + 1; j < NUM2; j++) {
        z[j] = z[j - 1] + step_z;
    }

    z[NUM2 / 2 - 1] = -c * step_z;
    for (k = NUM2 / 2 - 2; k >= M; k--) {
        z[k] = z[k + 1] - step_z;
    }
////////////////////////////////////////////////////////////

    //printf("%lf\t %lf\t %lf\t %lf\t %lf\t\n", x[0], x[NUM - 1], z[M], z[NUM2 / 2], z[NUM2 - 1]);

    prepare_functions();

    count_steps = 0;

    calculate_E();
    calculate_grad_E();

    printf("E			\t E_kinetic		\t E_potential		\t max_grad_E		\t\n");
    printf("%.14lf\t	%.14lf\t	%.14lf\t	%.14le\t\n", E, E_kinetic, E_potential, max_grad_E);

    max_grad_E = grad_E[0][0][M];
    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM; j++) {
            for (k = M; k < NUM2; k++) {
                if (max_grad_E < fabs(grad_E[i][j][k]))
                    max_grad_E = fabs(grad_E[i][j][k]);
            }
        }
    }

    normalize_to_1();
    printf("%.14lf\t	%.14lf\t	%.14lf\t	%.14le\t\n", E, E_kinetic, E_potential, max_grad_E);

    while (max_grad_E > 1.e-10) {

        for (i = 0; i < NUM; i++) {
            for (j = 0; j < NUM; j++) {
                for (k = M; k < NUM2; k++) {

                    vel[i][j][k] = G * (vel[i][j][k] + grad_E[i][j][k] * dt);
                    //уже учтено, что grad_E на самом деле (-1)*ГРАДИЕНТ

                    psi[i][j][k] = psi[i][j][k] + vel[i][j][k] * dt;
                }
            }
        }

        calculate_E();
        calculate_grad_E();

        max_grad_E = grad_E[0][0][M];
        for (i = 0; i < NUM; i++) {
            for (j = 0; j < NUM; j++) {
                for (k = M; k < NUM2; k++) {
                    if (max_grad_E < fabs(grad_E[i][j][k]))
                        max_grad_E = fabs(grad_E[i][j][k]);
                }
            }
        }

        if (count_steps == 100) {
            normalize_to_1();
            printf("%.14lf\t	%.14lf\t	%.14lf\t	%.14le\t\n",
                   E, E_kinetic, E_potential, max_grad_E);
            count_steps = 0;
        }
        count_steps++;
    }

    printf("%.14lf\t	%.14lf\t	%.14lf\t	%.14le\t\n", E, E_kinetic, E_potential, max_grad_E);
    normalize_to_1();
    printf("%.14lf\t	%.14lf\t	%.14lf\t	%.14le\t\n", E, E_kinetic, E_potential, max_grad_E);

    time_t tstop = time(NULL);
    sec = tstop - tstart;
    printf("spent_time = %lld\n\a", sec);

    //getchar();

    return 0;
}