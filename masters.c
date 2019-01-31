#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "output.txt"

#define H 10.   //характерный размер квантовой системы
//(пространственная бесконечность)

//#define a 3.    //расстояние, на которое развели ядра
//#define L 0.   // lambda parameter

#define NUM 100
#define NUM2 (2*NUM)

//шаги в левую и правую сторону от ядра по оси x соответственно
#define step1_x_left ((a / 2.) / (i0_x + ((i0_x-2) * (i0_x- 1) * step_x_left) / 2.))
#define step1_x_right ((H - a/2.) / (((NUM - i0_x) - 0.5) + ( ((NUM - i0_x)-2) * step_x_right)* ((NUM - i0_x) - 1) / 2.))

#define step1_y (H / ((NUM - 0.5) + ((NUM-2) * (NUM - 1) * step_y) / 2.))  // расстояние от нуля до первой точки
//(использована формула суммы арифметической прогрессии) (сумма слагаемых дает H)
#define c 0.699 //коэффициент определяющий размер до ближайшей к атому точку вдоль оси z

//#define z0 5.   ///на сколько поднимаем плоскость (относительно её начального положения z = -a)
//физическое расстояние между атомами и плоскостью = (a-z0)

#define i0_x ((int)((a / 2.) * (NUM - 1.) / H))   //индекс x узла, в котором находится ядро
#define M ((int)(z0 / step_z))    //индекс, начиная с которого заполняется массив точек z

#define H_var   4.  //чем больше параметр(переменного шага),
// тем меньше будут первые шаги вдоль соответстующей оси
// (график будет сильнее прогибаться)

#define step_x_left (H_var / (i0_x - 2))
#define step_x_right (H_var / ((NUM - i0_x) - 2))
#define step_y (H_var/(NUM - 2))
#define step_z (H/(NUM-1+c))

//#define G 0.995
#define dt 2.e-1

//##### GLOBAL VARIABLES    #####
double a;    //расстояние, на которое развели ядра
double L;   // lambda parameter
double z0;   ///на сколько поднимаем плоскость (относительно её начального положения z = -a)
//физическое расстояние между атомами и плоскостью = (a-z0)

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
                V[i][j][k] = -1. / (sqrt((x[i] + a / 2.) * (x[i] + a / 2.) + y[j] * y[j] + z[k] * z[k])) - \
                                1. / (sqrt((x[i] - a / 2.) * (x[i] - a / 2.) + y[j] * y[j] + z[k] * z[k]));

                psi[i][j][k] = exp(-sqrt((x[i] - a/2) * (x[i] - a/2) + y[j] * y[j] + z[k] * z[k]));
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

    E_kinetic = (sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 \
 + sum9 + sum10 + sum11 + sum12) / 8.;
    E_potential = (sum13 + sum14 + sum15 + sum16 + sum17 + sum18 \
 + sum19 + sum20) / 8.;
    surface_term = L*(sum21 + sum22 + sum23 + sum24) / 8.;
    numerator = E_kinetic + E_potential + surface_term;
    denumerator = (sum25 + sum26 + sum27 + sum28 + sum29 + sum30 + sum31 + sum32) / 8.;

    E = numerator / denumerator;

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
            grad_E[i][j][M] -= L * psi[i][j][M] * dx[i] * dy[j];
            grad_E[i + 1][j][M] -= L * psi[i + 1][j][M] * dx[i] * dy[j];
            grad_E[i][j + 1][M] -= L * psi[i][j + 1][M] * dx[i] * dy[j];
            grad_E[i + 1][j + 1][M] -= L * psi[i + 1][j + 1][M] * dx[i] * dy[j];

            for (k = M; k < NUM2 - 1; k++) {

                //слагаемые, входящие в производную числителя функционала энергии

                //КИНЕТИЧЕСКИЕ СЛАГАЕМЫЕ
                // += (i+1, j, k) // -= (i, j, k)
                d_v1 = (psi[i + 1][j][k] - psi[i][j][k]) * dy[j] * dz[k] / dx[i];;
                //помним, что записыываем (-1)*(градиент энергии)
                grad_E[i + 1][j][k] -= d_v1;
                grad_E[i][j][k] += d_v1;

                // += (i+1, j+1, k) // -= (i, j+1, k)
                d_v2 = (psi[i + 1][j + 1][k] - psi[i][j + 1][k]) * dy[j] * dz[k] / dx[i];
                grad_E[i + 1][j + 1][k] -= d_v2;
                grad_E[i][j + 1][k] += d_v2;

                // += (i+1, j, k+1) // -= (i, j, k+1)
                d_v3 = (psi[i + 1][j][k + 1] - psi[i][j][k + 1]) * dy[j] * dz[k] / dx[i];
                grad_E[i + 1][j][k + 1] -= d_v3;
                grad_E[i][j][k + 1] += d_v3;

                // += (i+1, j+1, k+1) // -= (i, j+1, k+1)
                d_v4 = (psi[i + 1][j + 1][k + 1] - psi[i][j + 1][k + 1]) * dy[j] * dz[k] / dx[i];
                grad_E[i + 1][j + 1][k + 1] -= d_v4;
                grad_E[i][j + 1][k + 1] += d_v4;
////////////////////////////////////////////////////////////////
                // += (i, j+1, k) // -= (i, j, k)
                d_v5 = (psi[i][j + 1][k] - psi[i][j][k]) * dx[i] * dz[k] / dy[j];
                grad_E[i][j + 1][k] -= d_v5;
                grad_E[i][j][k] += d_v5;

                // += (i+1, j+1, k) // -= (i+1, j, k)
                d_v6 = (psi[i + 1][j + 1][k] - psi[i + 1][j][k]) * dx[i] * dz[k] / dy[j];
                grad_E[i + 1][j + 1][k] -= d_v6;
                grad_E[i + 1][j][k] += d_v6;

                // += (i, j+1, k+1) // -= (i, j, k+1)
                d_v7 = (psi[i][j + 1][k + 1] - psi[i][j][k + 1]) * dx[i] * dz[k] / dy[j];
                grad_E[i][j + 1][k + 1] -= d_v7;
                grad_E[i][j][k + 1] += d_v7;

                // += (i+1, j+1, k+1) // -= (i+1, j, k+1)
                d_v8 = (psi[i + 1][j + 1][k + 1] - psi[i + 1][j][k + 1]) * dx[i] * dz[k] / dy[j];
                grad_E[i + 1][j + 1][k + 1] -= d_v8;
                grad_E[i + 1][j][k + 1] += d_v8;
/////////////////////////////////////////////////////////////////
                // += (i, j, k+1) // -= (i, j, k)
                d_v9 = (psi[i][j][k + 1] - psi[i][j][k]) * dx[i] * dy[j] / dz[k];
                grad_E[i][j][k + 1] -= d_v9;
                grad_E[i][j][k] += d_v9;

                // += (i+1, j, k+1) // -= (i+1, j, k)
                d_v10 = (psi[i + 1][j][k + 1] - psi[i + 1][j][k]) * dx[i] * dy[j] / dz[k];
                grad_E[i + 1][j][k + 1] -= d_v10;
                grad_E[i + 1][j][k] += d_v10;

                // += (i, j+1, k+1) // -= (i, j+1, k)
                d_v11 = (psi[i][j + 1][k + 1] - psi[i][j + 1][k]) * dx[i] * dy[j] / dz[k];
                grad_E[i][j + 1][k + 1] -= d_v11;
                grad_E[i][j + 1][k] += d_v11;

                // += (i+1, j+1, k+1) // -= (i+1, j+1, k)
                d_v12 = (psi[i + 1][j + 1][k + 1] - psi[i + 1][j + 1][k]) * dx[i] * dy[j] / dz[k];
                grad_E[i + 1][j + 1][k + 1] -= d_v12;
                grad_E[i + 1][j + 1][k] += d_v12;
////////////////////////////////////////////////////////////////
                //ПОТЕНЦИАЛЬНЫЕ СЛАГАЕМЫЕ
                // += (i, j, k)
                d_v13 = V[i][j][k] * psi[i][j][k] * dx[i] * dy[j] * dz[k];
                grad_E[i][j][k] -= d_v13;

                // += (i+1, j, k)
                d_v14 = V[i + 1][j][k] * psi[i + 1][j][k] * dx[i] * dy[j] * dz[k];
                grad_E[i + 1][j][k] -= d_v14;

                // += (i, j+1, k)
                d_v15 = V[i][j + 1][k] * psi[i][j + 1][k] * dx[i] * dy[j] * dz[k];
                grad_E[i][j + 1][k] -= d_v15;

                // += (i, j, k+1)
                d_v16 = V[i][j][k + 1] * psi[i][j][k + 1] * dx[i] * dy[j] * dz[k];
                grad_E[i][j][k + 1] -= d_v16;

                // += (i+1, j+1, k)
                d_v17 = V[i + 1][j + 1][k] * psi[i + 1][j + 1][k] * dx[i] * dy[j] * dz[k];
                grad_E[i + 1][j + 1][k] -= d_v17;

                // += (i+1, j, k+1)
                d_v18 = V[i + 1][j][k + 1] * psi[i + 1][j][k + 1] * dx[i] * dy[j] * dz[k];
                grad_E[i + 1][j][k + 1] -= d_v18;

                // += (i, j+1, k+1)
                d_v19 = V[i][j + 1][k + 1] * psi[i][j + 1][k + 1] * dx[i] * dy[j] * dz[k];
                grad_E[i][j + 1][k + 1] -= d_v19;

                // += (i+1, j+1, k+1)
                d_v20 = V[i + 1][j + 1][k + 1] * psi[i + 1][j + 1][k + 1] * dx[i] * dy[j] * dz[k];
                grad_E[i + 1][j + 1][k + 1] -= d_v20;
////////////////////////////////////////////////////////////////
                //слагаемые, входящие в производную знаменателя функционала энергии
                // += (i, j, k)
                d_w1 = psi[i][j][k] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i][j][k] += d_w1;

                // += (i+1, j, k)
                d_w2 = psi[i + 1][j][k] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i + 1][j][k] += d_w2;

                // += (i, j+1, k)
                d_w3 = psi[i][j + 1][k] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i][j + 1][k] += d_w3;

                // += (i, j, k+1)
                d_w4 = psi[i][j][k + 1] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i][j][k + 1] += d_w4;

                // += (i+1, j+1, k)
                d_w5 = psi[i + 1][j + 1][k] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i + 1][j + 1][k] += d_w5;

                // += (i+1, j, k+1)
                d_w6 = psi[i + 1][j][k + 1] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i + 1][j][k + 1] += d_w6;

                // += (i, j+1, k+1)
                d_w7 = psi[i][j + 1][k + 1] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i][j + 1][k + 1] += d_w7;

                // += (i+1, j+1, k+1)
                d_w8 = psi[i + 1][j + 1][k + 1] * dx[i] * dy[j] * dz[k] * numerator / denumerator;
                grad_E[i + 1][j + 1][k + 1] += d_w8;

            }
        }
    }


    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM; j++) {
            for (k = M; k < NUM2; k++) {
                grad_E[i][j][k] = grad_E[i][j][k] / 4. /
                                  denumerator;
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

//функция, считающая косинус между векторами и возвращающая коэффициент трения
double calc_G() {
    int i, j, k;
    double dot_ab, norm_a, norm_b, cos_ab, ans;
    dot_ab = norm_a = norm_b = 0;

    for (i = 0; i < NUM; i++) {
        for (j = 0; j < NUM; j++) {
            for (k = M; k < NUM2; k++) {
                dot_ab += vel[i][j][k] * grad_E[i][j][k];
                norm_a += vel[i][j][k] * vel[i][j][k];
                norm_b += grad_E[i][j][k] * grad_E[i][j][k];
            }
        }
    }

    norm_a = sqrt(norm_a);
    norm_b = sqrt(norm_b);

    if (norm_a * norm_b == 0) {
        cos_ab = 0;
    }
    else {
        cos_ab = dot_ab / (norm_a * norm_b);
    }
    ans = (1. + 0.1 * tanh((cos_ab - 0.5) * 2.) - 0.1 * tanh((1. - 0.5) * 2.));

   return ans;
}

int main() {
    int i, j, k, count_steps, count_prints;
    time_t sec;
    double G;
    //double F, v;

    FILE *f;
    f = fopen(OUTFILE, "w");
    fprintf(f, "L\t z0\t a\t	E		\tE_kinetic	\tE_potential	\tmax_grad_E\t\n");

    //Инициализация параметров задачи и тело программы
    L = 0.;
    z0 = 10.;
    while (z0 > step_z) {
        a = 10.;
        while (a >= 0.) {

            time_t tstart = time(NULL);
            //переменный шаг по x
            x[i0_x] = a / 2. + step1_x_right / 2.;   //делим пополам, чтобы наши точки были тождественны узлам решетки
            for (i = i0_x + 1; i < NUM; i++) {
                x[i] = x[i - 1] + step1_x_right * (1 + (i - (i0_x + 1)) * step_x_right);
            }

            if (i0_x > 2) {
                x[i0_x - 1] =
                        a / 2. - step1_x_left / 2.;   //делим пополам, чтобы наши точки были тождественны узлам решетки
                for (i = i0_x - 2; i >= 0; i--) {
                    x[i] = x[i + 1] - step1_x_left * (1 - (i - (i0_x - 2)) * step_x_left);
                }
            }
            else if (i0_x > 0) { //т.е. здесь обрабатывается максимум две точки
                x[i0_x - 1] = a / 2. - (a / 2.) / (i0_x) / 2.;   //делим пополам, чтобы наши точки были тождественны узлам решетки
                for (i = i0_x - 2; i >= 0; i--) {
                    x[i] = x[i + 1] - (a / 2.) / (i0_x);
                }

            }



            //переменный шаг по y
            y[0] = step1_y / 2.;   //делим пополам, чтобы наши точки были тождественны узлам решетки
            for (j = 1; j < NUM; j++) {
                y[j] = y[j - 1] + step1_y * (1 + (j - 1) * step_y);
            }

            //постоянный шаг по z
            z[NUM2 / 2] = c * step_z;
            for (k = NUM2 / 2 + 1; k < NUM2; k++) {
                z[k] = z[k - 1] + step_z;
            }

            z[NUM2 / 2 - 1] = -c * step_z;
            for (k = NUM2 / 2 - 2; k >= M; k--) {
                z[k] = z[k + 1] - step_z;
            }
            ////////////////////////////////////////////////////////////
            if(i0_x > 2) {
                printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t\n", x[0], step1_x_left / 2., step1_x_right / .2,
                       x[i0_x + 1],
                       x[NUM - i0_x - 1], x[NUM - 1]);
            }
            else{
                printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t\n", x[0], (a / 2.) / (i0_x) / 2., step1_x_right / .2,
                       x[i0_x + 1],
                       x[NUM - i0_x - 1], x[NUM - 1]);
            }

            prepare_functions();

            count_steps = 0;
            count_prints = 1;

            calculate_E();
            calculate_grad_E();
            G = calc_G();

            printf("N	\tE		\tE_kinetic	\tE_potential	\tmax_grad_E\t\n");
            printf("%d\t	%.12lf\t	%.12lf\t	%.12lf\t	%.12le\t\n", count_prints, E, E_kinetic, E_potential,
                   max_grad_E);
            count_prints++;

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
            printf("%d\t	%.12lf\t	%.12lf\t	%.12lf\t	%.12le\t\n", count_prints, E, E_kinetic, E_potential,
                   max_grad_E);
            count_prints++;

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
                G = calc_G();

            normalize_to_1();

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
                    //normalize_to_1();
                    printf("%d\t	%.12lf\t	%.12lf\t	%.12lf\t	%.12le\t\n", count_prints,
                           E, E_kinetic, E_potential, max_grad_E);
                    count_prints++;
                    count_steps = 0;
fflush(stdout);

                }
                count_steps++;
            }

            printf("%.12lf\t	%.12lf\t	%.12lf\t	%.12le\t\n", E, E_kinetic, E_potential, max_grad_E);
            normalize_to_1();
            printf("%.12lf\t	%.12lf\t	%.12lf\t	%.12le\t\n", E, E_kinetic, E_potential, max_grad_E);

            time_t tstop = time(NULL);
            sec = tstop - tstart;
            printf("spent_time = %lld\n\a", sec);


            fprintf(f, "%.1f\t %.1f\t %.1f\t	%.12lf\t	%.12lf\t	%.12lf\t	%.12le\t\n", L, z0, a, E, E_kinetic,
                    E_potential, max_grad_E);

fflush(stdout);
fflush(f);

            a -= 0.5;
        }
        z0 -= 0.5;
    }
    fclose(f);
    //getchar();

    return 0;
}