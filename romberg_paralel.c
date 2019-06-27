#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define MAX_STEPS 10
#define NTHREADS 5//(int)floor((MAX_STEPS+1)/2)
#define LIM_INF 0.001
#define LIM_SUP 10.001
#define ACCURACY 0.0001 // 1E-8

typedef struct{
    int row;
    int collum;
    double limInf;
} thread_args;

thread_args args[NTHREADS];
double R[MAX_STEPS][MAX_STEPS];
double h[MAX_STEPS]; // step size


void dump_row(size_t i, double *R){
   printf("R[%li] = ", i);
   for(size_t j = 0; j <= i; ++j){
      printf("%f ", R[j]);
   }
   printf("\n");
}

double function(double x){
    double result;
    result = x * log(x);
    result = pow(2.71828182846, result);
    result = log10(result)/log10(4);
    result = result * sin(x);
    return result;
}

double romberg(double (*f/* function to integrate */)(double), double /*lower limit*/ a, double /*upper limit*/ b, size_t max_steps, double /*desired accuracy*/ acc){
   double R1[max_steps], R2[max_steps]; //buffers
   double *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   Rp[0] = (f(a) + f(b))*h*.5; //first trapezoidal step

   //dump_row(0, Rp);

   for(size_t i = 1; i < max_steps; ++i){
      h /= 2.;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(size_t j = 1; j <= ep; ++j){
         c += f(a+(2*j-1)*h);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(size_t j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
      dump_row(i, Rc);

     /* if(i > 1 && fabs(Rp[i-1]-Rc[i]) < ACCURACY){
         return Rc[i-1];
      }*/

      //swap Rn and Rc as we only need the last row
      double *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[max_steps-1]; //return our best guess
}

void *romberg_collumZero(void *x){
    double c = 0.;
    double a;
    int row, collum, i, j, limInf;
    int id;
    id = *(int *)x;
    //printf("thread id: %d\n", id);
    
    row = args[id].row;
    collum = args[id].collum;
    limInf = args[id].limInf;
    //printf("R[%d][%d]\n", row, collum);

    size_t ep = 1 << (row-1); //2^(n-1)
    for(size_t j = 1; j <= ep; ++j){
        c += function(LIM_INF+(2*j-1)*h[row]);
    }

    R[row][collum] = h[row]*c + .5*R[row-1][collum];
    //printf("R[%d][%d] = %f\n", row, collum,R[row][collum]);
    return NULL;
}

void *romberg_collumN(void *x){
    int row, collum, j;
    double n_k;
    int id;
    id = *(int *)x;

    row = args[id].row;
    collum = args[id].collum;
    //printf("R[%d][%d]\n", row, collum);

    n_k = pow(4, collum);
    R[row][collum] = (n_k*R[row][collum-1] - R[row-1][collum-1])/(n_k-1); //compute R(i,j)
    //printf("R[%d][%d] = %f\n", row, collum,R[row][collum]);
    return NULL;
}

void *doNothing(void *x){
    //printf("");
    return NULL;
}

double romberg_paralel(double (*f )(double), double limInf, double limSup, size_t max_steps){
    pthread_t threads[NTHREADS];
    int tidx[NTHREADS];
    int rc, i, j, k, thread_idx;
    double result, erro;
    int index[max_steps];
    for(i = 0; i< NTHREADS; i++){
        args[i].limInf = LIM_INF;
    }

    // Inicia o controle (index) e preenche os H's
    j = 1;
    h[0] = (limSup - limInf);
    for(i=0; i<max_steps; i++){
        index[i] = j;
        //preenche os H's
        if (i > 0){h[i] = h[i-1]/2.0;} 
        j--;
    }

    // Calcula o primeiro elemento
    R[0][0] = (function(limInf) + function(limSup))*h[0]*.5;
    //printf("R[%d][%d] = %f\n", 0, 0,R[0][0]);


    // Calcula os elementos
    for(i=0; i<(max_steps-1)*2; i++){
        // percorrer os indexs
        //printf("\n\ninteracaoa %d\n", i);
        for(j=0; j<max_steps; j++){
            thread_idx = 0;
            if(index[j]==0){
                args[thread_idx].row = j;
                args[thread_idx].collum = index[j];
                tidx[thread_idx] = thread_idx;
                pthread_create(&threads[thread_idx], NULL, romberg_collumZero, &tidx[thread_idx]);
                //pthread_join(threads[thread_idx], NULL);
                thread_idx++;
                //calcula R[j][index[j]];
                index[j]++;
            }
            else if(index[j]<=j){
                if (index[j]>0){
                    args[thread_idx].row = j;
                    args[thread_idx].collum = index[j];
                    tidx[thread_idx] = thread_idx;
                    pthread_create(&threads[thread_idx], NULL, romberg_collumN, &tidx[thread_idx]);
                    //pthread_join(threads[thread_idx], NULL);
                    thread_idx++;
                    //calcula R[j][index[j]];

                    // Critério de parada
                    
                    /*if(j == index[j] && j>1){
                        printf("Calculando erro!\n");
                        erro = fabs((R[j-1][index[j]-1])-(R[j][index[j]]));
                        if(erro < ACCURACY){
                            printf("~~~~Parou em R[%d][%d]\n",j, index[j]);
                            return R[j][index[j]];
                        }
                    }*/
                }
                index[j]++;
            }
            else if(index[j] < 0){
                index[j]++;
            }
	    while(thread_idx < NTHREADS){
            pthread_create(&threads[thread_idx], NULL, doNothing, &tidx[thread_idx]);
            thread_idx++;
	    }
            for(k=0; k<NTHREADS; k++){
                pthread_join(threads[k], NULL);
            }
        }
    }
    result = R[MAX_STEPS-1][MAX_STEPS-1];
    return result;
}

int main(int argc, char *argv[]){
    double resultado1, resultado2;
    time_t seconds1, seconds2, seconds3;
    clock_t CPU_time_1, CPU_time_2, CPU_time_3, CPU_time_4;

    // Execução Sequencial
    seconds1 = time(NULL);
    CPU_time_1 = clock();
    resultado1 = romberg(function, LIM_INF, LIM_SUP, MAX_STEPS, 0.00001);
    CPU_time_2 = clock();
    // Execucao Paralela
    seconds2 = time(NULL);
    CPU_time_3 = clock();
    resultado2 = romberg_paralel(function, LIM_INF, LIM_SUP, MAX_STEPS);
    seconds3 = time(NULL);
    CPU_time_4 = clock();

    // Mostra resutados
    printf("resultado paralelo: %f\n",resultado2);
    printf("resultado simples: %f\n",resultado1);
    printf("Tempo total seq: %ld -- tempo total para: %ld\n", seconds2 - seconds1, seconds3 - seconds2);
    printf("Clocks seq: %ld -- Clocks para: %ld\n", CPU_time_2 - CPU_time_1, CPU_time_4 - CPU_time_3);
    printf("Ganho: %f", -1.0*((((float)CPU_time_4 - (float)CPU_time_3) * 100.0) / ((float)CPU_time_2 - (float)CPU_time_1)-100.0));
    return 1;
}
