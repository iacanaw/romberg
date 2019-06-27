#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define MAX_STEPS 11
#define NTHREADS (int)floor((MAX_STEPS+1)/2)
#define LIM_INF 0.001
#define LIM_SUP 10.001

double R[MAX_STEPS][MAX_STEPS];
double h = (LIM_SUP-LIM_INF); // step size

void *my_thread(void *x){
    int t_id;
    t_id = *((int *) x);
    printf("Thread %i inicializada!\n", t_id);
    return NULL;
}

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

   dump_row(0, Rp);

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

      /* if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc){
         return Rc[i-1];
      }*/

      //swap Rn and Rc as we only need the last row
      double *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[max_steps-1]; //return our best guess
}

double romberg_paralel(){
    pthread_t threads[NTHREADS];
    int thread_args[NTHREADS];
    int rc, idx, i, j, n_col;
    double result;
    int index[MAX_STEPS];

    j = 1;

    for(i=0; i<MAX_STEPS; i++){
        index[i] = j;
        printf("index[%d] = %d\n", i, index[i]);
        j--;
    }

    // Calcula o primeiro elemento
    R[0][0] = (function(LIM_INF) + function(LIM_SUP))*h*.5;

    // Calcula os elementos
    for(i=0; i<=MAX_STEPS*2; i++){
        // perorrer os indexs
        printf("\n\ninteracao %d\n", i+1);
        for(j=0; j<MAX_STEPS; j++){
            if(index[j]==0){
                //calcula R[j][index[j]];
                printf("R[%d][%d]\n", j, index[j]);
                index[j]++;
            }
            else if(index[j]<=j){
                if (index[j]>0){
                    //calcula R[j][index[j]];
                    printf("R[%d][%d]\n", j, index[j]);
                }
                index[j]++;
            }
            else if(index[j] < 0){
                index[j]++;
            }
        }
    }

    return result;
}

int main(int argc, char *argv[]){
    double resultado1, resultado2;

    resultado1 = romberg(function, 0.01, 50, 11, 0.00001);
    printf("resultado simples: %f\n",resultado1);
    resultado2 = romberg_paralel();
    printf("resultado paralelo: %f\n",resultado2);
    return 1;
}
