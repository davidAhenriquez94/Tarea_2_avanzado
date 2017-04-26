#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float energy(float pressure, float density, float velocity);
float F1(float w1, float w2, float w3);
float F2(float w1, float w2, float w3);
float F3(float w1, float w2, float w3);
void Lax_Wendroff(float* w1, float* w2, float* w3);
float time_step_calculator(float cfl,float* w1, float* w2, float* w3,float dx);

int N = 1000;
float dx  = 1.0/1000;

int main(void){
  float* w_1 = malloc((N+1)*sizeof(float));
  float* w_2 = malloc((N+1)*sizeof(float));
  float* w_3 = malloc((N+1)*sizeof(float));
     

  for(int i = 0; i<(N+1); i++){
    if(i<(N/2+1)){
      w_1[i]  = 1.0;
      w_2[i]  = 0.0;
      w_3[i]  = energy(1.0,1.0,0.0);
    }
    else{
      w_1[i] = 0.1;
      w_2[i]  = 0.0;
      w_3[i]  = energy(0.1,0.1,0.0);
    }
  }

  Lax_Wendroff(w_1,w_2,w_3);

  for(int i = 0; i<(N+1); i++){
    printf("%f\t%f\t%f\n",w_1[i],((2.0/3)*(w_3[i]-(w_2[i]*w_2[i])/(2*w_1[i]))),w_2[i]/w_1[i]);
  }

  free(w_1);
  free(w_2);
  free(w_3);

}

float energy(float pressure, float density, float velocity){
  return ((3.0/2)*pressure + (density/2.0)*velocity*velocity);
}

float F1(float w1, float w2, float w3){
  return w2;
}
float F2(float w1, float w2, float w3){
  return ((w2*w2)/w1) + (2.0/3)*(w3 - ((w2*w2)/(2.0*w1)));
}    
float F3(float w1, float w2, float w3){
  return (w3 + (2.0/3)*(w3 - ((w2*w2)/(2.0*w1))))*(w2/w1);
}
float max_calculator(float *arr){
  float max = arr[0];
  for(int i = 1; i<N; i++){
    if(arr[i]>max){
      max  = arr[i];
    }
  }
  return max;
}
float time_step_calculator(float cfl,float* w1, float* w2, float* w3,float dx){
  float w_dt[N]; 
  for(int i = 0; i<N; i++){
    float t = (float) fabs(w2[i]/w1[i]) + pow((10.0/9)*((w3[i]/w1[i]) - (w2[i]*w2[i])/(2*w1[i]*w1[i])),0.5);
    w_dt[i] = t; 
  }
  return cfl*(dx/max_calculator(w_dt));
}

void Lax_Wendroff(float* w1, float* w2, float* w3){
  float* temp_1 = malloc((N-1)*sizeof(float));
  float w_medios_1_f  = 0;
  float w_medios_1_p = 0;
  float* temp_2 = malloc((N-1)*sizeof(float));
  float w_medios_2_f = 0;
  float w_medios_2_p  = 0;
  float* temp_3 = malloc((N-1)*sizeof(float));
  float w_medios_3_f  = 0;
  float w_medios_3_p = 0;
  
  int n = 1;
  while( n < 5850 ){
    
    float dt = time_step_calculator(0.1,w1,w2,w3,dx); 

    for(int j = 0; j<(N-1); j++ ){
      w_medios_1_p = ((w1[j]+w1[j+1])/2.0) - (1.0/2)*(dt/dx)*(F1(w1[j+1],w2[j+1],w3[j+1]) - F1(w1[j],w2[j],w3[j]));
      w_medios_1_f = ((w1[j+1]+w1[j+2])/2.0) - (1.0/2)*(dt/dx)*(F1(w1[j+2],w2[j+2],w3[j+2]) - F1(w1[j+1],w2[j+1],w3[j+1]));
      
      w_medios_2_p = ((w2[j]+w2[j+1])/2.0) - (1.0/2)*(dt/dx)*(F2(w1[j+1],w2[j+1],w3[j+1]) - F2(w1[j],w2[j],w3[j]));
      w_medios_2_f = ((w2[j+1]+w2[j+2])/2.0) - (1.0/2)*(dt/dx)*(F2(w1[j+2],w2[j+2],w3[j+2]) - F2(w1[j+1],w2[j+1],w3[j+1]));
      
      w_medios_3_p = ((w3[j]+w3[j+1])/2.0) - (1.0/2)*(dt/dx)*(F3(w1[j+1],w2[j+1],w3[j+1]) - F3(w1[j],w2[j],w3[j]));
      w_medios_3_f = ((w3[j+1]+w3[j+2])/2.0) - (1.0/2)*(dt/dx)*(F3(w1[j+2],w2[j+2],w3[j+2]) - F3(w1[j+1],w2[j+1],w3[j+1]));

      temp_1[j] = w1[j+1] - (dt/dx)*(F1(w_medios_1_f,w_medios_2_f,w_medios_3_f) - F1(w_medios_1_p,w_medios_2_p,w_medios_3_p));
      temp_2[j]= w2[j+1] - (dt/dx)*(F2(w_medios_1_f,w_medios_2_f,w_medios_3_f) - F2(w_medios_1_p,w_medios_2_p,w_medios_3_p));
      temp_3[j]= w3[j+1] - (dt/dx)*(F3(w_medios_1_f,w_medios_2_f,w_medios_3_f) - F3(w_medios_1_p,w_medios_2_p,w_medios_3_p));
    }
    for(int j = 0; j<(N-1); j++){
      w1[j+1] = temp_1[j];
      w2[j+1] = temp_2[j];
      w3[j+1] = temp_3[j];
    }
    n+=1;
  }
  free(temp_1);
  free(temp_2);
  free(temp_3);
}
