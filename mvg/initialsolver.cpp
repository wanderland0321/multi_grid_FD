#include <cmath>
#include "constant.h"

void initialsolver(int Nxy, double *xy, double *vxy, double *vxyexact,double *fv){

    int i, j, ixy;
    for(i=0;i<Nxy;i++){
      for(j=0;j<Nxy;j++){

        ixy           = j*Nxy+i;
        vxy[ixy]      = 0.0;
        vxyexact[ixy] = sin(10*PI*xy[i])*sin(20*PI*xy[j]);   
        fv[ixy]       = 100*PI*PI*sin(10*PI*xy[i])*sin(20*PI*xy[j]) + 400*PI*PI*sin(10*PI*xy[i])*sin(20*PI*xy[j]);    

        if(i==0){
            vxy[ixy]      = 0;
            vxyexact[ixy] = 0;
            fv[ixy]       = 0;
        }

        if(i==(Nxy-1)){
            vxy[ixy]      = 0;
            vxyexact[ixy] = 0;
            fv[ixy]       = 0;
        }

        if(j==0){
            vxy[ixy]      = 0;
            vxyexact[ixy] = 0;
            fv[ixy]       = 0;
        }

        if(j==(Nxy-1)){
            vxy[ixy]      = 0;          
            vxyexact[ixy] = 0;    
            fv[ixy]       = 0;    
        }
                                                                       
      } 
    }


}
