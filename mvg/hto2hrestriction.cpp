#include <iostream>
#include <cmath>
#include "constant.h"
using namespace std;

// hto2hrestriction( Nxy, residualmatrix[i], fvmatrix[i+1], 2);
void hto2hrestriction( int Nxyi, double *vri, double *fvii, int nfactor){

     int Ntopii = (Nxyi-1)/2;
     int Nxyii  = Ntopii+1;
     int i, j, ixyi, ixyii, ixyupi,ixyupli, ixyupri, ixydowni, ixydownli, ixydownri, ixylefti, ixyrighti;

     for(i=1;i<(Nxyii-1);i++){
       for(j=1;j<(Nxyii-1);j++){

          ixyii         =  j*Nxyii + i;   
          ixyi          =  j*nfactor*Nxyi  + i*nfactor;      
          ixyupi        =  ixyi   + Nxyi;
          ixyupli       =  ixyupi - 1;
          ixyupri       =  ixyupi + 1;
          ixydowni      =  ixyi   - Nxyi;
          ixydownli     =  ixydowni - 1;
          ixydownri     =  ixydowni + 1;
          ixylefti      =  ixyi - 1;
          ixyrighti     =  ixyi + 1; 
          fvii[ixyii]   =  0.25*vri[ixyi] + 0.0625*(vri[ixyupli]+vri[ixyupri]+vri[ixydownli]+vri[ixydownri])+0.125*(vri[ixyupi]+vri[ixydowni]+vri[ixylefti]+vri[ixyrighti]);     
          
       }
     }
     

}


