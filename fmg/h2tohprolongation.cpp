#include <iostream>
#include <cmath>
#include "constant.h"
using namespace std;

//  h2tohprolongation( Nxy, vxymatrix[i], vxymatrix[i-1] );
//  h2tohprolongation( N/2+1, fv2, fv, 2);    //   for test;
void h2tohprolongation( int Nxyii, double *vxyii, double *vxyi, int nfactor){
  
     int Ntopii  = Nxyii - 1;
     int Ntopi   = nfactor*Ntopii;
     int Nxyi    = Ntopi + 1;
     int i, j, ixyi, ixyii, ixyupi,ixyupli, ixyupri, ixydowni, ixydownli, ixydownri, ixylefti, ixyrighti; 

     for(i=1;i<(Nxyii-1);i++){
        for(j=1;j<(Nxyii-1);j++){

          ixyii         =  j*Nxyii + i;   
          ixyi          =  j*nfactor*Nxyi  + i*nfactor; 
     
          ixyupi        =  ixyi   + Nxyi;
          ixyupli       =  ixyupi - 1;
          ixyupri       =  ixyupi + 1;

          ixydowni      =  ixyi     - Nxyi;
          ixydownli     =  ixydowni - 1;
          ixydownri     =  ixydowni + 1;

          ixylefti      =  ixyi - 1;
          ixyrighti     =  ixyi + 1; 
          
          vxyi[ixyi]      =  vxyi[ixyi] + vxyii[ixyii];

          vxyi[ixyupi]    =  vxyi[ixyupi]    + 0.5*vxyii[ixyii];
          vxyi[ixydowni]  =  vxyi[ixydowni]  + 0.5*vxyii[ixyii];
          vxyi[ixylefti]  =  vxyi[ixylefti]  + 0.5*vxyii[ixyii];
          vxyi[ixyrighti] =  vxyi[ixyrighti] + 0.5*vxyii[ixyii];

          vxyi[ixyupli]     =  vxyi[ixyupli]     + 0.25*vxyii[ixyii];
          vxyi[ixyupri]     =  vxyi[ixyupri]     + 0.25*vxyii[ixyii];
          vxyi[ixydownli]   =  vxyi[ixydownli]   + 0.25*vxyii[ixyii];
          vxyi[ixydownri]   =  vxyi[ixydownri]   + 0.25*vxyii[ixyii];

        }
     }


}



/*


          
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
          
          vxyi[ixyi]      =  vxyi[ixyii] + vxyii[ixyii];

          vxyi[ixyupi]    =  vxyi[ixyupi]    + 0.5*vxyii[ixyii];
          vxyi[ixydowni]  =  vxyi[ixydowni]  + 0.5*vxyii[ixyii];
          vxyi[ixylefti]  =  vxyi[ixylefti]  + 0.5*vxyii[ixyii];
          vxyi[ixyrighti] =  vxyi[ixyrighti] + 0.5*vxyii[ixyii];

          vxyi[ixyupli]     =  vxyi[ixyupli]     + 0.25*vxyii[ixyii];
          vxyi[ixyupri]     =  vxyi[ixyupri]     + 0.25*vxyii[ixyii];
          vxyi[ixydownli]   =  vxyi[ixydownli]   + 0.25*vxyii[ixyii];
          vxyi[ixydownri]   =  vxyi[ixydownri]   + 0.25*vxyii[ixyii];


*/
