#include <iostream>
#include <cmath>
#include "constant.h"
using namespace std;

//  
void prolongation( int Nxyii, double *vxyii, double *vxyi, int nfactor){
  
     int Ntopii  = Nxyii - 1;
     int Ntopi   = nfactor*Ntopii;
     int Nxyi    = Ntopi + 1;
     int i, j, ixyi, ixyii, ixyupi,ixyupli, ixyupri, ixydowni, ixydownli, ixydownri, ixylefti, ixyrighti; 
     double  *vtemp     = new double [Nxyi*Nxyi];  

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


          vtemp[ixyi]      =     vtemp[ixyi] + vxyii[ixyii];    //    
 
          vtemp[ixyupi]    =     vtemp[ixyupi]      + 0.5*vxyii[ixyii];    //   
          vtemp[ixydowni]  =     vtemp[ixydowni]    + 0.5*vxyii[ixyii];    //    
          vtemp[ixylefti]  =     vtemp[ixylefti]    + 0.5*vxyii[ixyii];    //    
          vtemp[ixyrighti] =     vtemp[ixyrighti]   + 0.5*vxyii[ixyii];    //    

          vtemp[ixyupli]     =   vtemp[ixyupli]     + 0.25*vxyii[ixyii];  //  
          vtemp[ixyupri]     =   vtemp[ixyupri]     + 0.25*vxyii[ixyii];  //   
          vtemp[ixydownli]   =   vtemp[ixydownli]   + 0.25*vxyii[ixyii];  //    
          vtemp[ixydownri]   =   vtemp[ixydownri]   + 0.25*vxyii[ixyii];  //   

       }
     }



     for(i=1;i<(Nxyi-1);i++){
        for(j=1;j<(Nxyi-1);j++){ 
          ixyi            =  j*Nxyi  + i; 
          vxyi[ixyi]      =  vtemp[ixyi];     //         //     vxyi[ixyi] + vxyii[ixyii];    //    

         }
      }


   delete[] vtemp;

}

