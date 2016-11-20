# include <iostream>
# include <cmath>
# include <new>
# include <iomanip>
using namespace std;

void jacobsolverexact( int Nxy, double dxy, double *vxy, double *vxynew, double *fv, double *vresidual, double tolfactor){

  int     i, j, k;
  int     ixy;
  int     ixyup, ixydown, ixyleft, ixyright;
  double  *temp;
  double  residual = 1.0, residual2 = 1.0, residualinit = 0.0;
  double  w = 2.0/3.0;  
  double  fcoef = 1.0/(dxy*dxy);  

  k = 0;
while(residual2>tolfactor){
//  while(k<20){
    k = k + 1;
    residual2 = 0.0;

    for(i=1;i<(Nxy-1);i++){
      for(j=1;j<(Nxy-1);j++){
         ixy          =  j*Nxy + i;
         ixyup        =  ixy + Nxy;  
         ixydown      =  ixy - Nxy;   
         ixyleft      =  ixy - 1;     
         ixyright     =  ixy + 1;     
         residual     =  0.25*( vxy[ixyup]+vxy[ixydown]+vxy[ixyleft]+vxy[ixyright]+dxy*dxy*fv[ixy] ) - vxy[ixy];    //     
         vxynew[ixy]  =  vxy[ixy]  + w*residual;
         residual2    =  residual2 + residual*residual;
      }
    }

    residual2 = sqrt(residual2/((Nxy-1)*(Nxy-1)));   

    // get the initial residual;
    if(k==1)
       residualinit = residual2;

    temp      = vxy;
    vxy       = vxynew;
    vxynew    = temp;   

  } // while;


     //  to compute residual;
    for(i=1;i<(Nxy-1);i++){
      for(j=1;j<(Nxy-1);j++){
         ixy            =  j*Nxy + i;
         ixyup          =  ixy + Nxy;  
         ixydown        =  ixy - Nxy;   
         ixyleft        =  ixy - 1;     
         ixyright       =  ixy + 1;     
         residual       =  -fcoef*( vxy[ixyup]+vxy[ixydown]+vxy[ixyleft]+vxy[ixyright] - 4.0*vxy[ixy] ) - fv[ixy];         
         vresidual[ixy] =  -residual;      
      }
    }


}


