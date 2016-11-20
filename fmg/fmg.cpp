#include <iostream>
#include <cmath>
#include "constant.h"
#include "multigrid_header.h"
using namespace std;

// fmg( N, dxy, vxy, fv, nlevel); 
void fmg( int Ntop, double dxytop, double *vxy, double *fv, int nlevel){

     int i, j, k, ntemp, counter, ixy, ixyup, ixydown, ixyleft, ixyright; 
     int Nxy;
     double dxy, L2residual;
     int ntotal = (int) log2(Ntop);

     int      *vNxy        = new int     [nlevel];
     double   *vdxy        = new double  [nlevel];
     double **fvmatrix     = new double* [nlevel];
     double **vxymatrix    = new double* [nlevel];
     double **vxynewmatrix = new double* [nlevel];
     double **vresmatrix   = new double* [nlevel];
 

     for(i=0;i<nlevel;i++){
        vNxy[i]         = (int) pow(2.0,ntotal-i) + 1;
        vdxy[i]         = pow(2.0,i)*dxytop;               
        fvmatrix[i]     = new double [vNxy[i]*vNxy[i]];
        vxymatrix[i]    = new double [vNxy[i]*vNxy[i]];
        vxynewmatrix[i] = new double [vNxy[i]*vNxy[i]];
        vresmatrix[i]   = new double [vNxy[i]*vNxy[i]];
     }         

     fvmatrix[0]  = fv;
     vxymatrix[0] = vxy;   

    // restrict f vectors;
    for(i=1;i<nlevel;i++){      
       restriction(vNxy[i-1], fvmatrix[i-1], fvmatrix[i], 2);
    }


   for(i=(nlevel-1);i>0;i--){
     
      Nxy = vNxy[i];
      dxy = vdxy[i]; 
      if(i==(nlevel-1)){
         jacobsolverexact( Nxy, dxy, vxymatrix[i], vxynewmatrix[i], fvmatrix[i], vresmatrix[i], 0.00000000001);
         prolongation( Nxy, vxymatrix[i], vxymatrix[i-1], 2);
       }
      else{
          vcycle( Nxy-1, dxy, vxymatrix[i], vxynewmatrix[i], fvmatrix[i], nlevel-i);
          prolongation( Nxy, vxymatrix[i], vxymatrix[i-1], 2);
       }

    } 

    // the finest grid level; 
    vcycle( vNxy[0]-1, vdxy[0], vxymatrix[0], vxynewmatrix[0], fvmatrix[0], nlevel );


}


