#include <iostream>
#include <cmath>
#include "constant.h"
#include "multigrid_header.h"
using namespace std;

void vcycle( int Ntop, double dxy, double *vxy, double *vxynew, double *fv, int nlevel){


     int i, j, k, Nxy, ixy, ntemp;
     int ixyup, ixydown, ixyleft, ixyright;
     int ntotal      =  (int) log2(Ntop);
    // int nlevel      =  5;  //  
     int ithlevel;  
     double L2error;   
     double tolfactor;


     // create array of pointers: the vxy for each level;
     double **vxymatrix        = new double* [nlevel]; 
     double **vxynewmatrix     = new double* [nlevel];
     double **fvmatrix         = new double* [nlevel];
     double **residualmatrix   = new double* [nlevel];
     double *vdxy              = new double  [nlevel];


     for(i=0;i<nlevel;i++){
         Nxy               = (int) pow(2.0,ntotal-i) + 1;
         vxymatrix[i]      = new double [Nxy*Nxy];
         vxynewmatrix[i]   = new double [Nxy*Nxy];
         fvmatrix[i]       = new double [Nxy*Nxy];
         residualmatrix[i] = new double [Nxy*Nxy];
         vdxy[i]           = pow(2.0,i)*dxy;
      }


      //***** pointing to vxy, vxynew, fv;
      vxymatrix[0]    = vxy;
      vxynewmatrix[0] = vxynew;
      fvmatrix[0]     = fv;


      //***** restriction;
      for(i=0;i<(nlevel-1);i++){ 
          Nxy      = (int) pow(2.0,ntotal-i) + 1;
          jacobsolver( Nxy, vdxy[i], vxymatrix[i], vxynewmatrix[i], fvmatrix[i], residualmatrix[i], 0.01);
          hto2hrestriction( Nxy, residualmatrix[i], fvmatrix[i+1], 2);       //       restriction; 

       }



      //***** solve for the last level;
      ntemp = nlevel-1;
      Nxy   = (int) pow(2.0,ntotal-ntemp) + 1;
      jacobsolver( Nxy, vdxy[ntemp], vxymatrix[ntemp], vxynewmatrix[ntemp], fvmatrix[ntemp], residualmatrix[ntemp], 0.00000000000001);  

      //***** prolongation;
      for(i=ntemp;i>0;i--){ 
          Nxy      = (int) pow(2.0,ntotal-i) + 1;
          h2tohprolongation( Nxy, vxymatrix[i], vxymatrix[i-1], 2 );         //       prolongation; 
      }


     //***** MAKE SURE: the dynamically allocated memory needs to be free, OTHERWISE, you will see the computer memory usage is increasing, AND,
     //***** eventually, the code running is slowed down;
     for(i=1;i<nlevel;i++){
         delete[] vxymatrix[i];
         delete[] vxynewmatrix[i]; 
         delete[] fvmatrix[i];
         delete[] residualmatrix[i];  
      }
     delete[] residualmatrix[0];
     delete[] vdxy;  
     delete[] vxymatrix;
     delete[] vxynewmatrix;
     delete[] fvmatrix;
     delete[] residualmatrix;

}


















