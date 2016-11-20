
void initialsolver(int Nxy, double *xy, double *vxy, double *vxyexact,double *fv);
void vcycle( int Ntop, double dxy, double *vxy, double *vxynew, double *fv, int nfactor);
void jacobsolver( int Nxy, double dx, double *vxy, double *vxynew, double *fv, double *vresidual, double tolfactor);
void hto2hrestriction( int Nxy, double *residuali, double *fvii, int nfactor);
void h2tohprolongation( int Nxyii, double *vxyii, double *vxyi, int nfactor);
void newtonsolver( int Nxy, double dxy, double *vxy, double *vxynew, double *fv, double *vresidual, double tolfactor);
void fmg( int Ntop, double dxy, double *vxy, double *fv, int nlevel);
void jacobsolverexact( int Nxy, double dxy, double *vxy, double *vxynew, double *fv, double *vresidual, double tolfactor);

void restriction( int Nxyi, double *vri, double *fvii, int nfactor);
void prolongation( int Nxyii, double *vxyii, double *vxyi, int nfactor);
