/*-------------------------------------------------------------*/
/*  CLB - Cudne LB - Stencil Version                           */
/*     CUDA based Adjoint Lattice Boltzmann Solver             */
/*     Author: Lukasz Laniewski-Wollk                          */
/*     Developed at: Warsaw University of Technology - 2012    */
/*-------------------------------------------------------------*/

/*
Model created by Travis Mitchell 10-03-2016. Purpose of model is
to give an introduction into model creation in TCLB, tutorial 
file is available upon request.

Model solves d2q9 files via applying the single relaxation time
BGK-lattice Boltzmann method
*/

#include <math.h>
#define PI 3.141592653

CudaDeviceFunction float2 Color() {
// used for graphics - can usually ignore function
 /*       float2 ret;
        vector_t u = getU();
        ret.x = sqrt(u.x*u.x + u.y*u.y);
        if (NodeType == NODE_Solid){
                ret.y = 0;
        } else {
                ret.y = 1;
        }
        return ret;*/
}

CudaDeviceFunction void Init(){
        real_t u[2] = {Velocity_x, Velocity_y};
	PressureF = 0.0; // To be updated in next stage.
        real_t d    = Density;

	if (Rad > 0) {
	real_t drop = Rad -  sqrt((X-cX)*(X-cX) + (Y-cY)*(Y-cY))- Rad*cos(45.0*PI/180.0);

       	if( drop > dLow ) {
		d  = drop;
	} else {
		d  = dLow;
	}
	}
/*       	if (Rad > 0){
		real_t Ri = sqrt( (X - cX)*(X - cX) + (Y - cY)*(Y - cY) );
		d = 0.5 * (dHigh + dLow)
		  - 0.5 * (dHigh - dLow) * tanh(2.0*(Ri - Rad)/5);
	}*/

	SetEquilibrium(d,u);

	DensityF = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8];
	// if (X==256 && Y==256) printf("Init: Density=%.4f, Pressure=%.4f \n", DensityF, PressureF);
}

CudaDeviceFunction void calcDensityF() {
	DensityF = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8];
}

CudaDeviceFunction void calcPressureF() {
DensityF = DensityF(0,0);

real_t dhdx =  ( DensityF(1,0) - DensityF(-1,0) )/2.0;
real_t dhdy =  ( DensityF(0,1) - DensityF(0,-1) )/2.0;
real_t d2h  =  ( DensityF(1,0) - 2.0 * DensityF(0,0) + DensityF(-1,0) )/((1.0 + dhdx*dhdx)*sqrt(1.0 + dhdx*dhdx))
            +  ( DensityF(0,1) - 2.0 * DensityF(0,0) + DensityF(0,-1) )/((1.0 + dhdy*dhdy)*sqrt(1.0 + dhdy*dhdy));

real_t kappa = (n-1.0)*(m-1.0)/((n-m)* hstar);

real_t Angle = ContactAngle * PI/180.0;

PressureF   =  - SurfaceTension * d2h - SurfaceTension * Angle * kappa * ( pow( (hstar/(DensityF(0,0)+hc)), n)
	       - pow( (hstar/(DensityF(0,0)+hc)), m ) ) ;

	// if (X==256 && Y==256) printf("calcP: Density=%.4f, Pressure=%.4f \n", DensityF, PressureF);

}

CudaDeviceFunction void SetEquilibrium(real_t d, real_t u[2])
{
real_t g = Gravitation;
real_t t1 = u[0]*u[0] + u[1]*u[1];

f[0] = d - (5.0/6.0)*g*d*d - (2.0/3.0)*d*t1;
f[1] =  wf[1]*d*(1.5*g*d + 3.0*(d2q9_ex[1]*u[0] + d2q9_ey[1]*u[1])+ 4.5*(d2q9_ex[1]*u[0] + d2q9_ey[1]*u[1])*(d2q9_ex[1]*u[0] + d2q9_ey[1]*u[1]) - 1.5*t1 );
f[2] =  wf[2]*d*(1.5*g*d + 3.0*(d2q9_ex[2]*u[0] + d2q9_ey[2]*u[1])+ 4.5*(d2q9_ex[2]*u[0] + d2q9_ey[2]*u[1])*(d2q9_ex[2]*u[0] + d2q9_ey[2]*u[1]) - 1.5*t1 );
f[3] =  wf[3]*d*(1.5*g*d + 3.0*(d2q9_ex[3]*u[0] + d2q9_ey[3]*u[1])+ 4.5*(d2q9_ex[3]*u[0] + d2q9_ey[3]*u[1])*(d2q9_ex[3]*u[0] + d2q9_ey[3]*u[1]) - 1.5*t1 );
f[4] =  wf[4]*d*(1.5*g*d + 3.0*(d2q9_ex[4]*u[0] + d2q9_ey[4]*u[1])+ 4.5*(d2q9_ex[4]*u[0] + d2q9_ey[4]*u[1])*(d2q9_ex[4]*u[0] + d2q9_ey[4]*u[1]) - 1.5*t1 );
f[5] =  wf[5]*d*(1.5*g*d + 3.0*(d2q9_ex[5]*u[0] + d2q9_ey[5]*u[1])+ 4.5*(d2q9_ex[5]*u[0] + d2q9_ey[5]*u[1])*(d2q9_ex[5]*u[0] + d2q9_ey[5]*u[1]) - 1.5*t1 );
f[6] =  wf[6]*d*(1.5*g*d + 3.0*(d2q9_ex[6]*u[0] + d2q9_ey[6]*u[1])+ 4.5*(d2q9_ex[6]*u[0] + d2q9_ey[6]*u[1])*(d2q9_ex[6]*u[0] + d2q9_ey[6]*u[1]) - 1.5*t1 );
f[7] =  wf[7]*d*(1.5*g*d + 3.0*(d2q9_ex[7]*u[0] + d2q9_ey[7]*u[1])+ 4.5*(d2q9_ex[7]*u[0] + d2q9_ey[7]*u[1])*(d2q9_ex[7]*u[0] + d2q9_ey[7]*u[1]) - 1.5*t1 );
f[8] =  wf[8]*d*(1.5*g*d + 3.0*(d2q9_ex[8]*u[0] + d2q9_ey[8]*u[1])+ 4.5*(d2q9_ex[8]*u[0] + d2q9_ey[8]*u[1])*(d2q9_ex[8]*u[0] + d2q9_ey[8]*u[1]) - 1.5*t1 );
}

CudaDeviceFunction void Run() {
	DensityF = DensityF(0,0);
	AddToMass(DensityF);
	PressureF= PressureF(0,0);
	// if (X==256 && Y==256) printf("Run: Density=%.4f, Pressure=%.4f \n", DensityF, PressureF);

// This defines the dynamics that we run at each node in the domain.
    switch (NodeType & NODE_BOUNDARY) {
	case NODE_Solid:
	case NODE_Wall:
		BounceBack();
		break;
    }
	if (NodeType & NODE_BGK)
	{
		CollisionBGK();
	}
}

CudaDeviceFunction void CollisionBGK() {
	// if (X==256 && Y==256) printf("Collision1: Density=%.4f, Pressure=%.4f \n", DensityF, PressureF);
	// if (X==256 && Y==256) printf("Collision1: f0=%.4f, f1=%.4f, f2=%.4f \n", f[0], f[1], f[2]);

real_t f_temp[9];
for (int i=0; i<9; i++) {
	f_temp[i] = f[i];
	}

real_t d = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8];
real_t u[2];

u[0]  = (f[1]-f[3]+f[5]-f[6]-f[7]+f[8])/d;
u[1]  = (f[2]-f[4]+f[5]+f[6]-f[7]-f[8])/d;

SetEquilibrium(d,u);
	// if (X==256 && Y==256) printf("Collision2: Density=%.4f, Pressure=%.4f \n", d, PressureF);
	// if (X==256 && Y==256) printf("Collision2: f0=%.4f, f1=%.4f, f2=%.4f \n", f[0], f[1], f[2]);

DensityF = DensityF(0,0);
PressureF= PressureF(0,0);

/***************************** Centered scheme proposed by Zhou ******************************/
real_t s_p1p1    =   0.25 * ( DensityF(1,1) + DensityF(1,0) + DensityF(0,1) + DensityF(0,0)  );
real_t s_m1p1    =   0.25 * ( DensityF(-1,1) + DensityF(-1,0) + DensityF(0,0)  + DensityF(0,1) );
real_t s_m1m1    =   0.25 * ( DensityF(-1,-1) + DensityF(-1,0) + DensityF(0,0)  + DensityF(0,-1) );
real_t s_p1m1    =   0.25 * ( DensityF(0,-1) + DensityF(1,-1) + DensityF(1,0) + DensityF(0,0)  );
/**************************** Pressure forcing based on the centered scheme ******************************/
/********************* For the first 4 velocities it is just a forward difference ****************/
real_t Fx_p1p0   = - 0.5 * ( DensityF(1,0) + DensityF(0,0) ) * ( PressureF(1,0) - PressureF(0,0)  );
real_t Fx_m1m0   = - 0.5 * ( DensityF(-1,0) + DensityF(0,0) ) * ( PressureF(0,0)  - PressureF(-1,0) );
real_t Fy_p0p1   = - 0.5 * ( DensityF(0,1) + DensityF(0,0) ) * ( PressureF(0,1) - PressureF(0,0)  );
real_t Fy_m0m1   = - 0.5 * ( DensityF(0,-1) + DensityF(0,0) ) * ( PressureF(0,0)  - PressureF(0,-1) );

/* For velocities 5-8 it is the average of the contained two forward differences */
/* Keep in mind that for x these are different differences then for y */
real_t Fx_p1p1   = - s_p1p1 * ( PressureF(1,1) + PressureF(1,0) - PressureF(0,1) - PressureF(0,0)    ) * 0.5 ;
real_t Fy_p1p1   = - s_p1p1 * ( PressureF(1,1) - PressureF(1,0) + PressureF(0,1) - PressureF(0,0)    ) * 0.5 ;

real_t Fx_m1p1   = - s_m1p1 * ( PressureF(0,1)  + PressureF(0,0)  - PressureF(-1,0) - PressureF(-1,1)  ) * 0.5 ;
real_t Fy_m1p1   = - s_m1p1 * ( PressureF(0,1)  - PressureF(0,0)  - PressureF(-1,0) + PressureF(-1,1)  ) * 0.5 ;

real_t Fx_m1m1   = - s_m1m1 * ( PressureF(0,0)   + PressureF(0,-1)  - PressureF(-1,0) - PressureF(-1,-1) ) * 0.5 ;
real_t Fy_m1m1   = - s_m1m1 * ( PressureF(0,0)   - PressureF(0,-1)  + PressureF(-1,0) - PressureF(-1,-1) ) * 0.5 ;

real_t Fx_p1m1   = - s_p1m1 * ( PressureF(1,0) + PressureF(1,-1)  - PressureF(0,-1) - PressureF(0,0)   ) * 0.5 ;
real_t Fy_p1m1   = - s_p1m1 * ( PressureF(1,0) - PressureF(1,-1)  - PressureF(0,-1) + PressureF(0,0)   ) * 0.5 ;

real_t shear_1, shear_2, shear_3, shear_4;
if (Substrate == 1){
	shear_1 = -3*viscosity*u[0]/(DensityF+3*SlipLength+1.5*SlipLength*SlipLength/DensityF);
	shear_2 = -3*viscosity*u[1]/(DensityF+3*SlipLength+1.5*SlipLength*SlipLength/DensityF);
	shear_3 = 0.0;
	shear_4 = 0.0;
} else {
	shear_1 = 0;
	shear_2 = 0;
}
/*if( rho[i][j] > 0.052 && l > 50000 )
		f_hoz 	  = 0.0001;
else
		f_hoz	  = 0;
*/
/*if(l<=forceme){
	  dT_1 = ( temp[idx(ip,j)] - temp[idx(i,j)] );
	  dT_2 = ( temp[idx(i,jp)] - temp[idx(i,j)] );
	  dT_3 = ( temp[idx(i,j)] - temp[idx(in,j)] );
	  dT_4 = ( temp[idx(i,j)] - temp[idx(i,jn)] );
}else{
	  dT_1 = 0; dT_2 = 0; dT_3 = 0; dT_4 = 0;
}*/
/********************************* Shear stress at the bottem ***********************************/
/*fx_1	  = - g * (z_bed[idx(ip,j)] - z_bed[idx(i,j)])  * (rho[idx(i,j)] + rho[idx(ip,j)]) * 0.5  + shear_1 + f_hoz;
fx_0	  = - g * (z_bed[idx(i,j)] - z_bed[idx(in,j)])  * (rho[idx(i,j)] + rho[idx(in,j)]) * 0.5  + shear_1 + f_hoz;
fy_1	  = - g * (z_bed[idx(i,jp)] - z_bed[idx(i,j)])  * (rho[idx(i,j)] + rho[idx(i,jp)]) * 0.5  + shear_2;
fy_0	  = - g * (z_bed[idx(i,j)] - z_bed[idx(i,jn)])  * (rho[idx(i,j)] + rho[idx(i,jn)]) * 0.5  + shear_2;
*/
/************************************** Collision ***********************************************/
f[0]    =  omega * f[0] + (1 - omega)*f_temp[0];
f[1]    =  omega * f[1] + (1 - omega)*f_temp[1] + (1.0/6.0) * ( wf[1] * (shear_1 + Fx_p1p0));	// + 1.5 * dT_1));
f[2]    =  omega * f[2] + (1 - omega)*f_temp[2] + (1.0/6.0) * ( wf[2] * (shear_2 + Fy_p0p1));	// + 1.5 * dT_2));
f[3]    =  omega * f[3] + (1 - omega)*f_temp[3] + (1.0/6.0) * ( wf[3] * (-shear_1 - Fx_m1m0));	// - 1.5 * dT_3));
f[4]    =  omega * f[4] + (1 - omega)*f_temp[4] + (1.0/6.0) * ( wf[4] * (-shear_2 - Fy_m0m1));	// - 1.5 * dT_4));
f[5]    =  omega * f[5] + (1 - omega)*f_temp[5] + (1.0/6.0) * ( wf[5] * (shear_1 + shear_2 + Fx_p1p1 + Fy_p1p1));	// + 1.5*(dT_1+dT_2)));
f[6]    =  omega * f[6] + (1 - omega)*f_temp[6] + (1.0/6.0) * ( wf[6] * (-shear_1 + shear_2 - Fx_m1p1 + Fy_m1p1));	// + 1.5*(-dT_3+dT_2)));
f[7]    =  omega * f[7] + (1 - omega)*f_temp[7] + (1.0/6.0) * ( wf[7] * (-shear_1 - shear_2 - Fx_m1m1 - Fy_m1m1));	// + 1.5*(-dT_3-dT_4)));
f[8]    =  omega * f[8] + (1 - omega)*f_temp[8] + (1.0/6.0) * ( wf[8] * (shear_1 - shear_2 + Fx_p1m1 - Fy_p1m1));	// + 1.5*(dT_1-dT_4)));

	// if (X==256 && Y==256) printf("Collision*: Density=%.4f \n",f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8] );
	// if (X==256 && Y==256) printf("Collision*: f0=%.4f, f1=%.4f, f2=%.4f \n", f[0], f[1], f[2]);


}



CudaDeviceFunction void BounceBack() {
// Method to reverse distribution functions along the bounding nodes.
     	real_t uf;
	uf = f[3];
	f[3] = f[1];
	f[1] = uf;
	uf = f[4];
	f[4] = f[2];
	f[2] = uf;
	uf = f[7];
	f[7] = f[5];
	f[5] = uf;
	uf = f[8];
	f[8] = f[6];
	f[6] = uf;
}


CudaDeviceFunction real_t getRho() {
// This function defines the macroscopic density at the current node.
	DensityF = DensityF(0,0);
	PressureF= PressureF(0,0);
	// if (X==256 && Y==256) printf("getRho Density=%.4f, Pressure=%.4f \n", DensityF, PressureF);
	return DensityF;
}

CudaDeviceFunction real_t getP(){
	DensityF = DensityF(0,0);
	PressureF= PressureF(0,0);
	return PressureF;
}

CudaDeviceFunction vector_t getU() {
// This function defines the macroscopic velocity at the current node.
	DensityF = DensityF(0,0);
	PressureF= PressureF(0,0);
	// if (X==256 && Y==256) printf("getU: Density=%.4f, Pressure=%.4f \n", DensityF, PressureF);
	vector_t u;
	u.x = (( f[8]-f[7]-f[6]+f[5]-f[3]+f[1] )/DensityF );
	u.y = ((-f[8]-f[7]+f[6]+f[5]-f[4]+f[2] )/DensityF );
	u.z = 0;
	return u;
}
