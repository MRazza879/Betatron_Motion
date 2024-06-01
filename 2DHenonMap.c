//C environment for the simulation of a 2D Hénon map
//Author: Marta Razza
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//This code computes the tracking of particles in a 2D Hénon map and the Lyapunov error.
//The aim is to study the stability of particle orbits along a circular accelerator modelled by a Hénon map. 
//With the help of the tracking it is possible to catch the main features of the particle motion, 
//but as largely explained in the report LE  allows to obtain information about particles without the need of tracking them.
//More details about theory and formulas are given in the report.


//The modulation is due to unavoidable external effects 
//such as the coupling between synchrotron and betatron motion 
//(respectively motion on longitudinal and transversal plane) 
//or power supply ripple and can be modelled by a set of non-linear oscillations



//Definition of functions: 
//the first one computes the Hénon map, 
//the second one the tangent Hénon map.

void henon(double* x,double* p, double nu,int i, double eps);
void tangent(double* etax, double* etap,double x,double p, double nu,int i, double eps);

int main(int argc, char *argv[])
{

	if (argc != 5) {
		printf("Usage: %s nu eps N_steps N_iterations\n", argv[0]);
		return 1;
	}

	//input parameters
	double nu = atof(argv[1]); //the tune gives the oscillation frequency of betatronic motion
	double eps = atof(argv[2]); //the modulation 
	int N_steps = atoi(argv[3]); //number of steps to explore lattice
	int N_iterations = atoi(argv[4]); //number of iterations to compute LE and REM

	char filename [50];
	sprintf(filename, "tracking_nu%s_eps%s.txt", argv[1], argv[2]);
	FILE *track = fopen(filename, "w");
	if (track == NULL) {
    	perror("Error opening file");
    	return 1;
	}

	
	char filename2 [50];
	sprintf(filename2, "lyapunov_error2D_nu%s__eps%s_%s.txt", argv[1],argv[2],argv[4]);
	FILE *lyap = fopen(filename2, "w");
	if (lyap == NULL) {
		perror("Error opening file");
		return 1;
	}
	
	double x0,p0,x,p,eta1x,eta1p,eta2x,eta2p,ERLyap;
	int i,k,j;

	
	//Particle tracking
	//In this way we select 30 different initial conditions evolving them using Hénon Map. 
	//This allows to catch the main features of the phase-space

	for(i=1;i<=30;i++)
	{
		x0=(double)i*(1./50.);
		p0=x0;
		x=x0;
		p=p0;

		for(k=1;k<=1000;k++)
		{
			henon(&x,&p,nu,k,eps);
			fprintf(track,"%lf\t%lf\n",x,p);
		}
	}

	//End of particle tracking

	//Here the Lyapunov error is computed 
	//With the same tune and the same modulation used to extract results with the tracking, 
	//LE gives a more evident visualization of stable and unstable regions.
	//A comparison is shown in the report.

	for (i=(-1)*N_steps;i<N_steps;i++) //Cycle on x
	{
		for (j=(-1)*N_steps;j<N_steps;j++) //Cycle on p
		{
			x0=(double)i*(1.)/(double)N_steps;
			p0=(double)j*(1.)/(double)N_steps;
			eta1x=1;
			eta1p=0;
			eta2x=0;
			eta2p=1;
			x=x0;
			p=p0;
		
            for(k=1;k<=N_iterations;k++)
            {
            	henon(&x,&p,nu,k,eps);
            	tangent(&eta1x,&eta1p,x,p,nu,k,eps);
            	tangent(&eta2x,&eta2p,x,p,nu,k,eps); 
                ERLyap=(sqrt(eta1x*eta1x+eta1p*eta1p+eta2x*eta2x+eta2p*eta2p));
			}
		
        	fprintf(lyap,"%lf\t%lf\t%lf\n",x0,p0,ERLyap);

		}
	}
	fclose(track);
	fclose(lyap);

			
}

//this function computes the Hénon map evolving the initial conditions x0 and p0
void henon(double* x,double* p,double nu,int i,double eps)
{
	double omega;
	double xold,pold;
  	double omega1 = 2.*M_PI/868.12; 
	double tk[7]={omega1,2.*omega1,3.*omega1,6.*omega1,7.*omega1,10.*omega1,12.*omega1};    
  	double epsilon[7]={1.E-4,0.218E-4,0.708E-4,0.254E-4,0.100E-4,0.078E-4,0.218E-4};    
  	double t0=2.*3.14159265358*nu;
	omega=t0*(1.+eps*(epsilon[0]*cos(tk[0]*i)+epsilon[1]*cos(tk[1]*i)+epsilon[2]*cos(tk[2]*i)+epsilon[3]*cos(tk[3]*i)+epsilon[4]*cos(tk[4]*i)+epsilon[5]*cos(tk[5]*i)+epsilon[6]*cos(tk[6]*i)));
	xold=*x;
	pold=*p;
	*x=cos(omega)*xold+sin(omega)*(pold+xold*xold);
	*p=-sin(omega)*xold+cos(omega)*(pold+xold*xold);
}


//this function computes the tangent Hénon map that evolves the tangent vector (etax,etap) starting from the initial conditions x0 and p0 as given by the definition of Lyapunov error
void tangent(double* etax, double* etap,double x,double p, double nu, int i, double eps)
{
	double eta_x_old=*etax;
	double eta_p_old=*etap;
	double omega;
	double omega1 = 2.*M_PI/868.12; 
	double t0=2.*3.14159265358*nu;
	double tk[7]={omega1,2.*omega1,3.*omega1,6.*omega1,7.*omega1,10.*omega1,12.*omega1};   
  	double epsilon[7]={1.E-4,0.218E-4,0.708E-4,0.254E-4,0.100E-4,0.078E-4,0.218E-4};     
    omega=t0*(1.+eps*(epsilon[0]*cos(tk[0]*i)+epsilon[1]*cos(tk[1]*i)+epsilon[2]*cos(tk[2]*i)+epsilon[3]*cos(tk[3]*i)+epsilon[4]*cos(tk[4]*i)+epsilon[5]*cos(tk[5]*i)+epsilon[6]*cos(tk[6]*i)));
    *etax=(cos(omega)+2.*sin(omega)*x)*eta_x_old+sin(omega)*eta_p_old;
    *etap=(-sin(omega)+2.*cos(omega)*x)*eta_x_old+cos(omega)*eta_p_old;

}



