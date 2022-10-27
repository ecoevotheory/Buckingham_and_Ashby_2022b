/***********************************************************************************************************
 * [t,S,I,R,EQFLAG] = Simulation_earlybTradeoff_function(t_max,a,beta,b0,aR,bR,c1b,c2b,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total)
 ***********************************************************************************************************/

/* Compile in Matlab using mex Simulation_earlybTradeoff_function.c */

#include <mex.h>
#include <math.h>


/***********************************
 * Constant parameter values
 ***********************************/
#define MAXSTEPS 1e6 /* Maximum number of steps for ODE solver */
#define INTERVAL 1e2 /* Check if the system is close to equilibrium */
#define EPS 1e-6 /* ODE solver tolerance */
#define TINY 1e-30 /* Constant value for solver */
#define TINY2 1e-30 /* Constant value for solver to avoid tolerance issues */
/* RK solver parameters */
#define b21 0.2
#define b31 3.0/40.0
#define b32 9.0/40.0
#define b41 0.3
#define b42 -0.9
#define b43 1.2
#define b51 -11.0/54.0
#define b52 2.5
#define b53 -70.0/27.0
#define b54 35.0/27.0
#define b61 1631.0/55296
#define b62 175.0/512.0
#define b63 575.0/13824.0
#define b64 44275.0/110592
#define b65 253.0/4096.0
#define c1 37.0/378.0
#define c3 250.0/621.0
#define c4 125.0/594.0
#define c6 512.0/1771.0
#define dc5 -277.00/14336

/*************************************
 * Define structure for model parameters
 *************************************/
struct PARAM{
    double t_max;
    double a;
    double beta;
    double b0;
    double aR;
    double bR;
    double c1b;
    double c2b;
    double q;
    double f;
    double alpha;
    double gamma;
    double eqtol;
    int strain_total;
};

/*************************************
 * Function prototypes
 *************************************/
int my_rungkut (double *T, double *S_out, double *I_out, double *R_out, double *EQFLAG, double *init_pop, double *zeta, struct PARAM *p);
void rkqs(double *S, double *I, double *R, double *DSDT, double *DIDT, double *DRDT,  double *oldh, double *hnext, double *S_SCALE, double *I_SCALE, double *R_SCALE, double *zeta, struct PARAM *p);
void rkck(double *S, double *I,  double *R, double *DSDT, double *DIDT,  double *DRDT, double *Sout, double *Iout,  double *Rout, double *Serr, double *Ierr,  double *Rerr, double oldh, double *zeta, struct PARAM *p);
void dynamic(double *S, double *I, double *R, double *DSDT, double *DIDT, double *DRDT, double *zeta, struct PARAM *p);
double FMAX(double, double);
double FMIN(double, double);

/*************************************
 * Main function
 *************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *T, *S, *I, *R, *EQFLAG, *init_pop, *zeta, *parameter;
    double *t_temp, *S_temp, *I_temp, *R_temp;
    int i, j, k, colLen, maxsteps;
    struct PARAM p;


   /* Allocate inputs */
    if(nrhs!=16){
       mexErrMsgTxt("Incorrect number of input arguments!\n");
    }

    else{
        parameter= mxGetPr(prhs[0]);
        p.t_max= *parameter;
        parameter= mxGetPr(prhs[1]);
        p.a= *parameter;
        parameter= mxGetPr(prhs[2]);
        p.beta= *parameter;
        parameter= mxGetPr(prhs[3]);
        p.b0= *parameter;
	parameter= mxGetPr(prhs[4]);
        p.aR= *parameter;
	parameter= mxGetPr(prhs[5]);
        p.bR= *parameter;
	parameter= mxGetPr(prhs[6]);
        p.c1b= *parameter;
	parameter= mxGetPr(prhs[7]);
        p.c2b= *parameter;
        zeta= mxGetPr(prhs[8]);
        parameter= mxGetPr(prhs[9]);
        p.q= *parameter;
        parameter= mxGetPr(prhs[10]);
        p.f= *parameter;
        parameter= mxGetPr(prhs[11]);
        p.alpha= *parameter;
        parameter= mxGetPr(prhs[12]);
        p.gamma= *parameter;
        parameter= mxGetPr(prhs[13]);
        p.eqtol= *parameter;
        init_pop= mxGetPr(prhs[14]);    
        parameter= mxGetPr(prhs[15]);
        p.strain_total= (int)*parameter;
    }

    maxsteps = (int)MAXSTEPS;
    
    /* Allocate memory */
    t_temp = malloc(maxsteps*sizeof(double));
    S_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    I_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    R_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    
    /* Initialise this output */           
    plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
    EQFLAG = mxGetPr(plhs[4]);
    EQFLAG[0] = 0;

    /* Call ODE solver */
    colLen = my_rungkut(t_temp, S_temp, I_temp, R_temp, EQFLAG, init_pop, zeta, &p); 

    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    T = mxGetPr(plhs[0]);
    S = mxGetPr(plhs[1]);
    I = mxGetPr(plhs[2]);
    R = mxGetPr(plhs[3]);

    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        T[i] = t_temp[i];
        for (j=0;j<p.strain_total;j++) {
	    S[i+j*colLen]=S_temp[i+j*maxsteps];
	    I[i+j*colLen]=I_temp[i+j*maxsteps];
	    R[i+j*colLen]=R_temp[i+j*maxsteps];
        }
    }
    
    /* Free memory */
    free(t_temp);
    free(S_temp);
    free(I_temp); 
    free(R_temp);
    
    return;

}

/*****************************************
 * ODE solver
 ****************************************/
int my_rungkut (double *T, double *S_out, double *I_out, double *R_out, double *EQFLAG, double *init_pop, double *zeta, struct PARAM *p){
    
    double *S, *I,  *R, *DSDT, *DIDT, *DRDT, *S_SCALE, *I_SCALE, *R_SCALE;
    double *SMIN, *SMAX, *IMIN, *IMAX, *RMIN, *RMAX, hnext[1], oldh[1];
    double t, nextcheck;
    int i, j, exitflag, count, maxsteps;
    
    /* Allocate memory */
    S = malloc(p->strain_total*sizeof(double));
    I = malloc(p->strain_total*sizeof(double));
    R = malloc(p->strain_total*sizeof(double));
    DSDT = malloc(p->strain_total*sizeof(double));
    DIDT = malloc(p->strain_total*sizeof(double));
    DRDT = malloc(p->strain_total*sizeof(double));
    S_SCALE = malloc(p->strain_total*sizeof(double));
    I_SCALE = malloc(p->strain_total*sizeof(double));
    R_SCALE = malloc(p->strain_total*sizeof(double));
    SMIN = malloc(p->strain_total*sizeof(double)); 
    SMAX = malloc(p->strain_total*sizeof(double)); 
    IMIN = malloc(p->strain_total*sizeof(double)); 
    IMAX = malloc(p->strain_total*sizeof(double)); 
    RMIN = malloc(p->strain_total*sizeof(double)); 
    RMAX = malloc(p->strain_total*sizeof(double)); 
    
    /* Other parameters */
    exitflag = 1;
    count=0;
    /* k=1;*/
    oldh[0] = 1e-3;
    hnext[0] = 1e-3;
    t=0;
    nextcheck = INTERVAL;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations */
    for(i=0;i<p->strain_total;i++){
        S[i] = init_pop[3*i];
        I[i] = init_pop[3*i+1];
        R[i] = init_pop[3*i+2];
    }
    
    /* Initialise equilibrium arrays */
    for(i=0;i<p->strain_total;i++){
        SMIN[i] = S[i];
        SMAX[i] = S[i];
        IMIN[i] = I[i];
        IMAX[i] = I[i];
        RMIN[i] = R[i];
        RMAX[i] = R[i];
    }
    
    /* Update output */
    T[0]=t;
    for (i=0; i<p->strain_total; i++) {
        S_out[i*maxsteps] = S[i];
        I_out[i*maxsteps] = I[i];
        R_out[i*maxsteps] = R[i];
    }

    
    /* Main loop: */
    do{
        /* This ensures the final step lands us on the final time point */
        if(1.1*hnext[0]>(p->t_max-t)){
            hnext[0] = p->t_max-t;
            oldh[0] = p->t_max-t;
            t=p->t_max;
            exitflag=0;
        }
        else{
            oldh[0] = hnext[0];
            t+=oldh[0];
        }
        if(t>=p->t_max) {
            t=p->t_max;
            exitflag=0;
        }
        /* This is where the equations are first solved */

        dynamic(S, I, R, DSDT, DIDT,  DRDT, zeta, p);

        
        /* Adjust the step size to maintain accuracy */
        for (i=0; i<p->strain_total; i++){
            S[i] = FMAX(S[i],0);
            S_SCALE[i]=fabs(S[i])+fabs(DSDT[i]*(*oldh))+TINY;
            I[i] = FMAX(I[i],0);
            I_SCALE[i]=fabs(I[i])+fabs(DIDT[i]*(*oldh))+TINY;
            R[i] = FMAX(R[i],0);
            R_SCALE[i]=fabs(R[i])+fabs(DRDT[i]*(*oldh))+TINY;
        }
        
        /* RK solver & adaptive step-size */

        rkqs(S, I, R, DSDT, DIDT, DRDT, oldh, hnext, S_SCALE, I_SCALE, R_SCALE, zeta, p);

        /* Make sure nothin has gone negative */
        for (i=0; i<p->strain_total; i++){
            S[i] = FMAX(S[i],0);
            I[i] = FMAX(I[i],0);
            R[i] = FMAX(R[i],0);
        }
        
        /* Update output */
        count++;
        T[count] = t;
        for (i=0; i<p->strain_total; i++) {
            S_out[count + i*maxsteps] = S[i];
            I_out[count + i*maxsteps] = I[i];
            R_out[count + i*maxsteps] = R[i];
        }


        /* For equilibrium check */
        for (i=0; i<p->strain_total; i++){
            SMIN[i] = FMIN(SMIN[i],S[i]);
            SMAX[i] = FMAX(SMAX[i],S[i]);
            IMIN[i] = FMIN(IMIN[i],I[i]);
            IMAX[i] = FMAX(IMAX[i],I[i]);
            RMIN[i] = FMIN(RMIN[i],R[i]);
            RMAX[i] = FMAX(RMAX[i],R[i]);
        }


        /* Check if we're close to equilibrium */
        if(t>nextcheck){
            exitflag = 0;
            for (i=0; i<p->strain_total; i++){
                if(fabs(SMAX[i]-SMIN[i])>p->eqtol || fabs(IMAX[i]-IMIN[i])>p->eqtol || fabs(RMAX[i]-RMIN[i])>p->eqtol){
                    exitflag = 1;
                    break;
		}
            }
            /* If close to equilibrium, then break */
            if(exitflag==0){
                t=p->t_max; 
                T[count] = t; 
                EQFLAG[0] = 1; 
                break; 
            } 
            
            /* If not, then reset min/max values for each class */
            nextcheck+=INTERVAL;
            for (i=0; i<p->strain_total; i++){
                    SMIN[i] = S[i];
                    SMAX[i] = S[i];
                    IMIN[i] = I[i];
                    IMAX[i] = I[i];
                    RMIN[i] = R[i];
                    RMAX[i] = R[i];
            }
        }
    }while(count<(maxsteps-1) && t<=p->t_max && exitflag);
    count++;

    
    /* Free memory */
    free(S);
    free(I);
    free(R);
    free(DSDT);
    free(DIDT);
    free(DRDT);
    free(S_SCALE);
    free(I_SCALE);
    free(R_SCALE);
    free(SMIN);
    free(IMIN);
    free(RMIN);
    free(SMAX);
    free(IMAX);
    free(RMAX);
    
    return count;
}

/***************************************
 * This generates the adaptive step-size
 **************************************/
void rkqs(double *S, double *I, double *R, double *DSDT, double *DIDT,  double *DRDT, double *oldh, double *hnext, double *S_SCALE, double *I_SCALE, double *R_SCALE, double *zeta, struct PARAM *p)
{
    double *S_temp, *I_temp,  *R_temp, *S_err, *I_err, *R_err;
    double htemp, errmax;
    int i, j, count;
    
    /* Allocate memory */
    S_temp = malloc(p->strain_total*sizeof(double));
    I_temp = malloc(p->strain_total*sizeof(double));
    R_temp = malloc(p->strain_total*sizeof(double));
    S_err = malloc(p->strain_total*sizeof(double));
    I_err = malloc(p->strain_total*sizeof(double));
    R_err = malloc(p->strain_total*sizeof(double));
    
    count = 0;
    for(;;)
    {
        rkck(S, I, R, DSDT, DIDT, DRDT, S_temp, I_temp,  R_temp, S_err, I_err,  R_err, *oldh, zeta, p); 
        
        errmax= 0.0;
        for(i=0;i<p->strain_total;i++){
            errmax= FMAX(errmax, fabs(S_err[i]/(S_SCALE[i]))); 
            errmax= FMAX(errmax, fabs(I_err[i]/(I_SCALE[i])));   
            errmax= FMAX(errmax, fabs(R_err[i]/(R_SCALE[i])));     
        }

        errmax/= EPS;

        if(errmax<=1.0) break;
        htemp= 0.9*(*oldh)*pow(errmax, -0.25);
        *oldh= (*oldh>=0.0 ? FMAX(htemp, 0.1*(*oldh)) : FMIN(htemp, 0.1*(*oldh)));
        count++;
            

        if(count>1e4){
            printf("%f\n",errmax);
            mexErrMsgTxt("stuck in loop!\n");
            break;
        }
    }    
    if(errmax > 1.89E-4) {
        *hnext= 0.9*(*oldh)*pow(errmax, -0.2);
    }
    else {
        *hnext= 5.0*(*oldh);
    }    
    *hnext = FMAX(*hnext, p->t_max/MAXSTEPS);
    
    for(i=0;i<p->strain_total;i++){
        S[i] = S_temp[i];
        I[i] = I_temp[i];
        R[i] = R_temp[i];
    }
    
    /* Free memory */
    free(S_temp);
    free(I_temp);
    free(S_err);
    free(I_err);
    free(R_temp);
    free(R_err);
}

/**************************************
 * Standard RK solver
 **************************************/
void rkck(double *S, double *I, double *R, double *DSDT, double *DIDT, double *DRDT, double *Sout, double *Iout, double *Rout, double *Serr, double *Ierr, double *Rerr, double oldh, double *zeta, struct PARAM *p){
    int i, j;
    double *Sk1, *Sk2, *Sk3, *Sk4, *Sk5, *Sk6, *Stemp;
    double *Ik1, *Ik2, *Ik3, *Ik4, *Ik5, *Ik6, *Itemp;
    double *Rk1, *Rk2, *Rk3, *Rk4, *Rk5, *Rk6, *Rtemp;
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, 
            dc6=c6-0.25;
    
    /* Allocate memory */    
    Sk1 = malloc(p->strain_total*sizeof(double));
    Sk2 = malloc(p->strain_total*sizeof(double));
    Sk3 = malloc(p->strain_total*sizeof(double));
    Sk4 = malloc(p->strain_total*sizeof(double));
    Sk5 = malloc(p->strain_total*sizeof(double));
    Sk6 = malloc(p->strain_total*sizeof(double));
    Stemp = malloc(p->strain_total*sizeof(double));

    Ik1 = malloc(p->strain_total*sizeof(double));
    Ik2 = malloc(p->strain_total*sizeof(double));
    Ik3 = malloc(p->strain_total*sizeof(double));
    Ik4 = malloc(p->strain_total*sizeof(double));
    Ik5 = malloc(p->strain_total*sizeof(double));
    Ik6 = malloc(p->strain_total*sizeof(double));
    Itemp = malloc(p->strain_total*sizeof(double));

    Rk1 = malloc(p->strain_total*sizeof(double));
    Rk2 = malloc(p->strain_total*sizeof(double));
    Rk3 = malloc(p->strain_total*sizeof(double));
    Rk4 = malloc(p->strain_total*sizeof(double));
    Rk5 = malloc(p->strain_total*sizeof(double));
    Rk6 = malloc(p->strain_total*sizeof(double));
    Rtemp = malloc(p->strain_total*sizeof(double));
    
    for(i=0;i<p->strain_total;i++){
        Stemp[i] = S[i] + b21*oldh*DSDT[i];
        Itemp[i] = I[i] + b21*oldh*DIDT[i];
        Rtemp[i] = R[i] + b21*oldh*DRDT[i];
    }
    dynamic(Stemp, Itemp, Rtemp, Sk2, Ik2, Rk2, zeta, p);


    
    for(i=0;i<p->strain_total;i++){
        Stemp[i] = S[i]+oldh*(b31*DSDT[i]+b32*Sk2[i]);
	Itemp[i] = I[i]+oldh*(b31*DIDT[i]+b32*Ik2[i]);
        Rtemp[i] = R[i]+oldh*(b31*DRDT[i]+b32*Rk2[i]);
    }    
    dynamic(Stemp, Itemp, Rtemp, Sk3, Ik3, Rk3, zeta, p);
    
    for(i=0;i<p->strain_total;i++){
        Stemp[i] = S[i]+oldh*(b41*DSDT[i]+b42*Sk2[i]+b43*Sk3[i]);
	Itemp[i] = I[i]+oldh*(b41*DIDT[i]+b42*Ik2[i]+b43*Ik3[i]);
        Rtemp[i] = R[i]+oldh*(b41*DRDT[i]+b42*Rk2[i]+b43*Rk3[i]);
    }


    dynamic(Stemp, Itemp, Rtemp, Sk4, Ik4, Rk4, zeta, p);
    
    for(i=0;i<p->strain_total;i++){
        Stemp[i] = S[i]+oldh*(b51*DSDT[i]+b52*Sk2[i]+b53*Sk3[i]+b54*Sk4[i]);
	Itemp[i] = I[i]+oldh*(b51*DIDT[i]+b52*Ik2[i]+b53*Ik3[i]+b54*Ik4[i]);
        Rtemp[i] = R[i]+oldh*(b51*DRDT[i]+b52*Rk2[i]+b53*Rk3[i]+b54*Rk4[i]);
    }



    dynamic(Stemp, Itemp, Rtemp, Sk5, Ik5, Rk5, zeta, p); 
    
    for(i=0;i<p->strain_total;i++){
        Stemp[i] = S[i]+oldh*(b61*DSDT[i]+b62*Sk2[i]+b63*Sk3[i]+b64*Sk4[i]+b65*Sk5[i]);
	Itemp[i] = I[i]+oldh*(b61*DIDT[i]+b62*Ik2[i]+b63*Ik3[i]+b64*Ik4[i]+b65*Ik5[i]);
        Rtemp[i] = R[i]+oldh*(b61*DRDT[i]+b62*Rk2[i]+b63*Rk3[i]+b64*Rk4[i]+b65*Rk5[i]);
    }

    dynamic(Stemp, Itemp, Rtemp, Sk6, Ik6, Rk6, zeta, p);
    
    for(i=0;i<p->strain_total;i++){
        Sout[i]= S[i]+oldh*(c1*DSDT[i]+c3*Sk3[i]+c4*Sk4[i]+c6*Sk6[i]);
        Serr[i]= oldh*(dc1*DSDT[i]+dc3*Sk3[i]+dc4*Sk4[i]+dc5*Sk5[i]+dc6*Sk6[i]);
        Iout[i]= I[i]+oldh*(c1*DIDT[i]+c3*Ik3[i]+c4*Ik4[i]+c6*Ik6[i]);
        Ierr[i]= oldh*(dc1*DIDT[i]+dc3*Ik3[i]+dc4*Ik4[i]+dc5*Ik5[i]+dc6*Ik6[i]);
        Rout[i]= R[i]+oldh*(c1*DRDT[i]+c3*Rk3[i]+c4*Rk4[i]+c6*Rk6[i]);
        Rerr[i]= oldh*(dc1*DRDT[i]+dc3*Rk3[i]+dc4*Rk4[i]+dc5*Rk5[i]+dc6*Rk6[i]);
    }

    /* Free memory */
    free(Sk1);
    free(Sk2);
    free(Sk3);
    free(Sk4);
    free(Sk5);
    free(Sk6);
    free(Stemp);
    free(Ik1);
    free(Ik2);
    free(Ik3);
    free(Ik4);
    free(Ik5);
    free(Ik6);
    free(Itemp);
    free(Rk1);
    free(Rk2);
    free(Rk3);
    free(Rk4);
    free(Rk5);
    free(Rk6);
    free(Rtemp);
}

/**************************************
 * Population and evolutionary dynamics
 **************************************/
void dynamic(double *S, double *I, double *R, double *DSDT, double *DIDT, double *DRDT, double *zeta, struct PARAM *p){
    
    int i, j;
    double N, allinfecteds;

    /* Population sums */
    N = 0;
    for(i=0;i<p->strain_total;i++){ 
	N = N + S[i] + I[i] +R [i];
    }

    allinfecteds = 0;
    for(i=0;i<p->strain_total;i++){ 
	allinfecteds = allinfecteds + I[i];
    }
    

    /* ODEs */
    for(i=0;i<p->strain_total;i++){
	DSDT[i] = p->a*(1 - p->q*N)*(S[i] + p->f*I[i]) + p->aR*(1-p->q*N)*R[i] - p->beta*S[i]*allinfecteds + p->gamma*I[i] - (p->b0*(1+((p->c1b*(1-exp(-p->c2b*zeta[i])))/(1-exp(-p->c2b)))))*S[i] - zeta[i]*S[i];
        DIDT[i] = p->beta*S[i]*allinfecteds - p->gamma*I[i] - (p->b0*(1+((p->c1b*(1-exp(-p->c2b*zeta[i])))/(1-exp(-p->c2b)))))*I[i]*(1+p->alpha);
	DRDT[i] = -p->bR*R[i] + zeta[i]*S[i];
    }
}

/***************************************
 * Return maximum of two inputs
 ***************************************/
double FMAX(double l, double r)
{
    if(l>r)return l;
    else   return r;
}

/***************************************
 * Return minimum of two inputs
 ***************************************/
double FMIN(double l, double r)
{
    if(l<r)return l;
    else   return r;
}
