#include <R.h> 
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>	// RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#define SWAP(a,b)   { int t; t=a; a=b; b=t; }  // Macro for swapping

// Cost functions  

double mll_var(int *n, double *sumstat, int end, int start, int seglen, double shape){
    double x3=*(sumstat+*n + *n +2+end)-*(sumstat+*n + *n +2+start);
    if(x3<=0){x3=0.00000000001;}
    return(seglen*(log(2*M_PI)+log(x3/seglen)+1)); /* M_PI is in Rmath.h  */
  } 

double mll_meanvar(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x2=*(sumstat+*n+1+end)-*(sumstat+*n+1+start); // this relies on the R code doing things in the correct order!
    double x=*(sumstat+end)-*(sumstat+start);
  double sigsq=(x2-((x*x)/seglen))/seglen;
  if(sigsq<=0){sigsq=0.00000000001;}
  return(seglen*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
}


double mll_mean(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x2=*(sumstat+*n+1+end)-*(sumstat+*n+1+start); // this relies on the R code doing things in the correct order!
    double x=*(sumstat+end)-*(sumstat+start);
  return(x2-(x*x)/seglen);
}

double mll_meanvar_exp(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x=*(sumstat+end)-*(sumstat+start);
  return(2*seglen*(log(x)-log(seglen)));
}

double mll_meanvar_gamma(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x=*(sumstat+end)-*(sumstat+start);
  return(2*seglen*shape*(log(x)-log(seglen*shape)));
}

double mll_meanvar_poisson(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x=*(sumstat+end)-*(sumstat+start);
  if(x==0){return(0);}
  else{return(2*x*(log(seglen)-log(x)));}
}

double mbic_var(int *n, double *sumstat, int end, int start, int seglen, double shape){
    double x3=*(sumstat+*n + *n +2+end)-*(sumstat+*n + *n +2+start);
    if(x3<=0){x3=0.00000000001;}
    return(seglen*(log(2*M_PI)+log(x3/seglen)+1)+log(seglen)); /* M_PI is in Rmath.h  */
  } 

double mbic_meanvar(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x2=*(sumstat+*n+1+end)-*(sumstat+*n+1+start); // this relies on the R code doing things in the correct order!
    double x=*(sumstat+end)-*(sumstat+start);
  double sigsq=(x2-((x*x)/seglen))/seglen;
  if(sigsq<=0){sigsq=0.00000000001;}
  return(seglen*(log(2*M_PI)+log(sigsq)+1)+log(seglen)); /* M_PI is in Rmath.h  */
}


double mbic_mean(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x2=*(sumstat+*n+1+end)-*(sumstat+*n+1+start); // this relies on the R code doing things in the correct order!
    double x=*(sumstat+end)-*(sumstat+start);
  return(x2-(x*x)/seglen+log(seglen));
}

double mbic_meanvar_exp(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x=*(sumstat+end)-*(sumstat+start);
  return(2*seglen*(log(x)-log(seglen))+log(seglen));
}

double mbic_meanvar_gamma(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x=*(sumstat+end)-*(sumstat+start);
  return(2*seglen*shape*(log(x)-log(seglen*shape))+log(seglen));
}

double mbic_meanvar_poisson(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x=*(sumstat+end)-*(sumstat+start);
  if(x==0){return(0);}
  else{return(2*x*(log(seglen)-log(x))+log(seglen));}
}

// code to choose cost function  
  
  const static struct {
    char *name;
    double (*func)(int *n, double *sumstat, int end, int start, int seglen, double shape);
  } function_map [] = {
{"var.norm", mll_var},
{"mean.norm", mll_mean},
{"meanvar.norm", mll_meanvar},
{"meanvar.exp", mll_meanvar_exp},
{"meanvar.gamma", mll_meanvar_gamma},
{"meanvar.poisson", mll_meanvar_poisson},
{"var.norm.mbic", mbic_var},
{"mean.norm.mbic", mbic_mean},
{"meanvar.norm.mbic", mbic_meanvar},
{"meanvar.exp.mbic", mbic_meanvar_exp},
{"meanvar.gamma.mbic", mbic_meanvar_gamma},
{"meanvar.poisson.mbic", mbic_meanvar_poisson},
  };

double call_function(const char *name,int *n, double *sumstat, int end, int start, int seglen, double shape)
{
  int k;
  
  for (k = 0; k <= (sizeof(function_map) / sizeof(function_map[0])); k++) {
    if (!strcmp(function_map[k].name, name) && function_map[k].func) {
      return(function_map[k].func(n,sumstat,end,start, seglen, shape));
    }
  }
	return(INFINITY);
}



void min_which(double *array,int n,double *minout,int *whichout){
  // Function to find minimum of an array with n elements that is put in min 
  *minout=*array;
	*whichout=0;
	int i;
	for(i=1;i<n;i++){
		if(*(array+i)< *minout){
			*minout= *(array+i); 
			*whichout=i;
		}
	}
}

void max_which(double *array,int n,double *maxout,int *whichout){
	// Function to find maximum of an array with n elements that is put in max 
	*maxout=*array;
	*whichout=0;
	int i;
	for(i=1;i<n;i++){
		if(*(array+i)> *maxout){
			*maxout= *(array+i);
			*whichout=i;
		}
	}
}

void order_vec( int a[], int n ){   
	int i, j;
	for(i = 0; i < n; i++){         // Make a pass through the array for each element
	  for(j = 1; j < (n-i); j++){			// Go through the array beginning to end
			if(a[j-1] > a[j])       // If the the first number is greater, swap it 
				SWAP(a[j-1],a[j]);   
		}
	}
}
