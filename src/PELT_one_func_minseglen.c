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
//#include "cost_general_functions.c" // commented out as already imported in BinSeg_one_func_minseglen.c

//static int *lastchangecpts;
//static double *lastchangelike;
static int *checklist;
static double *tmplike;
static int *tmpt;


void FreePELT(error)
int *error; // Error code from PELT C function, non-zero => error 	
{
  if(*error==0){
    free((void *)checklist);
    free((void *)tmplike);
    free((void *)tmpt);
  }
}

void PELT(cost_function,sumstat,n,pen,cptsout,error, shape, minseglen, lastchangelike, lastchangecpts, numchangecpts)
     char **cost_function; //Descibe the cost function used i.e. norm.mean.cost (change in mean in normal distributed data)  
     double *sumstat;  //array of summary statistics of the time series  
     int *n;			// Length of the time series 
     double *pen;  // Penalty used to decide if a changepoint is significant 
     int *cptsout;    // Vector of identified changepoint locations 
     int *error;   // 0 by default, nonzero indicates error in code 
     double *shape; // only used when cost_func is the gamma likelihood 
     int *minseglen; //minimum segment length 
     double *lastchangelike; // stores likelihood up to that time using optimal changepoint locations up to that time 
     int *lastchangecpts; // stores last changepoint locations 
     int *numchangecpts; //stores the current number of changepoints 
{
  
  
  
  //int checklist[*n]; 
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 1;
    goto err1;
  }
  
  int nchecklist;
  double minout;
  
  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 2;
    goto err2;
  }
  
  //int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 3;
    goto err3;
  }
  
  int tstar,i,whichout,nchecktmp, ind,j;
  void min_which();
  double call_function();

 
  lastchangelike[0]= -*pen;  // null (last changepoint at 0) 
  lastchangecpts[0]=0; 
   
  for(j=*minseglen;j<(2*(*minseglen));j++){
      lastchangelike[j] = call_function(*cost_function,n,sumstat,j,0,j,*shape); 
    }
    
  
  for(j=1;j<(2*(*minseglen));j++){ 
    lastchangecpts[j] = 0;
  }

  for(j=*minseglen;j<(2*(*minseglen));j++){ 
    numchangecpts[j] =1;
  }
  nchecklist=2;
  checklist[0]=0;
  checklist[1]=*minseglen;
 
  
  for(tstar=2*(*minseglen);tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); // checks if user has interrupted the R session and quits if true 
    
   if ((lastchangelike[tstar]) == 0){ 
     for(i=0;i<(nchecklist+1);i++){
        tmplike[i]=lastchangelike[checklist[i]] + call_function(*cost_function,n,sumstat,tstar,checklist[i],tstar-checklist[i], *shape) + *pen;
      }
    
   min_which(tmplike,nchecklist,&minout,&whichout); //updates minout and whichout with min and which element 
    lastchangelike[tstar]=minout;
    lastchangecpts[tstar]=checklist[whichout]; 
    numchangecpts[tstar]=numchangecpts[lastchangecpts[tstar]]+1; 
    
    // Update checklist for next iteration, first element is next tau 
      nchecktmp=0;
    for(i=0;i<nchecklist;i++){
      if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
        *(checklist+nchecktmp)=checklist[i];
        nchecktmp+=1;
      }
    }
   }
   
    *(checklist+nchecklist)=tstar-(*minseglen-1);  // atleast 1 obs per seg
    nchecklist+=1;
   // end taustar
}


// put final set of changepoints together
   int ncpts=0;
  int last=*n;
  while(last!=0){
    *(cptsout + ncpts) = last; 
    last=lastchangecpts[last];
    ncpts+=1;
  }
  free(tmpt);
 err3:  free(tmplike);
 err2:  free(checklist);
 err1:  return;
}

// Cost functions  
/*

double mll_var(int *n, double *sumstat, int end, int start, int seglen, double shape){
    double x3=*(sumstat+*n + *n +2+end)-*(sumstat+*n + *n +2+start);
    if(x3<=0){x3=0.00000000001;}
    return(seglen*(log(2*M_PI)+log(x3/seglen)+1)); // M_PI is in Rmath.h  
  } 

double mll_meanvar(int *n, double *sumstat, int end, int start, int seglen, double shape){
  double x2=*(sumstat+*n+1+end)-*(sumstat+*n+1+start); // this relies on the R code doing things in the correct order!
    double x=*(sumstat+end)-*(sumstat+start);
  double sigsq=(x2-((x*x)/seglen))/seglen;
  if(sigsq<=0){sigsq=0.00000000001;}
  return(seglen*(log(2*M_PI)+log(sigsq)+1)); // M_PI is in Rmath.h  
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

// code to choose cost function  
  
  const static struct {
    char *name;
    double (*func)(int *n, double *sumstat, int end, int start, int seglen, double shape);
  } function_map [] = {
{ "norm.var", mll_var},
{"norm.mean", mll_mean},
{"norm.meanvar", mll_meanvar},
{"exp", mll_meanvar_exp},
{"gamma", mll_meanvar_gamma},
{"poisson", mll_meanvar_poisson},
  };

double call_function(const char *name,int *n, double *sumstat, int end, int start, int seglen, double shape)
{
  int k;
  
  for (k = 0; k <= (sizeof(function_map) / sizeof(function_map[0])); k++) {
    if (!strcmp(function_map[k].name, name) && function_map[k].func) {
      return function_map[k].func(n,sumstat,end,start, seglen, shape);
    }
  }
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
*/
