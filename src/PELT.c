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

static int *lastchangecpts;
static double *lastchangelike;
static int *checklist;
static double *tmplike;
static int *tmpt;

void FreePELT(error)
	int *error; /* Error code from PELT C function, non-zero => error */	
	{
	if(*error==0){
		free((void *)lastchangecpts);
	  free((void *)lastchangelike);
	  free((void *)checklist);
	  free((void *)tmplike);
	  free((void *)tmpt);
	}
}

void PELT_var_norm(y2,n,pen,cptsout,error)
  double *y2;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  int *error;   /* 0 by default, nonzero indicates error in code */
  {
	// R code does know.mean and fills mu if necessary

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
  int *lastchangecpts;
  lastchangecpts = (int *)calloc((*n+1)*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  double *lastchangelike;
  lastchangelike = (double *)calloc((*n+1),sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc((*n+1),sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc((*n+1),sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc((*n+1),sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_var();
	void min_which();
  
	lastchangelike[0]= -*pen; /* null (last changepoint at 0) */
	lastchangecpts[0]=0; lastchangecpts[*n+0]=0;	
	lastchangelike[1]=mll_var(*(y2+1),1);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
	lastchangelike[2]=mll_var(*(y2+2),2);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=2;
	lastchangelike[3]=mll_var(*(y2+3),3);
	lastchangecpts[3]=0; lastchangecpts[*n+3]=3;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_var(*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar]=minout;
		lastchangecpts[tstar]=checklist[whichout]; lastchangecpts[*n+tstar]=tstar;

		/* Update checklist for next iteration, last element is next tau */
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}
  free(tmpt);
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}

void binseg_var_norm(y2,n,pen,Q,cptsout,likeout,op_cps)
  double *y2;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
	int *Q;			/* Max number of changepoints */
  int *cptsout;    /* Q length vector of identified changepoint locations */
	double *likeout;		/* Q length vector of likelihood ratio values for changepoints in cptsout */
	int *op_cps;		/* Optimal number of changepoint for pen supplied */
  {
	// R code does know.mean and fills mu if necessary
	// must -0.5*mll_var to get same as R code

	double oldmax=1000,null,lambda[*n],maxout;
	int q,p,i,j,whichout,st,end;
	int tau[*Q+2]; // max ncpts is Q, +2 is for 0 and n
	tau[0]=0;
	tau[1]= *n;

	double mll_var();
	void max_which();
	void order_vec();

  for(q=0;q<*Q;q++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */
    for(p=0;p<*n;p++){lambda[p]=0;}
    i=1;
    st=tau[0]+1;
		end=tau[1];
    null= (-0.5) * mll_var(*(y2+end)-*(y2+st-1),end-st+1);
    for(j=2;j<(*n-2);j++){
      if(j==end){
        st=end+1;
				i=i+1;
				end=tau[i];
        null= (-0.5) * mll_var(*(y2+end)-*(y2+st-1),end-st+1);
      }
			else{
				if(((j-st)>1)&&((end-j)>1)){
	        lambda[j]= ((-0.5) * mll_var(*(y2+j)-*(y2+st-1),j-st+1)) + ((-0.5) * mll_var(*(y2+end)-*(y2+j),end-j)) - null;
				}
      }
    }
    max_which(lambda,*n,&maxout,&whichout);

    *(cptsout+q)=whichout;
		*(likeout+q)= (oldmax<=maxout) ? oldmax : maxout ;
    oldmax= *(likeout+q);
		tau[q+2]=whichout;
    order_vec(tau,q+3);
  }

	int stop=0;
	*op_cps=0;
	while((stop==0)&&(*op_cps < *Q)){
		if((2* *(likeout+*op_cps)) >= *pen){ (*op_cps)++; }
		else{ stop=1; }
	}
}


double mll_var(double x, int n){
	if(x<=0){x=0.00000000001;}
	return(n*(log(2*M_PI)+log(x/n)+1)); /* M_PI is in Rmath.h  */
}

void PELT_meanvar_norm(y2,y,n,pen,cptsout,error)
  double *y2;    /* Summary statistic for the time series */
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  int *error;   /* 0 by default, nonzero indicates error in code */
  {
	// R code does know.mean and fills mu if necessary

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
  int *lastchangecpts;
  lastchangecpts = (int *)calloc((*n+1)*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  double *lastchangelike;
  lastchangelike = (double *)calloc(*n+1,sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_meanvar();
	void min_which();
	
	lastchangelike[0]= -*pen;
	lastchangecpts[0]=0; lastchangecpts[*n+0]=0;
	lastchangelike[1]=mll_meanvar(*(y+1),*(y2+1),1);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
	lastchangelike[2]=mll_meanvar(*(y+2),*(y2+2),2);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=2;
	lastchangelike[3]=mll_meanvar(*(y+3),*(y2+3),3);
	lastchangecpts[3]=0; lastchangecpts[*n+3]=3;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar(*(y+tstar)-*(y+checklist[i]),*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar]=minout;
		lastchangecpts[tstar]=checklist[whichout]; lastchangecpts[*n+tstar]=tstar;

		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}
  free(tmpt);
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}



void binseg_meanvar_norm(y2,y,n,pen,Q,cptsout,likeout,op_cps)
  double *y2;    /* Summary statistic for the time series */
	double *y;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
	int *Q;			/* Max number of changepoints */
  int *cptsout;    /* Q length vector of identified changepoint locations */
	double *likeout;		/* Q length vector of likelihood ratio values for changepoints in cptsout */
	int *op_cps;		/* Optimal number of changepoint for pen supplied */
  {
	// R code does know.mean and fills mu if necessary
	// must -0.5*mll_var to get same as R code

	double oldmax=1000,null,lambda[*n],maxout;
	int q,p,i,j,whichout,st,end;
	int tau[*Q+2]; // max ncpts is Q, +2 is for 0 and n
	tau[0]=0;
	tau[1]= *n;

	double mll_meanvar();
	void max_which();
	void order_vec();

  for(q=0;q<*Q;q++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(p=0;p<*n;p++){lambda[p]=0;}
    i=1;
    st=tau[0]+1;
		end=tau[1];
    null= (-0.5) * mll_meanvar(*(y+end)-*(y+st-1),*(y2+end)-*(y2+st-1),end-st+1);
    for(j=2;j<(*n-2);j++){
      if(j==end){
        st=end+1;
				i=i+1;
				end=tau[i];
        null= (-0.5) * mll_meanvar(*(y+end)-*(y+st-1),*(y2+end)-*(y2+st-1),end-st+1);
      }
			else{
				if(((j-st)>1)&&((end-j)>1)){
	        lambda[j]= ((-0.5) * mll_meanvar(*(y+j)-*(y+st-1),*(y2+j)-*(y2+st-1),j-st+1)) + ((-0.5) * mll_meanvar(*(y+end)-*(y+j),*(y2+end)-*(y2+j),end-j)) - null;
				}
      }
    }
    max_which(lambda,*n,&maxout,&whichout);

    *(cptsout+q)=whichout;
		*(likeout+q)= (oldmax<=maxout) ? oldmax : maxout ;
    oldmax= *(likeout+q);
		tau[q+2]=whichout;
    order_vec(tau,q+3);
  }

	int stop=0;
	*op_cps=0;
	while((stop==0)&&(*op_cps < *Q)){
		if((2* *(likeout+*op_cps)) >= *pen){ (*op_cps)++; }
		else{ stop=1; }
	}
}

double mll_meanvar(double x, double x2, int n){
	double sigsq=(x2-((x*x)/n))/n;
	if(sigsq<=0){sigsq=0.00000000001;}
	return(n*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
}

void PELT_mean_norm(y2,y,n,pen,cptsout,error)
  double *y2;    /* Summary statistic for the time series */
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  int *error;   /* 0 by default, nonzero indicates error in code */
  {
	// R code does know.mean and fills mu if necessary

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
  int *lastchangecpts;
  lastchangecpts = (int *)calloc((*n+1)*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  double *lastchangelike;
  lastchangelike = (double *)calloc(*n+1,sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_mean();
	void min_which();
	
	lastchangelike[0]= -*pen;
	lastchangecpts[0]=0; lastchangecpts[*n+0]=0;
	lastchangelike[1]=mll_mean(*(y+1),*(y2+1),1);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
	lastchangelike[2]=mll_mean(*(y+2),*(y2+2),2);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=2;
	lastchangelike[3]=mll_mean(*(y+3),*(y2+3),3);
	lastchangecpts[3]=0; lastchangecpts[*n+3]=3;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_mean(*(y+tstar)-*(y+checklist[i]),*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar]=minout;
		lastchangecpts[tstar]=checklist[whichout]; lastchangecpts[*n+tstar]=tstar;

		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}
  free(tmpt);
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}


void binseg_mean_norm(y2,y,n,pen,Q,cptsout,likeout,op_cps)
  double *y2;    /* Summary statistic for the time series */
	double *y;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
	int *Q;			/* Max number of changepoints */
  int *cptsout;    /* Q length vector of identified changepoint locations */
	double *likeout;		/* Q length vector of likelihood ratio values for changepoints in cptsout */
	int *op_cps;		/* Optimal number of changepoint for pen supplied */
  {
	// R code does know.mean and fills mu if necessary
	// must -0.5*mll_var to get same as R code

	double oldmax=1000,null,lambda[*n],maxout;
	int q,p,i,j,whichout,st,end;
	int tau[*Q+2]; // max ncpts is Q, +2 is for 0 and n
	tau[0]=0;
	tau[1]= *n;

	double mll_mean();
	void max_which();
	void order_vec();

  for(q=0;q<*Q;q++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(p=0;p<*n;p++){lambda[p]=0;}
    i=1;
    st=tau[0]+1;
		end=tau[1];
    null= (-0.5) * mll_mean(*(y+end)-*(y+st-1),*(y2+end)-*(y2+st-1),end-st+1);
    for(j=2;j<(*n-2);j++){
      if(j==end){
        st=end+1;
				i=i+1;
				end=tau[i];
        null= (-0.5) * mll_mean(*(y+end)-*(y+st-1),*(y2+end)-*(y2+st-1),end-st+1);
      }
			else{
				if(((j-st)>1)&&((end-j)>1)){
	        lambda[j]= ((-0.5) * mll_mean(*(y+j)-*(y+st-1),*(y2+j)-*(y2+st-1),j-st+1)) + ((-0.5) * mll_mean(*(y+end)-*(y+j),*(y2+end)-*(y2+j),end-j)) - null;
				}
      }
    }
    max_which(lambda,*n,&maxout,&whichout);

    *(cptsout+q)=whichout;
		*(likeout+q)= (oldmax<=maxout) ? oldmax : maxout ;
    oldmax= *(likeout+q);
		tau[q+2]=whichout;
    order_vec(tau,q+3);
  }

	int stop=0;
	*op_cps=0;
	while((stop==0)&&(*op_cps < *Q)){
		if((2* *(likeout+*op_cps)) >= *pen){ (*op_cps)++; }
		else{ stop=1; }
	}
}

double mll_mean(double x, double x2, int n){
	return(x2-(x*x)/n);
}




void PELT_meanvar_exp(y,n,pen,cptsout,error)
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  int *error;   /* 0 by default, nonzero indicates error in code */
  {
	// R code does know.mean and fills mu if necessary

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
  int *lastchangecpts;
  lastchangecpts = (int *)calloc((*n+1)*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  double *lastchangelike;
  lastchangelike = (double *)calloc(*n+1,sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_meanvar_exp();
	void min_which();
	
	lastchangelike[0]= -*pen;
	lastchangecpts[0]=0; lastchangecpts[*n+0]=0;
	lastchangelike[1]=mll_meanvar_exp(*(y+1),1);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
	lastchangelike[2]=mll_meanvar_exp(*(y+2),2);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=2;
	lastchangelike[3]=mll_meanvar_exp(*(y+3),3);
	lastchangecpts[3]=0; lastchangecpts[*n+3]=3;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar_exp(*(y+tstar)-*(y+checklist[i]),tstar-checklist[i])+*pen;
		}
		min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar]=minout;
		lastchangecpts[tstar]=checklist[whichout]; lastchangecpts[*n+tstar]=tstar;

		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}
  free(tmpt);
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}


void binseg_meanvar_exp(y,n,pen,Q,cptsout,likeout,op_cps)
	double *y;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
	int *Q;			/* Max number of changepoints */
  int *cptsout;    /* Q length vector of identified changepoint locations */
	double *likeout;		/* Q length vector of likelihood ratio values for changepoints in cptsout */
	int *op_cps;		/* Optimal number of changepoint for pen supplied */
  {
	// R code does know.mean and fills mu if necessary
	// must -0.5*mll_var to get same as R code

	double oldmax=1000,null,lambda[*n],maxout;
	int q,p,i,j,whichout,st,end;
	int tau[*Q+2]; // max ncpts is Q, +2 is for 0 and n
	tau[0]=0;
	tau[1]= *n;

	double mll_meanvar_exp();
	void max_which();
	void order_vec();

  for(q=0;q<*Q;q++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(p=0;p<*n;p++){lambda[p]=0;}
    i=1;
    st=tau[0]+1;
		end=tau[1];
    null= (-0.5) * mll_meanvar_exp(*(y+end)-*(y+st-1),end-st+1);
    for(j=2;j<(*n-2);j++){
      if(j==end){
        st=end+1;
				i=i+1;
				end=tau[i];
        null= (-0.5) * mll_meanvar_exp(*(y+end)-*(y+st-1),end-st+1);
      }
			else{
				if(((j-st)>1)&&((end-j)>1)){
	        lambda[j]= ((-0.5) * mll_meanvar_exp(*(y+j)-*(y+st-1),j-st+1)) + ((-0.5) * mll_meanvar_exp(*(y+end)-*(y+j),end-j)) - null;
				}
      }
    }
    max_which(lambda,*n,&maxout,&whichout);

    *(cptsout+q)=whichout;
		*(likeout+q)= (oldmax<=maxout) ? oldmax : maxout ;
    oldmax= *(likeout+q);
		tau[q+2]=whichout;
    order_vec(tau,q+3);
  }

	int stop=0;
	*op_cps=0;
	while((stop==0)&&(*op_cps < *Q)){
		if((2* *(likeout+*op_cps)) >= *pen){ (*op_cps)++; }
		else{ stop=1; }
	}
}


double mll_meanvar_exp(double x, int n){
	return(2*n*(log(x)-log(n)));
}


void PELT_meanvar_gamma(y,n,pen,cptsout,shape,error)
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
	double *shape;		/* Shape parameter (fixed) */
  int *error;   /* 0 by default, nonzero indicates error in code */
  {
	// R code does know.mean and fills mu if necessary

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
  int *lastchangecpts;
  lastchangecpts = (int *)calloc((*n+1)*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  double *lastchangelike;
  lastchangelike = (double *)calloc(*n+1,sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_meanvar_gamma();
	void min_which();
	
	lastchangelike[0]= -*pen;
	lastchangecpts[0]=0; lastchangecpts[*n+0]=0;
	lastchangelike[1]=mll_meanvar_gamma(*(y+1),1,*shape);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
	lastchangelike[2]=mll_meanvar_gamma(*(y+2),2,*shape);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=2;
	lastchangelike[3]=mll_meanvar_gamma(*(y+3),3,*shape);
	lastchangecpts[3]=0; lastchangecpts[*n+3]=3;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar_gamma(*(y+tstar)-*(y+checklist[i]),tstar-checklist[i],*shape)+*pen;
		}
		min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar]=minout;
		lastchangecpts[tstar]=checklist[whichout]; lastchangecpts[*n+tstar]=tstar;

		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}
  free(tmpt);
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}

void binseg_meanvar_gamma(y,n,pen,Q,cptsout,likeout,op_cps,shape)
	double *y;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
	int *Q;			/* Max number of changepoints */
  int *cptsout;    /* Q length vector of identified changepoint locations */
	double *likeout;		/* Q length vector of likelihood ratio values for changepoints in cptsout */
	int *op_cps;		/* Optimal number of changepoint for pen supplied */
	double *shape;		/* Shape parameter (fixed) */
  {
	// R code does know.mean and fills mu if necessary
	// must -0.5*mll_var to get same as R code

	double oldmax=1000,null,lambda[*n],maxout;
	int q,p,i,j,whichout,st,end;
	int tau[*Q+2]; // max ncpts is Q, +2 is for 0 and n
	tau[0]=0;
	tau[1]= *n;

	double mll_meanvar_gamma();
	void max_which();
	void order_vec();

  for(q=0;q<*Q;q++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(p=0;p<*n;p++){lambda[p]=0;}
    i=1;
    st=tau[0]+1;
		end=tau[1];
    null= (-0.5) * mll_meanvar_gamma(*(y+end)-*(y+st-1),end-st+1,*shape);
    for(j=2;j<(*n-2);j++){
      if(j==end){
        st=end+1;
				i=i+1;
				end=tau[i];
        null= (-0.5) * mll_meanvar_gamma(*(y+end)-*(y+st-1),end-st+1,*shape);
      }
			else{
				if(((j-st)>1)&&((end-j)>1)){
	        lambda[j]= ((-0.5) * mll_meanvar_gamma(*(y+j)-*(y+st-1),j-st+1,*shape)) + ((-0.5) * mll_meanvar_gamma(*(y+end)-*(y+j),end-j,*shape)) - null;
				}
      }
    }
    max_which(lambda,*n,&maxout,&whichout);

    *(cptsout+q)=whichout;
		*(likeout+q)= (oldmax<=maxout) ? oldmax : maxout ;
    oldmax= *(likeout+q);
		tau[q+2]=whichout;
    order_vec(tau,q+3);
  }

	int stop=0;
	*op_cps=0;
	while((stop==0)&&(*op_cps < *Q)){
		if((2* *(likeout+*op_cps)) >= *pen){ (*op_cps)++; }
		else{ stop=1; }
	}
}

double mll_meanvar_gamma(double x, int n,double shape){
	return(2*n*shape*(log(x)-log(n*shape)));
}




void PELT_meanvar_poisson(y,n,pen,cptsout,error)
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  int *error;   /* 0 by default, nonzero indicates error in code */
  {
	// R code does know.mean and fills mu if necessary

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
  int *lastchangecpts;
  lastchangecpts = (int *)calloc((*n+1)*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  double *lastchangelike;
  lastchangelike = (double *)calloc(*n+1,sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_meanvar_poisson();
	void min_which();
	
	lastchangelike[0]= -*pen;
	lastchangecpts[0]=0; lastchangecpts[*n+0]=0;
	lastchangelike[1]=mll_meanvar_poisson(*(y+1),1);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
	lastchangelike[2]=mll_meanvar_poisson(*(y+2),2);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=2;
	lastchangelike[3]=mll_meanvar_poisson(*(y+3),3);
	lastchangecpts[3]=0; lastchangecpts[*n+3]=3;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar_poisson(*(y+tstar)-*(y+checklist[i]),tstar-checklist[i])+*pen;
		}
		min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar]=minout;
		lastchangecpts[tstar]=checklist[whichout]; lastchangecpts[*n+tstar]=tstar;

		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}
  free(tmpt);
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}


void binseg_meanvar_poisson(y,n,pen,Q,cptsout,likeout,op_cps)
	double *y;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
	int *Q;			/* Max number of changepoints */
  int *cptsout;    /* Q length vector of identified changepoint locations */
	double *likeout;		/* Q length vector of likelihood ratio values for changepoints in cptsout */
	int *op_cps;		/* Optimal number of changepoint for pen supplied */
  {
	// must -0.5*mll_var to get same as R code

	double oldmax=1000,null,lambda[*n],maxout;
	int q,p,i,j,whichout,st,end;
	int tau[*Q+2]; // max ncpts is Q, +2 is for 0 and n
	tau[0]=0;
	tau[1]= *n;

	double mll_meanvar_poisson();
	void max_which();
	void order_vec();

  for(q=0;q<*Q;q++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

    for(p=0;p<*n;p++){lambda[p]=0;}
    i=1;
    st=tau[0]+1;
		end=tau[1];
    null= (-0.5) * mll_meanvar_poisson(*(y+end)-*(y+st-1),end-st+1);
    for(j=2;j<(*n-2);j++){
      if(j==end){
        st=end+1;
				i=i+1;
				end=tau[i];
        null= (-0.5) * mll_meanvar_poisson(*(y+end)-*(y+st-1),end-st+1);
      }
			else{
				if(((j-st)>1)&&((end-j)>1)){
	        lambda[j]= ((-0.5) * mll_meanvar_poisson(*(y+j)-*(y+st-1),j-st+1)) + ((-0.5) * mll_meanvar_poisson(*(y+end)-*(y+j),end-j)) - null;
				}
      }
    }
    max_which(lambda,*n,&maxout,&whichout);

    *(cptsout+q)=whichout;
		*(likeout+q)= (oldmax<=maxout) ? oldmax : maxout ;
    oldmax= *(likeout+q);
		tau[q+2]=whichout;
    order_vec(tau,q+3);
  }

	int stop=0;
	*op_cps=0;
	while((stop==0)&&(*op_cps < *Q)){
		if((2* *(likeout+*op_cps)) >= *pen){ (*op_cps)++; }
		else{ stop=1; }
	}
}


double mll_meanvar_poisson(double x, int n){
	if(x==0){return(INFINITY);}
	else{return(2*x*(log(n)-log(x)));}
}

// IN PROGRESS
/*void PELT_ar_norm(data,max.lag,n,pen,cptsout,error)
/*  double *data;    /* the original time series */
/*	int *max.lag;			/* Maximum lag to consider */
/*	int *n;			/* Length of the time series */
/*  double *pen;  /* Penalty used to decide if a changepoint is significant */
/*  int *cptsout;    /* Vector of identified changepoint locations */
/*  int *error;   /* 0 by default, nonzero indicates error in code */
/*  {

	//int lastchangecpts[*n][2]; /* stores last changepoint locations   */
/*  int *lastchangecpts;
  lastchangecpts = (int *)calloc(*n*2,sizeof(int));
  if (lastchangecpts==NULL)   {
    *error = 1;
    goto err1;
  }

  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
/*  double *lastchangelike;
  lastchangelike = (double *)calloc(*n,sizeof(double));
  if (lastchangelike==NULL)   {
    *error = 2;
    goto err2;
  }

  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n,sizeof(int));
  if (checklist==NULL)   {
    *error = 3;
    goto err3;
  }

  int nchecklist;
	double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n,sizeof(double));
  if (tmplike==NULL)   {
    *error = 4;
    goto err4;
  }

	//int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n,sizeof(int));
  if (tmpt==NULL)   {
    *error = 5;
    goto err5;
  }

	int tstar,i,whichout,nchecktmp;

	double mll_ar();
	void min_which();
  
	/* create summary statistics from original data */
/*	sumstat=(double *)calloc((*n+1)*(*max.lag+1),sizeof(double));
	if(sumstat==NULL) {
		*error=6;
		goto err6;
	}

	for(i=0;i<(max.lag+1);i++){
		double tmp=0.0;
		for(j=0;j<(*n+1);j++){
			if((j-i)>=0){
				tmp+= (*(data+j)* *(data+j-i)) /* replaces the cumulative sum in R*/
/*				*(sumstat+(i*(*max.lag+1))+j)=tmp;
			}
		}
	}

/*	lastchangelike[0]=mll_ar(*(y2+1),*max.lag,1);
	lastchangecpts[0]=0; lastchangecpts[*n+0]=1; // no longer required */ 
/*	lastchangelike[1]=mll_ar(*(sumsatat+2*(*max.lag+1)),*max.lag,2);
	lastchangecpts[1]=0; lastchangecpts[*n+1]=2;  
	lastchangelike[2]=mll_ar(*(sumstat+3*(*max.lag+1)),*max.lag,3);
	lastchangecpts[2]=0; lastchangecpts[*n+2]=3;

	nchecklist=1;
	checklist[0]=2;
	for(tstar=4;tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */

/*    for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]-1] + mll_ar(*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		tmplike[nchecklist]=mll_var(*(y2+tstar),tstar);
		min_which(tmplike,nchecklist+1,&minout,&whichout); /*updates minout and whichout with min and which element */
/*		lastchangelike[tstar-1]=minout;
		if(whichout==nchecklist){		lastchangecpts[tstar-1]=0; lastchangecpts[*n+tstar-1]=tstar; } /* Null */
/*		else{			lastchangecpts[tstar-1]=checklist[whichout]; lastchangecpts[*n+tstar-1]=tstar; } /* Alt */
		
		/* Update checklist for next iteration, first element is next tau */
/*		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar-1]+*pen)){
				*(checklist+nchecktmp)=checklist[i];
				nchecktmp+=1;
			}
		}
		*(checklist+nchecktmp)=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar
	
	// put final set of changepoints together
	int ncpts=0;
	int last=*n;
	while(last!=0){
		*(cptsout+ncpts)=lastchangecpts[*n+last-1];
		last=lastchangecpts[last-1];
		ncpts+=1;
	}

	free(sumstat);
err6:  free(tmpt);	
err5:  free(tmplike);
err4:  free(checklist);
err3:  free(lastchangelike);
err2:  free(lastchangecpts);
err1:  return;
}


double mll_ar(double *x, int max.lag, int n){
		/* x is a vector containing the unscaled covariances starting with lag0 and ending with max.lag */
/*		if(n<max.lag){max.lag=n-1;}
		if(sum(x==0)>0){x[x==0]=1e-10;} /* change this to be a for loop */
/*		double mdl[max.lag];
		mdl[0]=(n/2)*log(2*M_PI* (*x)/n);
		if(max.lag==0){p=0;}
		else{
			mdl[1]=(n/2)*log(2*M_PI*(*x/n - (*(x+1)* *(x+1))/(*x *n)));
			if(mdl[0]<=mdl[1]){p=0;}
			else if(max.lag==1){p=1;}
			else{
				ind=0;
				p=2;
        while((ind==0) & (p<=max.lag)){
	        mdl[p]=(n/2)*(log(2*M_PI)+log(*x/n-(x[2:p]/n)%*%solve(toeplitz(x[1:(p-1)]/n),x[2:p]/n)))
          if(mdl[p]<mdl[p-1]){p=p+1;}
 	        else{ind=1;p=p-1;}
      	}
    	}
  	}
		if(p==0){ return(mdl[0]); }
		else{ return(mdl[p-1]); }
}

*/

void min_which(double *array,int n,double *minout,int *whichout){
	/* Function to find minimum of an array with n elements that is put in min */
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
	/* Function to find maximum of an array with n elements that is put in max */
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

