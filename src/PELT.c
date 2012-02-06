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

void PELT_var_norm(y2,n,pen,cptsout)
  double *y2;    /* Summary statistic for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  {
	// R code does know.mean and fills mu if necessary

	int lastchangecpts[*n][2]; /* stores last changepoint locations   */
	double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
	int checklist[*n],nchecklist;
	double tmplike[*n],minout;
	int tmpt[*n];
	int tstar,i,whichout,nchecktmp;

	double mll_var();
	void min_which();
	
	lastchangelike[0]=mll_var(*(y2+1),1);
	lastchangecpts[0][0]=0; lastchangecpts[0][1]=1;
	lastchangelike[1]=mll_var(*(y2+2),2);
	lastchangecpts[1][0]=0; lastchangecpts[1][1]=2;
	lastchangelike[2]=mll_var(*(y2+3),3);
	lastchangecpts[2][0]=0; lastchangecpts[2][1]=3;

	nchecklist=1;
	checklist[0]=2;
	for(tstar=4;tstar<(*n+1);tstar++){
		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]-1] + mll_var(*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		tmplike[nchecklist]=mll_var(*(y2+tstar),tstar);
		min_which(tmplike,nchecklist+1,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar-1]=minout;
		if(whichout==nchecklist){		lastchangecpts[tstar-1][0]=0; lastchangecpts[tstar-1][1]=tstar; } /* Null */
		else{			lastchangecpts[tstar-1][0]=checklist[whichout]; lastchangecpts[tstar-1][1]=tstar; } /* Alt */
		
		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
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
		*(cptsout+ncpts)=lastchangecpts[last-1][1];
		last=lastchangecpts[last-1][0];
		ncpts+=1;
	}
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
	if(x<0){x=0.00000000001;}
	return(n*(log(2*M_PI)+log(x/n)+1)); /* M_PI is in Rmath.h  */
}

void PELT_meanvar_norm(y2,y,n,pen,cptsout)
  double *y2;    /* Summary statistic for the time series */
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  {
	// R code does know.mean and fills mu if necessary

	int lastchangecpts[*n][2]; /* stores last changepoint locations   */
	double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
	int checklist[*n],nchecklist;
	double tmplike[*n],minout;
	int tmpt[*n];
	int tstar,i,whichout,nchecktmp;

	double mll_meanvar();
	void min_which();
	
	lastchangelike[0]=mll_meanvar(*(y+1),*(y2+1),1);
	lastchangecpts[0][0]=0; lastchangecpts[0][1]=1;
	lastchangelike[1]=mll_meanvar(*(y+2),*(y2+2),2);
	lastchangecpts[1][0]=0; lastchangecpts[1][1]=2;
	lastchangelike[2]=mll_meanvar(*(y+3),*(y2+3),3);
	lastchangecpts[2][0]=0; lastchangecpts[2][1]=3;

	nchecklist=1;
	checklist[0]=2;
	for(tstar=4;tstar<(*n+1);tstar++){
		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]-1] + mll_meanvar(*(y+tstar)-*(y+checklist[i]),*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		tmplike[nchecklist]=mll_meanvar(*(y+tstar),*(y2+tstar),tstar);
		min_which(tmplike,nchecklist+1,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar-1]=minout;
		if(whichout==nchecklist){		lastchangecpts[tstar-1][0]=0; lastchangecpts[tstar-1][1]=tstar; } /* Null */
		else{			lastchangecpts[tstar-1][0]=checklist[whichout]; lastchangecpts[tstar-1][1]=tstar; } /* Alt */
		
		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
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
		*(cptsout+ncpts)=lastchangecpts[last-1][1];
		last=lastchangecpts[last-1][0];
		ncpts+=1;
	}
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
	return(n*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
}

void PELT_mean_norm(y2,y,n,pen,cptsout)
  double *y2;    /* Summary statistic for the time series */
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  {
	// R code does know.mean and fills mu if necessary

	int lastchangecpts[*n][2]; /* stores last changepoint locations   */
	double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
	int checklist[*n],nchecklist;
	double tmplike[*n],minout;
	int tmpt[*n];
	int tstar,i,whichout,nchecktmp;

	double mll_mean();
	void min_which();
	
	lastchangelike[0]=mll_mean(*(y+1),*(y2+1),1);
	lastchangecpts[0][0]=0; lastchangecpts[0][1]=1;
	lastchangelike[1]=mll_mean(*(y+2),*(y2+2),2);
	lastchangecpts[1][0]=0; lastchangecpts[1][1]=2;
	lastchangelike[2]=mll_mean(*(y+3),*(y2+3),3);
	lastchangecpts[2][0]=0; lastchangecpts[2][1]=3;

	nchecklist=1;
	checklist[0]=2;
	for(tstar=4;tstar<(*n+1);tstar++){
		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]-1] + mll_mean(*(y+tstar)-*(y+checklist[i]),*(y2+tstar)-*(y2+checklist[i]),tstar-checklist[i])+*pen;
		}
		tmplike[nchecklist]=mll_mean(*(y+tstar),*(y2+tstar),tstar);
		min_which(tmplike,nchecklist+1,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar-1]=minout;
		if(whichout==nchecklist){		lastchangecpts[tstar-1][0]=0; lastchangecpts[tstar-1][1]=tstar; } /* Null */
		else{			lastchangecpts[tstar-1][0]=checklist[whichout]; lastchangecpts[tstar-1][1]=tstar; } /* Alt */
		
		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
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
		*(cptsout+ncpts)=lastchangecpts[last-1][1];
		last=lastchangecpts[last-1][0];
		ncpts+=1;
	}
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




void PELT_meanvar_exp(y,n,pen,cptsout)
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
  {
	// R code does know.mean and fills mu if necessary

	int lastchangecpts[*n][2]; /* stores last changepoint locations   */
	double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
	int checklist[*n],nchecklist;
	double tmplike[*n],minout;
	int tmpt[*n];
	int tstar,i,whichout,nchecktmp;

	double mll_meanvar_exp();
	void min_which();
	
	lastchangelike[0]=mll_meanvar_exp(*(y+1),1);
	lastchangecpts[0][0]=0; lastchangecpts[0][1]=1;
	lastchangelike[1]=mll_meanvar_exp(*(y+2),2);
	lastchangecpts[1][0]=0; lastchangecpts[1][1]=2;
	lastchangelike[2]=mll_meanvar_exp(*(y+3),3);
	lastchangecpts[2][0]=0; lastchangecpts[2][1]=3;

	nchecklist=1;
	checklist[0]=2;
	for(tstar=4;tstar<(*n+1);tstar++){
		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]-1] + mll_meanvar_exp(*(y+tstar)-*(y+checklist[i]),tstar-checklist[i])+*pen;
		}
		tmplike[nchecklist]=mll_meanvar_exp(*(y+tstar),tstar);
		min_which(tmplike,nchecklist+1,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar-1]=minout;
		if(whichout==nchecklist){		lastchangecpts[tstar-1][0]=0; lastchangecpts[tstar-1][1]=tstar; } /* Null */
		else{			lastchangecpts[tstar-1][0]=checklist[whichout]; lastchangecpts[tstar-1][1]=tstar; } /* Alt */
		
		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
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
		*(cptsout+ncpts)=lastchangecpts[last-1][1];
		last=lastchangecpts[last-1][0];
		ncpts+=1;
	}
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


void PELT_meanvar_gamma(y,n,pen,cptsout,shape)
	double *y;		/*Summary statistics for the time series */
	int *n;			/* Length of the time series */
  double *pen;  /* Penalty used to decide if a changepoint is significant */
  int *cptsout;    /* Vector of identified changepoint locations */
	double *shape;		/* Shape parameter (fixed) */
  {
	// R code does know.mean and fills mu if necessary

	int lastchangecpts[*n][2]; /* stores last changepoint locations   */
	double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
	int checklist[*n],nchecklist;
	double tmplike[*n],minout;
	int tmpt[*n];
	int tstar,i,whichout,nchecktmp;

	double mll_meanvar_gamma();
	void min_which();
	
	lastchangelike[0]=mll_meanvar_gamma(*(y+1),1,*shape);
	lastchangecpts[0][0]=0; lastchangecpts[0][1]=1;
	lastchangelike[1]=mll_meanvar_gamma(*(y+2),2,*shape);
	lastchangecpts[1][0]=0; lastchangecpts[1][1]=2;
	lastchangelike[2]=mll_meanvar_gamma(*(y+3),3,*shape);
	lastchangecpts[2][0]=0; lastchangecpts[2][1]=3;

	nchecklist=1;
	checklist[0]=2;
	for(tstar=4;tstar<(*n+1);tstar++){
		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]-1] + mll_meanvar_gamma(*(y+tstar)-*(y+checklist[i]),tstar-checklist[i],*shape)+*pen;
		}
		tmplike[nchecklist]=mll_meanvar_gamma(*(y+tstar),tstar,*shape);
		min_which(tmplike,nchecklist+1,&minout,&whichout); /*updates minout and whichout with min and which element */
		lastchangelike[tstar-1]=minout;
		if(whichout==nchecklist){		lastchangecpts[tstar-1][0]=0; lastchangecpts[tstar-1][1]=tstar; } /* Null */
		else{			lastchangecpts[tstar-1][0]=checklist[whichout]; lastchangecpts[tstar-1][1]=tstar; } /* Alt */
		
		/* Update checklist for next iteration, first element is next tau */
		nchecktmp=0;
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
		*(cptsout+ncpts)=lastchangecpts[last-1][1];
		last=lastchangecpts[last-1][0];
		ncpts+=1;
	}
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

