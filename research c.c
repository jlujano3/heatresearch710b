/* This is John Lujano's research with Dr. Tausch for the paper entitled
An Optimization Method for Moving Interface Problems Goverend by the Heat Equation*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ONETHRD 0.333333333333333333333333333333
#define M_SQPI 1.77245385090552
#define M_SQPINV 0.564189583547756287

double rKnown(double t){
  double g;
  g=t+1;
  return g;
}

double dRKnown(double t){
  double g;
  g=1;
  return g;
}

double u(double x,double xo,double t,double to){
  double g;
  g=exp(-(x-xo)*(x-xo)/(4.0*(t+to)))/sqrt(t+to);
  return g;
}

double du(double x,double xo,double t,double to){
  double g;
  g = exp(-(x-xo)*(x-xo)/(4.0*(t+to)))/sqrt(t+to);
  g *= -(x-xo)/(2*(t+to));
  return g;
}

double uo(double x,double xo,double to){
  double g;
  //g=u(x,xo,0.0,to);
  //g=0;
  return g;
}

double qKnown(double xo,double t,double to){
  double r,r1,q;
  r  = rKnown(t);
  r1 = dRKnown(t);
  q = du(r, xo, t, to) + u(r, xo, t, to)*r1;

  return q;
}

double kv(double x,double xo,double t,double to){
  double g;
  g=0.5*M_SQPINV*exp(-(x-xo)*(x-xo)/(4.0*(t-to)));
  return g;
}

/*Used to test the single layer potential with trivial kv*/

double kv2(double x,double xo,double t,double to){
  double g;
  g=1.0;
  return g;
}

double kvtntn(){
  double g;
  g=0.5*M_SQPINV;
  return g;
}

/*Used to test the single layer potential with trivial kv*/

double kvtntn2(){
  double g;
  g=1.0;
  return g;
}

double kd(double x,double xo,double t,double to){
  double g;
  g=exp(-(x-xo)*(x-xo)/(4.0*(t-to)))*(x-xo)*0.25*M_SQPINV/((t-to));
  return g;
}

double kdtntn(double t){
  double g;
  g=dRKnown(t)*0.25*M_SQPINV;
  return g;
}

double sLay(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv(r[n],r[0],tn,0)*q[0]/sqrt(tn-0);
  for (i=1; i<n; i++){
    sum += kv(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  return sum*h;
}

double sLay2(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv2(r[n],r[0],tn,0)*q[0]/sqrt(tn-0);
  for (i=1; i<n; i++){
    sum += kv2(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  return sum*h;
}

double dLay(int n,double h,double *r,double *u){
  double tn=n*h;
  int i;
  double sum=0.5*kd(r[n],r[0],tn,0)*u[0]/sqrt(tn-0);

  for (i=1; i<n; i++){
    sum += kd(r[n],r[i],tn,i*h)*u[i]/sqrt(tn-i*h);
  }
  return sum*h;
}

double mun(int n,double h){
  double tn=n*h;
  int i;
  double sum=0.5*h/sqrt(tn-0);
  for (i=1; i<n; i++){
    sum += h/sqrt(tn-i*h);
  }
  return 2.0*sqrt(tn)-sum;
}

double initPot(int n, double h, double xo,double to,double *r){
  int numInts=10000;
  double b=10.0, z;
  double tn=(n)*h, s4t=sqrt(4.0*tn);
  double a=r[n]/s4t;

  if (a<b){
    a=a;
  }
  else{
    a=b;
  }

  double delta=(b+a)/numInts;
  double sum=0.5*(exp(-a*a)*uo(r[n]-a*s4t,xo,to)+exp(-b*b)*uo(r[n]+b*s4t,xo,to));
  int i;
  for(i=1; i<numInts; i++){
    z = -a+delta*i;
    sum += exp(-z*z)*uo(r[n]+z*s4t,xo,to);
  }
  return sum*delta*M_SQPINV;
}

double qTestExact(double t){
  double g;
  g=exp(t)*sqrt(M_PI)*erf(sqrt(t));
  return g;
}

double qEval(double t){
  double ans=(1+t*t+t*t*t*t);
  return ans;
}

double q2Exact(double t){
  double ans=M_PI*(1+3*t*t/8+35*t*t*t*t/128);
  return ans;
}

double sLayAdj(int n, double h, double *r, double *qt){
  double ktntjsum=0;
  double tn=n*h;
  int i;
  for (i=1;i<n-1;i++){
    ktntjsum+=h*kv(r[n],r[i],tn,i*h)*qt[i]/(sqrt(tn-i*h));
  }
  return ktntjsum;
}

double sLayAdj2(int n, double h, double *r, double *qt){
  double ktntjsum=0;
  double tn=n*h;
  int i;
  for (i=1;i<n-1;i++){
    ktntjsum+=h*kv2(r[n],r[i],tn,i*h)*qt[i]/(sqrt(tn-i*h));
  }
  return ktntjsum;
}

double qto(int n, double h){
  double sum=M_PI/2;
  double tn=n*h;
  int i;
  for (i=1;i<n-1;i++){
    sum-=h*sqrt(tn-i*h)/(tn*sqrt(i*h));
  }
  return sum;
}

double munAdj(int n, double h){
  double tn=n*h;
  double sum=0.5*M_PI*tn;
  int i;
  for (i=1;i<n-1;i++){
    sum-=h*sqrt(i*h)/(sqrt(tn-i*h));
  }
  return sum/sqrt(tn);
}

int main(int nargs, char *argv[]){

  double xo=2.0,to=2.0,tmax=1.0,dt;
  int i,numInts=4096,verbose=0;
  double *r, *t, *qKnownList, *q;
  double *qTest, *uPasser, *sLays;
  double *dLays,*munValsD,*initPote,*denom,*qtList;
  double *sLaysAdj,*dLaysAdj,*munValsAdj,*initPoteAdj;
  double *denomAdj,*qAdj;
  double *qEvalList, *qtoList;
 
  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
      case 'v': verbose = atoi( argv[i]+3 );
        break;
      case 'T': tmax = atof( argv[i]+3 );
        break;
      case 'N': numInts = atoi( argv[i]+3 );
        break;
      }

  dt=tmax/numInts;
  
  r=(double*)malloc( (numInts+1)*sizeof(double) );
  t=(double*)malloc( (numInts+1)*sizeof(double) );
  qKnownList=(double*)malloc( (numInts+1)*sizeof(double) );
  uPasser=(double*)malloc( (numInts+1)*sizeof(double) );
  //For the SLay Test
  qTest=(double*)malloc( (numInts+1)*sizeof(double) );
  //These are the qs for the second integration strategy testing module
  qEvalList=(double*)malloc( (numInts+1)*sizeof(double) );
  for(i=0;i<numInts+1;i++){
    t[i]=0+i*dt;
    r[i]=rKnown(t[i]);
    qKnownList[i]=qKnown(xo,t[i],to);
    uPasser[i]=u(r[i],xo,i*dt,to);
    qTest[i]=exp(t[i]);
    qEvalList[i]=qEval(i*dt)/sqrt(i*dt);
  }
  q=(double*)malloc( (numInts+1)*sizeof(double) );
  q[0]=qKnown(xo,t[0],to);
  qAdj=(double*)malloc( (numInts+1)*sizeof(double) );
  qAdj[0]=sqrt(M_PI)*(u(rKnown(0),xo,0,to)-uo(0,xo,to));
  //qAdj[0]=0.03;
  sLays=(double*)malloc( (numInts+1)*sizeof(double) );
  sLays[0]=0;
  dLays=(double*)malloc( (numInts+1)*sizeof(double) );
  dLays[0]=0;
  munValsD=(double*)malloc( (numInts+1)*sizeof(double) );
  munValsD[0]=0;
  initPote=(double*)malloc( (numInts+1)*sizeof(double) );
  initPote[0]=0;
  denom=(double*)malloc( (numInts+1)*sizeof(double) );
  denom[0]=0;
  denomAdj=(double*)malloc( (numInts+1)*sizeof(double) );
  denomAdj[0]=0;  
  sLaysAdj=(double*)malloc( (numInts+1)*sizeof(double) );
  sLaysAdj[0]=0;
  munValsAdj=(double*)malloc( (numInts+1)*sizeof(double) );
  munValsAdj[0]=0;
  initPoteAdj=(double*)malloc( (numInts+1)*sizeof(double) );
  initPoteAdj[0]=0;
  qtoList=(double*)malloc( (numInts+1)*sizeof(double) );
  qtoList[0]=0;
  
  int j;
  for(j=1;j<numInts+1;j++){
    sLays[j]=sLay(j,dt,r,q);
    sLaysAdj[j]=sLayAdj(j,dt,r,qAdj);
    dLays[j]=dLay(j,dt,r,uPasser);
    munValsD[j]=kdtntn(j*dt)*mun(j,dt)*uPasser[j];
    initPote[j]=initPot(j,dt,xo,to,r);
    denom[j]=mun(j,dt)*kvtntn();
    denomAdj[j]=munAdj(j,dt)*kvtntn();
    qtoList[j]=qto(j,dt)*kv(r[j],r[0],j*dt,0)*qAdj[0];
    q[j]=(-0.5*uPasser[j]-sLays[j]+dLays[j]+munValsD[j]+initPote[j])/denom[j];
    qAdj[j]=(-0.5*uPasser[j]-sLaysAdj[j]-qtoList[j]+dLays[j]+munValsD[j]+initPote[j])/denomAdj[j];
    printf("%s  -> %s\n","Adjoint","Original");
    //Percent Error
    //printf("%f  %f\n",(qAdj[j]-qKnownList[j])/qKnownList[j]*100,(q[j]-qKnownList[j])/qKnownList[j]*100);
    //Adjoint and Regular Error
    printf("%f  %f\n",qAdj[j]-qKnownList[j],q[j]-qKnownList[j]);
    //Print Adjoint
    //printf("%10.9f\n",qAdj[j]);
    //Print Non-Adjoint q
    //printf("%10.9f\n",q[j]);
  }

#if 0
  printf("%s\n","u(rKnown(0),xo,0,to)");
  printf("%10.9f\n",u(rKnown(0),xo,0,to));
  printf("%s\n","uo(0,xo,to)");
  printf("%10.9f\n",uo(0,xo,to));
#endif

#if 0
  /*Use this module to test the new integration technique*/
  //I am feeding the q values into qEval List
  printf("%s\n","Integration Test With 2 Singularities");
  printf("%s\n","Numerical");
  printf("%10.9f\n",munAdj(numInts,dt)*kvtntn2()*qEvalList[numInts]+sLayAdj2(numInts,dt,r,qEvalList)+qto(numInts,dt)*qEval(0)*kv2(r[numInts],r[0],(numInts+1)*dt,0));
  printf("%s\n","Analytical");
  printf("%10.9f\n",q2Exact(tmax));
  printf("%s\n","Difference");
  printf("%10.9f\n",munAdj(numInts,dt)*kvtntn2()*qEvalList[numInts]+sLayAdj2(numInts,dt,r,qEvalList)+qto(numInts,dt)*qEval(0)*kv2(r[numInts],r[0],(numInts)*dt,0)-q2Exact(tmax));
#endif  
  
#if 0
  /*Build An Initial Potential Testing Module*/
  double ans2=initPot(10,dt,xo,to,r);
  double analyt=0.622853;
  printf("%s\n","Initial Potential Tester");
  printf("%s\n","Numerical");
  printf("%f\n",ans2);
  printf("%s\n","Analytical");
  printf("%f\n",analyt);
#endif 

#if 0
  /*This is proof that if we feed the single layer scheme the correct q value at each arrival, then the scheme is accurate with a trivial kv*/
  int testpoint=10;
  double ans3=sLay2(testpoint,dt,r,qTest)+mun(testpoint,dt)*qTest[testpoint];
  double analyt2=qTestExact(testpoint*dt);
  printf("%s\n","SLay Test");
  printf("%s\n","Numerical");
  printf("%10.9f\n",ans3);
  printf("%s\n","Analytical");
  printf("%10.9f\n",analyt2);
#endif  
}