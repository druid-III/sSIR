#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define BUFF 128

struct sDY
{
  double             *a_RB;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sDY DY;

struct sSA
{
  double             *a_RB;
  double              sum;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sSA SA;

struct sVA
{
  double             *a_RB;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sVA VA;

struct sDELAY
{
  double             *a_RB;
  double              sum;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sDELAY DELAY;

struct sCORR
{
  double             *a_RB;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sCORR CORR;

struct sFIR
{
  double             *a_RB;
  double             *a_H;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sFIR FIR;

struct sIIRdir
{
  double             *a_RB1;
  double             *a_RB2;
  double             *a_H1;
  double             *a_H2;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sIIRdir IIRdir;

struct sIIRcan
{
  double             *a_RB;
  double             *a_H1;
  double             *a_H2;
  const unsigned int  L;
  unsigned int        p;
};
typedef struct sIIRcan IIRcan;

typedef unsigned int UI;
typedef const unsigned int const CUIC;

//............................................................................//

//............................................................................//

bool fn_DELAY(double *, UI, DY *);

bool fn_SA(double *, UI, SA *);

bool fn_VARIANCE(double *, UI, VA *);
bool fn_CORR(double *, double *, UI, CORR *);

bool fn_FIR(double *, UI, FIR *);
bool fn_IIRdir(double *, UI, IIRdir *);
bool fn_IIRcan(double *, UI, IIRcan *);


/*

bool fn_AGC(double *, UI, AGC *);
/**/


// ========================================================================== //
bool fn_DELAY(double *ps, UI ns, DY *O)
{                                                                             // printf("DELAY procession\n");
  double temp;

  for(UI n=0; n<ns; n++)
  {
    temp = O->a_RB[O->p];                                                     // oldest 2 buff

    O->a_RB[O->p] = ps[n]; ps[n] = temp;                                      // new - input? buff 2 output
    
    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
} 
// ========================================================================== //
bool fn_SA(double *ps, UI ns, SA *O)
{                                                                             // printf("SA procession\n");
  for(UI n=0; n<ns; n++)
  {
    O->sum -= O->a_RB[O->p];                                                  // oldest out

    O->a_RB[O->p] = ps[n];   ps[n] = 0.;                                      // new - input
    
    O->sum += O->a_RB[O->p]; ps[n] = O->sum/(double)O->L;                     // sum all of they

    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
} 
// ========================================================================== //
bool fn_VARIANCE(double *ps, UI ns, VA *O)
{                                                                             // printf("VARIANCE procession\n");
  double sum1, sum2 = 0.;

  UI aEnd = O->L - 1;                                                         // max sequential number beff cell 4 back step

  for(UI n=0; n<ns; n++)
  {
    O->a_RB[aEnd - O->p] = ps[n];

    for(UI c=0; c<O->L; c++)
    {
      sum1 += O->a_RB[((aEnd - O->p) + c)%O->L];
      sum2 += 
           O->a_RB[((aEnd - O->p) + c)%O->L]*O->a_RB[((aEnd - O->p) + c)%O->L];
    }

    ps[n] = ((sum2 - sum1*(sum1/(double)O->L))/((double)O->L - 1.));

    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
}
// ========================================================================== //
bool fn_CORR(double *ps, double *cf, UI ns, CORR *O)
{                                                                             // printf("CORRELATOR procession\n");
  UI aEnd = O->L - 1;                                                         // max sequential number beff cell 4 back step

  for(UI n=0; n<ns; n++)
  {
    O->a_RB[aEnd - O->p] = ps[n]; ps[n] = 0.;

    for(UI c=0; c<O->L; c++)
    {
      ps[n] += O->a_RB[((aEnd - O->p) + c)%O->L]*cf[c];
    }

    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
}
// ========================================================================== //
bool fn_FIR(double *ps, UI ns, FIR *O)
{                                                                             // printf("FIR procession\n");
  UI aEnd = O->L - 1;                                                         // max sequential number beff cell 4 back step

  for(UI n=0; n<ns; n++)
  {
    O->a_RB[aEnd - O->p] = ps[n]; ps[n] = 0.;

    for(UI c=0; c<O->L; c++)
    {
      ps[n] += O->a_RB[((aEnd - O->p) + c)%O->L]*O->a_H[c];
    }

    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
}
// ========================================================================== //
bool fn_IIRdir(double *ps, UI ns, IIRdir *O)
{                                                                             // printf("direct IIR procession\n");
  UI aEnd = O->L - 1;                                                         // max sequential number beff cell 4 back step

  for(UI n=0; n<ns; n++)
  {  
    O->a_RB1[aEnd - O->p] = ps[n];

    ps[n] = O->a_RB1[aEnd - O->p]*O->a_H1[0];

    for(UI c=1; c<O->L; c++)
    {
      ps[n] += O->a_RB1[((aEnd - O->p) + c)%O->L]*O->a_H1[c];
    }

    ps[n] += O->a_RB2[O->p];                                                  // feedBack from F(n - 1)
        
    for(UI c=1; c<O->L; c++)
    {
      ps[n] += O->a_RB2[((aEnd - O->p) + c)%O->L]*O->a_H2[c];
    }

    O->a_RB2[aEnd - O->p] = ps[n];                                            // y -> F(n - 1)

    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
}

// ========================================================================== //
bool fn_IIRcan(double *ps, UI ns, IIRcan *O)
{                                                                             // printf("canonical IIR procession\n");
  UI aEnd = O->L - 1;                                                         // max sequential number beff cell 4 back step

  for(UI n=0; n<ns; n++)
  {    
    for(UI c=1; c<O->L; c++)
    {
      ps[n] -= O->a_RB[((aEnd - O->p) + c)%O->L]*O->a_H2[c];                  // feedBack from F(n - 1)
    }

    O->a_RB[aEnd - O->p] = ps[n];                                             // y -> F(n - 1)

    ps[n] = O->a_RB[aEnd - O->p]*O->a_H1[0];

    for(UI c=1; c<O->L; c++)
    {
      ps[n] += O->a_RB[((aEnd - O->p) + c)%O->L]*O->a_H1[c];
    }
                                              
    O->p = (O->p + 1)%O->L;                                                   // ring buffer pointer step up
  }
  return true;
}

// -------------------------------------------------------------------------- //
/**/
int main()
{
// -------------------------------------------------------------------------- //
  printf ("// -DELAY-processing------------------------------------------->\n");

  double a_d[4]; for(UI c=0; c<4; c++) a_d[c] = 0.;

  double pa_4D[BUFF]; for(UI c=0; c<BUFF; c++) pa_4D[c] = 0.; 
  pa_4D[16]  = 1.; 

  pa_4D[110] = 1.; 
  pa_4D[111] = 2.; 
  pa_4D[112] = 3.;

  pa_4D[64] =  1.; 
  pa_4D[65] =  1.; 
  pa_4D[66] =  1.;
  pa_4D[67] =  1.; 
  pa_4D[68] = -1.; 
  pa_4D[69] = -1.;
  pa_4D[70] = -1.; 
  pa_4D[71] = -1.; 


  DY O_DY =
  {
    .L    = 4,

    .a_RB = a_d,
    .p    = 0
  };

  fn_DELAY(pa_4D, BUFF, &O_DY);

  for(UI c=0; c<BUFF; c++) printf ("%-4u sampl is %e\n", c, pa_4D[c]);
// -------------------------------------------------------------------------- //
  printf ("// -SA-processing---------------------------------------------->\n");

  double a_s[4]; for(UI c=0; c<4; c++) a_s[c] = 0.;

  double pa_4SA[BUFF]; for(UI c=0; c<BUFF; c++) pa_4SA[c] = 0.; 
  pa_4SA[16]  = 1.; 

  pa_4SA[110] = 1.; 
  pa_4SA[111] = 2.; 
  pa_4SA[112] = 3.;

  pa_4SA[64] =  1.; 
  pa_4SA[65] =  1.; 
  pa_4SA[66] =  1.;
  pa_4SA[67] =  1.; 
  pa_4SA[68] = -1.; 
  pa_4SA[69] = -1.;
  pa_4SA[70] = -1.; 
  pa_4SA[71] = -1.; 


  SA O_SA =
  {
    .L    = 4,

    .a_RB = a_s,
    .sum  = 0.,
    .p    = 0
  };

  fn_SA(pa_4SA, BUFF, &O_SA);

  for(UI c=0; c<BUFF; c++) printf ("%-4u sampl is %e\n", c, pa_4SA[c]);
// -------------------------------------------------------------------------- //
  printf ("// -FIR-processing--------------------------------------------->\n");

  double a_b[3]; for(UI c=0; c<3; c++) a_b[c] = 0.;
  double a_k[3]; for(UI c=0; c<3; c++) a_k[c] = (double)(c+1);

  double pa_4FIR[BUFF]; for(UI c=0; c<BUFF; c++) pa_4FIR[c] = 0.; 
  pa_4FIR[16]  = 1.; 

  pa_4FIR[110] = 1.; 
  pa_4FIR[111] = 2.; 
  pa_4FIR[112] = 3.;

  FIR O_FIR =
  {
    .L   = 3,

    .a_RB = a_b,
    .a_H  = a_k,
    .p    = 0
  };

  fn_FIR(pa_4FIR, BUFF, &O_FIR);

  for(UI c=0; c<BUFF; c++) printf ("%-4u sampl is %e\n",c , pa_4FIR[c]);

// -------------------------------------------------------------------------- //
  printf ("// -IIR-processing--------------------------------------------->\n");

  double a_b1[2]; for(UI c=0; c<2; c++) a_b1[c] = 0.;
  double a_b2[2]; for(UI c=0; c<2; c++) a_b2[c] = 0.;

  double a_k1[2] = 
    {
      -3.0,
       2.0,
    };

  double a_k2[2] = 
   {

       NAN,
      -0.5,
   };

  double pa_4dIIR[BUFF]; for(UI c=0; c<BUFF; c++) pa_4dIIR[c] = 0.; 
  pa_4dIIR[16]  = 1.; 

  for(UI c=64; c<BUFF; c++) pa_4dIIR[c] = 1.; 

  IIRdir O_IIRdir =
  {
    .L     = 2,

    .a_RB1 = a_b1,
    .a_RB2 = a_b2,
    .a_H1  = a_k1,
    .a_H2  = a_k2,
    .p     = 0
  };

  fn_IIRdir(pa_4dIIR, BUFF, &O_IIRdir);

  for(UI c=0; c<BUFF; c++) printf ("%-4u sampl is %e\n",c , pa_4dIIR[c]);
// -------------------------------------------------------------------------- //
  printf ("// -IIR-processing--------------------------------------------->\n");

  double a_B[2]; for(UI c=0; c<2; c++) a_B[c] = 0.;

  double a_K1[2] = 
    {
      -3.0,
       2.0,
    };

  double a_K2[2] = 
   {

       NAN,
      -0.5,
   };

  double pa_4cIIR[BUFF]; for(UI c=0; c<BUFF; c++) pa_4cIIR[c] = 0.; 
  pa_4cIIR[16]  = 1.; 

  for(UI c=64; c<BUFF; c++) pa_4cIIR[c] = 1.; 

  IIRcan O_IIRcan =
  {
    .L    = 2,

    .a_RB = a_B,
    .a_H1 = a_k1,
    .a_H2 = a_k2,
    .p    = 0
  };

  fn_IIRcan(pa_4cIIR, BUFF, &O_IIRcan);

  for(UI c=0; c<BUFF; c++) printf ("%-4u sampl is %e\n",c , pa_4cIIR[c]);
// -------------------------------------------------------------------------- //

  return 1;
}
