#include <mex.h>
#include <stdio.h>
#include <math.h>

static void usage()
{
  char str[1024];
  sprintf(str, "usage: gmm_distance(uy, SSigma, logweights)\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int k, d, N, K;
  const double* uy     ;
  const double* SSigma;
  const mwSize* sy     ;
  const mwSize* sm     ;
  const double* uyj;
  const double *logwts;
  double* energy;
  double* energyk;
  double invSSigma_tab[64];
  double sumlogSSigma;

  int i, j;
  double diff, tmp;

  if (nrhs < 3|| nlhs > 1)
    {
      usage();
      return;
    }
  plhs = plhs;

  for (k = 0; k < 3; ++k)
    if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]))
      {
	usage();
	return;
      }
  uy     = mxGetData(prhs[0]);
  SSigma = mxGetData(prhs[1]);
  logwts = mxGetData(prhs[2]);
  sy     = mxGetDimensions(prhs[0]);
  sm     = mxGetDimensions(prhs[1]);

  d = sy[0];
  N = sy[1];
  K = sm[1];

  plhs[0] = mxCreateDoubleMatrix(N, K, mxREAL);
  energy = (double*) mxGetPr(plhs[0]);

  //printf("%d %d %d\n", K, N, d);
  for (k = 0; k < K; ++k)
    {
      sumlogSSigma = logwts[k];
      for (i = 0; i < d; ++i)
	{
	  tmp              = SSigma[k * d + i];
	  invSSigma_tab[i] = tmp;
	  sumlogSSigma    += -log(tmp);
	}
      energyk = energy + k * N;
      for (j = 0; j < N; ++j)
	{
	  tmp = sumlogSSigma;
	  uyj = uy + j * d;
	  for (i = 0; i < d; ++i)
	    {
	      diff  = uyj[i];
	      tmp  += diff * diff * invSSigma_tab[i];
	    }
	  energyk[j] = tmp;
	}
    }
}
