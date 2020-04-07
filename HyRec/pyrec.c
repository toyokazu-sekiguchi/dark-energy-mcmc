#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#include "pyrec.h"

HRATEEFF rate_table;
TWO_PHOTON_PARAMS twog_params;
REC_COSMOPARAMS param;
int firstTime = 0;
double *xe_output, *tm_output, **Dfnu_hist, *Dfminus_Ly_hist[3];
double logstart;

void hyrec_init() {

  /* Build effective rate table */
  //char *buffer = (char *) malloc (1024);
  //getcwd (buffer, 1024);
  //chdir(HYRECPATH);
  rate_table.logTR_tab = create_1D_array(NTR);
  rate_table.TM_TR_tab = create_1D_array(NTM);
  rate_table.logAlpha_tab[0] = create_2D_array(NTM, NTR);
  rate_table.logAlpha_tab[1] = create_2D_array(NTM, NTR);
  rate_table.logR2p2s_tab = create_1D_array(NTR);
  read_rates(&rate_table);

  /* Read two-photon rate tables */
  read_twog_params(&twog_params);
  //chdir(buffer);
  //free(buffer);
}

void rec_build_history_wrap(double tcmb, double obh2, double odmh2,
                            double okh2, double odeh2, double w0, double wa,
                            double yp, double nnu, double mnu[3]){
  int j;
  
  param.T0 = tcmb;
  param.obh2 = obh2;
  param.omh2 = obh2+odmh2;
  param.okh2 = okh2;
  param.odeh2 = odeh2;
  param.w0 = w0;
  param.wa = wa;
  param.Y = yp;
  param.Nnueff = nnu;
  for(j=0; j<3; j++) {
    param.mnu[j] = mnu[j];
  };
  if (firstTime ==0) {
    //char *buffer = (char *) malloc (1024);
    //getcwd (buffer, 1024);
    //chdir(HYRECPATH);
    init_mnu();
    //chdir(buffer);
    //free(buffer);
  }
  param.fsR = param.meR = 1.;
  rec_set_derived_params(&param);

  if (firstTime == 0) {
    hyrec_init();
    firstTime=1;
    logstart = -log(1.+ZSTART);
    xe_output          = create_1D_array(param.nz);
    tm_output          = create_1D_array(param.nz);
  }
  Dfnu_hist          = create_2D_array(NVIRT, param.nzrt);
  Dfminus_Ly_hist[0] = create_1D_array(param.nzrt);
  Dfminus_Ly_hist[1] = create_1D_array(param.nzrt);
  Dfminus_Ly_hist[2] = create_1D_array(param.nzrt);
  rec_build_history(&param, &rate_table, &twog_params, xe_output, tm_output,Dfnu_hist, Dfminus_Ly_hist);
  free_2D_array(Dfnu_hist, NVIRT);
  free(Dfminus_Ly_hist[0]);
  free(Dfminus_Ly_hist[1]);
  free(Dfminus_Ly_hist[2]);

}

double hyrec_xe(double a){
  double loga = log(a);
  if (loga < logstart) {
    return xe_output[0];
  }
  return rec_interp1d(logstart, DLNA, xe_output, param.nz, loga);
}

double hyrec_tm(double a){
  double loga = log(a);
  if (loga < logstart) {
    return tm_output[0];
  }
  return rec_interp1d(logstart, DLNA, tm_output, param.nz, loga);
}
