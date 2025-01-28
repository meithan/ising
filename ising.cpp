#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include "IsingModel.h"
#include "utils.h"
using namespace std;

/*===================================*/

/* MODEL PARAMETERS */

// Size of grid (NGRID x NGRID)
const int NGRID = 100;

// Temperature (in units of J/k, so that Tc=2.2693706, dimensionless)
// >> TEMPERATURE PASSED AS SOLE COMMAND LINE ARGUMENT
double TEMP;

// Number of generations to simulate per run
const int NUM_GENS = 10000;

// Number of runs to simulate
const int NUM_RUNS = 1;

// Initial magnetization -- determines how the initial spin states are set
// Accepted values for INIT_MAGN_MODE: INIT_MAGN_AUTO or INIT_MAGN_MANUAL
// If INIT_MAGN_AUTO, the initial magnetization will be set to the exact
// equilibrium value, so as to shorten the initial transient: this is M = 0
// when T >= Tc, or an expression below Tc.
// If instead INIT_MAGN_MANUAL is used, the user specifies the initial
// magnetization through INIT_MAGN (must be in [-1, 1]).
// In both cases, the spins will be initially set randomly but weighted so
// that the resulting initial magnetization is near (but not exactly) the
// desired value.
const int INIT_MAGN_MODE = INIT_MAGN_AUTO;
const float INIT_MAGN = 0.0;

// Data directory -- trailing slash optional
const char datadir[] = ".";

// Generations between full grid dumps
// Set this value to zero for no grid dumps
const int DUMP_GRID_EVERY = 1000;

/*===================================*/

int main(int argc, char* argv[]) {

  int gen, run;
  clock_t sclock, rclock, lclock;
  time_t ltime;
  double elapsed;
  double _init_magn;
  char fname[192];
  char datadir2[128];
  char tempstr[7];
  ofstream seriesfile, gridsfile;

  sclock = clock();

  // Read temperature from command line
  if (argc < 2) {
    cerr << "Must provide temperature as first argument!" << endl;
    return 1;
  } else {
    TEMP = atof(argv[1]);
  }
  sprintf(tempstr, "T%.3f", TEMP);

  // Remove slash to datadir if present
  if (datadir[strlen(datadir)-1] == '/') {
    strncat(datadir2, datadir, strlen(datadir)-1);
  } else {
    strcpy(datadir2, datadir);
  }
  // Create model
  IsingModel model(NGRID, TEMP);

  // Determine initial magnetization
  if (INIT_MAGN_MODE == INIT_MAGN_AUTO) {
    if (TEMP < TEMP_CRIT) {
      _init_magn = pow(1 - pow(sinh(2/TEMP), -4), 0.125);
    } else {
      _init_magn = 0.0;
    }
  } else if ((INIT_MAGN_MODE == INIT_MAGN_MANUAL)) {
    _init_magn = INIT_MAGN;
  } else {
    printf("INIT_MAGN_MODE must be either INIT_MAGN_AUTO or INIT_MAGN_MANUAL. Aborting.\n");
    return 1;
  }

  printf("Temperature T=%f\n", TEMP);
  printf("%i x %i Ising model\n", NGRID, NGRID);
  printf("%i run%s\n", NUM_RUNS, NUM_RUNS > 1 ? "s" : "");
  printf("%i generations\n", NUM_GENS);
  printf("Datadir is %s/\n", datadir2);

  // Loop over runs
  for (run = 0; run < NUM_RUNS; run++) {

    printf("\n=== Starting run %i/%i ===\n", run+1, NUM_RUNS);
    rclock = clock();
    ltime = time(NULL);
    cout << asctime(localtime(&ltime));

    // Reset model
    model.reset_stats();
    model.cur_gen = 0;
    model.set_magnetization(_init_magn);
    model.update_energy();
    model.update_magnetization();

    // Open series file for this run and write header
    if (NUM_RUNS == 1) {
      sprintf(fname, "%s/%s_series.dat", datadir2, tempstr);
    } else {
      sprintf(fname, "%s/%s_r%03i_series.dat", datadir2, tempstr, run);
    }
    seriesfile.open(fname);
    printf("Recording time series in file %s\n",fname);
    seriesfile << "# " << asctime(localtime(&ltime));
    seriesfile << "# Temperature = " << fixed << TEMP << "\n";
    seriesfile << "# " << NGRID << " x " << NGRID << " grid\n";
    seriesfile << "# Columns: Magnetization, Energy\n";

    // Open grid file for this run and write header
    if (DUMP_GRID_EVERY > 0) {
      if (NUM_RUNS == 1) {
        sprintf(fname, "%s/%s_grids.dat", datadir2, tempstr);
      } else {
        sprintf(fname, "%s/%s_r%03i_grids.dat", datadir2, tempstr, run);
      }
      printf("Recording grids in file %s\n",fname);
      gridsfile.open(fname);
      gridsfile << "# " << asctime(localtime(&ltime));
      gridsfile << "# Temperature = " << fixed << TEMP << "\n";
      gridsfile << "# " << NGRID << " x " << NGRID << " grid\n";
    }

    printf("Initial magnetization M=%f\n", model.global_magnetization);
    printf("Simulating %i generations ...\n", NUM_GENS);

    // Dump state and grid of start state
    seriesfile << scientific << model.global_magnetization;
    seriesfile << " " << (double)(model.global_energy)/model.NCELLS;
    seriesfile << endl;
    if (DUMP_GRID_EVERY > 0) {
      gridsfile << "# GEN 0" << endl;
      for (int i = 0; i < NGRID; i++) {
        for (int j = 0; j < NGRID; j++) {
          if (model.grid[i][j] == 1) gridsfile << 1;
          else gridsfile << 0;
          // if (j < NGRID - 1) gridsfile << " ";
        }
        gridsfile << endl;
      }
    }
    elapsed = (double)(clock()-rclock)/CLOCKS_PER_SEC;
    printf("[%.3f] gen 0 | M = %f | E = %f\n", elapsed, model.global_magnetization, ((double)model.global_energy)/model.NCELLS);

    // Simulate NUM_GENS generations
    for (gen = 1; gen <= NUM_GENS; gen++) {
      model.doGeneration();
      seriesfile << scientific << model.global_magnetization;
      seriesfile << " " << (double)(model.global_energy)/model.NCELLS;
      seriesfile << endl;
      if (DUMP_GRID_EVERY > 0 && gen % DUMP_GRID_EVERY == 0) {
        gridsfile << "# GEN " << gen << endl;
        for (int i = 0; i < NGRID; i++) {
          for (int j = 0; j < NGRID; j++) {
            if (model.grid[i][j] == 1) gridsfile << 1;
            else gridsfile << 0;
            // if (j < NGRID - 1) gridsfile << " ";
          }
          gridsfile << endl;
        }
      }
      if (gen % (NUM_GENS/10) == 0) {
        elapsed = (double)(clock()-rclock)/CLOCKS_PER_SEC;
        printf("[%.3f] gen %i | M = %f | E = %f\n", elapsed, gen, model.global_magnetization, ((double)model.global_energy)/model.NCELLS);
      }
    }

    ltime = time(NULL);
    elapsed = (double)(clock()-rclock)/CLOCKS_PER_SEC;
    seriesfile << "# Finished " << asctime(localtime(&ltime));
    seriesfile << "# Elapsed " << elapsed << " s";
    seriesfile.close();
    if (DUMP_GRID_EVERY > 0) {
      gridsfile << "# Finished " << asctime(localtime(&ltime));
      gridsfile << "# Elapsed " << elapsed << " s";
      gridsfile.close();
    }
    printf("%s", asctime(localtime(&ltime)));
    printf("Run completed in %.3f s\n", elapsed);
    printf("=== Run %i/%i complete ===\n", run+1, NUM_RUNS);

  }

  if (NUM_RUNS > 1) {
    printf("\n=== All runs complete! ===\n");
    ltime = time(NULL);
    elapsed = (double)(clock()-sclock)/CLOCKS_PER_SEC;
    printf("Finished: %s", asctime(localtime(&ltime)));
    printf("Total elapsed: %.1f s\n", elapsed);
  }

  return 0;

}
