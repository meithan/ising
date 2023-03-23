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

// Initial magnetization (-1 <= INIT_MAGN <= 1)
// Set INIT_MAGN = 0 for random
const float INIT_MAGN = 0.0;

// Data directory -- trailing slash optional
const char datadir[] = ".";

// Generations between full grid dumps
const int DUMP_GRID_EVERY = 1000;

/*===================================*/

int main(int argc, char* argv[]) {

  int gen, run;
  clock_t sclock, rclock, lclock;
  time_t ltime;
  double elapsed;
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
    model.set_magnetization(INIT_MAGN);
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
    seriesfile << "# Grid N = " << NGRID << "\n";
    seriesfile << "# Data columns\n";
    seriesfile << "# 1: generation number\n";
    seriesfile << "# 2: magnetization for this gen\n";
    seriesfile << "# 3: energy per site for this gen\n";

    // Open grid file for this run and write header
    if (NUM_RUNS == 1) {
      sprintf(fname, "%s/%s_grids.dat", datadir2, tempstr);
    } else {
      sprintf(fname, "%s/%s_r%03i_grids.dat", datadir2, tempstr, run);
    }
    printf("Recording grids in file %s\n",fname);
    gridsfile.open(fname);
    gridsfile << "# " << asctime(localtime(&ltime));
    gridsfile << "# Temperature = " << fixed << TEMP << "\n";
    gridsfile << "# Grid N = " << NGRID << "\n";

    printf("Initial magnetization M=%f\n", model.global_magnetization);
    printf("Simulating %i generations ...\n", NUM_GENS);

    // Dump state and grid of start state
    seriesfile << 0;
    seriesfile << " " << scientific << model.global_magnetization;
    seriesfile << " " << (double)(model.global_energy)/model.NCELLS;
    seriesfile << endl;
    gridsfile << "# GEN 0" << endl;
    for (int i = 0; i < NGRID; i++) {
      for (int j = 0; j < NGRID; j++) {
        if (model.grid[i][j] == 1) gridsfile << 1;
        else gridsfile << 0;
        // if (j < NGRID - 1) gridsfile << " ";
      }
      gridsfile << endl;
    }
    elapsed = (double)(clock()-rclock)/CLOCKS_PER_SEC;
    printf("[%.3f] gen 0 | E = %f | M = %f\n", elapsed, ((double)model.global_energy)/model.NCELLS, model.global_magnetization);

    // Simulate NUM_GENS generations
    for (gen = 1; gen <= NUM_GENS; gen++) {
      model.doGeneration();
      seriesfile << gen;
      seriesfile << " " << scientific << model.global_magnetization;
      seriesfile << " " << (double)(model.global_energy)/model.NCELLS;
      seriesfile << endl;
      if (gen % DUMP_GRID_EVERY == 0) {
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
        printf("[%.3f] gen %i | E = %f | M = %f\n", elapsed, gen, ((double)model.global_energy)/model.NCELLS, model.global_magnetization);
      }
    }

    ltime = time(NULL);
    elapsed = (double)(clock()-rclock)/CLOCKS_PER_SEC;
    printf("%s", asctime(localtime(&ltime)));
    printf("Run completed in %.3f s\n", elapsed);
    printf("=== Run %i/%i complete ===\n", run+1, NUM_RUNS);
    seriesfile.close();
    gridsfile.close();

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
