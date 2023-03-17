#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "IsingModel.h"
#include "utils.h"

/*============================================================================*/

/*==================================\\
|| Ising Model class implementation ||
\\==================================*/

/*============================================================================*/

// CONSTRUCTORS

// Basic constructor
// Only grid size and temperature must be provided.
// No sample stats nor moving averages are tracked.
IsingModel::IsingModel (int p_NGRID, double p_TEMP) {

  // Set model parameters
  TEMP = p_TEMP;
  NGRID = p_NGRID;
  NCELLS = NGRID*NGRID;
  NUM_SAMPLES = 0;
  SAMPLE_MIN = 0;
  SAMPLE_MAX = 0;
  START_GEN = 1;
  NUM_DATA = 0;
  sample_magn = NULL;
  sample_mean = NULL;
  sample_var = NULL;
  sample_M2 = NULL;
  sample_npts = NULL;
  sample_size = NULL;
  sample_cells = NULL;
  track_samples = false;

  // Do common tasks
  common_constructor();

}

/*============================================================================*/

// Constructor for advanced model: tracks sample stats and moving averages

IsingModel::IsingModel (int p_NGRID, double p_TEMP, int p_NUM_SAMPLES, int p_SAMPLE_MIN, int p_SAMPLE_MAX, int p_START_GEN, int p_NUM_DATA) {

  double a;

  // Set model parameters
  TEMP = p_TEMP;
  NGRID = p_NGRID;
  NCELLS = NGRID*NGRID;
  NUM_SAMPLES = p_NUM_SAMPLES;
  SAMPLE_MIN = p_SAMPLE_MIN;
  SAMPLE_MAX = p_SAMPLE_MAX;
  START_GEN = p_START_GEN;
  NUM_DATA = p_NUM_DATA;
  track_samples = true;

  // Allocate sample stats arrays
  sample_magn = (double*) malloc(NUM_SAMPLES*sizeof(double));
  sample_mean = (double*) malloc(NUM_SAMPLES*sizeof(double));
  sample_var = (double*) malloc(NUM_SAMPLES*sizeof(double));
  sample_M2 = (double*) malloc(NUM_SAMPLES*sizeof(double));
  sample_npts = (int*) malloc(NUM_SAMPLES*sizeof(int));

  // Determine sample sizes and allocate sample cell lists
  sample_size = (int*) malloc(NUM_SAMPLES*sizeof(int));
  sample_cells = (int**) malloc(NUM_SAMPLES*sizeof(int*));
  a = pow((double)SAMPLE_MAX/SAMPLE_MIN, 1.0/(NUM_SAMPLES-1));
  for (int s = 0; s < NUM_SAMPLES; s++) {
    sample_size[s] = round(SAMPLE_MIN*pow(a,s));
    sample_cells[s] = (int*) malloc(sample_size[s]*sizeof(int));
  }

  // Randomly pick cells to be sampled
  pickSamples();

  // Allocate running mean array
  rundata = (double*) malloc(NUM_SAMPLES*sizeof(double));

  // Do common tasks
  common_constructor();

}

/*============================================================================*/

// Common tasks done by all constructors

void IsingModel::common_constructor () {

  int i;

  // Default flip strategy. User must change after class instantiation.
  flip_strategy = STRATEGY_SHUFFLE;

  // Default dynamics. User must change after class instantiation
  trans_dynamics = DYNAMICS_METROPOLIS;

  // Allocate spin grid (initialize grid_copy to null ptr)
  grid = (int**) malloc(NGRID*sizeof(int*));
  for (i = 0; i < NGRID; i++) {
    grid[i] = (int*) malloc(NGRID*sizeof(int));
  }
  grid_copy = NULL;

  // Dead cells -- turned OFF by default
  dead_cells = NULL;
  useDeadCells = false;

  // Allocate and initialize flip_order
  flip_order = (int*) malloc(NCELLS*sizeof(int));
  for (i = 0; i < NCELLS; i++) {
    flip_order[i] = i;
  }

  // Reset all stats
  reset_stats();
  cur_gen = 0;

  // Seed RNG
  srand(time(NULL));

}

/*============================================================================*/

// Reset all stats
void IsingModel::reset_stats () {
  global_energy = 0;
  global_magnetization = 0.0;
  global_mean = 0.0;
  global_variance = 0.0;
  global_npoints = 0;
  global_M2 = 0;
  if (track_samples) {
    for (int s = 0; s < NUM_SAMPLES; s++) {
      sample_mean[s] = 0.0;
      sample_var[s] = 0.0;
      sample_npts[s] = 0;
      sample_M2[s] = 0.0;
    }
  }
  run_mean = 0.0;
  run_var = 0.0;
  nextdata = 0;
}

/*============================================================================*/

// Computes the current global magnetization of the grid
void IsingModel::update_energy () {
  global_energy = 0;
  for (int i = 0; i < NGRID; i++) {
    for (int j = 0; j < NGRID; j++) {
      global_energy += compute_energy_cell(i, j, false);
    }
  }
}

/*============================================================================*/

// Computes the current global magnetization of the grid
void IsingModel::update_magnetization () {
  int i, j, sum;
  sum = 0;
  for (i = 0; i < NGRID; i++) {
    for (j = 0; j < NGRID; j++) {
      sum += grid[i][j];
    }
  }
  global_magnetization = sum/(double)(NCELLS);
}

/*============================================================================*/

// Computes the magnetization in a sample of the grid
void IsingModel::update_sample_magn (int sample) {
  int i, x, y, sum;
  sum = 0;
  for (i = 0; i < sample_size[sample]; i++) {
    getCellCoords(sample_cells[sample][i], x, y);
    sum += grid[x][y];
  }
  sample_magn[sample] = sum/(double)(sample_size[sample]);
}

/*============================================================================*/

// Randomize spins everywhere
void IsingModel::randomize () {
  int i, j, s;
  for (i = 0; i < NGRID; i++) {
    for (j = 0; j < NGRID; j++) {
      if (rand()%2 == 0){
        grid[i][j] = -1;
      } else {
        grid[i][j] = +1;
      }
    }
  }
  update_magnetization();
  for (s = 0; s < NUM_SAMPLES; s++) {
    update_sample_magn(s);
  }
}

/*============================================================================*/

// Sets spin to get a global magnetization close to the given value
// This will probablistically change the spins with a distribution that
// should yield a global magnetization close to the desired value.
void IsingModel::set_magnetization (double magn) {
  int i, j, s;
  double p = (magn+1)/2.0;
  for (i = 0; i < NGRID; i++) {
    for (j = 0; j < NGRID; j++) {
      if (rand()/(double)RAND_MAX<=p){
        grid[i][j] = +1;
      } else {
        grid[i][j] = -1;
      }
    }
  }
  update_magnetization();
  for (s = 0; s < NUM_SAMPLES; s++) {
    update_sample_magn(s);
  }
}

/*============================================================================*/

// ASCII representation of the current state of the grid
void IsingModel::display () {
  int i, j;
  for (i = 0; i < NGRID; i++) {
    for (j = 0; j < NGRID; j++) {
      if (grid[i][j]==+1) {
        printf("+");
      } else {
        printf("-");
      }
      if (j<(NGRID-1)) printf(" ");
    }
    printf ("\n");
  }
}

/*============================================================================*/

// Activates dead cells
// This will allocate the dead_cells array (if not already allocated) and
// turn on the useDeadCells flag
void IsingModel::activateDeadCells() {

  int i, j;

  // Allocate dead_cells array if not allocated
  if (!dead_cells) {
    dead_cells = (bool**) malloc(NGRID*sizeof(bool*));
    for (i = 0; i < NGRID; i++) {
      dead_cells[i] = (bool*) malloc(NGRID*sizeof(bool));
    }
  }

  // Initialize to all false
  for (i = 0; i < NGRID; i++) {
    for (j = 0; j < NGRID; j++) {
       dead_cells[i][j] = false;
    }
  }

  useDeadCells = true;

}

/*============================================================================*/

// Randomizes dead cells
// This will turn cells dead at random with probability density, which must
// be a number in the range [0,1]
void IsingModel::randomizeDead(double density) {

  int i, j;
  for (i = 0; i < NGRID; i++) {
    for (j = 0; j < NGRID; j++) {
      if (rand()/(double)RAND_MAX <= density){
        dead_cells[i][j] = true;
      } else {
        dead_cells[i][j] = false;
      }
    }
  }

}
/*============================================================================*/

// Advances the grid by one "generation"
// One generation is defined as having attempted a flip for *all* cells.
// The order in which the cells are temptatively flipped depends on the
// defined flip strategy, which is one of the following:
// STRATEGY_RANDOM: flips are done completely at random, stopping after NCELLS
//                  flips have been attempted.
// STRATEGY_SHUFFLED: before each generation, the order of the flips is shuffled
//                    (by means of the Fisher-Yates shuffle).
// STRATEGY_SEQUENTIAL: the flips are done sequentially from left to right and
//                      from top to bottom
// STRATEGY_PEANO: does a sequential flip following two converging Peano-like
//                 curves to try to minimize direction bias
// STRATEGY_COPY: the grid is copied, and the flips are done sequentially but
//                using the copied grid for neighbor information.
// Note that in all strategies except the last no copy of the grid is made,
// so later flips may depend on the results of previous ones.
void IsingModel::doGeneration () {

  int i, j, x, y, tmp;
  int i1, j1, i2, j2, count, next, d1, d2;

  switch (flip_strategy) {

  case STRATEGY_SHUFFLE:

    // Shuffle flip order
    for (i = 0; i < NCELLS; i++) {
      x = rand() % (NCELLS-i) + i;
      tmp = flip_order[i];
      flip_order[i] = flip_order[x];
      flip_order[x] = tmp;
    }
    // Attempt flip for all cells
    for (i = 0; i < NCELLS; i++) {
      getCellCoords(flip_order[i], x, y);
      tryCellFlip(x,y,false);
    }
    break;

  case STRATEGY_RANDOM:

    // Completely random flips. Stops after NCELLS flips.
    for (tmp = 1; tmp <= NCELLS; tmp++) {
      i = rand() % NGRID;
      j = rand() % NGRID;
      tryCellFlip(i,j,false);
    }
    break;

  case STRATEGY_SEQUENTIAL:

    // Do flips in sequential (index) order
    for (i = 0; i < NGRID; i++) {
      for (j = 0; j < NGRID; j++) {
        tryCellFlip(i,j,false);
      }
    }
    break;

  case STRATEGY_PEANO:

    // The flips are done sequentially following two Peano-like curves,
    // one that starts at (0,0) corner and the other in the opposite corner.
    // This is done to attenuate the directional bias of the normal sequential
    // order.
    i1 = 0;
    j1 = 0;
    i2 = NGRID-1;
    j2 = NGRID-1;
    count = 0;
    next = 0;
    d1 = +1;
    d2 = -1;
    while (count < NCELLS) {
      if (next == 0) {
        if ((d1 == +1 && i1 == NGRID-1) || (d1 == -1 && i1 == 0)){
          j1 += 1;
          d1 *= -1;
        } else {
          i1 += d1;
        }
        tryCellFlip(i1,j1,false);
        next = 1;
      } else if (next == 1) {
        if ((d2 == +1 && i2 == NGRID-1) || (d2 == -1 && i2 == 0)){
          j2 -= 1;
          d2 *= -1;
        } else {
          i2 += d2;
        }
        tryCellFlip(i2,j2,false);
        next = 0;
      }
      count++;
    }
    break;

  case STRATEGY_COPY:

    // In this strategy the current grid state is first copied. Then, the flips
    // are done sequentially but using the unchanging copied state to determine
    // neighboring spins. This eliminates the possible dependence on previous
    // flips in the same generation, effectively doing all flips simultaneously

    // Allocate grid_copy array if not allocated
    if (!grid_copy) {
      grid_copy = (int**) malloc(NGRID*sizeof(int*));
      for (i = 0; i < NGRID; i++) {
        grid_copy[i] = (int*) malloc(NGRID*sizeof(int));
      }
    }

    // Copy grid
    for (i = 0; i < NGRID; i++) {
      for (j = 0; j < NGRID; j++) {
        grid_copy[i][j] = grid[i][j];
      }
    }

    // Do flips (energy computed using grid copy)
    for (i = 0; i < NGRID; i++) {
      for (j = 0; j < NGRID; j++) {
        tryCellFlip(i,j,true);
      }
    }
    break;

  }

  // Update stats (and sample stats, if applicable)
  cur_gen++;
  if (cur_gen>=START_GEN) {
    update_stats();
    if (track_samples) update_sample_stats();
  }

}

/*============================================================================*/

// Returns the energy of a cell
// The grid wraps around at the edges (toroidal symmetry)
// The from_copy boolean determines if the neighbor information is pulled from
// the current state of the grid or from a copy of the previous generation's
// grid.
int IsingModel::compute_energy_cell (int i, int j, bool from_copy) {

  int ip, im, jp, jm;
  int energy, neigh_sum;
  int** _grid;

  ip = i+1;
  im = i-1;
  jp = j+1;
  jm = j-1;
  if      (i == 0) im = NGRID-1;
  else if (i == NGRID-1) ip = 0;
  if      (j == 0) jm = NGRID-1;
  else if (j == NGRID-1) jp = 0;

  if (from_copy) {
    _grid = grid_copy;
  } else {
    _grid = grid;
  }

  neigh_sum = 0;
  if (!useDeadCells || !dead_cells[ip][j])
    neigh_sum += _grid[ip][j];
  if (!useDeadCells || !dead_cells[im][j])
    neigh_sum += _grid[im][j];
  if (!useDeadCells || !dead_cells[i][jp])
    neigh_sum += _grid[i][jp];
  if (!useDeadCells || !dead_cells[i][jm])
    neigh_sum += _grid[i][jm];
  energy = -_grid[i][j] * neigh_sum;
  
  return energy;

}

/*============================================================================*/

// Attempts to flip cell (i,j)
// This will compute the change in energy the flip would produce
// and then apply the Metropolis algorithm:
// 1) if the new energy is lower or unchanged, always flip
// 2) if the new energy is higher, flip with probability e^(-deltaE/T)
// The from_copy boolean determines if the neighbor information is pulled from
// the current state of the grid or from a copy of the previous generation's
// grid.
void IsingModel::tryCellFlip (int i, int j, bool from_copy) {

  int old_E, new_E, deltaE;
  int ID, s;
  double prob;
  bool do_flip;
  int ip, im, jp, jm;

  old_E = compute_energy_cell(i, j, from_copy);
  new_E = -old_E;   // Always true since E_i = s_i*(sum_neighs s_n)
  deltaE = new_E - old_E;

  if (trans_dynamics == DYNAMICS_METROPOLIS) {

    // Metropolis dynamics
    if (new_E <= old_E) {
      // New energy is lower or unchanged: always flip cell
      prob = 1.0;
    } else {
      // New energy is higher: try thermal flip
      prob = exp(-deltaE/TEMP);
    }

  } else if (trans_dynamics == DYNAMICS_GLAUBER) {

    // Glauber dynamics
    prob = 1/(1 + exp(deltaE/TEMP));

  }

  // Roll the "die"
  if ((rand()/((double)RAND_MAX)) <= prob) {
    do_flip = true;
  } else {
    do_flip = false;
  }

  if (do_flip) {

    // Flip cell
    grid[i][j] *= -1;

    // Update global energy and magnetization
    global_magnetization += grid[i][j]*2/(double)(NCELLS);
    global_energy += deltaE;

    // Update magnetization of sample if cell in list
    if (track_samples) {
      getCellID(i, j, ID);
      for (s = 0; s < NUM_SAMPLES; s++) {
        if (inSample(ID,s)) {
          sample_magn[s] += grid[i][j]*2/(double)(sample_size[s]);
        }
      }
    }

  }

}

/*============================================================================*/

// Determines whether a cell is in the cell list of sample number s
// Assumes the cell list is sorted!
bool IsingModel::inSample(int ID, int s) {
  int idx;
  idx = binary_search(sample_cells[s], sample_size[s], ID);
  if (idx != -1) return true;
  else return false;
}

/*============================================================================*/

// Updates the global mean and variance of the magnetization
// Uses the current value of the magnetization as the "new" data point
void IsingModel::update_stats () {
  double delta;
  global_npoints++;
  delta = global_magnetization - global_mean;
  global_mean = global_mean + delta/global_npoints;
  global_M2 = global_M2 + delta*(global_magnetization - global_mean);
  if (global_npoints==1) {
    global_variance = 0.0;
  } else {
    global_variance = global_M2/(global_npoints-1);
  }
}

/*============================================================================*/

// Updates the sample mean and variance of the sample magnetization
// Uses the current value of the sample magnetization as the "new" data point
void IsingModel::update_sample_stats () {
  double delta;
  for (int s = 0; s < NUM_SAMPLES; s++) {
    sample_npts[s]++;
    delta = sample_magn[s] - sample_mean[s];
    sample_mean[s] = sample_mean[s] + delta/sample_npts[s];
    sample_M2[s] = sample_M2[s] + delta*(sample_magn[s] - sample_mean[s]);
    if (sample_npts[s]==1) {
      sample_var[s] = 0.0;
    } else {
      sample_var[s] = sample_M2[s]/(sample_npts[s]-1);
    }
  }
}

/*============================================================================*/

// Remembers last NUM_DATA magnetization values for running mean/variance
void IsingModel::update_data () {
  rundata[nextdata] = global_magnetization;
  if (nextdata>=NUM_DATA) {
    nextdata = 0;
  } else {
    nextdata++;
  }
}

/*============================================================================*/

// Computes the running mean and variance using the stored values
void IsingModel::running_stats () {
  int i;
  run_mean = 0.0;
  run_var = 0.0;
  for (i = 0; i < NUM_DATA; i++) {
    run_mean += rundata[i];
  }
  run_mean /= (double)NUM_DATA;
  for (i = 0; i < NUM_DATA; i++) {
    run_var += (rundata[i]-run_mean)*(rundata[i]-run_mean);
  }
  run_var /= (double)(NUM_DATA-1);
}

/*============================================================================*/

// Randomly selects samples of cells for statistical followup
// The arrays sample_cells are filled with cell IDs, and then each is sorted
// for quick cell lookup.
void IsingModel::pickSamples () {
  int i, s, x;
  int* pool = (int*) malloc(NCELLS*sizeof(int));
  for (s = 0; s < NUM_SAMPLES; s++) {
    for (i = 0; i < NCELLS; i++) {
      pool[i] = i;
    }
    for (i = 0; i < sample_size[s]; i++) {
      x = rand() % (NCELLS-i) + i;
      sample_cells[s][i] = pool[x];
      pool[x] = pool[i];
    }
    quicksort(sample_cells[s], sample_size[s]);
  }
}

/*============================================================================*/

// Converts a cell ID to an (x,y) position
// Coords and IDs start at zero.
void IsingModel::getCellCoords(int ID, int &x, int &y) {
  x = ID % NGRID;
  y = ID / NGRID;
}

/*============================================================================*/

// Converts an (x,y) position into a cell ID
// Coords and IDs start at zero.
void IsingModel::getCellID(int x, int y, int &ID) {
  ID = x + y*NGRID;
}

/*============================================================================*/
