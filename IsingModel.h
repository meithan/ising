#ifndef ISING_H
#define ISING_H

/*===============================\\
|| Ising Model class declaration ||
\\===============================*/

class IsingModel {

  public:

  /*==========================================================================*/

  /* STATIC CONSTANTS */

  // The critical temperature in units of J/k
  // Tc = 2/ln(1+sqrt(2))
  static const double TEMP_CRIT = 2.26918531421302;

  /* MEMBER VARIABLES */

  // Temperature
  // Given in units of J/k, where J is the coupling constant of the interaction
  // and k is Boltzmann's constant.
  double TEMP;

  // Size of grid (NGRID x NGRID)
  int NGRID;

  // Total number of cells, equal to NGRID*NGRID
  int NCELLS;

  // 2D grid for spin states
  // grid[NGRID][NGRID]
  int** grid;

  // Copy of the grid (not allocated if not needed)
  // grid_copy[NGRID][NGRID]
  int** grid_copy;

  // Dead cells
  bool** dead_cells;
  bool useDeadCells;
  double DEAD_DENS;

  // Flip strategy.
  // See the doGeneration class documentation for information on valid options.
  int flip_strategy;
  static const int STRATEGY_SHUFFLE = 0;
  static const int STRATEGY_RANDOM = 1;
  static const int STRATEGY_SEQUENTIAL = 2;
  static const int STRATEGY_PEANO = 3;
  static const int STRATEGY_COPY = 4;

  // Dynamics
  int trans_dynamics;
  static const int DYNAMICS_METROPOLIS = 0;
  static const int DYNAMICS_GLAUBER = 1;

  // List of cell IDs for randomized flipping order
  // flip_order[NCELLS]
  int* flip_order;

  // Current generation (will never reset)
  int cur_gen;

  // Global statistics
  int global_energy;
  double global_magnetization;
  double global_mean;
  double global_variance;
  double global_M2;
  int global_npoints;

  // Track statistics in samples?
  bool track_samples;

  // Number of samples to monitor (logarithmic spacing)
  int NUM_SAMPLES;

  // Minimum sample size
  int SAMPLE_MIN;

  // Maximum sample size
  int SAMPLE_MAX;

  // Generation in which to start recording stats
  int START_GEN;

  // Sample statistics -- optional
  // These all have size NUM_SAMPLES
  double* sample_magn;
  double* sample_mean;
  double* sample_var;
  double* sample_M2;
  int* sample_npts;

  // Array of sample sizes (number of cells in each sample)
  // sample_size[NUM_SAMPLES]
  int* sample_size;

  // List of cells to use for sample statistics
  // sample_cells[NUM_SAMPLES][<number of cells in this sample>]
  int** sample_cells;

  // Number of points to remember for running mean/variance
  int NUM_DATA;

  // Running mean/variance
  double* rundata;
  double run_mean, run_var;
  int nextdata;

  /*==========================================================================*/

  /* MEMBER FUNCTIONS */

  IsingModel(int, double);
  IsingModel(int, double, int, int, int, int, int);
  void common_constructor();
  int compute_energy_cell(int, int, bool);
  void randomize();
  void set_magnetization(double);
  void reset_stats();
  void display();
  void update_energy();
  void update_magnetization();
  void update_sample_magn(int);
  void activateDeadCells();
  void randomizeDead(double);
  void doGeneration();
  void tryCellFlip(int,int,bool);
  void update_stats();
  void update_sample_stats();
  void update_data();
  void running_stats();
  void getCellCoords(int, int&, int&);
  void getCellID(int, int, int&);
  bool inSample(int,int);
  void pickSamples();

};

#endif // ISING_H
