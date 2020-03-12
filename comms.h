/*=========================================================/
/  Comms header, PX425 2018 assignment 4.                  /
/=========================================================*/

/* Number of processors in the communicator */
extern int p;

/* Rank of the current MPI processor/task */
extern int my_rank;

/* Coordinates of my_rank within a Cartesian grid */
extern int my_rank_coords[2];

/* Neighbouring MPI processors left, right, down, up */
extern int my_rank_neighbours[4];


/* Function prototypes */
void comms_initialise(int argc, char **argv);
void comms_processor_map();
void comms_get_global_mag(double local_mag,double *global_mag);
void comms_halo_swaps();
void comms_get_global_grid();
void comms_finalise();

