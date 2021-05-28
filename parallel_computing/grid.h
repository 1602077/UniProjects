/*=========================================================/
/  Grid header, PX425 2018 assignment 4.                   /
/=========================================================*/

/* Integer constants defining coordinates and directions */
extern const int x,y,left,right,down,up;

/* Size of local grid */
extern int grid_domain_size;
extern int grid_domain_start[2];
extern int grid_domain_end[2];

/* Local grid of spins */
extern int **grid_spin;

/* Local random field */
extern double **grid_field;

/* Global grid of spins */
extern int **global_grid_spin;

/* Grid halo */
extern int **grid_halo;

/* Subdomain info */
extern int grid_subdom_size;
extern int grid_subdom_start[4][2];
extern int grid_subdom_end[4][2];

/* This is a struct of four pointers which point to the */
/* up, down, left and right neigbour spins.            */
typedef struct {
  int *up;
  int *down;
  int *left;
  int *right;
} neighbour_ptr_t;

/* Pointers to neighbours */
extern neighbour_ptr_t **grid_neighbours;

/* Function prototypes */
void grid_initialise_local(int Ngrid,int p,int my_rank,int *my_rank_coords);
void grid_initialise_global(int Ngrid,int my_rank,int p);
void grid_shutdown_local();
void grid_shutdown_global(int Ngrid,int my_rank,int p);
void grid_set_neighbour_spins();
