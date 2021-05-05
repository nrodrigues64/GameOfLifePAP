
#include "easypap.h"
#include "rle_lexer.h"

#include <omp.h>
#include <stdbool.h>
#include <string.h>
#include <sys/mman.h>
#include <stdint.h>
#include <unistd.h>

static unsigned color = 0xFFFF00FF; // Living cells have the yellow color

typedef unsigned cell_t;

int iter = 1;
int *tmp;
static cell_t *restrict _table = NULL, *restrict _alternate_table = NULL;

static inline cell_t *table_cell (cell_t *restrict i, int y, int x)
{
  return i + y * DIM + x;
}

// This kernel does not directly work on cur_img/next_img.
// Instead, we use 2D arrays of boolean values, not colors
#define cur_table(y, x) (*table_cell (_table, (y), (x)))
#define next_table(y, x) (*table_cell (_alternate_table, (y), (x)))

void life_init (void)
{
  // life_init may be (indirectly) called several times so we check if data were
  // already allocated
  if (_table == NULL) {
    const unsigned size = DIM * DIM * sizeof (cell_t);

    PRINT_DEBUG ('u', "Memory footprint = 2 x %d bytes\n", size);

    _table = mmap (NULL, size, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    _alternate_table = mmap (NULL, size, PROT_READ | PROT_WRITE,
                             MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  }
}

void life_finalize (void)
{
  const unsigned size = DIM * DIM * sizeof (cell_t);

  munmap (_table, size);
  munmap (_alternate_table, size);
}

// This function is called whenever the graphical window needs to be refreshed
void life_refresh_img (void)
{
  for (int i = 0; i < DIM; i++)
    for (int j = 0; j < DIM; j++)
      cur_img (i, j) = cur_table (i, j) * color;
}

static inline void swap_tables (void)
{
  cell_t *tmp = _table;

  _table           = _alternate_table;
  _alternate_table = tmp;
}

///////////////////////////// Sequential version (seq)
static int compute_new_state (int y, int x)
{
  unsigned n      = 0;
  unsigned me     = cur_table (y, x) != 0;
  unsigned change = 0;

  if (x > 0 && x < DIM - 1 && y > 0 && y < DIM - 1) {

    for (int i = y - 1; i <= y + 1; i++)
      for (int j = x - 1; j <= x + 1; j++)
        n += cur_table (i, j);

    n = (n == 3 + me) | (n == 3);
    if (n != me)
      change |= 1;

    next_table (y, x) = n;
  }

  return change;
}


unsigned life_compute_seq (unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++) {
    int change = 0;

    monitoring_start_tile (0);

    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        change |= compute_new_state (i, j);

    monitoring_end_tile (0, 0, DIM, DIM, 0);

    swap_tables ();

    if (!change)
      return it;
  }

  return 0;
}

///////////////////////////// Tiled sequential version (tiled)

// Tile inner computation
static int do_tile_reg (int x, int y, int width, int height)
{
  int change = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
      change |= compute_new_state (i, j);

  return change;
}

static int do_tile (int x, int y, int width, int height, int who)
{
  int r;

  monitoring_start_tile (who);

  r = do_tile_reg (x, y, width, height);

  monitoring_end_tile (x, y, width, height, who);
  
  return r;
  
}

unsigned life_compute_tiled (unsigned nb_iter)
{
  unsigned res = 0;

  for (unsigned it = 1; it <= nb_iter; it++) {
    unsigned change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        change |= do_tile (x, y, TILE_W, TILE_H, omp_get_thread_num());

    swap_tables ();

    if (!change) { // we stop when all cells are stable
      res = it;
      break;
    }
  }

  return res;
}




///////////////////////////// Tiled omp version (tiled)
unsigned life_compute_omp (unsigned nb_iter)
{
  
  unsigned res = 0;
  
  for (unsigned it = 1; it <= nb_iter; it++) {
    unsigned change = 0;
  #pragma omp parallel for schedule(static)
    for (int y = 0; y < DIM; y += TILE_H) 
      for (int x = 0; x < DIM; x += TILE_W)
        #pragma omp critical 
        {
        change |= do_tile (x, y, TILE_W, TILE_H, omp_get_thread_num());
        }
      
    
    swap_tables ();
    //printf("change : %d\n", change);
    if (!change) { // we stop when all cells are stable
      res = it;
      break;
    }
  }

  return res;
}
///////////////////////////// Initial configs

void life_draw_guns (void);

static inline void set_cell (int y, int x)
{
  cur_table (y, x) = 1;
  if (opencl_used)
    cur_img (y, x) = 1;
}

static inline int get_cell (int y, int x)
{
  return cur_table (y, x);
}

static void inline life_rle_parse (char *filename, int x, int y,
                                   int orientation)
{
  rle_lexer_parse (filename, x, y, set_cell, orientation);
}

static void inline life_rle_generate (char *filename, int x, int y, int width,
                                      int height)
{
  rle_generate (x, y, width, height, get_cell, filename);
}

void life_draw (char *param)
{
  if (param && (access (param, R_OK) != -1)) {
    // The parameter is a filename, so we guess it's a RLE-encoded file
    life_rle_parse (param, 1, 1, RLE_ORIENTATION_NORMAL);
  } else
    // Call function ${kernel}_draw_${param}, or default function (second
    // parameter) if symbol not found
    hooks_draw_helper (param, life_draw_guns);
}

static void otca_autoswitch (char *name, int x, int y)
{
  life_rle_parse (name, x, y, RLE_ORIENTATION_NORMAL);
  life_rle_parse ("data/rle/autoswitch-ctrl.rle", x + 123, y + 1396,
                  RLE_ORIENTATION_NORMAL);
}

static void otca_life (char *name, int x, int y)
{
  life_rle_parse (name, x, y, RLE_ORIENTATION_NORMAL);
  life_rle_parse ("data/rle/b3-s23-ctrl.rle", x + 123, y + 1396,
                  RLE_ORIENTATION_NORMAL);
}

static void at_the_four_corners (char *filename, int distance)
{
  life_rle_parse (filename, distance, distance, RLE_ORIENTATION_NORMAL);
  life_rle_parse (filename, distance, distance, RLE_ORIENTATION_HINVERT);
  life_rle_parse (filename, distance, distance, RLE_ORIENTATION_VINVERT);
  life_rle_parse (filename, distance, distance,
                  RLE_ORIENTATION_HINVERT | RLE_ORIENTATION_VINVERT);
}

// Suggested cmdline: ./run -k life -s 2176 -a otca_off -ts 64 -r 10 -si
void life_draw_otca_off (void)
{
  if (DIM < 2176)
    exit_with_error ("DIM should be at least %d", 2176);

  otca_autoswitch ("data/rle/otca-off.rle", 1, 1);
}

// Suggested cmdline: ./run -k life -s 2176 -a otca_on -ts 64 -r 10 -si
void life_draw_otca_on (void)
{
  if (DIM < 2176)
    exit_with_error ("DIM should be at least %d", 2176);

  otca_autoswitch ("data/rle/otca-on.rle", 1, 1);
}

// Suggested cmdline: ./run -k life -s 6208 -a meta3x3 -ts 64 -r 50 -si
void life_draw_meta3x3 (void)
{
  if (DIM < 6208)
    exit_with_error ("DIM should be at least %d", 6208);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      otca_life (j == 1 ? "data/rle/otca-on.rle" : "data/rle/otca-off.rle",
                 1 + j * (2058 - 10), 1 + i * (2058 - 10));
}

// Suggested cmdline: ./run -k life -a bugs -ts 64
void life_draw_bugs (void)
{
  for (int y = 16; y < DIM / 2; y += 32) {
    life_rle_parse ("data/rle/tagalong.rle", y + 1, y + 8,
                    RLE_ORIENTATION_NORMAL);
    life_rle_parse ("data/rle/tagalong.rle", y + 1, (DIM - 32 - y) + 8,
                    RLE_ORIENTATION_NORMAL);
  }
}

// Suggested cmdline: ./run -k life -v omp -a ship -s 512 -m -ts 16
void life_draw_ship (void)
{
  for (int y = 16; y < DIM / 2; y += 32) {
    life_rle_parse ("data/rle/tagalong.rle", y + 1, y + 8,
                    RLE_ORIENTATION_NORMAL);
    life_rle_parse ("data/rle/tagalong.rle", y + 1, (DIM - 32 - y) + 8,
                    RLE_ORIENTATION_NORMAL);
  }

  for (int y = 43; y < DIM - 134; y += 148) {
    life_rle_parse ("data/rle/greyship.rle", DIM - 100, y,
                    RLE_ORIENTATION_NORMAL);
  }
}

void life_draw_stable (void)
{
  for (int i = 1; i < DIM - 2; i += 4)
    for (int j = 1; j < DIM - 2; j += 4) {
      set_cell (i, j);
      set_cell (i, j + 1);
      set_cell (i + 1, j);
      set_cell (i + 1, j + 1);
    }
}

void life_draw_guns (void)
{
  at_the_four_corners ("data/rle/gun.rle", 1);
}

void life_draw_random (void)
{
  for (int i = 1; i < DIM - 1; i++)
    for (int j = 1; j < DIM - 1; j++)
      if (random () & 1)
        set_cell (i, j);
}

// Suggested cmdline: ./run -k life -a clown -s 256 -i 110
void life_draw_clown (void)
{
  life_rle_parse ("data/rle/clown-seed.rle", DIM / 2, DIM / 2,
                  RLE_ORIENTATION_NORMAL);
}

void life_draw_diehard (void)
{
  life_rle_parse ("data/rle/diehard.rle", DIM / 2, DIM / 2,
                  RLE_ORIENTATION_NORMAL);
}
#include <immintrin.h>
#if defined(ENABLE_VECTO) && (AVX2 == 1)


void display_vec (__m256i v)
{
  char tmp [32];

  _mm256_store_si256 ((__m256i *)tmp, v);

  for (int i = 0; i < 32; i++){
    printf ("%u ", tmp [i]);
    if(i%32==31)
      printf("\n");
  }
}


/*static int compute_new_state (int y, int x)
{
  unsigned n      = 0;
  unsigned me     = cur_table (y, x) != 0;
  unsigned change = 0;

  if (x > 0 && x < DIM - 1 && y > 0 && y < DIM - 1) {

    for (int i = y - 1; i <= y + 1; i++)
      for (int j = x - 1; j <= x + 1; j++)
        n += cur_table (i, j);

    n = (n == 3 + me) | (n == 3);
    if (n != me)
      change |= 1;

    next_table (y, x) = n;
  }

  return change;
}*/

//./run -k life -v vec -s 32 -a stable -i 1 -n
static int compute_new_state_vec (int y, int x,int height, int width)
{
  const __m256i zero = _mm256_set1_epi8(0);
  __m256i change = _mm256_set1_epi8(0);
  char x_c[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  
  const __m256i dim = _mm256_set1_epi8(DIM);

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j+=32) {
      //Pour les x
      
      __m256i tmp_x = _mm256_set1_epi8( x );
      __m256i x_vec = _mm256_load_si256((__m256i *)x_c);
      x_vec = _mm256_add_epi8(tmp_x,x_vec); 
      
      
      __m256i y_vec = _mm256_set1_epi8( i );
      
      
      

      __m256i mask_x_more_than_zero = _mm256_cmpgt_epi8(x_vec,zero);
      __m256i mask_x_less_than_DIM = _mm256_cmpgt_epi8(dim,x_vec);
      

      __m256i mask_y_more_than_zero = _mm256_cmpgt_epi8(y_vec,zero);
      __m256i mask_y_less_than_DIM = _mm256_cmpgt_epi8(dim,y_vec);
      
      


      __m256i x_y = _mm256_loadu_si256((__m256i *)&cur_table(y,x)); 
      __m256i mask_me = _mm256_cmpeq_epi8(x_y,zero);
      __m256i me_if = _mm256_set1_epi8(0);
      __m256i me_else = _mm256_set1_epi8(1);
      __m256i me = _mm256_blendv_epi8(me_else,me_if,mask_me);

      __m256i x_y_right = _mm256_set1_epi8(0);
      //x_y_right = _mm256_blendv_epi8(zero,(__m256i *)&cur_table(y,x+1));
    

      __m256i x_y_left = _mm256_set1_epi8(0);
      x_y_left = _mm256_blendv_epi8(zero,_mm256_loadu_si256((__m256i *)&cur_table(y,x-1)),mask_x_more_than_zero);
      printf("x_y_left :");
      display_vec(zero);
      display_vec(_mm256_loadu_si256((__m256i *)&cur_table(y,x)));
      display_vec(x_y_left);

      /*__m256i mask_x_y_up = _mm256_cmpgt_epi8(_mm256_load_si256(&y_vec),_mm256_set1_epi8(0));
      
      __m256i x_y_up = _mm256_blendv_epi8(_mm256_set1_epi8(0),_mm256_loadu_si256((__m256i *)&cur_table(y-1,x)),mask_x_y_up);*/
      __m256i x_y_up = _mm256_set1_epi8(0);
      if(y == 0) {
        x_y_up = _mm256_set1_epi8(0);
      } else {
        x_y_up = _mm256_loadu_si256((__m256i *)&cur_table(y-1,x));
      }

      __m256i x_y_down = _mm256_set1_epi8(0);
      if(y == DIM-1) {
        x_y_down = _mm256_set1_epi8(0);
      } else {
        x_y_down = _mm256_loadu_si256((__m256i *)&cur_table(y+1,x));
      }
      
      __m256i x_y_up_right = _mm256_set1_epi8(0);
      if(y == 0 || x+32 == DIM-1) {
        x_y_up_right = _mm256_set1_epi8(0);
      } else {
        x_y_up_right = _mm256_loadu_si256((__m256i *)&cur_table(y-1,x+1));
      }

      __m256i x_y_up_left = _mm256_set1_epi8(0);
      if(y == 0 || x == 0) {
        x_y_up_left = _mm256_set1_epi8(0);
      } else {
        x_y_up_left = _mm256_loadu_si256((__m256i *)&cur_table(y-1,x-1));
      }

      __m256i x_y_down_right =  _mm256_set1_epi8(0);
      if(y == DIM-1 || x+32 == DIM-1) {
        x_y_down_right =  _mm256_set1_epi8(0);
      } else {
        x_y_down_right = _mm256_loadu_si256((__m256i *)&cur_table(y+1,x+1));
      }
      
      __m256i x_y_down_left = _mm256_set1_epi8(0);
      if(y == DIM-1 || x == 0) {
        x_y_down_left = _mm256_set1_epi8(0);
      } else {
        x_y_down_left = _mm256_loadu_si256((__m256i *)&cur_table(y+1,x-1));
      }
      
      __m256i n = _mm256_add_epi8(x_y,x_y_right);
      
      n = _mm256_add_epi8(n, x_y_left);
      
      n = _mm256_add_epi8(n, x_y_up);
      
      n = _mm256_add_epi8(n, x_y_down);
      
      n = _mm256_add_epi8(n, x_y_up_right);
      n = _mm256_add_epi8(n, x_y_up_left);
      n = _mm256_add_epi8(n, x_y_down_right);
      n = _mm256_add_epi8(n, x_y_down_left);
      
      __m256i me_3 = _mm256_add_epi8(me,_mm256_set1_epi8(3));
      __m256i mask_me_3 = _mm256_cmpeq_epi8(n, me_3); 
      __m256i mask_3 = _mm256_cmpeq_epi8(n, _mm256_set1_epi8(3));
      __m256i mask_n = _mm256_add_epi8(mask_me_3,mask_3);
      __m256i n_if = _mm256_set1_epi8(1);
      __m256i n_else = _mm256_set1_epi8(0);
      n = _mm256_blendv_epi8(n_else,n_if,mask_n);
      
      __m256i mask_me_n = _mm256_cmpeq_epi8(n,me);
      change = _mm256_set1_epi8(0);
      __m256i change_else = _mm256_set1_epi8(1);
      change = _mm256_blendv_epi8(change_else,change,mask_me_n); 

    }
  
  return !_mm256_testz_si256(change,_mm256_set1_epi8(1));
  
  //
  
  
  //int * a = (int*) &change;
  //return 0;
  
  
}

static int do_tile_vec (int x, int y, int width, int height)
{
  
  __m256i change = _mm256_set1_epi8(0);
  __m256i one = _mm256_set1_epi8(1); 
  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++) {
      __m256i tmp = _mm256_set1_epi8(compute_new_state (i, j));
      __m256i mask_change = _mm256_cmpeq_epi8(tmp,one);
      change = _mm256_blendv_epi8(change,tmp,mask_change);
    }
  int * a = (int*)&change;
  return a[0];
}

//////////////////////////// Lazy vect version (tiled)
// Tile inner computation

static int do_tile_rec (int x, int y, int width, int height)
{
  int change = 0;
  
  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
      change |= compute_new_state (i, j);

  return change;
}

static int do_tile_v (int x, int y, int width, int height, int who)
{
  int r;

  monitoring_start_tile (who);

  r = compute_new_state_vec (x, y, width, height);

  monitoring_end_tile (x, y, width, height, who);
  
  return r;
  
}

unsigned life_compute_lazy_vec (unsigned nb_iter)
{
  compute_new_state_vec(0,0,32,32);
  unsigned res = 0;
  //compute_new_state_vec(0,0);
  for (unsigned it = 1; it <= nb_iter; it++) {
    unsigned change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W) {
        change |= do_tile (x, y, TILE_W, TILE_H, 0);
        
      }
        
    swap_tables ();

    if (!change) { // we stop when all cells are stable
      res = it;
      break;
    }
  }
  return res;
}

unsigned life_compute_lazy_vec_omp (unsigned nb_iter)
{
  unsigned res = 0;

  for (unsigned it = 1; it <= nb_iter; it++) {
    unsigned change = 0;
  #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W) {
        change |= do_tile_v (x, y, TILE_W, TILE_H, omp_get_thread_num());
      }
        
    swap_tables ();

    if (!change) { // we stop when all cells are stable
      res = it;
      break;
    }
  }
  return res;
}

#endif