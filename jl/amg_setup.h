#ifndef AMG_SETUP
#define AMG_SETUP

/* Main function build a structure crs_data required to solve AMG */
void amg_setup(uint n, const ulong *id, uint nz_unassembled, const uint *Ai, 
    const uint* Aj, const double *Av,  struct crs_data *data);

/*******************************************************************************
* AMG functions
*******************************************************************************/
void coarsen(double *vc, struct csr_mat *A, double ctol);

void interpolation(struct csr_mat *W, struct csr_mat *Af, struct csr_mat *Ac, 
    struct csr_mat *Ar, double gamma2, double tol);

/*******************************************************************************
* Algebraic functions
*******************************************************************************/
static void mat_max(double *y, struct csr_mat *A, double *f, double *x, 
    double tol);

uint lanczos(double **lambda, struct csr_mat *A);

static void tdeig(double *lambda, double *y, double *d, const double *v,
                  const int n);
static double sec_root(double *y, const double *d, const double *v,
                       const int ri, const int n);
static double rat_root(const double a, const double b, const double c,
                       const double sign);
static double sum_3(const double a, const double b, const double c);

void chebsim(double *m, double *c, double rho, double tol);

void sparsify(double *S, struct csr_mat *A, double tol);

uint pcg(double *v, struct csr_mat *A, double *r, double *M, double tol);

void min_skel(struct csr_mat *W_skel, struct csr_mat *R);

void solve_weights(struct csr_mat *W, struct csr_mat *W0, double *lam, 
    struct csr_mat *W_skel, struct csr_mat *Af, struct csr_mat *Ar, uint rnc,
    double *alpha, double *u, double *v, double tol);

void interp(struct csr_mat *X, struct csr_mat *A, struct csr_mat *B, double *u, 
    double *lambda);

static void mv_utt(double *y, uint n, const double *U, const double *x);

static void mv_ut(double *y, uint n, const double *U, const double *x);

static void sp_restrict_unsorted(double *y, uint yn, const uint *map_to_y,
    uint xn, const uint *xi, const double *x);

static void sp_restrict_sorted(double *y, uint Rn, const uint *Ri, uint xn, 
    const uint *xi, const double *x);

/*******************************************************************************
* Exctract a sub-matrix in csr format
*******************************************************************************/
// subA = A(vr, vc), where vr and vc are vectors of 0 and 1 
// It is assumed that vr has (at least) size A->rn and vc has size A->cn
void sub_mat(struct csr_mat *subA, struct csr_mat *A, double* vr, double *vc);

/*******************************************************************************
* Exctract a sub-vector
*******************************************************************************/
// a = b(v), where v is vector made of 0 and 1 
void sub_vec(double *a, double *b, double* v, uint n);
// idem but for slong type (used for gs_id)
void sub_slong(slong *a, slong *b, double* v, uint n);

/*******************************************************************************
* Vector-vector operations
*******************************************************************************/
// a[i] = a[i] (op) b[i] for i = 0,n-1
enum vv_ops {plus, minus, ewmult}; //+, -, element-wise multiplication (.* in 
                                   // Matlab)
void vv_op(double *a, double *b, uint n, enum vv_ops op);

// Dot product between two vectors
double vv_dot(double *a, double *b, uint n);

/*******************************************************************************
* Binary operations
*******************************************************************************/
// For 'and', 'or' and  'xor' operations:
//     mask[i] = 1 if (mask[i] (op) a[i]) is true
//             = 0 otherwise
//     i = 0, n-1
//
// For 'not' operation: 
//     mask[i] = not(a[i])
//     i = 0, n-1
enum bin_ops {and_op, or_op, xor_op, not_op};
void bin_op(double *mask, double *a, uint n, enum bin_ops op);

/*******************************************************************************
* Mask operations
*******************************************************************************/
// mask[i] = 1 if (a[i] (op) trigger) is true
//         = 0 otherwise
// i = 0, n-1
enum mask_ops {gt, lt, ge, le, eq}; //>, <, >=, <=, =
void mask_op(double *mask, double *a, uint n, double trigger, enum mask_ops op);

/*******************************************************************************
* Extremum operations
*******************************************************************************/
// *extr = op(a[i]) (op = max or min)
// *idx  = index of extr
// i = 0, n-1
enum extr_ops {max, min};
void extr_op(double *extr, uint *idx, double *a, uint n, enum extr_ops op);

/*******************************************************************************
* Array operations
*******************************************************************************/
// a[i] = op(a[i])
// i = 0, n-1
enum array_ops {abs_op, sqrt_op, minv_op, sqr_op, sum_op, norm2_op}; 
// absolute value, sqrt, multiplicative inverse, square, sum, 2-norm
// If op = absolute value, sqrt, multiplicative inverse, square:
//      function returns 0
// If op = sum, 2-norm: 
//      function returns corresponding value
double array_op(double *a, uint n, enum array_ops op);

// a[i] = v
// i = 0, n-1
void init_array(double *a, uint n, double v);

// Operations between an array and a scalar
// a[i] = a[i] op scal
// i = 0, n-1
enum ar_scal_ops {mult_op}; // a[i] = a[i] * scal
void ar_scal_op(double *a, double scal, uint n, enum ar_scal_ops op);

/*******************************************************************************
* Diagonal operations
*******************************************************************************/
// Extract diagonal: D = diag(A)
void diag(double *D, struct csr_mat *A);

// Operations between csr and diagonal matrices
// A = A (op) D
enum diagcsr_ops {dplus, dminus, dmult, multd}; // A = {A+D, A-D, D*A, A*D}
void diagcsr_op(struct csr_mat *A, double *D, enum diagcsr_ops op);

/*******************************************************************************
* Others
*******************************************************************************/
// Copy csr matrix B <-- A
void copy_csr(struct csr_mat *B, struct csr_mat *A);
void csr_free(struct csr_mat **A);

/*******************************************************************************
* Functions used to build sparse matrix and id array for gs
*******************************************************************************/
// Matrix under coordinate list format (used as intermediate step to build csr)
typedef struct {uint i, j; double v;} coo_mat; 
typedef struct {coo_mat coo_A; uint dest;} coo_mat_dest; 

/* Build matrix using csr format */
void build_csr(struct csr_mat *A, coo_mat *coo_A, uint nnz);

/* Sorting functions */
int comp_coo_v (const void * a, const void * b);
int comp_coo_ij (const void * a, const void * b);
int comp_coo_ji (const void * a, const void * b);
int comp_uint (const void * a, const void * b);
int comp_gs_id(const void * a, const void * b);
/* Unused functions
int comp_coo_i (const void * a, const void * b);
int comp_coo_j (const void * a, const void * b);
*/

/* Remove duplicates in a sorted list */
static uint remdup(uint *array, uint size);

/* Function to build sparse matrix and gs_id */
void build_setup_data(struct csr_mat *A, uint n, const ulong *id,  
    uint nz_unassembled, const uint *Ai, const uint* Aj, const double *Av,   
    struct crs_data *data);

#endif
