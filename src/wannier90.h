#ifndef WANNIER90_H
#define WANNIER90_H

#include <ISO_Fortran_binding.h>
#include <string.h>

/* CFI_cdesc_t provides a mechanism for passing variable length char* arrays */
// see https://community.intel.com/t5/Intel-Fortran-Compiler/C-interoperablilty-and-character-strings/td-p/1084167

/* the internals of the fortran data structure are not visible from C
 * provide a type wrapping a void* pointer to memory allocated by fortran
 */
typedef struct {
  void* cptr = NULL;
} w90_data;

void w90_set_option_double2d_c(w90_data, CFI_cdesc_t*, void*, int, int);
void w90_set_option_double2d(w90_data blob, char* key, void* data, int x, int y) {
        CFI_cdesc_t desc;
        CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
        w90_set_option_double2d_c(blob, &desc, data, x, y);
};
void w90_set_option_int2d_c(w90_data, CFI_cdesc_t*, void*, int, int);
void w90_set_option_int2d(w90_data blob, char* key, void* data, int x, int y) {
        CFI_cdesc_t desc;
        CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
        w90_set_option_int2d_c(blob, &desc, data, x, y);
};
void w90_set_option_int1d_c(w90_data, CFI_cdesc_t*, void*, int);
void w90_set_option_int1d(w90_data blob, char* key, void* data, int x) {
        CFI_cdesc_t desc;
        CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
        w90_set_option_int1d_c(blob, &desc, data, x);
};
void w90_set_option_int_c(w90_data, CFI_cdesc_t*, int);
void w90_set_option_int(w90_data blob, char* key, int data) {
        CFI_cdesc_t desc;
        CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
        w90_set_option_int_c(blob, &desc, data);
};
void w90_input_setopt_c(w90_data blob, CFI_cdesc_t*, int&, int);
void w90_input_setopt(w90_data blob, char* seed, int& ierr) {
        CFI_cdesc_t desc;
        CFI_establish(&desc, seed, CFI_attribute_other, CFI_type_char, strlen(seed), 0, NULL);
        ierr = 0;
        w90_input_setopt_c(blob, &desc, ierr, 0);
}
void w90_create(w90_data&);
void w90_delete(w90_data&);
void w90_disentangle_c(w90_data, int&);
void w90_get_centres_c(w90_data, void*);
void w90_get_gkpb_c(w90_data, void*);
void w90_get_nn_c(w90_data, void*);
void w90_get_nnkp_c(w90_data, void*);
void w90_get_spreads_c(w90_data, void*);
void w90_input_reader_c(w90_data, int&);
void w90_overlaps_c(w90_data, int&);
void w90_project_overlap_c(w90_data, int&);
void w90_set_eigval_c(w90_data, void*);
void w90_set_m_local_c(w90_data, void*);
void w90_set_u_matrix_c(w90_data, void*);
void w90_set_u_opt_c(w90_data, void*);
void w90_wannierise_c(w90_data, int&);
#ifdef MPI_VERSION
#include <mpi.h>
void w90_set_comm(void*, int);
void w90_set_comm_c(w90_data blob, MPI_Comm comm){
        int fcomm = MPI_Comm_c2f(comm);
        w90_set_comm(blob%cptr, fcomm);
}
#endif
#endif
