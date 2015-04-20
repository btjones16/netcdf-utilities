#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <netcdf.h>
#include "log.h"
#include "netcdf_utilities.h"

void load_xy(const char * restrict filename, float ** x, float ** y, 
	     size_t * nnodes) {
  /* Open the input file. */
  int ncid, ncerr;
  if((ncerr = nc_open(filename, NC_NOWRITE, &ncid)) != NC_NOERR)
    printf("%s\n", nc_strerror(ncerr));
  /* Read the length of x and y. */
  nc_inq_dimlen_name(ncid, "node", nnodes);
  /* Allocate memory for x and y. */
  if((*x = malloc(*nnodes * sizeof(float))) == NULL)
    printf("Malloc error.\n");
  if((*y = malloc(*nnodes * sizeof(float))) == NULL)
    printf("Malloc error.\n");
  /* Read the x and y variables. */
  nc_get_var_float_name(ncid, "x", *x);
  nc_get_var_float_name(ncid, "y", *y);
  /* Close the input file. */
  nc_close(ncid);
}

void choose_nodes_to_keep(size_t n_nodes, float * x, float * y,
			  const float x_min, const float x_max,
			  const float y_min, const float y_max,
			  size_t * n_keeps, int ** keeps) {
  *n_keeps = 0;
  for(size_t i = 0; i < n_nodes; ++i)
    if((x[i] >= x_min) & (x[i] <= x_max) & (y[i] >= y_min) & (y[i] <= y_max))
      ++(*n_keeps);
  log_info("Keeping %lu nodes.", *n_keeps);
  *keeps = malloc(*n_keeps * sizeof(int));
  if(*keeps == NULL)
    log_error("Malloc error allocating keeps");
  int idx = 0;
  for(size_t i = 0; i < n_nodes; ++i)
    if((x[i] >= x_min) & (x[i] <= x_max) & (y[i] >= y_min) & (y[i] <= y_max))
      (*keeps)[idx++] = (int) i;
}

void load_elems(const int ncid, size_t * n_elems, int ** elems) {
  nc_inq_dimlen_name(ncid, "nele", n_elems);
  *elems = malloc(*n_elems * 3 * sizeof(int));
  nc_get_var_int_name(ncid, "nv", *elems);
}

/** \todo Check if the nodes and nv are 0 or 1 indexed. */

void choose_elems_to_keep(const size_t n_elems,
			  const int * restrict elems,
			  const size_t n_nodes_keep,
			  const int * restrict nodes_keep,
			  size_t * n_elems_keep,
			  int ** elems_keep) {
  size_t n = 0;
  int * keeps = malloc(n_elems * sizeof(int));
  for(size_t i = 0; i < n_elems; ++i) {
    bool matches[3] = {false, false, false};
    for(size_t j = 0; j < n_nodes_keep; ++j) {
      if(elems[i] == nodes_keep[j])
	matches[0] = true;
      if(elems[n_elems + i] == nodes_keep[j])
	matches[1] = true;
      if(elems[2 * n_elems + i] == nodes_keep[j])
	matches[2] = true;
    }
    if(matches[0] & matches[1] & matches[2]) {
      keeps[i] = true;
      ++n;
    } else {
      keeps[i] = false;
    }
  }
  /* Pass 2: reduce the array. */
  *n_elems_keep = n;
  *elems_keep = malloc(n * sizeof(int));
  int idx = 0;
  for(size_t i = 0; i < n_elems; ++i)
    if(keeps[i]) 
      (*elems_keep)[idx++] = (int) i;
  log_info("Keeping %lu elements.", n);
  /* Cleanup. */
  free(keeps);
}

int create_output_file(const int ncid_in, const char * restrict filename, 
		       const size_t n_nodes_keep, const size_t n_elems_keep,
		       const int * restrict nodes_keep, 
		       const int * restrict elems_keep) {
  /* Create a mapping in the other direction. */
  size_t n_nodes = 0, n_elems = 0;
  nc_inq_dimlen_name(ncid_in, "node", &n_nodes);
  int * nodes = malloc(n_nodes * sizeof(int));
  for(size_t i = 0; i < n_nodes_keep; ++i)
    nodes[nodes_keep[i]] = i - 1;
  nc_inq_dimlen_name(ncid_in, "nele", &n_elems);
  int * elems = malloc(n_elems * sizeof(int));
  for(size_t i = 0; i < n_elems_keep; ++i)
    elems[elems_keep[i]] = i - 1;
  int ncid;
  /* Create the new file. */
  nc_check(nc_create(filename, NC_64BIT_OFFSET, &ncid));
  log_info("Created %s", filename);
  /* Read the number of dimensions in the input file. */
  int n_dims;
  nc_check(nc_inq_ndims(ncid_in, &n_dims));
  log_info("Transferring %d dimensions.", n_dims);
  int dimid_node, dimid_nele;
  /* Copy the dimensions. For the node and nele dimensions, resize them. */
  for(int i = 0; i < n_dims; ++i) {
    int dimid;
    size_t dimlen;
    char dimname[1024];
    nc_check(nc_inq_dimname(ncid_in, i, dimname));
    if(!strcasecmp(dimname, "node")) {
      dimlen = n_nodes_keep;
      dimid_node = i;
    } else if(!strcasecmp(dimname, "nele")) {
      dimlen = n_elems_keep;
      dimid_nele = i;
    } else {
      nc_check(nc_inq_dimlen(ncid_in, i, &dimlen));
    }
    nc_check(nc_def_dim(ncid, dimname, dimlen, &dimid));
  }
  /* Read the number of variables in the input file. */
  int n_vars;
  nc_check(nc_inq_nvars(ncid_in, &n_vars));
  log_info("Transferring %d variables.", n_vars);
  /* Copy the variables. */
  for(int i = 0; i < n_vars; ++i) {
    int varid;
    int vartype;
    nc_check(nc_inq_vartype(ncid_in, i, &vartype));
    char varname[1024];
    nc_check(nc_inq_varname(ncid_in, i, varname));
    int ndims;
    nc_check(nc_inq_varndims(ncid_in, i, &ndims));
    int dimids[ndims];
    nc_check(nc_inq_vardimid(ncid_in, i, dimids));
    log_info("Creating variable %s.", varname);
    nc_check(nc_def_var(ncid, varname, vartype, ndims, dimids, &varid));
  }
  /* Exit define mode. */
  log_info("All variables created, exiting define mode.");
  nc_check(nc_enddef(ncid));
  /* Transfer the data. */
  for(int i = 0; i < n_vars; ++i) {
    nc_check(nc_sync(ncid));
    log_info("Transferring variable %d", i);
    /* Read the variable info. */
    int vartype;
    nc_check(nc_inq_vartype(ncid, i, &vartype));
    int ndims;
    nc_check(nc_inq_varndims(ncid, i, &ndims));
    int dimids[ndims];
    nc_check(nc_inq_vardimid(ncid, i, dimids));
    int xtype;
    nc_check(nc_inq_vartype(ncid, i, &xtype));
    /* 
     * Transfer the variable. If the node or nele variable is included, then
     * iterate along it for the transfer. Otherwise just do a block transfer. 
     */
    int iter_dim = -1;
    for(int j = 0; j < ndims; ++j)
      if((dimids[j] == dimid_node) | (dimids[j] == dimid_nele))
	iter_dim = j;
    /* Compute the data size. */
    size_t data_size = 4;
    for(int j = 0; j < ndims; ++j) {
      size_t dimlen;
      nc_check(nc_inq_dimlen(ncid, dimids[j], &dimlen));
      data_size *= dimlen;
    }
    if(iter_dim == -1) { /* No node or nele, raw data copy. */
      log_info("Raw data transfer.");
      void * data = malloc(data_size);
      switch(xtype) {
      case NC_CHAR:
	nc_check(nc_get_var_text(ncid_in, i, data));
	nc_check(nc_put_var_text(ncid, i, data));
	break;
      case NC_FLOAT:
	nc_check(nc_get_var_float(ncid_in, i, data));
	nc_check(nc_put_var_float(ncid, i, data));	
	break;
      case NC_INT:
	nc_check(nc_get_var_int(ncid_in, i, data));
	nc_check(nc_put_var_int(ncid, i, data));
	break;
      default:
	log_error("Unsupported variable type.");
	break;
      }
      free(data);
    } else if(dimids[iter_dim] == dimid_node) {
      log_info("Nodal data transfer.");
      size_t start[ndims];
      size_t count[ndims];
      for(int j = 0; j < ndims; ++j) {
	start[j] = 0;
	nc_check(nc_inq_dimlen(ncid, dimids[j], &(count[j])));
      }
      data_size /= count[iter_dim];
      void * data = malloc(data_size);
      count[iter_dim] = 1;
      for(size_t j = 0; j < n_nodes_keep; ++j) {
	start[iter_dim] = nodes_keep[j];
	switch(xtype) {
	case NC_CHAR:
	  nc_check(nc_get_vara_text(ncid_in, i, start, count, data));
	  start[iter_dim] = j;
	  nc_check(nc_put_vara_text(ncid, i, start, count, data));
	  break;
	case NC_FLOAT:
	  nc_check(nc_get_vara_float(ncid_in, i, start, count, data));
	  start[iter_dim] = j;
	  nc_check(nc_put_vara_float(ncid, i, start, count, data));
	  break;
	case NC_INT:
	  nc_check(nc_get_vara_int(ncid_in, i, start, count, data));
	  char varname[1024];
	  nc_check(nc_inq_varname(ncid_in, i, varname));
	  if(!strcasecmp(varname, "nbsn") |
	     !strcasecmp(varname, "wet_nodes_time")) {
	    int n_data_points = data_size / 4;
	    for(int k = 0; k < n_data_points; ++k)
	      ((int*) data)[k] = nodes[((int*) data)[k]] + 1;
	  } else if(!strcasecmp(varname, "nbve") |
		    !strcasecmp(varname, "wet_nodes_prev_int")) {
	    int n_data_points = data_size / 4;
	    for(int k = 0; k < n_data_points; ++k)
	      ((int*) data)[k] = elems[((int*) data)[k]] + 1;
	  }
	  start[iter_dim] = j;
	  nc_check(nc_put_vara_int(ncid, i, start, count, data));
	  break;
	default:
	  log_error("Unsupported data type.");
	}
      }
      free(data);
    } else if(dimids[iter_dim] == dimid_nele) {
      log_info("Elemental data transfer.");
      size_t start[ndims];
      size_t count[ndims];
      for(int j = 0; j < ndims; ++j) {
	start[j] = 0;
	nc_check(nc_inq_dimlen(ncid, dimids[j], &(count[j])));
      }
      data_size /= count[iter_dim];
      void * data = malloc(data_size);
      count[iter_dim] = 1;
      for(size_t j = 0; j < n_elems_keep; ++j) {
	start[iter_dim] = elems_keep[j];
	switch(xtype) {
	case NC_CHAR:
	  nc_check(nc_get_vara_text(ncid_in, i, start, count, data));
	  start[iter_dim] = j;
	  nc_check(nc_put_vara_text(ncid, i, start, count, data));
	  break;
	case NC_FLOAT:
	  nc_check(nc_get_vara_float(ncid_in, i, start, count, data));
	  start[iter_dim] = j;
	  nc_check(nc_put_vara_float(ncid, i, start, count, data));
	  break;
	case NC_INT:
	  nc_check(nc_get_vara_int(ncid_in, i, start, count, data));
	  char varname[1024];
	  nc_check(nc_inq_varname(ncid_in, i, varname));
	  if(!strcasecmp(varname, "nv")) {
	    int n_data_points = data_size / 4;
	    for(int k = 0; k < n_data_points; ++k)
	      ((int*) data)[k] = nodes[((int*) data)[k]] + 1;
	  } else if(!strcasecmp(varname, "nbe") |
		    !strcasecmp(varname, "wet_cells") |
		    !strcasecmp(varname, "wet_cells_prev_int") |
		    !strcasecmp(varname, "wet_cells_prev_ext")) {
	    int n_data_points = data_size / 4;
	    for(int k = 0; k < n_data_points; ++k)
	      ((int*) data)[k] = elems[((int*) data)[k]] + 1;
	  }
	  start[iter_dim] = j;
	  nc_check(nc_put_vara_int(ncid, i, start, count, data));
	  break;
	default:
	  log_error("Unsupported data type.");
	}
      }
      free(data);
    } else {
      log_error("Invalid dimension, not transferring data.");
    }
  }
  /* Close the file. */
  free(nodes);
  free(elems);
  log_info("Transfer complete, closing the file.");
  nc_close(ncid);
  return ncid;
}

void transfer_data(const int ncid_in, const int ncid_out, 
		   const size_t n_nodes_keep, const int * restrict nodes_keep,
		   const size_t n_elems_keep, const int * restrict elems_keep) {
  /** \todo */
}

int main(int argc, char * argv[]) {
  /* Read the command line arguments. */
  const char * input_file = argv[1];
  const char * output_file = argv[2];
  const float x_min = atof(argv[3]);
  const float x_max = atof(argv[4]);
  const float y_min = atof(argv[5]);
  const float y_max = atof(argv[6]);
  log_info("Subsetting %s into %s using %f < x < %f and %f < y < %f.",
	   input_file, output_file, x_min, x_max, y_min, y_max);
  /* Open the input file. */
  int ncid_in, ncerr;
  if((ncerr = nc_open(input_file, NC_NOWRITE, &ncid_in)) != NC_NOERR)
    log_error("Could not open input file with error %s", nc_strerror(ncerr));
  /* Open the input file and load x and y. These pointers need to be freed. */
  float * x, * y; 
  size_t n_nodes;
  load_xy(input_file, &x, &y, &n_nodes);
  /* Determine which nodes to keep. */
  size_t n_keeps;
  int * keeps;
  choose_nodes_to_keep(n_nodes, x, y, x_min, x_max, y_min, y_max, 
		       &n_keeps, &keeps);
  /* Load the elements. */
  size_t n_elems;
  int * elems;
  load_elems(ncid_in, &n_elems, &elems);
  /* Determine which elements to keep. */
  size_t n_elems_keep;
  int * elems_keep;
  choose_elems_to_keep(n_elems, elems, n_keeps, keeps, 
		       &n_elems_keep, &elems_keep);
  free(x);
  free(y);
  /* Create the new output file and transfer the data. */
  create_output_file(ncid_in, output_file, n_keeps, n_elems_keep,
		     keeps, elems_keep);
  /* Cleanup. */
  nc_close(ncid_in);
  free(elems_keep);
  free(elems);
  free(keeps);
  return 0;
}
