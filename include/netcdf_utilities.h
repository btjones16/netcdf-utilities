#ifndef __NEMO_SRC_NETCDF_UTILITIES_H__
#define __NEMO_SRC_NETCDF_UTILITIES_H__

#include <stdbool.h>

/** 
 * \defgroup NetCDF NetCDF Utilities 
 *  \{ 
 */

/** 
 *  Execute stmt and checks if the return value is NC_NOERR. If so, no further
 *  action is taken. If not, a string version of the error is logged and the
 *  code passed in as the ellipses is executed. This is declared as a macro, so
 *  the error correction code executes with access to the variables in the 
 *  calling scope. This macro uses a variable named _nc_error, so any portion of
 *  the error correcting code that accesses or modifies a variable with the
 *  same name is likely to have unexpected consequences. 
 */
#define nc_check(stmt, ...)                                             \
  {                                                                     \
    int _nc_error;                                                      \
    if((_nc_error = stmt) != NC_NOERR) {                                \
      log_warning("%s.", nc_strerror(_nc_error));                       \
      __VA_ARGS__;							\
    }                                                                   \
  }

/** \defgroup nc_get_var nc_get_var_type_name
 * The group of macros named nc_get_var_type_name pair calls to nc_inq_varid() 
 * and nc_get_var_type(). Each macro accepts variable length arguments, which 
 * may contain code that executes if either nc_inq_varid() or nc_get_var_type()
 * fails. The macros use variables named _varid and _error internally, so any
 * portion of the error code that accesses variables with the same names may
 * have unintended effects. The macros call nc_check() internally, so the 
 * caller does not need to. See the documentation for the NetCDF C library
 * for additional details on allowable parameter types.
 * \param ncid NetCDF ID, from a previous call to nc_open or nc_create.
 * \param var Variable name.
 * \param fp Pointer to the location into which the data value is read. See the
 *           NetCDF C library documentation for additional details.
 * \param ... Arbitrary code to execute if either nc_inq_varid() or 
 *            nc_get_var_type() returns an error.
 * \{
 */

#define nc_get_var_float_name(ncid, var, fp, ...)        		      \
  {                                                                           \
    int _varid;							  	      \
    int _error = false;                                                       \
    nc_check(nc_inq_varid(ncid, var, &_varid), _error = 1; __VA_ARGS__); \
    if(!_error) nc_check(nc_get_var_float(ncid, _varid, fp), ##__VA_ARGS__);  \
  }

#define nc_get_var_int_name(ncid, var, ip, ...)        		              \
  {                                                                           \
    int _varid;							  	      \
    int _error = false;                                                       \
    nc_check(nc_inq_varid(ncid, var, &_varid), _error = 1; __VA_ARGS__); \
    if(!_error) nc_check(nc_get_var_int(ncid, _varid, ip), ##__VA_ARGS__); \
  }

/** \} */

/** 
 * \defgroup nc_get_vara nc_get_vara_type_name
 * Equivalent to the nc_get_var_type_name group of macros but for 
 * nc_get_vara_type NetCDF calls. Only the additional parameters from the 
 * nc_get_var_type_name are documented here. For addtional details, see the 
 * documentation for that group or for the NetCDF C libraries.
 * \param start An array of start indices for each dimension of the variable.
 * \param count An array of counts of the number of indices to read along each
 *              dimension.
 * 
 * \{
 */
#define nc_get_vara_float_name(ncid, var, start, count, fp, ...)              \
  {                                                                           \
    int _varid;							  	      \
    int _error = false;                                                       \
    nc_check(nc_inq_varid(ncid, var, &_varid), _error = 1; __VA_ARGS__); \
    if(!_error) nc_check(nc_get_vara_float(ncid, _varid, start, count, fp),   \
                         ##__VA_ARGS__);				      \
  }

#define nc_get_vara_int_name(ncid, var, start, count, ip, ...)                \
  {                                                                           \
    int _varid;							  	      \
    int _error = false;                                                       \
    nc_check(nc_inq_varid(ncid, var, &_varid), _error = 1; __VA_ARGS__); \
    if(!_error) nc_check(nc_get_vara_int(ncid, _varid, start, count, fp),     \
                         ##__VA_ARGS__);				      \
  }



/** \} */

#define nc_inq_dimlen_name(ncid, name, lengthp, ...)                          \
  {                                                                           \
    int _dimid;                                                               \
    int _error = false;                                                       \
    nc_check(nc_inq_dimid(ncid, name, &_dimid), _error = 1; __VA_ARGS__); \
    if(!_error) nc_check(nc_inq_dimlen(ncid, _dimid, lengthp), ##__VA_ARGS__); \
  }


/** \} */

#endif
