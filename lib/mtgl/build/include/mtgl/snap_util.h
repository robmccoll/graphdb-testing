/*  _________________________________________________________________________
 *
 *  MTGL: The MultiThreaded Graph Library
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top MTGL directory.
 *  _________________________________________________________________________
 */

/****************************************************************************/
/*! \file snap_util.h

    \brief Serial implementation of snapshot functions for use on Linux/Mac.

    \author Karen Devine (kddevin@sandia.gov)

    \date 5/6/2009
*/
/****************************************************************************/

#ifndef MTGL_SNAP_UTIL_H
#define MTGL_SNAP_UTIL_H

// Don't use these stubs if actually working on the MTA/XMT.
#ifdef __MTA__

#include <luc/luc_exported.h>
#include <snapshot/client.h>

#else

#include <cstdio>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <iostream>

// Return types and data types used by SNAP
#include <sys/types.h>
#include <sys/stat.h>

#include <mtgl/types.hpp>

typedef int luc_error_t;
typedef uint64_t luc_endpoint_id_t;
#ifdef _WIN32
typedef struct _stat snap_stat_buf;
#else
typedef struct stat snap_stat_buf;
#endif

// ERROR CODES RETURNED BY SNAP
#define LUC_ERR_OK 1
#define SNAP_ERR_OK 1
#define LUC_ERR_TIMEOUT -13
#define SNAP_ERR_BAD_PARAMETER 2
#define SNAP_ERR_OPEN_FAILED 3
#define SNAP_ERR_READ_FAILED 4
#define SNAP_ERR_RESOURCE_FAILURE 5
#define SNAP_ERR_RESTORE_FAILED 6
#define SNAP_ERR_SNAPSHOT_FAILED 7
#define SNAP_ERR_WRITE_FAILED 8

// Define this to some dummy value if not defined
#ifndef SNAP_ANY_SW
#define SNAP_ANY_SW 0
#endif

namespace mtgl {

/// Initializes the snapshot system; we don't have to do anything.
inline luc_error_t snap_init()
{
  return SNAP_ERR_OK;
}

/// Writes contents of buffer buf to file fname.
inline luc_error_t snap_snapshot(
  char fname[],
  void *buf,
  size_t count,
  int64_t *err)
{
  FILE *fp = fopen(fname, "wb");

  if (fp == NULL)
  {
    fprintf(stderr, "ERROR in snap_snapshot:  File %s not found\n", fname);

    if (err != NULL) *err = -1;

    return SNAP_ERR_SNAPSHOT_FAILED;
  }

  fwrite(buf, 1, count, fp);
  fclose(fp);

  if (err != NULL) *err = 0;

  return SNAP_ERR_OK;
}

/// Reads contents of file fname into buffer buf.
inline luc_error_t snap_restore(
  char fname[],
  void *buf,
  size_t count,
  int64_t *err)
{
  FILE *fp = fopen(fname, "rb");

  if (fp == NULL)
  {
    fprintf(stderr, "ERROR in snap_restore:  File %s not found\n", fname);

    if (err != NULL) *err = -1;

    return SNAP_ERR_RESTORE_FAILED;
  }

  if (fread(buf, 1, count, fp) != count)
  {
    fprintf(stderr, "ERROR in snap_restore:  fread() and stat() file sizes "
                    "don't match.\n");

    return SNAP_ERR_RESTORE_FAILED;
  }

  fclose(fp);

  if (err != NULL) *err = 0;

  return SNAP_ERR_OK;
}

/// Stub not yet implemented since it is not yet used in MTGL.
inline luc_error_t snap_pwrite(
  char fname[],
  luc_endpoint_id_t swEP,
  const void *buf,
  size_t count,
  off_t offset,
  int64_t *err)
{
  fprintf(stderr, "ERROR:  snap_pwrite not yet implemented.\n");

  if (err != NULL) *err = -1;

  return SNAP_ERR_WRITE_FAILED;
}

/// Stub not yet implemented since it is not yet used in MTGL.
inline luc_error_t snap_pread(
  char fname[],
  luc_endpoint_id_t swEP,
  void *buf,
  size_t count,
  off_t offset,
  int64_t *err)
{
  fprintf(stderr, "ERROR:  snap_pread not yet implemented.\n");

  if (err != NULL) *err = -1;

  return SNAP_ERR_READ_FAILED;
}

/// \brief Returns info about a file.
///
/// We typedef snap_stat_buf to be of type struct stat; the two structs have
/// similar fields and, in particular, share the st_size field used in MTGL.
inline luc_error_t snap_stat(
  char  fname[],
  luc_endpoint_id_t swEP,
  snap_stat_buf *statBuf,
  int64_t *err)
{
#ifdef _WIN32
  int tmp = _stat(fname, statBuf);
#else
  int tmp = stat(fname, statBuf);
#endif

  if (err != NULL) *err = tmp;

  if (tmp == 0) {
      return SNAP_ERR_OK;
  } else {
      return SNAP_ERR_OPEN_FAILED;
  }
}

}

#endif

#endif
