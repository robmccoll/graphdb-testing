
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * STINGER CONFIGURATION
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * This file contains all of the parts of STINGER that cannot be configured at 
 * run time.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ * 
 * VERTEX Configuration
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* Enable a per vertex key-value store that stores arbitrary data keys 
 * and values.  Uses up maximum number of vertices * sizeof(struct *)
 * aditional data when enabled.
 */
// #define STINGER_VERTEX_KEY_VALUE_STORE


/* Define a custom vertex type type (normally int64_t is used) 
 */
// #define STINGER_VERTEX_CUSTOM_TYPE
// typedef int64_t vtype_t;

/* Define a custom vertex weight type
 * NOTE: you must provide the vweight_fetch_add and vweight_atomic_fetch_add
 * listed below.  The atomic function must guarantee an atomic add if running
 * with any form of parallelism. By default, this is simply an alias 
 * for an atomic_fetch_add intrinsic
 */
// #define STINGER_VERTEX_CUSTOM_WEIGHT
// typedef int64_t vweight_t;
// vweight_t stinger_vweight_atomic_fetch_add(vweight_t * p, vweight_t val);
// vweight_t stinger_vweight_fetch_add(vweight_t * p, vweight_t val);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * EDGE Configuration
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//TODO currently this is still in stinger-defs.h
//#define STINGER_EDGEBLOCK_SIZE 32

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * PHYSMAP Configuration
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* This will reenable the old completely independent tree-based physical 
 * mapper.  It is highly memory inefficient and recommended only if you have
 * legacy code that makes use of it
 */
//#define STINGER_FORCE_OLD_MAP 1
