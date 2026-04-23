/*
 * Lagrange Physics Library
 * Header-only C library for physics simulation
 * 
 * Version: 1.0.0
 * License: MIT
 * 
 * USAGE:
 *   Include this header in your source files:
 *     #include <lagrange.h>
 *   
 *   In exactly ONE source file, define the implementation:
 *     #define LAGRANGE_IMPLEMENTATION
 *     #include <lagrange.h>
 * 
 * DEPENDENCIES:
 *   - C11 or later
 *   - Standard C library (stdlib, math, string)
 *   
 * COMPILER FLAGS:
 *   -lm (link math library on Unix)
 */

#ifndef LAGRANGE_H
#define LAGRANGE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Version info */
#define LAGRANGE_VERSION_MAJOR 1
#define LAGRANGE_VERSION_MINOR 0
#define LAGRANGE_VERSION_PATCH 0

/*============================================================================
 * Core Modules
 *===========================================================================*/

#include "lagrange/types.h"
#include "lagrange/math.h"
#include "lagrange/transform.h"
#include "lagrange/body.h"
#include "lagrange/collider.h"
#include "lagrange/integrator.h"
#include "lagrange/gravity.h"
#include "lagrange/world.h"
#include "lagrange/sim.h"
#include "lagrange/transfer.h"
#include "lagrange/attitude_control.h"
#include "lagrange/koopman.h"
#include "lagrange/wavelet.h"
#include "lagrange/adaptive_integrator.h"
#include "lagrange/regularization.h"
#include "lagrange/autodiff.h"

/*============================================================================
 * Implementation
 *===========================================================================*/

#ifdef LAGRANGE_IMPLEMENTATION

/* All implementations are in the header files above as static inline functions.
 * No separate implementation file is needed for the core library.
 * This section is kept for future expansion if needed.
 */

#endif /* LAGRANGE_IMPLEMENTATION */

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_H */
