#pragma once

// Optional: MATLAB Engine wrapper for the windowed GloptiPoly SOS solver.
//
// This is provided as a *bridge* so that a C++ pipeline can still invoke the
// exact same moment–SOS formulation you used in MATLAB (GloptiPoly 3).
//
// It is disabled by default. Enable by configuring CMake with:
//   -DMAGCAL_ENABLE_MATLAB_ENGINE=ON
// and then adding the MATLAB include/lib paths to your build.
//
// Notes:
// - MATLAB Engine C++ API differs between MATLAB versions.
// - This header intentionally contains only scaffolding and a clear call point.
// - A full, dependency-free SOS SDP implementation is nontrivial and is not
//   included in this lightweight C++ package.

#include "magcal/types.hpp"

namespace magcal {

#ifdef MAGCAL_ENABLE_MATLAB_ENGINE

// Pseudocode signature:
//
// Solve a small window of samples using GloptiPoly SOS in MATLAB and return (T,h).
//
// Mat3 solve_window_sos_matlab(const Mat3X& Yw, Mat3& T_out, Vec3& h_out);
//
// Implement this in your project by:
// - Starting MATLAB Engine
// - Transferring Yw into MATLAB
// - Calling an m-file that wraps your GloptiPoly formulation
// - Fetching T,h

#endif

} // namespace magcal
