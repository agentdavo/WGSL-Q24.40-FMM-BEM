# WGSL Q24.40 FMM BEM Framework

## Overview

This project provides a WebGPU Shading Language (WGSL) library for performing high-precision fixed-point arithmetic using the **Q24.40** format (24 integer bits, 40 fractional bits). This format is emulated using standard WGSL `i32` and `u32` types. The library includes core arithmetic operations, complex number support, basic 3D point operations, and helper functions.

Additionally, this repository contains a structural framework for implementing the Fast Multipole Method (FMM) accelerated Boundary Element Method (BEM) in WGSL, specifically targeting 3D acoustics problems potentially using the Burton-Miller formulation. It defines the necessary WGSL kernel entry points and data flow for the standard FMM stages (P2M, M2M, M2L, L2L, L2P, P2P).

**Note:** While the core Q24.40 arithmetic library structure is present and basic operations compile, the FMM and BEM specific mathematical functions within the kernels (e.g., spherical harmonics evaluation, Bessel/Hankel function evaluation, FMM translation operators, BEM integral evaluations including regularization) are currently **placeholders or simplified first attempts**. The focus of the FMM part has been the **kernel structure, data layout, integration with the Q24.40 library, and basic data flow.** Significant work remains to implement the correct physics and robust numerical methods.

## Features

### Core Q24.40 Library (`bemMath.wgsl`)

* **Q24.40 Representation:** Emulates 64-bit fixed-point numbers (`Q24_40` struct) using `i32` (high part, signed) and `u32` (low part, unsigned).
* **Basic Arithmetic:**
    * Addition (`q24_40_add`)
    * Subtraction (`q24_40_sub`)
    * Multiplication (`q24_40_mul` - includes robust 128-bit intermediate calculation)
    * Division (`q24_40_div` - Newton-Raphson based, requires careful handling near zero)
    * Square Root (`q24_40_sqrt` - Babylonian method based)
    * Negation (`q24_40_negate`)
    * Absolute Value (`q24_40_abs` - Placeholder, currently uses float conversion)
    * Comparison (`q24_40_lt`, `gt`, `le`, `ge`, `eq`, `neq` - Placeholders assumed)
    * Min/Max (`q24_40_min`, `q24_40_max` - Placeholders assumed)
* **Complex Number Support:**
    * `ComplexQ24_40` struct.
    * Complex arithmetic (`add`, `sub`, `mul`, `div`, `conj`, `abs_sq`, `abs`, `scale`, `scale_inv`, `axpy`).
* **3D Point Operations:**
    * `Point3D` struct.
    * Vector subtraction (`point3d_sub`), dot product (`point3d_dot`), squared norm (`point3d_norm_sq`), norm (`point3d_norm`).
* **Coordinate Transformations:**
    * Cartesian to Spherical (`cart2sph` - returns `SphericalCoord` struct).
    * Spherical to Cartesian (`sph2cart` - accepts `SphericalCoord` struct).
    * *Note: Underlying trig functions (`sin`, `cos`, `acos`, `atan2`) and `exp` currently use `f32` conversion and are placeholders needing native Q24.40 implementation.*
* **Special Functions (Placeholders/Simplified):**
    * Associated Legendre Polynomial calculation (`get_Ynm`) - Needs verification.
    * Spherical Bessel function `j_n` (`get_jn`) - Simplified implementation, uses `f32` math internally.
    * Spherical Hankel function `h_n` (`get_hn`) - Simplified implementation, uses `f32` math internally.
    * Spherical Harmonic evaluation (`eval_spherical_harmonic`) - Relies on `get_Ynm` and placeholder trig functions.
* **Helpers:** Packing/unpacking functions, 64-bit integer helpers, bit counting.

### FMM Kernel Framework (`fmmKernelsSphereBEM.wgsl`)

* **Structural Implementation:** Provides WGSL compute shader entry points (`_main` functions) and basic calculation function signatures (`fmm_p2m`, `fmm_m2m`, etc.) for standard FMM stages.
* **Data Structures:** Defines basic WGSL structs (`Panel`, `MultipoleExpansion`, `LocalExpansion`) compatible with FMM concepts.
* **Q24.40 Integration:** Demonstrates how to unpack data, call core Q24.40 math functions (though often using simplified FMM physics), and pack results within the kernel structure. Calls to `get_jn`/`get_hn` now pass the scalar `k*r` argument correctly wrapped as a complex number.
* **Accumulation Strategy:** Implements result accumulation for L2L, L2P, and P2P stages using separate partial result buffers and final summation kernels (`combine_l2l_main`, `final_sum_main`) to avoid reliance on complex atomic operations for custom structs.

## BEM Context

The library and framework are designed with 3D Boundary Element Methods, particularly for acoustics (Helmholtz equation), in mind.

* **Kernel Evaluation:** Includes functions for evaluating potential (`kernel_helmholtz_potential`) and force (`kernel_helmholtz_force`) kernels at a point, and basic quadrature integration (`eval_gauss_quadrature`, `eval_dGdn`).
* **Singularity Handling:**
    * Basic Gaussian quadrature is implemented in `eval_gauss_quadrature`.
    * `eval_dGdn` includes an *attempted implementation* of Brandão's finite part approach for the hypersingular kernel integral. **This implementation requires careful verification against references and robust handling of potential divisions by zero.**
    * The self-term case (source point = target point) in P2P interactions requires a dedicated, robust singularity subtraction or analytical method, which is currently **not correctly implemented**. The simple distance check in `eval_dGdn` is likely incorrect for the true self-term limit.

## WGSL Implementation Details

* The library relies on WGSL structs to represent Q24.40, ComplexQ24_40, Point3D, etc.
* Data is transferred to/from the GPU using packed formats (`vec2<u32>` for Q24_40, `vec4<u32>` for ComplexQ24_40) via storage buffers. `pack_` and `unpack_` helper functions are provided.
* Panel geometry and FMM coefficients are assumed to be stored in tightly packed `array<u32>` buffers, requiring careful matching between JavaScript packing and WGSL unpacking logic.
* FMM kernels often require parameters (box indices, element counts, constants like wavenumber `k`, Burton-Miller `alpha`) passed via uniform buffers.

## Usage / Integration

1.  **Combine WGSL:** Concatenate `bemMath.wgsl` and `fmmKernelsSphereBEM.wgsl` (if separated) into a single string for shader module creation.
2.  **JavaScript Driver:** An external JavaScript application is required to:
    * Initialize WebGPU (`setup.js`).
    * Load and compile the WGSL code (`setup.js`).
    * Manage the BEM mesh data (panels, vertices, normals) and octree structure.
    * Precompute FMM interaction lists (neighbor lists, M2L lists).
    * Create and populate GPU buffers for geometry, the current solution vector (`lambda`), FMM coefficients (`M`, `L`), neighbor lists, uniforms, and partial/final result vectors.
    * Orchestrate the FMM process by dispatching the appropriate WGSL kernels (P2M, M2M, M2L, L2L, Combine L2L, L2P, P2P, Final Sum) in the correct order, managing dependencies between stages.
    * Integrate the FMM matrix-vector product calculation within an iterative solver like GMRES or FGMRES.
    * Read back final results.
3.  **Testing:** Use `tests.js` for testing the core Q24.40 library functions and potentially `testFMM.js` (if developed) for testing the structure and data flow of individual FMM kernels.

## Current Status & Limitations

* **Core Q24.40 Library:** Compiles successfully. Basic arithmetic operations (`add`, `sub`, `mul`, `div`, `sqrt`, `negate`) and complex operations are structurally present. Accuracy and edge-case handling (especially for `div`, `sqrt`) depend heavily on the underlying `f32` conversions in placeholder functions. **Needs native Q24.40 implementations for trig, exp, comparisons, etc.**
* **FMM Kernel Structure:** The WGSL entry points, bindings, and basic loop structures for all FMM stages are implemented. Calls to helper functions like `get_jn`/`get_hn` use correctly formatted arguments (scalar `kr` wrapped as complex). The accumulation strategy uses staging buffers and summation kernels.
* **PLACEHOLDERS / NEEDS VERIFICATION:** The following crucial components **do not contain correct/complete mathematical implementations** and require significant work:
    * **Q24.40 Math:** Trigonometric (`sin`, `cos`, `acos`, `atan2`), exponential (`exp`), comparison (`lt`, `gt`, etc.), min/max functions need native implementations.
    * **FMM Physics:**
        * Evaluation of spherical harmonics (`get_Ynm`, `eval_spherical_harmonic`) needs verification.
        * Evaluation of spherical Bessel/Hankel functions (`get_jn`, `get_hn`) uses simplified, `f32`-based placeholders.
        * FMM Translation Operators (M2M, M2L, L2L) within the `fmm_...` functions are placeholders and need correct implementation based on formulas (e.g., using Wigner D-matrices or other standard methods).
        * L2P/M2P evaluation logic needs verification (e.g., use of Ynm vs its conjugate).
    * **BEM Physics / Quadrature:**
        * The Brandão finite part implementation in `eval_dGdn` needs verification against literature and robust division-by-zero handling.
        * Singularity handling for the self-interaction (P2P, `kernel_helmholtz_bem` when source=target) is not correctly implemented.
* **Atomic Operations:** The current accumulation strategy avoids WGSL atomics.
* **Optimization:** The WGSL code is not yet optimized for performance.

## Future Work / Roadmap

1.  **Implement Q24.40 Math:** Implement accurate, native fixed-point versions of `sin`, `cos`, `acos`, `atan2`, `exp`, comparisons, min/max (e.g., using CORDIC or Taylor series).
2.  **Implement FMM Physics:**
    * Verify/implement correct Spherical Harmonic calculations.
    * Implement accurate Q24.40 versions of Spherical Bessel/Hankel functions.
    * Fill in the FMM translation operator formulas (M2M, M2L, L2L) using appropriate coefficients (e.g., Clebsch-Gordan, Wigner symbols).
    * Verify L2P/M2P summation formulas (e.g., conjugation of Ynm).
3.  **Implement BEM Physics & Singularity Handling:**
    * Verify/correct the Brandão finite part implementation in `eval_dGdn` or replace with a different verified method.
    * Develop robust numerical integration and analytical/regularization techniques for the self-interaction (P2P) kernel.
4.  **Develop Full JS Driver:** Create the complete JavaScript FMM driver to manage the octree and orchestrate kernel launches.
5.  **Integrate with Solver:** Integrate the FMM MatVec into a WGSL/JS implementation of GMRES or FGMRES.
6.  **Testing & Validation:** Create comprehensive tests comparing results against analytical solutions or validated libraries. Implement the `cpuVerify...` functions in `testFMM.js`.
7.  **Optimization:** Profile and optimize WGSL kernels. Explore workgroup memory usage.
8.  **Atomics (Optional):** Re-evaluate using atomics if the `atomics` WGSL extension becomes widely available *and* an efficient, correct method for atomic Q24.40 addition can be devised.

## Dependencies

* WebGPU compatible browser and hardware.
* (Potentially) `packed_4x8_integer_dot_product` WGSL extension.
