# WGSL Q24.40 FMM BMM Framework

## Overview

This project provides a WebGPU Shading Language (WGSL) library for performing high-precision fixed-point arithmetic using the **Q24.40** format (24 integer bits, 40 fractional bits). This format is emulated using standard WGSL `i32` and `u32` types. The library includes core arithmetic operations, complex number support, basic 3D point operations, and helper functions.

Additionally, this repository contains a structural framework for implementing the Fast Multipole Method (FMM) accelerated Boundary Element Method (BEM) in WGSL, specifically targeting 3D acoustics problems potentially using the Burton-Miller formulation. It defines the necessary WGSL kernel entry points and data flow for the standard FMM stages (P2M, M2M, M2L, L2L, L2P, P2P).

**Note:** While the core Q24.40 arithmetic library is functional and tested, the FMM and BEM specific mathematical functions within the FMM kernels (e.g., spherical harmonics, Bessel/Hankel functions, translation operators, BEM integral evaluations) are currently **placeholders**. The focus of the FMM part is the **kernel structure, data layout, and integration with the Q24.40 library**.

## Features

### Core Q24.40 Library (`bemMath.wgsl`)

*   **Q24.40 Representation:** Emulates 64-bit fixed-point numbers (`Q24_40` struct) using `i32` (high part, signed) and `u32` (low part, unsigned).
*   **Basic Arithmetic:**
    *   Addition (`q24_40_add`)
    *   Subtraction (`q24_40_sub`)
    *   Multiplication (`q24_40_mul` - includes robust 128-bit intermediate calculation)
    *   Division (`q24_40_div` - Newton-Raphson based)
    *   Square Root (`q24_40_sqrt` - Babylonian method based)
    *   Negation (`q24_40_negate`)
    *   Absolute Value (`q24_40_abs`)
    *   Comparison (`q24_40_lt`, `gt`, `le`, `ge`, `eq`, `neq`)
    *   Min/Max (`q24_40_min`, `q24_40_max`)
*   **Complex Number Support:**
    *   `ComplexQ24_40` struct.
    *   Complex arithmetic (`add`, `sub`, `mul`, `div`, `conj`, `abs_sq`, `abs`, `scale`, `scale_inv`, `axpy`).
*   **3D Point Operations:**
    *   `Point3D` struct.
    *   Vector subtraction (`point3d_sub`), dot product (`point3d_dot`), squared norm (`point3d_norm_sq`), norm (`point3d_norm`).
*   **Coordinate Transformations:**
    *   Cartesian to Spherical (`cart2sph` - returns `SphericalCoords` struct).
    *   Spherical to Cartesian (`sph2cart` - accepts `SphericalCoords` struct).
    *   *Note: Underlying trig functions (`sin`, `cos`, `acos`, `atan2`) are currently placeholders.*
*   **Helpers:** Packing/unpacking functions, 64-bit integer helpers, bit counting.

### FMM Kernel Framework (`fmmKernelsSphereBEM.wgsl`)

*   **Structural Implementation:** Provides WGSL compute shader entry points (`_main` functions) and core calculation function signatures for standard FMM stages:
    *   `P2M` (Panel-to-Multipole)
    *   `M2M` (Multipole-to-Multipole)
    *   `M2L` (Multipole-to-Local)
    *   `L2L` (Local-to-Local)
    *   `L2P` (Local-to-Panel/Point)
    *   `P2P` (Panel-to-Panel / Direct Near-Field)
*   **Data Structures:** Defines basic WGSL structs (`PanelData`, `MultipoleExpansion`, `LocalExpansion`) compatible with FMM concepts.
*   **Q24.40 Integration:** Demonstrates how to unpack data, call core Q24.40 math functions, and pack results within the kernel structure.
*   **Accumulation Strategy:** Implements result accumulation for L2L, L2P, and P2P stages using separate partial result buffers and final summation kernels (`combine_l2l_main`, `final_sum_main`) to avoid reliance on complex atomic operations for custom structs.

## BEM Context

The library and framework are designed with 3D Boundary Element Methods, particularly for acoustics (Helmholtz equation), in mind. The FMM kernel structure includes placeholders anticipating the needs of the Burton-Miller formulation, which involves combining multiple integral kernel types (G, dG/dn(y), dG/dn(x), d²G/dn(x)dn(y)).

*   Placeholders for BEM kernel evaluation (`evaluate_inner_bem`, `evaluate_outer_bem`, `evaluate_direct_bem_integral`) are included but require specific implementation based on the chosen BEM formulation and numerical integration techniques.
*   Singularity handling for the self-term (`P2P` where source=target) is crucial for BEM and is marked as a placeholder needing careful implementation.

## WGSL Implementation Details

*   The library relies on WGSL structs to represent Q24.40, ComplexQ24_40, Point3D, etc.
*   Data is transferred to/from the GPU using packed formats (`vec2<u32>` for Q24_40, `vec4<u32>` for ComplexQ24_40) via storage buffers. `pack_` and `unpack_` helper functions are provided.
*   Panel geometry and FMM coefficients are assumed to be stored in tightly packed `array<u32>` buffers, requiring careful matching between JavaScript packing and WGSL unpacking logic (`unpack_panel_data`, `unpack_vec3_q24_40`).
*   FMM kernels often require parameters (box indices, element counts, constants like wavenumber `k`, Burton-Miller `alpha`) passed via uniform buffers.

## Usage / Integration

1.  **Combine WGSL:** Concatenate `bemMath.wgsl` and `fmmKernelsSphereBEM.wgsl` into a single string for shader module creation.
2.  **JavaScript Driver:** An external JavaScript application is required to:
    *   Initialize WebGPU (`setup.js`).
    *   Load and compile the WGSL code (`setup.js`).
    *   Manage the BEM mesh data (panels, vertices, normals) and octree structure.
    *   Precompute FMM interaction lists (neighbor lists, M2L lists).
    *   Create and populate GPU buffers for geometry, the current solution vector (`lambda`), FMM coefficients (`M`, `L`), neighbor lists, uniforms, and partial/final result vectors.
    *   Orchestrate the FMM process by dispatching the appropriate WGSL kernels (P2M, M2M, M2L, L2L, Combine L2L, L2P, P2P, Final Sum) in the correct order, managing dependencies between stages.
    *   Integrate the FMM matrix-vector product calculation within an iterative solver like GMRES or FGMRES.
    *   Read back final results.
3.  **Testing:** Use `tests.js` for testing the core Q24.40 library functions and `testFMM.js` for testing the structure and data flow of individual FMM kernels.

## Current Status & Limitations

*   **Core Q24.40 Library:** Functional and passes basic arithmetic and complex number tests. Performance is dependent on WGSL compiler/GPU.
*   **FMM Kernel Structure:** The WGSL entry points, bindings, and basic loop structures for all FMM stages are implemented. The accumulation strategy uses staging buffers and summation kernels.
*   **PLACEHOLDERS:** The following crucial components are **placeholders** and **do not contain correct mathematical implementations**:
    *   Trigonometric functions for Q24.40 (`sin`, `cos`, `acos`, `atan2`).
    *   Evaluation of spherical harmonics (`evaluate_Ylm`), associated Legendre polynomials.
    *   Evaluation of spherical Bessel and Hankel functions (`evaluate_inner_bem`, `evaluate_outer_bem`).
    *   FMM Translation Operators (`apply_m2m_operator`, `apply_m2l_operator`, `apply_l2l_operator`).
    *   BEM Kernel Evaluations, including derivatives and hypersingular terms (`evaluate_direct_bem_integral`).
    *   Numerical quadrature and singularity handling within `calculate_p2p_interaction`.
*   **Atomic Operations:** The current accumulation strategy avoids WGSL atomics due to their limitations with custom structs and the complexity of atomic carry propagation for Q24.40.
*   **Optimization:** The WGSL code is not yet optimized for performance (e.g., workgroup memory usage, vectorization).

## Future Work / Roadmap

1.  **Implement Q24.40 Math:** Implement accurate fixed-point versions of trigonometric functions, inverse trig functions, and potentially other special functions needed for FMM/BEM (e.g., using CORDIC or Taylor series).
2.  **Implement FMM/BEM Physics:** Fill in the placeholder sections in the `calculate_...` functions with the correct mathematical formulas for Helmholtz FMM translations and Burton-Miller BEM kernel evaluations.
3.  **Implement Singularity Handling:** Develop robust numerical integration and analytical/regularization techniques for the self-interaction (P2P) kernel.
4.  **Develop Full JS Driver:** Create the complete JavaScript FMM driver to manage the octree and orchestrate kernel launches.
5.  **Integrate with Solver:** Integrate the FMM MatVec into a WGSL/JS implementation of GMRES or FGMRES.
6.  **Testing & Validation:** Create comprehensive tests comparing results against analytical solutions or validated libraries. Implement the `cpuVerify...` functions in `testFMM.js`.
7.  **Optimization:** Profile and optimize WGSL kernels for specific GPU architectures. Explore workgroup memory usage.
8.  **Atomics (Optional):** Re-evaluate using atomics if the `atomics` WGSL extension becomes widely available *and* an efficient, correct method for atomic Q24.40 addition (handling carries) can be devised.

## Dependencies

*   WebGPU compatible browser and hardware.
*   (Potentially) `packed_4x8_integer_dot_product` WGSL extension (for the core math library, though not directly used by the FMM structure). Check `navigator.gpu.wgslLanguageFeatures`.
