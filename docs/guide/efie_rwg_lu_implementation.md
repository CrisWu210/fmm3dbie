# EFIE (Galerkin + RWG) implementation guide based on fmm3dbie

This guide shows a minimal-rewrite path to build a new EFIE solver while reusing the existing RWG geometry, quadrature-correction, and operator-evaluation infrastructure already in the repository.

## 1) Reuse existing building blocks

- Geometry/RWG conversion: `open_go3_4_rwg` in `src/maxwell/fix_tri.f90`.
- RWG <-> point-value projection helpers: `rwg2vr`, `vr2rwg_vector`, `vr2rwg_scalar` in `src/maxwell/fix_tri.f90`.
- Near quadrature generation + FMM-accelerated operator application: `getnearquad_em_cfie_rwg_pec` and `lpcomp_em_cfie_rwg_pec_addsub` in `src/maxwell/em_cfie_rwg_pec.f90`.
- Existing end-to-end RWG example and data flow: `examples/maxwell/example_em_efie_rwg_pec.f90` and `examples/maxwell/pec/efie_rwg_open.f90`.
- Dense LU wrappers (complex): `zgausselim` in `src/common/lapack_wrap.f90`.

## 2) New solver skeleton (dense matrix + LU)

Create a new routine, e.g. `em_efie_rwg_lu_solver`, following the input/output style of `em_cfie_rwg_solver` but replacing GMRES by:

1. Assemble full dense matrix `A(nrwg,nrwg)`.
2. Assemble RHS by reusing existing RHS generation routines.
3. Solve `A x = b` with `zgausselim` (or direct LAPACK call if preferred).

A practical matrix-assembly strategy with maximal reuse:

- For each basis index `j = 1..nrwg`, set test vector `e_j` (Kronecker basis).
- Call existing operator matvec (`lpcomp_em_cfie_rwg_pec_addsub`) to obtain `y = K e_j`.
- Write `A(:,j) = y` (plus identity/scaling term if your EFIE formulation requires it).

This keeps the difficult parts (singular/near quadrature, oversampling, RWG projections, FMM far interactions) unchanged.

## 3) Suggested file organization

- New example driver: `examples/maxwell/example_em_efie_rwg_lu_pec.f90`
- New build recipe: `examples/maxwell/example_em_efie_rwg_lu_pec.make`
- Optional new core source: `src/maxwell/em_efie_rwg_lu_pec.f90`

Prefer to mirror naming, argument lists, and initialization steps from `example_em_efie_rwg_pec.f90` to reduce integration risk.

## 4) Validation path

1. Start from the same sphere mesh/input as existing RWG example.
2. Compare LU solution versus existing iterative solution at low/moderate DoF.
3. Reuse existing accuracy test utilities to compare scattered fields.
4. Record matrix condition estimate (`dcond`) returned by `zgausselim`.

## 5) Complexity notes

- Dense LU is `O(nrwg^3)` time and `O(nrwg^2)` memory.
- Recommended only for small-to-medium problems or as a reference solver.
- For larger systems, keep current FMM + GMRES path.

## 6) Minimal migration checklist

- [ ] Copy geometry/RWG setup block from `example_em_efie_rwg_pec.f90`.
- [ ] Reuse RHS generation (`get_rhs_em_cfie_rwg_pec` or `_sphere`).
- [ ] Reuse near/far setup from `em_cfie_rwg_solver` (once per matrix build).
- [ ] Build matrix by repeated calls to existing matvec routine.
- [ ] Call `zgausselim` and post-process solution.
- [ ] Reuse existing field/accuracy test routines for regression checks.

## 7) What to do next (hands-on sequence)

If you are starting right now, execute the following in order.

### Step A. Clone an existing RWG example into a new LU example

1. Copy `examples/maxwell/example_em_efie_rwg_pec.f90` to
   `examples/maxwell/example_em_efie_rwg_lu_pec.f90`.
2. Keep the geometry setup and RHS generation parts unchanged.
3. Replace only the solver call section with your new `em_efie_rwg_lu_solver`.

### Step B. Implement `em_efie_rwg_lu_solver`

Suggested structure:

1. Initialize geometry-near/far data exactly like `em_cfie_rwg_solver`.
2. Allocate `A(nrwg,nrwg)`, `rhs(nrwg)`, `soln(nrwg)`.
3. For `j=1..nrwg`:
   - set `e_j`;
   - call `lpcomp_em_cfie_rwg_pec_addsub(..., e_j, ..., y, ...)`;
   - write `A(:,j)=y` (+ identity/scale for your EFIE convention).
4. Call `zgausselim(nrwg, A, rhs, info, soln, dcond)`.
5. Reuse existing field/accuracy routines.

### Step C. Add an example make target

1. Copy `examples/maxwell/example_em_efie_rwg_pec.make` to
   `examples/maxwell/example_em_efie_rwg_lu_pec.make`.
2. Rename executable/object names consistently.
3. Keep the same library links (`-lfmm3dbie`, `-lfmm3d`).

### Step D. Validate from small to larger systems

1. Start with a small sphere mesh and loose tolerance.
2. Compare LU current against iterative current from the existing example.
3. Tighten tolerance and increase DoF gradually.
4. Monitor memory because dense LU scales as `O(nrwg^2)` storage.
