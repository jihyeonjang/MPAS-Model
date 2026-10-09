# RRTMGP integration design (photochem-dev)

## Status (2026-10-09)
This branch is **preparatory**, not a functional RRTMGP integration. The RRTMG behavior is unchanged except for CAM interface-array bounds fixes. Do not select RRTMGP in the namelist until an adapter, solver calls, coefficient initialization and diagnostics are implemented, compiled and tested.

## Source comparison
- Official https://github.com/earth-system-radiation/rte-rrtmgp : numerical RTE+RRTMGP engine, independent of GFS/CCPP. Prefer this source as an immutable external checkout pinned to a tested SHA.
- https://github.com/NCAR/ccpp-physics/tree/main/physics/Radiation/RRTMGP : GFS host wrappers, clouds, aerosol optics, overlap and SW/LW drivers; the CCPP repository itself pins a separate NCAR/rte-rrtmgp submodule. Treat CCPP sources as an interface *reference*, not a drop-in MPAS implementation.
- https://github.com/NOAA-EMC/ufsatm : the `ccpp/physics` entry is a separate submodule, and `ccpp/suites/suite_FV3_GFS_v17_p8_rrtmgp.xml` orders setup/pre, surface, cloud microphysics, overlap, diagnostics, aerosol optics, SW, LW, post.

## Recommended implementation
1. Add official RTE+RRTMGP to `src/core_atmosphere/Externals.cfg` at a *verified stable version* (currently official CMake reports 1.9.2), with `local_path = ./physics_rrtmgp/rte-rrtmgp`. Never edit this checkout. Audit `manage_externals` support for hash-pinning before writing the external.
2. In MPAS, introduce a new `mpas_atmphys_rrtmgp_adapter.F90` that owns coefficient data loading, gas-concentration objects, initialization/finalization, and 2D (ncol,nlay) chunks. The official core library does **not** do netCDF coefficient IO; implement MPI rank-root reading and rank-local initialization in MPAS, or reuse its *example* IO helper externally without editing it.
3. Add SW/LW calls through `mpas_atmphys_driver_radiation_sw.F` and `mpas_atmphys_driver_radiation_lw.F`, including allocation, unpacking and output mapping, but keep current RRTMG options unchanged.
4. Only after fully implemented, add `rrtmgp_sw`/`rrtmgp_lw` choices in `mpas_atmphys_control.F` and `Registry.xml`; config must fail fast for missing coefficients, incompatible cloud/aerosol options, absent gas fields, or half-enabled schemes.
5. Add CMake integration via the official library's CMake targets (rte and rrtmgp); use an out-of-tree Fortran module directory to prevent `machine.mod` conflicts already present in MPAS/CCPP photochemistry. Legacy Make needs separate archives/module paths and explicit dependency order. Build must work with both MPAS precision modes or explicitly prohibit unsupported ones.
6. Map MPAS bottom-to-top arrays to upstream convention after checking `top_at_1`. Pass layer/interface pressure [Pa], temperatures [K], dry-air molar gas mixing ratios, gas coefficient files, SW solar zenith cosine/TOA incident spectrum and broadband/bandwise surface direct/diffuse albedo, LW emissivity and skin temperature.
7. Retain prognostic `qo3` from `config_o3_source='prognostic'`; ozone **units** (moist/dry mass basis versus molar VMR) must be confirmed at the actual MPAS state source, and SW+LW must use identical O3. With climatological ozone, interpolate once to layers. Gas `ty_gas_concs%set_vmr` expects dry-air VMR, not mass fraction.
8. Clouds: convert MPAS condensate mixing ratios and layer mass to LWP/IWP [g m-2], handle cloud fraction and McICA masks deterministically per column/time step, and preserve Thompson cloud effective radii [micron] after validation. Do not simply treat fractional clouds as grid-box mean optical paths. Include snow/precipitation optical effects only when supported and validated; RRTMGP's native cloud model may not match legacy RRTMG.
9. Aerosols: RRTMG's existing fixed SW bands/Thompson optics **cannot** be inserted unchanged into RRTMGP spectral bands/g-points. Implement explicit spectral mapping or scientifically consistent aerosol optics; document a clear-sky/no-aerosol limitation for the initial prototype.
10. Derive heating from **flux divergence** on MPAS layer pressure differences and check sign, Pa/hPa conversions and W kg-1 -> K s-1. Reuse the existing radiation tendency/diagnostic pathways (SW/LW heating, OLR, TOA and surface SW/LW, clear sky, net/total flux).
11. Initialize coefficients once per MPI process, cache them, and process columns in chunks, minimizing temporary copies and redundant SW/LW preprocessing. Ensure thread-safety and deterministic cloud sampling under OpenMP.
12. Verify with upstream RFMIP clear-sky reference tests, isolated single-column clear/cloudy cases, night SW=0, pressure/energy closure, conservation and tracer ozone sensitivity; then compile and run MPAS legacy Make and CMake with MPI/OMP and make RRTMG-vs-RRTMGP climate and performance comparisons. **Differences are expected** and not automatically errors.

## Bugs / warnings found
- MPAS: CAM LW/SW drivers allocate `emstot_p` and `abstot_p` up to `kme` but write through `kte+1`. Fixed on this branch to allocate to `kme+1`. Also removed redundant LW down-flux allocation checks. Run bounds-checked regression.
- MPAS: examine `vinterp_ozn(1,ncols,... p2d(its:ite,...),o3clim_p(1,...)` in LW driver if the code ever uses non-unit first column indexing; do not classify as proven error without runtime evidence.
- Upstream unresolved GPU correctness reports: https://github.com/earth-system-radiation/rte-rrtmgp/issues/405 and /issues/406 . Pin/test the intended release; no upstream edits here.

## Acceptance gates
(A) untouched RRTMG numerical regression and CAM bounds-checked regression;
(B) independent state-to-optics conversion unit tests;
(C) RRTMGP clear-sky flux validation;
(D) all-sky cloud/aerosol validation;
(E) MPAS atmospheric radiation + ozone coupling;
(F) MPI+OpenMP determinism, CMake/Make and measured CPU scaling.

A production RRTMGP selection must NOT be exposed until gates A through E are passed.
