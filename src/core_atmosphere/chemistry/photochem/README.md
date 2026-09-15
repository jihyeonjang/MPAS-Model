# NRL photochemistry component

This MPAS-owned component calls the NRL ozone and stratospheric-water schemes
from the pristine `ccpp-physics` external. The external tracks
`production/GFS.v17`; its observed review-time commit was
`8007f2c3e32e20a4aea56576ccf50aaca7d18f44` (traceability only).

Direct MPAS coupling preserves surface-to-TOA order. `k=1` is the lowest model
layer and is passed as upstream `k=1`; no generic physics-state reversal is
applied. Tracers are mapped directly and component tendencies are additive.
