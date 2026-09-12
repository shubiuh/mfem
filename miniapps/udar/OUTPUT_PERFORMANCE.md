# UDAR derived-field output measurements

Measured on 2026-09-12 with the existing Windows Release oneAPI/MPI/MUMPS
build, one MPI rank, CPU execution, binary ParaView output and the production
writer's default one-subdivision sampling. CUDA, HIP, OCCA and libCEED are
disabled in this build. These are output-only measurements, not complete
Maxwell survey timings or accelerator performance measurements.
The installed hardware reports an Intel Arc 140T GPU; neither a CUDA/HIP
toolchain nor an enabled MFEM backend for this GPU is available in the current
build environment. GPU absence is not assumed.

## Method

The `udar_output_benchmark` CMake target calls `WriteUDARVTKFields`, the same
function used by the forward solver. It projects manufactured real/imaginary
Nedelec fields for all three transmitter orientations and computes RT magnetic
fields with the discrete curl. Conductivity is spatially varying and
anisotropic. Each size has three alternating raw/full output repetitions plus
one separate instrumented full-output pass. Raw contains the twelve E/H vector
fields; full also contains six derived fields per transmitter orientation:
E/H magnitudes, electric/magnetic energy, ohmic loss and Poynting vector.

Projection and curl preparation are outside the timed region. Save includes
mesh and field evaluation, binary serialization and file writing. Per-call
instrumentation separately times derived-coefficient Eval, including FE
interpolation, conductivity evaluation and coefficient algebra. Instrumented
timings include clock overhead and must be interpreted with the unprofiled
repetitions. No instrumented measurements are used in the raw/full medians.
Filesystem cache and host load contribute to run-to-run variability.

## Results

| Hex elements | ND order | Complex unknowns | Raw Save median (s) | Full Save median (s) | Profiled Save (s) | Derived Eval (s) |
|---:|---:|---:|---:|---:|---:|---:|
| 4,096 | 2 | 104,544 | 0.5264 | 1.4049 | 1.6895 | 0.7865 |
| 32,768 | 2 | 811,200 | 4.7101 | 13.5845 | 13.4218 | 5.9226 |

The larger case performs 4,718,592 derived-coefficient evaluations per station.
Derived evaluation occupies about 44% of the instrumented Save time, and full
output takes about 8.87 seconds more than raw output at the respective medians.
That difference includes extra serialization/writing as well as evaluation;
it is not an estimate of GPU kernel speedup. The larger full-output unprofiled
runs ranged from 13.28 to 16.02 seconds. The measurements establish material
host evaluation cost within output at a substantial FE size. They do not
establish its share of a real production solve.

The profiled and unprofiled full VTU files are byte-identical at both sizes:
profiling preserves output values, array names, ordering and topology.

Raw CSVs from this run are `build/udar-output-n16/timing.csv` and
`build/udar-output-n32/timing.csv`. Reproduce after initializing oneAPI:

```text
cmake --build build/oneapi-mpi-mumps --target udar_output_benchmark
build/oneapi-mpi-mumps/miniapps/udar/udar_output_benchmark.exe -n 16 -o 2 -r 3 -out build/udar-output-n16
build/oneapi-mpi-mumps/miniapps/udar/udar_output_benchmark.exe -n 32 -o 2 -r 3 -out build/udar-output-n32
```

## Device implementation decision

The output-cost trigger is met: host derived evaluation is material in this
large-mesh output workload. Production output remains on the host pending a
device implementation and validation; no GPU speedup is claimed.

Two concrete integration constraints remain:

1. The available MFEM build has no accelerator backend. CPU measurements cannot
   establish a transfer-inclusive GPU improvement or device memory correctness.
2. In this checkout `QuadratureInterpolator::SupportsFESpace` accepts VALUE,
   INTEGRAL and H_DIV mappings, but not H_CURL. The existing ND electric fields
   cannot simply be passed to its device physical-value interpolation. A
   dedicated ND sampling path is needed, including edge orientations,
   nonconforming constraints, physical transformations and exact ParaView
   sampling-point ordering. Moving only scalar algebra would leave the repeated
   host FE interpolation and introduce transfers.

The next implementation should batch E/H interpolation at the existing output
points, compute all derived arrays from the shared samples, and transfer final
arrays once for the host serializer. Retain host fallback for unsupported
meshes and devices. Before selecting it automatically, compare all derived
arrays and measure interpolation, kernel work, transfers and total Save on an
accelerator-enabled build. Use `-vtk-profile` on a representative forward survey
to assess the end-to-end benefit alongside this isolated output benchmark.
