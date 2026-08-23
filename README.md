# SMASH Afterburner Analysis — Modular

Modular C++17 rewrite of the HBT/femtoscopy analysis that consumes SMASH
Afterburner output, prepares validated particles and traverses physical pairs
for the required primitive HBT channels.

The current executable performs configuration, HBT startup resolution,
subevent-local particle preparation, deterministic accepted-particle shuffling,
numerical-rejection reporting, physical pair traversal, pair kT/mT calculation
and validation, origin/slice routing, Lab-to-LCMS and LCMS-to-PRF
transformations, pair-frame observables, eager raw histogram accumulation, and
complete-sample post-processing. Outer events may be processed in parallel with
worker-private integer histogram state and deterministic reduction. After all
workers have joined, it selects statistical regions, evaluates normalized
exact-bin models, runs Gaussian and Gaussian-plus-exponential fits with ROOT
Minuit2, calculates signed delta-t moments, and writes canonical production CSV
output.

## Build requirements

### Required software

Required software:

- CMake 3.16 or newer.
- A C++17 compiler.
- `yaml-cpp` with a supported exported CMake target.
- ROOT with the Minuit2 component and the `ROOT::Minuit2` target.
- Doxygen with Graphviz `dot` available.

### Offline dependency policy

The project has a strict offline dependency rule:

> The build never downloads third-party packages or source code.

CMake resolves all third-party dependencies locally:

```cmake
find_package(Doxygen REQUIRED COMPONENTS dot)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Minuit2)
```

Both exported target names are supported:

```text
yaml-cpp::yaml-cpp
yaml-cpp
```

If a required local installation is unavailable, CMake configuration fails
immediately. ROOT must export `ROOT::Minuit2`; no exact ROOT version is pinned.
There is no `FetchContent`, vendoring, network fallback or automatic package
download.

### Compiler diagnostics

All project targets are compiled with:

```text
-Wall -Wextra -Wpedantic -Werror
```

## Build and test

### Default toolchain

Default compiler:

```bash
rm -rf build
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Different compilers must use different build directories so that CMake does not
reuse a cache created for another toolchain.

### Linux with GCC

```bash
rm -rf build-gcc docs/doxygen docs/doxygen-warnings.log && \
CC=gcc CXX=g++ cmake -S . -B build-gcc && \
cmake --build build-gcc && \
ctest --test-dir build-gcc --output-on-failure && \
./build-gcc/event_preparation_real_fixture_test config/main.yaml && \
cmake --build build-gcc --target doxygen

echo $?
```

### Linux with Clang

```bash
rm -rf build-clang docs/doxygen docs/doxygen-warnings.log && \
CC=clang CXX=clang++ cmake -S . -B build-clang && \
cmake --build build-clang && \
ctest --test-dir build-clang --output-on-failure && \
./build-clang/event_preparation_real_fixture_test config/main.yaml && \
cmake --build build-clang --target doxygen

echo $?
```

### macOS with AppleClang

The default Apple toolchain can be used with:

```bash
rm -rf build docs/doxygen docs/doxygen-warnings.log && \
cmake -S . -B build && \
cmake --build build && \
ctest --test-dir build --output-on-failure && \
./build/event_preparation_real_fixture_test config/main.yaml && \
cmake --build build --target doxygen

echo $?
```

On macOS, the commands `clang++`, `g++`, and `c++` may all resolve to
AppleClang. Always check the reported compiler version before assuming that a
different compiler has been tested:

```bash
clang++ --version
g++ --version
c++ --version
```

### Standard CTest suite

The standard suite currently contains 51 CTest tests. It includes the
`doxygen_contract`, `naming_catalog_test`, and deterministic
`subevent_shuffle_test` gates. The direct real-fixture regression is executed
separately by the release gate.

## Run the executable

The executable requires exactly one argument: the global run-configuration
file.

```bash
./build/smash-afterburner-analysis-modular config/main.yaml
/usr/bin/time -p ./build/smash-afterburner-analysis-modular config/main.yaml
```

`main()` only validates the command-line contract, delegates the run to
`app::AnalysisRunner`, and delegates result serialization to the `output`
module.

### Global progress display

Enabled-HBT production runs report one global progress display to `stderr`.
The completed-work unit is one fully processed subevent, and the denominator is
`number_of_events * number_of_subevents`. The observer is notified only after a
subevent has completed and its results have been accumulated; no terminal
updates occur in the per-particle or per-pair hot paths.

On an interactive terminal, the display uses ANSI bright green (`92`) and
re-queries the current terminal width before every redraw. The bar therefore
adapts to terminal resizing while the run is active. On narrow terminals the
bar is shortened first; subevent and event detail lines are compacted or
omitted when they no longer fit, so the renderer does not deliberately wrap a
progress line.

A typical interactive display is:

```text
HBT [=======================>------]  76.4%
    1528 / 2000 subevents
    event 1 / 1
```

The percentage describes completed subevents only. After all subevents are
processed, post-sample statistical analysis and production-output writing are
shown as textual stages rather than assigning invented percentages to fitting
or serialization work.

When `stderr` is redirected, the renderer emits plain text with no ANSI control
sequences. Progress is throttled to coarse 10% milestones, followed by explicit
stage records for post-processing, analysis completion, production-output
writing and final completion.

The progress observer is application-level orchestration only. It does not
participate in particle selection, pair processing, histogramming or fitting
decisions, and HBT scientific modules do not depend on the console renderer.

## Outer-event parallelism

Phase 8 parallelizes only the outer-event loop. Work is claimed dynamically so
workers can continue with the next event when event costs differ. Every worker
owns its mutable scientific accumulation state, including a private full raw
histogram state; workers do not share histogram bins.

After all workers have joined successfully, the main thread reconstructs
observable summaries in canonical outer-event/subevent order and reduces the
worker histogram states with checked integer additions. Phase 7, ROOT Minuit2
fitting, and production output then run serially on the final reduced state.

Explicitly classified recoverable failures are handled at the smallest safe
scope. A subevent-local failure is reported as
`SkippedDueToSubeventFailure`, contributes zero, and processing continues with
the next subevent. An event-local input failure is reported as
`SkippedDueToEventFailure`; the complete event is transactionally discarded,
including any preparation, pair, and raw-histogram contribution produced before
the failure, and processing continues with other outer events.

A fatal worker failure requests cancellation, joins all workers, and prevents
the scientific merge, Phase 7, and successful production output. Unclassified
exceptions are never converted into skipped states.

The determinism contract is exact within the same build and environment:
`threads: 1` and `threads: N` must produce the same integer scientific state,
ordered diagnostics, and production files. The thread scheduler must not become
an observable part of the analysis.

## Generate Doxygen documentation

Doxygen and Graphviz `dot` must already be installed locally. They are required
at CMake configuration time and are never downloaded automatically by the
project. `doxygen_contract` is also part of the standard CTest suite.

From a clean repository state:

```bash
rm -rf build docs/doxygen && \
cmake -S . -B build && \
cmake --build build --target doxygen
```

The generated HTML documentation is written to:

```text
docs/doxygen/
```

Open `docs/doxygen/index.html` to browse it. The Doxygen target scans
`README.md`, `src/`, and `tests/`; documentation warnings make the target fail.
Graphviz is enabled for class/collaboration, include/dependency, directory,
call and caller graphs.

## Repository architecture

### Source layout

```text
src/
├── app/
│   ├── analysis_progress.h
│   ├── analysis_runner.*
│   ├── console_progress_bar.*
│   └── hbt_config_gate.*
├── common/
│   ├── four_vector.h
│   ├── kinematics.*
│   └── kinematics_validation.*
├── config/
│   ├── run_config.h
│   └── run_config_loader.*
├── hbt/
│   ├── boosts/
│   ├── channels/
│   ├── config/
│   ├── event/
│   ├── fits/
│   ├── histograms/
│   ├── pair/
│   ├── reporting/
│   ├── selection/
│   ├── species/
│   └── startup/
├── input/
│   ├── afterburner_reader.*
│   ├── emission_point_resolver.*
│   └── sampler_reader.*
├── output/
│   ├── analysis_output_writer.*
│   └── hbt_production_output.*
└── main.cpp
```

### Responsibility boundaries

`common`
: Reusable mathematical calculations and non-fatal numerical validation.
  It has no HBT policy, event context or output behavior.

`config`
: Global run configuration and path resolution.

`input`
: Parsing Afterburner/Sampler input and resolving emission positions.

`hbt`
: HBT domain configuration, species/channels, selection, event preparation,
  pair traversal/counting, pair kinematics and validation, HBT frame boosts,
  OSL pair-frame observables, kinetic slice routing/counting, eager raw
  histogram allocation/accumulation, post-sample statistical regions,
  normalization, exact-bin models, Minuit2 fitting, startup derivation,
  pair-processing accounting and in-memory rejection reporting.

`app`
: Orchestration. It joins configuration, readers, HBT preparation and pair
  processing without owning pair mathematics or result formatting.

`output`
: The only result-serialization boundary. Scientific modules provide data;
  they do not choose output formats or write analysis-result files.

This output boundary is a strict architectural rule.

## Global run configuration

### Example

Example:

```yaml
events_path: "../../"
output_path: "../output"
number_of_events: 1
number_of_subevents: 2000
threads: 1
hbt_enabled: true
hbt_config_path: "hbt.yaml"
```

### Paths and subevent cardinality

`events_path` and `output_path` are mandatory non-empty scalar paths. Relative
paths are resolved relative to the directory containing the global run YAML.
`output_path` is a run-level filesystem setting, not an HBT scientific setting.
The production writer requires its root to be absent or empty so stale output
cannot be mistaken for the current run.

`number_of_subevents` is the strict expected number of subevents for each outer
event, not a global total. A Sampler/Afterburner cardinality mismatch confined
to one outer event is reported as `SkippedDueToEventFailure`: the event is not
silently truncated or padded, its complete scientific contribution is rolled
back, and processing may continue with other outer events.

`threads` is an optional run-level execution control and defaults to `1`. A
value of `1` uses the serial outer-event path, `0` requests automatic hardware
concurrency with a fallback of one worker, and values greater than `1` request
at most that many workers. The effective worker count never exceeds
`number_of_events`. This setting changes execution only; it does not change HBT
scientific configuration or post-sample analysis.

### Conditional HBT configuration loading

When HBT is disabled, the HBT configuration file is not loaded.

When HBT is enabled, `hbt_config_path` is mandatory and is resolved relative to
the directory containing the global run configuration.

## HBT configuration

The current HBT configuration contains five required scientific blocks:

- enabled analysis products/channels;
- particle acceptance;
- pair-slicing configuration;
- raw-histogram binning;
- origin mode.

There are no implicit scientific defaults.

### Enabled channels and products

`hbt_enabled_channels` describes requested final analysis products assembled
from canonical primitive channels.

Product configuration uses canonical ASCII channel tokens. Unicode aliases are
not part of the public configuration grammar. Parsing preserves the
primitive-channel order written by the user and also stores the configured
ASCII product expression, after trimming only outer horizontal whitespace, for
run-level traceability metadata.

A primitive channel may not be repeated inside one final product, and final
products may not be repeated semantically: `A+B` and `B+A` are the same product
for duplicate detection. Startup independently validates the same invariants.
Derivation of startup requirements deduplicates primitive channels globally so
a physical primitive channel is prepared only once even when several products
share it.

Each primitive-channel definition also owns its canonical text name. Parsing
and output therefore share the same channel-domain naming source instead of
maintaining duplicate token maps. The configured spelling is metadata only; it
does not change canonical channel identity, routing, fan-out or physics.

Required species are similarly derived once from the unique required primitive
channels.

### Particle acceptance

The longitudinal variable is explicitly configured as either:

```text
rapidity
pseudorapidity
```

Each species belongs to one configured cut group:

- pions;
- kaons;
- nucleons;
- sigmas;
- lambdas.

Particle-level cuts use open boundaries:

```text
|longitudinal variable| < longitudinal_abs_max
pt_min_gev < pT < pt_max_gev
```

Before any particle-level observable is calculated, the final Afterburner
four-momentum must be finite and its energy component must be strictly positive.
A non-positive energy is a recoverable numerical rejection and never reaches pT,
rapidity/pseudorapidity, origin classification, invariant-mass reconstruction,
Sampler lookup or emission-position resolution.

The positive-energy requirement is common to both longitudinal-variable modes.
Rapidity-specific validation reuses that guarantee instead of checking energy
again. Pseudorapidity does not use energy in its formula, but particles with
`E <= 0` have already been removed from the HBT flow.

A valid particle outside the configured physical cuts increments the ordinary
particle acceptance rejection counter. It is not a numerical anomaly.

### Pair slicing

`hbt_pair_slicing` independently configures `kt` and `mt` routing axes.

For each axis:

```yaml
kt:
  enabled: false
```

or:

```yaml
mt:
  enabled: true
  bin_edges_gev:
    - 0.50
    - 0.70
    - 0.90
```

Rules:

- `enabled` is mandatory;
- an enabled axis requires at least two bin edges;
- every edge must be finite and non-negative;
- edges must be strictly increasing;
- a disabled axis may omit `bin_edges_gev`;
- a disabled axis may retain a fully validated edge sequence;
- disabled axes do not route or filter pairs;
- consecutive edges define half-open slices `[edge_i, edge_(i+1))`.

The configuration supports no kinetic slicing, kT only, mT only, or the
Cartesian combination of both. Physical pair construction, pair counting,
kT/mT calculation and numerical validation, and pair-origin routing are now
part of the executable.

`PairSliceRouting` is connected to production pair processing. It receives
already calculated and validated pair kinematics and resolves at most one
destination from the enabled slicing axes. A value outside an enabled axis
range produces no slice route; this is not a numerical rejection and does not
change formed, valid or origin-routed pair counts.

When both axes are enabled, one physical pair maps to at most one Cartesian
`kT x mT` cell. Cells use deterministic kT-major ordering. When both axes are
disabled, no fictitious inclusive slice is created because the existing
`routed_P`, `routed_PR`, and `routed_PRD` counts already provide the inclusive
origin accounting.

The slice lookup is performed once for a physical pair. The resulting
`PairSliceRoute` is then reusable for every compatible requested origin route.
Under `OriginMode::All`, for example, a `P-P` pair uses the same resolved slice
for P, PR and PRD rather than repeating the slice lookup. The route also stores
its deterministic `flat_slice_index`, calculated once during routing and reused
by downstream slice consumers.

`PairSliceCountAccumulator` consumes the already resolved slice route, the
already resolved pair-origin routes and the primitive-channel index. It uses
`flat_slice_index` directly and increments counts by slice, origin and channel
without recalculating kT, mT or the slice layout. The hot-path bin lookup uses
the configuration already validated at startup and `std::upper_bound`; it does
not revalidate all bin edges for each pair.

`pair_processor` performs the slice lookup after kT/mT validation and
pair-origin routing. It commits per-slice counts only after every numerical
gate required for that pair succeeds. The runner accumulates the
local slice summaries into run totals. The output boundary serializes the
run-total slice summaries together with the exact validated `PairSlicingConfig`
retained in `AnalysisRunSummary.startup.hbt_config`; it does not define, copy or
reconstruct kinetic bins independently. Raw histogram accumulation consumes the
already resolved `flat_slice_index` directly rather than reconstructing kinetic
bins or reading the slice-count summaries.

For every origin route and primitive channel, production enforces:

```text
sum(slice counts) <= routed_origin_count
```

Under `OriginMode::All`, every populated slice also satisfies:

```text
slice_P <= slice_PR <= slice_PRD
```

The three origin counts overlap and must not be summed as a partition.

### Raw histogram binning

`hbt_histograms` configures three independent uniform raw-histogram families:

```yaml
hbt_histograms:
  osl:
    nbins: 3000
    min_fm: 0.0
    max_fm: 300.0
  radial:
    nbins: 3000
    min_fm: 0.0
    max_fm: 300.0
  delta_t:
    nbins: 6000
    min_fm_c: -300.0
    max_fm_c: 300.0
```

All keys are explicit. `nbins >= 1`, all boundaries must be finite, and each
maximum must be strictly greater than its minimum. OSL and radial require an
exact zero minimum. Delta-t retains its sign and requires `min < 0 < max`.
Ranges are half-open `[min,max)`. Values below the minimum increment underflow;
values at or above the maximum increment overflow. No clamp, epsilon, bound
swap or other numerical repair is applied.

The reciprocal uniform bin width is resolved once during configuration loading
and reused by accumulation. Histogram bins are independent of kT/mT slice
boundaries.

### Origin mode

Allowed values are:

```text
primordial
primordial_rescattering
primordial_rescattering_decay
all
```

Origin flags are inclusive. A primordial particle therefore belongs to all
three nested origin selections.

`PairOriginRouting` implements pair-origin membership without forming or
recalculating a physical pair. The exact P/R/D membership contract is:

```text
origin A  origin B  primordial  prim+rescatt  prim+rescatt+decay
P         P              yes           yes                  yes
P         R               no           yes                  yes
R         P               no           yes                  yes
R         R               no           yes                  yes
P         D               no            no                  yes
D         P               no            no                  yes
R         D               no            no                  yes
D         R               no            no                  yes
D         D               no            no                  yes
```

Here `P`, `R`, and `D` denote primordial, rescattering, and decay contributions.
The pair memberships are obtained by intersecting the two particles' inclusive
origin memberships.

`all` requests the three nested selections in the same run; it is not a fourth
physical origin category. Under `all`, one valid physical pair is routed to all
and only the compatible selections in the table above. A `P-P` pair therefore
routes to all three selections, while any pair containing `D` routes only to the
widest selection.

The routing component is connected to `pair_processor` in production. Routing
is resolved only after pair kT/mT has passed numerical validation. Route counts
are committed only after every numerical gate required for that pair succeeds,
so a rejected pair increments no origin route. The same physical pair and its
already calculated kT/mT are reused; no second pair traversal or kinematic
recalculation is introduced.

For an individual origin mode, only the requested route is active. Under `all`,
the nested production counts satisfy, for every primitive channel:

```text
routed_P <= routed_PR <= routed_PRD = valid
```

The three routed counts overlap under `all` and must not be summed as a
partition of `valid`.

## Current event-preparation pipeline

### Processing order

For each outer event, Sampler and Afterburner input are read from the same
outer-event directory. Afterburner subevent IDs define the independent
subevent units.

For every raw Afterburner particle, the current order is:

```text
raw record
  -> identify species from charge + PDG
  -> required species?
  -> validate finite Afterburner four-momentum
  -> require strictly positive Afterburner energy
  -> calculate pT once and validate it
  -> validate selected y/eta input domain
  -> calculate selected y/eta once and validate it
  -> apply configured pT and y/eta acceptance
  -> classify origin
  -> requested origin eligible?
  -> calculate invariant mass squared from Afterburner four-momentum
  -> validate invariant mass squared
  -> calculate and validate invariant mass
  -> resolve emission position
  -> construct final hbt::Particle
  -> append to the current-subevent accepted-particle staging vector

after the complete subevent has been prepared
  -> deterministic shuffle seeded by outer event + subevent ID
  -> move shuffled particles into EventBuffers
```

### Deterministic subevent shuffle

SMASH particle print order can carry non-physical ordering correlations. Pair
construction assigns first/second roles through container order, so signed
anti-symmetric observables such as `delta_t = t1 - t2` must not inherit that
input ordering.

The runner therefore stages the complete accepted and required particle set for
each subevent, then performs exactly one `std::shuffle` before species grouping
and pair formation. The shuffle acts on complete `hbt::Particle` values, so all
stored scientific fields remain attached to the same particle. Its local
`std::mt19937_64` seed depends only on the one-based outer-event number and the
Afterburner subevent ID. It does not depend on worker identity, scheduling, or
configured thread count.

Selection cuts, requested-origin filtering, invariant-mass validation and
emission-position resolution all occur before the shuffle. Pair construction
and origin/product/slice fan-out occur after it.

### Subevent pair traversal

After the current subevent has been fully prepared, `EventBuffers` is traversed
once for each unique required primitive channel while the buffers are still
alive.

Pair enumeration follows these rules:

- identical species use `i < j`, producing `N(N-1)/2` physical pairs;
- different species use the full Cartesian product `A x B`;
- required primitive channels are traversed in startup order;
- duplicate or invalid channel lists are structural errors;
- no particles are copied into a global pair vector;
- no pair is formed across different subevents.

The iterator returns ordered per-channel counts for the current subevent.
Those counts retain their original meaning: every physical pair formed before
pair-level numerical rejection or kinetic slicing.

### Pair kinematics and recoverable numerical rejection

Each formed physical pair is processed once while the current `EventBuffers`
remain alive. `PairKinematics` is calculated exactly once from the final
Afterburner momenta and the validated invariant masses already stored in the
two accepted particles. The total Lab-frame pair four-momentum is formed once
and retained for later frame transformations:

```text
P = p1 + p2
Kx = Px / 2
Ky = Py / 2
kT = sqrt(Kx^2 + Ky^2)

mavg = (m1 + m2) / 2
mT = sqrt(mavg^2 + kT^2)
```

The mT definition is fixed to the average particle mass; there is no runtime
mT-definition mode and the raw OSCAR mass column is not used.

Calculated kT is validated first, followed by mT. A non-finite result rejects
only that pair from later pair-observable processing. The two particles remain
valid and may participate in other pairs. There is no clamp, repair or
fallback.

Every rejected pair remains part of the formed-pair count and is stored in a
complete in-memory `RejectedPairReport`. Current reasons are:

```text
non_finite_kt
non_finite_mt
non_finite_delta_t_lab
non_finite_delta_t_lcms
non_finite_delta_t_prf
non_finite_r_out_lcms
non_finite_r_out_prf
non_finite_r_side
non_finite_r_long
non_finite_r_radial_lcms
non_finite_r_radial_prf
```

Each record retains outer-event number, subevent ID, primitive channel,
one-based deterministic pair ordinal within that channel, both particle
kinematic snapshots, calculated kT/mT, and the exact rejection reason.

For every primitive channel, both locally and in run totals, pair processing
must satisfy exactly:

```text
formed = valid + numerical_rejected
```

The rejection records must also reproduce the numerical-rejection counts
exactly. Any mismatch is a structural internal failure rather than a tolerated
analysis result. Ordered local counts are accumulated into run totals with
exact channel identity/order checks and explicit integer-overflow detection.

After kT/mT validation, each surviving pair is passed once to origin routing
and performs at most one kinetic-slice lookup. Those resolved routes are reused
without recalculating kT, mT, origin membership, or slice identity. When the
pair's path requires frame observables, final valid/origin/slice counts are
committed only after frame finiteness validation succeeds.

A pair outside the configured slice window remains valid and origin routed but
increments no slice count and does not enter frame calculation. With both
slicing axes disabled, no dummy slice is created and the slice summaries remain
empty.

### Pair-frame transformations and observables

Frame observables are calculated only after pair kinematics, kT/mT validation,
origin routing and the single slice lookup have completed. With at least one
slicing axis enabled, a pair outside the configured slice domain does not enter
the frame transformations. With both axes disabled, every pair that passes
kT/mT validation enters frame calculation. Whenever frames are required, their
finiteness gate is part of the pair's final numerical-validity decision; the
slice-admission lookup itself never creates a numerical rejection.

For every admitted pair, the relative Lab separation is formed once:

```text
delta_x_lab = x1 - x2
```

The HBT-specific boosts live in `src/hbt/boosts/`. Lab to LCMS uses the total
pair four-momentum already stored in `PairKinematics`:

```text
beta_LCMS = Pz / P0
```

For exact `beta_LCMS == 0.0`, Lab already is the pair LCMS and the relative
separation is reused without calculating a Lorentz factor. Any representable
non-zero beta, however small, executes the full longitudinal boost. No legacy
epsilon such as `1e-30` is used. The same locally calculated LCMS boost factor
is reused to derive the later out-direction PRF boost speed; individual
particle momenta are not boosted or re-summed.

In LCMS, the signed OSL basis uses the beam direction as long. For `kT > 0`,
out is defined by `K_T / kT` using the already calculated Kx, Ky and kT. For
exact `kT == 0.0`, qx, qy and qT are calculated lazily from the two particle
momenta and qT defines out. Accepted particles already satisfy `pT > 0`, so
this exact branch has `qT > 0`; there is no second degeneracy fallback and no
forced zero for out or side. The side direction is `ez x out`.

LCMS to PRF is a one-dimensional boost in the OSL out direction. It transforms
only `delta_t_lcms` and `r_out_lcms`; side and long are reused unchanged. For
exact `beta_out == 0.0`, PRF already equals LCMS and no Lorentz factor is
calculated. Non-zero beta values use the full boost without an epsilon.

`PairFrameObservables` retains exactly the transient values needed by later HBT
analysis:

```text
delta_t_lab
delta_t_lcms
delta_t_prf
r_out_lcms
r_out_prf
r_side
r_long
r_radial_lcms
r_radial_prf
```

The three temporal separations remain distinct for later frame comparisons.
Side and long are stored once because the LCMS-to-PRF boost leaves them
unchanged. Both radial distances use only spatial OSL components; time is not
included. The complete object is calculated once per admitted physical pair and
validated for finiteness before the downstream consumer is called. The Phase-6
raw-histogram accumulator consumes the same object synchronously by const
reference and does not retain per-pair frame state. Phase 7 later consumes only
the completed raw histogram counts, not the per-pair object.

### Raw histogram accumulation

Phase 6 materializes all raw-histogram state before the first subevent. State
exists only for requested final products and requested origins. With active
slicing, every origin owns one global destination plus the configured slice
destinations; the global destination is the union of those configured slices.
With slicing disabled, no dummy slice state is allocated.

Each destination contains nine logical raw histograms:

```text
OSL absolute marginals: |r_out_lcms|, |r_out_prf|, |r_side|, |r_long|
radial:                 r_radial_lcms, r_radial_prf
signed relative time:   delta_t_lab, delta_t_lcms, delta_t_prf
```

The four OSL histograms share one binning, the two radial histograms share one
binning and the three delta-t histograms share one binning. Mutable bins use
`std::uint64_t`, with independent underflow and overflow counters for every
logical histogram. Counter overflow is detected before wraparound.

A primitive-channel pair reaches the histogram boundary once. Absolute OSL
values and the nine bin classifications are calculated once and then reused
for every compatible final product, origin, global destination and routed
slice. Product fan-out is resolved before event processing and indexed directly
with the existing `channel_index`; `PairOriginRoutes` and `flat_slice_index` are
reused without route reconstruction. No histogram normalization or fit state is
stored in Phase 6.

### Subevent isolation

`EventBuffers` exists only for the current subevent. Pair traversal finishes
before the buffers are destroyed and before the next subevent begins.
Particles from different subevents are never mixed.

## Post-sample statistical analysis (Phase 7)

`hbt::analyze_histograms()` is called exactly once after all configured outer
events have been processed. It consumes the completed Phase-6 raw histogram
state and the validated HBT configuration. It does not modify raw counts,
revisit pairs, recompute kinematics/frames, rebuild routing or serialize files.

### Statistical regions and normalized distributions

For OSL and radial shape histograms, the selected region starts at bin zero,
retains leading empty bins, includes the modal peak and stops immediately before
the first empty bin to the right of that modal plateau. Later disconnected
islands are excluded. An all-zero histogram has no selected region.

For signed `delta_t`, the modal plateau is located first and the selected region
extends to the first empty bin on each side. No minimum-entry threshold or
silent tail repair is applied.

Presentation distributions are normalized only over the selected region. For a
uniform bin width `dx`, each selected raw count `n_i` is written as
`pdf_i = n_i / (N_selected * dx)` with counting error
`d_pdf_i = sqrt(n_i) / (N_selected * dx)`.

### Exact-bin shape models and likelihood

Shape fits use exact analytic bin-edge integrals; model values are never sampled
at bin centers. The implemented component shapes are:

```text
OSL Gaussian:       exp(-x^2 / (4 R^2))
radial Gaussian:    r^2 exp(-r^2 / (4 R^2))
OSL exponential:    exp(-x / R_tail)
radial exponential: r^2 exp(-r / R_tail)
```

The full statistical region keeps the historical contiguous-region semantics and
is used for presentation normalization and for the mixed model. The pure
Gaussian has a separate compact core region. A non-increasing PAVA estimate is
used only to choose the upper Gaussian edge; raw counts are never smoothed in
the likelihood. Radial core selection starts at the right edge of the modal
plateau and excludes the first PAVA bin at or below 10% of the raw modal
maximum. OSL core selection starts at `x = 0` and excludes the first PAVA bin at
or below 10% of the PAVA level at the origin. In both geometries the existing
first-empty-bin boundary is the safety limit and the available full-region edge
is the fallback when the 10% level is not reached.

Each component is independently normalized to unit probability over the region
used by that fit. The mixed probability is
`f_core * Gaussian + (1 - f_core) * tail`. There is no free amplitude and all
three mixed estimators use the same expected counts
`mu_i = N_selected * p_i`. The default production estimator is the binned
Poisson deviance, which retains observed zero-count bins. Two additional
independent mixed fits are evaluated with Neyman chi-square
`sum_(n_i>0) (n_i-mu_i)^2/n_i` and Pearson chi-square
`sum_i (n_i-mu_i)^2/mu_i`. Neyman deliberately omits observed zero-count bins,
matching the historical weighting used for the R_core(mT) comparison; Pearson
retains them. The probability vector must sum to one within an allowance
derived only from `double` roundoff and the selected-bin count. There is no
scientific epsilon, post-hoc renormalization, or free fit amplitude.

### ROOT/Minuit2 fitting

The pure Gaussian fit has one physical parameter, `R > 0`, represented through
a log-radius parameter. Poisson, Neyman and Pearson are fitted independently on
the same compact 10%-core region. Each estimator runs two deterministic starts:
the moment-derived Gaussian radius and `R_HM`. `R_HM` is obtained from a
shape-constrained PAVA estimate before linearly interpolating the half-maximum
crossing(s): OSL uses a non-increasing envelope from the origin, while radial
histograms use a least-squares unimodal regression with a non-decreasing branch
followed by a non-increasing branch. This prevents isolated radial endpoint bins
from defining the core mode. The resulting FWHM is converted to the model
radius; raw counts, `N_selected`, and fit regions are unchanged. Among valid
MIGRAD minima the smallest objective value is selected for MINOS. The mixed
model remains fitted over the full statistical region.

The mixed model fits `R_core > 0`, `R_tail > 0` and `f_core` in `[0,1]`. No
ordering between `R_core` and `R_tail` is imposed or used as a validity
criterion. Each estimator runs 36 deterministic starts, the Cartesian product
of `R_core = {R_G, 0.5 R_HM, R_HM, 2 R_HM}`,
`R_tail = {0.5, 1, 2} R_tail,mom` and
`f_core = {0.25, 0.50, 0.75}`. The `R_G` seed always comes from the pure
Gaussian fit using the same estimator. Poisson, Neyman and Pearson therefore
remain completely independent.

Numerically valid mixed MIGRAD minima are first grouped into numerical basins
using the final `(log(R_core), log(R_tail), f_core)` coordinates. The physical
Gaussian-core basin is identified by the observed half-maximum scale: for each
basin, the arithmetic mean of `log(R_core)` is compared with `log(R_HM)`, and
the basin minimizing `|mean(log(R_core)) - log(R_HM)|` is selected. This is a
relative-scale criterion and imposes no ordering between `R_core` and `R_tail`.
Within that selected basin, the valid minimum with the smallest objective value
is used for MINOS. Basin multiplicity is retained only as a stability diagnostic
and is not an acceptance veto. For diagnostic studies, the terminal physical
`R_core`, `R_tail`, and `f_core` coordinates of every one of the 36 starts are
also serialized, independently of whether that start is ultimately accepted.
These endpoint columns do not participate in fit selection. Poisson remains the
default estimator used by backward-compatible plotting columns.

The one-dimensional radial-mT count-threshold machinery is retained, but its
current value is `N_selected >= 0`, so no non-empty slice is vetoed by this
criterion. `N_selected` itself and its histogram-selection definition are
unchanged. Invalid MIGRAD/MINOS or covariance states, degenerate core fractions,
missing Gaussian anchors, invalid half-maximum seeds, invalid objectives and
insufficient regions are preserved explicitly. An invalid fit does not fabricate
parameter estimates or fit curves.

### Signed delta-t statistics

Signed `delta_t` does not use Minuit2. The selected raw counts define a weighted
bin-center mean, population sigma and
`sigma_error = sigma / sqrt(2 * (N_selected - 1))` when `N_selected > 1`.
Negative or non-finite computed variance is reported explicitly; it is not
clamped or repaired.

## Single-particle kinematics

### Calculation and validation responsibilities

`common::kinematics` is calculation-only. It does not validate values, apply
cuts, reject particles, print diagnostics or repair data.

`common::kinematics_validation` is validation-only. Its operations return
validity and have no side effects. They do not abort, throw, print, clamp or
apply HBT cuts.

The orchestration follows the pattern:

```text
validate required input domain
  -> calculate once
  -> validate calculated result
  -> apply analysis policy
```

### Positive-energy domain

After finite momentum has been established, particle acceptance requires the
final Afterburner energy `E` to be strictly positive. `E <= 0` is rejected once
at this common particle boundary before pT or the selected longitudinal
variable is calculated. The value of `E` is retained as the diagnostic scalar
for that rejection.

This guarantee is reused downstream. Rapidity validation does not repeat the
positive-energy check, and `EmissionPointResolver` receives only particles that
already satisfy `E > 0`.

### Invariant mass

The invariant mass used by HBT is reconstructed from the final Afterburner
four-momentum:

```text
m^2 = E^2 - px^2 - py^2 - pz^2
m   = sqrt(m^2)
```

The squared value must be finite and strictly positive. No tolerance, clamp or
replacement is applied.

The invariant mass is calculated once during particle preparation, before any
Sampler lookup or emission-position calculation that may be required for that
particle. A valid mass is stored directly in `hbt::Particle` for reuse by later
pair processing.

The raw OSCAR mass column is retained only as diagnostic input data; it is not
the HBT invariant-mass calculation.

## Recoverable numerical particle rejection

### Rejection policy

A numerical failure local to one particle does not terminate the run. That
particle is removed from the HBT flow and recorded in
`hbt::RejectedParticleReport`.

### Rejection reasons

Current numerical rejection reasons are:

```text
non-finite momentum
non-positive energy
non-finite transverse momentum
invalid rapidity input
non-finite rapidity
invalid pseudorapidity input
non-finite pseudorapidity
non-finite invariant mass squared
non-positive invariant mass squared
non-finite invariant mass
non-finite Sampler emission position
non-finite propagated emission position
non-finite Afterburner emission position
```

The stable serialized token for the energy-domain rejection is
`non_positive_energy`. Its diagnostic scalar is the rejected Afterburner energy.

No numerical rejection triggers a silent fallback, clamp or repaired value.

### Complete rejection record

Every rejected particle record stores:

```text
outer event number
subevent ID
particle ID
raw PDG
raw electric charge
identified SpeciesId
raw Afterburner four-momentum
raw Afterburner position
raw OSCAR mass column
ncoll
time_last_coll
exact rejection reason
optional diagnostic scalar
optional diagnostic position
```

`outer event number + subevent ID + particle ID` uniquely identifies the
particle within the processed production.

`RejectedParticleReport` is only an in-memory store. It owns complete records
and exact counts by rejection reason. It performs no I/O and makes no scientific
decisions.

## Failure scope and recovery contract

Recovery is classified by the smallest scope that can be discarded safely. It
is separate from particle-local numerical rejection and from valid empty
events/subevents.

- A missing mandatory Sampler match is a recoverable subevent-local failure.
  The current subevent is reported as `SkippedDueToSubeventFailure`, its staged
  accepted particles are discarded, it forms no pairs and contributes no raw
  histogram entries. The remaining declared rows are consumed so the reader
  stays aligned, then processing continues with the next subevent.
- Runtime failures from Sampler/Afterburner event input operations, including
  event-local open/read/structure failures and strict configured subevent
  cardinality mismatch, are recoverable at outer-event scope. The event is
  reported as `SkippedDueToEventFailure`; its preparation and pair summaries
  are reset and its event-local raw histograms are not merged. Other outer
  events may continue.
- Fatal conditions remain fatal when safe recovery cannot be guaranteed. These
  include invalid or missing global configuration, overflow, broken internal
  invariants/logic errors, failures from scientific pair or histogram
  processing, progress-infrastructure failures, and unclassified exceptions.
  In the parallel path they request cancellation, all workers are joined, and
  scientific merge, Phase 7, and successful production output are skipped.

`SkippedDueToEventEmpty` and `SkippedDueToSubeventEmpty` remain valid
zero-contribution states rather than failures. The final analysis summary
serializes aggregate counts plus diagnostics: event failures include outer-event
number and reason; subevent failures include outer-event number, subevent ID and
reason.

## Emission-position contract

Emission-position selection keeps a fixed precedence.

`EmissionPointResolver` is reached only after particle acceptance has
established strictly positive Afterburner energy. `E <= 0` is therefore not an
emission-position fallback condition: such a particle has already been rejected
and never reaches Sampler lookup, propagation or raw Afterburner selection.

### Mandatory Sampler branch

For:

```text
primordial && ncoll == 0
```

Sampler lookup by `(subevent_id, ID, PDG)` is mandatory.

A missing mandatory lookup is `SkippedDueToSubeventFailure`, not a fallback to
raw Afterburner position. The complete current subevent contributes zero; after
its remaining declared rows are consumed, processing continues with the next
subevent. The diagnostic records the outer event, subevent, particle ID and PDG.

A selected non-finite Sampler position is instead a recoverable numerical
rejection of that particle. It does not switch to another source and does not
invalidate the complete subevent.

### Propagation branch

When the propagation inputs are valid, the emission point is propagated back
from the Afterburner state using the particle velocity.

A non-finite propagated result rejects that particle numerically. It is not
clamped and does not fall back to raw Afterburner position.

### Raw Afterburner branch

When the propagation branch is not applicable, the raw Afterburner position is
selected.

A non-finite selected raw position rejects that particle numerically.

The final momentum always comes from Afterburner. Sampler may determine the
emission position but never replaces the momentum.

## Final accepted particle

`hbt::Particle` currently stores:

```text
SpeciesId
resolved emission position
final Afterburner four-momentum
validated invariant mass
inclusive OriginFlags
raw PDG
raw charge
```

It contains no pair observables, histogram state, product IDs or output-format
state.

## Analysis summary and output

### Run summary

`app::AnalysisRunSummary` owns the resolved startup state and, when HBT is
enabled, completed event-preparation and pair-processing summaries, the
complete-sample raw histogram state and the post-sample histogram analysis.

The post-sample result is stored separately from the raw histogram state and is
aligned with the same product/origin/global-or-slice layout. Raw `std::uint64_t`
histogram counts remain unchanged after accumulation.

The HBT preparation summary contains:

- processed outer-event and subevent counts;
- raw particle count;
- unsupported-species count;
- supported but unrequired-species count;
- ordinary physical acceptance rejection count;
- complete numerical rejection report;
- origin rejection count;
- accepted-particle count;
- emission-source counts;
- accepted counts by required species;
- accepted counts by subevent.

The HBT pair-processing summary contains:

- ordered run-total formed-pair counts for every required primitive channel;
- ordered run-total numerically valid pair counts after applicable gates;
- ordered run-total numerical pair-rejection counts;
- complete numerical pair-rejection records and reason counts;
- the configured pair-origin routing mode;
- ordered run-total `routed_P`, `routed_PR`, and `routed_PRD` counts;
- run-total slice counts indexed by slice, origin route and primitive channel;
- the same formed/valid/rejected, origin-route and slice counts for every
  subevent;
- outer-event number and subevent ID for each local summary.

### Output boundary

`src/output/analysis_output_writer.*` and `src/output/hbt_production_output.*`
are serialization-only boundaries. Scientific modules return data and never
write result files. The executable writes the diagnostic run summary to stdout
and writes Phase-7 production files under the resolved mandatory `output_path`.

The stdout summary includes particle/pair numerical rejections, formed/valid
pair counts, P/PR/PRD routing, slice counts and aggregate raw-histogram range
warnings. Pair names come from the canonical channel catalog.

Production output writes one run-level traceability table at the output root:

```text
<output_path>/product_catalog.csv
```

`product_catalog.csv` records `product_index`, the preserved
`configured_expression`, the catalog-derived `canonical_expression`,
`channel_index`, the canonical primitive-channel token, and both canonical
species tokens. One row is written per primitive channel in each configured
product, preserving configured channel order.

Scientific result files use this hierarchy:

```text
<output_path>/product_<index>/global/<P|PR|PRD>/
  LAB/dt/delta_t/
  LCMS/osl/r_out/
  LCMS/osl/r_side/
  LCMS/osl/r_long/
  LCMS/radial/r_radial/
  LCMS/dt/delta_t/
  PRF/osl/r_out/
  PRF/radial/r_radial/
  PRF/dt/delta_t/

<output_path>/product_<index>/<configured-slice-token>/<P|PR|PRD>/
  LCMS/radial/r_radial/
  LCMS/dt/delta_t/
  PRF/radial/r_radial/
  PRF/dt/delta_t/
```

Slice directory names are derived from the actual validated user-configured
kinetic edges, never from a hard-coded width. With mT edges `0.50, 0.70, 0.90`,
for example, the directories are `mT_slice0_0.5-0.7` and
`mT_slice1_0.7-0.9`. kT uses the analogous `kT_slice...` token. If both axes are
enabled, their tokens are joined with `__` in the same kT-major flat ordering
used by pair routing. OSL result directories are global-only; kinetic slices are
used for the radial `R_core` dependence (and the existing signed delta-t slice
outputs).

Shape directories always contain `fit_parameters.csv`; `distribution.csv` is
written when a statistical region exists. Distribution columns contain the
normalized `pdf` and `d_pdf`, plus Gaussian and/or mixed fitted PDF columns only
when the corresponding fit is fully valid. Gaussian fitted values are written
only inside the compact Gaussian-core region. The backward-compatible
`gaussian_fit_pdf` and `mixed_fit_pdf` columns are the Poisson curves; valid
alternatives are written as `gaussian_fit_pdf_neyman`,
`gaussian_fit_pdf_pearson`, `mixed_fit_pdf_neyman` and
`mixed_fit_pdf_pearson`. All mixed curves span the full statistical region. The
parameter table contains an explicit `estimator` column. It records the exact
region used by each fit, independent Poisson/Neyman/Pearson Gaussian `R_G_core`,
and independent mixed `R_core`, `R_tail`, and `f_core` values, asymmetric MINOS
errors, objective minima, covariance state, all 36 mixed-start diagnostics,
the terminal `R_core`, `R_tail`, and `f_core` coordinates for every mixed start,
and the selected-minimum basin multiplicity diagnostic.
Delta-t directories contain
`statistics.csv` with status, `N_selected`, mean, sigma and sigma error and,
when a region exists, `distribution.csv`.

The production writer consumes the already-derived Phase-7 state and performs no
scientific recalculation. Its output root must be absent or empty.

The diagnostic stdout pair-slice serialization uses the exact validated
`PairSlicingConfig` preserved in `AnalysisRunSummary.startup.hbt_config`. The
output module has no independent bin definition and does not reconstruct slice
limits. Before serialization it requires the pair-summary slicing configuration
to match the startup configuration exactly and requires the slice-count origin
mode to match the run-total origin-route mode. A pair summary without the
corresponding startup HBT configuration is rejected. Raw-histogram diagnostics
likewise require the exact configured product/origin/slice/family layout, and
configured histogram range boundaries are serialized with 17-digit double
precision.

The serialized pair-origin and pair-slice fields are:

```text
pair_origin_route_mode
routed_P_pair_counts_by_primitive_channel
routed_PR_pair_counts_by_primitive_channel
routed_PRD_pair_counts_by_primitive_channel
pair_slice_kt_enabled
pair_slice_kt_bin_edges_gev
pair_slice_mt_enabled
pair_slice_mt_bin_edges_gev
pair_slice_counts
```

Each `pair_slice_counts` entry stores zero-based `kt_slice_index` and
`mt_slice_index` values together with its P/PR/PRD primitive-channel counts. A
disabled axis is serialized with the stable index token `none`; retained
validated bin edges are still serialized exactly from the run configuration.

The radial-mT `N_selected` threshold machinery remains in place but is currently
set to zero. This allows all non-empty radial mT slices to attempt the fits while
preserving the same code path for future threshold studies such as 1000, 5000 or
10000.

`histogram_range_warnings` contains one aggregate entry only for each logical
histogram with non-zero underflow or overflow. Entries identify the final
product, origin, global or `flat_slice_index` destination, observable,
configured `[min,max)` range, and both counters. Individual pair values and
percentages are not serialized by this diagnostic.

Scientific modules do not directly serialize or write analysis outputs.

## Plot exported production distributions

`scripts/plot_results.py` renders the production CSV tree with Python and
Matplotlib. It is a presentation-only tool: it reads exported `pdf`, `d_pdf`,
the Gaussian and mixed fitted-PDF columns for the selected estimator and does
not refit, renormalize, interpolate, or reconstruct scientific results. Use
`--estimator {poisson,neyman,pearson}` to select both independently fitted
Gaussian and mixed estimators; `poisson` remains the default for backward
compatibility.

```bash
python scripts/plot_results.py /path/to/results
python scripts/plot_results.py /path/to/results --estimator neyman
python scripts/plot_results.py /path/to/results --estimator pearson --plots-dir /path/to/figures
```

`scripts/plot_paper_figures.py` accepts the same `--estimator` option. It uses
the corresponding Gaussian and mixed curves for Figures 2--4 and 7--8 and the
matching mixed `R_core` rows in `fit_parameters.csv` for Figure 9. Figure styles
and layouts are estimator-independent.

If `--plots-dir` is omitted, the script writes to a sibling directory named
`<results_root>_plots`. Matplotlib must already be installed in the Python
environment.

## Tests

### Standard CTest suite

The standard CTest suite contains 51 tests covering:

- Doxygen documentation with Graphviz navigation as a CTest contract;
- canonical species/channel naming and legacy-token rejection;
- species metadata and identification;
- four-vector data handling;
- common kinematic calculations;
- non-fatal kinematic validation;
- Afterburner parsing;
- Sampler parsing and lookup;
- emission-position resolution;
- particle acceptance decisions;
- event buffers;
- deterministic per-subevent accepted-particle shuffling, including stable
  event/subevent seeding and preservation of complete particle records;
- primitive-channel pair iteration;
- pair-count construction and accumulation;
- pair kT/mT calculation, including reusable total pair four-momentum;
- Lab-to-LCMS relative-separation transformation and exact zero-beta behavior;
- LCMS-to-PRF out-direction transformation and exact zero-beta behavior;
- composed Lab/LCMS/PRF pair-frame observables, including OSL orientation and
  the exact kT-zero qT fallback;
- pair-frame-observable finiteness validation before downstream consumption;
- pair kT/mT numerical validation;
- rejected-pair storage and reason accounting;
- formed/valid/rejected pair-processing partitioning;
- pair-origin membership and requested-origin routing, including all nine
  P/R/D combinations and `OriginMode::All`;
- production pair-origin routing through the pair processor and runner;
- kT-only, mT-only and Cartesian pair-slice routing, including half-open
  boundaries, out-of-range behavior and reusable flat slice indices;
- per-slice pair counting by origin route and primitive channel;
- production slice integration through the pair processor and runner,
  including out-of-range pairs and disabled slicing;
- rejected-particle storage;
- origin selection;
- channel catalog and requirement derivation;
- HBT configuration parsing and validation, including preservation of the
  configured ASCII product expression as traceability metadata;
- run configuration loading;
- HBT startup gating and derivation;
- startup-resolved primitive-channel-to-final-product histogram fan-out;
- eager raw-histogram allocation, half-open binning, marginal underflow/overflow
  accounting, checked counter overflow, origin/product/slice fan-out,
  test-instrumented proof of exactly nine bin classifications per consumed pair
  independent of destination fan-out, and complete-sample persistence;
- pair-processing-to-raw-histogram integration across independent subevents;
- statistical-region selection for shape and signed delta-t histograms;
- presentation normalization and counting-error densities;
- exact-edge Gaussian/exponential OSL and radial model integrals;
- unit-normalized model probabilities plus binned Poisson deviance, Neyman
  chi-square and Pearson chi-square, including their documented zero-count
  semantics and explicit probability-normalization validation;
- compact-core independent Poisson/Neyman/Pearson Gaussian fitting with moment
  and half-maximum starts, plus independent 36-start mixed Minuit2 fits that
  identify the R_HM-anchored Gaussian-core basin, select the smallest objective
  inside that basin, retain basin multiplicity and every start endpoint
  diagnostically, and validate MIGRAD/covariance/MINOS states and asymmetric
  physical errors;
- post-sample F6-to-F7 analysis without pair re-traversal or raw-count mutation;
- canonical production-output hierarchy, run-level `product_catalog.csv`
  traceability metadata, fit/statistics CSVs and omission of invalid fit curves;
- event-preparation integration;
- analysis-runner behavior;
- Phase-8 outer-event scheduling, worker-private accumulation, ordered merge,
  cancellation/join behavior and thread-count-independent results;
- output serialization, including production pair-slicing bins, run-total
  slice counts, kT-only/mT-only disabled-axis `none` indices, aggregate
  histogram-range diagnostics and inconsistent startup/summary routing-state
  rejection;
- executable startup contract.

### Numerical-rejection integration coverage

The integration suite explicitly verifies that non-positive-energy rejection
occurs before mandatory Sampler access, numerical mass rejection occurs before
Sampler access, and non-finite propagation rejects only the affected particle
while the run continues. It also verifies that pair processing crosses the
runner boundary with exact formed, valid, numerical-rejection, origin-route and
slice counts. Focused pair tests additionally cover the two HBT boosts, the
composed OSL/frame-observable path, exact zero-beta branches, the exact kT-zero
qT fallback, all eleven pair-level numerical rejection reasons, complete
diagnostic records, production origin routing, production slicing, slice-count
invariants, all nine frame-observable rejection mappings before histogram
consumption, stable rejection-token serialization, and raw histogram
fan-out/persistence.

## Real-fixture regression

### Regression executable

The real-fixture executable is always built:

```text
event_preparation_real_fixture_test
```

It is intentionally not part of the default CTest suite unless an explicit
configuration path is supplied at CMake configuration time.

### Direct execution

Direct execution:

```bash
./build/event_preparation_real_fixture_test config/main.yaml
echo $?
```

A successful regression returns exit code `0`.

### Optional CTest registration

It can also be registered in CTest explicitly:

```bash
cmake -S . -B build-real \
  -DHBT_REAL_REGRESSION_MAIN_CONFIG=config/main.yaml
cmake --build build-real
ctest --test-dir build-real --output-on-failure
```

### Frozen aggregate counts

The current frozen event-preparation regression is:

```text
outer events processed              1
subevents processed              2000
raw particles                  791392
unsupported species             41459
unrequired species             292628
physical acceptance rejections 391285
numerical rejections                0
origin rejections                   0
accepted particles              66020
emission positions: Sampler     11759
emission positions: propagation 54261
emission positions: Afterburner     0
```

### Frozen primitive-channel pair counts

The same real fixture now also freezes the run-total number of physical pairs
formed for every required primitive channel:

```text
pi_plus_pi_plus      178440
pi_minus_pi_minus    221678
p_p                     686
p_bar_p_bar              335
pi_plus_p              22503
pi_minus_p_bar         16654
k_plus_p                2864
k_minus_p_bar           1852
k_plus_k_plus           3070
k_minus_k_minus         2759
```

These totals were independently reproduced with the legacy pair loops for the
inclusive `primordial_rescattering_decay` selection. All 10 primitive-channel
formed-pair counts match the modular implementation exactly. The comparison is
before pair-level kT/mT cuts or slicing, matching the modular formed-pair
contract.

For this real fixture, every formed pair also has finite kT and mT. The frozen
regression therefore requires the valid-pair count to equal the formed-pair
count for all 10 channels and requires zero numerical pair rejections.

### Frozen pair-origin route counts

The same regression runs with `OriginMode::All` and freezes the three nested
legacy origin selections for every required primitive channel:

```text
channel                 routed_P  routed_PR  routed_PRD
pi_plus_pi_plus             7713       8754      178440
pi_minus_pi_minus          10241      11774      221678
p_p                           17         20         686
p_bar_p_bar                    6          8         335
pi_plus_p                    693        828       22503
pi_minus_p_bar               559        682       16654
k_plus_p                     164        188        2864
k_minus_p_bar                125        143        1852
k_plus_k_plus                495        569        3070
k_minus_k_minus              459        494        2759
```

These 30 routed counts reproduce the legacy nested primordial,
primordial+rescattering, and primordial+rescattering+decay diagnostics exactly.
Because this fixture has zero numerical pair rejections, `routed_PRD` also
equals both `valid` and `formed` for every channel. This equality is
fixture-specific; in a run with a non-finite pair, the general contract remains
`formed = valid + numerical_rejected` and, under `all`,
`routed_PRD = valid`.

## Current validation status

### Standard suite

The Phase-8 implementation plus deterministic subevent shuffle, before the
later recoverable event/subevent failure classification described above, passed
the full official release gate on all three required toolchains:

```text
macOS / AppleClang 21.0.0.21000101   51/51 PASS
Linux / GCC 11.5.0                    51/51 PASS
Linux / Clang 21.1.8                  51/51 PASS
```

Each gate completed the full build, all 51 CTest tests, the direct real-fixture
regression, the Doxygen target, and exited with status zero. This includes the
deterministic `subevent_shuffle_test` together with the existing architecture,
numerical, integration, ROOT/Minuit2, production-output and Doxygen coverage.

The recoverable event/subevent failure implementation was added after those gate
logs. The current tree therefore requires a fresh AppleClang/GCC/Clang gate
before those PASS results can be attributed to this later recovery patch.
Focused integration coverage now exercises a missing mandatory Sampler match as
`SkippedDueToSubeventFailure` and a late event-local cardinality failure as
`SkippedDueToEventFailure`, including rollback of already accumulated event
science.

### Thread-count determinism and timing

A production determinism check using 10 outer events and 2000 subevents per
event compared `threads: 1` with `threads: 10`. After normalizing only the
different output-root directory names, all 271 production files were
byte-identical. The common normalized tree fingerprint was:

```text
1cecfaafb10e9dc062ada6cda353b488022b6419ece04ab809a7aa10b828def8
```

The same workload completed in 87.47 s wall time with one thread and 29.09 s
with ten threads, corresponding to about 3.01x speedup. This is an illustrative
timing check, not the final Phase-8 scaling benchmark. High-statistics
scientific validation of the corrected signed delta-t distributions is still in
progress and is not claimed by these build/determinism gates.

### Real-fixture regression

The direct real-fixture regression on the current tree with traceability and
progress reporting returns exit code `0` on AppleClang, Linux GCC and Linux
Clang. The same frozen formed-pair counts, valid-pair equality and zero
numerical pair rejections pass on all three toolchains. All 30 frozen P/PR/PRD
route counts pass as well.

This validates production `PairOriginRouting` against the legacy nested origin
diagnostics while preserving the frozen pair totals. Production slice routing
and slice counting are also traversed with both slicing axes disabled by the
current fixture configuration. The regression requires empty slice summaries,
confirming that disabled slicing creates no dummy slice and preserves all frozen
production results. With slicing disabled, the same production pair-processing
path reaches frame-observable calculation for every pair passing kT/mT
validation; the fixture's finite frames then reach the Phase-6 raw-histogram
accumulator without changing the frozen pair accounting. The runner also
completes Phase-7 post-processing; focused unit/integration tests validate the
fit and production-output contents in detail.

## Development rules

Current implementation rules:

- C++17.
- `-Wall -Wextra -Wpedantic -Werror`.
- No compiler warnings accepted.
- No automatic dependency downloads.
- No implicit scientific defaults.
- No silent numerical repairs, clamps or source fallbacks.
- No duplicate expensive physical calculations when a validated value can be
  reused.
- No cross-subevent HBT pair formation.
- Scientific modules do not write result outputs.
- `main()` only validates its CLI and delegates work.
- Every function, class, struct, method and interface must be documented in
  Doxygen format, including parameters, return values, exceptions and
  preconditions when applicable.
- New behavior requires focused tests.
- Module crossings require integration coverage.
- Final validation always runs the complete standard CTest suite; focused tests
  do not replace the full-suite gate.
- Every implementation requires clean AppleClang, Linux GCC and Linux Clang
  validation plus the direct real-fixture regression.
- The README/roadmap is updated only after all required toolchains and the
  real-fixture regression pass with zero warnings and zero errors.
