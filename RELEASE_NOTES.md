## Release d3.0 - Architecture Modernization & Bounding Rectangles Algorithm

### Overview

This release continues development on the **degenerate_cases** branch, building on the full
degenerate-case support shipped in d2.0. The focus of this release is a major architectural
modernization of the codebase together with the addition of a new bounding-rectangles-based
intersection algorithm.

---

### New Features

#### Bounding Rectangles Intersection Algorithm (`-a32`)
A new algorithm that finds all pairs of segments whose **axis-aligned bounding rectangles**
overlap. It uses a sweep-line approach over X and a segment tree / list structure for Y, giving
an efficient filter that can serve as a pre-pass or a standalone intersection detector. Select
it with `-a32` (or include `32` in the bitmask for `-a`). The algorithm is fully supported for
integer-coordinate segments (`-si`).

#### Thread Count Option (`-Tn`)
The number of threads used by the parallel Balaban algorithm (`-a16`) can now be set explicitly
with `-Tn` (default 8). Previously the thread count was hard-coded.

#### Benchmark Repeat Option (`-Rn`)
A new `-Rn` flag controls the minimum number of timing repetitions performed in release builds
(default 3). Useful for reducing measurement noise on fast runs.

#### Windows Thread Affinity and Priority (`-An`, `-Pp`)
On Windows builds it is now possible to:
- Pin the main thread to a specific logical CPU with `-An` (0-based index; `-A-1` disables
  affinity, which is the default).
- Set process/thread priority with `-Pp`:
  - `p=0` — no change (default)
  - `p=1` — ABOVE_NORMAL process + THREAD_PRIORITY_ABOVE_NORMAL
  - `p=2` — HIGH process + THREAD_PRIORITY_HIGHEST
  - `p=3` — REALTIME process + THREAD_PRIORITY_TIME_CRITICAL

---

### Architecture Improvements

#### C++20 Concepts for Segment Collections
Formal C++20 `concept` constraints (`SegmentCollection`, `SegmentCollectionWithSentinels`, …)
are now used throughout `new_int_finder.h`, `fast_finder.h`, and `optimal_finder.h`. This
catches interface mismatches at compile time and serves as machine-checked documentation for
implementors of custom segment collections.

#### CRTP `CollectionBase` - Unified Segment Collections
All four concrete segment-collection classes (`CLine1SegmentCollection`,
`CLine2SegmentCollection`, `CArcSegmentCollection`, `CGraphSegmentCollection`) now inherit
from the new `CollectionBase<Derived, Registrator>` CRTP template. Shared bookkeeping
(endpoint arrays, rank mappings, clone/unclone logic, InsDel dispatch, bubble-sort helpers)
lives in the base, eliminating duplicated code.

#### Generic `CRemaper`
The remapper that splits collinear overlapping segments is now templated on the underlying
collection type, so it can pre-process any collection that satisfies the collection concept —
not just `CIntegerSegmentCollection`.

#### `uint4` Indices Throughout
All segment indices and counts have been changed from signed `int4` to unsigned `uint4`. This
removes a class of signed/unsigned comparison warnings and makes the API contracts clearer.

#### Refactored `InsDel` and Endpoint Bubbling Interface
The `InsDel` logic and the endpoint-bubbling interface have been cleaned up to provide a
single, well-defined entry point for inserting and deleting segments at stripe boundaries,
making it easier to follow the degenerate-case handling path.

#### Sentinel Segments for `CIntegerSegmentCollection`
Sentinel (guard) segments are now added automatically by `CIntegerSegmentCollection`,
consistent with the sentinel contract described in the wiki. This removes the need for callers
to inject sentinels manually.

#### Optional Parallel Sorting of Endpoints
Endpoint sorting during collection preparation can now run in parallel (controlled by a
compile-time switch), providing a measurable speedup for very large inputs.

#### Improved In-Source Documentation
`fast_finder.h` and `optimal_finder.h` now carry detailed comments explaining the strip
subdivision strategy, the role of `FindR` / `FindRNoChecks` / `SISFindR`, and the
`MultipleCutting` dispatch.

---

### Notes

- This release requires a **C++20**-capable compiler (MSVC 2019 16.11+, GCC 10+, Clang 12+).
- The bounding-rectangles algorithm (`-a32`) replaces the old `bentley_ottmann` slot (which
  was removed in d1.0) in the algorithm bitmask.
- For standard (non-degenerate) use cases continue to use the `master` branch releases (`v3.x`).
