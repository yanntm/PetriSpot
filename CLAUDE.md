# PetriSpot — Working Notes for Claude

## Project layout

PetriSpot is a C++23 command-line tool for P/T Petri nets in PNML: invariant
computation (P/T flows and semi-flows, mature, state of the art) and, under
construction, an explicit heuristic walk engine for reachability. It is a
subprocess of ITS-Tools (Java, `~/git/ITStools`), not a replacement for it.

* `Petri/` : the CMake project (`Petri/CMakeLists.txt`). Build with
  `cmake -S Petri -B build && cmake --build build`; `buildPetriSpot.sh` at
  root does the full static build with libexpat as the CI does. Binaries:
  `build/petri32|64|128` (integer width via `-DVAL`) and `kersconv`.
* `Petri/src/` : header-only templates on the integer type `T`, one `Petri.cpp`
  driver. Being reorganised into folders (see `WALK_PLAN.md` section 3):
  `core/` sparse vectors, matrices, net; `parse/` PNML and property parsers;
  `expr/` property AST; `invariants/` the flow solver; `io/` exporters;
  `walk/` the explorer.
* `Petri/examples/` : a few small nets to show capabilities and for quick tests.
  Benchmarks and MCC models do **not** live in this repo (see below).
* `Petri/test/` : test and bench scripts, logs.
* `PetriSage/` : python, out of scope, never touch.
* `README.md`, `KERS.md` : user docs (KERS is the binary sparse matrix format
  for program-to-program exchange). `WALK_PLAN.md` : the walk engine design and
  plan of attack; the reference until folder `algorithm.md` files take over.

Reference code bases, read-only from this project; never edit them from here:

* `~/git/ITStools` : the Java side. `fr.lip6.move.gal.structural.walk`
  (`RandomExplorer`, `WalkUtils`) is the ancestor of `Walker.h`; the property
  parsers and MCC output conventions live there too.
* `~/git/libHSC` : vendors copies of our `core/` files under
  `include/hsc/petri/`; its `surface/` is the s-expression syntax we may mirror.
* `~/git/PetriVizu/public/syntax.md` : the human-friendly property syntax.
* `~/git/pnmcc-models-20xx` : MCC model corpora (benchmark inputs, run on the
  cluster, never copied here).
* Branch `origin/er/link-spot` (2021): source of the `expr/` tree and the MCC
  property parser. Read it with `git show origin/er/link-spot:<path>`; do not
  check it out or merge it.

Every folder bears a `README.md` (source map, services) and/or `algorithm.md`
(abstract description of the algorithm). Always read the README before using or
editing files in a folder. Maintain these docs in sync with code.

## Discipline (mandatory)

- Do not spawn agents.
- We are solo: master only, no branches. Commit directly to master.
- One commit per file by preference when editing an existing core file; bulk
  commits are fine for one logical change, moved code, new datasets.
- Use `git commit -F -` with a quoted heredoc and a terse message. No backticks
  in a double-quoted `-m`.
- `git add/rm/mv` close to `git commit`; don't leave the index open. Remember
  `git mv` / `git rm` already populate the index: commit before adding more.
- **Pushing**: don't ask, don't do it. The user pushes.
- **Never rewrite history**: no reset/rebase/amend, even unpushed. Fix forward.
- **NO persistent "memory" files**: the user does not use them (not inspectable
  in the repo, they bloat unseen). Capture anything durable in documentation in
  the repo. Never edit this file unless prompted by the user.
- Work in the folder: avoid `/tmp` and the scratchpad. If a session crashes,
  all traces of work must be *in* the repo.
- Keep files roughly under 500 LOC, single responsibility, organised in folders.
- Comments are context free: describe the code in that file, its current
  state, not its callers, not its history (git has the history).
- Doc first for new code: write the folder `README.md` / `algorithm.md`, then
  transcribe it into code.
- Test often, but let the user guide how much; too many unit tests are churn on
  a fast-moving code base. Probes longer than about 5 lines go in `Petri/test/`
  as a file (rm after use if one-shot, else commit).
- Large model files (hundreds of MB PNML) are never committed.

## Working style (how the user wants me to operate)

- **Look around, then propose.** New directions start with a written design
  in the repo (markdown) that the user peruses and comments; code follows
  validation. Present intermediate results, do not start a new direction
  without user validation.
- **Diagnostics self-bound, about 15s per example.** A blown timeout is a
  finding, report it. One input per invocation; redirect long output to
  `Petri/test/logs/`, never pipe long runs to `tail`.
- **No manual process management.** Never `kill`/`pkill`, no `&`/`nohup`.
  For long self-terminating runs use Bash `run_in_background` and wait on
  completion.
- **Sparse is the substrate.** Transitions touch a handful of places out of
  tens of thousands; every data structure in the hot path is sparse
  (`SparseArray`, `MatrixCol`). Do not introduce dense per-place or
  per-transition work in the walk loop without measuring and asking.
- **Design for threads.** Hot state is thread-local; shared knowledge goes
  through an explicit, documented structure.
- **Honest results.** Distinguish what we failed at from what a reference tool
  cannot handle; never misrepresent data. We are scientists.
- C++23, no compiler-specific extensions beyond what the build already uses
  (`__uint128_t` for `petri128`). Explicit signatures, `T` template on the
  integer width kept end to end.
