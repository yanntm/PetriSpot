# `test/probes/` — one-off analysis scripts kept for reference

* `project_places.py` — over-approximating projection of a PNML net onto a
  set of places (regex): the other places are removed with their arcs, so
  every firing sequence of the original stays fireable and the reachable
  markings of the kept places are a superset. Any bound proved on the
  projection holds on the original. Writes `model.pnml` and an
  `UpperBounds.xml` with one `place-bound` per kept place.

  Finding (2026-09-04), BridgeAndVehicles-PT-V80P50N10 (identical to the COL
  instance): the semi-flow `SUR_PONT_A + CAPACITE + SUR_PONT_B = 50` gives 50
  for the bridge places while the walk reaches 10 and never more; the two
  UpperBounds on `SUR_PONT_B` are the ones no tool answered in 2026. The
  counter is one-hot (`COMPTEUR_0..10`): a decision increments it before each
  admission, no decision exists at 10, the switch goes through a drain phase
  whose exit (`basculement`) reads `CAPACITE >= 50` (bridge empty). Projecting
  on `CAPACITE, SUR_PONT_*, CONTROLEUR_*, CHOIX_*, VIDANGE_*, COMPTEUR_*`
  (20 of 188 places, 50 of 2108 transitions after merging) gives a net with
  264 reachable states on which ITS-Tools proves `SUR_PONT_A = SUR_PONT_B = 10`
  exhaustively; with the walk's 10 on the original, the bound is exactly 10.

  ```
  python3 Petri/test/probes/project_places.py bench/models/BridgeAndVehicles-PT-V80P50N10/model.pnml \
      bench/models/BridgeSlice-V80P50N10 '^(CAPACITE|SUR_PONT_[AB]|CONTROLEUR_[12]|CHOIX_[12]|VIDANGE_[12]|COMPTEUR_[0-9]+)$'
  its-tools -pnfolder bench/models/BridgeSlice-V80P50N10 -examination UpperBounds -timeout 120
  ```
