# Guide: the CI chain

Two repositories, two CIs, one dependency between them. Everything the harness
installs is a CI artifact, never a local build (`BENCH.md`), so a change is
benchmarked only once it has travelled this chain.

## PetriSpot

`.github/workflows/{linux,osx,windows}.yml`, on every push to `master`.
`linux.yml` runs `buildPetriSpot.sh` (static libexpat, stripped
`petri32/64/128` and `kersconv` into `website/`) and deploys `website/` to the
branch `Inv-Linux` (`Inv-OSX`, `Inv-Windows` likewise). The binary is then at

```
https://github.com/yanntm/PetriSpot/raw/Inv-Linux/petri64
```

Check what is published:

```
git fetch origin Inv-Linux && git log -1 --format='%h %ci %s' origin/Inv-Linux   # "Deploying to Inv-Linux from @ yanntm/PetriSpot@<sha>"
wget -q https://github.com/yanntm/PetriSpot/raw/Inv-Linux/petri64 -O /tmp/p && chmod +x /tmp/p && /tmp/p -h | grep -c ctlSteps
```

Run status without `gh`:

```
curl -s "https://api.github.com/repos/yanntm/PetriSpot/actions/runs?per_page=4" | python3 -c "
import sys,json
for w in json.load(sys.stdin)['workflow_runs']: print(w['created_at'], w['status'], w['conclusion'], w['head_sha'][:8], w['name'])"
```

## ITS-Tools

`~/git/ITStools`, origin `lip6/ITStools`, `.github/workflows/build.yml` on
push to `master`: Maven with Tycho from `fr.lip6.move.gal.parent`, the
command line product, the GraalVM native image built from it
(`ITS-commandline/native/build-native.sh`), and `website/` deployed to
`gh-pages`. Published files:

```
https://lip6.github.io/ITSTools/fr.lip6.move.gal.itscl.product-linux.gtk.x86_64.zip
https://lip6.github.io/ITSTools/its-tools-native
```

The build of the plugin `fr.lip6.petrispot.binaries` **downloads `petri64`
from `Inv-Linux` at build time** (Maven `get` task; the same for Spot's
`spotutil`, LTSmin, GreatSPN). Hence the order: push PetriSpot, wait for
its Linux build to deploy, then push ITS-Tools. A product built before the
deploy bundles the previous binary, and a new Java flag meets an old binary:
CLI11 rejects it and every walk dies silently ("PetriSpot walker: 0/n solved,
exit 109"). Check the pair:

```
curl -sIL https://lip6.github.io/ITSTools/fr.lip6.move.gal.itscl.product-linux.gtk.x86_64.zip | grep -i last-modified
cd /data/ythierry/MCC26deploy/products/itstools-ci-check && rm -rf * && wget -q https://lip6.github.io/ITSTools/fr.lip6.move.gal.itscl.product-linux.gtk.x86_64.zip && unzip -q *.zip
sha256sum plugins/fr.lip6.petrispot.binaries_*/bin/petri64 /tmp/p     # same hash as the Inv-Linux binary
ls plugins/fr.lip6.move.gal.application.pnmcc_*.jar                     # the product's build stamp, yyyyMMddHHmm
```

Re-triggering a CI run without `gh` or a token: push a commit (a doc line
does). The ITS-Tools run takes about 10 minutes, the pages deployment a few
more.

## The harness picks it up

`MCC-drivers/install_itstools.sh` clones `yanntm/ITS-Tools-MCC` into
`itstools/`, whose own `install_itstools.sh` downloads the product zip and the
native image (dropped when the URL fails); `MCC-drivers/petrispot/install.sh`
downloads `petri64` from `Inv-Linux` (`PETRISPOT_BIN=<path>` installs a
local build instead). Then the deploy sequence of `docs/CLUSTER.md`.

## MCC-analysis pages

`~/git/MCC-analysis` builds the result pages from the collected `csv/`
folders; its `campaign/example.json` names the sets (`ITS-Tools latest`, the
dated ones); pushed pages appear on its gh-pages.
