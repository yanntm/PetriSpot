#!/usr/bin/env python3
"""Project a P/T PNML onto a set of places: every other place is removed with
its arcs, transitions left without arcs are dropped, transitions that became
identical are merged. Removing places only removes constraints and effects,
so the result over-approximates the reachable markings of the kept places:
any bound proved on the projection holds on the original.

Usage: project_places.py model.pnml out_dir PLACE_REGEX
Writes out_dir/model.pnml and out_dir/UpperBounds.xml (one place-bound per kept place)."""
import re, sys, os, html
src, out, rx = sys.argv[1], sys.argv[2], re.compile(sys.argv[3])
x = open(src).read()
places = {}
for m in re.finditer(r'<place id="([^"]*)">(.*?)</place>', x, re.S):
    im = re.search(r'<initialMarking>.*?<text>(\d+)</text>', m.group(2), re.S)
    places[m.group(1)] = int(im.group(1)) if im else 0
kept = [p for p in places if rx.search(p)]
keptset = set(kept)
def w(body):
    m = re.search(r'<inscription>.*?<text>(\d+)</text>', body, re.S); return int(m.group(1)) if m else 1
trans = {}  # name -> (pre dict, post dict)
for t in re.findall(r'<transition id="([^"]*)"', x): trans[t] = ({}, {})
for s, t, body in re.findall(r'<arc id="[^"]*" source="([^"]*)" target="([^"]*)"(.*?)</arc>', x, re.S):
    if s in places:
        if s in keptset: trans[t][0][s] = w(body)
    elif t in keptset:
        trans[s][1][t] = w(body)
merged = {}
for t, (pre, post) in trans.items():
    if not pre and not post: continue
    key = (tuple(sorted(pre.items())), tuple(sorted(post.items())))
    merged.setdefault(key, t)
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'model.pnml'), 'w') as f:
    f.write('<?xml version="1.0"?>\n<pnml xmlns="http://www.pnml.org/version-2009/grammar/pnml"><net id="slice" type="http://www.pnml.org/version-2009/grammar/ptnet"><page id="page0">\n')
    for p in kept:
        f.write('<place id="%s"><name><text>%s</text></name>' % (p, p))
        if places[p]: f.write('<initialMarking><text>%d</text></initialMarking>' % places[p])
        f.write('</place>\n')
    n = 0
    for (pre, post), t in merged.items():
        f.write('<transition id="%s"><name><text>%s</text></name></transition>\n' % (t, t))
        for p, v in pre:
            f.write('<arc id="a%d" source="%s" target="%s"><inscription><text>%d</text></inscription></arc>\n' % (n, p, t, v)); n += 1
        for p, v in post:
            f.write('<arc id="a%d" source="%s" target="%s"><inscription><text>%d</text></inscription></arc>\n' % (n, t, p, v)); n += 1
    f.write('</page></net></pnml>\n')
with open(os.path.join(out, 'UpperBounds.xml'), 'w') as f:
    f.write('<?xml version="1.0"?>\n<property-set xmlns="http://mcc.lip6.fr/">\n')
    for i, p in enumerate(kept):
        f.write('<property><id>slice-UpperBounds-%02d-%s</id><description>max %s</description><formula><place-bound><place>%s</place></place-bound></formula></property>\n' % (i, html.escape(p), p, p))
    f.write('</property-set>\n')
print("kept %d of %d places, %d of %d transitions after merging" % (len(kept), len(places), len(merged), len(trans)))
