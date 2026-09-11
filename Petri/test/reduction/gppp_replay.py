"""One-shot arbitrary-precision replay of the GPPP CTL counterexample."""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

root = ET.parse(Path(sys.argv[1]) / 'model.pnml').getroot()
for element in root.iter():
    element.tag = element.tag.rsplit('}', 1)[-1]
marking = {p.attrib['id']: int(p.findtext('initialMarking/text', '0'))
           for p in root.iter('place')}
names = {t.findtext('name/text', t.attrib['id']): t.attrib['id']
         for t in root.iter('transition')}
arcs = [(a.attrib['source'], a.attrib['target'], int(a.findtext('inscription/text', '1')))
        for a in root.iter('arc')]
trace = Path(sys.argv[2]).read_text().split('EU: path to the right side : ', 1)[1].splitlines()[0].split()
for label, path in [('C++ evidence', trace)]:
    m = marking.copy()
    assert m['_3PG'] < 2467342475
    largest = max(m.values())
    for step, name in enumerate(path):
        t = names[name]
        pre = [(s, w) for s, target, w in arcs if target == t]
        post = [(target, w) for s, target, w in arcs if s == t]
        assert all(m[p] >= w for p, w in pre), (label, step, name, pre)
        for p, w in pre:
            m[p] -= w
        for p, w in post:
            m[p] += w
        largest = max(largest, max(m.values()))
    assert m['Ru5P'] < m['Pyr']
    assert m['E4P'] < 228211673
    assert m['FBP'] < m['b2'] and m['GAP'] < 2707223165
    print(label, len(path), 'enabled firings; maximum marking', largest)
    print({p: m[p] for p in ['Ru5P', 'Pyr', 'E4P', 'FBP', 'b2', 'GAP']})
    print('AG body false at endpoint; outer until left side false initially: property FALSE')
