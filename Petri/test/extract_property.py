#! /usr/bin/env python3
"""Build an MCC model folder holding one property, for a reproduction.

    extract_property.py <model folder> <examination> <id or suffix> <out folder>
    extract_property.py bench/models/Foo-PT-01 LTLCardinality 03 bench/models/Foo-f03

Copies model.pnml and writes <examination>.xml with the single matching
property, in the MCC namespace the parsers expect. The suffix matches the end
of the property id, so "03" picks <model>-<examination>-03.
"""
import os, shutil, sys
import xml.etree.ElementTree as ET

MCC = "http://mcc.lip6.fr/"


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    src, exam, want, dst = sys.argv[1:]
    ET.register_namespace("", MCC)
    tree = ET.parse(os.path.join(src, exam + ".xml"))
    root = tree.getroot()
    props = root.findall(f"{{{MCC}}}property")
    keep = [p for p in props if p.findtext(f"{{{MCC}}}id", "").endswith(want)]
    if len(keep) != 1:
        sys.exit(f"{len(keep)} properties match {want!r} among {len(props)}: "
                 + ", ".join(p.findtext(f'{{{MCC}}}id', '') for p in keep))
    for p in props:
        if p is not keep[0]:
            root.remove(p)
    os.makedirs(dst, exist_ok=True)
    shutil.copy(os.path.join(src, "model.pnml"), os.path.join(dst, "model.pnml"))
    tree.write(os.path.join(dst, exam + ".xml"), xml_declaration=True, encoding="UTF-8")
    print(f"{dst}: {keep[0].findtext(f'{{{MCC}}}id')} of {len(props)}")


if __name__ == "__main__":
    main()
