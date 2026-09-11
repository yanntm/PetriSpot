"""Remove JAR signing metadata from an isolated flat/native product copy."""
from pathlib import Path
import sys
import zipfile


def strip(product: Path) -> None:
    for path in product.rglob("*.jar"):
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            signatures = {name for name in names if name.upper().startswith("META-INF/")
                          and name.upper().endswith((".SF", ".RSA", ".DSA", ".EC"))}
            if not signatures:
                continue
            temporary = path.with_suffix(".jar.unsigned")
            with zipfile.ZipFile(temporary, "w") as output:
                for entry in archive.infolist():
                    if entry.filename not in signatures:
                        output.writestr(entry, archive.read(entry))
        temporary.replace(path)
        print(path)


if __name__ == "__main__":
    strip(Path(sys.argv[1]))
