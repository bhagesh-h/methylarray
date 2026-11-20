#!/usr/bin/env python3
"""Assemble the sources of the documentation website (MkDocs, mkdocs.yml).

The pages stay single-sourced: README.md becomes the home page, docs/*.md the
guides, CITATIONS.md and CHANGELOG.md their own pages. Links between these
files are rewritten to the site's page names; links to any other file in the
repository (LICENSE, nextflow_schema.json, assets/...) point to GitHub.

Usage: python3 dev/docs/build_site.py [output folder, default .site/src]
Then:  mkdocs build --strict   (or: dev/nfdev docs)
"""
import os
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / ".site" / "src"
GITHUB = "https://github.com/bhagesh-h/methylarray/blob/main/"

# repository file -> page in the site
PAGES = {
    "README.md": "index.md",
    "CITATIONS.md": "citations.md",
    "CHANGELOG.md": "changelog.md",
    "docs/usage.md": "usage.md",
    "docs/output.md": "output.md",
    "docs/CONTRIBUTING.md": "contributing.md",
    "docs/developer.md": "developer.md",
    "docs/development/lint_exceptions.md": "development/lint_exceptions.md",
}
# files copied as they are (images)
ASSETS = {"docs/images/metro_map.svg": "images/metro_map.svg"}

LINK = re.compile(r"(\]\()([^)\s]+)(\))")


def rewrite(source: str, text: str) -> str:
    """Rewrite the relative links of one page, read from repository file `source`."""
    base = Path(source).parent
    page_dir = Path(PAGES[source]).parent

    def fix(match: re.Match) -> str:
        target = match.group(2)
        if re.match(r"^(https?:|mailto:|#)", target):
            return match.group(0)
        path, _, anchor = target.partition("#")
        resolved = os.path.normpath((base / path).as_posix()) if path else source
        suffix = f"#{anchor}" if anchor else ""
        page = PAGES.get(resolved) or ASSETS.get(resolved)
        if page:
            return f"{match.group(1)}{os.path.relpath(page, page_dir)}{suffix}{match.group(3)}"
        return f"{match.group(1)}{GITHUB}{resolved}{suffix}{match.group(3)}"

    return LINK.sub(fix, text)


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    for source, page in PAGES.items():
        target = OUT / page
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(rewrite(source, (REPO / source).read_text(encoding="utf-8")), encoding="utf-8")
    for source, asset in ASSETS.items():
        target = OUT / asset
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / source, target)
    print(f"{len(PAGES)} pages and {len(ASSETS)} images in {OUT}")


if __name__ == "__main__":
    main()
