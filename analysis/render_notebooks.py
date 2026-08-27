#!/usr/bin/env python3
"""Render the figure notebooks to self-contained HTML, in the style of the
comparative-splicing and simulation bundles.

Quarto reads YAML front matter from a *raw cell* at the top of an .ipynb, so
this inserts (or refreshes) that cell before rendering. The cell is inert --
raw cells are never executed -- and re-running this script is idempotent.

Quarto does NOT re-execute a notebook: it renders the outputs already stored in
the file. That is the point (the heavy cells stay run once) but it also means a
notebook whose outputs live only in an unsaved editor renders with no figures.
This script refuses to render such a notebook rather than emit a silently empty
page. Save the notebook in Jupyter/VSCode first.

    python3 analysis/render_notebooks.py            # render all
    python3 analysis/render_notebooks.py Figure4    # render one
"""
import json, os, shutil, subprocess, sys

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(BASE), 'docs')
AUTHOR = "Carlos F Buen Abad Najar"
DATE = "2026-08-26"

NOTEBOOKS = [
    ("Figure2/Figure2.ipynb",
     "Figure 2: differential unproductive splicing across human tissues"),
    ("Figure2/Figure2_heatmap.ipynb",
     "Figure 2c: unproductive splicing with concordant differential expression"),
    ("Figure4/Figure4.ipynb",
     "Figure 4: genetic basis of variation in unproductive splicing"),
    ("Figure5/Figure5_Python.ipynb",
     "Figure 5: unproductive splicing in Alzheimer's disease (Python panel)"),
    ("Figure5/Figure5_R.ipynb",
     "Figure 5: unproductive splicing in Alzheimer's disease (R panels)"),
]

FRONT_MATTER = """---
title: "{title}"
author: "{author}"
date: "{date}"
format:
  html:
    code-fold: true
    code-tools: true
    toc: true
    toc-depth: 3
    theme: cosmo
    embed-resources: true
execute:
  echo: true
  warning: false
---
"""


def set_front_matter(path, title):
    nb = json.load(open(path))
    body = FRONT_MATTER.format(title=title, author=AUTHOR, date=DATE)
    lines = body.split('\n')
    cell = {"cell_type": "raw", "metadata": {},
            "source": [l + '\n' for l in lines[:-1]]}
    first = nb['cells'][0] if nb['cells'] else None
    if first and first['cell_type'] == 'raw' and \
            ''.join(first['source']).lstrip().startswith('---'):
        nb['cells'][0] = cell          # refresh ours
    else:
        nb['cells'].insert(0, cell)
    json.dump(nb, open(path, 'w'), indent=1)


IMG = ('image/png', 'image/jpeg', 'image/svg+xml')


def output_census(path):
    """(code cells, cells with any output, cells with an image output)."""
    nb = json.load(open(path))
    code = [c for c in nb['cells'] if c['cell_type'] == 'code']
    any_out = sum(1 for c in code if c.get('outputs'))
    imgs = sum(1 for c in code
               if any(k in o.get('data', {}) for o in c.get('outputs', [])
                      for k in IMG))
    return len(code), any_out, imgs


def main(argv):
    quarto = shutil.which('quarto') or os.path.expanduser('~/.local/bin/quarto')
    if not os.path.exists(quarto):
        sys.exit("quarto not found; install it or put it on PATH")
    os.makedirs(OUT, exist_ok=True)

    todo = [(n, t) for n, t in NOTEBOOKS
            if not argv or any(a in n for a in argv)]
    if not todo:
        sys.exit(f"no notebook matches {argv}")

    bad = []
    for n, _ in todo:
        code, any_out, imgs = output_census(os.path.join(BASE, n))
        # A figure notebook that has been run and saved has an image in most of
        # its plotting cells. No images at all means the outputs are still only
        # in an unsaved editor.
        if imgs == 0 or any_out < code / 2:
            bad.append((n, code, any_out, imgs))
    if bad:
        print("REFUSING -- these notebooks have no stored figures, so they would\n"
              "render with prose and code but no plots. Save them in Jupyter or\n"
              "VSCode first (quarto renders stored outputs; it does not re-run).\n")
        for n, code, any_out, imgs in bad:
            print(f"    {n:34s} {code:3d} code cells, {any_out:3d} with output, "
                  f"{imgs:3d} with a figure")
        return 1

    for nbpath, title in todo:
        src = os.path.join(BASE, nbpath)
        code, any_out, imgs = output_census(src)
        print(f"--- {nbpath}  ({imgs} figures in {any_out}/{code} cells)")
        set_front_matter(src, title)
        subprocess.run([quarto, "render", os.path.basename(src), "--to", "html",
                        "--output-dir", os.path.relpath(OUT, os.path.dirname(src))],
                       cwd=os.path.dirname(src), check=True)

    print("\nrendered into", OUT)
    for f in sorted(os.listdir(OUT)):
        if f.endswith('.html'):
            mb = os.path.getsize(os.path.join(OUT, f)) / 1e6
            print(f"   {f:34s} {mb:6.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
