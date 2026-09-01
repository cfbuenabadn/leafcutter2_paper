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

    python3 analysis/render_notebooks.py                  # render all
    python3 analysis/render_notebooks.py Figure4          # render one
    python3 analysis/render_notebooks.py --check          # report only, render nothing
    python3 analysis/render_notebooks.py --sync-external  # re-copy drifted external pages
"""
import json, os, shutil, subprocess, sys, tempfile

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(BASE), 'docs')
AUTHOR = "Carlos F Buen Abad Najar"
DATE = "2026-08-26"

# (path, page title, produces inline figures)
# The sashimi notebook prepares inputs for an external tool and draws nothing,
# so the "no figures means you forgot to save" check must not apply to it.
NOTEBOOKS = [
    ("Figure1/Figure1.ipynb",
     "Figure 1: classification of unproductive splice junctions", True),
    ("Figure1/Figure1_fig1i.ipynb",
     "Figure 1i: unproductive splicing vs intron length and expression", True),
    ("Figure2/Figure2.ipynb",
     "Figure 2: differential unproductive splicing across human tissues", True),
    ("Figure2/Figure2_heatmap.ipynb",
     "Figure 2c: unproductive splicing with concordant differential expression", True),
    ("Figure2/Figure2_prepare_sashimi.ipynb",
     "Figure 2d: GABBR1 sashimi plot inputs", False),
    ("Figure4/Figure4.ipynb",
     "Figure 4: genetic basis of variation in unproductive splicing", True),
    ("Figure5/Figure5_Python.ipynb",
     "Figure 5: unproductive splicing in Alzheimer's disease (Python panel)", True),
    ("Figure5/Figure5_R.ipynb",
     "Figure 5: unproductive splicing in Alzheimer's disease (R panels)", True),
]

# Pages copied in from the companion repositories. They are rendered there, not
# here, so the only check available is whether our copy still matches theirs.
_BF = '/project/yangili1/cfbuenabadn'
EXTERNAL = {
    'ComparativeSplicingFigures.html':
        f'{_BF}/20260825_comparativesplicing_paper/docs/ComparativeSplicingFigures.html',
    'NMD_GroupingDiscrepancy.html':
        f'{_BF}/20260825_comparativesplicing_paper/docs/NMD_GroupingDiscrepancy.html',
    'SimulationBenchmarkFigures.html':
        f'{_BF}/20260825_leaf2simulation_paper/docs/SimulationBenchmarkFigures.html',
}


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


def drop_duplicate_figures(nb):
    """Remove an execute_result image that duplicates a display_data image.

    A cell whose last expression returns a Figure/Axes stores the plot twice:
    once as execute_result (the returned object's repr) and once as
    display_data (the inline backend). Both carry the same PNG, so the
    rendered page embeds every panel twice. Dropping the execute_result copy
    changes nothing visible and halves the file.
    """
    dropped = 0
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        outs = c.get('outputs', [])
        shown = {o['data'][k] for o in outs if o.get('output_type') == 'display_data'
                 for k in IMG if k in o.get('data', {})}
        if not shown:
            continue
        keep = []
        for o in outs:
            if (o.get('output_type') == 'execute_result'
                    and any(o.get('data', {}).get(k) in shown for k in IMG)):
                dropped += 1
                continue
            keep.append(o)
        c['outputs'] = keep
    return dropped


def output_census(path):
    """(code cells, cells with any output, cells with an image output)."""
    nb = json.load(open(path))
    code = [c for c in nb['cells'] if c['cell_type'] == 'code']
    any_out = sum(1 for c in code if c.get('outputs'))
    imgs = sum(1 for c in code
               if any(k in o.get('data', {}) for o in c.get('outputs', [])
                      for k in IMG))
    return len(code), any_out, imgs



def stale_pages():
    """Notebooks whose rendered page is older than the notebook itself.

    Nothing forces docs/ to agree with the notebooks -- rendering is manual --
    so a page can sit wrong indefinitely. This is the check that catches it.
    """
    out = []
    for nbpath, _, _ in NOTEBOOKS:
        src = os.path.join(BASE, nbpath)
        html = os.path.join(OUT, os.path.basename(nbpath).replace('.ipynb', '.html'))
        if not os.path.exists(src):
            continue
        if not os.path.exists(html):
            out.append((nbpath, 'never rendered'))
        elif os.path.getmtime(src) > os.path.getmtime(html) + 5:
            out.append((nbpath, 'notebook is newer than its page'))
    return out


def diverged_external():
    """Copied-in pages that no longer match the repository they came from.

    Skipped silently when the source repository is not on this machine, since
    the copy is still perfectly serviceable without it.
    """
    out = []
    for name, origin in EXTERNAL.items():
        ours = os.path.join(OUT, name)
        if not os.path.exists(ours):
            out.append((name, 'missing from docs/'))
        elif not os.path.exists(origin):
            out.append((name, 'source repository not on this machine, cannot check'))
        elif open(ours, 'rb').read() != open(origin, 'rb').read():
            out.append((name, 'differs from the source repository'))
    return out


def report(sync=False):
    """Print both checks. With sync=True, re-copy any external page that drifted."""
    stale = stale_pages()
    print('\nrendered pages vs their notebooks:')
    if stale:
        for nb, why in stale:
            print(f'   STALE  {nb:<38} {why}')
        print('   fix: python3 analysis/render_notebooks.py')
    else:
        print(f'   all {len(NOTEBOOKS)} in sync')

    ext = diverged_external()
    print('copied-in pages vs their source repositories:')
    actionable = [(n, w) for n, w in ext if 'cannot check' not in w]
    if not ext:
        print(f'   all {len(EXTERNAL)} identical')
    for name, why in ext:
        tag = 'note ' if 'cannot check' in why else 'DIFF '
        print(f'   {tag} {name:<38} {why}')
        if sync and 'differs' in why:
            shutil.copy2(EXTERNAL[name], os.path.join(OUT, name))
            print(f'          re-copied from {EXTERNAL[name]}')
    if not actionable and ext:
        print('   (nothing actionable)')
    return len(stale) + len(actionable)


def main(argv):
    sync = '--sync-external' in argv
    if '--check' in argv:
        return 1 if report(sync=sync) else 0
    argv = [a for a in argv if not a.startswith('--')]

    quarto = shutil.which('quarto') or os.path.expanduser('~/.local/bin/quarto')
    if not os.path.exists(quarto):
        sys.exit("quarto not found; install it or put it on PATH")
    os.makedirs(OUT, exist_ok=True)

    todo = [e for e in NOTEBOOKS
            if not argv or any(a in e[0] for a in argv)]
    if not todo:
        sys.exit(f"no notebook matches {argv}")

    bad = []
    for n, _, wants_figs in todo:
        if not wants_figs:
            continue
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

    for nbpath, title, _ in todo:
        src = os.path.join(BASE, nbpath)
        code, any_out, imgs = output_census(src)
        print(f"--- {nbpath}  ({imgs} figures in {any_out}/{code} cells)")
        set_front_matter(src, title)

        # Render a de-duplicated copy so the saved notebook is left alone.
        with tempfile.TemporaryDirectory() as td:
            nb = json.load(open(src))
            n = drop_duplicate_figures(nb)
            if n:
                print(f"    dropped {n} duplicate figure output(s)")
            tmp = os.path.join(td, os.path.basename(src))
            json.dump(nb, open(tmp, 'w'), indent=1)
            subprocess.run([quarto, "render", os.path.basename(tmp), "--to", "html"],
                           cwd=td, check=True, stdout=subprocess.DEVNULL)
            html = tmp[:-len('.ipynb')] + '.html'
            shutil.copy2(html, os.path.join(OUT, os.path.basename(html)))

    report(sync=sync)
    print("\nrendered into", OUT)
    for f in sorted(os.listdir(OUT)):
        if f.endswith('.html'):
            mb = os.path.getsize(os.path.join(OUT, f)) / 1e6
            print(f"   {f:34s} {mb:6.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
