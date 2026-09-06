# Web report of a campaign

Superseded by `~/git/MCC-analysis/campaign/` (read its `README.md`): local
pages built from result sets, ours from the harness logs in
`/data/ythierry/MCC26run/`, the contest's from `raw-result-analysis.csv`,
cross-referenced against the field and against each other, served on
localhost by `serve.py` and read through an SSH tunnel. The markdown
`report.py` here stays as the plain text reading of a campaign folder; the
total examinations (`total-runs.csv`) are not in the pages yet.

The first design, one self-contained page per campaign folder with
DataTables and Plotly, is what the prototype does; the ideas kept for later
are the closing curves (t25..t100 per run) and the family and engine share
views of the total examinations, and a baseline mode over two campaign
folders, which the pairwise comparison of result sets now covers.
