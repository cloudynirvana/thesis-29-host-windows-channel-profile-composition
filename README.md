# Host infection windows that re-rank channel–profile composition without merging into Θ

**Thesis #29.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #22 (observation-channel Disease Profile composition) and Thesis #23 (host-ranked immunometabolic hypotheses).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Can declared host infection windows that re-rank immunometabolic hypotheses also change which named observation-channel compositions remain legal Disease Profiles, without promoting either the windows or the channel coefficients into Θ?

They can. Sixteen composition files meet one structural rule list and one schedule clause. The structural list is the disjoint-union contract. It does not read the schedule. Seven files pass it. On the ledger schedule the accept set has five files. On the host schedule, which adds the declared windows [1.5, 2.5] and [10, 12], the accept set has one file, `windowed`. Four files leave. None enters. The JSON Schema accept set has eleven files on both schedules. SHA-256 of the canonical kinetic object is `cf5be3430e25f75d9f87810c789195a3225388b0a9973de7970ec3890a70313e` at every checkpoint, including after a refused promotion that would have written `k_inf = 1/2` and `k_host = 1/9` into Θ.

The rank control still moves first place from U2 to U3 when the windows are attached. That order was recomputed here. It was not read from Thesis #23's results file. The Fisher diagonals of Thesis #22 were not recomputed.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `schema/windowed_composition.schema.json` | JSON Schema `1.3.0-windowed-compose` (schedule-blind) |
| `profiles/*.profile.yaml` | The sixteen composition files |
| `sim/windowed_compose.py` | Seeded accept/refuse sets and the rank control (seed 20260929) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/crossref_snapshot.json` | DOI snapshot consulted on 21 September 2026 |
| `sim/figures/` | Accept/refuse matrix, set movement, rank control |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/windowed_compose.py
python3 build_pdf.py
```

NumPy, Matplotlib, PyYAML, and jsonschema are required for the catalogue. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json`, the profile files, and `sim/figures/`.

## Cite

Ogbonna KE. Host infection windows that re-rank channel–profile composition without merging into Θ [Internet]. Thesis #29 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-29-host-windows-channel-profile-composition

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

Depends on: [Thesis #22](https://github.com/cloudynirvana/thesis-22-observation-channel-profile-composition), [Thesis #23](https://github.com/cloudynirvana/thesis-23-host-ranked-immunometabolic-hypotheses).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
