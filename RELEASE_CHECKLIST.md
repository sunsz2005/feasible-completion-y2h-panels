# Public-release checklist

The Git package is complete for checking recorded results. Before a full score replay or a versioned release:

1. The authors' original code has an MIT license. Confirm that any reused third-party code, if identified later, is excluded or separately credited; the MIT license does not apply to upstream datasets.
2. The YDIP and YBIOGRID snapshots are excluded from Git. Restore and hash-check them using `python tools/restore_yeast_inputs.py` before the full score replay. Do not add them to the public repository unless redistribution terms are confirmed.
3. Create a stable release and archive that release in a service that assigns a persistent DOI, if available.
4. Record the public repository URL, release tag, DOI or archive identifier, access conditions, and final license in the manuscript's **Availability of data and materials** statement. The current submission text still says that a public location has not been assigned.
5. Verify the uploaded tree using `python tools/verify_release.py`. For the full score replay, use Python 3.12.14 and the pinned dependencies described in `README.md`.

The archived lock files include historical absolute Windows workspace paths. These are provenance records, not paths a reader must recreate; the audit script rebases them. The data are interaction records and candidate identifiers, not participant-level clinical data. This checklist does not replace the journal's own submission declarations.

