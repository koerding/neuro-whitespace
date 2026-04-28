#!/usr/bin/env python3
"""
Brain region PubMed prevalence analysis.

For each region in REGIONS, this script computes two numbers from PubMed:

  inclusive_count
      Papers in 2015-2025 whose title/abstract matches the region.
      Match = any canonical phrase ("hippocampus") OR any disambiguated
      abbreviation (e.g. "VTA"[tiab] AND tegmental[tiab]).

  exclusive_count
      Inclusive count NOT (term of any descendant region).
      Only meaningful for umbrella terms. Example: exclusive
      "prefrontal cortex" excludes papers that focus on mPFC, dlPFC,
      vmPFC, or OFC. The inclusive/exclusive gap shows how much of an
      umbrella's literature is really about a named subdivision.

Why "inclusive" doesn't quite mean what it says: the regional term is
matched in title/abstract only, so a paper that names a subregion in the
abstract but never the umbrella will count for the subregion but not the
umbrella. Inclusive is "papers that mention this label", not "all papers
about anything within this region".

Inputs:
    REGIONS  hard-coded below. Each entry is a tuple
             (name, category, parent, phrases, ambigs):
               phrases   list[str]     unambiguous phrases, searched
                                       as "phrase"[tiab]
               ambigs    list[(str,str)]  (abbrev, ctx_OR_string).
                                       The abbrev is searched as
                                       ("ABBREV"[tiab] AND ctx[tiab]).
               parent    str | None    name of the broader region for
                                       the hierarchy used to compute
                                       exclusive counts.

Outputs (written to --out-dir):
    queries.json     every constructed query, saved BEFORE any fetch
                     so we have a record even if NCBI calls fail
    regions.csv      flat table: queries, counts, URLs, log values
    regions.json     same data as JSON
    prevalence.png   category-grouped lollipop chart on a log axis

Usage:
    python3 run.py --email YOU@example.com
    python3 run.py --dry-run         # print queries; no NCBI calls
    python3 run.py --api-key XYZ     # higher rate limit if you have one

Caveats: see README.md. Counts are a bibliometric proxy for term
prevalence, not anatomical ground truth or scientific importance.
"""

from __future__ import annotations
import argparse, csv, datetime as dt, json, math, sys, time
import urllib.parse, urllib.request, urllib.error
from collections import defaultdict
from pathlib import Path

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
DATE_FILTER = '("2015"[PDAT] : "2025"[PDAT])'
TOOL_NAME = "brain_region_prevalence"

# Hippocampal-subfield anchors used to disambiguate bare CA1 / CA2 / CA3.
# A naive "CAn"[tiab] AND hippocamp* leaks calcium-signaling papers,
# because PubMed normalizes Ca2+, Ca(2+) tokens to overlap with CA2.
# Phrase-only matching is the other extreme and undercounts legitimate
# subfield papers that say "Ca2+ imaging in CA2 of hippocampus" without
# matching one of our exact phrases. The compromise: require the bare
# abbreviation AND at least one piece of subfield-architecture vocabulary
# that essentially never appears in unrelated calcium-signaling papers
# (specific layers, projections, fibre systems, place cells, social
# memory). Verified empirically: this lifts CA2 from 212 (phrase-only)
# to ~1,000, vs. ~3,100 for the contaminated bare-AND-hippocamp* form.
CA_SUBFIELD_ANCHORS = (
    '"place cell" OR "place cells" OR "social memory" '
    'OR "pyramidal layer" OR "stratum pyramidale" '
    'OR "stratum radiatum" OR "stratum oriens" '
    'OR "dentate gyrus" OR "schaffer collateral" OR "schaffer collaterals" '
    'OR "mossy fiber" OR "mossy fibers" OR "perforant path" '
    'OR "hippocampal subfield" OR "hippocampal subfields"'
)

# Each region: (name, category, parent, phrases, ambigs)
#   phrases : list[str]  unambiguous phrases, searched as "phrase"[tiab]
#   ambigs  : list[(abbr, ctx_OR)]  abbreviation searched as
#             ("ABBR"[tiab] AND (ctx_term1[tiab] OR ctx_term2[tiab] ...))
#   parent  : name of broader region, or None
REGIONS = [
    # ---- Frontal / prefrontal cortex ----
    # Hierarchy: frontal cortex > {prefrontal cortex, motor cortex};
    # prefrontal cortex > {mPFC, dlPFC, vmPFC, OFC};
    # motor cortex > {premotor cortex, supplementary motor area}.
    ("prefrontal cortex", "frontal cortex", "frontal cortex",
        ["prefrontal cortex"],
        [("PFC", "prefrontal")]),
    ("medial prefrontal cortex", "frontal cortex", "prefrontal cortex",
        ["medial prefrontal cortex"],
        [("mPFC", "prefrontal")]),
    ("dorsolateral prefrontal cortex", "frontal cortex", "prefrontal cortex",
        ["dorsolateral prefrontal cortex"],
        [("dlPFC", "prefrontal"), ("DLPFC", "prefrontal")]),
    ("ventromedial prefrontal cortex", "frontal cortex", "prefrontal cortex",
        ["ventromedial prefrontal cortex"],
        [("vmPFC", "prefrontal"), ("VMPFC", "prefrontal")]),
    ("orbitofrontal cortex", "frontal cortex", "prefrontal cortex",
        ["orbitofrontal cortex", "orbital frontal cortex"], []),
    ("frontal cortex", "frontal cortex", None, ["frontal cortex"], []),
    ("motor cortex", "frontal cortex", "frontal cortex",
        ["motor cortex", "primary motor cortex"], []),
    ("premotor cortex", "frontal cortex", "motor cortex",
        ["premotor cortex"], []),
    ("supplementary motor area", "frontal cortex", "motor cortex",
        ["supplementary motor area"], []),

    # ---- Sensory / parietal / temporal / occipital ----
    # Hierarchy: occipital cortex > visual cortex > primary visual cortex.
    ("somatosensory cortex", "sensory cortex", None,
        ["somatosensory cortex", "primary somatosensory cortex"], []),
    ("parietal cortex", "sensory cortex", None, ["parietal cortex"], []),
    ("posterior parietal cortex", "sensory cortex", "parietal cortex",
        ["posterior parietal cortex"], [("PPC", "parietal")]),
    ("auditory cortex", "sensory cortex", None,
        ["auditory cortex", "primary auditory cortex"], []),
    ("visual cortex", "sensory cortex", "occipital cortex",
        ["visual cortex"], []),
    ("primary visual cortex", "sensory cortex", "visual cortex",
        ["primary visual cortex", "striate cortex"], []),
    ("occipital cortex", "sensory cortex", None, ["occipital cortex"], []),
    ("temporal cortex", "sensory cortex", None, ["temporal cortex"], []),
    ("inferior temporal cortex", "sensory cortex", "temporal cortex",
        ["inferior temporal cortex", "inferotemporal cortex"], []),
    ("superior temporal gyrus", "sensory cortex", "temporal cortex",
        ["superior temporal gyrus"], []),
    ("fusiform gyrus", "sensory cortex", None, ["fusiform gyrus"], []),
    ("angular gyrus", "sensory cortex", None, ["angular gyrus"], []),
    ("supramarginal gyrus", "sensory cortex", None, ["supramarginal gyrus"], []),

    # ---- Limbic / cingulate / insular cortex ----
    # Hierarchy: cingulate cortex > {ACC, PCC}.
    ("anterior cingulate cortex", "limbic cortex", "cingulate cortex",
        ["anterior cingulate cortex"], [("ACC", "cingulate")]),
    ("posterior cingulate cortex", "limbic cortex", "cingulate cortex",
        ["posterior cingulate cortex"], [("PCC", "cingulate")]),
    ("cingulate cortex", "limbic cortex", None, ["cingulate cortex"], []),
    ("retrosplenial cortex", "limbic cortex", None,
        ["retrosplenial cortex"], [("RSC", "retrosplenial")]),
    ("insular cortex", "limbic cortex", None, ["insular cortex", "insula"], []),
    ("entorhinal cortex", "limbic cortex", None,
        ["entorhinal cortex"], [("EC", "entorhinal")]),
    ("perirhinal cortex", "limbic cortex", None, ["perirhinal cortex"], []),
    ("parahippocampal cortex", "limbic cortex", None, ["parahippocampal cortex"], []),
    ("piriform cortex", "limbic cortex", None, ["piriform cortex"], []),

    # ---- Hippocampal formation ----
    # CA1 / CA2 / CA3: phrase forms PLUS the bare abbreviation guarded by
    # subfield-architecture anchors (see CA_SUBFIELD_ANCHORS). The bare
    # term alone leaks Ca2+ papers; phrase-only undercounts legitimate
    # papers that say "Ca2+ imaging in CA2 of mouse hippocampus" without
    # matching one of our canned phrases.
    ("hippocampus", "hippocampal", None,
        ["hippocampus", "hippocampal formation"], []),
    ("dentate gyrus", "hippocampal", "hippocampus", ["dentate gyrus"], []),
    ("CA1", "hippocampal", "hippocampus",
        ["hippocampal CA1", "CA1 region", "CA1 pyramidal",
         "CA1 neurons", "CA1 of the hippocampus", "CA1 area", "CA1 subfield"],
        [("CA1", CA_SUBFIELD_ANCHORS)]),
    ("CA2", "hippocampal", "hippocampus",
        ["hippocampal CA2", "CA2 region", "CA2 of the hippocampus",
         "CA2 pyramidal", "CA2 neurons", "CA2 area", "CA2 subfield"],
        [("CA2", CA_SUBFIELD_ANCHORS)]),
    ("CA3", "hippocampal", "hippocampus",
        ["hippocampal CA3", "CA3 region", "CA3 pyramidal",
         "CA3 neurons", "CA3 of the hippocampus", "CA3 area", "CA3 subfield"],
        [("CA3", CA_SUBFIELD_ANCHORS)]),
    ("subiculum", "hippocampal", "hippocampus", ["subiculum"], []),

    # ---- Amygdala ----
    ("amygdala", "amygdala", None, ["amygdala"], []),
    ("basolateral amygdala", "amygdala", "amygdala",
        ["basolateral amygdala"], [("BLA", "amygdala OR basolateral")]),
    ("central amygdala", "amygdala", "amygdala",
        ["central amygdala", "central nucleus of the amygdala"],
        [("CeA", "amygdala OR central")]),

    # ---- Septum / BNST / basal forebrain / claustrum ----
    ("bed nucleus of the stria terminalis", "limbic subcortical", None,
        ["bed nucleus of the stria terminalis"],
        [("BNST", "stria OR terminalis OR amygdala")]),
    ("lateral septum", "limbic subcortical", None, ["lateral septum"], []),
    ("medial septum", "limbic subcortical", None, ["medial septum"], []),
    ("septal nuclei", "limbic subcortical", None,
        ["septal nuclei", "septal nucleus"], []),
    ("basal forebrain", "limbic subcortical", None, ["basal forebrain"], []),
    ("diagonal band", "limbic subcortical", None,
        ["diagonal band of Broca", "nucleus of the diagonal band"], []),
    ("claustrum", "limbic subcortical", None, ["claustrum"], []),

    # ---- Basal ganglia ----
    ("striatum", "basal ganglia", None, ["striatum", "corpus striatum"], []),
    ("dorsal striatum", "basal ganglia", "striatum", ["dorsal striatum"], []),
    ("ventral striatum", "basal ganglia", "striatum", ["ventral striatum"], []),
    ("nucleus accumbens", "basal ganglia", "striatum",
        ["nucleus accumbens"], [("NAc", "accumbens OR striat* OR reward")]),
    ("caudate nucleus", "basal ganglia", "striatum", ["caudate nucleus"], []),
    ("putamen", "basal ganglia", "striatum", ["putamen"], []),
    ("globus pallidus", "basal ganglia", None,
        ["globus pallidus"],
        [("GPi", "pallidus OR pallidal"), ("GPe", "pallidus OR pallidal")]),
    ("ventral pallidum", "basal ganglia", None, ["ventral pallidum"], []),
    ("subthalamic nucleus", "basal ganglia", None,
        ["subthalamic nucleus"], [("STN", "subthalamic OR basal ganglia")]),

    # ---- Thalamus ----
    ("thalamus", "thalamus", None, ["thalamus"], []),
    ("mediodorsal thalamus", "thalamus", "thalamus",
        ["mediodorsal thalamus", "mediodorsal nucleus of the thalamus",
         "mediodorsal thalamic nucleus"], []),
    ("anterior thalamic nuclei", "thalamus", "thalamus",
        ["anterior thalamic nuclei", "anterior thalamic nucleus"], []),
    ("pulvinar", "thalamus", "thalamus", ["pulvinar"], []),
    ("lateral geniculate nucleus", "thalamus", "thalamus",
        ["lateral geniculate nucleus"],
        [("LGN", "geniculate OR visual OR thalam*")]),
    ("medial geniculate nucleus", "thalamus", "thalamus",
        ["medial geniculate nucleus"],
        [("MGN", "geniculate OR auditory OR thalam*")]),
    ("reticular thalamic nucleus", "thalamus", "thalamus",
        ["reticular thalamic nucleus", "thalamic reticular nucleus"],
        [("TRN", "thalam* OR reticular")]),

    # ---- Hypothalamus ----
    ("hypothalamus", "hypothalamus", None, ["hypothalamus"], []),
    ("arcuate nucleus", "hypothalamus", "hypothalamus", ["arcuate nucleus"], []),
    ("paraventricular nucleus", "hypothalamus", "hypothalamus",
        ["paraventricular nucleus"],
        [("PVN", "hypothalam* OR paraventricular")]),
    ("suprachiasmatic nucleus", "hypothalamus", "hypothalamus",
        ["suprachiasmatic nucleus"],
        [("SCN", "circadian OR suprachiasmatic")]),
    ("mammillary bodies", "hypothalamus", "hypothalamus",
        ["mammillary bodies", "mammillary body", "mamillary bodies"], []),

    # ---- Epithalamus / pineal ----
    ("habenula", "epithalamus", None, ["habenula", "habenular nuclei"], []),
    ("pineal gland", "epithalamus", None, ["pineal gland"], []),

    # ---- Midbrain ----
    # midbrain is itself a child of brainstem, so brainstem's exclusive
    # count transitively excludes substantia nigra, VTA, etc.
    ("midbrain", "midbrain", "brainstem", ["midbrain", "mesencephalon"], []),
    ("substantia nigra", "midbrain", "midbrain",
        ["substantia nigra"],
        [("SNc", "nigra OR dopamin*"), ("SNr", "nigra OR reticulata")]),
    ("ventral tegmental area", "midbrain", "midbrain",
        ["ventral tegmental area"],
        [("VTA", "tegmental OR midbrain OR dopamin* OR mesolimbic")]),
    ("superior colliculus", "midbrain", "midbrain",
        ["superior colliculus"], [("SC", "colliculus OR tectum")]),
    ("inferior colliculus", "midbrain", "midbrain", ["inferior colliculus"], []),
    ("red nucleus", "midbrain", "midbrain", ["red nucleus", "nucleus ruber"], []),
    ("periaqueductal gray", "midbrain", "midbrain",
        ["periaqueductal gray", "periaqueductal grey"],
        [("PAG", "periaqueductal")]),
    ("interpeduncular nucleus", "midbrain", "midbrain",
        ["interpeduncular nucleus"], []),

    # ---- Pons / medulla / brainstem ----
    ("brainstem", "brainstem", None, ["brainstem", "brain stem"], []),
    ("pons", "brainstem", "brainstem", ["pons"], []),
    ("pontine nuclei", "brainstem", "pons",
        ["pontine nuclei", "pontine nucleus"], []),
    ("medulla oblongata", "brainstem", "brainstem",
        ["medulla oblongata"], []),
    ("nucleus tractus solitarius", "brainstem", "medulla oblongata",
        ["nucleus tractus solitarius", "nucleus of the solitary tract"],
        [("NTS", "solitar* OR brainstem")]),
    ("area postrema", "brainstem", "medulla oblongata", ["area postrema"], []),
    ("parabrachial nucleus", "brainstem", "brainstem",
        ["parabrachial nucleus"], []),
    ("raphe nuclei", "brainstem", "brainstem",
        ["raphe nuclei", "raphe nucleus"], []),
    ("dorsal raphe nucleus", "brainstem", "raphe nuclei",
        ["dorsal raphe nucleus", "dorsal raphe"], []),
    ("locus coeruleus", "brainstem", "brainstem",
        ["locus coeruleus", "locus ceruleus"], []),

    # ---- Cerebellum ----
    ("cerebellum", "cerebellum", None, ["cerebellum"], []),
    ("cerebellar cortex", "cerebellum", "cerebellum", ["cerebellar cortex"], []),
    ("cerebellar vermis", "cerebellum", "cerebellum", ["cerebellar vermis"], []),
    ("cerebellar nuclei", "cerebellum", "cerebellum",
        ["cerebellar nuclei", "deep cerebellar nuclei"], []),
    ("dentate nucleus (cerebellar)", "cerebellum", "cerebellar nuclei",
        ["cerebellar dentate nucleus", "dentate nucleus of the cerebellum"], []),
    ("fastigial nucleus", "cerebellum", "cerebellar nuclei",
        ["fastigial nucleus"], []),
    ("interposed nucleus", "cerebellum", "cerebellar nuclei",
        ["interposed nucleus", "nucleus interpositus"], []),
    ("flocculus", "cerebellum", "cerebellum", ["flocculus"], []),

    # ---- Olfactory ----
    ("olfactory bulb", "olfactory", None, ["olfactory bulb"], []),
    ("olfactory tubercle", "olfactory", None, ["olfactory tubercle"], []),
]


# ---------- Query construction ----------

def context_to_query(ctx: str) -> str:
    """Wrap an 'A OR B OR C' shorthand into a tagged PubMed expression.

    Each context term is tagged with [tiab] so PubMed restricts the
    match to titles and abstracts. Wildcards (e.g. 'hippocamp*') flow
    through unchanged.

    >>> context_to_query("cortex OR cortical OR prefrontal")
    '(cortex[tiab] OR cortical[tiab] OR prefrontal[tiab])'
    """
    terms = [t.strip() for t in ctx.split(" OR ") if t.strip()]
    return "(" + " OR ".join(f"{t}[tiab]" for t in terms) + ")"


def region_term(region) -> str:
    """Build the OR-joined PubMed term for one region.

    No date filter is applied here — that gets ANDed in by build_queries.
    The output is a single parenthesized expression suitable for use
    inside a larger query (e.g. as the LHS of a NOT clause).

    For each canonical phrase: ``"phrase"[tiab]``.
    For each ambiguous abbreviation: ``("ABBR"[tiab] AND ctx[tiab])`` so
    the abbreviation only counts when accompanied by a domain term.
    """
    _, _, _, phrases, ambigs = region
    parts = [f'"{p}"[tiab]' for p in phrases]
    for abbr, ctx in ambigs:
        if ctx:
            parts.append(f'("{abbr}"[tiab] AND {context_to_query(ctx)})')
        else:
            parts.append(f'"{abbr}"[tiab]')
    return "(" + " OR ".join(parts) + ")"


def descendants(name: str, by_parent: dict) -> list:
    """Return every transitive child of `name`, depth-first.

    `by_parent` maps a parent region name to the list of region tuples
    whose `parent` field equals that name. Used to assemble the NOT
    clause for an umbrella term: every descendant must be subtracted,
    not just direct children, so brainstem's exclusive count also
    excludes substantia nigra (a grandchild via midbrain).

    A `visited` set guards against cycles created by accidental
    parent/child loops in REGIONS — without it a typo could make the
    function spin forever.
    """
    out, stack = [], list(by_parent.get(name, []))
    visited = {name}
    while stack:
        r = stack.pop()
        if r[0] in visited:
            continue
        visited.add(r[0])
        out.append(r)
        stack.extend(by_parent.get(r[0], []))
    return out


def validate_regions(regions):
    """Sanity-check the REGIONS table at startup.

    Catches the two failure modes that produce silent, wrong results
    instead of an obvious crash: duplicated region names (which collide
    in the by_parent dict), and parent strings that don't reference any
    real region (which silently disable the exclusive logic for that
    branch).
    """
    names = [r[0] for r in regions]
    seen = set()
    dups = []
    for n in names:
        if n in seen:
            dups.append(n)
        seen.add(n)
    if dups:
        raise ValueError(f"duplicate region names: {dups}")
    nameset = set(names)
    bad = [(r[0], r[2]) for r in regions if r[2] is not None and r[2] not in nameset]
    if bad:
        raise ValueError(f"unknown parent references: {bad}")


def build_queries(regions):
    """Construct inclusive and exclusive PubMed query strings per region.

    Returns a list of dicts (no counts yet) so we can dump
    ``queries.json`` before any network call. If the live fetch later
    blows up, we still have a complete record of what would have been
    asked.

    For leaf regions (no descendants), the exclusive query is set equal
    to the inclusive query, so downstream code can treat the two columns
    uniformly without special-casing.
    """
    validate_regions(regions)
    by_parent = defaultdict(list)
    for r in regions:
        if r[2]:
            by_parent[r[2]].append(r)
    queries = []
    for r in regions:
        name, cat, parent, *_ = r
        base = region_term(r)
        inclusive = f"({base}) AND {DATE_FILTER}"
        descs = descendants(name, by_parent)
        if descs:
            child_term = "(" + " OR ".join(region_term(c) for c in descs) + ")"
            exclusive = f"(({base}) NOT {child_term}) AND {DATE_FILTER}"
        else:
            exclusive = inclusive
        queries.append({
            "region": name, "category": cat, "parent": parent or "",
            "n_descendants": len(descs),
            "inclusive_query": inclusive,
            "exclusive_query": exclusive,
        })
    return queries


# ---------- NCBI fetch ----------

def esearch_count(term, email=None, api_key=None, retries=4, sleep_s=0.34):
    """Call NCBI ESearch and return ``(count, full_url)``.

    Uses ``rettype=count`` so the response is small (no PMID list).
    Sleeps ``sleep_s`` after each successful request — NCBI's default
    cap is ~3 req/s without an API key and ~10 req/s with one. The
    caller, not this function, is responsible for choosing a polite
    sleep value.

    Errors are retried with linear backoff. Anything that looks like a
    transient network or parsing failure is retried up to ``retries``
    times; a payload that lacks ``count`` is also treated as a retry.
    """
    params = {"db": "pubmed", "term": term, "retmode": "json",
              "rettype": "count", "tool": TOOL_NAME}
    if email: params["email"] = email
    if api_key: params["api_key"] = api_key
    url = NCBI_ESEARCH + "?" + urllib.parse.urlencode(params)
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            time.sleep(sleep_s)
            res = payload.get("esearchresult", {})
            if "count" not in res:
                raise RuntimeError(f"no count in payload: {payload}")
            return int(res["count"]), url
        except (urllib.error.URLError, urllib.error.HTTPError,
                json.JSONDecodeError, RuntimeError) as e:
            last_err = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"ESearch failed after {retries} retries for {term[:120]!r}: {last_err}")


def fetch_all(queries, email, api_key, sleep_s):
    """Run ESearch for every entry in ``queries`` and return rows.

    For regions with at least one descendant the function makes two
    calls (inclusive, then exclusive). For leaf regions it makes one
    call and reuses the result for both columns. Progress is streamed
    to stderr so a long run is observable.
    """
    rows = []
    today = dt.date.today().isoformat()
    n = len(queries)
    for i, q in enumerate(queries, 1):
        inc, inc_url = esearch_count(q["inclusive_query"], email, api_key, sleep_s=sleep_s)
        if q["n_descendants"] > 0:
            exc, exc_url = esearch_count(q["exclusive_query"], email, api_key, sleep_s=sleep_s)
        else:
            exc, exc_url = inc, inc_url
        row = dict(q)
        row.update({
            "inclusive_count": inc, "exclusive_count": exc,
            "inclusive_url": inc_url, "exclusive_url": exc_url,
            "query_date": today,
        })
        rows.append(row)
        sys.stderr.write(f"[{i:>3}/{n}] {q['region']:<35} inc={inc:>7,}  exc={exc:>7,}\n")
    return rows


# ---------- Normalization & output ----------

def normalize(rows, key):
    """Add ``log10_{key}_p1`` and ``norm_{key}`` columns in place.

    ``log10(x + 1)`` because counts can be zero, and PubMed counts span
    several orders of magnitude (~40 to ~70k); a linear scale would let
    a couple of giants flatten everything else. Min-max into [0, 1] is
    purely for downstream visualization or thresholding — it is sample
    relative, not an absolute prevalence score.
    """
    logs = [math.log10(r[key] + 1) for r in rows]
    mn, mx = min(logs), max(logs)
    span = (mx - mn) or 1.0
    for r, l in zip(rows, logs):
        r[f"log10_{key}_p1"] = round(l, 4)
        r[f"norm_{key}"] = round((l - mn) / span, 4)
    return rows


def write_csv(rows, path):
    fields = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def make_lollipop(rows, out_png):
    """Render a horizontal lollipop chart to ``out_png``.

    Each region gets one row. The horizontal stem runs from x=1 to the
    inclusive count on a symlog axis. A circle marks the inclusive
    count; a diamond marks the exclusive count when it differs (only
    umbrella regions). Categories are ordered as first encountered in
    ``rows``, with regions sorted ascending by inclusive count within
    each category and a faint colored band behind them.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cat_order = []
    seen = set()
    for r in rows:
        if r["category"] not in seen:
            cat_order.append(r["category"])
            seen.add(r["category"])

    by_cat = defaultdict(list)
    for r in rows:
        by_cat[r["category"]].append(r)
    for c in by_cat:
        by_cat[c].sort(key=lambda r: r["inclusive_count"])  # ascending; will invert y

    display = []
    for c in cat_order:
        display.extend(by_cat[c])

    cmap = plt.get_cmap("tab20")
    cat_colors = {c: cmap(i % 20) for i, c in enumerate(cat_order)}

    fig_h = max(8, len(display) * 0.24)
    fig, ax = plt.subplots(figsize=(13, fig_h))

    ys = list(range(len(display)))
    for y, r in zip(ys, display):
        col = cat_colors[r["category"]]
        inc = max(r["inclusive_count"], 1)
        exc = max(r["exclusive_count"], 1)
        ax.hlines(y, 1, inc, color=col, alpha=0.5, lw=1.6)
        ax.plot(inc, y, "o", color=col, ms=8, mec="black", mew=0.4,
                label="inclusive" if y == 0 else None)
        if r["n_descendants"] > 0 and r["exclusive_count"] != r["inclusive_count"]:
            ax.plot(exc, y, "D", color=col, ms=5.5, mec="black", mew=0.3,
                    alpha=0.95)

    ax.set_yticks(ys)
    ax.set_yticklabels([r["region"] for r in display], fontsize=8.5)
    ax.set_xscale("symlog")
    ax.set_xlabel("PubMed papers, 2015-2025 (title/abstract, log scale)",
                  fontsize=10)
    ax.set_title(
        "Brain region prevalence in the PubMed literature, 2015-2025\n"
        "Circle = inclusive count.   Diamond = exclusive of subregion mentions.",
        fontsize=12)

    # Category bands and labels on the right
    prev_cat, band_start = None, 0
    xmax_for_label = max(r["inclusive_count"] for r in display) * 2
    for i, r in enumerate(display + [None]):
        cur = r["category"] if r else None
        if cur != prev_cat and prev_cat is not None:
            ax.axhspan(band_start - 0.5, i - 0.5,
                       color=cat_colors[prev_cat], alpha=0.07, zorder=0)
            mid = (band_start + i - 1) / 2
            ax.text(xmax_for_label, mid, prev_cat, fontsize=9,
                    color=cat_colors[prev_cat], ha="left", va="center",
                    fontweight="bold")
            band_start = i
        prev_cat = cur

    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3, ls=":")
    ax.set_xlim(left=1)
    fig.subplots_adjust(left=0.22, right=0.82, top=0.97, bottom=0.04)
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", default=None,
                    help="contact email sent to NCBI (recommended)")
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--sleep", type=float, default=0.34,
                    help="seconds between calls; 0.34 ~3 req/s without API key")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--dry-run", action="store_true",
                    help="print queries; do not call NCBI")
    args = ap.parse_args()

    queries = build_queries(REGIONS)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    with open(out / "queries.json", "w", encoding="utf-8") as f:
        json.dump(queries, f, indent=2)

    if args.dry_run:
        for q in queries:
            print(f"\n=== {q['region']} ({q['category']}) ===")
            print("INC:", q["inclusive_query"])
            if q["n_descendants"]:
                print("EXC:", q["exclusive_query"])
        return

    rows = fetch_all(queries, args.email, args.api_key, args.sleep)
    rows = normalize(rows, "inclusive_count")
    rows = normalize(rows, "exclusive_count")

    write_csv(rows, out / "regions.csv")
    with open(out / "regions.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)
    make_lollipop(rows, out / "prevalence.png")

    print(f"\nWrote {out/'regions.csv'}")
    print(f"Wrote {out/'regions.json'}")
    print(f"Wrote {out/'prevalence.png'}")

    top = sorted(rows, key=lambda r: r["inclusive_count"], reverse=True)[:15]
    print("\nTop 15 by inclusive count (2015-2025):")
    for r in top:
        print(f"  {r['region']:<35} {r['inclusive_count']:>8,}  "
              f"(exclusive {r['exclusive_count']:,})")

if __name__ == "__main__":
    main()
