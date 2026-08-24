"""Ancestor-closure computation and fast slim mapping (replaces goatools.mapslim).

mapslim() re-walks the DAG on every call and can't be memoized across a
changing slim set. Instead we compute, once, the full set of ancestors (self
included) for every GO term ever observed in the GAF, then mapping a term to
a slim becomes plain set algebra - fast enough to redo on every UI click.
"""

_MAP_MEMO = {}


def build_ancestor_closure(godag, seed_ids, use_part_of=True):
    """Return {go_id: frozenset(go_id and all its is_a/part_of ancestors)}
    for every id in seed_ids. Iterative post-order DFS with memoization so
    shared ancestors are only walked once."""
    closure = {}

    def parents_of(term):
        # is_a parents, plus part_of relationship targets if enabled
        parents = list(term.parents)
        if use_part_of:
            rel = getattr(term, "relationship", None) or {}
            parents += list(rel.get("part_of", []))
        return parents

    def visit(go_id):
        if go_id in closure:
            return closure[go_id]
        term = godag.get(go_id)
        if term is None:
            closure[go_id] = frozenset({go_id})
            return closure[go_id]
        acc = {go_id}
        # guard against cycles (shouldn't exist in a real DAG, but be safe)
        closure[go_id] = frozenset(acc)
        for parent in parents_of(term):
            acc |= visit(parent.id)
        closure[go_id] = frozenset(acc)
        return closure[go_id]

    for seed in seed_ids:
        visit(seed)
    return closure


def build_depth_map(godag, seed_ids):
    """Return {go_id: depth} using goatools' own depth attribute, falling
    back to 0 for ids missing from the release-matched DAG."""
    depth = {}
    for go_id in seed_ids:
        term = godag.get(go_id)
        depth[go_id] = term.depth if term is not None else 0
    return depth


def map_term_to_slim(go_id, slim_ids, anc):
    """Direct (minimal) slim ancestors of go_id: intersect its ancestor
    closure with the slim set, then drop any candidate that is itself an
    ancestor of another candidate (keep only the most specific hits)."""
    key = (go_id, slim_ids)
    cached = _MAP_MEMO.get(key)
    if cached is not None:
        return cached

    hits = anc.get(go_id, frozenset()) & slim_ids
    if len(hits) > 1:
        hits = {c for c in hits if not any(c != o and c in anc.get(o, ()) for o in hits)}

    _MAP_MEMO[key] = hits
    return hits
