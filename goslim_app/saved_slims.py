"""Read/write/list saved working-slim lists as human-readable JSON files."""

import json
import re
from pathlib import Path


def slugify(name):
    """Turn a display name into a filesystem-safe slug."""
    slug = re.sub(r"[^a-z0-9]+", "-", name.strip().lower()).strip("-")
    return slug or "untitled"


def _iter_saved(saved_dir):
    """Yield (mtime, stem, payload) for every readable saved slim file."""
    saved_dir = Path(saved_dir)
    if not saved_dir.exists():
        return
    for path in saved_dir.glob("*.json"):
        try:
            payload = read_slim(path)
        except (json.JSONDecodeError, OSError):
            continue
        yield path.stat().st_mtime, path.stem, payload


def list_saved(saved_dir):
    """Return {slug: display_label} for every saved slim, newest first."""
    entries = sorted(_iter_saved(saved_dir), reverse=True)
    labels = {}
    for _mtime, stem, payload in entries:
        stats = payload.get("stats", {})
        labels[stem] = f"{payload.get('name', stem)} ({stats.get('n_terms', '?')} terms)"
    return labels


def list_saved_for_aspect(saved_dir, aspect, key_prefix=""):
    """Return {key_prefix+slug: display_label} for saved slims matching one
    GO aspect - used to offer saved slims as full-fledged 'base slim' choices.
    key_prefix disambiguates these keys from the static base-slim choices."""
    entries = sorted(_iter_saved(saved_dir), reverse=True)
    labels = {}
    for _mtime, stem, payload in entries:
        if payload.get("aspect") != aspect:
            continue
        stats = payload.get("stats", {})
        labels[f"{key_prefix}{stem}"] = f"Saved: {payload.get('name', stem)} ({stats.get('n_terms', '?')} terms)"
    return labels


def read_slim(path):
    """Load one saved slim JSON file."""
    with open(path) as f:
        return json.load(f)


def write_slim(saved_dir, payload):
    """Write a new saved slim, appending -2/-3/... on filename collision.
    Never silently overwrites an existing save."""
    saved_dir = Path(saved_dir)
    saved_dir.mkdir(parents=True, exist_ok=True)
    base_slug = slugify(payload["name"])

    slug = base_slug
    n = 2
    while (saved_dir / f"{slug}.json").exists():
        slug = f"{base_slug}-{n}"
        n += 1

    path = saved_dir / f"{slug}.json"
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    return path
