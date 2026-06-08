#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cloudflare "Manifest path must be URI encoded" fixer.

Run this from the ROOT of your repo (the folder that contains `public/`).

What it does:
  1. DIAGNOSE  : lists every file under public/ whose path contains
                 characters Cloudflare's manifest validator rejects
                 (non-ASCII, spaces, and a set of reserved/special chars),
                 and flags macOS NFD (decomposed Hangul) filenames.
  2. FIX       : two independent fixes you can choose:
       --normalize : rename NFD -> NFC (keeps Korean names, just fixes the
                     broken byte form). Often enough on its own.
       --asciify   : transliterate/replace non-ASCII names with safe ASCII
                     and rewrite href/src references inside html/css/js.

Usage:
    python3 cf_fix.py                 # diagnose only (no changes)
    python3 cf_fix.py --normalize     # apply NFC normalization
    python3 cf_fix.py --asciify       # rename to ASCII + fix links
    python3 cf_fix.py --normalize --asciify
"""

import os
import re
import sys
import unicodedata
import argparse

PUBLIC = "public"
TEXT_EXT = {".html", ".htm", ".css", ".js", ".json", ".svg", ".xml", ".txt"}

# Characters that are safe unencoded in a path segment.
SAFE = re.compile(r"^[A-Za-z0-9._~/\-]+$")


def is_problematic(rel_path: str) -> bool:
    """True if the path would need URI encoding (Cloudflare rejects it)."""
    return not SAFE.match(rel_path)


def is_nfd(name: str) -> bool:
    """True if the filename is not in NFC (e.g. macOS decomposed Hangul)."""
    return unicodedata.normalize("NFC", name) != name


def walk_files(root):
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            full = os.path.join(dirpath, f)
            rel = os.path.relpath(full, root).replace(os.sep, "/")
            yield full, rel


# ---------------------------------------------------------------------------
def diagnose():
    problems = []
    for full, rel in walk_files(PUBLIC):
        bad = is_problematic(rel)
        nfd = is_nfd(rel)
        if bad or nfd:
            tags = []
            if nfd:
                tags.append("NFD")
            if bad:
                # what kind of bad?
                if re.search(r"[^\x00-\x7F]", rel):
                    tags.append("non-ASCII")
                if " " in rel:
                    tags.append("space")
                specials = set(re.findall(r"[#%?\[\]{}&+'\"()@!$,;=]", rel))
                if specials:
                    tags.append("special:" + "".join(sorted(specials)))
            problems.append((rel, tags))

    print(f"Scanned files under '{PUBLIC}/'.")
    print(f"Problematic paths: {len(problems)}\n")
    for rel, tags in problems:
        print(f"  [{', '.join(tags)}]  {rel}")
    if not problems:
        print("  (none found — the issue may be elsewhere)")
    print()
    return [p[0] for p in problems]


# ---------------------------------------------------------------------------
def normalize_nfc():
    """Rename any NFD filename/dir to its NFC form (Korean names preserved)."""
    renamed = 0
    # rename deepest first so parent dir renames don't break child paths
    for dirpath, dirs, files in os.walk(PUBLIC, topdown=False):
        for name in files + dirs:
            nfc = unicodedata.normalize("NFC", name)
            if nfc != name:
                src = os.path.join(dirpath, name)
                dst = os.path.join(dirpath, nfc)
                os.rename(src, dst)
                renamed += 1
                print(f"  NFC: {name}  ->  {nfc}")
    print(f"\nNormalized {renamed} name(s) to NFC.")
    return renamed


# ---------------------------------------------------------------------------
def slugify(name: str) -> str:
    """Make an ASCII-safe filename, preserving the extension."""
    base, ext = os.path.splitext(name)
    # transliterate what we can, drop the rest
    base_nfkd = unicodedata.normalize("NFKD", base)
    ascii_base = base_nfkd.encode("ascii", "ignore").decode("ascii")
    ascii_base = re.sub(r"[^A-Za-z0-9._\-]+", "_", ascii_base).strip("_")
    if not ascii_base:
        # pure-Hangul name with nothing transliterable -> stable hash tag
        ascii_base = "file_" + format(abs(hash(name)) % 10_000_000, "07d")
    return ascii_base + ext


def asciify():
    """Rename non-ASCII files to ASCII and rewrite references in text files."""
    mapping = {}  # old_rel -> new_rel
    for dirpath, dirs, files in os.walk(PUBLIC, topdown=False):
        for name in files:
            if is_problematic(name) or is_nfd(name):
                new = slugify(name)
                src_rel = os.path.relpath(os.path.join(dirpath, name), PUBLIC)
                # avoid collisions
                target = os.path.join(dirpath, new)
                k = 1
                while os.path.exists(target) and os.path.basename(target) != name:
                    b, e = os.path.splitext(new)
                    target = os.path.join(dirpath, f"{b}_{k}{e}")
                    k += 1
                os.rename(os.path.join(dirpath, name), target)
                dst_rel = os.path.relpath(target, PUBLIC)
                mapping[name] = os.path.basename(target)
                print(f"  ASCII: {src_rel}  ->  {dst_rel}")

    if not mapping:
        print("  (no non-ASCII files to rename)")
        return mapping

    # rewrite references (match by basename, safest for varied link forms)
    print("\nUpdating references in html/css/js ...")
    changed = 0
    for full, rel in walk_files(PUBLIC):
        if os.path.splitext(full)[1].lower() not in TEXT_EXT:
            continue
        with open(full, "r", encoding="utf-8", errors="ignore") as fh:
            txt = fh.read()
        new_txt = txt
        for old, new in mapping.items():
            if old in new_txt:
                new_txt = new_txt.replace(old, new)
        if new_txt != txt:
            with open(full, "w", encoding="utf-8") as fh:
                fh.write(new_txt)
            changed += 1
            print(f"  updated refs in {rel}")
    print(f"\nRewrote references in {changed} file(s).")
    return mapping


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--normalize", action="store_true",
                    help="rename NFD filenames to NFC (keeps Korean)")
    ap.add_argument("--asciify", action="store_true",
                    help="rename non-ASCII files to ASCII and fix links")
    args = ap.parse_args()

    if not os.path.isdir(PUBLIC):
        print(f"ERROR: '{PUBLIC}/' not found. Run this from the repo root.")
        sys.exit(1)

    print("=== DIAGNOSE ===")
    diagnose()

    if args.normalize:
        print("=== NORMALIZE (NFD -> NFC) ===")
        normalize_nfc()
        print()
    if args.asciify:
        print("=== ASCIIFY ===")
        asciify()
        print()

    if args.normalize or args.asciify:
        print("Re-checking after fixes:")
        remaining = diagnose()
        if not remaining:
            print("All clear. Now commit & push:")
            print('  git add -A && git commit -m "Fix asset filenames for Cloudflare" && git push')


if __name__ == "__main__":
    main()