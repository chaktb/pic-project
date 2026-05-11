#!/usr/bin/env python3
"""Auto-build hook for pic-project.

Triggered by Claude Code Stop hook. Steps:
  1. Re-run stale n_*.py plotting scripts (those containing 'savefig').
  2. Sync any n_*.png from repo root into public/.
  3. Auto-create public/<topic>/index.html for new n_*.png topics.
  4. Regenerate landing-page cards between AUTO_CARDS_START / END markers.
  5. git add + commit + push (main branch only).

Exits 0 on success or partial-failure so the hook never blocks Claude.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
PUBLIC = REPO / "public"
LOG = []


def log(msg: str) -> None:
    LOG.append(f"[auto-build] {msg}")


def find_python() -> str:
    cand = Path.home() / ".nix-profile/bin/python"
    if cand.exists():
        return str(cand)
    for name in ("python3", "python"):
        p = shutil.which(name)
        if p:
            return p
    return "python3"


def newer(a: Path, b: Path) -> bool:
    if not b.exists():
        return True
    return a.stat().st_mtime > b.stat().st_mtime


def run_plotting_scripts(py_bin: str) -> None:
    for py in sorted(REPO.glob("n_*.py")):
        try:
            text = py.read_text()
        except Exception:
            continue
        if "savefig" not in text:
            continue
        png = py.with_suffix(".png")
        if not newer(py, png):
            continue
        log(f"running {py.name}")
        env = {**os.environ, "MPLBACKEND": "Agg"}
        r = subprocess.run(
            [py_bin, str(py)],
            capture_output=True,
            env=env,
            cwd=str(REPO),
            timeout=120,
        )
        if r.returncode != 0:
            err = r.stderr.decode(errors="replace")[:200]
            log(f"  failed: {err}")


def sync_pngs() -> None:
    PUBLIC.mkdir(exist_ok=True)
    for png in REPO.glob("n_*.png"):
        target = PUBLIC / png.name
        if newer(png, target):
            target.write_bytes(png.read_bytes())
            log(f"synced {png.name} -> public/")


PAGE_TEMPLATE = """<!-- AUTO_GENERATED_PAGE -->
<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>{title}</title>
<style>
:root {{ color-scheme: light dark; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; max-width: 880px; margin: 2rem auto; padding: 0 1rem; line-height: 1.6; }}
.back {{ display: inline-block; color: #666; text-decoration: none; margin-bottom: 1rem; }}
.back:hover {{ text-decoration: underline; }}
img {{ width: 100%; height: auto; border: 1px solid #ddd; border-radius: 8px; }}
@media (prefers-color-scheme: dark) {{ img {{ border-color: #333; }} }}
</style>
</head>
<body>
<a class="back" href="../">← Home</a>
<h1>{title}</h1>
<img src="../{png}" alt="{topic}" />
</body>
</html>
"""


def autogen_pages() -> None:
    for png in REPO.glob("n_*.png"):
        topic = png.stem[2:]  # strip "n_"
        page_dir = PUBLIC / topic
        page = page_dir / "index.html"
        if page.exists():
            continue
        page_dir.mkdir(parents=True, exist_ok=True)
        title = topic.replace("-", " ").replace("_", " ").title()
        page.write_text(PAGE_TEMPLATE.format(topic=topic, png=png.name, title=title))
        log(f"created {page.relative_to(REPO)}")


def update_landing_cards() -> None:
    landing = PUBLIC / "index.html"
    if not landing.exists():
        return
    content = landing.read_text()
    if "AUTO_CARDS_START" not in content or "AUTO_CARDS_END" not in content:
        return  # markers absent, nothing to do

    cards = []
    for d in sorted(PUBLIC.iterdir()):
        if not d.is_dir():
            continue
        idx = d / "index.html"
        if not idx.exists():
            continue
        try:
            html = idx.read_text()
        except Exception:
            continue
        m = re.search(r"<title>([^<]+)</title>", html)
        title = m.group(1).strip() if m else d.name.replace("-", " ").title()
        label = "Topic"
        # First sentence of first <p> (if any) becomes the card description
        pm = re.search(r"<p[^>]*>([^<]+)</p>", html)
        desc = pm.group(1).strip() if pm else f"See {d.name}."
        if len(desc) > 140:
            desc = desc[:137] + "..."
        # Big heading: uppercase short label from title
        heading = title.upper()
        cards.append(
            f'      <a class="card" href="./{d.name}/">'
            f'<div class="label">{label}</div>'
            f'<h2>{heading}</h2>'
            f'<p>{desc}</p></a>'
        )

    new_block = (
        "<!-- AUTO_CARDS_START -->\n"
        + "\n".join(cards)
        + "\n      <!-- AUTO_CARDS_END -->"
    )
    new_content = re.sub(
        r"<!-- AUTO_CARDS_START -->.*?<!-- AUTO_CARDS_END -->",
        new_block,
        content,
        flags=re.DOTALL,
    )
    if new_content != content:
        landing.write_text(new_content)
        log("updated landing cards")


def commit_and_push() -> None:
    branch = subprocess.run(
        ["git", "branch", "--show-current"],
        capture_output=True, text=True, cwd=str(REPO),
    ).stdout.strip()
    if branch != "main":
        log(f"branch is '{branch}', skipping commit/push")
        return

    status = subprocess.run(
        ["git", "status", "--porcelain"],
        capture_output=True, text=True, cwd=str(REPO),
    ).stdout.strip()
    if not status:
        return

    # Stage tracked changes + new files in safe paths only
    subprocess.run(["git", "add", "-A"], check=True, cwd=str(REPO))
    msg = f"auto-build: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}"
    r = subprocess.run(
        ["git", "commit", "-m", msg],
        capture_output=True, text=True, cwd=str(REPO),
    )
    if r.returncode != 0:
        # Possibly nothing actually staged (e.g. all changes ignored). Bail quietly.
        return
    log(f"committed: {msg}")
    p = subprocess.run(
        ["git", "push"], capture_output=True, text=True, cwd=str(REPO),
    )
    if p.returncode != 0:
        log(f"push failed: {p.stderr.strip()[:200]}")
    else:
        log("pushed to origin/main")


def main() -> int:
    # Drain stdin (hook payload) but we don't need it
    try:
        sys.stdin.read()
    except Exception:
        pass

    py_bin = find_python()
    try:
        run_plotting_scripts(py_bin)
        sync_pngs()
        autogen_pages()
        update_landing_cards()
        commit_and_push()
    except Exception as e:
        log(f"FATAL: {e!r}")

    if LOG:
        msg = "\n".join(LOG)
        # Print as systemMessage so it appears in the UI compactly
        print(json.dumps({"systemMessage": msg, "suppressOutput": True}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
