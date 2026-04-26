#!/usr/bin/env python3
import argparse
import errno
import os
import shutil
import sys
import secrets

from pathlib import Path
from typing import List

MAX_LEN = 51  # path length limit

ENV_VARS = ["TMPDIR"]

def real_len(p):
    try:
        return len(p.resolve(strict=True).as_posix())
    except FileNotFoundError:
        return len(p.as_posix())

def ensure_dir(p):
    p.mkdir(parents=True, exist_ok=True)

def same_filesystem(a, b):
    try:
        return os.stat(a).st_dev == os.stat(b).st_dev
    except FileNotFoundError:
        return False

def candidate_roots():
    roots: List[Path] = []
    for ev in ENV_VARS:
        v = os.environ.get(ev)
        if v:
            roots.append(Path(v))
    roots.append(Path.cwd())
    return [r for r in roots if r.exists() and os.access(r, os.W_OK)]

def ancestor_roots(src):
    roots = []
    seen = set()
    cur = src.resolve().parent  # start at parent
    while True:
        if cur is None:
            break
        sp = cur.as_posix()
        if sp not in seen and cur.exists() and os.access(cur, os.W_OK):
            roots.append(cur)
            seen.add(sp)
        if cur.parent == cur:  # reached root
            break
        cur = cur.parent
    return roots


def make_very_short(root, upper = 5):
    root_abs = root.resolve()
    root_len = len(root_abs.as_posix())
    prefixes = ("i.", "x.", ".i.", ".x.")

    for pref in prefixes:
        # remaining length for the random suffix
        budget = MAX_LEN - root_len - 1 - len(pref)
        if budget < 1:
            continue
        n = min(upper, budget)

        suffix = secrets.token_hex((n + 1) // 2)[:n]
        cand = root_abs / (pref + suffix)

        if len(cand.as_posix()) <= MAX_LEN and not cand.exists():
            return cand

    return None


def _mapped_symlink_target(src_link, src_root, dst_root, dst_path):
    raw = os.readlink(src_link)
    target_abs = Path(raw) if os.path.isabs(raw) else (src_link.parent / raw).resolve(strict=False)
    try:
        rel = target_abs.relative_to(src_root)
        mapped_abs = dst_root / rel
        return os.path.relpath(mapped_abs, start=dst_path.parent)
    except ValueError:
        return raw

def _replace_with_symlink(link_target, link_path):
    try:
        link_path.unlink()
    except FileNotFoundError:
        pass
    except IsADirectoryError:
        try:
            link_path.rmdir()
        except OSError as e:
            if e.errno != errno.ENOTEMPTY:
                raise
            raise RuntimeError(f"Refusing to replace non-empty directory: {link_path}") from e
    os.symlink(link_target, link_path)

def hardcopy(src, dst, fallback_copy=True):
    src = Path(src).resolve()
    dst = Path(dst).resolve()
    dst.mkdir(parents=True, exist_ok=True)
    copy_type = "hardlink"
    for root, dirnames, filenames in os.walk(src, topdown=True, followlinks=False):
        root = Path(root)
        rel = root.relative_to(src)
        dest_root = dst / rel
        dest_root.mkdir(parents=True, exist_ok=True)

        # Create entries for directories at this level
        for name in dirnames:
            s = root / name
            d = dest_root / name
            if s.is_symlink():
                target = _mapped_symlink_target(s, src, dst, d)
                _replace_with_symlink(target, d)
            else:
                d.mkdir(exist_ok=True)

        # Files (regular + symlinks-to-files)
        for name in filenames:
            s = root / name
            d = dest_root / name
            if s.is_symlink():
                target = _mapped_symlink_target(s, src, dst, d)
                _replace_with_symlink(target, d)
            else:
                try:
                    try:
                        d.unlink()
                    except FileNotFoundError:
                        pass
                    os.link(s, d)
                except OSError:
                    if fallback_copy:
                        shutil.copy2(s, d)
                        copy_type = "copy"
    return copy_type

def link_or_copy_tree(src, dst):
    """
    Attempt to hardlink src/* into dst/.
    Returns "hardlink" if linked, otherwise "copy".
    """
    # If not same FS, we’ll copy everything
    if not same_filesystem(src, dst.resolve().parent):
        shutil.copytree(src, dst, dirs_exist_ok=True)
        return "copy"
    else:
        return hardcopy(src, dst)


def pick_target(src):
    """
    Decide target and initial mode.
    - If src path already short: return (src, "original")
    - Else: try ancestors first (same FS), then env/common roots.
    For non-original, we’ll materialize and then report "hardlink" or "copy".
    """

    src_real = src.resolve()
    if real_len(src_real) < MAX_LEN:
        return src_real, "original"

    # 1) Try ancestors (same FS most likely)
    for r in ancestor_roots(src_real):
        d = make_very_short(r)
        if d:
            return d, "materialize"

    # 2) Try env/common roots
    for r in candidate_roots():
        d = make_very_short(r)
        if d:
            return d, "materialize"

    raise RuntimeError("No short writable directory (<57 chars) could be created")

def main():
    ap = argparse.ArgumentParser(
        description="Ensure IEDB path is short; hardlink when possible, else copy."
    )
    ap.add_argument("--src", required=True, help="Path to existing IEDB directory")
    args = ap.parse_args()
    src = Path(args.src)
    if not src.is_dir():
        print(f"smart_link_iedb.py: source not a directory: {src}", file=sys.stderr)
        sys.exit(1)

    target, mode = pick_target(src)

    if mode == "original":
        print(target.as_posix())
        print("original")
        return

    # Materialize (link or copy)
    try:
        final_mode = link_or_copy_tree(src, target)
    except Exception as e:
        print(
            "smart_link_iedb.py: materialization failed: "
            f"{e}. Provide a writable short path or bind-mount a short path.",
            file=sys.stderr,
        )
        sys.exit(2)

    print(target.as_posix())
    print(final_mode)  # "hardlink" or "copy"

if __name__ == "__main__":
    main()
