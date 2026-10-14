"""
check_refs.py — every DOI in a .bib must resolve, and to the paper the entry describes.

    python check_refs.py paper/refs.bib

For each entry with a doi field, Crossref is queried and the registered metadata is
compared with the entry: year, volume, first page (or article number) and the first
author's family name. Any unresolved DOI or mismatch is listed; exit status 1 if any.
Entries without a doi are listed as such (not a failure).
"""

import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.crossref.org/works/"
UA = "pulse-ai-check-refs/1.0"
_LATEX = re.compile(r"\\[a-zA-Z]+\s*|\\.|[{}]")


def entries(text):
    """(key, {field: value}) for each @type{key, ...} entry."""
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)\n\}", text, re.S):
        fields = {k.lower(): v.strip()
                  for k, v in re.findall(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", m[3])}
        yield m[2], fields


def norm(s):
    """Plain lowercase ASCII letters/digits: strips LaTeX, accents, punctuation."""
    s = _LATEX.sub("", s or "")
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]", "", s.lower())


def crossref(doi):
    url = API + urllib.parse.quote(doi, safe="/")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)["message"]
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == 2:
                raise
        except urllib.error.URLError:
            if attempt == 2:
                raise
        time.sleep(2 * (attempt + 1))


def check(key, f):
    """List of problems for one entry ([] = OK)."""
    m = crossref(f["doi"])
    if m is None:
        return ["DOI does not resolve (Crossref 404)"]
    problems = []
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    if "year" in f and year and int(f["year"]) != year:
        problems.append(f"year {f['year']} != Crossref {year}")
    if "volume" in f and m.get("volume") and norm(f["volume"]) != norm(m["volume"]):
        problems.append(f"volume {f['volume']} != Crossref {m['volume']}")
    if "pages" in f:
        first = re.split(r"-+|–", f["pages"])[0]
        ref = (m.get("page") or m.get("article-number") or "")
        ref_first = re.split(r"-", ref)[0]
        if ref_first and norm(first) != norm(ref_first):
            problems.append(f"first page {first} != Crossref {ref_first}")
    authors = m.get("author") or []
    if "author" in f and authors:
        a0 = authors[0]
        if f["author"].startswith("{"):
            # corporate author: Crossref often splits the name across given/family
            bib_first = f["author"]
            cr_first = a0.get("name") or f"{a0.get('given', '')} {a0.get('family', '')}"
            ok = norm(f"{a0.get('family', '')}") in norm(bib_first) and \
                norm(cr_first)[:20] == norm(bib_first)[:20]
        else:
            bib_first = f["author"].split(" and ")[0].split(",")[0]
            cr_first = a0.get("family") or a0.get("name", "")
            ok = norm(bib_first) == norm(cr_first)
        if not ok:
            problems.append(f"first author {bib_first!r} != Crossref {cr_first!r}")
    return problems


def main(path):
    text = open(path, encoding="utf-8").read()
    failed = 0
    for key, f in entries(text):
        if "doi" not in f:
            print(f"  --    {key}: no doi")
            continue
        try:
            problems = check(key, f)
        except Exception as exc:          # noqa: BLE001 - network failure is a failure
            problems = [f"lookup error: {type(exc).__name__}: {exc}"]
        print(f"  {'OK  ' if not problems else 'FAIL'}  {key}  {f['doi']}")
        for p in problems:
            print(f"          {p}")
        failed += bool(problems)
    print(f"\n{failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "paper/refs.bib"))
