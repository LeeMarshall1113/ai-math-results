"""Novelty queries for 10-face quadrilateral torus, run 2026-10-29. Prints query strings and hits."""
import json, urllib.request, urllib.parse, re, time
def get(u):
    req = urllib.request.Request(u, headers={"User-Agent": "novelty-check"})
    return urllib.request.urlopen(req, timeout=60).read().decode()
for q in ['Seamons torus quadrilateral', '"quadrilateral faces" torus', 'toroidal polyhedron quadrilateral', 'creators.name:Seamons',
          'Szilassi', '"polyhedral torus"', '"{4,4}" torus realization', 'torus "ten" quadrilateral faces', '"equivelar" torus']:
    u = "https://zenodo.org/api/records?" + urllib.parse.urlencode({"q": q, "size": 25})
    print("== Zenodo", u)
    try:
        d = json.loads(get(u))
        for h in d["hits"]["hits"]:
            print("  ", h["metadata"].get("publication_date"), "|", h["metadata"]["title"][:110], "|", h.get("doi"))
    except Exception as e:
        print("  ERR", e)
    time.sleep(2)
for q in ['abs:toroidal AND abs:quadrilateral AND abs:faces', 'abs:Szilassi', 'abs:"polyhedral torus"', 'abs:equivelar AND abs:torus', 'abs:"{4,4}" AND abs:polyhedral', 'abs:"polyhedral realization" AND abs:torus']:
    u = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": q, "max_results": 100, "sortBy": "submittedDate"})
    print("== arXiv", u)
    try:
        x = get(u)
        for e in re.findall(r"<entry>(.*?)</entry>", x, re.S):
            t = re.search(r"<title>(.*?)</title>", e, re.S).group(1).split()
            p = re.search(r"<published>(.*?)</published>", e).group(1)[:10]
            print("  ", p, "|", " ".join(t)[:120])
    except Exception as e:
        print("  ERR", e)
    time.sleep(3)
for q in ["toroidal polyhedron","polyhedral torus","Szilassi"]:
    u = "https://api.osf.io/v2/preprints/?" + urllib.parse.urlencode({"filter[title]": q})
    print("== OSF preprints", u)
    try:
        d = json.loads(get(u)); [print("  ", r["attributes"]["title"]) for r in d["data"]]
    except Exception as e: print("  ERR", e)
    u = "https://api.osf.io/v2/nodes/?" + urllib.parse.urlencode({"filter[title]": q})
    print("== OSF nodes", u)
    try:
        d = json.loads(get(u)); [print("  ", r["attributes"]["title"]) for r in d["data"]]
    except Exception as e: print("  ERR", e)
