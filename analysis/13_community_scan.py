"""Track M: has the hair-loss community already been talking about this? Searches YouTube for pirfenidone, Kitoscell and
fibrosis in hair loss, and records title, channel, date, views and description for each hit. Frozen in
data/community/youtube_scan.json on the first run (search results change daily), so reruns are exact.
Reddit's API and the Cloudflare-protected forums (HairLossTalk, Ray Peat Forum) blocked scripted access on Oct 4; what was
found there by web search is written up in RESEARCH_LOG.md, Track M."""
import json
import re
import time
import urllib.parse
import urllib.request

import pandas as pd

from common import RESULTS, ROOT, save_summary

FROZEN = ROOT / "data" / "community" / "youtube_scan.json"
QUERIES = ["pirfenidone hair loss", "pirfenidone hair", "kitoscell hair", "topical pirfenidone scalp", "perifollicular fibrosis hair loss",
           "fibrosis hair loss", "CXCL12 hair loss", "more plates more dates pirfenidone", "more plates more dates fibrosis"]
UA = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/129.0 Safari/537.36",
      "Accept-Language": "en-US,en;q=0.9"}
TOPIC = re.compile(r"(?i)pirfenidone|kitoscell|pirfesca|fibros|cxcl12")


def page(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60).read().decode("utf-8", "ignore")


def search(q):
    m = re.search(r"var ytInitialData = (\{.*?\});</script>", page("https://www.youtube.com/results?search_query=" + urllib.parse.quote_plus(q)))
    out = []

    def walk(o):
        if isinstance(o, dict):
            if "videoRenderer" in o:
                v = o["videoRenderer"]
                out.append({"video_id": v.get("videoId"), "title": "".join(r.get("text", "") for r in v.get("title", {}).get("runs", [])),
                            "channel": "".join(r.get("text", "") for r in v.get("ownerText", {}).get("runs", []))})
            for x in o.values():
                walk(x)
        elif isinstance(o, list):
            for x in o:
                walk(x)
    if m:
        walk(json.loads(m.group(1)))
    return out[:20]


def details(vid):
    m = re.search(r"var ytInitialPlayerResponse = (\{.*?\});(?:var|</script>)", page("https://www.youtube.com/watch?v=" + vid))
    if not m:
        return {}
    d = json.loads(m.group(1))
    vd, mf = d.get("videoDetails", {}), d.get("microformat", {}).get("playerMicroformatRenderer", {})
    return {"published": (mf.get("publishDate") or "")[:10], "views": int(vd.get("viewCount") or 0), "minutes": round(int(vd.get("lengthSeconds") or 0) / 60, 1),
            "description": (vd.get("shortDescription") or "")[:3000]}


if FROZEN.exists():
    scan = json.loads(FROZEN.read_text())
else:
    hits = {}
    for q in QUERIES:
        for v in search(q):
            if TOPIC.search(v["title"]):
                hits.setdefault(v["video_id"], {**v, "queries": []})["queries"].append(q)
        time.sleep(1)
    for vid, v in hits.items():
        v.update(details(vid))
        time.sleep(1)
    scan = {"scanned_on": time.strftime("%Y-%m-%d"), "queries": QUERIES, "videos": list(hits.values())}
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps(scan, indent=1))

df = pd.DataFrame(scan["videos"])
df["mentions_pirfenidone"] = (df.title + " " + df.description.fillna("")).str.contains("(?i)pirfenidone|kitoscell|pirfesca")
df = df.sort_values("published")
df.drop(columns=["description"]).to_csv(RESULTS / "13_community_videos.csv", index=False)
pd.set_option("display.width", 220)
print(f"YouTube scan from {scan['scanned_on']}: {len(df)} on-topic videos")
print(df[["published", "channel", "views", "mentions_pirfenidone", "title"]].to_string(index=False))
first = df[df.mentions_pirfenidone].head(1)
save_summary("13", {"scanned_on": scan["scanned_on"], "videos": int(len(df)), "videos_on_pirfenidone": int(df.mentions_pirfenidone.sum()),
                    "earliest_pirfenidone_video": first[["published", "channel", "title"]].to_dict(orient="records"),
                    "channels": df.groupby("channel").size().sort_values(ascending=False).to_dict()})
