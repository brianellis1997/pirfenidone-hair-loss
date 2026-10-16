# An Audited LLM Evidence Agent for Drug Repurposing: Testing Pirfenidone for Pattern Hair Loss

Code for my High-risk Project in A I 380H, AI in Healthcare (UT Austin, Fall 2026). Brian Ellis.

Pirfenidone is an approved drug for lung scarring (idiopathic pulmonary fibrosis). Male pattern hair loss has a
scarring side next to the hormonal one, and nobody has published a test of whether pirfenidone reaches it. This project
checks the idea against public gene expression, genetics, single-cell, FDA adverse event and literature data, and builds
an LLM evidence agent whose grades are checked against blind reference labels.

The report, slides and recorded presentation are turned in on Canvas. This repo is the code.

## Start here

Open [`HRP_Code_Walkthrough.ipynb`](HRP_Code_Walkthrough.ipynb). It's already executed, so GitHub shows every output.
It goes track by track: the question, the core code pulled straight from each script, the saved results and the figure.
Three numbers are recomputed live from the intermediate tables (the FDA comparator ratio, the evidence agent's grades and
the dose arithmetic) to show the saved results come from the code.

## What's where

| Path | What's there |
|---|---|
| `analysis/00_download_data.sh` | Downloads every public dataset into `data/` (about 800 MB, mostly the baldness GWAS) |
| `analysis/01_*.py` to `15_*.py` | One script per research track, run in order. Each writes `results/NN_*` and `results/NN_summary.json` |
| `analysis/run_all.sh` | Rebuilds everything from scratch |
| `src/` | Dataset loaders, differential expression and drug scoring |
| `results/` | Every table and summary, with run logs in `results/logs/` |
| `figs/report/` | The report figures (PDF and PNG) and their captions |
| `data/llm_agent/` | The frozen PubMed set and every LLM reply, so the evidence agent replays exactly |
| `data/community/` | The frozen YouTube scan for the community track |
| `DATA_SOURCES.md` | Every dataset and API, with accession numbers |
| `build_walkthrough.py` | Generates the walkthrough notebook |

## Rerunning it

```bash
pip install -r requirements.txt   # Python 3.11
zsh analysis/run_all.sh           # downloads the data, then runs tracks 01 to 15 (about 10 minutes after downloads)
```

The evidence agent replays every LLM reply from `data/llm_agent/`, so a rebuild needs no API key. Only a fresh run with
new prompts or new abstracts calls OpenRouter, and that needs `OPENROUTER_API_KEY` set. PubMed counts and openFDA reports
keep growing, so a rebuild months from now can drift a little from the saved numbers in tracks 06 and 07.

All data is public. Sources and accession numbers are in [`DATA_SOURCES.md`](DATA_SOURCES.md).

## AI disclosure

I wrote and tested the code together with Claude (Opus 5.5), and it also helped me come up with ideas and draft the
notebook and this README. The question and the calls on where to take it were mine, and I went through the results
myself before anything went in. The evidence agent itself runs on Llama 3.3 70B, DeepSeek V3.1, Gemini 2.5 Flash and
GPT-4o-mini through OpenRouter.
