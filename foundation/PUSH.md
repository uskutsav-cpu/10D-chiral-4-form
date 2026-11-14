# Getting this onto GitHub

I can't create the repo for you (no access to your GitHub account), so:

1. On github.com click **New repository**, name it `selfdual-5form-invariants`,
   leave it empty (no README, no .gitignore).

2. Then from the unzipped folder:

```bash
git init
git add .
git commit -m "verified core: 6D reproduces 5 invariants, 10D order 4 gives 1"
git branch -M main
git remote add origin https://github.com/<you>/selfdual-5form-invariants.git
git push -u origin main
```

3. Add David as a collaborator (Settings -> Collaborators).

## Working agreement worth setting now

- **Branches and PRs from here on.** No more direct pushes to main. Five
  commits for a whole project means neither of you can review or bisect.
- **`pytest tests/` must pass before any merge.** That suite exists because
  the same overflow bug produced confident wrong answers three times.
- **Split the work at a clean interface:** one of you owns candidate
  generation (`graphs.py`), the other owns targets and validation (the
  Hilbert-series degree list, and the rank checks that say whether a level
  is complete). Neither blocks the other.
