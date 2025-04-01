# Git repository cleanup guide

This guide explains how to inspect and clean a GitHub repository by identifying large files, removing them with `git filter-repo`, deleting leftover empty folders, and force-pushing the cleaned repo.


## 1. Print the heaviest files in the repository history

Run this command to print the 20 largest files ever committed (size in MB):

```bash
git rev-list --objects --all | \
  git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' | \
  grep '^blob' | \
  sort -k3 -n -r | \
  head -n 20 | \
  awk '{printf "%.2f MB\t%s\n", $3/1024/1024, $4}'
```


## 2. Remove unwanted files or folders using `git filter-repo`

General syntax to remove files or directories from the full Git history:

```bash
pip install git-filter-repo
git filter-repo \
  --force \
  --invert-paths \
  --path <exact/path/to/remove> \
  --path-glob '<pattern/to/match>'
```

Examples:
- `--path _static/` removes a specific directory.
- `--path-glob '*.mp4'` removes all `.mp4` files.
- You can specify multiple `--path` and `--path-glob` options.

> ⚠️ This rewrites history. Make a backup or work in a disposable clone if needed.

> ⚠️ If there are issues about empty folders not being deleted you can run the following command in a separate terminal to speed up the process

```bash
while true; do find .git/objects -type d -empty ! -name info ! -name pack -print -delete; sleep 0.2; done
```

---

## 4. Push the cleaned repository back to GitHub

After filtering the repository, you need to reconfigure the remote and force-push:

```bash
# Remove and re-add origin (if needed)
git remote add origin https://github.com/your-username/your-repo.git

# Force-push cleaned history
git push --force origin main

# Force-push all tags
git push --force origin --tags
```


## 5. Prevent re-adding build artifacts

Add a `.gitignore` file to exclude unwanted files in the future:
