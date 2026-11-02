# Keeping this repository traceable

The initial publication uses the repository's existing default branch, `Dosage-Compensation`. STAR and Telescope source snapshots are versioned separately. Keep historical snapshots intact; place changed scientific workers/configuration in a clearly named new version and document why recounting or remapping is required.

For later local maintenance, the command sequence is:

```bash
git clone https://github.com/Lubna-Younas/Drosophila-TE-landscape.git
cd Drosophila-TE-landscape
git switch Dosage-Compensation
git pull --ff-only
git switch -c maintenance/describe-the-change

# Edit the relevant code/docs; preserve original execution snapshots.
python3 scripts/validate_repository.py
python3 scripts/run_synthetic_checks.py
git diff --check
git diff --stat
git status --short

# Select reviewed files explicitly; do not add raw data or private manifests.
git add README.md docs/ scripts/ metadata/ hpc/ .github/ .gitignore
git diff --cached --check
git diff --cached --stat
git commit -m "Describe the concrete workflow or documentation change"
git push -u origin maintenance/describe-the-change
```

Open a pull request into `Dosage-Compensation`. Summarize the changed behavior, affected annotation/assembly identities, and validation. A changed execution snapshot needs new source manifest entries with both original and published hashes; do not simply overwrite the recorded hashes to conceal a modification. Validate links and check staged content for personal paths, account/host details, credentials, sequencing data, and generated outputs before publishing.

The initial release was uploaded through the authenticated GitHub connector, not by pretending these local `git push` commands had already run. GitHub Actions performs static/metadata checks only; the repository does not automatically submit HPC jobs or rerun scientific analyses.
