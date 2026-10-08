#!/bin/bash
# Assemble the curated GitHub release repo, the full-provenance bundle,
# and the downloadable zip. Idempotent: wipes and rebuilds github_release
# content (except README.md / UPSTREAM.md / .gitignore which are authored).
set -euo pipefail

ROOT=/home/z/my-project
REL=$ROOT/github_release
DL=$ROOT/download

# ---- 1. Curated release content -------------------------------------------
mkdir -p $REL/papers $REL/research_notes $REL/sources/flagship $REL/sources/companion $REL/scripts

# Final PDFs (v1 and v2 — v2 is the audit revision; v1 kept, never overwritten)
cp $DL/lonely_runner_type_mismatch_paper.pdf                       $REL/papers/
cp $DL/lonely_runner_pair_sum_lattices_conditional_theorem_paper.pdf $REL/papers/
cp $DL/lonely_runner_type_mismatch_paper_v2.pdf                   $REL/papers/
cp $DL/lonely_runner_pair_sum_lattices_conditional_theorem_paper_v2.pdf $REL/papers/
cp $DL/lonely_runner_type_mismatch_paper_v3.pdf                   $REL/papers/
cp $DL/lonely_runner_pair_sum_lattices_conditional_theorem_paper_v3.pdf $REL/papers/
cp $DL/flagship_audit_adjudication.md                              $REL/papers/
cp $DL/audit_adjudication_addendum_m10.md                           $REL/papers/
cp $DL/companion_audit_adjudication.md                             $REL/papers/

# Research notes (the 20 dated .md notes, not the README)
for f in $DL/*.md; do
  b=$(basename "$f")
  [ "$b" = "README.md" ] || cp "$f" $REL/research_notes/
done

# Flagship sources (v1 + v2)
cp -r $DL/paper2_sources/. $REL/sources/flagship/
cp -r $DL/paper2_sources_v2/. $REL/sources/flagship_v2/
cp -r $DL/paper2_sources_v3/. $REL/sources/flagship_v3/
cp $DL/paper2_sources_v3/merge_cover.py $REL/sources/flagship_v3/ 2>/dev/null || true
cp $ROOT/scripts/paper2/merge_cover.py $REL/sources/flagship/ 2>/dev/null || true
cp $ROOT/scripts/paper2_v2/merge_cover.py $REL/sources/flagship_v2/ 2>/dev/null || true

# Companion sources (v1 + v2)
cp -r $DL/paper_sources/. $REL/sources/companion/
cp -r $DL/paper_sources_v2/. $REL/sources/companion_v2/
cp -r $DL/paper_sources_v3/. $REL/sources/companion_v3/
cp $DL/paper_sources_v3/merge_cover.py $REL/sources/companion_v3/ 2>/dev/null || true
cp $ROOT/scripts/paper/merge_cover.py $REL/sources/companion/ 2>/dev/null || true
cp $ROOT/scripts/paper_v2/merge_cover.py $REL/sources/companion_v2/ 2>/dev/null || true

# The adjudication verification run
cp $ROOT/scripts/audit_verify_v2.py $REL/scripts/
cp $ROOT/scripts/out_audit_verify_v2.json $REL/scripts/

# Research scripts + persisted outputs (exclude build dirs + pycache)
rsync -a --exclude 'paper/' --exclude 'paper2/' --exclude '__pycache__/' \
      --exclude '*.pyc' $ROOT/scripts/ $REL/scripts/ 2>/dev/null || {
  # rsync fallback with cp
  (cd $ROOT/scripts && find . -maxdepth 1 -type d ! -name paper ! -name paper2 \
    ! -name __pycache__ ! -name . -exec cp -r {} $REL/scripts/ \;)
  (cd $ROOT/scripts && find . -maxdepth 1 -type f ! -name '*.pyc' -exec cp {} $REL/scripts/ \;)
}

# License from upstream repo
cp $ROOT/lonely-runner/LICENSE $REL/LICENSE

# ---- 2. Init the release repo (single clean commit) ------------------------
cd $REL
if [ -d .git ]; then rm -rf .git; fi
git init -q -b main
git config user.email "z@container"; git config user.name "Z User"
git add -A
git commit -q -m "Curated release: LRC pair-sum-lattices program (flagship 20pp + companion 64pp, sources, scripts, provenance outputs)"
echo "release repo: $(git rev-parse --short HEAD)  files: $(git ls-files | wc -l)  size: $(du -sh --exclude=.git . | cut -f1)"

# ---- 3. Full-provenance bundle of the working repo --------------------------
cd $ROOT
git bundle create $DL/lrc_full_provenance.bundle main --all
echo "bundle: $(du -sh $DL/lrc_full_provenance.bundle | cut -f1)"

# ---- 4. Downloadable zip of the release -------------------------------------
cd $(dirname $REL)
rm -f $DL/lrc_github_release.zip
zip -qr $DL/lrc_github_release.zip $(basename $REL) -x "*/.git/*"
echo "zip: $(du -sh $DL/lrc_github_release.zip | cut -f1)"
