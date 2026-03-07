# Changelog

## 2026-03-07 - Mentorship Documentation Pass

### Why we made this pass
- We wanted the docs to teach the workflow, not just describe it.
- We rewrote explanations so junior developers and wet-lab collaborators can troubleshoot confidently.
- We aligned all high-level guidance to the same five-station mental model.

### What changed
- changing names to clarify processes
- Rewrote `README.md` around:
  - `Scanner`
  - `The Starting Line`
  - `The Race Begins`
  - `Security Check`
  - `The Winning Bunch`
- Added empathetic troubleshooting tips with plain-language fixes.
- Updated inline comments/docstrings in core modules to explain **why** each step exists:
  - `src/sequence_generator.py`
  - `src/binding_scorer.py`
  - `src/target_analyzer.py`
  - `src/filter_rank.py`
  - `src/pipeline.py`

### Behavior impact
- No intentional algorithm changes from this documentation pass.
- Primary logic remains empirical: real count ingestion, CPM normalization, enrichment scoring, hard-stop target verification, and ranked shortlist output.
