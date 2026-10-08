# Decision records

This folder records development decisions for bmm: what we decided, why, and
what follows from it. The folder is excluded from the package build
(`^.dev$` in `.Rbuildignore`).

The public roadmap lives in a pinned GitHub issue and in the milestones. These
records hold the reasoning behind it.

## Format

One file per decision, numbered in order: `NNNN-short-title.md`. Each record
has:

- **Status**: proposed, accepted, or superseded by NNNN
- **Date**
- **Context**: the situation that forced a decision
- **Decision**
- **Consequences**: what changes, what it costs, and open points

Records are not rewritten after acceptance. A changed decision gets a new
record, and the old one is marked superseded with a link to the new one.

## Index

| No. | Decision | Status |
|---|---|---|
| [0001](0001-release-scope-1.4-to-1.6.md) | Release scope for 1.4.0, 1.5.0 and 1.6.0 | accepted 2026-10-04 |
| [0002](0002-release-integration-branches.md) | Integration branches for the next two releases | accepted 2026-10-04 |
| [0003](0003-model-fit-files.md) | Large model-fit files: keep the history, stop the growth | accepted 2026-10-04 |
