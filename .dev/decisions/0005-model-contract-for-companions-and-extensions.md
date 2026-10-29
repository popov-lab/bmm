# 0005: The model contract for companion and extension packages

- **Status:** proposed
- **Date:** 2026-10-29
- **Issues:** #409 (item 7), #442, #554; bmmtools #10 and #22
- **Amends:** [0001](0001-release-scope-1.4-to-1.6.md) § 1.6.0 (the
  simulation generic)

## Context

Two packages outside bmm depend on its models today. bmmtools (parameter
recovery, simulation-based calibration, validation reports) and bmm-explorer
(an interactive page per model) both need to know what a model is: its
parameters, links, fixed values, default priors, the columns it reads, and
how to simulate data from it. More readers are coming: #442 lets other
packages define bmm models, and those models must be described and simulated
the same way.

Measured on `develop` (f04fd0fb and 4ae8bfff, 2026-10-29):

- Every model object carries the metadata the readers need (`resp_vars`,
  `other_vars`, `domain`, `task`, `name`, `citation`, `version`,
  `requirements`, `parameters`, `links`, `fixed_parameters`,
  `default_priors`; `init_ranges` in 12 of the 14 constructors,
  `alternative_parameters` in 2). No exported function returns it.
  `parameter_info()` shows part of it; `init_ranges` and
  `alternative_parameters` are reachable only with `$`.
- bmmtools reads the object with `$` (`links`, `fixed_parameters`,
  `alternative_parameters`, `parameters`, about 90 reads of `resp_vars` and
  `other_vars`), keeps its own copy of the twelve-link table because
  `link_transform()` is internal, hard-codes `mu1`, and carries one adapter
  per model that maps the model's parameters to the arguments of its `r*()`
  generator. bmm 1.4.1's `criterion`/`bias` alternative parameters broke its
  `sdt_yn` adapter.
- bmm-explorer hard-codes each page's parameter list and prior ranges. Its
  fixture generator calls 25 internal functions; that generator is a test
  oracle for the numerics and stays outside this record.
- The 16 `r*()` generators follow three argument conventions: natural-scale
  names (`p_mem`, `p_nt`) where the model says `thetat`, `thetant`; `n`
  against `n_trials`; and parameter lists (`rm3()`, `rmpt()`). `rimm()` and
  `dimm()` default `c` and `a` differently. No test spans them.
- The lifecycle rules in CONTRIBUTING ("Deprecating and removing features")
  cover exported functions and arguments, parameter names, documented
  defaults, the data layout and saved fits. They say nothing about model-object
  fields or companion packages, and "internal functions are not covered".
- Three changes already scheduled break a reader that does not follow them:
  `mu1` to `mu` and `pnt` to `ptarget` in the rebuilt circular mixture models
  (1.5.0), the `log` default of the density functions (#553, 1.6.0), and the
  removal of `supported_models()` and `parameters()` (1.6.0).
- Nothing in bmm's CI builds or tests a companion package.

Parts of a fix were designed in five places that do not reference each other:
#409 item 7 (a documented accessor, promised for 1.5.0), #442 phases 4 and 7
(metadata in the object, a versioned contract), #554 (one vocabulary for
`domain`, `task` and `name`), the bmmtools requests behind bmmtools #10 and
#22 (export `link_transform()`, a simulation generic), and 0001 § 1.6.0 ("a
simulation generic (new issue)", never filed). This record joins them.

## Decision

### One contract, two readers

bmm defines the model contract once. It has two readers:

- **Companion packages** read it: the accessor, the simulation generic, the
  exported densities and generators, and `link_transform()`.
- **Extension packages** (#442) write to it: a model they register must give
  the accessor everything it returns for a built-in model, and must have a
  method for the simulation generic.

A built-in model and a registered model are held to the same contract.

### The accessor

`model_spec()` is exported and returns a model's specification as a plain
list. The name `model_info()` is taken: insight exports that generic and bmm
registers `model_info.bmmfit()` for it.

### The stable surface

The contract is exactly what `model_spec()` returns:

- name, package, version, domain and task (the #554 vocabulary);
- parameters, each with its description, link, fixed value (if any) and
  alternative group (if any);
- default priors (`main`, `effects`, `sd`);
- init ranges;
- response variables and other variables;
- data-column roles;
- the names of the model's density and generator functions;
- the contract version.

The lifecycle rules apply to these items as they apply to exported functions.
Reading the model object with `$`, its attributes, and internal functions stay
outside the contract. Companion packages move to the accessor.

bmm ships no JSON writer for the specification. A reader that needs JSON
serialises the list itself; bmm-explorer already runs R in its fixture
oracle.

### The simulation generic

`simulate_data(model, pars, design = NULL, scale = c("link", "native"), ...)`
is an exported S3 generic on `bmmodel`. It takes parameter values under the
model's own names, on the link or the native scale, and a design (trial
counts, set sizes, non-target locations, conditions), and returns a data frame
that `bmm()` accepts with the same model. Methods wrap the existing `r*()`
functions, as `rm3()` and `rmpt()` already take parameter lists. The `r*()`
signatures do not change.

A method for `stats::simulate()` was considered. Its signature is
`(object, nsim, seed, ...)`, with no place for parameters or a design, so the
whole interface would pass through `...`. dRiftDM and ggdag export a
`simulate_data()` generic of their own; when one of them is attached after
bmm it masks bmm's, so the documentation calls it as `bmm::simulate_data()`
where that matters.

Release: the generic lands in 1.5.0 if the circular mixture stack (#434) has
merged into `develop` by 1 Dec 2026, otherwise in 1.6.0 on `dev-1.6.0`. This
amends 0001 § 1.6.0, which placed it in 1.6.0 without condition.

### `link_transform()`

`link_transform()` is exported with its twelve element-wise links. `softmax`
is not one of them: it needs the whole vector of a parameter group, and stays
on its own path.

### The contract version

`model_spec()` reports a contract version. It becomes `1.0` with the 1.6.0
release, where #442 phase 7 pins it (`bmm_contract_version()` and a test of
the exported generics and their arguments). Before that release the contract
may still change, under the announcement rule below.

### The announcement rule

A pull request that changes or removes a contract item:

1. adds a NEWS bullet under a heading `### Companion packages` in the section
   of its release;
2. files one issue in each companion repository (GidonFrischkorn/bmmtools,
   GidonFrischkorn/bmm-explorer), from a template, with the old and the new
   form side by side, before the release that contains the change.

This extends rule 3 of CONTRIBUTING ("Announce breaking changes") from what
users type to the contract items. The first applications are `mu1` to `mu`
and `pnt` to `ptarget` (1.5.0), the `log` default (#553, 1.6.0), and the
removal of `supported_models()` and `parameters()` (1.6.0).

### The downstream job

A workflow in popov-lab/bmm checks out GidonFrischkorn/bmmtools and
GidonFrischkorn/bmm-explorer and runs their checks against the bmm that was
pushed. It runs on pushes to `develop` and `dev-*`, weekly, and on manual
dispatch; not on pull requests, because the stacked model PRs would start it
many times. It never blocks a merge. The other maintainer confirms the
reference to repositories outside the organisation before the workflow
merges. The companions' own weekly jobs against `develop` stay in place.

### Recovery evidence

bmmtools validation cards are the one home for parameter-recovery results
(#409 item 6). The website's "Parameter recovery" link points to the card site
once it is published, and the old recovery site is archived. A new-model pull
request links a card or attaches a `validation_plan()` file; CONTRIBUTING, the
pull request template and `CLAUDE.md` say so.

## Consequences

- Open model pull requests add a `simulate_data()` method at their next base
  merge once the generic is on their base, and take the #554 strings for
  `domain`, `task` and `name`: the circular mixture levels (#583 to #590),
  TCC (#591 to #593, #582), the psychometric function (#507), ezcdm (#524) and,
  on `dev-1.6.0`, the race models (#500 to #502, #565). The convention test
  names what is missing, as the constructor-dots test of #552 does.
- `use_model_template()` scaffolds a `simulate_data()` method, and the
  model-cases registry gains a `simulate` case kind: one more convention test
  (ADR 0004, R4).
- bmmtools removes its adapters and its copy of the link table once the
  installed bmm has the generic and the export; until bmm 1.5.0 is on CRAN it
  keeps both paths behind a version check. bmm-explorer reads each page's
  parameters and prior ranges from the specification instead of typing them.
- The downstream job costs runner time on every push to `develop` and
  `dev-*`, and someone has to read its failures. A failure that comes from an
  intended contract change is answered by the announcement issue; any other
  failure is a bug in bmm or in the companion.
- `model_spec()` and the generic join the exported API, so their pull
  requests take the pipeline-or-API review tier.
- Open points: whether `model_spec()` also accepts a fitted model (decided
  in the accessor's pull request); the 1 Dec gate for the generic; the other
  maintainer's confirmation of the downstream job; the `package` and `version`
  fields for registered models, which #442's registry fills when `develop` is
  forward-merged into `dev-1.6.0`.
