# Contributing to bmm

The `bmm` package is designed to be community driven and thus we
appreciate your input! We attempt to make contributing to `bmm` as easy
and transparent as possible. These contributor guidelines are designed
to clarify different types of contributions and how they will be
acknowledge in the package publication.

The following contributions will be acknowledge in the `NEWS` for each
release of the package:

- Reporting a bug
- Submitting a bug fix
- Discussing and proposing improvements on the current state of the
  code, in particular:
- improving stability of model estimation
- speeding up model estimation
- generalizing functions across multiple already implemented `bmmodels`

More extensive contributions will be acknowledge by being listed as a
contributor in the package documentation. These entail:

- adding a new model to `bmm`
- adding new function that ease the use of multiple models implemented
  in `bmm`
- adding functions that implement the communication of `bmm` with other
  R packages, such as `emmeans`, `tidybayes`, `bayesplot`, etc.

## Package Development on Github

We use Github to host all code, track issues and feature requests, as
well as accept pull requests. Detailed info on the development process
can be found in the `bmm` [Developer
Notes](https://venpopov.com/bmm/dev/dev-notes/index.html)

All changes and additions to code have to be submitted via pull
requests. We consider all pull requests as propositions for changes to
the codebase. Thus, we reserve the right to not merge changes, if we
feel that they are in conflict with general principles implemented in
the `bmm` package. To avoid the rejection of pull requests, please
consider contacting us before committing changes to `bmm`, especially if
they involve changing a large number of files.

If you are interested in contributing to `bmm`, please follow the
following steps.

1.  Fork the repo and create your branch from `develop`.
2.  If you’ve added code that should be tested, add tests.
3.  Update the documentation for the changes you implemented.
4.  Ensure that all unit tests passed.
5.  Ensure the R CMD Checks passed.
6.  Issue that pull request!

## Use a Consistent Coding Style

Please follow the general coding style used throughout `bmm`. This
entails:

- labeling variables and functions using `snake_case`
- avoid upper case labels in variable and function names

The `bmm` [Developer
Notes](https://venpopov.com/bmm/dev/dev-notes/index.html) provide an
introduction into the file organization of the package. Please follow
the guidelines where to put functions associated with the different
steps in fitting `bmmodels`. If you have questions or an unsure about
where to add code, feel free to ask us. There is a dedicated
[Discussion](https://github.com/venpopov/bmm/discussions) page for
informal chats and questions.

## Any contributions you make will be under the GPL-2 Software License

In short, when you submit code changes, your submissions are understood
to be under the same
[GPL-2](https://choosealicense.com/licenses/gpl-2.0/) that covers the
project. Feel free to contact us if that’s a concern.

## Report bugs using Github’s [issues](https://github.com/venpopov/bmm/issues)

We use GitHub issues to track public bugs. Report a bug by [opening a
new issue](https://github.com/venpopov/bmm/issues); it’s that easy!

## Write bug reports with detail, background, and sample code

When reporting a bug, please use the provided template.

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can.
- What you expected would happen
- What actually happens
- Notes (including why you think this might be happening, or stuff you
  tried that didn’t work)

We appreciate thorough bug reports *a lot*.
