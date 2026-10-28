# Contributing to GORILLA

## Questions and bugs

For questions about any aspect of the software, contact the main developer
Michael Eder (<eder@tugraz.at>) or open an issue on the
[issue tracker](https://github.com/itpplasma/GORILLA/issues). The tracker is
also the place to report bugs, propose features or discuss the design.

## Code contributions

Contributions are made through GitHub pull requests. External contributors
should use the ["fork and pull"](https://guides.github.com/activities/forking/)
workflow; core developers use branches in the main repository.

1. Open an issue to discuss the proposed contribution.
2. Implement the change in your fork or branch.
3. Open a pull request against `main`, where the change is reviewed before it
   is accepted.

## Continuous integration and merging

- Every pull request, including each layer of a stacked pull request, runs the
  fast CI tier. Add the label `full-ci` to run the full tier before queueing.
- Review a stack layer by layer. When it is ready, retarget the top layer to
  `main` so that the stack is merged as one pull request; the lower pull
  requests then show as merged.
- Merge with "Merge when ready" (or `gh pr merge --merge --auto`). The pull
  request enters the merge queue, the full CI tier runs on the exact merge
  result, and on success it lands on `main` as a merge commit. `main` is not
  re-tested on push.
- A branch does not need to be up to date with `main`. Merge `main` into a
  pull request only when it conflicts.
