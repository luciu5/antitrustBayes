# Release dependency pins

`release-dependencies.json` is the automatic dependency source for `release/**`
pushes and pull requests. Its schema is a flat JSON object with exactly three
keys: `antitrust`, `coordination`, and `trade`. Each value must be the full
40-character commit SHA for that package. The workflow validates the key set
and every SHA before fetching; it checks out each dependency outside the package
build tree and records the resolved SHA.

The checked-in values are previous-manifest placeholders. Replace them with the
final candidate commit SHAs after those candidates have been committed. A
release certification started with `workflow_dispatch` remains separate and
requires all three full SHA inputs; those inputs are checked against the
checked-out commits. Ordinary development pushes and pull requests keep their
existing `master`/`main` bootstrap refs.

## Structural dependency metadata regression

Release-pin CI runs [`test-verified-dependency-metadata.R`](test-verified-dependency-metadata.R)
after installing the selected dependency checkouts. The script extracts both
R verification heredocs from `qa.yml` and executes them against those source
checkouts and the installed packages. It covers the `antitrust` source and
installed metadata with no `Remotes` field, wrong `antitrust`/`coordination`
pins in `coordination`/`trade`, and source/installed version mismatches. To run
the same check locally, select a library containing the installed package set
and pass the three source directories followed by the selected full SHAs:

```sh
Rscript --vanilla ci/test-verified-dependency-metadata.R \
  /path/to/antitrust /path/to/coordination /path/to/trade \
  ANTITRUST_SHA COORDINATION_SHA
```
