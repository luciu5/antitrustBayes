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
