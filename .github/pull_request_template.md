## What and why

<!-- What does this change, why, and how did you check it? -->

## Changelog entry

<!--
One line for CHANGELOG.md, written for a user of the library, or "None" for a
change no user can see (tests, CI, internal refactoring). The maintainers
edit and place it.
-->

## Compatibility

<!--
Your best guess; the maintainers make the final classification
(RELEASING.md). Tick one.
-->

- [ ] No user-visible change
- [ ] Compatible fix or improvement (patch)
- [ ] Compatible addition or deprecation (minor)
- [ ] Needs downstream code or tuning to change (major)

Can this move a filter's output (estimates, covariance, innovation, or which
steps succeed)? If so, which filters and configurations, by how much, and how
was it checked?

<!-- "No", or the details. -->

## Checks

- [ ] Tests cover the change, in single and double precision where it matters
- [ ] User-visible changes are documented (README, usage guide or header comment)
