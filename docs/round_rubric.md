# This project's map, for an adversarial round's attacker

**Declared to the shared harness as `[paths] rubric` in `round.toml`, and spliced into every attacker's
brief.** Rules 2 and 3 of the harness's rubric make the argument: check the literature this project
holds before asserting anything about prior work, and check a statistical claim against this project's
own standards. They name no paths, because the harness serves more than one repository. This file says
where ours are. Editing it changes every brief generated afterwards, so do it when no round is in flight.

## The literature (rule 2)

`refs/catalog_index.json` lists the papers this project draws on: slug, citation and DOI. It is
rendered from a separate literature library and **scoped to the shelves this project cites**
(genome alignment, bioinformatics tooling, statistical model selection and inference, software
testing). A paper missing from it may still be held by the library under another shelf; if a claim
needs one, say so rather than assuming it is absent. `refs/LIBRARY_VERSION` is the library revision
the index was rendered from.

Each paper has an adversarially fact-checked summary in the library. Read it **at the pinned
revision**, never at whatever the library's working tree has checked out, so another repository's
state cannot change this round's verdict:

    git -C ../library show "$(tail -1 refs/LIBRARY_VERSION)":papers_md/summaries/<slug>.md

When a summary is silent on the number you need, open the full text. `../library/papers_md/<slug>.md`
is a text conversion of the held PDF; it is not versioned, so there is nothing to pin.

**The library is not public.** Without a checkout at `../library`, the summaries and full texts are
unavailable: say in your findings that you worked from the catalog index alone, and do not assert what
a paper says from memory.

`manuscript/references.bib` is what the manuscript itself cites. A held paper in the index whose result
would change a claim here, and which the manuscript does not cite, is a prior-art finding.

## The method cards (rule 3)

This project keeps **no method cards** yet. Say so in your findings rather than searching for them, and
check statistical claims against the cited papers instead.

## The claims already adjudicated

`scripts/round history "<subject>"` searches every attack, adjudication and completeness critic across
every closed round and prints the verdict where one is readable. Check it before you file.
