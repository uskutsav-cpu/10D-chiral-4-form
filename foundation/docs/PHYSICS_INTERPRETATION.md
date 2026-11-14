# Physics interpretation

**STATUS: REQUIRES COAUTHOR CONFIRMATION.** The algebra below is certified; the
physical reading of it is not, and the manuscript marks it accordingly.

## The setting

A chiral four-form gauge field in ten dimensions has a self-dual five-form field
strength. A non-linear theory of such a field is specified by a Lagrangian, or
equivalently by a flow equation, built from Lorentz invariants of `F`. In four
and six dimensions the available invariants are few enough that a completely
general deformation can be written down term by term. In ten dimensions there
are 81 functionally independent invariants, so that route is closed, and progress
has come instead from specific constructions --- above all stress-tensor flows.

## The question this answers

A construction generates *some* invariants. Which ones, and what does it miss?

Until now that question could not be posed precisely at degree ten, because the
degree-ten space itself was not known. Twelve explicit structures were available,
with no statement that they were independent, and no statement that they were
complete.

## The answer, at degree ten

    dim A10 = 14      the whole degree-ten space
    dim D10 = 11      what the stress flow reaches
    dim Q10 =  3      what it does not

Three independent directions of degree-ten local structure lie outside the
reachable subspace, and we exhibit representatives of all three. That is the
concrete content: not that the flow is incomplete in some vague sense, but that
its degree-ten deficit is exactly three-dimensional and explicitly presentable.

A second, sharper fact: `P10 ⊂ D10`. Products of lower-degree invariants are all
reachable. So the three missing directions are not products of anything simpler
--- they are genuinely new degree-ten structure, and no amount of composing
lower-degree results will produce them.

## What a theorist can now do that they could not before

Concretely, and this is the answer to "what new calculation becomes possible":

1. **Expand any degree-ten invariant in a certified basis.** Given a candidate
   term arising from any source --- an effective-action computation, an amplitude
   matching exercise, a deformation ansatz --- its coordinates in `A10` can be
   computed exactly, and its class in `Q10` read off. Previously a candidate could
   only be compared against the twelve published expressions, and a failure to
   match was uninformative because that list was not known to be complete.

2. **Decide reachability.** Whether a proposed degree-ten term can arise from the
   stress flow at all is now a finite exact computation: does its class in `Q10`
   vanish?

3. **Parameterise the deficit.** A theorist wanting to extend the flow, or to
   consider a more general construction, now knows exactly how many new
   parameters are needed at degree ten (three) and has explicit representatives
   to attach them to.

4. **Avoid a specific error.** The published twelve are not a product complement
   and their non-product content is eleven, not twelve. Anyone counting
   "primitive" degree-ten structures from the published list would be off by one.

## The arithmetic caveat that attaches to items 2 and 3

Items 1 and 4 rest on `dim A10`, `dim B10`, `dim G10` and `dim P10`, none of
which is exposed to an exceptional prime: each space is spanned by exactly as
many explicit invariants as its modular rank, which pins the dimension over `Q`
with no prime excluded.

Items 2 and 3 rest on `D10`, which is not in that position. It is assembled by
admitting a direction only when it raises the rank modulo `p`, so `dim D10 = 11`
is a lower bound over `Q` and `dim Q10 = 3` is an upper bound. Read carefully,
that means:

- **Item 2 is safe in the direction that matters.** A term whose class in `Q10`
  is non-zero is genuinely unreachable, because a bad prime could only make the
  flow reach *more*, and a class that survives a larger `D10` survives. A term
  whose class vanishes is reachable modulo `p`; that is the weaker reading.
- **Item 3's count is "three, or fewer".** A theorist extending the flow needs at
  most three new degree-ten parameters, and exactly three if 32717 and 32749 are
  both good primes.

This is stated so that nobody builds on "exactly three" where "at most three" is
what has been established. Discharging it needs exact evaluation over `Z`; see
PO-09.

## What is deliberately NOT claimed

- No Type IIB coefficient.
- No supersymmetric completion.
- No causality or unitarity statement.
- No amplitude.
- No uniqueness statement.
- No all-orders theorem.
- Nothing about degree twelve beyond that it is unclassified.

None of these follows from the algebra computed here, and the manuscript says so
in its limitations section.

## Honest assessment of significance

This is a structural and computational result, not a new physical prediction. Its
value is that it converts a question that was previously unanswerable --- "is this
degree-ten term reachable, and if not, what is missing?" --- into a finite exact
computation. Whether that is interesting enough for the intended venue is a
judgement for the coauthors, and the manuscript does not assert it.
