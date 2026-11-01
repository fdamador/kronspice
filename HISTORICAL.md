# Gabriel Kron and universal engineering

kronspice is named after Gabriel Kron (1901–1968), and its `kron` package
implements his methods. This file records where those methods came from and
what Kron meant to do with them.

"Universal engineering" appears to be Kron's name for the aim behind his work:
one method for every engineering system, built on geometry and networks. The
method started with electrical machines, then spread to networks, fields,
elastic structures, nuclear reactors and quantum mechanics. It ended in an
unfinished theory of multi-dimensional networks.

**How far the evidence goes.** The phrase is Kron's: it appears in the titles
of two of his 1956 articles. We have not read either article, so the account of
the idea below is reconstructed from his other writing and from people who knew
him. [Sources](#sources) separates what was read directly from what is known
only from citations or excerpts.

## Timeline

| Year | Event |
|---|---|
| 1901 | Born in Nagybánya, Hungary (later Baia Mare, Romania) |
| 1921 | Reaches New York with his brother Joseph; both enter the University of Michigan that September |
| 1926–28 | Travels around the world, largely on foot, studying mathematics along the way |
| 1930 | First paper, "Generalized Theory of Electrical Machinery" |
| 1931–33 | Warner Brothers contract; his department closes and he studies tensor analysis in Romania while still being paid |
| 1934 | "Non-Riemannian Dynamics of Rotating Electrical Machinery" published; joins General Electric in May |
| 1935 | Montefiore Prize, University of Liège |
| 1939 | *Tensor Analysis of Networks* |
| 1942 | *A Short Course in Tensor Analysis for Electrical Engineers* (*Tensors for Circuits*, Dover, 1959) |
| 1944–45 | Equivalent circuits for Maxwell's equations, differential equations and the Schrödinger equation |
| 1951 | *Equivalent Circuits of Electric Machinery* |
| 1953–55 | First tearing papers |
| 1956 | The two "universal engineering" articles |
| 1957–59 | *Diakoptics* serial in *The Electrical Journal* (London); book 1963 |
| 1959–65 | Papers on multi-dimensional space filters, self-organizing automata, and networks beyond graphs |
| 1966 | Retires from General Electric in July |
| 1968 | Dies in Schenectady on March 25 |
| 1969 | Union College symposium on October 14; papers published in 1973 as *Gabriel Kron and Systems Theory* |

## Life

Kron grew up in Transylvania, the youngest of eight children of a dry goods
merchant. After the Romanian takeover of the region closed his Gymnasium in
1918, he decided his future lay abroad. He tutored his brother Joseph through
the full Gymnasium course in under a year, and the two reached New York in
January 1921 (Alger 1973).

After Michigan, Kron designed induction motors in Los Angeles, then moved to
Robbins and Myers in Ohio. Once his citizenship came through in 1926, he left on the
trip he had planned, through the Pacific, Australia and Asia to Romania. He
read differential equations and vector analysis along the way. He later wrote
that on the walks through Queensland and Asia, engineering structures began to
appear to him as

> a collection of numerous multi-dimensional spaces connected together into
> one unit

in the way that countries are linked by roads, customs and laws. If the links
vanish, he wrote, nothing physical is lost, only what turns the separate parts
into a single whole. Only later did he find that mathematicians had built such
a calculus under the name of tensor analysis (Kron, *Life and Times*, p. 161,
quoted in Alger 1973).

General Electric took him on in 1934, after a conference that Philip Alger
arranged to air his ideas. He stayed for 32 years under fifteen managers,
moving between machines, power systems, turbine stress analysis, electronics
with Simon Ramo, and reactor piles. Most of his theoretical writing was done
at night and on weekends (Alger 1973).

## How the method developed

### The generalized machine (1930–38)

Kron's starting point was that every type of electric machine must be a
special case of one general machine, and that understanding it would lead to
new types (Alger 1973). His 1934 paper recast rotating machinery in the
language of non-Riemannian geometry.

Banesh Hoffmann summarized the result in 1958. Kron builds prototype
"primitive machines", and obtains all other machines from them by tensor
transformation. His transformations are often singular and sometimes
non-holonomic, which is where much of the criticism fell (*Tensors for
Circuits*, pp. ix–x). P. Le Corbeiller called it a method "in which one and
the same tensor equation applies to every conceivable type of machine" (quoted
by Hoffmann, pp. viii–ix).

### Networks as geometry (1935–42)

*Tensor Analysis of Networks* applied the same idea to static networks. Start
from the primitive network, with every branch unconnected. Interconnection is
then a change of reference frame, given by a connection tensor `C`, and the
equations of the connected network follow from the transformation laws:

```
i = C·i'    Z' = Cᵀ·Z·C    e' = Cᵀ·e
```

In the early 1930s Kron had also looked for the dual of Maxwell's mesh method.
This produced the junction-pair concept and the idea of the "open mesh", now
the open path (Happ 1971, pp. xi–xii).

### Equivalent circuits for fields (1944–54)

Kron then turned the method on continuous systems, modeling them as circuits
that could be set up on a network analyzer:

- electromagnetic fields, for antennas, waveguides and cavities (1944)
- ordinary and partial differential equations, and the elastic field (1945)
- the Schrödinger equation (1945)
- the nuclear reactor (1954)

### Tearing, or diakoptics (1953–63)

A problem of losses in interconnected power systems led Kron to a way of
solving large systems piecewise (Happ 1973). The method tears a system into
subsystems, solves each one, then corrects the solutions for their
interconnection, and the result is exact for the whole. It appeared in the
*Journal of Applied Physics* (1953), *Proceedings of the IRE* (1954) and
elsewhere, including an inversion of a 256 × 256 matrix by sections
(*Engineering*, listed under 1955). Kron then applied it beyond circuits, to
elastic structures (1955) and to nonlinear plastic structures (1956).

The foundation came partly from others. J. Paul Roth gave tearing an
algebraic-topology basis (1955, 1959), and Happ built the contour theory and
multilevel algorithms, because, in Happ's words, the theoretical foundation
"still remained to be developed" (1971, p. xii). Kron did not accept that
textbook topology was enough. In 1962 he called algebraic topology "a false alarm
for electrical networks", said he believed tearing was absent from topology texts, and told
electrical engineers to do their own topological research.

### Polyhedra and n-networks (1958–68)

Work on forecasting and optimization led Kron toward networks of higher
dimension, which his 1968 obituary calls "seemingly far-afield". His clearest
statement that we have read is "Camouflaging Electrical 1-Networks as Graphs"
(1962), which makes two connected claims.

**An electric network is not a graph.** Before excitation, Kron wrote, a network
may be viewed as a polyhedron: nodes, branches, triangular planes and so on, all
interconnected. Once excited, only the branches count, so the proper model "is
not a graph, nor a planar polyhedron, but a 1-network", a configuration of
branches with no nodes or planes. He argued that Kirchhoff's current law rules
out nodes and his voltage law rules out planes. He also separated the connection
matrix from the incidence matrix, and held that the first cannot be derived from
the second.

**Higher-dimensional networks carry waves, not currents.** Every result for a
1-network carries over to a p-network of any dimension, by what he called a
"generalization postulate". A sequence of n + 1 p-networks, from a 0-network to
an n-network, can be connected into an n-dimensional polyhedron, with incidence
matrices linking neighboring networks. An electric current cannot cross such a
structure, but an electromagnetic wave can, and those wave quantities satisfy
Stokes' theorem at the boundaries. He proposed these excited polyhedra for
multi-dimensional curve fitting, partial differential equations, Fourier
transforms and information-theoretical problems.

**The oscillatory polyhedron as a "self-organizing automaton".** The
companion paper, "Multidimensional Curve-Fitting with Self-Organizing Automata"
(*J. Math. Anal. Appl.*, 1962), says what the excited polyhedron was for. A
polyhedron and its dual carry a sequence of multidimensional waves; immersed in
a stationary or moving magnetohydrodynamic plasma, "the amorphous field
crystallizes into a sequence of 2k sets of transmission networks", coupled and
energized by k-dimensional magnetohydrodynamic generators, which he equates
with "generalized" rotating electrical machines (abstract, p. 46). Even without
motion the structure "may assume a self-adaptive 'oscillatory' state" in which
one polyhedron represents any number of functions of k variables together with
their divided differences, "without needing any adjustment or interference by
the analyst" — that is what he meant by a self-organizing automaton. The
worked examples are least-squares regression of several functions of two
variables at four non-uniformly spaced points; the oscillatory model fits six
arbitrary functions at once, while non-oscillatory polyhedra fit only one or a
small class (p. 46). It even fits a set of random numbers, of which he wrote
"the significance of this is not yet understood" (p. 69). In the introduction he
placed the whole program himself within the integration theory of exterior
differential forms of Cartan, de Rham, Hodge, Whitney and Steenrod, "tensors in
the large" on the polyhedral networks of algebraic topology (p. 47).

**The electromagnetic structure of the model** (pp. 53–54, 66–68). Each
q-network is built from his 1944 Maxwell network (Fig. 1, "satisfying eight
Maxwell's equations simultaneously": four field and four constitutive
equations). Because even- and odd-dimensional spaces behave differently, one
complete wave with e, h, b and d needs two networks of neighbouring dimension,
so a k-dimensional polyhedron carries only k/2 full waves; all even-dimensional
q-networks are magnetic (inductors and resistances), all odd-dimensional ones
dielectric (capacitors and conductances), and a dual polyhedron with the roles
exchanged is interconnected with the first both conductively and inductively,
supplying the ideal transformers. The total power input to the network is the
error-square of the regression, so the best fit is the state of least stored
electromagnetic energy, and a perfect fit stores none, which is why both
polyhedra must be made oscillatory before excitation (p. 67). Kron called the
physical interpretation "only tentative" and wrote that "some aspects of the
physics are still puzzling, but the numerical results speak for themselves"
(p. 66). The "k-dimensional space" of the paper is the space of the k
independent variables of a data matrix (p. 47), not physical space: the
polyhedron is an estimating device, "one engineering tool" (p. 50). The one
bridge to nature is the closing section "Crystals as Self-Organizing Automata"
(pp. 67–68): a crystal's atoms as the vertices of a polyhedron, its X-ray-excited
dipole waves and diffracted field as a self-sustaining system of which the
polyhedral waves are "surprisingly complete analogues", hence the crystal
computer, and a remark that a crystal "may also act as a model for some classes
of neural phenomena".

**Kron's own words, 1958–65, at second hand.** Balasubramanian's 1968
Liverpool thesis (pp. 7–18) quotes at length from six of the late papers we
have not obtained, and is our only text for them. From "Multidimensional space
filters" (1958): the "sequence of alternating magnetic and dielectric networks
may be called a 'Wave-Model'"; the configuration of 1- to k-dimensional volume
elements carrying the waves is a "polyhedron"; two or more wave models,
interconnected or inductively coupled, are a "Multi Dimensional Space Filter";
the waves "are always accompanied by free and bound electric charges as well as
by magnetic poles"; and the dielectric lines and magnetic planes of his 1944
Maxwell circuit form the first complete sequence (thesis pp. 8–9). From "Basic
concepts" (1959): odd-dimensional spaces conduct the dielectric quantities
(solenoidal d, e and lamellar D, E), even-dimensional spaces the magnetic ones
(pp. 10–11). From the General Electric report "Building-blocks of
self-organizing automata – I" (60GL164, 1960): "The stage for all cognitive
activities is an underlying n-dimensional field of magnetohydrodynamic plasma.
Two spatially orthogonal polyhedra are immersed in the fluid in order to
crystallise the amorphous non linear field into a host of interconnected
n-dimensional magnetohydrodynamic generators", which "form the elementary
'neurons' of the projected 'artificial brain'" (pp. 15–17). From "The
misapplication of graph theory" (1962): "electricity demands its own highly
specialised underlying topological structure"; both the polyhedron and its
dual "must be present to propagate an electromagnetic wave", so the complete
1-network is "a 1-network with an environment", and in the two-phase model
the eight electromagnetic parameters appear as the four tensors F_αβ, H_αβ,
s_α, s^α (pp. 11–12). From "The frustrating search" (1963): the discovery of
differential topology's "tensors in-the-large" (Cartan, de Rham, Hodge,
Whitney, Eilenberg, Steenrod, Chern and others); incidence matrices satisfy
Stokes' theorem while connection matrices satisfy Kirchhoff's laws; and the
crystal as "artificial brain" (pp. 12–15). From "Invisible dual
(n-1)-networks" (1965): a region of n-space is represented not by n
one-dimensional axes but by two orthogonal reference axes that are a
p-hyperplane and an (n−p)-hyperplane, so a continuous region has "a 'primal'
p-network plus a 'dual' (n-p) network as its abstract reference frame", a
generalisation of the two-phase frame of rotating machines (pp. 17–18). These
passages, not the 1962 curve-fitting paper, are where Kron states the physical
reading (register entry H20 in CLAUDE.md). The thesis's own verdict: it
"searched, in vain, for some supporting work along the lines of Kron's wave
structure", found that "many details have not been disclosed by Kron himself,
even where he has published numerical results", and "does not reach as far as
Kron's automaton" (pp. 19–20).

Later titles continue the argument, such as "Invisible Dual (n-1)-Networks
Induced by Electric 1-Networks" (1965). The record of that abstract says the
visible and invisible p-networks around the branches together form a
non-Riemannian space rather than a graph or a polyhedron.

This work was never finished. Kron himself wrote in 1962 that "much
organizational and development work still lies ahead", and Happ called the
polyhedral theories "still unexplored" in 1973. Kron saw crystals as the way to
build such structures, with a molecule's atoms as a polyhedron's vertices, and
described "crystal computers" as "within the realm of possibility" (1962, p.
173). After retiring, he meant to write more books, especially on the crystal
computer he thought would be his crowning achievement (Alger 1973). Tolan recalls
being shown what Kron considered proof of its superiority. Our sources document
no working device or independent test.

## What "universal engineering" meant

The two articles carrying the phrase are:

- "Electrical power engineering as a spearhead of 'universal' engineering,"
  *Bulletin of Electrical Engineering Education*, December 1956.
- "Diakoptics — a gateway into universal engineering," *Electrical Journal*
  (London) 157, 1940–1945, 1956.

Their titles suggest a sequence: power engineering leads the way and diakoptics
opens the door to something larger. Five other sources indicate what the larger
thing was.

**Kron's 1939 argument.** Most sciences have equations of identical form, and
"That universal language is geometry". Results stated in that language can be
pooled across sciences. Engineering was worse placed, he wrote, having grown in
isolation into closed specialties, "with all its attendant inconveniences, like
having a different theory for each machine" (*Tensor Analysis of Networks*,
Ch. VIII, pp. 196–197).

**The range of his applications.** Kron carried the same tools (primitive
systems, transformations, equivalent circuits, and later tearing) across
machines, networks, fields, structures and reactors, and wrote on linear
programming in the same journal series.

**The polyhedra.** In 1962 Kron wrote that such structures "can be used for
modelling complex non-electrical problems that have a large number of independent
variables" (p. 162). This is the closest we have to the one structure the
title implies, though he left mechanical, thermodynamic and other parameters "for
future considerations" (p. 167).

**Hoffmann's reading.** "For Kron uses tensors to unify great classes of
physical systems" (*Tensors for Circuits*, p. ix).

**The people around him.** The Japanese Research Association of Applied
Geometry, which Kron influenced, published its memoirs as the *Unifying Study
of the Basic Problems in Engineering Sciences by Means of Geometry*. His
obituary traces to his influence a "trend toward the teaching of engineering
as a unified discipline". Keith Bowden later described his search as one for
"an ontology of engineering" (via Wikipedia, *Diakoptics*).

Together these support a working definition: **universal engineering treats
every engineering system as an instance of one geometric structure of
interconnected networks, so that a method or result carries from one field to
another by transformation.** This is our synthesis, not Kron's wording. His
1956 articles should be checked against it.

## Reception

The reaction was hostile for decades. Mathematicians objected that he
redefined terms and broke established rules (Alger 1973). J. Slepian likened
the work to fruit salad, and Duschek and Hochrainer dismissed it in their tensor
textbook. As Hoffmann observed, critics said both that it was wrong and that it
was right but not new (*Tensors for Circuits*, p. vii).

Some of the criticism was fair, and his defenders said so. Hoffmann conceded
faults of rigor and exposition, including name-dropping such as *Unified Field
Theory*, but still judged him "an innovator of major importance" (p. viii).
Alger reported that some of Kron's proofs were faulty though his conclusions
were almost always right, and called him "a pioneer, not an educator".

Kron gave as hard as he got. His 1962 paper accuses graph theorists of
"sleight-of-hand" and heads one section "Transforming graphs into instruments of
torture". On this point he lost: graph-based network theory became the standard
treatment.

Support came from Alger, Hoffmann, Paul Langevin (who, Hoffmann says, saw to
the Montefiore Prize), and groups in Britain and Japan: the Tensor Club of Great
Britain, its *Matrix and Tensor Quarterly*, and the Research Association of
Applied Geometry in Tokyo.

## Legacy

**What became standard.**

- Computers: his obituary credits him with the "logical and procedural
  foundations" for programming a computer to compile and solve network
  equations automatically.
- Power systems: Kirchmayer's economic dispatch books (1958, 1959) rest on his
  work, as does Happ's *Diakoptics and Networks* (1971). Kron reduction, the
  Schur complement of a network's admittance matrix, still carries his name
  (Dörfler and Bullo).
- Numerical methods: tearing has since reappeared as domain decomposition in
  parallel computing (Lai 1994), and the transmission-line matrix method cites
  his 1944 paper (Russer).
- Electromagnetic compatibility: tensorial analysis of networks is still in use
  (Leman et al. 2009; Maurice et al. 2014).

**What did not.** Universal engineering never became an organized discipline
under that name, and the polyhedral program was carried on only by a few. In
1973 Edwin Tolan knew of no one who had taken up the work where Kron left off,
but in Britain some had: N. V. Balasubramanian's Liverpool thesis (1968, under
J. W. Lynn) treated Kron's multidimensional wave automaton with differential
forms; H. Nicholson at Sheffield wrote two reports on the polyhedron model and
the scattering problem (1973, 1974); and Lynn and Russell published "Kron's Wave
Automaton" in 1974. Tolan suspected that Kron was working toward a general
theory of the universe as a whole system, but he offered that as his own
belief, and it should be read that way. Later work on Kron's polyhedra and dual
networks exists (Popkov 2007; Petrov 2008), but we have seen it only in
excerpts.

## In kronspice

| Kron's idea | Where |
|---|---|
| Primitive network and connection tensors | [kron/primitive.go](kron/primitive.go), [kron/connection.go](kron/connection.go) |
| Open and closed paths (orthogonal network) | [kron/orthogonal.go](kron/orthogonal.go) |
| Tearing | [kron/diakoptics.go](kron/diakoptics.go), [kron/happ.go](kron/happ.go) |
| Generalized machine | [kron/machine.go](kron/machine.go) |
| Non-electrical domains as networks | [doc/domains.md](doc/domains.md) |
| Schrödinger equation as a network | [examples/kron_schrodinger/main.go](examples/kron_schrodinger/main.go) |

kronspice does not yet implement Kron's polyhedra, p-networks or crystal
computer. The physical reading of the 1962 polyhedral waves is register entry
H20 in [CLAUDE.md](CLAUDE.md#hypotheses-under-test), with a build plan for
`kron/polyhedron.go`.

## Where the sources disagree

| Item | Used here | Other versions |
|---|---|---|
| Death | March 25, 1968 (IEEE obituary, Alger, Happ) | October 25, 1968 (Wikipedia). The obituary appeared in September 1968, so October cannot be right |
| Birth | July 23, 1901 (Alger) | December 1, 1901 (Wikipedia) |
| Arrival in the US and entry to Michigan | January 1921, September 1921 (Alger; the obituary also gives 1921) | December 1920 departure, fall 1922 entry, 1925 graduation (Wikipedia) |
| Non-Riemannian paper written | During the year in Romania, about 1932–33 (Alger) | 1931 (Maurice et al.) |
| Montefiore Prize | 1935 (obituary, Alger, Happ) | 1933 (Maurice et al.) |
| Paper's journal volume | *J. Math. Phys.* 13 (Kron's publication list; Kron 1962) | 12 (Roth 1959) |
| *A Short Course* | 1942 (publication list) | 1940 (Alger) |
| Michigan honor, 1936 | Honorary M.S. (obituary, Happ) | "Doctor of Engineering" (Anders website) |
| End of GE career | Retired July 1966 (obituary, Happ) | Worked there until his death (AIP catalog) |
| *Diakoptics* serial | 20 chapters, June 7, 1957 – February 13, 1959 (Kron 1962; publication list) | 21 installments (Wikipedia). The GE report edition has 21 chapters, which may be the source of the confusion |
| "Camouflaging" paper | *Quarterly of Applied Mathematics* 20 (the paper's own page headers) | Vol. 10 (publication list) |
| Number of papers | More than 100 (Alger) | More than 50 (obituary) |

We left out claims online that link Kron to free energy from "negative
resistors", wartime secret projects or similar. We found them only in
speculative sources, and they are not supported by the sources listed here.

## Sources

### Read for this file

Only the parts named were read, unless noted otherwise.

- G. Kron, *Tensor Analysis of Networks*, Wiley, 1939; Macdonald, 1965, with a
  new introduction: Ch. VIII §1 and the bibliography of the 1939 text. In
  `reference/`.
- G. Kron, *Tensors for Circuits*, Dover, 1959: the introduction by B. Hoffmann,
  dated June 1958, pp. vii–xi. In `reference/`.
- "Gabriel Kron, 1901–1968," *IEEE Transactions on Circuit Theory* 15(3), 174,
  September 1968. [Copy](https://punainen.org/~biotek/kron/kron.pdf).
- P. L. Alger, "The Evolution of an Engineering Scientist," pp. 5–17; H. H.
  Happ, Introduction; E. K. Tolan, Foreword; and "Publications of Gabriel Kron
  1930–1968," pp. 165–172, all in H. H. Happ (ed.), *Gabriel Kron and Systems
  Theory*, Union College Press, 1973. Read as transcribed on
  [U. Anders's website](http://www.quantum-chemistry-history.com/Kron_Dat/KronGabriel1.htm),
  not from the printed book. Page numbers come from the book's contents list.
  Publication details in this file come from that list unless noted.
- G. Kron, "Camouflaging electrical 1-networks as graphs," *Quarterly of Applied
  Mathematics* 20(2), 161–174, 1962: read in full, as extracted text (equations
  and figures are garbled). Its bibliography misprints *Tensor Analysis of
  Networks* as 1931. In `reference/`.
- N. V. Balasubramanian, *Differential structures in electromagnetic field and
  circuit theory*, PhD thesis, University of Liverpool, July 1968, supervised
  by J. W. Lynn: Chapter 1 (introduction, pp. 6–22, with its long quotations
  from Kron's 1958–65 papers), the opening of Chapter 2 and the reference
  list, as extracted text. The rest (exterior forms, Kron's Maxwell network in
  forms, waveguides, machine frames) is unread.
  [Open copy](https://livrepository.liverpool.ac.uk/3173930/1/642931.pdf); in
  `reference/` (`642931.pdf`, 361 pages, with a text layer).
- G. Kron, "Multidimensional curve-fitting with self-organizing automata,"
  *Journal of Mathematical Analysis and Applications* 5, 46–69, 1962 (submitted
  by Richard Bellman): abstract, introduction and closing summary, as extracted
  text. In `reference/` (`PII_ 0022-247X(62)90005-7.pdf`).
- H. H. Happ, *Diakoptics and Networks*, Academic Press, 1971: preface and
  acknowledgments. In `reference/`.
- J. P. Roth, "An application of algebraic topology: Kron's method of tearing,"
  *Quarterly of Applied Mathematics* 17, 1–24, 1959: introduction and references.
  In `reference/`.
- O. Maurice et al., "Kron's method and cell complexes for magnetomotive and
  electromotive forces," 2014, [arXiv:1412.0187](https://arxiv.org/abs/1412.0187):
  history section. In `reference/`.
- S. Leman et al., "Use of the circuit approach to solve large EMC problems,"
  *C. R. Physique* 10, 70–82, 2009: abstract and introduction. In `reference/`.
- P. Russer, "The Transmission Line Matrix Method," TU München: reference list
  only. In `reference/`.

### Known only from citations, catalogs or excerpts

We have not read these. They should be checked before being relied on, and the
first group is where to start.

**Most useful to obtain:**

- G. Kron, the two 1956 articles listed under
  [What "universal engineering" meant](#what-universal-engineering-meant).
  Their only citation we found is the publication list.
- G. Kron, "Tearing, tensors and topological models," *American Scientist* 45,
  401–413, 1957.
- G. Kron, "The frustrating search for a geometrical model of electrodynamic
  networks," *Tensor* (n.s.) 13, 111–128, 1963.
- G. Kron, the 1965 introduction to *Tensor Analysis of Networks*. Tolan points
  to it as the best evidence for Kron's late views.
- P. L. Alger (ed.), *The Life and Times of Gabriel Kron*, Mohawk Development
  Service, Schenectady, 1969. Contains Kron's autobiographical writing.
- The rest of *Gabriel Kron and Systems Theory*: B. Hoffmann, "On G. Kron's
  methods and achievements"; H. H. Happ, "The development of diakoptics,"
  pp. 83–120; A. Brameller and D. W. Mortifee, "The influence of G. Kron in the
  United Kingdom," pp. 143–151; K. Kondo, "An Oriental expansion of Kron's
  science beyond electrical engineering," pp. 153–164.

**Kron's other papers:**

- "Generalized theory of electrical machinery," *AIEE Transactions* 49,
  666–683, 1930.
- "Non-Riemannian dynamics of rotating electrical machinery," *J. Math. Phys.*
  13, 103–194, 1934.
- "Quasi-holonomic dynamical systems," *Physics* 7, 143–152, 1936; and
  "Invariant form of the Maxwell–Lorentz field equations for accelerated
  systems," *J. Appl. Phys.* 9, 196–208, 1938 (Balasubramanian refs. 40–41;
  the second bears on §3l of CLAUDE.md).
- "Equivalent circuits of compressible and incompressible fluid flow fields,"
  *J. Aeronautical Sciences* 12, 221–231, 1945 (Balasubramanian ref. 55).
- "Equivalent circuit of the field equations of Maxwell — I," *Proc. IRE* 32,
  289–299, 1944.
- "A set of principles to interconnect the solutions of physical systems,"
  *J. Appl. Phys.* 24, 965–980, 1953.
- "A method of solving very large physical systems in easy stages," *Proc. IRE*
  42, 680–686, 1954.
- "Electric circuit models of the nuclear reactor," *AIEE Transactions* 73,
  259–265, 1954.
- "Inverting a 256 x 256 matrix," *Engineering* 178, 309–312, listed under 1955.
- "Solving highly complex elastic structures in easy stages," *ASME Journal of
  Applied Mechanics* 77, 235–244, 1955.
- "Tearing and interconnecting as a form of transformation," *Quarterly of
  Applied Mathematics* 13, 147–159, 1955. In `reference/`
  (`S0033-569X-1955-73431-7.pdf`); not yet read for this file. Its reference
  list cites N. Coburn, "A new approach to Kron's work," *J. Math. Phys.* 17,
  112–122, 1938, and a paper by Hermann Weyl on current distribution in a
  conducting network (*Revista Matemática Hispano-Americana*).
- "Solution of complex nonlinear plastic structures by the method of tearing,"
  *J. Aeronautical Sciences* 23, 1956.
- "Multiple substitution of basic vectors in linear programming," *Matrix and
  Tensor Quarterly* 7, 1956.
- "A generalization of the calculus of finite differences to nonuniformly
  spaced variables," *AIEE Transactions* Part I, 77, 539–544, 1958; and
  "Multidimensional space filters," *Matrix and Tensor Quarterly* 9(2), 1958.
  Both cited in the 1962 *J. Math. Anal. Appl.* paper as the groundwork of the
  curve-fitting polyhedron; Balasubramanian gives the first as vol. IX,
  pp. 40–46, December 1958, and quotes both (pp. 8–11).
- "Building-blocks of self-organizing automata – I," General Electric Report
  60GL164, August 1960, 24 pp. Quoted by Balasubramanian (pp. 15–17); the
  "dynamo-type automata" paper in *Matrix and Tensor Quarterly* 11, 1960, is
  presumably its published form.
- "Basic concepts of multi-dimensional space filters," *AIEE Transactions* 78,
  554–561, 1959.
- "Graphs as illegitimate models of electrical networks," *Matrix and Tensor
  Quarterly* 12, 1–10, 1961.
- "The misapplication of graph theory to electrical networks," *AIEE
  Transactions* 81, 257–267, 1962.
- "Invisible dual (n-1)-networks induced by electric 1-networks," *IEEE
  Transactions on Circuit Theory* 12(4), 464–470, 1965.

Copies of Kron's 1945 papers, "Electric circuit models of the Schrödinger
equation" and "Numerical solution of ordinary and partial differential equations
by means of equivalent circuits," and of *Equivalent Circuits of Electric
Machinery* are in `reference/`. We did not use them for this history beyond
their citations.

**Others:**

- B. Hoffmann, "Kron's method of subspaces," *Quarterly of Applied Mathematics*
  2, 218–231, 1944; "Kron's non-Riemannian electrodynamics," *Reviews of Modern
  Physics* 21(3), 1949.
- J. P. Roth, "An application of algebraic topology to numerical analysis: on
  the existence of a solution to the network problem," *PNAS* 41(7), 518–521,
  1955. In `reference/` (`roth-1955-an-application-of-algebraic-topology-to-
  numerical-analysis-on-the-existence-of-a-solution-to-the-network.pdf`); not
  yet read for this file.
- H. Nicholson, "The Structure of Kron's Polyhedron Model and the Scattering
  Problem," ACSE Research Report 14, University of Sheffield, 1973
  ([open copy](https://eprints.whiterose.ac.uk/id/eprint/86081/)); journal
  version *Proc. IEE* 120(12), 1545–1557, 1973.
- H. Nicholson, "System Concepts in Kron's Polyhedron Model and the Scattering
  Problem," ACSE Research Report 26, University of Sheffield, 1974
  ([open copy](https://eprints.whiterose.ac.uk/id/eprint/86082/)); journal
  version *Int. J. Control* 20(4), 529–555, 1974. In `reference/`
  (`acse research report 26.pdf`, a scan without a text layer); not yet read.
- J. W. Lynn and R. A. Russell, "Kron's Wave Automaton," in J. J. van Dixhoorn
  and F. J. Evans (eds.), *Physical Structure in Systems Theory*, Academic
  Press, 1974, p. 131.
- W. Givens, review of *Tensor Analysis of Networks*, *Bulletin of the AMS* 47(7),
  536–538, 1941.
- L. K. Kirchmayer, *Economic Operation of Power Systems*, Wiley, 1958; *Economic
  Control of Interconnected Systems*, Wiley, 1959.
- K. Bowden, "Physical computation and parallelism," *International Journal of
  General Systems* 27(1–3), 93–103, 1998.
- C. H. Lai, "Diakoptics, domain decomposition and parallel computing," *The
  Computer Journal* 37(10), 840–846, 1994.
- M. C. Duffy, §9.12 in I. Grattan-Guinness (ed.), *Companion Encyclopedia of the
  History and Philosophy of the Mathematical Sciences*, vol. 2, Routledge, 1994,
  pp. 1228–1234.
- F. Dörfler, F. Bullo, "Kron reduction of graphs with applications to electrical
  networks," [arXiv:1102.2950](https://arxiv.org/abs/1102.2950).
- V. Popkov, "Gabriel Kron's universal engineering,"
  [working paper](https://www.researchgate.net/publication/313407243_GABRIEL_KRON'S_UNIVERSAL_ENGINEERING),
  2007. Not peer reviewed as far as we can tell; search excerpts only.
- Petrov, "Tensor method and dual networks in electrical engineering," *Russian
  Electrical Engineering*, 2008,
  [doi:10.3103/S1068371208120018](https://doi.org/10.3103/S1068371208120018).
  Abstract excerpt only.

### Reference pages

- [Gabriel Kron](https://en.wikipedia.org/wiki/Gabriel_Kron) and
  [Diakoptics](https://en.wikipedia.org/wiki/Diakoptics), Wikipedia.
- [AIP International Catalog of Sources](https://history.aip.org/catalog/icos/44175.html):
  Gabriel Kron papers, 1963–1969.
- [Penn State catalog record](https://catalog.libraries.psu.edu/catalog/714021) for
  *Gabriel Kron and Systems Theory*.
- [Online Books Page](https://onlinebooks.library.upenn.edu/webbin/book/lookupid?key=olbp55036):
  the General Electric report edition of *Diakoptics*, 1957–58, with page images at
  HathiTrust.
