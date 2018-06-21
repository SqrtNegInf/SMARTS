## SMARTS is like a regular expression for chemical compounds

SMARTS is a tool for structure searching, the process of finding a particular pattern
(a *subgraph*) in a molecule (a *graph*). Structural searches are
used in virtually every computer-based chemistry application.

Initially inspired by UNIX implementation of regular expressions, 
SMARTS evolved into somethings that can be viewed as a 
domain-specific grammar for chemical compounds, combining 
the functions of traversing a graph and matching a pattern.
An overview of the syntax can be in the [SMARTS summary](./cheat-sheet.txt)

### Defining terms, an example

* SMILES = __S__ implified __M__ olecular __I__ nput __L__ ine __E__ diting __S__ pecification
* SMARTS = __SM__ ILES __Ar__ bitrary __T__ arget __S__ pecification

A SMILES defines a specific chemical compound, here the pesticide DDT:  

*  Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl

A SMARTS is an expression that looks for a match of a particular arrangement 
of atoms/bonds in that compound, here di-aryl ethane with any 3 halogens:

* C(c)(c)C([Br,Cl,F])([Br,Cl,F])([Br,Cl,F])

This matches DDT, but three other compounds:

* 2,2-DIPHENYL-1,1,1-TRICHLOROETHANE
* 1,1,1-TRIFLUORO-2,2-DI(P-METHOXYPHENYL)ETHANE
* METHOXYCHLOR

Here's an image of the [DDT search results](./ddt.png)

### Scaffold dictionary

Precisely targeted structural searches are often non-trivial, so it makes sense
to develop a standard library for re-use. Here is a compendium of nearly 
[700 SMARTS targets](./scaffold-dictionary.txt).

### Scaffold screens

Another good example of the utility of SMARTS is shown by the scaffold screening code in this 
project [screen.for](./screen.for)

The scaffold screens consist of a mixture of custom code and SMARTS targets. The trade-off
is that while SMARTS targets are relatively easy to write, they are slow at runtime. Custom
code is generally 10x or more faster, but significantly harder to write and maintain.

To reduce the runtime penalty for SMARTS, a quick analysis is done for the input compound,
tallying atoms, bonds, rings (sizes and counts), basic structural features, etc.  This is
used to avoid running a lengthy SMARTS search when it cannot possibly succeed (say the target
contains two nitrogens, but the compound has only one).

A good comparison between the two approaches can be seen in the current SMARTS-based
'phenothiazine' code (around line 500) and the previously-used custom code, left commented
just after.
