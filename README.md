# SMARTS is like a regular expression for chemical structures

SMARTS is a tool for structure searching, the process of finding a particular pattern
(a *subgraph*) in a molecule (a *graph*). Sub- (and super-) structure searches are
used in virtually every computer-based chemistry application.

## Defining the terms

* SMILES = __S__ implified __M__ olecular __I__ nput __L__ ine __E__ diting __S__ pecification

* SMARTS = __SM__ ILES __Ar__ bitrary __T__ arget __S__ pecification

A SMILES defines a specific chemical compound.  A SMARTS is an expression that looks 
for a match of a particular arrangement of atoms/bonds in that compound.

Inspired by xxx

## Scaffold screens

A good example of the utility of SMARTS is shown by the scaffold screening code in this 
project ('screen.for')

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
