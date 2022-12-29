# Description

This programm is a gillespie algorithm that simulate the aging/equilibration of an associative polymer gel.

## The System

the system is a polymer chain that is anchored at one extremity, the other being dangling. Additional binding occurs as the simulation progress. A state of the system is described by a list of sub-polymer chain that represent the free strands in between two binding points and called a Loop. the state of the system can be denoted:
$$ \mathcal{C} = \{ \mathbf{R}_i,\{ \mathbf{r}_n\}_{n \in [0,k_i]},\ell_i \}_{i \in [0,N]}$$
Where $\mathbf{R}_i$ is a 3D vector that represent the position of the $i$ th anchoring point. $k_i$ is the number of binding sites in the vicinity of the loop $i$, and $\mathbf{r}_n$ their position. $\ell_i$ is the length of the loop $i$. The final dangling part of the polymer is also fully determined with this notation.

## Time evolution

At each step, a list of possible move is built (see the moves section), then a move is picked with a probability proportional to its rate :
$$P_i = \frac{r_i}{\sum_j r_j}$$
where $P_i$ is the probability of picking the rate $i$ and $r_i$ its rate.

## Parameters of the simulation
- $\rho_0$ is either the volume fraction of crosslinkers in the system, or the initial volume fraction when the system evolve.
- $\ell_{tot}$ the length of the total loop.
- $a$ is the size of the crosslinkers, $b$ is the size of the monomers. We assume that they are equal and define the unit length so far.
- $\tau_0$ is the characteristic time for the rates. It defines the unit time of the simulation.
- $E_b$ is the binding energy, $>0$. It defines the energy scale.
- $\beta = 1/(k_B T)$ is the inverse temperature.

<!---
## initialization

We start with a single polymer chain bound at its two extremities. The initialization consist in subsequent binding of the chain until every fragment had a chance to bind. Here is the loop :
> - i=0
>    - Select the chain i  
>    - Compute the probability of the chain i to have a binding point according to the formula :  
>   - $$P_\text{binding} = V_\text{free} \rho_\text{stickers}$$  
>   - with $V_\text{free}$ the free volume occupied by the ellipse of the chain, and $\rho_\text{stickers}$ the density of sticker in the system (which is a free parameter).  
>   - draw a random number : RAND between 0 and 1.
>       - if RAND < $P_\text{binding}$ then the chain is
 splitted into two sub chains  
>       - if RAND > $P_\text{binding}$ the process stops for this sub-chain and the initialization continue for the next chain.

> - repeat the process for the two sub-chains, and the next chains  

The initialization stops when all the chains / sub-chains etc... answer False to any additional binding.
-->

## Moves

We do not simulate the polymer dynamic, one of the main assumption is that the polymer dynamic occurs on smaller timescale than the one of binding / unbinding dynamic. Therefore, the only time consuming moves are the binding/unbinding ones, that are thus assumed to dominate the dynamic of the system.

We start by creating a list of possible moves.

<ins>

### Binding Moves 
</ins>

a binding move occurs in two steps:
1. select a loop
2. select a crosslinker and  a length for the loop.

the binding rate of loop is:
$$r_b(\ell_i,\mathbf{R}_i,\mathbf{R}_{i+1}) = \int \frac{\text{d}\ell}{b}\frac{1}{\tau_0}\sum_{\mathbf{r}_n} \Omega(\mathbf{R}_i,\mathbf{r}_n,\ell)\Omega(\mathbf{r}_n,\mathbf{R}_{i+1},\ell_i-\ell)$$

We then choose a binding event of the previously selected loop, according to the following rate :
$$r_b(\ell_i,\mathbf{R}_i,\mathbf{R}_{i+1},\mathbf{r}_n,\ell) = \frac{1}{\tau_0} \Omega(\mathbf{R}_i,\mathbf{r}_n,\ell)\Omega(\mathbf{r}_n,\mathbf{R}_{i+1},\ell_i-\ell)$$
With
$$\Omega(\mathbf{R}_i,\mathbf{r}_n,\ell) = (4\pi)^{\frac{\ell}{a}} \left(\frac{3}{2\pi\ell_i a} \right)^{3/2} \exp \left[-\frac{3}{2}\frac{(\mathbf{R}_i-\mathbf{r}_n)^2}{\ell_i a} \right]$$
We compute this binding rate for every portion of monomer that compose the loop in order to select a linker and a length according to the probability defined previously.

Notice that the binding of a loop to a crosslinker always affect two loops, thus the total rate of creating a bound in $\mathbf{r}$ and length $\ell_i$ in a loop bound in $\mathbf{R}_0$ and $\mathbf{R}_1$ at its extremities and total length $\ell$ reads :
$$\tau_0 r_b(\mathbf{R}_0,\mathbf{R}_1,\ell,\mathbf{r},\ell_i) = \left(\frac{3 \ell}{2 \pi \ell_i(\ell-\ell_i)} \right)^{3/2} \exp\left[-\frac{3}{2}\left( \frac{(\mathbf{R}_0-\mathbf{r})^2}{\ell_i}+\frac{(\mathbf{R}_1-\mathbf{r})^2}{\ell-\ell_i} - \frac{\mathbf{R}_0-\mathbf{R}_1)^2}{\ell}\right) \right]$$

for the dangling end of the polymer, it reads :
$$\tau_0 r(\mathbf{R}_0,\ell,\mathbf{r},\ell_i) = \left( \frac{3}{2\pi \ell_i}\right)^{3/2}\exp\left[ -\frac{3}{2}\frac{(\mathbf{R}_0-\mathbf{r})^2}{\ell_i}\right]$$

<ins>

### Unbinding moves
</ins>

Every anchoring point can unbind according to the following rate :
$$r_{ub} = \frac{e^{-\beta E_b}}{\tau_0}$$

### Slide linkers

Each linker can slide along the polymer by a length a (=1). Thus cum_rates has a length : loop.size()+(bounded_linker.size()-1) (1 is the (0,0,0) linker). The sliding rate is computed by each loop, that propose to slide its right linker by a step $\pm 1$. 

## Linkers management

A linker is a proper object, its main parameter is its position : three double x,y,z. It also owns a reference to every strand that have this linker in its vicinity. Every time a loop is modified a single linker changes its state from bound to unbound. This implies a change in the rate of the processes that concern this linker. To find all the strand that are in the vicinity of the concerned linker, we ask this linker to provide a reference to every strand in its vicinity. Based on the new status of the concerned linker, these strands are thus remade, and the rate of their associated processes are recomputed.

# Python and C++ objects
:warning:Not updated !:warning:
## System
The private variable of the system is:
- A **set<loop*>** of all the *loop* of the system. Think about allocating the memory we call it *loops*. plus 


The function of the system is:

- *evolve* function that
  - recompute *cum_rates* by iterating over *loops* and accessing the associated total rate.
  - A **vector<double>** of cumulative rates. the $i$ element is just the sum for $n$ from 0 to $i$ of *loop[n].get_rate()*. We call it *cum_rates*.
  - draw a random number *rand* and select the corresponding index $i$ so that:
  $$ cum\_rates[i]<rand<cum\_rates[i+1]$$
  - the corresponding loop is *loops[i]*
  - This loop is then asked for selecting a given crosslinker at a certain length, and returns it using *select_link_length*.
  - Remove the concerned *loop* from *loops*.
  - **delete** the *loop** pointer.
  - create the **new** *loop**s.
  - add the newly created *loop**s to *loops*.
  - returns a time increment.
- *output*${r_n}_{n \in [0,k]} : an array<array<3>,k> of the position of all binding points called r

  - sort the loop by their value of $\mathbf{R}$.
  - return three vectors that corresponds to $\mathcal{C}$.

## Loop

The main part of the system is a **vector<loop>** of loop. Each loop $i$ contains the following private variables:
- $\mathbf{R}$ : is a **array<double,3>** the right crosslinker. is called *R*
- $\ell$ : the length of the loop is a double called *ell_loop*
- $\{r_n\}_{n \in [0,k]}$ : an **vector<array<3>,k>** of the position of all binding points called *r*

The following variables accessible by accessor:

- an **array<double,k*ell/a>** of the binding rates at any crosslinker in length $\ell$ called *binding_rates*
- an unbinding rates called *unbinding_rates* and is a **double**

The following functions:

- a function the compute the binding rate of the loop to a crosslinker located to $\mathbf{r}_i$ at the length $\ell$.
- a function that select a crosslinker and a length (according to their respective rates) and returns the position of the crosslinker and the length selected. we call it *select_link_length*
- The constructor:
  - take a right crosslinker position *R*, a length *ell*.
  - generate a **array<array<3>,k>** *r* by drawing a $k$ from the Poisson distribution

## Dangling

The end of the polymer. It is a single object owned by the System. Very similar to a loop. It always evolves following an alternative path

## Function

- a function that compute the square of a 3D vector.  

:warning: End of the not updated section :warning:

# Remarks

## Miscellaneous remarks

- the loops are sorted by curvilinear coordinates along the polymer length in the set
- The system know a reference to all the crosslinkers. However it is ask to each loop and dangling to create their own linkers, at a distance set in their own object.
- We arange the linkers in the map with : $key = x+y\ell=z\ell^2$

## *Strands*, *Dangling* and *Loops*

### crosslinkers :

crosslinkers are generated in a volume around each loops and dangling bond, but also stored in the system. each time a loop or a strand is generated, we draw $N_\text{link}$ linkers from a Poisson distribution of mean $\rho V$ where V is the volume in which the crosslinkers are generated. We do not generate crosslinkers in the whole volume which would be very large, while the system only occupies a negligible part of the total volume.

We generate crosslinkers in :
- a cube of side $R = 2 \sqrt{\ell}$ for the dangling bond
*warning*
- a cube of side $R = 2 \sqrt{\ell}$ for the loops with anchoring point $\mathbf{R}_0$ and $\mathbf{R}_1$ are close so that : $|\mathbf{R}_0-\mathbf{R}_1| < 0.1 \ell$.
- a rectangle of large side : $\ell$ and small axis $2\sqrt{\ell}$ and small side : $R2 \sqrt{\ell}$ for loops where : $|\mathbf{R}_0-\mathbf{R}_1| > 0.1 \ell$.
*warning*

- for now we will just generate crosslinkers in a cube of side $R= |\mathbf{R}_0-\mathbf{R}_1|+2\sqrt{\ell}$

### Strands

strands is the parent class of *Loop* and *Dangling* it basically encapsulate most of what loop and dangling do. Few specificity are handle by the child class. Mostly geometrical properties (related to the space they occupy). The following function are overwritten by the child class :
- get_volume_limit : get the slicing limits to generate crosslinkers in the vicinity of the strand.
- compute_rate : compute the binding rate of the strand.
- random_in_volume : generate coordinates to create crosslinkers in the surrounding of the strand.

*warning* the constructor of strands  does not call *generate_binding_sites* and *compute_all_rates* despite these function being exclusively owned by the *Strand* class. This is because other geometrical parameters need to be set before we call them. Because *Strand* constructor is called before *Loop* or *Dangling* one we have to call these function at the end of the constructor.

### LoopLinkWraps and loop_link
In the file Container.h we create a class of container that manage the interaction between loop and linkers. The goal of this container is to manage the following dependencies :

- Each strand own a list of crosslinker in its vicinity. When a strand bind somewhere, one of the crosslinker is not free anymore. This implies to change the rate of all the strand that has this crosslinker in its vicinity.
- LoopLinkWrap thus also own a map that make the link between every linkers and the multiple strands that own this linker.
- To quickly access the list of linkers in a vicinity of a strand, the LoopLinkWrap contain a map3D which is a custom container that can perform fast slicing. The vicinity is thus defined as cube around an area.

The LoopLinkWrap must thus manage all the interaction between the three containers :
- set<Strand*,LessLoop> loops : list of the strands.
- map3d free_linkers : container of the free linkers for quick slicing.
- map<array<double,3>,Strand*> linker_to_strand : link between a linker and a strand.

