# Description

This programm is a gillespie algorithm that simulate the aging/equilibration of an associative polymer gel.

## The System

the system is a polymer chain that is anchored at its two extremities. Additional binding occurs as the simulation progress. A state of the system is described by a list of sub-polymer chain that represent the free strands in between two binding points. the state of the system can be denoted:
$$ \mathcal{C} = \{ \mathbf{R}_i,\{ \mathbf{r}_n\}_{n \in [0,k_i]},\ell_i \}_{i \in [0,N]}$$
Where $\mathbf{R}_i$ is a 3D vector that represent the position of the $i$th anchoring point. $k_i$ is the number of binding sites in the vicinity of the loop $i$, and $\mathbf{r}_n$ their position. $\ell_i$ is the length of the loop $i$.
## Time evolution

At each step, a list of possible move is built (see the moves section), then a move is picked with a probability proportional to its rate :
$$P_i = \frac{r_i}{\sum_j r_j}$$
where $P_i$ is the probability of picking the rate $i$ and $r_i$ its rate.

## Parameters of the simulation
- $\rho_0$ is either the volume fraction of crosslinkers in the system, or the initial volume fraction when the system evolve.
- $D$ the distance between the first, and final anchoring point.
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
>       - if RAND < $P_\text{binding}$ then the chain is splitted into two sub chains  
>       - if RAND > $P_\text{binding}$ the process stops for this sub-chain and the initialization continue for the next chain.

> - repeat the process for the two sub-chains, and the next chains  

The initialization stops when all the chains / sub-chains etc... answer False to any additional binding.
-->

## Moves

We do not simulate the polymer dynamic, one of the main assumption is that the polymer dynamic occurs on smaller timescale than the one of binding / unbinding dynamic. Therefore, the only time consuming moves are the binding/unbinding ones, that are thus assumed to dominate the dynamic of the system.

We start by creating a list of possible moves.

#### Binding Moves
a binding move occurs in two steps:
1. select a loop
2. select a crosslinker and  a length for the loop.

the binding rate of loop is:
$$r_b(\ell_i,\mathbf{R}_i,\mathbf{R}_{i+1}) = \int \frac{\text{d}\ell}{b}\frac{1}{\tau_0}\sum_{\mathbf{r}_n} \Omega(\mathbf{R}_i,\mathbf{r}_n,\ell)\Omega(\mathbf{r}_n,\mathbf{R}_{i+1},\ell_i-\ell)$$

We then choose a binding event of the previously selected loop, according to the following rate :
$$r_b(\ell_i,\mathbf{R}_i,\mathbf{R}_{i+1},\mathbf{r}_n,\ell) = \frac{1}{\tau_0} \Omega(\mathbf{R}_i,\mathbf{r}_n,\ell)\Omega(\mathbf{r}_n,\mathbf{R}_{i+1},\ell_i-\ell)$$
With
$$\Omega(\mathbf{R}_i,\mathbf{r}_n,\ell) = (4\pi)^{\frac{\ell}{a}} \left(\frac{3}{2\pi\ell_i a} \right)^{3/2} \exp \left[-\frac{3}{2}\frac{(\mathbf{R}_i-\mathbf{r}_n)^2}{\ell_i a} \right]$$
We compute this binding rate for every portion of monomer that compose the loop.

#### Unbinding moves
Every anchoring point can unbind according to the following rate :
$$r_{ub} = \frac{e^{-\beta E_b}}{\tau_0}$$

<!-- We choose a bond to unbind according to its unbinding rate :  
$$r_\text{unbind} = 1/\tau_0 e^{-\beta E_\text{bind}+\Delta S}$$
where $\Delta S$ is the difference of polymer entropy between the bound and unbound state, $E_\text{bind}$ is the binding energy, and $\beta$ the temperature.

#### Rebinding

We then re-do a similar re-binding chain of events similarly to the initialization step, but only with the bond affected by the unbinding event. Each rebinding event consist in drawing a random binding point coordinates : $\overrightarrow{r}$ and a random $\ell$ that corresponds of the length of the polymer after binding on the left-hand side.
-->
# Python and C++ objects

## System
The private variable of the system is:
- A **vector<loop*>** of all the *loop* of the system. Think about allocating the memory we call it *loops*.


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

## Function

- a function that compute the square of a 3D vector.  
