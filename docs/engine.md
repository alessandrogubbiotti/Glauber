---
layout: default
title: The Gillespie Algorithm
---

[Back to the index](index.html)

# The Gillespie Algorithm

The dynamics is simulated event by event, in continuous time, with no time step. The engine is [`src/interface_state.c`](https://github.com/alessandrogubbiotti/Glauber/blob/main/src/interface_state.c), and the programs of [the programs page](programs.html) are built on it.

---

## 🔹 The State

The engine keeps the **interfaces**, not the spins. An interface is a bond $b$ of the ring such that the spins at $b$ and at $b+1$ are different. The interfaces are kept in an array sorted by position, each one with the spin on its left.

Between two consecutive interfaces there is a **domain**, and its **width** is the number of its sites. On the ring the number of interfaces is even.

---

## 🔹 The Events

In **code units** the hop rate is $1$ per direction.

- **Annihilation:** a domain of width $1$ has rate $\mathrm{ann}$. Its spin flips and its two interfaces disappear.

- **Movement:** a domain of width at least $2$ has rate $2$. One of its two extreme spins flips, the left one or the right one with probability $1/2$, and the interface on that side moves by one bond into the domain. So an interface hops with rate $1$ per direction, into empty bonds only.

- **Creation:** a spin equal to both its neighbours flips with rate $\mathrm{create}$, and two interfaces appear at its sides. A domain of width $w \geq 2$ has $w - 2$ such spins, and when there are no interfaces all the spins of the ring are such. The creation is one event for the whole ring, with rate $\mathrm{create}$ times the number of these spins, and the spin that flips is chosen uniformly among them.

Every event flips exactly one spin.

---

## 🔹 The Clocks

Each domain carries the time of its next event, exponential with parameter the rate of the domain. One more time is carried for the creation. At each step the engine

1. takes the smallest of these times, among the domains and the creation;
2. moves the clock to it and executes the event;
3. updates the times that the event has changed.

After a movement the domain that has fired draws a new time, and the neighbour that has grown draws a new one only if its rate has changed, that is if its width has passed from $1$ to $2$. After an annihilation the merged domain draws a new time. When the rate of the creation changes, the time that is left to it is rescaled by the ratio between the old rate and the new one. After a creation, and when an interface crosses the bond that closes the ring, the times of all the domains are drawn again.

Since the times are exponential, drawing them again or rescaling them does not change the law of the process.

The macroscopic time is $\tau = \text{micro time}/N^2$, see [the model](model.html).

---

## 🔹 No Interfaces

When the last two interfaces annihilate the ring has one colour, and the engine remembers which one. A creation from this state gives a domain of the opposite colour inside it.

---

## 🔹 The Cost

An event costs a scan of the times of the domains, so it is proportional to the number of interfaces. The number of events in a unit of macroscopic time is about $2N^2$ times the number of interfaces. With $N = 150$, $L = 4$ and $\nu = 1$, up to the macroscopic time $4$, a trajectory is about $725\,000$ events.

The branch `heap-engine` keeps the times in a heap. It gives the same trajectories and, in the cases we have timed, it is not faster.

---

## 🔹 The Functions

They are declared in [`src/interface_state.h`](https://github.com/alessandrogubbiotti/Glauber/blob/main/src/interface_state.h).

- `istate_alloc(spins, l, N, ann, create, seed)` builds the state from an array of `l` spins, $+1$ or $-1$. Each state has its own generator, seeded with `seed`.
- `istate_advance_to(s, t)` runs the dynamics up to the macroscopic time `t`.
- `istate_to_spins(s, spins)` writes the spins of the current configuration.
- `istate_n_interfaces(s)` is the number of interfaces.
- `istate_track_flips(s)` and `istate_never_flipped_fraction(s)`: the fraction of the sites whose spin has never flipped since the first call, that is the persistence.
- `istate_free(s)` frees the state.

```c
#include "interface_state.h"

int spins[600];                       /* the initial configuration, +1 or -1 */
/* ... fill spins ... */
InterfaceState *s = istate_alloc(spins, 600, 150, ann, create, 12345UL);
for (int t = 1; t <= 400; t++) {
    istate_advance_to(s, t / 100.0);  /* macroscopic times 0.01, 0.02, ... */
    istate_to_spins(s, spins);
}
istate_free(s);
```
