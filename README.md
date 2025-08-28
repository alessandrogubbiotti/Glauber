# Glauber Dynamics in the Zero-Temperature Continuum Limit

This repository explores the **zero-temperature continuum limit** of the **1-dimensional Glauber dynamics**.  
The model considers sites spaced by  $1/N$ with

![equation](https://latex.codecogs.com/svg.latex?\dpi{150}\color{white}N%20=%20\frac{e^{2\beta}}{1%20-%20e^{\beta}},) 

and time is **diffusively rescaled**.

---
## 🔹 Pages

For a rendered version of the project documentation, see the [📄 GitHub Pages site](https://alessandrogubbiotti.github.io/Glauber/).

## 🔹 Description of the Dynamics

The dynamics can be described in terms of **interfaces (particles):**

- **Single interface motion:**  
  Before taking the diffusive limit, an interface moves to the right or left with rate  
  ![equation](https://latex.codecogs.com/svg.latex?\dpi{150}\color{white}\large%201),  
  provided no other interface is present.

  ![diffusion image](docs/images/diffusion.png)

- **Annihilation of neighboring particles:**  
  Two neighboring particles are annihilated with rate  
  ![equation](https://latex.codecogs.com/svg.latex?\dpi{150}\color{white}\large%202%20-%20\frac{1}{N^a}).

  ![annihilation image](docs/images/annihilation.png)

- **Creation of particles:**  
  At each site, two particles are created with rate 
  ![equation](https://latex.codecogs.com/svg.latex?\dpi{150}\color{white}\large%20\frac{1}{N^a}).

  ![creation image](docs/images/creation.png)

---

## 🔹 Invariant Measure

The invariant measure is that of a **Poisson point process** of intensity one:

- The number of particles follows a **Poisson distribution**.  
- Particle positions are **uniformly distributed**.

Example simulation output:

![simulation image](docs/images/simulation_32.png)

---


