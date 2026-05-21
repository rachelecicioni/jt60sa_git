# Integrated modelling for VUV diagnostic for JT-60SA

**Author:** R. Cicioni  (rachele.cicioni@igi.cnr.it)  
**Last update:** 2026-05-21  

---

## Info

This document provides the necessary context to independently use the codes developed to reproduce line-of-sight VUV signals within the VUV project for JT-60SA, coordinated by L. Carraro.
Both this document and the associated Git repository aim to clearly document the workflow and track the development steps in an organized and reproducible way, ensuring that the work can be easily understood, reused, and extended in the future.

## Context

For the JT-60SA tokamak in Japan, a future transition to tungsten (W) wall components is planned. This project aims to provide density and temperature profiles along arbitrary lines of sight for the VUV diagnostic system, for both argon (Ar) and tungsten (W) impurities. The synthetic data have been generated using integrated modelling with the JINTRAC suite of codes, specifically the COCONUT module. This framework simulates the plasma in coupled 2D maps covering both the core and the scrape-off layer (SOL), allowing for a self-consistent description of the entire plasma region.  
The modelling work was supervised by L. Garzotti in his role as scientific coordinator of the Transport and Confinement (TC) group.

---

## Index

1. [JINTRAC simulations](#1-jintrac-simulations)  
   1.1 [What is JINTRAC](#11-what-is-jintrac)  
   1.2 [Impurities modelling in JINTRAC](#12-impurities-modelling-in-jintrac)  
   1.3 [Simulations list](#2-simulations-list)  
2. [VUV LOS code](#3-vuv-los-code)  
   2.1 [How to use it](#31-how-to-use-it)  
   2.2 [IDL compatibility](#32-idl-compatibility)  
3. [Bibliography](#4-bibliography)

---

# 1. JINTRAC simulations

## 1.1 What is JINTRAC? - Briefly -

The JINTRAC [1] framework allows the simultaneous simulation of core, edge, and SOL regions, accounting for their physical interconnection in a self-consistent manner. At the core of JINTRAC is the COCONUT [2] coupling scheme, which couples the 1.5D core transport code JETTO [3] with the 2D SOL transport code EDGE2D [4]. In the simulations used in this collaboration, an integrated core+edge+SOL modelling framework (COCONUT) was used to evolve the plasma self-consistently across the entire domain (core and SOL).

---

## 1.2 Impurities modelling in JINTRAC

Impurities in JINTRAC are modelled using the SANCO impurity transport solver, which evolves the density of each ionization stage of a given element. In full detail, this would require solving up to Z+1 coupled equations per impurity species. To reduce computational cost in integrated modelling, a charge-state bundling (***superstage***) approach is often used. In this approach, multiple adjacent ionization stages are grouped into a single “superstage". The number of evolved variables is therefore reduced to: N_superstages ≪ Z. Within each superstage, ionization states are assumed to be in internal ***coronal equilibrium***, meaning that fast atomic processes between bundled charge states are not resolved explicitly [5].  
Example: for W (tungsten) it was used the 6-stages boundling. Instead of solving 74 ionization stages, they are grouped in 6 superstages only.  

***ADAS database coupling***  
Atomic data are taken from the ADAS database provided in JDC server:  
`/usr/local/depot/sim/data/adas/adas_v1.2/adf11`  
JINTRAC uses a hybrid ADAS-year structure, meaning that different physical processes may come from different ADAS “years” to ensure both accuracy and numerical stability. The processes considered are:
- ACD: Coll.-diel. effective recombination coefficients  
- SCD: Coll.-diel. effective ionisation coefficients  
- CCD: Coll.-rad. effective charge exchange coefficients  
- PLT: Coll.-rad. line power driven by excitation of dominant ions  
- PRB: Coll.-diel. continuum and line power driven by recombination and Bremsstrahlung of dominant ions  
- PRC: Coll.-rad. effective charge exchange recombination power loss rate coefficients  
- PLS: Coll.-rad. line power from selected transitions of dominant ions (obsolete?)  
- ZCD: Superstate averaged Z  
- YCD: Superstate averaged Z^2  
- ECD: Superstate ionisation potential (eV)  

Therefore, “ADAS YEAR” in JINTRAC is not a single consistent dataset, but a composed model.  
The parameter "ADAS LABEL" defines the bundling scheme applied to the impurity charge states. Bundling files are located with the coefficients for each processes, year, label are in:
`adf11/plt<year>/plt<year>_<element>02<label>.dat`  
Example: plt42_w_02_6.dat cointaints the coefficients for W for the PLT process, year 42, label '02_6' which corresponds to 6 superstages.  
The correspondence between the label and the number of superstages can be found in this page of JINTRAC wiki: https://users.euro-fusion.org/pages/data-cmg/wiki/ADAS_database_for_EDGE2D.html


---

## 1.3 Simulations list

List of simulations with full paths in JDC (JET Data Center) cluster:

- COCONUT, Ar+W case: /home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/dec0525/seq#1 

List of username and reference person in JDC cluster:
- jnv7243: R. Cicioni
- gm7685: S. Gabriellini
- lgarzot: L. Garzotti

List of username and reference person in EFGW (EUROfusion Gateway) cluster:
- g2rcicio: R. Cicioni
---

# 2. VUV LOS code

Python script used for VUV line-of-sight reconstruction:
`vuv_los_coconut.py`

Purpose:
- reconstruction of VUV LOS signals
- post-processing of simulation outputs

---

## 2.1 How to use it

Instructions for running the script and required inputs.

Example:
```bash
python vuv_los_coconut.py input_file
```

---

## 2.2 IDL compatibility

Notes on compatibility with existing IDL workflows and possible limitations.

---

# 3. Bibliography and wikipages

[1] JINTRAC Wiki: https://users.euro-fusion.org/pages/data-cmg/wiki/index.html    
[2] Taroni A et al. 1996 Proceedings of the 16th iaea international conference on fusion energy Proc. 16th IAEA Int. Conf. on Fusion Energy vol 2 (Montreal, Canada: IAEA) p 477  
[3] Cenacchi G and Taroni A 1988 Jetto a free boundary plasma transport code Tech. rep. ENEA   
[4] Simonini R, Corrigan G, Radford G, Spence J and Taroni A 1994 Contributions to Plasma Physics 34 368–373   
[5] Bundled description of heavy impurities in JETTO/SANCO and EDGE2D, Laura Lauro-Taroni    
https://users.euro-fusion.org/pages/data-cmg/wiki/files/JETTO_bundledimpurities.pdf 
