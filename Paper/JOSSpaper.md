---
title: 'STEDY: Software for TEnsegrity DYnamics'

tags:
  - tensegrity
  - Lagrangian
  - dynamics
  - MATLAB

authors:
  - name: Vaishnav Tadiparthi
    affiliation: 1

  - name: Shao-Chen Hsu
    affiliation: 1

  - name: Raktim Bhattacharya
    affiliation: 2
affiliations:
 - name: Graduate Research Assistant, Texas A&M University
   index: 1

 - name: Associate Professor, Texas A&M University
   index: 2

date: 13 September 2018

---

# Summary

A tensegrity system is an arrangement of axially-loaded elements (no element bends, even though the overall structure bends), that we loosely characterize as a network of bars and cables. The bars take compressive axial loads and the cables handle tensile loads. Since failure due to axial stresses happens at higher loads than at bending, a tensegrity structure has a higher strength-to-weight ratio. The famous architect Buckminster Fuller, in the 60's, coined the term tensegrity, combining the words tensile and integrity. Since then, tensegrity principles have found applications in diverse domains like architecture, biological modeling as well as civil engineering design. Tensegrity structures, through use of pre-stresses in the bars and cables, can also achieve controlled stiffness in the structure, which makes it attractive in applications such as soft-robotics and prosthetics. In essence, tensegrity principles can be applied in the design of any structure where mass is premium, a high strength-to-weight ratio is critical, and structural stiffness needs to be tailored in both space and time. These include several applications from various engineering sectors such as aerospace (morphing airframes), energy (wind turbine blades, off-shore structures) as well as biomedical engineering (stents, minimally invasive surgical tools) and many more. Clearly, a framework is required that can efficiently model the dynamics of tensegrity structures directly from the topology of bars and cables.

The dynamics of tensegrity systems is governed by multi-body dynamics, given by a set of ordinary differential equations. We have developed a Lagrangian formulation for deriving these differential equations directly from the given topology of members (bars and strings), and their mass and geometric properties. Three key features of classical tensegrity systems are a) actuations only occur via cables, b) bar-to-bar connections are pin joints, and c) the bars do not spin about their longitudinal axis. These properties are exploited to simplify the equations of motion. However, the Lagrangian framework presented here is general enough to allow modeling of general multi-body systems with actuated joints.

`STEDY` is a MATLAB package for conducting numerically accurate tensegrity dynamics. The scientific computing community's familiarity with MATLAB makes it a favored choice for performing simulations of this nature. The equations of motion are developed in Cartesian coordinates. Although the choice of non-minimum coordinates is fairly controversial in the robotics community on account of their configuration space having been reduced by joints, adopting Cartesian coordinates precludes any singularities developed in the mass matrix and facilitates derivations of elegant differential-algebraic equations governing the multibody dynamics. However, usage of non-minimum coordinates is likely to result in the solution drifting away from the constraint space because of integration errors. The software package includes an implementation of the direct correction scheme that minimises numerically induced constraint violations and energy conservation errors. We have demonstrated superiority of the method in terms of accuracy over the commercially available dynamics simulator, Simscape Multibody.

`STEDY` was designed to be used by researchers looking to simulate tensegrity dynamics without having to derive the differential equations for any structure in particular. The inputs that are required from the user are simply the parameters that uniquely define the tensegrity structure (initial nodal configuration, connectivity), the material and geometric properties, and simulation environment properties (inertially fixed nodes, external forces, duration, ODE solver options).

The software documentation and the paper cited below (contact us for a copy) comprehensively cover the functionality and the theory behind the Lagrangian formulation and the novel constraint correction method implemented in the package. We shall seek to incorporate techniques introducing control of the dynamic tensegrity system in future versions of the software.

# Acknowledgements

This work was supported by NSF IUSE/PFE: RED: REvolutionizing Diversity Of Engineering (REDO-E)Award Number:1730693; and NASA NIAC Phase II grant, on Tensegrity Approaches to In-Space Construction of a 1g Growable Habitat.

# References
1. Shao-Chen Hsu, Vaishnav Tadiparthi and Raktim Bhattacharya, "A Lagrangian Formulation for Constrained Multibody Dynamics in Tensegrity Systems", Manuscript submitted for publication.

    Please contact the authors at addyhsu@tamu.edu or vaishnavtv@tamu.edu for a copy of the submitted paper.
