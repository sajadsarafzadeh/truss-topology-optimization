# Truss Topology Optimization (MATLAB)

Topology optimization of a 2D truss structure implemented **from scratch in MATLAB**, including finite element modeling, sensitivity analysis, and a gradient-based optimization algorithm (CONLIN).

---

## 📌 Overview

The objective of this project is to **minimize structural compliance (maximize stiffness)** under a **volume constraint** by optimizing the cross-sectional areas of truss elements.

The implementation builds the full pipeline:
- Finite Element Analysis (FEA)
- Assembly of the global stiffness matrix
- Sensitivity (gradient) computation
- Constrained optimization using CONLIN

---

## ⚙️ Features

- 2D Truss Finite Element Solver
- Custom stiffness matrix assembly (`buildStiffnessMatrix.m`)
- Structural response computation (`computeStructure.m`)
- Stress evaluation (`computeStress.m`)
- Volume calculation (`computeVolume.m`)
- Sensitivity analysis:
  - Objective gradient (`gradOf.m`)
  - Constraint gradient (`gradConstr.m`)
- Optimization loop (`main_training_005.m`)
- Multiple design scenarios with different volume fractions

---

## 🧮 Problem Formulation

**Objective:**
Minimize compliance

C(x) = Fᵀu(x)

**Subject to:**
- Volume constraint:  
  V(x) ≤ Vmax
- Design variables:  
  Cross-sectional areas of truss elements

Where:
- x → design variables (areas)
- u(x) → displacement vector
- F → external force vector

---

## 📊 Results

The optimization was performed for different volume fractions:

- α = 0.05
- α = 0.10
- α = 0.15
- α = 0.20

### Key Observations

- Increasing α → **lower compliance → higher stiffness**
- Material is distributed along **main load paths**
- Critical elements (e.g., 64 and 82) dominate load transfer
- Less critical members shrink in low-volume scenarios

This demonstrates effective **automatic material redistribution** by the optimization algorithm.

---

## 📂 Project Structure
