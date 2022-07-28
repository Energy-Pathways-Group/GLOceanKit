---
layout: default
title: Amt
parent: WVTransform
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  Amt

negative wave coefficients at time t


---

## Description
Complex valued state variable with dimensions $$(k,l,j)$$ and units of $$m/s$$.

## Discussion

These are the *time dependent* coefficients of the internal gravity wave and inertial oscillation portion of the flow, denoted  $$A_- e^{-i \omega t} $$ in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995).

Unlike `Am`, these coefficients do not have their phases wound to time $$t_0$$ (`t0`), and **do** change for linear dynamics.

