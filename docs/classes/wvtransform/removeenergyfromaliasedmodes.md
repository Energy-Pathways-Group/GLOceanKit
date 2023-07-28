---
layout: default
title: removeEnergyFromAliasedModes
parent: WVTransform
grand_parent: Classes
nav_order: 147
mathjax: true
---

#  removeEnergyFromAliasedModes

remove all energy from aliased modes


---

## Declaration
```matlab
 removeEnergyFromAliasedModes(options)
```
## Parameters
+ `jFraction`  (optional) fraction of vertical mode to assume are not aliased (default 2/3)

## Discussion

  Removes the energy from modes that will alias with a quadratic
  multiplication. If running a nonlinear simulation, it is recommended that
  you perform this step before time-stepping.
      
