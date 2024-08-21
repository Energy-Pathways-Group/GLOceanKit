---
layout: default
title: solutionForModeAtIndex
parent: WVPrimaryFlowComponent
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  solutionForModeAtIndex

return the analytical solution for the mode at this index


---

## Declaration
```matlab
 solution = solutionForModeAtIndex(index)
```
## Parameters
+ `index`  non-negative integer less than nModes
+ `amplitude`  (optional) 'wvt' or 'random' (default)

## Returns
+ `solution`  an instance of WVAnalyticalSolution

## Discussion

  Returns WVAnalyticalSolution object for this index.
  The solution indices run from 1:nModes.
 
  The solution amplitude can be set to either 'wvt' or
  'random'. Setting the amplitude='wvt' will use the amplitude
  currently set in the wvt to initialize this solution.
  Otherwise an appropriate random amplitude will be created.
 
          
