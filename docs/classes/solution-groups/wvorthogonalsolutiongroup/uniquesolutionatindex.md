---
layout: default
title: uniqueSolutionAtIndex
parent: WVOrthogonalSolutionGroup
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  uniqueSolutionAtIndex

return the analytical solution at this solution index


---

## Declaration
```matlab
 solution = uniqueSolutionAtIndex(index)
```
## Parameters
+ `index`  non-negative integer less than nUniqueSolutions
+ `amplitude`  (optional) 'wvt' or 'random' (default)

## Returns
+ `solution`  an instance of WVAnalyticalSolution

## Discussion

  Returns WVAnalyticalSolution object for this solution index.
  The solution indices run from 1:nUniqueSolutions.
 
  The solution amplitude can be set to either 'wvt' or
  'random'. Setting the amplitude='wvt' will use the amplitude
  currently set in the wvt to initialize this solution.
  Otherwise an appropriate random amplitude will be created.
 
          
