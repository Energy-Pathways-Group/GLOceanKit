---
layout: default
title: addRandomFlow
parent: WVTransform
grand_parent: Classes
nav_order: 62
mathjax: true
---

#  addRandomFlow

add randomized flow to the existing state


---

## Declaration
```matlab
 addRandomFlow(flowComponentNames,options)
```
## Parameters
+ `flowComponentNames`  strings of flow component names names.
+ `uvMax`  (optional) maximum horizontal velocity

## Discussion

  Adds random amplitudes at all available modes. Optionally, you can
  specify which components of the flow should get initialized. For example,
 
  ```matlab
    wvt.addRandomFlow();
  ```
 
  will add noise at all modes, while
 
  ```matlab
    wvt.addRandomFlow('geostrophic','mda');
  ```
 
  will add random flow at only thegeostrophic and mean density anomaly flow
  components, while the wave and inertial oscillations components will
  remain untouched
 
        
