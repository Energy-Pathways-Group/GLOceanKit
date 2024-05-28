---
layout: default
title: initWithRandomFlow
parent: WVTransform
grand_parent: Classes
nav_order: 109
mathjax: true
---

#  initWithRandomFlow

initialize with a random flow state


---

## Declaration
```matlab
 initWithRandomFlow(flowComponentNames)
```
## Parameters
+ `flowComponentNames`  strings of flow component names names.
+ `uvMax`  (optional) maximum horizontal velocity

## Discussion

  Clears variables Ap,Am,A0 and then randomizes the flow by adding random
  amplitudes at all available modes. Optionally, you can specify which
  components of the flow should get initialized. For example,
 
  ```matlab
    wvt.initWithRandomFlow();
  ```
 
  will initialize all modes, while
 
  ```matlab
    wvt.initWithRandomFlow('geostrophic','mda');
  ```
 
  will initialize the flow with geostrophic and mean density anomaly flow
  components, while the wave and inertial oscillations components will be
  zero.
 
        
