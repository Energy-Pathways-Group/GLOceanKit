---
layout: default
title: addOperation
parent: WVTransform
grand_parent: Classes
nav_order: 64
mathjax: true
---

#  addOperation

add a WVOperation


---

## Discussion

  Several things happen when adding an operation.
  1. We check that dimensions exist for all output variables
  produced by this operation.
  2. We see if there are any existing output variables with the
  same name.
    2a. We remove the operation that produced the existing
    variables, if it exists.
  3. We map each new variable to this operation variableAnnotationNameMap
  4. Map each operation name to the operation
 
  In our revision,
    - The variableAnnotationNameMap will map the name to the
    variable annotation
    - The operationVariableNameMap will map the name to the
    operation
  
  
