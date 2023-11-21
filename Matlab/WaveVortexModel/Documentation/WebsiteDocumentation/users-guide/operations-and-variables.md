---
layout: default
title: Operations and variables
parent: Users guide
mathjax: true
nav_order: 4
---

#  Operations and variables

This document is meant to articulate the internal logic used for operations and variable annotations.

There are a couple of important points,
1. One operation can produce many variables.
2. Variables may be cached to prevent re-computation.

When the user adds an operation, the following takes place
1. The operation is associated with its name in the `operationNameMap'.
2. Each variable is added to the `variableAnnotationNameMap' and possibly the `timeDependentVariablesNameMap'.
3. Each variable is associated with this operation using the `operationVariableNameMap'.

When the user requests a variable, the following takes place
1. Check the variable cache to see if its already been computed (if so, return).
2. Use the `variableAnnotationNameMap', to find the variable annotation, and then retrieve the operation name
3. Use `operationNameMap' to retrieve the operation.
4. Run the operation, cache and return the results.

This seems convoluted. Why can we not just use `operationVariableNameMap'?

So now the user adds an operation which produces a variable that is already being used. We certainly could come up with a way to pick and choose which variables get computed. But in the short-term we need to just discard the other operation that produced the variables.
1. 
