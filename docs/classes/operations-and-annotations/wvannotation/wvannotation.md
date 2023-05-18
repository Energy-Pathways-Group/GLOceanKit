---
layout: default
title: WVAnnotation
parent: WVAnnotation
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVAnnotation

create a new instance of WVAnnotation


---

## Declaration
```matlab
 wvAnnotation = WVAnnotation(name,description,options)
```
## Parameters
+ `name`  name of the method, property, or variable
+ `description`  short description of the method, property, or variable
+ `detailedDescription`  (optional) a detailed description of the method, property, or variable

## Returns
+ `wvAnnotation`  a new instance of WVAnnotation

## Discussion

  Creates a new instance of WVAnnotation with a name,
  description and optional detailed description.
 
  If a markdown file of the same name is in the same directory
  or child directory, it will be loaded as the detailed
  description upon initialization.
 
            
