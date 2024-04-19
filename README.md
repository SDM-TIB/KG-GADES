# KG-GADES
This is an updated version of GADES in Python: https://github.com/SDM-TIB/GADES, the original version of GADES (in Java) is here: https://github.com/collarad/GADES/tree/master. 

The reference paper is: `Ignacio Traverso-Ribón, Maria-Esther Vidal, Benedikt Kämpgen and York Sure-Vetter. GADES: A Graph-based Semantic Similarity Measure. Semantics 2016, Leipzig.`

# Implementation Detail
![kg_gades](https://github.com/SDM-TIB/KG-GADES/blob/main/kg_gades.png)

# Example
We test this package with the `example.ttl` file. 
```python
from kg_gades import MyOWLOntology

MyOWLOntology("exmaple.ttl")
```
