# KG-GADES
This is an updated version of GADES in Python: https://github.com/SDM-TIB/GADES, the original version of GADES (in Java) is here: https://github.com/collarad/GADES/tree/master. 



# Implementation Detail
![kg_gades](https://github.com/SDM-TIB/KG-GADES/blob/main/kg_gades.png)

# Example
We test this package with the `example.ttl` file. 
```python
from kg_gades import MyOWLOntology

MyOWLOntology("exmaple.ttl")
```

# References
[1] Ignacio Traverso-Ribón, Maria-Esther Vidal, Benedikt Kämpgen and York Sure-Vetter. GADES: A Graph-based Semantic Similarity Measure. Semantics 2016, Leipzig.

[2] Ignacio Traverso-Ribón, Maria-Esther Vidal and Guillermo Palma. A Similarity Measure for Determining Relatedness Between Ontology Terms. 11th International Conference on Data Integration in the Life Sciences 2015 (DILS2015).

[3] Ignacio Traverso-Ribón, Maria-Esther Vidal. Exploiting Information Content and Semantics to Accurately Computing Similarity of GO-based Annotated Entities. 2015 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology

[4] Palma, G.; Vidal, M. E.; Haag, E.; Raschid, L.; Thor, A. Measuring Relatedness Between Scientific Entities in Annotation Datasets. ACM International Conference on Bioinformatics, Computational Biology, and Biomedical Informatics (BCB), 2013
