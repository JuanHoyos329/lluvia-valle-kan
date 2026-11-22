# **Lluvia en el Valle del Cauca con datos funcionales y redes KAN**

En el Valle del Cauca la lluvia tiene dos temporadas al año y varía con el ENSO y con el relieve. Las estaciones de la CVC la miden con precisión, pero son pocas y tienen huecos. El satélite cubre todo el departamento, aunque con sesgos en zonas de montaña.

Este repositorio reúne el trabajo para unir las dos fuentes entre 1981 y 2003. Primero se limpian e imputan los datos de las estaciones. Luego cada año de lluvia se representa como una curva suave, se describe su comportamiento y se mide cuánto se aleja el satélite de las estaciones. Por último, se entrena una red funcional con núcleo KAN que reconstruye la curva de una estación a partir de la del satélite, y se compara con un MLP, con un modelo lineal y con el satélite sin corregir.

El objetivo es tener una estimación de la lluvia más confiable donde no hay estaciones, y que además se pueda interpretar a la luz de la climatología del Valle.
