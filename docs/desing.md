# Diseño

En este documento se trabajará el el análisis del problema pasando por varias etapas antes del diseño del algoritmo propiamente, y antes de la implementación en código.

---

### Entender el problema

En el trabajo con bacteriofagos con capacidades líticas para el control de patógenos o cepas de interés se pueden obtener a través de diversas vías decenas, centenas y hasta miles de bacteriófagos diferentes. Estos pueden tener host-range(s), capacidad lítica, y hasta morfologías similares, pero finalmente se terminan diferenciando por leves diferencias ya sean fenotípicas (*caso1*) o genotípicas (*caso2*).

Sea cual sea el caso, realizar una selección que pretenda generar el mejor cóctel manualmente puede derivar en cesgos por la cantidad de datos manejados. En ese sentido es necesario desarrollar un software que facilite la visualización de las caracteristicas de cada fago aislado, y que genere candidatos para la integración de un cóctel.

Asimismo, es necesario contemplar que en el trabajo sistematizado de fagos, es necesario delimitar aquellos fagos a los cuales se les realice secuenciación, pues muchos de ellos terminarían siendo material no útil. Por lo mismo, dentro de este mismo programa, se plantean dos casos para la selección de fagos para la integración de un cóctel

**Caso 1:** El primer caso tiene por objetivo la delimitación del material biológico candidato a ser secuenciado para su posterior análisis. A través de este mismo caso también se generan cócteles candidatos para su selección posterior en el ´caso 2´.

**Caso 2:** El segundo caso surge una vez se cuenta con gran material biológico del que se tiene además información genómica. A través de este una selección mucho más fina puede ser lograda, con la llegada de nuevos parámetros a considerar. Se pretende que al realizar una selección con el *caso 2* se propongan versiones preliminares para su implementación comercial.

El siguiente matíz del problema radica en la necesidad de generar cócteles complementarios que abarquen la totalidad, o la mayoría de las cepas. Aquí hay dos enfoques, cada uno con su posible repertorio de ventajas y desventajas posibles. Complementación solapante, y no solapante. La solapante pretendería generar cócteles robustos con varias "balas" para un mismo objetivo. La no solapante, por otro lado, espera estar actuando sobre las cepas con mecanismos diferentes, ganando así flexibilidad y mayor posibilidad de abarcar cepas no contempladas.

---

### Definición de qué debe hacer el programa

Este programa está dividido en dos principales casos, en el primero que no se tiene información genética y que pretende delimitar el material a secuenciar, y el segundo caso que serviría para la integración aplicativa de un cóctel para uso.


***Caso 1***

El primer caso para ser útil debe parametrizar cada uno de los bacteriofagos aislados. Esto tendrá el fin último de generar un top 'n' fagos como candidatos para secuenciación.
Para ello, la información recibida por el programa será un archivo en formato .tsv que contenga los datos del host-range de los fagos a considerar, así como el perfil fenotípico de presencia de depolimeras (de haber dicho perfil). El tipo de dato de entrada en el host-range obtenido de spot-test será un rango del 0 - 4. Cada nivel del rango indica lo siguiente:

**0)** Actividad nula / no detectada.

**1)** Actividad mínima visible, fenotipo tipo *foggy*.

**2)** Actividad moderada: fenotipo *mod foggy* + *res prevalent*.

**3)** Actividad moderada: fenotipo *res prevalent*.

**4)** Actividad total, lisis completa.

Adicionalmente una de las columnas en la tabla podrá contener los datos del fenotipo de depolimerasa observado. (+) En caso de la presencia del fenotipo; (-) En caso de ausencia del fenotipo. De modo que el formato del `input` será:

|CEPA/FAGO | F1 | F2 | F3 | F4 | T |
|----------|----|----|----|----|---|
|C1        | 2  | 0  | 4  | 1  | 7 |
|C2        | 3  | 1  | 4  | 2  | 10|
|C3        | 3  | 0  | 3  | 2  | 8 |
|DP        | +  | -  | -  | -  |n/a|
|T         | 8  | 1  | 11 | 5  | 25|

El algoritmo por diseñar deberá ser capaz de leer los valores de las filas y las columnas generando puntajes para los fagos, así como para las cepas. Los valores tipo `output` generados serán:

#### Ranking de fagos:

**PhageID:** cada fago se le asigna el `ID` correspondiente a la entrada. Cada uno recibe las siguientes métricas.
**Raw-Score:** Corresponde a la sumatoria de los valores de su host-range.
**StrongCov:** Obtenido a partir de su fracción de cobertura con nivel "=4" (`fracscore4`).
**GoodCov:** Corresponde a la cobertura útil, con nivel superior a 3 (cov>=3).
**WideCov:** Métrica que corresponde a la amplitud 'útil', n° de cepas cubiertas con cov>=2.
**OhmPrev:** Esta métrica sirve para identificar la prevalencia de resistentes por fago; cov=<2.
**Depo:** (+) ó (-), presencia de depolimerasas.
**GlobalScore:** Ponderación configurable, a partir de este se genera el ranking de fagos.

#### Métricas p/cepa (suceptibility-score)(s-score)

> Las métricas por cepa tienen como objetivo principal identificar tanto a las cepas más suceptibles como posibles candidatas para ser cepas de producción de fagos, así como aquellas cepas que representen vulnerabilidad al tener poca o nula suceptibilidad.

**StrainID:** cada cepa se le asigna el `ID` correspondiente a la entrada. Cada uno recibe las siguientes métricas.
**Raw-Score:** Corresponde a la sumatoria de los valores de suceptibility-range.
**StrongSub:** Obtenido a partir de su fracción de suceptibilidad con nivel "=4" (`fracscore4`).
**ModSub:** Corresponde a la suceptibilidad significativa, con nivel superior a 3 (sub>=3). Estas últimas dos metricas representan buena suceptibilidad y posibilidad de fungir como cepas productoras de fagos.
**OhmPrev:** Esta métrica sirve para identificar la prevalencia de resistencia en la cepa; sub=<2.
**Vulnerability:** Esta mide la fracción cuya suceptibilidad sea sub=0. Estas últimas dos representan vulnerabilidad y necesidad de más aislamientos que las cubran.
**S-Score:** Ponderación configurable, a partir de este se genera un informe con las cepas más suceptibles y las más vulnerables.

#### Integración de cócteles ~ Índice de similitud de fagos

Una vez se tenga un listado de fagos prioritarios, se buscará generar cócteles, pero antes es necesario revisar los parametros necesarios para lo mismo. Aquí surgen los dos enfoques; con y sin complementación solapante. Solamente una vez que se tenga el listado de fagos prioritarios, el programa deberá generar índices de similitud con los fagos candidatos, de este modo ahorrando procesamiento innecesario. El programa tomará de manera iterativa uno de los fagos candidatos y se generarán índices de similitud con los demás fagos candidatos, se seleccionará uno de los fagos para cada caso (con y sin solapamiento), se tomará la unión de ambos conjuntos y se repite el ejercicio para 'n' fagos que se desee integren el cóctel. Los índices de similitud se no se imprimirán de manera explícita para el usuario, pero se guardarán si el mismo desea verlo. Sirven para generar los cócteles. Adicionalmente, de haber alguno, se seleccionarán fagos con evidencia empírica de depolimerasas, optando por el que mejor cumpla con el resto de parámetros. Un set de cócteles para cada caso será propuesto.


#### Sequence Priority

>Esto define la meta explícita del `caso1`. A partir de los datos recopilados, particularmente del ranking de fagos, así como de los cócteles obtenidos por índices de similitud se determinan cuáles son los fagos cuyo perfil les permite avanzar a la selección más fina, priorizando su análisis genómico.

**GlobalScore:** Tomado como el primer criterio de los fagos a secuenciar.


**NOTA:** Los cócteles propuestos en el `caso1` son meramente demostrativos y preliminares. Sirven para tener una guía y dinámica de trabajo y selección más fluida, no deben ser tomados de manera definitiva, sino como material de prueba, así como para priorizar secuenciación. A partir de los datos genómicos se generarán cócteles más finos lo que corresponde al `caso2`.

