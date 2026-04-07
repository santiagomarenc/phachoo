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

En formato de tabla se presenta cada Fago con su métrica correspondiente, en orden descenciente tomando como base **GlobalScore**.

**PhageID:** cada fago se le asigna el `ID` correspondiente a la entrada. Cada uno recibe las siguientes métricas. (Lectura directa)

**Raw-Score:** Sumatoria de los valores de host-range del fago sobre todas las cepas. Valor entero, sin normalizar. Métrica informativa, no entra al GlobalScore.

**FracCov:** Fracción de cepas con score >= 1 sobre el total de cepas. Valor entre 0 y 1. Métrica de soporte interno usada como denominador en otros cálculos; se reporta al usuario pero no entra al GlobalScore.

**WideCov:** Fracción de cepas con score >= 2 sobre el total de cepas. Valor entre 0 y 1. Captura la amplitud de cobertura útil mínima. Entra al GlobalScore.

**StrongCov:** Métrica de calidad de cobertura. Se calcula a partir de dos conteos internos (no reportados): `n4` = cepas con score >= 4, `n3` = cepas con score = 3. La fórmula combina una perspectiva global y una local:

```
StrongCov_global = (n4 × w4 + n3 × w3) / n_cepas_totales
StrongCov_local  = (n4 × w4 + n3 × w3) / n_cepas_cubiertas   (cepas con score >= 1)
StrongCov        = (StrongCov_global + StrongCov_local) / 2
```

Pesos fijos: `w4 = 1.0`, `w3 = 0.5`. Valor entre 0 y 1. Entra al GlobalScore.

**OhmPrev:** Fracción de cepas con score <= 2 sobre el total de cepas. Identifica la prevalencia de cobertura débil o resistentes en la colección completa. Valor entre 0 y 1. Entra al GlobalScore con peso negativo.

**Depo:** (+) o (-), presencia de depolimerasas. Lectura directa desde el input. Valor nulo si no hay datos. No entra al GlobalScore; actúa como criterio de desempate en el algoritmo de selección.

**RareCov:** Fracción de las cepas cubiertas por el fago (score >= 3) que son clasificadas como cepas difíciles (S-Score < 0.25). Valor entre 0 y 1. Métrica interna: no se reporta al usuario pero alimenta el algoritmo greedy como bonus para fagos que cubren cepas huérfanas. Requiere que las métricas por cepa estén calculadas previamente.

**GlobalScore:** Promedio ponderado de `WideCov`, `StrongCov` y `OhmPrev` (con peso negativo). Pesos fijos en el Caso 1: `Pw = 0.5`, `Ps = 1.0`, `Po = 0.25`.

```
GlobalScore = (WideCov×0.5 + StrongCov×1.0 + (1 - OhmPrev)×0.25) / 1.75
```

Valor entre 0 y 1. A partir de este se genera el ranking de fagos y se seleccionan los candidatos al cóctel.

#### Métricas p/cepa (susceptibility-score)(S-Score)

Salida en formato de tabla, de cada cepa con sus métricas, en orden descendente tomando como base **S-Score**.

> Las métricas por cepa tienen como objetivo principal identificar tanto las cepas más susceptibles como posibles candidatas para ser cepas de producción de fagos, así como aquellas cepas que representen vulnerabilidad al tener poca o nula susceptibilidad. El S-Score actúa como eje único del ranking: los valores altos señalan cepas productoras candidatas y los valores bajos señalan cepas vulnerables.

**StrainID:** cada cepa recibe el `ID` correspondiente a la entrada.

**Raw-Score:** Sumatoria de los valores de susceptibilidad de la cepa frente a todos los fagos. Valor entero, informativo, no entra al S-Score.

**S-Score:** Perfil ponderado de susceptibilidad. Captura la acumulación de niveles altos (3 y 4) sobre el total de fagos, con mayor peso a los 4s. Se calcula a partir de dos conteos internos: `n4` = fagos con score >= 4 frente a la cepa, `n3` = fagos con score = 3 frente a la cepa. Pesos fijos: `w4 = 1.0`, `w3 = 0.5`.

```
S-Score = (n4 × w4 + n3 × w3) / n_fagos_totales
```

Valor entre 0 y 1. Cepas con S-Score alto son candidatas a producción. Cepas con S-Score < 0.25 se clasifican automáticamente como **cepas difíciles**: tienen cobertura insuficiente, se reportan como vulnerables y se marcan como necesidad de más aislamientos.

**Vulnerability:** Fracción de fagos con score = 0 frente a la cepa, sobre el total de fagos. Valor entre 0 y 1. Métrica informativa complementaria al S-Score para identificar cepas sin ningún fago activo.

```
Vulnerability = n_score0 / n_fagos_totales
```

**NOTA:** Al finalizar el cálculo de métricas por cepa se genera una estructura interna `cepas_dificiles` con todas las cepas cuyo S-Score < 0.25. Esta estructura es consumida por el cálculo de `RareCov` en las métricas por fago, por lo que las métricas por cepa deben calcularse primero.

**Alertas en el reporte:** se emiten para cepas cuyo S-Score = 0 (ningún fago las cubre) y para cepas clasificadas como difíciles (S-Score < 0.25).

#### Integración de cócteles ~ Índice de similitud de fagos

Una vez se tenga un listado de fagos prioritarios, se generan dos cócteles candidatos mediante un **algoritmo greedy de cobertura marginal**, uno para cada enfoque de complementación. Los índices de similitud coseno se calculan sobre los fagos candidatos y se guardan opcionalmente si el usuario desea explorarlos, pero no son el motor de selección del cóctel.

**Índice de similitud coseno**

Cada fago se representa como un vector de scores frente a las cepas (excluyendo la fila `DP`, que es categórica). La similitud entre dos fagos A y B se define como:

```
cos(A, B) = (A · B) / (||A|| × ||B||)
```

El resultado es un valor entre 0 y 1. Valores cercanos a 1 indican fagos con patrones de host-range similares (redundantes). Valores cercanos a 0 indican fagos complementarios. Esta matriz se genera una sola vez sobre los fagos candidatos para evitar procesamiento innecesario.

**Criterios base**

1) Cobertura de cada una de las cepas, o de las máximas posibles.
2) Prioridad a fagos con evidencia empírica de depolimerasas como criterio de desempate final.

**Algoritmo greedy set cover**

Primero se toman los 'n' fagos candidatos con mayor `GlobalScore`. A partir de ellos se construyen dos cócteles de tamaño configurable mediante un proceso en dos fases.

**Estructura previa — cobertura de cepas difíciles**

Antes de correr el greedy se construye una estructura sobre **todos los fagos** (no solo los candidatos) que mapea cada cepa difícil a los fagos que la cubren con score >= 3, ordenados por score descendente:

```
cobertura_dificiles = {
    cepa_difícil: [fagos ordenados por score desc frente a esa cepa]
}
```

Esta estructura es consumida en la Fase 2 del algoritmo.

**Ganancia marginal ajustada**

Las cepas difíciles (S-Score < 0.25) tienen peso `W = 2` dentro del cálculo de ganancia; las cepas normales tienen `W = 1`. Esto favorece naturalmente a los fagos que cubren cepas huérfanas sin introducir parámetros externos al greedy.

**Fase 1 — Greedy base**

*Cóctel no solapante:* en cada iteración se selecciona el fago que maximiza la suma de pesos de cepas aún no cubiertas (score >= 2):

```
ganancia = Σ W(cepa)  para cada cepa nueva con score >= 2
           donde W = 2 si cepa es difícil, W = 1 si no
```

Las cepas ya cubiertas se descartan del cálculo siguiente.

*Cóctel solapante:* en cada iteración se selecciona el fago que maximiza la suma ponderada de cobertura, donde cada cepa difícil parte con peso doble y el peso de toda cepa decae ×0.5 por cada fago que la cubre:

```
ganancia = Σ pesos[cepa] × W(cepa)  para cada cepa con score >= 2
           donde W = 2 si cepa es difícil, W = 1 si no
           y pesos[cepa] arranca en 1.0, decae ×0.5 por cada fago que la cubre
```

**Jerarquía de desempate en Fase 1** (igual ganancia ajustada):
1. Mayor `GlobalScore`
2. Presencia de `Depo` (+)

**Fase 2 — Rescate de cepas difíciles**

Al terminar el greedy base, se verifica qué cepas difíciles quedaron sin cobertura adecuada (ningún fago del cóctel las cubre con score >= 2). Para cada una se consulta `cobertura_dificiles` y se agrega al cóctel el mejor fago disponible:

- *Cóctel no solapante:* se toma el fago de mayor score frente a esa cepa. En caso de empate, se resuelve por `GlobalScore`.
- *Cóctel solapante:* en caso de empate en score, se agregan ambos fagos al cóctel.

Si una cepa difícil no tiene ningún fago con score >= 3 en toda la colección, se reporta como irrescatable y se emite alerta de necesidad de más aislamientos.

Los fagos agregados en la Fase 2 pueden estar fuera del top-n de candidatos originales. El tamaño final del cóctel puede exceder el 'n' configurado, lo cual es esperado y comunicado al usuario.

**En cada ciclo:** Se proponen dos cócteles construidos a partir de los 'n' candidatos más los fagos de rescate si aplica. Ambos apuntan a los `criterios base`.


#### Sequence Priority

>Esto define la meta explícita del `caso1`. A partir de los datos recopilados, particularmente del ranking de fagos, así como de los cócteles obtenidos por índices de similitud se determinan cuáles son los fagos cuyo perfil les permite avanzar a la selección más fina, priorizando su análisis genómico.

**GlobalScore:** Tomado como el primer criterio de los fagos a secuenciar.


**NOTA:** Los cócteles propuestos en el `caso1` son meramente demostrativos y preliminares. Sirven para tener una guía y dinámica de trabajo y selección más fluida, no deben ser tomados de manera definitiva, sino como material de prueba, así como para priorizar secuenciación. A partir de los datos genómicos se generarán cócteles más finos lo que corresponde al `caso2`.




## Algoritmo

El siguiente flujo describe el procesamiento completo del Caso 1, desde la lectura del archivo de entrada hasta la generación del output. Cada etapa produce estructuras de datos que alimentan a las siguientes.

---

### Etapa 1 — Lectura y parsing del input

El programa recibe un archivo `.tsv`. Se identifican y separan los siguientes componentes:

- **Matriz de scores:** submatriz numérica con filas = cepas y columnas = fagos. Se excluyen la fila `DP`, la fila `T` y la columna `T`.
- **Fila DP:** fila categórica con valores (+) o (-) por fago. Se almacena como un diccionario separado. Si no existe en el archivo, se registra como ausente y `Depo` retorna nulo para todos los fagos.
- **Fila y columna T:** totales presentes en el input para referencia del usuario. Se ignoran en los cálculos; el programa genera sus propios totales.

Se extraen además las listas de IDs de cepas y de fagos tal como aparecen en el input. Estas son la base de todas las estructuras posteriores.

---

### Etapa 2 — Métricas por cepa

Se procesa cada cepa de la matriz de scores calculando `Raw-Score`, `S-Score` y `Vulnerability`. Al finalizar se generan dos estructuras:

- **Ranking de cepas:** lista ordenada de mayor a menor S-Score.
- **`cepas_dificiles`:** conjunto de IDs de cepas cuyo S-Score < 0.25. Se emiten alertas para cepas con S-Score = 0.

Estas estructuras quedan disponibles para la Etapa 3.

---

### Etapa 3 — Métricas por fago

Con `cepas_dificiles` disponible, se procesa cada fago calculando `Raw-Score`, `FracCov`, `WideCov`, `StrongCov`, `OhmPrev`, `Depo`, `RareCov` y `GlobalScore`. Internamente se calculan los conteos `n3` y `n4` necesarios para `StrongCov`.

Al finalizar se construye la estructura de rescate sobre **todos los fagos**:

- **`cobertura_dificiles`:** diccionario que mapea cada cepa difícil a la lista de fagos que la cubren con score >= 3, ordenados por score descendente. En caso de empate de score, se ordenan por GlobalScore descendente.
- **Ranking de fagos:** lista ordenada de mayor a menor GlobalScore.

---

### Etapa 4 — Selección de candidatos

Se toman los top-n fagos del ranking por GlobalScore. El valor de `n` es configurable por el usuario. Esta lista de candidatos es el input de la Etapa 5.

---

### Etapa 5 — Generación de cócteles

Se ejecuta el greedy set cover en dos fases (Fase 1 greedy base + Fase 2 rescate de cepas difíciles) para cada modo (solapante y no solapante), produciendo dos cócteles. El proceso está descrito en detalle en la sección anterior.

---

### Etapa 6 — Similitud coseno (opcional)

Si el usuario lo solicita, se calcula la matriz de similitud coseno sobre los fagos candidatos del top-n. Cada fago se representa como su vector de scores sobre todas las cepas (excluyendo DP). La matriz se guarda en un archivo separado.

---

### Etapa 7 — Output

Todos los reportes se escriben como archivos `.tsv` en un directorio de salida. No hay impresión en consola. Los archivos generados son:

- **ranking_fagos.tsv:** tabla con todas las métricas por fago, ordenada por GlobalScore descendente.
- **ranking_cepas.tsv:** tabla con S-Score y Vulnerability por cepa, ordenada por S-Score descendente. Incluye columna de alerta para cepas difíciles e irrescatables.
- **cocktail_non_overlapping.tsv:** lista de fagos del cóctel no solapante con su contribución de cobertura por cepa.
- **cocktail_overlapping.tsv:** lista de fagos del cóctel solapante con su contribución ponderada por cepa.
- **sequence_priority.tsv:** fagos del top-n ordenados por GlobalScore, más los fagos de rescate incorporados en la Etapa 5.
- **similarity_matrix.tsv** (opcional): matriz de similitud coseno entre candidatos, generada solo si el usuario lo solicitó.

---

## Estructura

```
main()
├── leer_input()
├── metricas_cepa()
├── metricas_fago()
├── seleccionar_candidatos()
├── generar_coctel()
│   ├── greedy(mode="non_overlapping")
│   └── greedy(mode="overlapping")
├── similitud_coseno()
└── guardar_resultados()
```

**`main()`** — Función orquestadora. Recibe los parámetros configurables del usuario (ruta del archivo, directorio de salida, n candidatos, tamaño del cóctel, umbral de cobertura, si guardar similitud), llama a cada función en orden y gestiona el flujo de datos entre ellas.

**`leer_input(ruta)`** — Lee el archivo `.tsv` y separa sus componentes: matriz de scores, fila DP (si existe) y listas de IDs de cepas y fagos. Retorna estas estructuras limpias para su uso posterior.

**`metricas_cepa(matriz, ids_cepas, ids_fagos)`** — Calcula Raw-Score, S-Score y Vulnerability para cada cepa. Retorna el ranking de cepas ordenado por S-Score y la estructura `cepas_dificiles`. Convención: si el denominador es 0, la métrica retorna 0.

**`metricas_fago(matriz, ids_fagos, ids_cepas, depo, cepas_dificiles)`** — Calcula Raw-Score, FracCov, WideCov, StrongCov, OhmPrev, Depo, RareCov y GlobalScore para cada fago. Construye `cobertura_dificiles`. Retorna el ranking de fagos y `cobertura_dificiles`. Convención: si el denominador es 0, StrongCov_local retorna 0.

**`seleccionar_candidatos(ranking_fagos, n)`** — Toma los top-n fagos del ranking por GlobalScore. Retorna la lista de candidatos.

**`greedy(candidatos, matriz, cepas_dificiles, cobertura_dificiles, mode)`** — Ejecuta el algoritmo greedy en dos fases (base + rescate) para el modo indicado. Retorna el cóctel resultante como lista de fagos.

**`generar_coctel(candidatos, matriz, cepas_dificiles, cobertura_dificiles)`** — Llama a `greedy()` en ambos modos y retorna los dos cócteles.

**`similitud_coseno(candidatos, matriz)`** — Calcula la matriz de similitud coseno entre los fagos candidatos. Solo se ejecuta si el usuario lo solicita. Retorna la matriz y la guarda como `similarity_matrix.tsv`.

**`guardar_resultados(ranking_fagos, ranking_cepas, coctel_no_solapante, coctel_solapante, directorio, similitud)`** — Escribe todos los reportes de salida como archivos `.tsv` en el directorio indicado: rankings, cócteles, alertas y lista de prioridad de secuenciación.
