# Pseudocódigo — Caso 1

---

## `leer_input(ruta)`

```
FUNCIÓN leer_input(ruta):

  1. Abrir el archivo .tsv en ruta
  2. Leer la primera fila → ids_fagos (excluir columna 0 y columna "T")
  3. Leer las filas restantes:
     PARA cada fila:
       id = valor columna 0
       SI id == "DP":
         depo ← { fago: valor para cada fago en ids_fagos }
         valores faltantes en DP → "-" por defecto
       SI id == "T":
         ignorar
       SINO:
         ids_cepas.append(id)
         matriz[id] ← { fago: int(valor) para cada fago en ids_fagos }
         valores faltantes en scores → 0 por defecto

  4. SI "DP" no fue encontrada:
     depo ← None

  RETORNAR matriz, ids_cepas, ids_fagos, depo
```

---

## `metricas_cepa(matriz, ids_cepas, ids_fagos)`

```
FUNCIÓN metricas_cepa(matriz, ids_cepas, ids_fagos):

  n_fagos ← len(ids_fagos)
  ranking_cepas ← []
  cepas_dificiles ← conjunto vacío

  PARA cada cepa en ids_cepas:

    scores ← [matriz[cepa][fago] para cada fago en ids_fagos]

    raw_score     ← suma(scores)
    n4            ← contar scores donde score >= 4
    n3            ← contar scores donde score == 3
    n_score0      ← contar scores donde score == 0

    SI n_fagos == 0:
      s_score       ← 0
      vulnerability ← 0
    SINO:
      s_score       ← (n4×1.0 + n3×0.5) / n_fagos
      vulnerability ← n_score0 / n_fagos

    alerta ← None
    SI s_score == 0:
      alerta ← "SIN COBERTURA"
    SINO SI s_score < 0.25:
      alerta ← "VULNERABLE"
      cepas_dificiles.add(cepa)

    ranking_cepas.append({
      StrainID:      cepa,
      RawScore:      raw_score,
      SScore:        s_score,
      Vulnerability: vulnerability,
      Alerta:        alerta
    })

  ranking_cepas ← ordenar por SScore descendente

  RETORNAR ranking_cepas, cepas_dificiles
```

---

## `metricas_fago(matriz, ids_fagos, ids_cepas, depo, cepas_dificiles)`

```
FUNCIÓN metricas_fago(matriz, ids_fagos, ids_cepas, depo, cepas_dificiles):

  n_cepas ← len(ids_cepas)
  ranking_fagos ← []
  cobertura_dificiles ← {}

  PARA cada fago en ids_fagos:

    scores ← [matriz[cepa][fago] para cada cepa en ids_cepas]

    raw_score   ← suma(scores)
    n_cubiertas ← contar scores donde score >= 1
    n_wide      ← contar scores donde score >= 2
    n4          ← contar scores donde score >= 4
    n3          ← contar scores donde score == 3
    n_ohm       ← contar scores donde score <= 2

    SI n_cepas == 0:
      frac_cov          ← 0
      wide_cov          ← 0
      strong_cov_global ← 0
      ohm_prev          ← 0
    SINO:
      frac_cov          ← n_cubiertas / n_cepas
      wide_cov          ← n_wide / n_cepas
      strong_cov_global ← (n4×1.0 + n3×0.5) / n_cepas
      ohm_prev          ← n_ohm / n_cepas

    SI n_cubiertas == 0:
      strong_cov_local ← 0
    SINO:
      strong_cov_local ← (n4×1.0 + n3×0.5) / n_cubiertas

    strong_cov ← (strong_cov_global + strong_cov_local) / 2

    SI depo == None:
      depo_fago ← None
    SINO:
      depo_fago ← depo[fago]

    cepas_fuertes ← [cepa para cada cepa en ids_cepas donde matriz[cepa][fago] >= 3]
    SI len(cepas_fuertes) == 0:
      rare_cov ← 0
    SINO:
      n_rare   ← contar cepas_fuertes donde cepa IN cepas_dificiles
      rare_cov ← n_rare / len(cepas_fuertes)

    global_score ← (wide_cov×0.5 + strong_cov×1.0 + (1 - ohm_prev)×0.25) / 1.75

    ranking_fagos.append({
      PhageID:     fago,
      RawScore:    raw_score,
      FracCov:     frac_cov,
      WideCov:     wide_cov,
      StrongCov:   strong_cov,
      OhmPrev:     ohm_prev,
      Depo:        depo_fago,
      RareCov:     rare_cov,
      GlobalScore: global_score
    })

  ranking_fagos ← ordenar por GlobalScore descendente

  PARA cada cepa en cepas_dificiles:
    fagos_rescue ← [fago para cada fago en ids_fagos
                    donde matriz[cepa][fago] >= 3]
    fagos_rescue ← ordenar por matriz[cepa][fago] desc,
                   luego por GlobalScore desc (desempate)
    cobertura_dificiles[cepa] ← fagos_rescue

  RETORNAR ranking_fagos, cobertura_dificiles
```

---

## `seleccionar_candidatos(ranking_fagos, n)`

```
FUNCIÓN seleccionar_candidatos(ranking_fagos, n):

  SI n > len(ranking_fagos):
    n ← len(ranking_fagos)

  candidatos ← ranking_fagos[0:n]

  RETORNAR candidatos
```

---

## `greedy(candidatos, matriz, ids_cepas, cepas_dificiles, cobertura_dificiles, mode, tamaño)`

```
FUNCIÓN greedy(candidatos, matriz, ids_cepas, cepas_dificiles,
               cobertura_dificiles, mode, tamaño):

  coctel ← []
  ids_candidatos ← [fago[PhageID] para cada fago en candidatos]

  SI mode == "non_overlapping":
    cubiertas ← conjunto vacío

  SI mode == "overlapping":
    pesos ← { cepa: 1.0 para cada cepa en ids_cepas }
    PARA cada cepa en cepas_dificiles:
      pesos[cepa] ← 2.0

  # ── FASE 1: Greedy base ──────────────────────────────────────

  REPETIR hasta len(coctel) == tamaño:

    mejor_fago     ← None
    mejor_ganancia ← -1

    PARA cada fago en ids_candidatos - coctel:

      SI mode == "non_overlapping":
        ganancia ← 0
        PARA cada cepa en ids_cepas:
          SI cepa NOT IN cubiertas AND matriz[cepa][fago] >= 2:
            SI cepa IN cepas_dificiles:
              ganancia += 2
            SINO:
              ganancia += 1

      SI mode == "overlapping":
        ganancia ← 0
        PARA cada cepa en ids_cepas:
          SI matriz[cepa][fago] >= 2:
            ganancia += pesos[cepa]

      SI ganancia > mejor_ganancia:
        mejor_ganancia ← ganancia
        mejor_fago ← fago

      SINO SI ganancia == mejor_ganancia:
        SI GlobalScore[fago] > GlobalScore[mejor_fago]:
          mejor_fago ← fago
        SINO SI GlobalScore[fago] == GlobalScore[mejor_fago]:
          SI Depo[fago] == "+" AND Depo[mejor_fago] != "+":
            mejor_fago ← fago

    coctel.append(mejor_fago)

    SI mode == "non_overlapping":
      PARA cada cepa en ids_cepas:
        SI matriz[cepa][mejor_fago] >= 2:
          cubiertas.add(cepa)

    SI mode == "overlapping":
      PARA cada cepa en ids_cepas:
        SI matriz[cepa][mejor_fago] >= 2:
          pesos[cepa] ← pesos[cepa] × 0.5

  # ── FASE 2: Rescate de cepas difíciles ───────────────────────

  PARA cada cepa en cepas_dificiles:

    cubierta ← FALSO
    PARA cada fago en coctel:
      SI matriz[cepa][fago] >= 2:
        cubierta ← VERDADERO
        BREAK

    SI NOT cubierta:
      SI cobertura_dificiles[cepa] está vacío:
        registrar cepa como irrescatable
        CONTINUAR

      SI mode == "non_overlapping":
        rescue_fago ← cobertura_dificiles[cepa][0]
        SI rescue_fago NOT IN coctel:
          coctel.append(rescue_fago)

      SI mode == "overlapping":
        mejor_score_rescue ← matriz[cepa][cobertura_dificiles[cepa][0]]
        PARA cada fago en cobertura_dificiles[cepa]:
          SI matriz[cepa][fago] == mejor_score_rescue AND fago NOT IN coctel:
            coctel.append(fago)
          SINO SI matriz[cepa][fago] < mejor_score_rescue:
            BREAK

  RETORNAR coctel
```

---

## `generar_coctel(candidatos, matriz, ids_cepas, cepas_dificiles, cobertura_dificiles, tamaño)`

```
FUNCIÓN generar_coctel(candidatos, matriz, ids_cepas, cepas_dificiles,
                        cobertura_dificiles, tamaño):

  coctel_no_solapante ← greedy(candidatos, matriz, ids_cepas,
                                cepas_dificiles, cobertura_dificiles,
                                mode="non_overlapping", tamaño)

  coctel_solapante    ← greedy(candidatos, matriz, ids_cepas,
                                cepas_dificiles, cobertura_dificiles,
                                mode="overlapping", tamaño)

  RETORNAR coctel_no_solapante, coctel_solapante
```

---

## `similitud_coseno(candidatos, matriz, ids_cepas)`

```
FUNCIÓN similitud_coseno(candidatos, matriz, ids_cepas):

  ids_candidatos ← [fago[PhageID] para cada fago en candidatos]
  n ← len(ids_candidatos)
  similitud ← matriz n×n inicializada en 0.0

  vectores ← {}
  PARA cada fago en ids_candidatos:
    vectores[fago] ← [matriz[cepa][fago] para cada cepa en ids_cepas]

  PARA i en 0..n-1:
    PARA j en i..n-1:

      A ← vectores[ids_candidatos[i]]
      B ← vectores[ids_candidatos[j]]

      dot    ← suma(A[k] × B[k] para cada k)
      norm_A ← raiz(suma(A[k]² para cada k))
      norm_B ← raiz(suma(B[k]² para cada k))

      SI norm_A == 0 OR norm_B == 0:
        cos ← 0.0
      SINO:
        cos ← dot / (norm_A × norm_B)

      similitud[i][j] ← cos
      similitud[j][i] ← cos

  RETORNAR similitud, ids_candidatos
```

---

## `guardar_resultados(ranking_fagos, ranking_cepas, coctel_no_solapante, coctel_solapante, cepas_irrescatables, directorio, similitud, ids_candidatos)`

```
FUNCIÓN guardar_resultados(ranking_fagos, ranking_cepas, coctel_no_solapante,
                            coctel_solapante, cepas_irrescatables, directorio,
                            similitud=None, ids_candidatos=None):

  # ── ranking_fagos.tsv ─────────────────────────────────────────
  escribir directorio/ranking_fagos.tsv:
    cabecera: PhageID, RawScore, FracCov, WideCov, StrongCov,
              OhmPrev, Depo, RareCov, GlobalScore
    PARA cada fago en ranking_fagos:
      escribir fila con sus valores

  # ── ranking_cepas.tsv ─────────────────────────────────────────
  escribir directorio/ranking_cepas.tsv:
    cabecera: StrainID, RawScore, SScore, Vulnerability, Alerta
    PARA cada cepa en ranking_cepas:
      escribir fila con sus valores
    SI cepas_irrescatables no está vacío:
      añadir al final: "CEPAS SIN FAGO DISPONIBLE (score>=3): [ids]"

  # ── cocktail_non_overlapping.tsv ──────────────────────────────
  escribir directorio/cocktail_non_overlapping.tsv:
    cabecera: PhageID, GlobalScore, StrongCov, WideCov, Depo, Origen
    PARA cada fago en coctel_no_solapante:
      origen ← "rescate" SI fago no estaba en candidatos, SINO "candidato"
      escribir fila con sus valores y origen

  # ── cocktail_overlapping.tsv ──────────────────────────────────
  escribir directorio/cocktail_overlapping.tsv:
    cabecera: PhageID, GlobalScore, StrongCov, WideCov, Depo, Origen
    PARA cada fago en coctel_solapante:
      origen ← "rescate" SI fago no estaba en candidatos, SINO "candidato"
      escribir fila con sus valores y origen

  # ── sequence_priority.tsv ─────────────────────────────────────
  escribir directorio/sequence_priority.tsv:
    cabecera: PhageID, GlobalScore, RareCov, Depo, Prioridad
    PARA cada fago en ranking_fagos que esté en candidatos o en rescate:
      prioridad ← "rescate" SI vino de Fase 2, SINO "candidato"
      escribir fila ordenada por GlobalScore desc

  # ── similarity_matrix.tsv (opcional) ──────────────────────────
  SI similitud != None:
    escribir directorio/similarity_matrix.tsv:
      cabecera: [vacío] + ids_candidatos
      PARA cada i en ids_candidatos:
        escribir fila: ids_candidatos[i] + [similitud[i][j] para cada j]
```

---

## `main(ruta, directorio, n, tamaño, guardar_similitud)`

```
FUNCIÓN main(ruta, directorio, n, tamaño, guardar_similitud=False):

  matriz, ids_cepas, ids_fagos, depo ← leer_input(ruta)

  ranking_cepas, cepas_dificiles ← metricas_cepa(matriz, ids_cepas, ids_fagos)

  ranking_fagos, cobertura_dificiles ← metricas_fago(matriz, ids_fagos,
                                                      ids_cepas, depo,
                                                      cepas_dificiles)

  candidatos ← seleccionar_candidatos(ranking_fagos, n)

  coctel_no_solapante, coctel_solapante ← generar_coctel(candidatos, matriz,
                                                          ids_cepas,
                                                          cepas_dificiles,
                                                          cobertura_dificiles,
                                                          tamaño)

  cepas_irrescatables ← cepas donde cobertura_dificiles[cepa] está vacío

  SI guardar_similitud:
    similitud, ids_candidatos ← similitud_coseno(candidatos, matriz, ids_cepas)
  SINO:
    similitud      ← None
    ids_candidatos ← None

  guardar_resultados(ranking_fagos, ranking_cepas,
                     coctel_no_solapante, coctel_solapante,
                     cepas_irrescatables, directorio,
                     similitud, ids_candidatos)
```
