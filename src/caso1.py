import csv
import math
import os


# ─────────────────────────────────────────────────────────────────────────────
# LECTURA
# ─────────────────────────────────────────────────────────────────────────────

def leer_input(ruta):
    matriz = {}
    ids_cepas = []
    ids_fagos = []
    depo = None

    with open(ruta, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        filas = list(reader)

    # Primera fila: IDs de fagos (excluir columna 0 y columna "T")
    cabecera = filas[0]
    ids_fagos = [col for col in cabecera[1:] if col.strip() != "T"]

    for fila in filas[1:]:
        if not fila:
            continue
        id_ = fila[0].strip()

        if id_ == "T":
            continue

        valores = fila[1:]
        # Alinear valores con ids_fagos (excluir columna T)
        valores_fagos = []
        j = 0
        for col_idx, col_name in enumerate(cabecera[1:]):
            if col_name.strip() == "T":
                j += 1
                continue
            val = valores[j].strip() if j < len(valores) else ""
            valores_fagos.append(val)
            j += 1

        if id_ == "DP":
            depo = {}
            for fago, val in zip(ids_fagos, valores_fagos):
                depo[fago] = val if val in ("+", "-") else "-"
        else:
            ids_cepas.append(id_)
            matriz[id_] = {}
            for fago, val in zip(ids_fagos, valores_fagos):
                try:
                    matriz[id_][fago] = int(val)
                except (ValueError, TypeError):
                    matriz[id_][fago] = 0

    return matriz, ids_cepas, ids_fagos, depo


# ─────────────────────────────────────────────────────────────────────────────
# MÉTRICAS POR CEPA
# ─────────────────────────────────────────────────────────────────────────────

def metricas_cepa(matriz, ids_cepas, ids_fagos):
    n_fagos = len(ids_fagos)
    ranking_cepas = []
    cepas_dificiles = set()

    for cepa in ids_cepas:
        scores = [matriz[cepa][fago] for fago in ids_fagos]

        raw_score = sum(scores)
        n4 = sum(1 for s in scores if s >= 4)
        n3 = sum(1 for s in scores if s == 3)
        n_score0 = sum(1 for s in scores if s == 0)

        if n_fagos == 0:
            s_score = 0.0
            vulnerability = 0.0
        else:
            s_score = (n4 * 1.0 + n3 * 0.5) / n_fagos
            vulnerability = n_score0 / n_fagos

        alerta = None
        if s_score == 0:
            alerta = "SIN COBERTURA"
        elif s_score < 0.25:
            alerta = "VULNERABLE"
            cepas_dificiles.add(cepa)

        ranking_cepas.append({
            "StrainID":      cepa,
            "RawScore":      raw_score,
            "SScore":        round(s_score, 4),
            "Vulnerability": round(vulnerability, 4),
            "Alerta":        alerta if alerta else ""
        })

    ranking_cepas.sort(key=lambda x: x["SScore"], reverse=True)
    return ranking_cepas, cepas_dificiles


# ─────────────────────────────────────────────────────────────────────────────
# MÉTRICAS POR FAGO
# ─────────────────────────────────────────────────────────────────────────────

def metricas_fago(matriz, ids_fagos, ids_cepas, depo, cepas_dificiles):
    n_cepas = len(ids_cepas)
    ranking_fagos = []
    cobertura_dificiles = {}

    for fago in ids_fagos:
        scores = [matriz[cepa][fago] for cepa in ids_cepas]

        raw_score   = sum(scores)
        n_cubiertas = sum(1 for s in scores if s >= 1)
        n_wide      = sum(1 for s in scores if s >= 2)
        n4          = sum(1 for s in scores if s >= 4)
        n3          = sum(1 for s in scores if s == 3)
        n_ohm       = sum(1 for s in scores if s <= 2)

        if n_cepas == 0:
            frac_cov          = 0.0
            wide_cov          = 0.0
            strong_cov_global = 0.0
            ohm_prev          = 0.0
        else:
            frac_cov          = n_cubiertas / n_cepas
            wide_cov          = n_wide / n_cepas
            strong_cov_global = (n4 * 1.0 + n3 * 0.5) / n_cepas
            ohm_prev          = n_ohm / n_cepas

        if n_cubiertas == 0:
            strong_cov_local = 0.0
        else:
            strong_cov_local = (n4 * 1.0 + n3 * 0.5) / n_cubiertas

        strong_cov = (strong_cov_global + strong_cov_local) / 2

        depo_fago = None if depo is None else depo.get(fago, None)

        cepas_fuertes = [cepa for cepa in ids_cepas if matriz[cepa][fago] >= 3]
        if not cepas_fuertes:
            rare_cov = 0.0
        else:
            n_rare = sum(1 for cepa in cepas_fuertes if cepa in cepas_dificiles)
            rare_cov = n_rare / len(cepas_fuertes)

        global_score = (wide_cov * 0.5 + strong_cov * 1.0 + (1 - ohm_prev) * 0.25) / 1.75

        ranking_fagos.append({
            "PhageID":     fago,
            "RawScore":    raw_score,
            "FracCov":     round(frac_cov, 4),
            "WideCov":     round(wide_cov, 4),
            "StrongCov":   round(strong_cov, 4),
            "OhmPrev":     round(ohm_prev, 4),
            "Depo":        depo_fago if depo_fago is not None else "N/A",
            "RareCov":     round(rare_cov, 4),
            "GlobalScore": round(global_score, 4)
        })

    ranking_fagos.sort(key=lambda x: x["GlobalScore"], reverse=True)

    # Construir cobertura_dificiles sobre todos los fagos
    gs_map = {r["PhageID"]: r["GlobalScore"] for r in ranking_fagos}
    for cepa in cepas_dificiles:
        fagos_rescue = [
            fago for fago in ids_fagos if matriz[cepa][fago] >= 3
        ]
        fagos_rescue.sort(
            key=lambda f: (matriz[cepa][f], gs_map.get(f, 0)),
            reverse=True
        )
        cobertura_dificiles[cepa] = fagos_rescue

    return ranking_fagos, cobertura_dificiles


# ─────────────────────────────────────────────────────────────────────────────
# SELECCIÓN DE CANDIDATOS
# ─────────────────────────────────────────────────────────────────────────────

def seleccionar_candidatos(ranking_fagos, n):
    n = min(n, len(ranking_fagos))
    return ranking_fagos[:n]


# ─────────────────────────────────────────────────────────────────────────────
# GREEDY SET COVER
# ─────────────────────────────────────────────────────────────────────────────

def greedy(candidatos, matriz, ids_cepas, cepas_dificiles,
           cobertura_dificiles, mode, tamaño):

    coctel = []
    ids_candidatos = [f["PhageID"] for f in candidatos]
    gs_map   = {f["PhageID"]: f["GlobalScore"] for f in candidatos}
    depo_map = {f["PhageID"]: f["Depo"] for f in candidatos}

    if mode == "non_overlapping":
        cubiertas = set()
    else:
        pesos = {cepa: 1.0 for cepa in ids_cepas}
        for cepa in cepas_dificiles:
            pesos[cepa] = 2.0

    # ── Fase 1: Greedy base ───────────────────────────────────────────────

    while len(coctel) < tamaño:
        mejor_fago     = None
        mejor_ganancia = -1

        for fago in ids_candidatos:
            if fago in coctel:
                continue

            if mode == "non_overlapping":
                ganancia = 0
                for cepa in ids_cepas:
                    if cepa not in cubiertas and matriz[cepa][fago] >= 2:
                        ganancia += 2 if cepa in cepas_dificiles else 1
            else:
                ganancia = sum(
                    pesos[cepa]
                    for cepa in ids_cepas
                    if matriz[cepa][fago] >= 2
                )

            if ganancia > mejor_ganancia:
                mejor_ganancia = ganancia
                mejor_fago = fago
            elif ganancia == mejor_ganancia and mejor_fago is not None:
                if gs_map.get(fago, 0) > gs_map.get(mejor_fago, 0):
                    mejor_fago = fago
                elif gs_map.get(fago, 0) == gs_map.get(mejor_fago, 0):
                    if depo_map.get(fago) == "+" and depo_map.get(mejor_fago) != "+":
                        mejor_fago = fago

        if mejor_fago is None:
            break

        coctel.append(mejor_fago)

        if mode == "non_overlapping":
            for cepa in ids_cepas:
                if matriz[cepa][mejor_fago] >= 2:
                    cubiertas.add(cepa)
        else:
            for cepa in ids_cepas:
                if matriz[cepa][mejor_fago] >= 2:
                    pesos[cepa] *= 0.5

    # ── Fase 2: Rescate de cepas difíciles ───────────────────────────────

    cepas_irrescatables = []

    for cepa in cepas_dificiles:
        cubierta = any(matriz[cepa][fago] >= 2 for fago in coctel)
        if cubierta:
            continue

        opciones = cobertura_dificiles.get(cepa, [])
        if not opciones:
            cepas_irrescatables.append(cepa)
            continue

        if mode == "non_overlapping":
            rescue = opciones[0]
            if rescue not in coctel:
                coctel.append(rescue)
        else:
            mejor_score = matriz[cepa][opciones[0]]
            for fago in opciones:
                if matriz[cepa][fago] == mejor_score:
                    if fago not in coctel:
                        coctel.append(fago)
                else:
                    break

    return coctel, cepas_irrescatables


# ─────────────────────────────────────────────────────────────────────────────
# GENERACIÓN DE CÓCTELES
# ─────────────────────────────────────────────────────────────────────────────

def generar_coctel(candidatos, matriz, ids_cepas, cepas_dificiles,
                   cobertura_dificiles, tamaño):

    coctel_no_solapante, irrescatables_no = greedy(
        candidatos, matriz, ids_cepas, cepas_dificiles,
        cobertura_dificiles, mode="non_overlapping", tamaño=tamaño
    )
    coctel_solapante, irrescatables_sol = greedy(
        candidatos, matriz, ids_cepas, cepas_dificiles,
        cobertura_dificiles, mode="overlapping", tamaño=tamaño
    )

    # Unión de irrescatables entre ambos modos
    cepas_irrescatables = list(set(irrescatables_no) | set(irrescatables_sol))

    return coctel_no_solapante, coctel_solapante, cepas_irrescatables


# ─────────────────────────────────────────────────────────────────────────────
# SIMILITUD COSENO
# ─────────────────────────────────────────────────────────────────────────────

def similitud_coseno(candidatos, matriz, ids_cepas):
    ids_candidatos = [f["PhageID"] for f in candidatos]
    n = len(ids_candidatos)

    vectores = {
        fago: [matriz[cepa][fago] for cepa in ids_cepas]
        for fago in ids_candidatos
    }

    sim = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i, n):
            A = vectores[ids_candidatos[i]]
            B = vectores[ids_candidatos[j]]

            dot    = sum(a * b for a, b in zip(A, B))
            norm_a = math.sqrt(sum(a ** 2 for a in A))
            norm_b = math.sqrt(sum(b ** 2 for b in B))

            cos = 0.0 if norm_a == 0 or norm_b == 0 else dot / (norm_a * norm_b)
            sim[i][j] = round(cos, 4)
            sim[j][i] = round(cos, 4)

    return sim, ids_candidatos


# ─────────────────────────────────────────────────────────────────────────────
# GUARDAR RESULTADOS
# ─────────────────────────────────────────────────────────────────────────────

def guardar_resultados(ranking_fagos, ranking_cepas, coctel_no_solapante,
                        coctel_solapante, cepas_irrescatables, directorio,
                        similitud=None, ids_candidatos_sim=None):

    os.makedirs(directorio, exist_ok=True)
    gs_map     = {f["PhageID"]: f for f in ranking_fagos}
    ids_cands  = {f["PhageID"] for f in gs_map.values()
                  if f["PhageID"] in [r["PhageID"] for r in ranking_fagos]}

    # ── ranking_fagos.tsv ─────────────────────────────────────────────────
    with open(os.path.join(directorio, "ranking_fagos.tsv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "PhageID", "RawScore", "FracCov", "WideCov", "StrongCov",
            "OhmPrev", "Depo", "RareCov", "GlobalScore"
        ], delimiter="\t")
        w.writeheader()
        w.writerows(ranking_fagos)

    # ── ranking_cepas.tsv ─────────────────────────────────────────────────
    with open(os.path.join(directorio, "ranking_cepas.tsv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "StrainID", "RawScore", "SScore", "Vulnerability", "Alerta"
        ], delimiter="\t")
        w.writeheader()
        w.writerows(ranking_cepas)
        if cepas_irrescatables:
            f.write(
                f"\n# CEPAS SIN FAGO DISPONIBLE (score>=3): "
                f"{', '.join(cepas_irrescatables)}\n"
            )

    # ── cócteles ──────────────────────────────────────────────────────────
    def escribir_coctel(nombre_archivo, coctel, ids_candidatos_set):
        with open(os.path.join(directorio, nombre_archivo), "w",
                  newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["PhageID", "GlobalScore", "StrongCov",
                        "WideCov", "Depo", "Origen"])
            for fago_id in coctel:
                datos = gs_map.get(fago_id, {})
                origen = "candidato" if fago_id in ids_candidatos_set else "rescate"
                w.writerow([
                    fago_id,
                    datos.get("GlobalScore", "N/A"),
                    datos.get("StrongCov", "N/A"),
                    datos.get("WideCov", "N/A"),
                    datos.get("Depo", "N/A"),
                    origen
                ])

    escribir_coctel("cocktail_non_overlapping.tsv",
                    coctel_no_solapante, ids_cands)
    escribir_coctel("cocktail_overlapping.tsv",
                    coctel_solapante, ids_cands)

    # ── sequence_priority.tsv ─────────────────────────────────────────────
    fagos_rescate = (
        set(coctel_no_solapante) | set(coctel_solapante)
    ) - ids_cands

    priority = []
    for f in ranking_fagos:
        if f["PhageID"] in ids_cands or f["PhageID"] in fagos_rescate:
            priority.append({
                "PhageID":     f["PhageID"],
                "GlobalScore": f["GlobalScore"],
                "RareCov":     f["RareCov"],
                "Depo":        f["Depo"],
                "Prioridad":   "rescate" if f["PhageID"] in fagos_rescate
                               else "candidato"
            })

    with open(os.path.join(directorio, "sequence_priority.tsv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "PhageID", "GlobalScore", "RareCov", "Depo", "Prioridad"
        ], delimiter="\t")
        w.writeheader()
        w.writerows(priority)

    # ── similarity_matrix.tsv (opcional) ──────────────────────────────────
    if similitud is not None and ids_candidatos_sim is not None:
        with open(os.path.join(directorio, "similarity_matrix.tsv"), "w",
                  newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow([""] + ids_candidatos_sim)
            for i, fago_id in enumerate(ids_candidatos_sim):
                w.writerow([fago_id] + similitud[i])


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main(ruta, directorio, n, tamaño, guardar_similitud=False):

    matriz, ids_cepas, ids_fagos, depo = leer_input(ruta)

    ranking_cepas, cepas_dificiles = metricas_cepa(
        matriz, ids_cepas, ids_fagos
    )

    ranking_fagos, cobertura_dificiles = metricas_fago(
        matriz, ids_fagos, ids_cepas, depo, cepas_dificiles
    )

    candidatos = seleccionar_candidatos(ranking_fagos, n)

    coctel_no_solapante, coctel_solapante, cepas_irrescatables = generar_coctel(
        candidatos, matriz, ids_cepas, cepas_dificiles,
        cobertura_dificiles, tamaño
    )

    sim = None
    ids_cands_sim = None
    if guardar_similitud:
        sim, ids_cands_sim = similitud_coseno(candidatos, matriz, ids_cepas)

    guardar_resultados(
        ranking_fagos, ranking_cepas,
        coctel_no_solapante, coctel_solapante,
        cepas_irrescatables, directorio,
        sim, ids_cands_sim
    )
