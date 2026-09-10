"""
SUPERPOSICION DE CASOS P1L3
===========================
P1L3 - Grupo 4

Parte C: con los casos G, Q, EX y EY se construye la combinacion
arbitraria

    R = lG*G + lQ*Q + lEX*EX + lEY*EY

primero sumando las respuestas unitarias (superposicion) y luego se
compara contra una corrida OpenSees explicita con las cargas combinadas
aplicadas de una sola vez. Se verifican desplazamiento, reaccion y
fuerza interna.
"""

import numpy as np
import openseespy.opensees as ops

from . import geometria as geo
from .casos import (Resultado, aplicar_sismo, configurar_y_analizar,
                    construir_edificio, correr_caso_gravedad,
                    seleccionar_elementos)
from .edificio import areas_tributarias, aplicar_cargas


def combinar(cfg, resultados, lambdas):
    """Suma lineal de las respuestas: R = sum(l*resultado).

    resultados: {caso: Resultado}; lambdas: {caso: float}.
    """
    desp = {lvl: np.zeros(6) for lvl in cfg["orden_niveles"]}
    reac = {}
    fuerzas = {}
    totales = np.zeros(6)
    for caso, lam in lambdas.items():
        r = resultados[caso]
        if abs(lam) < 1e-15:
            continue
        for lvl in cfg["orden_niveles"]:
            desp[lvl] += lam * r.desp[lvl]
        for nid, v in r.reac.items():
            reac[nid] = reac.get(nid, np.zeros(6)) + lam * v
        for nombre, v in r.fuerzas.items():
            fuerzas[nombre] = fuerzas.get(nombre, np.zeros_like(v)) + lam * v
        totales += lam * r.reacciones_totales
    comb = Resultado.__new__(Resultado)  # objeto liviano sin OPS
    comb.desp = desp
    comb.reac = reac
    comb.fuerzas = fuerzas
    comb._reac_tot = totales
    return comb


def correr_explicita(cfg, casos):
    """Corrida OpenSees con TODAS las cargas a la vez:
      gravedad   : lG*G + lQ*Q (patron 1)
      sismo      : lEX*EX + lEY*EY (patron 2)

    Devuelve (Resultado, info_cargas).
    """
    lam = casos["superposicion"]["lambdas"]
    qG, qQ = geo.carga_losa(cfg, casos)
    fuerzas = geo.fuerzas_laterales(cfg, casos)

    ed = construir_edificio(cfg)
    tags = seleccionar_elementos(ed)

    qc = {lvl: lam["G"] * qG[lvl] + lam["Q"] * qQ[lvl]
          for lvl in cfg["orden_niveles"]}
    w_map, resumen, area_piso_lvl, ax, ay, w_vol = \
        areas_tributarias(cfg, qc)
    total_grav = aplicar_cargas(cfg, ed, w_map, w_vol,
                                tag_patron=1, tag_series=1,
                                con_puntuales=False)

    fdir = {"EX": {lvl: lam["EX"] * fuerzas["EX"][lvl]
                   for lvl in cfg["orden_niveles"]},
            "EY": {lvl: lam["EY"] * fuerzas["EY"][lvl]
                   for lvl in cfg["orden_niveles"]}}
    aplicar_sismo(cfg, ed, fdir, escala=1.0, tag_patron=2, tag_series=2)

    configurar_y_analizar()

    info = {"gravedad_total_kN": total_grav,
            "lateral_EX_kN": sum(fdir["EX"].values()),
            "lateral_EY_kN": sum(fdir["EY"].values()),
            "w_map": w_map, "resumen": resumen}
    return Resultado(ed, tags), info


def comparar(comb, expl, cfg):
    """Compara superposicion vs corrida explicita.

    Verifica (al menos): desplazamiento (techo), reaccion (basal total
    y una reaccion puntual) y fuerza interna (elemento de control).

    Devuelve un dict con el resumen y una tabla lista para el informe.
    """
    rows = []

    def _err(a, b):
        den = max(np.max(np.abs(b)), np.max(np.abs(a)), 1e-12)
        return float(np.max(np.abs(np.asarray(a) - np.asarray(b))) / den)

    top = cfg["orden_niveles"][-1]
    # 1) desplazamiento del techo (master)
    d_c = comb.desp[top]
    d_e = expl.desp[top]
    rows.append({"variable": "Desplazamiento techo (master)",
                 "superposicion": d_c.tolist(),
                 "explicita": d_e.tolist(),
                 "error_rel": _err(d_c, d_e)})

    # 2) reaccion total basal
    r_c = comb.reacciones_totales
    r_e = expl.reacciones_totales
    rows.append({"variable": "Reacciones basales totales (global)",
                 "superposicion": r_c.tolist(),
                 "explicita": r_e.tolist(),
                 "error_rel": _err(r_c, r_e)})

    # 3) reaccion en un nodo base puntual
    nid = min(expl.reac.keys())
    n_c = comb.reac[nid]
    n_e = expl.reac[nid]
    rows.append({"variable": f"Reaccion nodo base {nid} (local)",
                 "superposicion": n_c.tolist(),
                 "explicita": n_e.tolist(),
                 "error_rel": _err(n_c, n_e)})

    # 4) fuerza interna (fuerzas basicas del elemento de control)
    for nombre, tag in expl.tags_control.items():
        if not tag:
            continue
        f_c = comb.fuerzas[nombre]
        f_e = expl.fuerzas[nombre]
        rows.append({"variable": f"Fuerza interna {nombre} (tag {tag})",
                     "superposicion": f_c.tolist(),
                     "explicita": f_e.tolist(),
                     "error_rel": _err(f_c, f_e)})

    max_err = max(r["error_rel"] for r in rows)
    return {"tabla": rows,
            "error_relativo_maximo": float(max_err),
            "ok": bool(max_err < 1e-6)}