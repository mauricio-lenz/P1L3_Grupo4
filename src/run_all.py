"""
P1L3 - Orquestador
==================
Ejecuta las cuatro partes del LAB P1L3 Grupo 4:

    python -m src.run_all

Parte A: caso de carga viva Q sobre la misma geometria tributaria de
         Semana 2 (verificacion de conservacion Q = q_Q * A).
Parte B: sismo pseudoestatico EX y EY con fuerza lateral en el CM de
         cada piso (F_i = 0.20 * W_i) y verificaciones.
Parte C: superposicion R = lG*G + lQ*Q + lEX*EX + lEY*EY (suma de
         respuestas unitarias) comparada contra una corrida OpenSees
         explicita con las cargas combinadas (desplazamiento, reaccion
         y fuerza interna).
Parte D: capacidad de la columna HA con Fiber Section - M-phi y
         primeros puntos P-M.

Todo lo configurable (intensidades, sismo, lambdas, seccion HA) esta en
data/casos_lab3.json. Los resultados y figuras van a results/.
"""

import json
from pathlib import Path

import numpy as np

from . import capacidad as cap
from . import casos as cs
from . import geometria as geo
from . import superposicion as sup
from .edificio import cargar_config

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results"
FIG = RES / "figures"
CSV = RES / "csv"
EXP = RES / "export"


def _guardar_csv(filename, rows, header):
    (CSV / filename).write_text(
        header + "\n" + "\n".join(",".join(str(x) for x in r) for r in rows)
        + "\n", encoding="utf-8")


def main():
    for d in (FIG, CSV, EXP):
        d.mkdir(parents=True, exist_ok=True)

    cfg = cargar_config(DATA / "edificio_config.json")
    casos = json.loads((DATA / "casos_lab3.json").read_text(encoding="utf-8"))
    cfg["cap_hormigon"] = casos["cap_hormigon"]
    qG, qQ = geo.carga_losa(cfg, casos)
    fuerzas = geo.fuerzas_laterales(cfg, casos)
    lambdas = casos["superposicion"]["lambdas"]

    res = {}

    # ---- Parte A/B: casos base G, Q, EX, EY (modelos independientes)
    for nombre, qmap, punt in (("G", qG, True), ("Q", qQ, False)):
        ed = cs.construir_edificio(cfg)
        tags = cs.seleccionar_elementos(ed)
        w_map, resumen, apisos, ax, ay, w_vol, total = \
            cs.correr_caso_gravedad(cfg, ed, qmap, con_puntuales=punt)
        res[nombre] = cs.registrar(cfg, ed, tags)
        fz = sum(r[2] for r in res[nombre].reac.values())
        res[nombre].meta = {"total_aplicado_kN": total,
                            "suma_resumen_kN": sum(r["carga_total"]
                                                   for r in resumen.values()),
                            "suma_reacciones_FZ": fz}
        if nombre == "Q":
            verifA = cs.verificar_conservacion_qv(cfg, qQ, resumen, apisos,
                                                  total, fz)
    for d in ("EX", "EY"):
        ed = cs.construir_edificio(cfg)
        tags = cs.seleccionar_elementos(ed)
        cs.correr_caso_sismo(cfg, ed, d, fuerzas[d])
        res[d] = cs.registrar(cfg, ed, tags)

    verifB = {d: cs.verificar_sismo(cfg, casos, fuerzas[d], res[d], d)
              for d in ("EX", "EY")}

    # ---- Parte C: superposicion vs corrida explicita combinada
    comb = sup.combinar(cfg, res, lambdas)
    expl, info_expl = sup.correr_explicita(cfg, casos)
    cmp_ = sup.comparar(comb, expl, cfg)

    # Demanda de la columna critica (combinacion) para marcar P-M
    f_col = expl.fuerzas["columna_critica"]
    demanda_pm = {"P": float(f_col[0]),
                  "M": float(np.hypot(f_col[4], f_col[5]))}

    # ---- Parte D: capacidad de la columna
    curvas_pm = cap.run_capacidad(cfg, FIG, demanda=demanda_pm)
    prop_sec = cap.propiedades_seccion(cfg)
    fibras, A_br = cap.discretizacion_fibras(cfg)
    n_fib_conf = cap.graficar_seccion(cfg, FIG / "fig_seccion_fibras.png")

    # ---- Salidas CSV
    dirs = ("G", "Q", "EX", "EY")
    hdr = "nivel," + ",".join("%s_%s" % (d, c)
                              for d in dirs
                              for c in "ux uy uz rx ry rz".split())
    _guardar_csv(
        "desplazamientos_master.csv",
        [[lvl, *(round(res[d].master(lvl)[i], 8)
                 for d in dirs for i in range(6))]
         for lvl in cfg["orden_niveles"]],
        hdr)

    tab = cmp_["tabla"]
    _guardar_csv(
        "superposicion.csv",
        [[r["variable"], r["error_rel"], r["superposicion"], r["explicita"]]
         for r in tab],
        "variable,error_relativo,superposicion,explicita")
    _guardar_csv(
        "cortes_basales.csv",
        [[d, verifB[d]["carga_lateral_total_kN"],
          verifB[d]["corte_basal_kN"],
          verifB[d]["desplazamiento_techo_m"],
          verifB[d]["sentido_deformada_correcto"]] for d in ("EX", "EY")],
        "caso,carga_lateral_kN,corte_basal_kN,desp_techo_m,sentido_ok")
    _guardar_csv(
        "pesos_sismicos.csv",
        [[lvl, round(ws["area_m2"], 2), round(ws["PP_kN"], 1),
          round(ws["Q_kN"], 1), round(ws["W_kN"], 1)]
         for lvl, ws in geo.pesos_sismicos(cfg, casos).items()],
        "nivel,area_m2,PP_kN,Q_kN,W_kN")
    body = []
    for P in curvas_pm["Ps"]:
        if P == curvas_pm["P_axial"]:     # punto de compresion axial pura (M=0)
            body.append([P, 0.0, 0.0, 0.0, 0.0])
            continue
        phis, Ms = curvas_pm["curvas"][P]
        i = int(np.argmax(Ms))
        body.append([P, Ms.max(), phis[i], Ms[-1], phis[-1]])
    _guardar_csv("cap_pm.csv", body, "P_kN,Mu_kN_m,phi_Mu,tail_M,tail_phi")
    _guardar_csv(
        "seccion.csv",
        [[k, v] for k, v in prop_sec.items()],
        "propiedad,valor")

    # ---- Informe en pantalla
    tot_area = geo.pesos_sismicos(cfg, casos)
    print("=" * 70)
    print("P1L3 - EDIFICIO U. DE LOS ANDES 2017_67 (Grupo 4)")
    print("=" * 70)
    print("Parte A (carga viva Q): Q total = %.1f kN; transferida = %.1f kN"
          % (verifA["q_total_kN"], verifA["carga_transferida_kN"]))
    print("   conservacion area: %s | reacciones: %s"
          % (verifA["ok_conservacion_area"], verifA["ok_conservacion_reacciones"]))
    for d in ("EX", "EY"):
        v = verifB[d]
        print("Parte B (%s): W_sismico total = %.1f kN, F = %.1f kN"
              % (d, sum(w["W_kN"] for w in tot_area.values()),
                 v["carga_lateral_total_kN"]))
        print("   corte basal = %.1f kN | desp techo = %.4f m | sentido ok: %s"
              % (v["corte_basal_kN"], v["desplazamiento_techo_m"],
                 v["sentido_deformada_correcto"]))
        print("   torsion techo Rz = %.2e rad" % v["torsion_por_piso"][-1]["Rz_rad"])
    print("Parte C: R = %.2f G + %.2f Q + %.2f EX + %.2f EY"
          % (lambdas["G"], lambdas["Q"], lambdas["EX"], lambdas["EY"]))
    print("   error relativo maximo (superposicion vs explicita): %.3e  OK=%s"
          % (cmp_["error_relativo_maximo"], cmp_["ok"]))
    print("   fuerza interna columna: %s" % [round(float(x), 2)
                                             for x in f_col])
    print("Parte D: seccion %.0fx%.0f, As=%.2f cm2, cuantia=%.4f"
          % (prop_sec["b_m"] * 100, prop_sec["h_m"] * 100,
             prop_sec["As_m2"] * 1e4, prop_sec["cuantia"]))
    print("   Puntos P-M (primeros):")
    for P, M in zip(curvas_pm["Ps"], curvas_pm["Mu"]):
        print("     P = %+7.0f kN   ->  Mu = %7.1f kN-m" % (P, M))
    print("   Demanda columna critica: P = %.1f kN, M = %.1f kN-m"
          % (demanda_pm["P"], demanda_pm["M"]))
    print("   Curvatura de fluencia estimada (d=%.2f m): %.5f 1/m"
          % (prop_sec["d_m"],
             prop_sec["fy_MPa"] * 1e3 / prop_sec["Es_kN_m2"]
             / (0.7 * prop_sec["d_m"])))
    print("Resumen seccion: %d fibras nucleo/recubrimiento, %d barras"
          % (n_fib_conf["n_fibras_concreto"], n_fib_conf["n_barras"]))


if __name__ == "__main__":
    main()