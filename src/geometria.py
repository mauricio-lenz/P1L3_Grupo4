"""
GEOMETRIA DE PISOS, MASAS SISMICAS Y FUERZAS LATERALES
=======================================================
P1L3 - Grupo 4

Funciones auxiliares para:

  * Carga viva        : intensidades q_Q por nivel (misma geometria
                        tributaria que el caso G de Semana 2).
  * Sismo pseudoestatico: centro de masa (CM) por piso, peso sismico
                        W_i = peso propio + 50% sobrecarga de uso y
                        fuerza lateral F_i = (fraccion_g) * W_i.

Unidades SI coherentes: m, kN, kN/m2, kN*m, m/s2.
"""


def carga_losa(cfg, casos):
    """Devuelve (qG, qQ): dicts {nivel: kN/m2} sin la clave 'nota'."""
    qG = {k: v for k, v in cfg["cargas"]["qG"].items() if k != "nota"}
    qQ = {k: v for k, v in casos["carga_viva_q_Q_kN_m2"].items()
          if k != "nota"}
    return qG, qQ


def centros_de_masa(cfg):
    """Centro de masa por piso considerando la losa real (rectangulo
    principal + voladizo eje J en Piso3/Piso4).

    Supuesto: la masa por unidad de area es uniforme en el piso, por lo
    que el CM coincide con el centroide del area de losa.

    Devuelve {nivel: (A_m2, x_cm_m, y_cm_m)}.
    """
    xc = [cfg["grilla_ejes"]["X"][k] for k in cfg["orden_x"]]
    yc = [cfg["grilla_ejes"]["Y"][k] for k in cfg["orden_y"]]
    x0, x1 = min(xc), max(xc)
    y0, y1 = min(yc), max(yc)
    a_main = (x1 - x0) * (y1 - y0)
    cx_main = (x0 + x1) / 2.0
    cy_main = (y0 + y1) / 2.0

    out = {}
    for lvl in cfg["orden_niveles"]:
        a = a_main
        cx = cx_main
        cy = cy_main
        for v in cfg.get("voladizos", {}).get("lista", []):
            if lvl not in v.get("niveles", []):
                continue
            xv0 = cfg["grilla_ejes"]["X"][v["desde_eje"]]
            xv1 = v["x_j_m"] + v.get("ancho_saliente_m", 0.0)
            yv = [cfg["grilla_ejes"]["Y"][e] for e in v["ejes_Y"]]
            yv0, yv1 = min(yv), max(yv)
            av = (xv1 - xv0) * (yv1 - yv0)
            cxv = (xv0 + xv1) / 2.0
            cyv = (yv0 + yv1) / 2.0
            cx = (cx * a + cxv * av) / (a + av)
            cy = (cy * a + cyv * av) / (a + av)
            a += av
        out[lvl] = (a, cx, cy)
    return out


def peso_estructura_por_piso(cfg):
    """Peso propio de la estructura (columnas + vigas + muros nucleo)
    tributario a cada piso.

    Determina el 'peso propio' estructural (ademas de la losa, que ya
    esta en qG) para la masa sismica. Por nivel:

      columnas: P-70x70 en todos los cruces de grilla, altura de piso.
      vigas   : V-60x80 en toda la grilla del nivel.
      muros   : muros del nucleo e=20 cm (los que pasan por el nivel).

    Devuelve {nivel: kN}.
    """
    xc = [cfg["grilla_ejes"]["X"][k] for k in cfg["orden_x"]]
    yc = [cfg["grilla_ejes"]["Y"][k] for k in cfg["orden_y"]]
    dx = [xc[i + 1] - xc[i] for i in range(len(xc) - 1)]
    dy = [yc[j + 1] - yc[j] for j in range(len(yc) - 1)]
    nx, ny = len(xc), len(yc)
    a_col = 0.70 * 0.70
    a_viga = 0.60 * 0.80
    gamma = cfg["materiales"][list(cfg["materiales"].keys())[0]]["gamma"]

    # longitudes de muros (solo los que participan en todos los niveles)
    mur_len = 0.0
    for m in cfg["muros"]["lista"]:
        mur_len += _largo_muro(cfg, m)

    out = {}
    for k, lvl in enumerate(cfg["orden_niveles"]):
        zprev = cfg["niveles"]["base"] if k == 0 else \
            cfg["niveles"][cfg["orden_niveles"][k - 1]]
        h = cfg["niveles"][lvl] - zprev
        vol_col = nx * ny * a_col * h
        vol_viga = (ny * sum(dx) + nx * sum(dy)) * a_viga
        vol_muro = mur_len * 0.20 * h
        out[lvl] = (vol_col + vol_viga + vol_muro) * gamma
    return out


def _largo_muro(cfg, m):
    if m["orientacion"] == "X":
        x0 = cfg["grilla_ejes"]["X"][m["desde_eje"]]
        x1 = cfg["grilla_ejes"]["X"][m["hasta_eje"]]
        return abs(x1 - x0)
    y0 = cfg["grilla_ejes"]["Y"][m["desde_eje"]]
    y1 = cfg["grilla_ejes"]["Y"][m["hasta_eje"]]
    return abs(y1 - y0)


def pesos_sismicos(cfg, casos):
    """Peso sismico por piso:  W_i = PP_i + 0.5 * Q_i.

      PP_i = qG_i * A_i (losa + terminaciones) + peso estructura
             (columnas/vigas/muros tributarios).
      Q_i  = qQ_i * A_i (sobrecarga de uso).

    Devolver {nivel: {"W_kN": ..., "PP_kN": ..., "Q_kN": ...}}.
    """
    qG, qQ = carga_losa(cfg, casos)
    cm = centros_de_masa(cfg)
    pp_est = peso_estructura_por_piso(cfg)
    frac = casos["sismo_pseudoestatico"]["fraccion_sobrecarga_masa"]
    out = {}
    for lvl in cfg["orden_niveles"]:
        a = cm[lvl][0]
        pp = qG[lvl] * a + pp_est[lvl]
        qside = qQ[lvl] * a
        out[lvl] = {"W_kN": pp + frac * qside,
                    "PP_kN": pp, "Q_kN": qside,
                    "area_m2": a, "qG": qG[lvl], "qQ": qQ[lvl]}
    return out


def fuerzas_laterales(cfg, casos):
    """Fuerza lateral en el CM de cada piso: F_i = fraccion_g * W_i.

    Admite sobreescritura manual por el profesor:
      sismo_pseudoestatico.fuerzas_manuales_EX / EY:
         [{"nivel": "...", "fuerza_kN": ...}, ...]

    Devuelve {"EX": {nivel: kN}, "EY": {nivel: kN}} (positivo segun
    +X y +Y globales).
    """
    frac = casos["sismo_pseudoestatico"]["fraccion_g"]
    ws = pesos_sismicos(cfg, casos)
    ex = {lvl: ws[lvl]["W_kN"] * frac for lvl in cfg["orden_niveles"]}
    ey = {lvl: ws[lvl]["W_kN"] * frac for lvl in cfg["orden_niveles"]}
    manual_x = casos["sismo_pseudoestatico"].get("fuerzas_manuales_EX")
    manual_y = casos["sismo_pseudoestatico"].get("fuerzas_manuales_EY")
    if manual_x:
        ex = {d["nivel"]: d["fuerza_kN"] for d in manual_x}
    if manual_y:
        ey = {d["nivel"]: d["fuerza_kN"] for d in manual_y}
    return {"EX": ex, "EY": ey}