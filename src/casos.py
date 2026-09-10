"""
CASOS BASE P1L3: G, Q, EX, EY
==============================
P1L3 - Grupo 4

Parte A (carga viva): misma geometria tributaria de Semana 2 con
intensidad q_Q. Verifica sum(Q transferida) = q_Q * A.

Parte B (sismo pseudoestatico): casos independientes EX y EY. La fuerza
lateral se aplica en el centro de masa de cada piso como fraccion de g
multiplicada por la masa del piso (W_i = peso propio + 50% sobrecarga).
Verifica carga lateral total, corte basal, sentido de la deformada y
torsion de piso.

El patron y los parametros son configurables desde data/casos_lab3.json
(sin tocar el codigo).
"""

import numpy as np
import openseespy.opensees as ops

from . import geometria as geo
from .edificio import Edificio, areas_tributarias, aplicar_cargas


# =====================================================================
# CONSTRUCCION Y ANALISIS
# =====================================================================

def construir_edificio(cfg):
    """Construye el modelo (nodos, apoyos, vigas, columnas, muros,
    diafragmas rigidos). Devuelve el objeto Edificio contenedor."""
    ed = Edificio(cfg)
    ed.construir()
    return ed


def configurar_y_analizar():
    """Configura el solucionador y ejecuta un paso estatico lineal."""
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError("El analisis estatico fallo (codigo {})".format(ok))
    ops.reactions()


def correr_caso_gravedad(cfg, ed, q_map, con_puntuales=False):
    """Aplica un caso de carga uniforme por piso (G o Q) mediante la
    misma geometria tributaria de Semana 2 (patron 1) y analiza.

    Devuelve (w_map, resumen, area_piso_lvl, area_x, area_y, w_vol,
    total_aplicado_kN).
    """
    w_map, resumen, area_piso_lvl, area_x, area_y, w_vol = \
        areas_tributarias(cfg, q_map)
    total = aplicar_cargas(cfg, ed, w_map, w_vol,
                           tag_patron=1, tag_series=1,
                           con_puntuales=con_puntuales)
    configurar_y_analizar()
    return w_map, resumen, area_piso_lvl, area_x, area_y, w_vol, total


def aplicar_sismo(cfg, ed, fuerzas_dir, escala=1.0,
                  tag_patron=2, tag_series=2):
    """Aplica fuerzas laterales en el CM de cada piso (patron 2).

    fuerzas_dir: {"EX": {nivel: F}, "EY": {nivel: F}} (F en kN).
    escala     : multiplicador (para la superposicion explicita).

    Por cada piso se aplica la fuerza horizontal F en la direccion
    solicitada y, como el CM no coincide con el nodo master del
    diafragma, el momento torsor equivalente (wrench) para reproducir
    la misma resultante:
        EX: (Fx, 0, 0) + Mz = -Fx * (y_cm - y_master)
        EY: (0, Fy, 0) + Mz =  Fy * (x_cm - x_master)
    """
    cm = geo.centros_de_masa(cfg)
    ops.timeSeries("Linear", tag_series)
    ops.pattern("Plain", tag_patron, tag_series)
    for d, fdict in fuerzas_dir.items():
        for lvl in ed.levels:
            f = fdict.get(lvl, 0.0) * escala
            if abs(f) < 1e-12:
                continue
            nm = ed.master_node[lvl]
            xm, ym, _ = ed.nodos[nm]
            _a, xc, yc = cm[lvl]
            dx_cm, dy_cm = xc - xm, yc - ym
            if d.upper() == "EX":
                ops.load(nm, f, 0.0, 0.0, 0.0, 0.0, -f * dy_cm)
            elif d.upper() == "EY":
                ops.load(nm, 0.0, f, 0.0, 0.0, 0.0, f * dx_cm)


def correr_caso_sismo(cfg, ed, direccion, fuerzas):
    """Aplica un caso sismico exclusivo (EX o EY) y analiza."""
    fuerzas_dir = {"EX": {}, "EY": {}}
    fuerzas_dir[direccion] = fuerzas
    aplicar_sismo(cfg, ed, fuerzas_dir)
    configurar_y_analizar()


# =====================================================================
# RESULTADOS
# =====================================================================

def seleccionar_elementos(ed):
    """Devuelve los tags de los elementos de control:
      columna_critica : columna del nivel Piso1 en el cruce del ultimo
                        eje X (I') y ultimo eje Y (Eje 1).
      viga_verifica   : primera viga X del nivel Piso2.
    """
    nx, ny = ed.nx, ed.ny
    col = None
    for c in ed.columnas:
        if (c["nivel"] == "Piso1" and c["ix"] == nx - 1
                and c["iy"] == ny - 1):
            col = c["elementTag"]
            break
    viga = None
    for tag, (n1, n2, tipo, ix, iy, k) in ed.coord_viga.items():
        if (tipo == "X" and ix == 0 and iy == 0 and k == 1):
            viga = tag
            break
    return {"columna_critica": col, "viga_verifica": viga}


class Resultado:
    """Resultados de un caso base: desplazamientos de los nodos master,
    reacciones basales y fuerzas internas basicas de los elementos de
    control."""

    def __init__(self, ed, tags_control):
        self.tags_control = tags_control
        self.desp = {lvl: np.array(ops.nodeDisp(ed.master_node[lvl]))
                     for lvl in ed.levels}
        self.reac = {}
        for nid, (x, y, z) in ed.nodos.items():
            if abs(z - ed.z_base) < 1e-9:
                self.reac[nid] = np.array(ops.nodeReaction(nid))
        self.fuerzas = {}
        for nombre, tag in tags_control.items():
            if tag is None:
                continue
            self.fuerzas[nombre] = np.array(ops.eleResponse(tag, "basicForce"))

    @property
    def reacciones_totales(self):
        if getattr(self, "_reac_tot", None) is not None:
            return self._reac_tot
        return sum(self.reac.values(), np.zeros(6))

    def master(self, lvl):
        return self.desp[lvl]


def registrar(cfg, ed, tags_control):
    return Resultado(ed, tags_control)


def verificar_conservacion_qv(cfg, qJ, resumen, area_piso_lvl,
                              total_aplicado, reac_fz, tol=1e-6):
    """Verificacion parte A: sum(Q transferida) = q_Q * A.

    qJ: {nivel: q_Q} intensidades de carga viva.
    """
    q_total_por_piso = {lvl: qJ[lvl] * area_piso_lvl[lvl]
                        for lvl in cfg["orden_niveles"]}
    q_total = sum(q_total_por_piso.values())
    carga_transferida = sum(r["carga_total"] for r in resumen.values())
    ok_area = abs(carga_transferida - q_total) / max(q_total, 1e-12) < tol
    ok_reacc = abs(reac_fz - q_total) / max(q_total, 1e-12) < tol
    return {
        "q_total_por_piso_kN": q_total_por_piso,
        "q_total_kN": q_total,
        "carga_transferida_kN": carga_transferida,
        "diferencia_relativa_transferida": abs(carga_transferida - q_total)
                                          / max(q_total, 1e-12),
        "suma_reacciones_FZ_kN": reac_fz,
        "diferencia_relativa_reacciones": abs(reac_fz - q_total)
                                         / max(q_total, 1e-12),
        "ok_conservacion_area": bool(ok_area),
        "ok_conservacion_reacciones": bool(ok_reacc),
    }


def verificar_sismo(cfg, casos, fuerzas, res, direccion):
    """Verificaciones parte B para un caso sismico EX o EY.

    - carga lateral total = suma de F_i
    - corte basal = suma de reacciones horizontales
    - sentido de la deformada (signo de desplazamiento del techo)
    - torsion de piso (rotacion Rz de cada piso con su excentricidad)
    """
    ws = geo.pesos_sismicos(cfg, casos)
    lateral = sum(fuerzas.values())
    if direccion == "EX":
        corte = sum(r[0] for r in res.reac.values())
        desp_techo = res.desp[cfg["orden_niveles"][-1]][0]
    else:
        corte = sum(r[1] for r in res.reac.values())
        desp_techo = res.desp[cfg["orden_niveles"][-1]][1]
    cm = geo.centros_de_masa(cfg)
    torsion = []
    for lvl in cfg["orden_niveles"]:
        rz = res.desp[lvl][5]
        _a, xc, yc = cm[lvl]
        torsion.append({"nivel": lvl, "Rz_rad": rz,
                        "CM_x": xc, "CM_y": yc})
    sentido_ok = (desp_techo >= 0) == (lateral >= 0) if abs(lateral) > 0 else True
    return {
        "direccion": direccion,
        "carga_lateral_total_kN": lateral,
        "corte_basal_kN": corte,
        "desplazamiento_techo_m": desp_techo,
        "sentido_deformada_correcto": bool(sentido_ok),
        "pesos_sismicos": ws,
        "torsion_por_piso": torsion,
        "deformada_m": {lvl: float(res.desp[lvl][0 if direccion == "EX" else 1])
                        for lvl in cfg["orden_niveles"]},
    }