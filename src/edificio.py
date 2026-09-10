"""
EDIFICIO INSTITUCIONAL - Analisis estatico de gravedad
=======================================================
Proyecto 1 - Grupo 4

Modelo 3D parametrico de un edificio de hormigon armado (5 niveles:
1 subterraneo + 4 pisos) con arriostramientos metalicos y muros de
corte, todo leido desde data/edificio_config.json (PROHIBIDO
hardcodear coordenadas/secciones/cargas en este script).

Requerimientos implementados:
  * Nodos en todas las intersecciones de grilla y niveles; apoyos
    empotrados (6 GDL) en la fundacion.
  * Columnas y vigas con elasticBeamColumn (analisis lineal elastico).
  * Muros de corte con metodo de "columna ancha" (wide column):
    columna equivalente en el centroide + brazos rigidos (rigidLink).
  * Diafragma rigido por nivel (rigidDiaphragm, constraints
    Transformation) amarrando GDL horizontales a un nodo master.
  * Arriostramientos metalicos como truss (solo axial).
  * Cargas de gravedad qG transferidas por AREAS TRIBUTARIAS a vigas
    via eleLoad -beamUniform (sin modelar la losa con shells).
  * Cargas puntuales excepcionales (equipos) sobre nodos o vigas.
  * Modulo QA/CSV: carga de losa por piso, suma de areas tributarias,
    reacciones basales, compatibilidad de diafragma.

Uso:
    .venv\\Scripts\\Activate.ps1
    python -m src.edificio
    python -m src.edificio data/edificio_config.json

NOTA: los resultados CSV se imprimen en consola (stdout) para poder
copiarlos y cruzarlos en una hoja de calculo.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")   # backend sin ventana: guarda la figura en archivo
import matplotlib.pyplot as plt
import numpy as np
import openseespy.opensees as ops

FIG_DIR = Path("results") / "figures"
FIG_EDIFICIO = "edificio_3d.png"
EXPORT_DIR = Path("results") / "export"


# =====================================================================
# 0. CARGAR CONFIGURACION
# =====================================================================

def cargar_config(ruta="data/edificio_config.json"):
    """Lee y valida el JSON externo con toda la parametrizacion."""
    with open(ruta, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    # Orden estable de ejes y niveles a partir de los dicts
    cfg["orden_x"] = list(cfg["grilla_ejes"]["X"].keys())
    cfg["orden_y"] = list(cfg["grilla_ejes"]["Y"].keys())
    cfg["orden_niveles"] = [k for k in cfg["niveles"].keys()
                            if k not in ("base", "nota")]
    return cfg


# =====================================================================
# 1. SECCIONES: FUNCIONES GENERICAS
# =====================================================================

def prop_rect(b, h, E, nu):
    """Seccion rectangular: area, inercias y torsion."""
    G = E / (2.0 * (1.0 + nu))
    A = b * h
    Iy = b * h ** 3 / 12.0
    Iz = h * b ** 3 / 12.0
    # Torsion aproximada de seccion rectangular
    a, c = max(b, h), min(b, h)
    J = 0.196 * (a * c ** 3) * (1.0 - 0.63 * c / a)
    return [A, E, G, J, Iy, Iz]


def dict_secciones(cfg):
    """Precomputa todas las secciones tipo del config.

    Devuelve {clase: {nombre: lista([A,E,G,J,Iy,Iz])}} y los
    materiales {nombre: {E,nu}}.
    """
    materiales = {}
    for nom, m in cfg["materiales"].items():
        materiales[nom] = m

    sec = {"vigas": {}, "columnas": {}, "muros": {}, "diagonales": {}}

    # Vigas
    for nom, s in cfg["secciones_tipo"]["vigas"].items():
        E = materiales[s["material"]]["E"]
        nu = materiales[s["material"]]["nu"]
        sec["vigas"][nom] = prop_rect(s["b"], s["h"], E, nu)

    # Columnas
    for nom, s in cfg["secciones_tipo"]["columnas"].items():
        E = materiales[s["material"]]["E"]
        nu = materiales[s["material"]]["nu"]
        sec["columnas"][nom] = prop_rect(s["b"], s["h"], E, nu)

    # Muros (columna ancha): seccion rectangular largo x espesor
    for nom, s in cfg["secciones_tipo"]["muros"].items():
        E = materiales[s["material"]]["E"]
        nu = materiales[s["material"]]["nu"]
        # largo L y espesor t se fijan al crear cada muro
        sec["muros"][nom] = {"E": E, "nu": nu, "t": s["t"]}

    # Diagonales (truss): solo area axial
    for nom, s in cfg["secciones_tipo"]["diagonales"].items():
        E = materiales[s["material"]]["E"]
        sec["diagonales"][nom] = {"A": s["A"], "E": E,
                                  "material": s["material"]}

    return sec, materiales


# =====================================================================
# 2. MODELO: NODOS, APOYOS, TRANSFORMACIONES Y ELEMENTOS
# =====================================================================

class Edificio:
    """Contenedor del modelo que acumula nodos/elementos y metadatos."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.sec, self.mat = dict_secciones(cfg)
        # Registro de nodos: id -> (x, y, z)
        self.nodos = {}
        self._next = 1
        # Nodos de la grilla por (ix, iy, k): id
        self.grid_node = {}
        self.nx = len(cfg["orden_x"])
        self.ny = len(cfg["orden_y"])
        self.levels = cfg["orden_niveles"]
        self.z_level = {k: cfg["niveles"][k] for k in self.levels}
        self.z_base = cfg["niveles"]["base"]
        # Contador GLOBAL de tags de elementos (evita colisiones)
        self._etag = 0
        # Estadisticas
        self.n_beams_by_level = {}
        self.area_tributaria_x = 0.0
        self.area_tributaria_y = 0.0
        self.master_node = {}
        self._transf = {}
        # Registros para exportar a Unity / viewer
        #   apoyos: {nid: [u1,u2,u3,r1,r2,r3]}  (1=restringido, 0=libre)
        self.apoyos = {}
        #   columnas: lista de dicts (tag, n1, n2, seccion, material, transf, L)
        self.columnas = []
        #   diagonales: lista de dicts (tag, n1, n2, seccion, material, tipo='truss', L)
        self.diagonales = []
        #   voladizo (eje J): vigas cantilever I'->J por nivel
        #   vol_vigas: {lvl: [(tag, n1, n2, L, iy)]}
        #   vol_coord:  {tag: (n1, n2, iy, k)}
        #   vol_saliente: {lvl: ancho_saliente_m}
        self.vol_vigas = {}
        self.vol_coord = {}
        self.vol_saliente = {}
        self.n_vol = 0

    def _nelem(self):
        """Devuelve un tag de elemento unico y global."""
        self._etag += 1
        return self._etag

    # -- helpers ------------------------------------------------------
    def _add_node(self, x, y, z):
        nid = self._next
        self._next += 1
        self.nodos[nid] = (x, y, z)
        return nid

    def _coord(self, name_x=None, name_y=None):
        """Convierte nombres de ejes a coordenadas X o Y."""
        x = self.cfg["grilla_ejes"]["X"][name_x] if name_x else None
        y = self.cfg["grilla_ejes"]["Y"][name_y] if name_y else None
        return x, y

    def _xcoords(self):
        return [self.cfg["grilla_ejes"]["X"][k] for k in self.cfg["orden_x"]]

    def _ycoords(self):
        return [self.cfg["grilla_ejes"]["Y"][k] for k in self.cfg["orden_y"]]

    def construir(self):
        """Construye todo el modelo (nodos, apoyos, elementos, muros,
        arriostramientos)."""
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        self._crear_materiales()
        self._crear_nodos()
        self._crear_transformaciones()
        self._crear_columnas()
        self._crear_vigas()
        self._crear_voladizos()
        self._crear_muros()
        self._crear_diafragmas()
        self._crear_arriostramientos()

    # -- nodos y apoyos ----------------------------------------------
    def _crear_nodos(self):
        """Nodo en cada interseccion de grilla y en cada nivel (+ base)."""
        xc = self._xcoords()
        yc = self._ycoords()
        # Base (fundacion)
        for iy in range(self.ny):
            for ix in range(self.nx):
                nid = self._add_node(xc[ix], yc[iy], self.z_base)
                self.grid_node[(ix, iy, -1)] = nid
                # Apoyo empotrado (6 GDL fijos)
                ops.node(nid, xc[ix], yc[iy], self.z_base)
                ops.fix(nid, 1, 1, 1, 1, 1, 1)
                self.apoyos[nid] = [1, 1, 1, 1, 1, 1]
        # Niveles
        for k, lvl in enumerate(self.levels):
            z = self.z_level[lvl]
            for iy in range(self.ny):
                for ix in range(self.nx):
                    nid = self._add_node(xc[ix], yc[iy], z)
                    self.grid_node[(ix, iy, k)] = nid
                    ops.node(nid, xc[ix], yc[iy], z)

    # -- transformaciones --------------------------------------------
    def _crear_materiales(self):
        """Define uniaxialMaterial Elastic para cada material del
        config. elasticBeamColumn y truss los referencian por tag."""
        self.mat_tag = {}
        for idx, (nom, m) in enumerate(self.mat.items(), start=1):
            ops.uniaxialMaterial("Elastic", idx, m["E"])
            self.mat_tag[nom] = idx

    def _crear_transformaciones(self):
        """Vecxz para columnas, vigas en X y vigas en Y.

        NOTA: ops.geomTransf no devuelve el tag (retorna None); el tag
        se fija como primer argumento, por eso se guardan literales."""
        ops.geomTransf("Linear", 1, 1, 0, 0)   # columnas (eje local x vertical)
        ops.geomTransf("Linear", 2, 0, 0, 1)   # vigas en X
        ops.geomTransf("Linear", 3, 0, 0, 1)   # vigas en Y
        self._transf["col"] = 1
        self._transf["viga_x"] = 2
        self._transf["viga_y"] = 3

    # -- columnas ----------------------------------------------------
    def _crear_columnas(self):
        """Columnas entre nivel y siguiente (la base incluye el tramo
        hacia el primer nivel, es decir el subterraneo)."""
        secciones = self.cfg["asignaciones"]["columnas_por_nivel"]
        n_cols = 0
        for iy in range(self.ny):
            for ix in range(self.nx):
                # Tramo base -> nivel 0 (subterraneo)
                n_prev = self.grid_node[(ix, iy, -1)]
                z_prev = self.z_base
                for k, lvl in enumerate(self.levels):
                    n_cur = self.grid_node[(ix, iy, k)]
                    z_cur = self.z_level[lvl]
                    seccion = secciones[lvl]
                    params = self.sec["columnas"][seccion]
                    material = self.cfg["secciones_tipo"]["columnas"][seccion]["material"]
                    tag = self._nelem()
                    n_cols += 1
                    ops.element("elasticBeamColumn", tag, n_prev, n_cur,
                                *params, self._transf["col"])
                    self.columnas.append({
                        "elementTag": tag, "tipo": "columna",
                        "n1": n_prev, "n2": n_cur,
                        "nivel": lvl, "seccion": seccion, "material": material,
                        "transf": self._transf["col"],
                        "L": float(z_cur - z_prev),
                        "ix": ix, "iy": iy, "k": k})
                    n_prev, z_prev = n_cur, z_cur
        self.n_cols = n_cols

    # -- vigas -------------------------------------------------------
    def _crear_vigas(self):
        """Vigas en la grilla de cada nivel. Guarda los tags de vigas
        por nivel para asignarles las cargas tributarias."""
        secciones = self.cfg["asignaciones"]["vigas_por_nivel"]
        xc = self._xcoords()
        yc = self._ycoords()
        self.vigas = {}          # lvl -> lista de (tag, n1, n2, L, tipo)
        self.coord_viga = {}     # tag -> (n1,n2,tipo,ix,iy,k)
        n_vigas = 0
        for k, lvl in enumerate(self.levels):
            seccion = secciones[lvl]
            params = self.sec["vigas"][seccion]
            lista = []
            # Vigas en X (paralelas a X, sobre cada eje Y)
            for iy in range(self.ny):
                for ix in range(self.nx - 1):
                    n1 = self.grid_node[(ix, iy, k)]
                    n2 = self.grid_node[(ix + 1, iy, k)]
                    tag = self._nelem()
                    n_vigas += 1
                    ops.element("elasticBeamColumn", tag, n1, n2,
                                *params, self._transf["viga_x"])
                    L = xc[ix + 1] - xc[ix]
                    lista.append((tag, n1, n2, L, "X"))
                    self.coord_viga[tag] = (n1, n2, "X", ix, iy, k)
            # Vigas en Y (paralelas a Y, sobre cada eje X)
            for ix in range(self.nx):
                for iy in range(self.ny - 1):
                    n1 = self.grid_node[(ix, iy, k)]
                    n2 = self.grid_node[(ix, iy + 1, k)]
                    tag = self._nelem()
                    n_vigas += 1
                    ops.element("elasticBeamColumn", tag, n1, n2,
                                *params, self._transf["viga_y"])
                    L = yc[iy + 1] - yc[iy]
                    lista.append((tag, n1, n2, L, "Y"))
                    self.coord_viga[tag] = (n1, n2, "Y", ix, iy, k)
            self.vigas[lvl] = lista
        self.n_vigas = n_vigas

    # -- voladizo (eje J) --------------------------------------------
    def _crear_voladizos(self):
        """Vigas cantilever del voladizo (eje J) en los niveles del
        config.

        Para cada voladizo se crea una viga cantilever en cada eje Y
        indicado, desde el nodo de la grilla en 'desde_eje' (I') hasta
        un nodo NUEVO en x_j_m (eje J), sin columnas bajo J (la losa
        cuelga de sus vigas). Se guardan aparte de coord_viga."""
        cfg = self.cfg
        yc = self._ycoords()
        n_vol = 0
        for v in cfg.get("voladizos", {}).get("lista", []):
            x0 = cfg["grilla_ejes"]["X"][v["desde_eje"]]
            x1 = v["x_j_m"]
            sal = v.get("ancho_saliente_m", 0.0)
            seccion = v["seccion_viga"]
            params = self.sec["vigas"][seccion]
            for lvl in v["niveles"]:
                k = self.levels.index(lvl)
                z = self.z_level[lvl]
                self.vol_saliente.setdefault(lvl, sal)
                for eje_y in v["ejes_Y"]:
                    iy = cfg["orden_y"].index(eje_y)
                    y = yc[iy]
                    # nodo raiz: en la grilla (eje 'desde_eje'), nivel k
                    n0 = self.grid_node[(cfg["orden_x"].index(v["desde_eje"]),
                                         iy, k)]
                    n1 = self._add_node(x1, y, z)
                    ops.node(n1, x1, y, z)
                    tag = self._nelem()
                    n_vol += 1
                    ops.element("elasticBeamColumn", tag, n0, n1,
                                *params, self._transf["viga_x"])
                    L = float(x1 - x0)
                    self.vol_vigas.setdefault(lvl, []).append(
                        (tag, n0, n1, L, iy))
                    self.vol_coord[tag] = (n0, n1, iy, k)
        self.n_vol = n_vol

    # -- muros (wide column) -----------------------------------------
    def _crear_muros(self):
        """Cada muro se modela como una columna ancha: columna
        equivalente en el centroide + rigidLink (brazo rigido) hacia
        los extremos en cada nivel."""
        xc = self._xcoords()
        yc = self._ycoords()
        X = self.cfg["grilla_ejes"]["X"]
        Y = self.cfg["grilla_ejes"]["Y"]
        self.muro_info = []
        nid_mur = 0
        for m in self.cfg["muros"]["lista"]:
            orient = m["orientacion"]
            nom_sec = m["seccion"]
            sec = self.sec["muros"][nom_sec]
            t = sec["t"]
            E = sec["E"]
            nu = sec["nu"]
            a = Y[m["linea"]] if orient == "X" else X[m["linea"]]
            # Recorrido de ejes segun orientacion
            if orient == "X":
                ejes = self.cfg["orden_x"]
                coords = xc
                idx0, idx1 = ejes.index(m["desde_eje"]), ejes.index(m["hasta_eje"])
                # linea es un eje Y constante
                fixed_y = a
                x0, x1 = coords[idx0], coords[idx1]
                L = abs(x1 - x0)
                centr = (x0 + x1) / 2.0
                fixed_x = centr
                # indices extremos en la grilla
                ex0, ex1 = idx0, idx1
            else:
                ejes = self.cfg["orden_y"]
                coords = yc
                idx0, idx1 = ejes.index(m["desde_eje"]), ejes.index(m["hasta_eje"])
                fixed_x = a
                y0, y1 = coords[idx0], coords[idx1]
                L = abs(y1 - y0)
                centr = (y0 + y1) / 2.0
                fixed_y = centr
                ex0, ex1 = idx0, idx1

            # seccion de la columna ancha (largo L x espesor t)
            A = L * t
            Iy = t * L ** 3 / 12.0
            Iz = L * t ** 3 / 12.0
            G = E / (2 * (1 + nu))
            a_, c_ = max(L, t), min(L, t)
            J = 0.196 * (a_ * c_ ** 3) * (1 - 0.63 * c_ / a_)
            params = [A, E, G, J, Iy, Iz]

            # indices de nivel
            ks = self.levels.index(m["desde_nivel"])
            ke = self.levels.index(m["hasta_nivel"])

            # Crear nodo de columna ancha por nivel (+ base)
            n_wall_prev = None
            z_prev = self.z_base
            for k in range(ks, ke + 1):
                z = self.z_level[self.levels[k]]
                if k == ks:
                    if ks == 0:
                        nbase = self._add_node(centr, fixed_y, self.z_base) \
                            if orient == "X" else self._add_node(fixed_x, centr, self.z_base)
                        ops.node(nbase, self._cx(orient, centr, fixed_y, fixed_x),
                                 self._cy(orient, centr, fixed_y, fixed_x), self.z_base)
                        ops.fix(nbase, 1, 1, 1, 1, 1, 1)
                        n_wall_prev = nbase
                        z_prev = self.z_base
                nw = self._add_node(self._cx(orient, centr, fixed_y, fixed_x),
                                    self._cy(orient, centr, fixed_y, fixed_x), z)
                ops.node(nw, self._cx(orient, centr, fixed_y, fixed_x),
                         self._cy(orient, centr, fixed_y, fixed_x), z)
                # elemento columna ancha
                if n_wall_prev:
                    tag = self._nelem()
                    nid_mur += 1
                    ops.element("elasticBeamColumn", tag,
                                n_wall_prev, nw, *params, self._transf["col"])
                # brazos rigidos a los extremos de la grilla en este nivel
                if orient == "X":
                    nendo1 = self.grid_node[(ex0, self.cfg["orden_y"].index(m["linea"]), k)]
                    nendo2 = self.grid_node[(ex1, self.cfg["orden_y"].index(m["linea"]), k)]
                else:
                    iyline = self.cfg["orden_x"].index(m["linea"])
                    nendo1 = self.grid_node[(iyline, ex0, k)]
                    nendo2 = self.grid_node[(iyline, ex1, k)]
                ops.rigidLink("beam", nw, nendo1)
                ops.rigidLink("beam", nw, nendo2)
                n_wall_prev = nw
                z_prev = z
            self.muro_info.append(m)

    @staticmethod
    def _cx(orient, centr, fixed_y, fixed_x):
        return centr if orient == "X" else fixed_x

    @staticmethod
    def _cy(orient, centr, fixed_y, fixed_x):
        return fixed_y if orient == "X" else centr

    # -- diafragmas rigidos ------------------------------------------
    def _crear_diafragmas(self):
        """Un diafragma rigido por nivel que amarra los GDL en planta
        (Ux, Uy, Rz) de los nodos de la losa a un nodo master.

        Se usa como master una COLUMNA CENTRAL (ya conectada a las
        columnas) para no crear nodos flotantes sin rigidez vertical."""
        self.master_node = {}
        self.diafragmas = {}
        # indice del nodo central de la grilla
        ixc = self.nx // 2
        iyc = self.ny // 2
        for k, lvl in enumerate(self.levels):
            nm = self.grid_node[(ixc, iyc, k)]
            self.master_node[lvl] = nm
            slaves = []
            for iy in range(self.ny):
                for ix in range(self.nx):
                    if (ix == ixc and iy == iyc):
                        continue
                    if (ix, iy, k) in self.grid_node:
                        slaves.append(self.grid_node[(ix, iy, k)])
            self.diafragmas[lvl] = {"master": nm, "slaves": list(slaves)}
            # Diafragma rigido perpendicular a Z (dir 3)
            ops.rigidDiaphragm(3, nm, *slaves)

    def _crear_arriostramientos(self):
        """Diagonales de acero (truss, solo axial) en los porticos
        exteriores indicados por el config, patron X."""
        ar = self.cfg["arriostramientos"]
        sec = self.sec["diagonales"][ar["seccion"]]
        A = sec["A"]
        E = sec["E"]
        mat_tag = self.mat_tag[sec["material"]]
        ks = self.levels.index(ar["desde_nivel"])
        ke = self.levels.index(ar["hasta_nivel"])
        xc = self._xcoords()
        yc = self._ycoords()
        n_diag = 0
        for fy in ar["frames_Y"]:
            iy = self.cfg["orden_y"].index(fy)
            for ix in range(self.nx - 1):
                for k in range(ks, ke + 1):
                    n_bot1 = self.grid_node[(ix, iy, k - 1)]
                    n_bot2 = self.grid_node[(ix + 1, iy, k - 1)]
                    n_top1 = self.grid_node[(ix, iy, k)]
                    n_top2 = self.grid_node[(ix + 1, iy, k)]
                    lvl = self.levels[k]
                    for diag in ((n_bot1, n_top2), (n_bot2, n_top1)):
                        tag = self._nelem()
                        n_diag += 1
                        ops.element("truss", tag, diag[0], diag[1], A, mat_tag)
                        c0 = np.array(self.nodos[diag[0]])
                        c1 = np.array(self.nodos[diag[1]])
                        self.diagonales.append({
                            "elementTag": tag, "tipo": "diagonal",
                            "n1": diag[0], "n2": diag[1],
                            "nivel": lvl, "seccion": ar["seccion"],
                            "material": sec["material"],
                            "transf": None, "patron": ar["patron"],
                            "frame": fy, "ix": ix,
                            "L": float(np.linalg.norm(c1 - c0))})
        self.n_diag = n_diag


# =====================================================================
# 3. AREAS TRIBUTARIAS Y CARGAS
# =====================================================================

def areas_tributarias(cfg, q_map=None):
    """Rutina algoritmica explicita de areas tributarias.

    Metodo simplificado de losas en dos direcciones:
      * El area de cada panel se reparte 50% a la direccion X y 50% a la
        direccion Y.
      * El 50% de cada direccion se divide a partes iguales entre las
        dos vigas paralelas (superior e inferior) que limitan el panel.

    El voladizo (eje J) se agrega como una franja de losa adicional en
    los niveles indicados: su area (desde_eje hasta x_j+saliente,
    sobre todo el ancho de ejes Y) se reparte 100% a la direccion X
    (vigas cantilever) para que el balance de areas y de carga siga
    siendo exacto por nivel.

    Parametros:
      q_map: {nivel: intensidad_kN_m2} opcional. Si no se entrega usa
             cfg["cargas"]["qG"] (caso G). Con la misma geometria
             tributaria se puede aplicar cualquier caso de carga
             uniforme por piso (p.ej. el caso Q con q_Q).

    Devuelve:
      w_map     {nivel: {tag_viga: w_kN/m}}   (solo grilla base)
      resumen   QA por nivel (agrega 'area_piso', 'carga_voladizo',
                'area_voladizo')
      area_piso_lvl  {nivel: m2} con el voladizo incluido
      area_x, area_y totales globales (m2)
      w_vol     {nivel: {iy: w_kN/m}}  carga en las vigas cantilever
    """
    xc = [cfg["grilla_ejes"]["X"][k] for k in cfg["orden_x"]]
    yc = [cfg["grilla_ejes"]["Y"][k] for k in cfg["orden_y"]]
    dx = [xc[i + 1] - xc[i] for i in range(len(xc) - 1)]
    dy = [yc[j + 1] - yc[j] for j in range(len(yc) - 1)]
    if q_map is None:
        qG = {k: v for k, v in cfg["cargas"]["qG"].items() if k != "nota"}
    else:
        qG = q_map
    area_piso_main = sum(dx) * sum(dy)
    area_piso_lvl = {lvl: area_piso_main for lvl in cfg["orden_niveles"]}

    # identificador de viga: (tipo, iX, iY, nivel) -> w acumulado
    w = {}
    area_x = 0.0
    area_y = 0.0
    resumen = {}
    for k, lvl in enumerate(cfg["orden_niveles"]):
        q = qG[lvl]
        area_x_lvl = 0.0
        area_y_lvl = 0.0
        for i in range(len(dx)):          # paneles en X
            for j in range(len(dy)):      # paneles en Y
                panel = dx[i] * dy[j]
                # 50% a direccion X -> repartido entre 2 vigas X
                w_x = q * (panel / 2.0) / 2.0 / dx[i]   # = q * dy[j] / 4
                # 50% a direccion Y -> entre 2 vigas Y
                w_y = q * (panel / 2.0) / 2.0 / dy[j]   # = q * dx[i] / 4
                # vigas X del panel: inferior (j) y superior (j+1)
                for jj in (j, j + 1):
                    key = ("X", i, jj, lvl)
                    w[key] = w.get(key, 0.0) + w_x
                    area_x += w_x * dx[i] / q if q else 0.0
                    area_x_lvl += w_x * dx[i] / q if q else 0.0
                # vigas Y del panel: izquierda (i) y derecha (i+1)
                for ii in (i, i + 1):
                    key = ("Y", ii, j, lvl)
                    w[key] = w.get(key, 0.0) + w_y
                    area_y += w_y * dy[j] / q if q else 0.0
                    area_y_lvl += w_y * dy[j] / q if q else 0.0
        total = sum(w[(t, i, j, lvl)] * (dx[i] if t == "X" else dy[j])
                    for (t, i, j, ll) in w if ll == lvl)
        resumen[lvl] = {"qG": qG[lvl], "carga_total": total,
                        "area_trib_X": area_x_lvl, "area_trib_Y": area_y_lvl,
                        "area_piso": area_piso_lvl[lvl],
                        "carga_voladizo": 0.0, "area_voladizo": 0.0,
                        "suma_w_x": sum(w[(t, i, j, lvl)] for (t, i, j, ll) in w
                                        if ll == lvl and t == "X"),
                        "suma_w_y": sum(w[(t, i, j, lvl)] for (t, i, j, ll) in w
                                        if ll == lvl and t == "Y")}

    # --- Voladizo (eje J): franja de losa cantilever ---------------
    w_vol = {}
    for v in cfg.get("voladizos", {}).get("lista", []):
        x0 = cfg["grilla_ejes"]["X"][v["desde_eje"]]
        x1 = v["x_j_m"]
        sal = v.get("ancho_saliente_m", 0.0)
        width_x = (x1 + sal) - x0
        yv = [cfg["grilla_ejes"]["Y"][e] for e in v["ejes_Y"]]
        hgt = max(yv) - min(yv)
        A_v = width_x * hgt
        L_total_vigas = (x1 - x0) * len(v["ejes_Y"])
        for lvl in v["niveles"]:
            q = qG[lvl]
            area_piso_lvl[lvl] += A_v
            area_x += A_v
            resumen[lvl]["area_piso"] = area_piso_lvl[lvl]
            resumen[lvl]["carga_total"] += q * A_v
            resumen[lvl]["area_trib_X"] += A_v   # 100% carga a vigas cantilever
            resumen[lvl]["carga_voladizo"] += q * A_v
            resumen[lvl]["area_voladizo"] += A_v
            w_m = q * A_v / L_total_vigas   # distribuida en las vigas J
            w_vol[lvl] = {cfg["orden_y"].index(e): w_m for e in v["ejes_Y"]}
    return w, resumen, area_piso_lvl, area_x, area_y, w_vol


def aplicar_cargas(cfg, edificio, w_map, w_vol=None, tag_patron=1,
                   tag_series=1, con_puntuales=True):
    """Aplica las cargas tributarias (eleLoad) y las puntuales.

    w_vol: {nivel: {iy: w_kN/m}} con la carga de las vigas cantilever
    del voladizo.

    Parametros extras (para construir casos independientes en el mismo
    modelo o en modelos reconstruidos):
      tag_patron/tag_series: etiquetas del load pattern y la serie.
      con_puntuales: False para casos donde no aplican las puntuales
      (p.ej. caso de carga viva Q).

    Devuelve la carga vertical total aplicada (suma de las cargas de
    losa sobre vigas + puntuales) para la verificacion de conservacion.
    """
    ops.timeSeries("Linear", tag_series)
    ops.pattern("Plain", tag_patron, tag_series)
    total = 0.0
    # Cargas tributarias sobre vigas
    for (tipo, i, j, lvl), w_val in w_map.items():
        # hallar el tag de viga
        tag = _tag_viga(edificio, tipo, i, j, lvl)
        # en coordenadas locales: carga vertical en z=-w (hacia abajo)
        ops.eleLoad("-ele", tag, "-type", "-beamUniform", 0.0, -w_val)
        # longitud del tramo para acumular carga aplicada
        n1, n2 = edificio.coord_viga[tag][0], edificio.coord_viga[tag][1]
        p1 = np.array(edificio.nodos[n1])
        p2 = np.array(edificio.nodos[n2])
        L = float(np.linalg.norm(p2 - p1))
        total += w_val * L
    # Cargas del voladizo (vigas cantilever I'->J)
    if w_vol:
        for lvl, beams in edificio.vol_vigas.items():
            w_lvl = w_vol.get(lvl, {})
            for (tag, n1, n2, L, iy) in beams:
                w_val = w_lvl.get(iy, 0.0)
                if w_val:
                    ops.eleLoad("-ele", tag, "-type",
                                "-beamUniform", 0.0, -w_val)
                    total += w_val * L
    # Cargas puntuales
    if con_puntuales:
        for p in cfg["cargas"]["puntuales"]:
            peso_kN = p["peso_kg"] * cfg["cargas"]["g"] / 1000.0 * p["factor"]
            if p["tipo"] == "nodo":
                nid = edificio.grid_node[
                    (cfg["orden_x"].index(p["eje_x"]),
                     cfg["orden_y"].index(p["eje_y"]),
                     edificio.levels.index(p["nivel"]))]
                ops.load(nid, 0.0, 0.0, -peso_kN, 0.0, 0.0, 0.0)
                total += peso_kN
            elif p["tipo"] == "viga":
                # aplicar como vector nodal consistente en los extremos de la viga
                total += _aplicar_puntual_viga(edificio, cfg, p)
    return total


def _tag_viga(edificio, tipo, i, j, lvl):
    """Devuelve el tag de la viga (tipo X/Y, indices i,j, nivel)."""
    for (tag, n1, n2, L, t) in edificio.vigas[lvl]:
        if t != tipo:
            continue
        # verificar por coordenadas/nodos
        if tipo == "X":
            if _indices_ok(edificio, tag, i, j, tipo):
                return tag
        else:
            if _indices_ok(edificio, tag, i, j, tipo):
                return tag
    raise KeyError(f"No se hallo viga {tipo} ({i},{j}) en {lvl}")


def _indices_ok(edificio, tag, i, j, tipo):
    """Compara nodos de la viga con la grilla para validar indices."""
    ni, nj = edificio.coord_viga[tag][0], edificio.coord_viga[tag][1]
    cfg = edificio.cfg
    if tipo == "X":
        # viga X: nodos en (ix, iy, k) y (ix+1, iy, k)
        return edificio.coord_viga[tag][3] == i and edificio.coord_viga[tag][4] == j
    else:
        return edificio.coord_viga[tag][3] == i and edificio.coord_viga[tag][4] == j


def _aplicar_puntual_viga(edificio, cfg, p):
    """Carga puntual en una viga: vector nodal consistente en extremos.

    Para una carga vertical P a distancia a del nodo i (viga de largo L),
    el vector de extremo fijo (consistent load) equivale a:
        V_i = P*(L-a)^2*(L+2a)/L^3
        V_j = P*a^2*(3L-2a)/L^3
        M_i = P*a*(L-a)^2/L^2
        M_j = -P*a^2*(L-a)/L^2
    (signos en coordenadas locales de OpenSees ajustados por el eje)."""
    P = p["peso_kg"] * cfg["cargas"]["g"] / 1000.0 * p["factor"]
    nivel = p["nivel"]
    tipo_deseado = p["direccion"]  # "X" o "Y"
    # localizar viga por nivel + eje
    tag, ni, nj = None, None, None
    for (t, n1, n2, L, tt) in edificio.vigas[nivel]:
        if tt == tipo_deseado and _es_viga_de(edificio, t, p):
            tag, ni, nj, L = t, n1, n2, L
            break
    if tag is None:
        print(f"[aviso] no se hallo viga para puntual: {p['descripcion']}")
        return 0.0
    a = L * p["posicion"]
    # nodos locales de la viga
    id_ini, id_fin = ni, nj
    # vector consistente (en global vertical)
    Vi = P * (L - a) ** 2 * (L + 2 * a) / L ** 3
    Vj = P * a ** 2 * (3 * L - 2 * a) / L ** 3
    Mi = P * a * (L - a) ** 2 / L ** 2
    Mj = -P * a ** 2 * (L - a) / L ** 2
    ops.load(id_ini, 0, 0, -Vi, 0, 0, 0)
    ops.load(id_fin, 0, 0, -Vj, 0, 0, 0)
    return P


def _es_viga_de(edificio, tag, p):
    """True si la viga 'tag' esta en el nivel y eje del puntual."""
    n1, n2, tipo, ix, iy, kk = edificio.coord_viga[tag]
    lvl = edificio.cfg["orden_niveles"][kk]
    if lvl != p["nivel"]:
        return False
    if tipo != p["direccion"]:
        return False
    # chequear que la viga yace sobre el eje_x o eje_y indicado
    if tipo == "X":
        return edificio.cfg["orden_y"][iy] == p["eje_y"]
    else:
        return edificio.cfg["orden_x"][ix] == p["eje_x"]


# =====================================================================
# 4. ANALISIS
# =====================================================================

def analizar():
    """Analisis estatico lineal de gravedad."""
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError("El analisis fallo (codigo {})".format(ok))
    ops.reactions()


# =====================================================================
# 5. VERIFICACION QA Y EXPORTACION CSV
# =====================================================================

def qa_verificacion(cfg, edificio, resumen, area_piso_lvl, area_x, area_y,
                    carga_aplicada):
    """Imprime en consola las verificaciones en formato CSV."""
    S = "=" * 78
    print("\n" + S)
    print("MODULO QA - VERIFICACION DE EQUILIBRIO GLOBAL (CSV)")
    print(S)

    # --- Carga de losa por piso y areas tributarias -----------------
    print("\n# nivel,qG_kN_m2,carga_losa_piso_kN,area_piso_m2,"
          "suma_w_vigasX_kN_m,suma_w_vigasY_kN_m")
    tot_losa = 0.0
    for lvl in cfg["orden_niveles"]:
        r = resumen[lvl]
        tot_losa += r["carga_total"]
        area = r["area_piso"]
        print(f"{lvl},{r['qG']:.3f},{r['carga_total']:.3f},{area:.3f},"
              f"{r['suma_w_x']:.3f},{r['suma_w_y']:.3f}")
    print(f"TOTAL_LOZAS,{sum(cfg['cargas']['qG'][l] for l in cfg['orden_niveles']):.3f},"
          f"{tot_losa:.3f},,,")

    print("\n# verificacion_areas_tributarias")
    print("area_piso_losa_teorica_m2,{:.6f}".format(
        resumen[cfg["orden_niveles"][0]]["area_piso"]))
    print("# nivel,area_tributaria_dirX_m2,area_tributaria_dirY_m2,total_m2,"
          "igual_a_area_piso")
    for lvl in cfg["orden_niveles"]:
        r = resumen[lvl]
        tot = r["area_trib_X"] + r["area_trib_Y"]
        ok = "SI" if abs(tot - r["area_piso"]) < 1e-6 else "NO"
        print(f"{lvl},{r['area_trib_X']:.6f},{r['area_trib_Y']:.6f},{tot:.6f},{ok}")
    print(f"TOTAL_DIRX_TODOS_PISOS_m2,{area_x:.6f}")
    print(f"TOTAL_DIRY_TODOS_PISOS_m2,{area_y:.6f}")
    print(f"niveles,{len(cfg['orden_niveles'])}")

    # --- Voladizo (eje J) -------------------------------------------
    has_vol = any(r["area_voladizo"] > 0 for r in resumen.values())
    if has_vol:
        print("\n# voladizo_eje_J")
        print("nivel,area_losa_voladizo_m2,carga_losa_voladizo_kN,"
              "porcentaje_area_piso_pct")
        for lvl in cfg["orden_niveles"]:
            r = resumen[lvl]
            if r["area_voladizo"] <= 0:
                continue
            pct = 100.0 * r["area_voladizo"] / r["area_piso"]
            print(f"{lvl},{r['area_voladizo']:.3f},"
                  f"{r['carga_voladizo']:.3f},{pct:.2f}")

    # --- Reacciones basales -----------------------------------------
    reacs = {}
    for nid, (x, y, z) in edificio.nodos.items():
        if abs(z - edificio.z_base) < 1e-9:
            r = ops.nodeReaction(nid)
            reacs[nid] = r
    sum_fx = sum(r[0] for r in reacs.values())
    sum_fy = sum(r[1] for r in reacs.values())
    sum_fz = sum(r[2] for r in reacs.values())
    print("\n# reacciones_basales")
    print("nodo,x_m,y_m,fx_kN,fy_kN,fz_kN,mx_kN_m,my_kN_m,mz_kN_m")
    for nid, r in reacs.items():
        x, y, z = edificio.nodos[nid]
        print(f"{nid},{x:.3f},{y:.3f},{r[0]:.4f},{r[1]:.4f},{r[2]:.4f},"
              f"{r[3]:.4f},{r[4]:.4f},{r[5]:.4f}")

    # --- Conservacion de carga --------------------------------------
    print("\n# conservacion_carga")
    print(f"carga_total_aplicada_kN,{carga_aplicada:.6f}")
    print(f"suma_reacciones_FZ_kN,{sum_fz:.6f}")
    print(f"diferencia_kN,{abs(carga_aplicada) - abs(sum_fz):.6e}")
    print(f"suma_reacciones_FX_kN,{sum_fx:.6f}")
    print(f"suma_reacciones_FY_kN,{sum_fy:.6f}")

    # --- Compatibilidad de diafragma --------------------------------
    print("\n# compatibilidad_diafragma (residuo rigid-body vs master)")
    print("nivel,ux_master_m,uy_master_m,rz_master_rad,max_residuo_Ux_m,max_residuo_Uy_m")
    for k, lvl in enumerate(edificio.levels):
        nm = edificio.master_node[lvl]
        d = ops.nodeDisp(nm)
        uxm, uym, rz = d[0], d[1], d[5]
        xm, ym, _ = edificio.nodos[nm]
        resid_x = 0.0
        resid_y = 0.0
        for iy in range(edificio.ny):
            for ix in range(edificio.nx):
                nid = edificio.grid_node[(ix, iy, k)]
                x, y, _ = edificio.nodos[nid]
                dd = ops.nodeDisp(nid)
                # desplazamiento esperado por cuerpo rigido
                ex = uxm - rz * (y - ym)
                ey = uym + rz * (x - xm)
                resid_x = max(resid_x, abs(dd[0] - ex))
                resid_y = max(resid_y, abs(dd[1] - ey))
        print(f"{lvl},{uxm:.6e},{uym:.6e},{rz:.6e},{resid_x:.6e},{resid_y:.6e}")

    print("\n" + S)
    print("FIN QA")
    print(S + "\n")


# =====================================================================
# 6. VISUALIZACION 3D
# =====================================================================

def _coord(edificio, nid):
    """Devuelve las coordenadas (x,y,z) de un nodo, con la deformada
    superpuesta (para el dibujo) si el modelo ya fue analizado y se pide
    escala. Sin deformada si edificio._vis_escala es None."""
    p = np.array(edificio.nodos[nid], float)
    if getattr(edificio, "_vis_escala", None):
        d = np.array(ops.nodeDisp(nid)[:3])
        return p + edificio._vis_escala * d
    return p


def visualizar(edificio, cfg, ruta=None):
    """Genera una vista 3D del edificio (columnas, vigas, diagonales y
    muro) con la deformada por gravedad amplificada y la guarda como PNG.

    No requiere datos extra: reconstruye las lineas desde la grilla
    (grid_node), las vigas (coord_viga) y la config de arriostramientos
    y muros. No modifica en nada el modelo analizado.
    """
    if ruta is None:
        ruta = FIG_DIR / FIG_EDIFICIO

    xc = [cfg["grilla_ejes"]["X"][k] for k in cfg["orden_x"]]
    yc = [cfg["grilla_ejes"]["Y"][k] for k in cfg["orden_y"]]
    nx, ny = len(xc), len(yc)
    niveles = list(edificio.levels)

    # Escala de la deformada (para que se vea amplificada)
    tag_master = edificio.master_node[niveles[-1]]
    u_max = np.linalg.norm(ops.nodeDisp(tag_master)[:3])
    edificio._vis_escala = (0.1 / u_max) if u_max > 1e-12 else 0.0

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(projection="3d")

    # ---- Columnas (rojo): tramo base -> cada nivel, por cada interseccion
    for iy in range(ny):
        for ix in range(nx):
            n_prev = edificio.grid_node[(ix, iy, -1)]
            for k in range(len(niveles)):
                n_cur = edificio.grid_node[(ix, iy, k)]
                p0 = _coord(edificio, n_prev)
                p1 = _coord(edificio, n_cur)
                ax.plot(*zip(p0, p1), color="tab:red", lw=3, alpha=0.7)
                n_prev = n_cur

    # ---- Vigas (azul): de coord_viga
    for tag, (n1, n2, _tipo, _ix, _iy, _k) in edificio.coord_viga.items():
        p0 = _coord(edificio, n1)
        p1 = _coord(edificio, n2)
        ax.plot(*zip(p0, p1), color="tab:blue", lw=2.5, alpha=0.7)

    # ---- Voladizo (eje J): vigas cantilever + losa saliente + parapeto
    for lvl, beams in edificio.vol_vigas.items():
        sal = edificio.vol_saliente.get(lvl, 0.0)
        for (tag, n1, n2, L, iy) in beams:
            x1, y1, z1 = _coord(edificio, n2)
            p0 = _coord(edificio, n1)
            ax.plot(*zip(p0, (x1, y1, z1)), color="tab:blue", lw=3.2,
                    alpha=0.9)
            # losa saliente: del eje J a J+sal
            ax.plot([x1, x1 + sal], [y1, y1], [z1, z1],
                    color="tab:pink", lw=2.2, alpha=0.9)
            # parapeto 0.75 m en el borde
            ax.plot([x1 + sal, x1 + sal], [y1, y1],
                    [z1, z1 + 0.75], color="tab:purple", lw=3, alpha=0.8)

    # ---- Diagonales (verde): patron X en porticos exteriores
    ar = cfg["arriostramientos"]
    ks = niveles.index(ar["desde_nivel"])
    ke = niveles.index(ar["hasta_nivel"])
    for fy in ar["frames_Y"]:
        iy = cfg["orden_y"].index(fy)
        for ix in range(nx - 1):
            for k in range(ks, ke + 1):
                n_bot1 = edificio.grid_node[(ix, iy, k - 1)]
                n_bot2 = edificio.grid_node[(ix + 1, iy, k - 1)]
                n_top1 = edificio.grid_node[(ix, iy, k)]
                n_top2 = edificio.grid_node[(ix + 1, iy, k)]
                ax.plot(*zip(_coord(edificio, n_bot1), _coord(edificio, n_top2)),
                        color="tab:green", lw=2, alpha=0.8)
                ax.plot(*zip(_coord(edificio, n_bot2), _coord(edificio, n_top1)),
                        color="tab:green", lw=2, alpha=0.8)

    # ---- Muro de corte (naranja): linea de muro en cada tramo de nivel
    for m in cfg["muros"]["lista"]:
        orient = m["orientacion"]
        if orient == "Y":                     # muro paralelo a Y, sobre un eje X
            xf = cfg["grilla_ejes"]["X"][m["linea"]]
            y0 = cfg["grilla_ejes"]["Y"][m["desde_eje"]]
            y1 = cfg["grilla_ejes"]["Y"][m["hasta_eje"]]
            pts = [(xf, y0), (xf, y1)]
        else:                                 # muro paralelo a X, sobre un eje Y
            yf = cfg["grilla_ejes"]["Y"][m["linea"]]
            x0 = cfg["grilla_ejes"]["X"][m["desde_eje"]]
            x1 = cfg["grilla_ejes"]["X"][m["hasta_eje"]]
            pts = [(x0, yf), (x1, yf)]
        i0 = niveles.index(m["desde_nivel"])
        i1 = niveles.index(m["hasta_nivel"])
        zbot = edificio.z_base
        for k in range(i0, i1 + 1):
            ztop = edificio.z_level[niveles[k]]
            # contorno rectangular del muro en este tramo (en su plano)
            ax.plot([pts[0][0], pts[1][0]], [pts[0][1], pts[1][1]], [zbot, zbot],
                    color="tab:orange", lw=4, alpha=0.6)
            ax.plot([pts[0][0], pts[1][0]], [pts[0][1], pts[1][1]], [ztop, ztop],
                    color="tab:orange", lw=4, alpha=0.6)
            ax.plot([pts[0][0], pts[0][0]], [pts[0][1], pts[0][1]], [zbot, ztop],
                    color="tab:orange", lw=4, alpha=0.6)
            ax.plot([pts[1][0], pts[1][0]], [pts[1][1], pts[1][1]], [zbot, ztop],
                    color="tab:orange", lw=4, alpha=0.6)
            zbot = ztop

    # ---- Nodos (puntos)
    for nid, (x, y, z) in edificio.nodos.items():
        p = _coord(edificio, nid)
        ax.scatter(*p, s=8, color="black", alpha=0.5)

    # ---- Apariencia
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    xm = (min(xc) + max(xc)) / 2.0
    ym = (min(yc) + max(yc)) / 2.0
    ax.set_box_aspect((max(xc) - min(xc), max(yc) - min(yc),
                       max(edificio.z_level.values()) - edificio.z_base))
    ax.view_init(elev=25, azim=-55)

    from matplotlib.lines import Line2D
    legend = [
        Line2D([0], [0], color="tab:red", lw=3, label="Columnas"),
        Line2D([0], [0], color="tab:blue", lw=2.5, label="Vigas"),
        Line2D([0], [0], color="tab:green", lw=2, label="Diagonales (X)"),
        Line2D([0], [0], color="tab:orange", lw=4, label="Muro de corte"),
        Line2D([0], [0], color="tab:pink", lw=2.2, label="Losa voladizo (eje J)"),
        Line2D([0], [0], color="tab:purple", lw=3, label="Parapeto"),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=9)
    ax.set_title("Edificio institucional - deformada por gravedad "
                 "(ampl. x{:.0f})".format(1.0 / edificio._vis_escala)
                 if edificio._vis_escala else
                 "Edificio institucional")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(ruta, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return ruta


# =====================================================================
# 6b. EXPORTACION JSON PARA EL VIEWER (UNITY)
# =====================================================================

def exportar_json(edificio, cfg, resumen, area_piso_lvl, area_x, area_y,
                  w_map, w_vol, carga_aplicada):
    """Exporta el modelo y sus verificaciones a JSON legibles por el
    viewer/Unity. No modifica el analisis ya realizado.

    Escribe en results/export/:
      nodos.json, elementos.json, diafragmas.json, apoyos.json,
      secciones.json, tributarias.json, verificaciones.json
    """
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. NODOS --------------------------------------------------
    nodos = []
    for nid, (x, y, z) in edificio.nodos.items():
        nodos.append({"tag": nid, "x": x, "y": y, "z": z})
    _write(EXPORT_DIR / "nodos.json", {"nodos": nodos})

    # ---- 1b. GRILLA (ejes/niveles por nombre, para el viewer) -----
    def _ejes_limpos(d):
        return {k: v for k, v in d.items() if k != "nota"}

    _write(EXPORT_DIR / "grilla.json", {
        "ejes_X": _ejes_limpos(cfg["grilla_ejes"].get("X", {})),
        "ejes_Y": _ejes_limpos(cfg["grilla_ejes"].get("Y", {})),
        "niveles": _ejes_limpos(cfg["niveles"])})

    # ---- 2. APOYOS ------------------------------------------------
    apoyos = [{"tag": nid, "ux": a[0], "uy": a[1], "uz": a[2],
               "rx": a[3], "ry": a[4], "rz": a[5]}
              for nid, a in edificio.apoyos.items()]
    _write(EXPORT_DIR / "apoyos.json", {"apoyos": apoyos, "n": len(apoyos)})

    # ---- 3. ELEMENTOS (columnas + vigas + diagonales + muros) ------
    secc_vigas = cfg["secciones_tipo"]["vigas"]
    asig_vigas = cfg["asignaciones"]["vigas_por_nivel"]
    elementos = [dict(c) for c in edificio.columnas]
    elementos += [dict(d) for d in edificio.diagonales]
    for tag, (n1, n2, tipo, ix, iy, k) in edificio.coord_viga.items():
        lvl = edificio.levels[k]
        seccion = asig_vigas[lvl]
        material = secc_vigas[seccion]["material"]
        c0 = np.array(edificio.nodos[n1])
        c1 = np.array(edificio.nodos[n2])
        elementos.append({
            "elementTag": tag, "tipo": "viga",
            "n1": n1, "n2": n2, "nivel": lvl,
            "seccion": seccion, "material": material,
            "transf": edificio._transf["viga_x"] if tipo == "X"
                      else edificio._transf["viga_y"],
            "orientacion": tipo,  # eje de la viga en planta (X o Y)
            "ix": ix, "iy": iy, "k": k,
            "L": float(np.linalg.norm(c1 - c0))})
    # Muros equivalentes (geometria de su plano, no como EF)
    for m in cfg["muros"]["lista"]:
        orient = m["orientacion"]
        seccion = m["seccion"]
        materia_m = cfg["secciones_tipo"]["muros"][seccion]["material"]
        elementos.append({
            "tipo": "muro", "nombre": m["nombre"], "orientacion": orient,
            "linea_eje": m["linea"], "desde_eje": m["desde_eje"],
            "hasta_eje": m["hasta_eje"], "desde_nivel": m["desde_nivel"],
            "hasta_nivel": m["hasta_nivel"],
            "seccion": seccion, "material": materia_m,
            "espesor": cfg["secciones_tipo"]["muros"][seccion]["t"]})
    # Vigas cantilever del voladizo (eje J)
    secc_vol = cfg["secciones_tipo"]["vigas"]
    for lvl, beams in edificio.vol_vigas.items():
        for (tag, n1, n2, L, iy) in beams:
            vv = next(v for v in cfg["voladizos"]["lista"]
                      if lvl in v["niveles"])
            seccion = vv["seccion_viga"]
            material = secc_vol[seccion]["material"]
            c0 = np.array(edificio.nodos[n1])
            c1 = np.array(edificio.nodos[n2])
            elementos.append({
                "elementTag": tag, "tipo": "viga", "es_voladizo": True,
                "n1": n1, "n2": n2, "nivel": lvl,
                "seccion": seccion, "material": material,
                "transf": edificio._transf["viga_x"],
                "orientacion": "X",
                "desde_eje": vv["desde_eje"], "hasta_eje": vv["hasta_eje"],
                "eje_Y": cfg["orden_y"][iy], "k": edificio.vol_coord[tag][3],
                "L": float(np.linalg.norm(c1 - c0))})
    _write(EXPORT_DIR / "elementos.json",
           {"nodos": list(edificio.nodos.keys()), "elementos": elementos})

    # ---- 4. DIAFRAGMAS --------------------------------------------
    diafragmas = [{"nivel": lvl, "master": ed["master"],
                   "slaves": ed["slaves"], "direccion": 3}
                  for lvl, ed in edificio.diafragmas.items()]
    _write(EXPORT_DIR / "diafragmas.json", {"diafragmas": diafragmas})

    # ---- 5. SECCIONES Y MATERIALES --------------------------------
    _write(EXPORT_DIR / "secciones.json",
           {"secciones_tipo": cfg["secciones_tipo"],
            "materiales": cfg["materiales"]})

    # ---- 6. AREAS TRIBUTARIAS Y CARGAS DE LOSA --------------------
    # tag_viga -> lista de (tag, w) acumulada. Las vigas X/Y se
    # identifican por (tipo, ix, iy, k) en w_map y por coord_viga.
    # Reconstruimos el id de viga usado en w_map.
    xc = [cfg["grilla_ejes"]["X"][k] for k in cfg["orden_x"]]
    yc = [cfg["grilla_ejes"]["Y"][k] for k in cfg["orden_y"]]
    dx = [xc[i + 1] - xc[i] for i in range(len(xc) - 1)]
    dy = [yc[j + 1] - yc[j] for j in range(len(yc) - 1)]
    # w_map: { (tipo, i, j, lvl): w }
    trib = []
    for tag, (n1, n2, tipo, ix, iy, k) in edificio.coord_viga.items():
        lvl = edificio.levels[k]
        w = w_map.get((tipo, ix, iy, lvl), 0.0)
        L = dx[ix] if tipo == "X" else dy[iy]
        p1 = np.array(edificio.nodos[n1])
        p2 = np.array(edificio.nodos[n2])
        Lr = float(np.linalg.norm(p2 - p1))
        trib.append({
            "elementTag": tag, "tipo": "viga", "nivel": lvl,
            "orientacion": tipo, "ix": ix, "iy": iy, "k": k,
            "w_kN_m": w, "longitud_m": Lr,
            "asa_x_m": dx[ix] if tipo == "X" else 0.0,
            "asa_y_m": dy[iy] if tipo == "Y" else 0.0,
            "carga_losa_kN": w * Lr})
    _write(EXPORT_DIR / "tributarias.json", {
        "qG": {k: v for k, v in cfg["cargas"]["qG"].items() if k != "nota"},
        "area_piso_m2": {lvl: area_piso_lvl[lvl]
                         for lvl in cfg["orden_niveles"]},
        "carga_total_losa_kN": sum(r["carga_total"] for r in resumen.values()),
        "carga_puntual_kN": carga_aplicada
                            - sum(r["carga_total"] for r in resumen.values()),
        "voladizo": {
            "descripcion": cfg["voladizos"].get("nota", ""),
            "por_nivel": [
                {"nivel": lvl, "area_losa_voladizo_m2": r["area_voladizo"],
                 "carga_losa_voladizo_kN": r["carga_voladizo"]}
                for lvl, r in resumen.items() if r["area_voladizo"] > 0],
            "vigas_cantilever": [
                {"elementTag": tag, "nivel": lvl,
                 "eje_Y": cfg["orden_y"][iy],
                 "desde_m": cfg["grilla_ejes"]["X"][
                     vv["desde_eje"]],
                 "hasta_m": vv["x_j_m"],
                 "saliente_m": edificio.vol_saliente.get(lvl, 0.0),
                 "longitud_m": L,
                 "w_kN_m": w_vol[lvl][iy] if w_vol and lvl in w_vol
                            and iy in w_vol[lvl] else 0.0}
                for lvl, beams in edificio.vol_vigas.items()
                for (tag, n1, n2, L, iy) in beams
                for vv in cfg["voladizos"]["lista"] if lvl in vv["niveles"]]},
        "vigas": trib})

    # ---- 7. VERIFICACIONES ----------------------------------------
    reacs = {}
    for nid, (x, y, z) in edificio.nodos.items():
        if abs(z - edificio.z_base) < 1e-9:
            reacs[nid] = ops.nodeReaction(nid)
    sum_fx = sum(r[0] for r in reacs.values())
    sum_fy = sum(r[1] for r in reacs.values())
    sum_fz = sum(r[2] for r in reacs.values())
    _write(EXPORT_DIR / "verificaciones.json", {
        "carga_por_piso": [
            {"nivel": lvl, "qG_kN_m2": resumen[lvl]["qG"],
             "carga_losa_kN": resumen[lvl]["carga_total"],
             "area_piso_m2": resumen[lvl]["area_piso"],
             "area_trib_X": resumen[lvl]["area_trib_X"],
             "area_trib_Y": resumen[lvl]["area_trib_Y"],
             "area_voladizo_m2": resumen[lvl]["area_voladizo"],
             "igual_a_area_piso":
                 abs(resumen[lvl]["area_trib_X"] + resumen[lvl]["area_trib_Y"]
                     - resumen[lvl]["area_piso"]) < 1e-6}
            for lvl in cfg["orden_niveles"]],
        "conservacion_carga": {
            "carga_total_aplicada_kN": carga_aplicada,
            "suma_reacciones_FZ_kN": sum_fz,
            "diferencia_kN": abs(carga_aplicada) - abs(sum_fz),
            "suma_reacciones_FX_kN": sum_fx,
            "suma_reacciones_FY_kN": sum_fy},
        "equilibrio_global": {
            "suma_FX_kN": sum_fx, "suma_FY_kN": sum_fy,
            "suma_FZ_kN": sum_fz},
        "reacciones_basales": [
            {"nodo": nid, "x": edificio.nodos[nid][0],
             "y": edificio.nodos[nid][1], "fx": r[0], "fy": r[1],
             "fz": r[2], "mx": r[3], "my": r[4], "mz": r[5]}
            for nid, r in reacs.items()]})
    return EXPORT_DIR


def _write(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


# =====================================================================
# 7. MAIN
# =====================================================================

def main():
    ruta = sys.argv[1] if len(sys.argv) > 1 else "data/edificio_config.json"
    print("=" * 78)
    print("EDIFICIO INSTITUCIONAL - ANALISIS DE GRAVEDAD")
    print("Proyecto 1 Grupo 4 | OpenSeesPy")
    print("=" * 78)

    cfg = cargar_config(ruta)
    print(f"\nCargando config: {ruta}")

    # Areas tributarias (antes de construir, para saber cargas)
    w_map, resumen, area_piso_lvl, area_x, area_y, w_vol = \
        areas_tributarias(cfg)

    # Construir modelo
    print("\n[1] Construyendo modelo...")
    ed = Edificio(cfg)
    ed.construir()
    print(f"    {len(ed.nodos)} nodos | {ed.n_cols} columnas | "
          f"{ed.n_vigas} vigas | {ed.n_vol} vigas-voladizo | "
          f"{ed.n_diag} diagonales | {len(ed.muro_info)} muros")

    # Cargas
    print("\n[2] Aplicando cargas (areas tributarias + puntuales)...")
    carga_aplicada = aplicar_cargas(cfg, ed, w_map, w_vol)
    print(f"    Carga total aplicada: {carga_aplicada:.4f} kN")

    # Analisis
    print("\n[3] Analisis estatico lineal de gravedad...")
    analizar()
    print("    OK")

    # QA
    print("\n[4] Verificacion QA / exportacion CSV:")
    qa_verificacion(cfg, ed, resumen, area_piso_lvl, area_x, area_y,
                    carga_aplicada)

    # Visualizacion 3D con deformada
    print("\n[5] Generando vista 3D con deformada por gravedad...")
    ruta_fig = visualizar(ed, cfg)
    print(f"    Figura guardada en: {ruta_fig}")

    # Exportacion JSON para el viewer/Unity
    print("\n[6] Exportando modelo/verificaciones a JSON (Unity viewer)...")
    ruta_exp = exportar_json(ed, cfg, resumen, area_piso_lvl, area_x, area_y,
                             w_map, w_vol, carga_aplicada)
    print(f"    Exportados en: {ruta_exp}")
    for f in sorted(ruta_exp.iterdir()):
        print(f"      {f.name}")


if __name__ == "__main__":
    main()
