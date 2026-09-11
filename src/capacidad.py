"""
CAPACIDAD DE LA COLUMNA HA (P1L3)
=================================
P1L3 - Grupo 4

Parte D: se construye una Fiber Section de la columna critica (0.70x0.70,
fc 25 MPa, 8xphi25, recubrimiento 5 cm) mediante el comando section Fiber
de OpenSees, se discretiza el nucleo y el recubrimiento con fibra de
concreto y el refuerzo con layer de acero. Sobre ella se corren analisis
de curvatura (M-phi) para distintos niveles de carga axial P y con esos
maximos se dibujan los primeros puntos de la curva de interaccion P-M.
"""

import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt


def cap(cfg):
    """Seccion HA de la parte D editada desde casos_lab3.json."""
    return cfg["cap_hormigon"]


def discretizacion_fibras(cfg):
    """Discretizacion de fibra equivalente a la del section Fiber.

    Devuelve una lista de (y, z, A, tag_material) para dibujar la figura:
    nucleo (mat 1), recubrimiento (mat 2) y barras (mat 3).
    """
    hch = cap(cfg)
    b, h = hch["b_m"], hch["h_m"]
    cover = hch["recubrimiento_libre_m"]
    nc = hch["discretizacion"]["n_divisiones_nucleo"]
    diam = hch["refuerzo_longitudinal"]["diametro_m"]

    y1, z1 = h / 2.0, b / 2.0
    cyl = np.linspace(-y1 + cover, y1 - cover, nc)
    czl = np.linspace(-z1 + cover, z1 - cover, nc)
    dy = 2 * (y1 - cover) / nc
    dz = 2 * (z1 - cover) / nc

    out = []
    for yy in cyl:
        for zz in czl:
            out.append((yy, zz, dy * dz, 1))
    t = cover
    ny = np.linspace(-y1 + t / 2, y1 - t / 2, nc)
    nz = np.linspace(-z1 + t / 2, z1 - t / 2, nc)
    for yy in ny:
        out.append((yy, z1 - t / 2, (2 * (y1 - t)) / nc * t, 2))
        out.append((yy, -z1 + t / 2, (2 * (y1 - t)) / nc * t, 2))
    for zz in nz:
        out.append((y1 - t / 2, zz, (2 * (z1 - t)) / nc * t, 2))
        out.append((-y1 + t / 2, zz, (2 * (z1 - t)) / nc * t, 2))
    A_br = np.pi * diam ** 2 / 4.0
    for zz in (z1 - cover, -z1 + cover):
        for yy in np.linspace(-y1 + cover, y1 - cover, 3):
            out.append((yy, zz, A_br, 3))
    for yy in (y1 - cover, -y1 + cover):
        out.append((yy, 0.0, A_br, 3))
    return out, A_br


def construir_seccion_fibra(cfg):
    """Construye la Fiber Section de la columna en OpenSeesPy (ndm=2).

    Crea los materiales Concrete01 (nucleo y recubrimiento) y Steel01, y el
    section Fiber con patch de nucleo, recubrimiento y layer longitudinal.

    Convencion OpenSees: compresion con valores negativos.
    """
    hch = cap(cfg)
    b, h = hch["b_m"], hch["h_m"]
    cover = hch["recubrimiento_libre_m"]
    nc = hch["discretizacion"]["n_divisiones_nucleo"]
    nb = hch["refuerzo_longitudinal"]["n_barras"]
    diam = hch["refuerzo_longitudinal"]["diametro_m"]
    nl = hch["materiales_no_lineales"]
    c01 = nl["Concrete01"]
    s01 = nl["Steel01"]

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.uniaxialMaterial("Concrete01", 1, c01["fpc_kN_m2"],
                         c01["epsc0"], c01["fpcu_kN_m2"], c01["epsu"])
    ops.uniaxialMaterial("Concrete01", 2, c01["fpc_kN_m2"],
                         c01["epsc0"], 0.0, -0.006)
    ops.uniaxialMaterial("Steel01", 3, s01["fy_kN_m2"], s01["Es_kN_m2"],
                         s01["b_ratio"])

    y1, z1 = h / 2.0, b / 2.0
    ops.section("Fiber", 1)
    ops.patch("rect", 1, nc, nc, cover - y1, cover - z1,
              y1 - cover, z1 - cover)
    ops.patch("rect", 2, nc, 1, -y1, z1 - cover, y1, z1)
    ops.patch("rect", 2, nc, 1, -y1, -z1, y1, cover - z1)
    ops.patch("rect", 2, 1, nc, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch("rect", 2, 1, nc, y1 - cover, cover - z1, y1, z1 - cover)

    As = np.pi * diam ** 2 / 4.0
    # 3+3 en caras superiores/inferiores, 1+1 en costados => 8 barras
    ops.layer("straight", 3, 3, As, cover - y1, cover - z1,
              y1 - cover, cover - z1)
    ops.layer("straight", 3, 3, As, cover - y1, z1 - cover,
              y1 - cover, z1 - cover)
    ops.layer("straight", 3, 1, As, cover - y1, 0.0, cover - y1, 0.0)
    ops.layer("straight", 3, 1, As, y1 - cover, 0.0, y1 - cover, 0.0)

    return {"sec_tag": 1, "As_bar": As, "n_fibras_steel": nb}


def moment_curvature(cfg, P, maxK=0.06, numIncr=1200, verbose=False):
    """Analisis M-phi de la fiber section bajo carga axial constante P.

    Convencion de OpenSees: P < 0 = compresion, P > 0 = traccion.
    Devuelve (phis, Ms) con unidades [1/m] y [kN-m].
    """
    construir_seccion_fibra(cfg)

    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    ops.element("zeroLengthSection", 1, 1, 2, 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, P, 0.0, 0.0)

    ops.integrator("LoadControl", 0.0)
    ops.system("SparseGeneral", "-piv")
    ops.test("NormUnbalance", 1e-9, 30, 0)
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.algorithm("Newton")
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError("no converge la etapa axial P=%.0f" % P)

    ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, 0.0, 1.0)

    dK = maxK / numIncr
    ops.integrator("DisplacementControl", 2, 3, dK, 1, dK, dK)
    phis, Ms = [], []
    for _ in range(numIncr):
        o = ops.analyze(1)
        if o != 0:
            if verbose:
                print("M-phi (P=%.0f): corte en paso %d (err %d)"
                      % (P, _, o))
            break
        phis.append(ops.nodeDisp(2)[2])
        ops.reactions()
        Ms.append(-ops.nodeReaction(1)[2])
    if not phis:
        raise RuntimeError("no converge el analisis M-phi con P=%.0f" % P)
    return np.array(phis), np.array(Ms)


def propiedades_seccion(cfg):
    """Propiedades para el informe (dimensiones, cuantia, EI, etc.)."""
    hch = cap(cfg)
    b, h = hch["b_m"], hch["h_m"]
    nb = hch["refuerzo_longitudinal"]["n_barras"]
    diam = hch["refuerzo_longitudinal"]["diametro_m"]
    fc_MPa = hch["fc_MPa"]
    fy_MPa = hch["fy_MPa"]
    Es = hch["Es_GPa"] * 1e6
    As = nb * np.pi * diam ** 2 / 4.0
    Ag = b * h
    Ec = 4700 * np.sqrt(fc_MPa) * 1e3          # MPa -> kN/m2 (ACI 318)
    return {
        "b_m": b, "h_m": h, "Ag_m2": Ag, "n_barras": nb, "diam_m": diam,
        "As_m2": As, "cuantia": As / Ag,
        "recubrimiento_m": hch["recubrimiento_libre_m"],
        "fc_MPa": fc_MPa, "fy_MPa": fy_MPa, "Es_kN_m2": Es,
        "d_m": h - hch["recubrimiento_libre_m"], "Ec_kN_m2": Ec,
        "EI_gross_kN_m2": Ec * b * h ** 3 / 12.0,
        "eps_core": hch["materiales_no_lineales"]["Concrete01"]["epsu"],
    }


def graficar_seccion(cfg, out_file):
    """Figura con la discretizacion de fibras y el refuerzo."""
    hch = cap(cfg)
    fibras, A_br = discretizacion_fibras(cfg)
    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    col = {1: "#d8d8d8", 2: "#f2a4a4", 3: "#222222"}
    seen = set()
    for f in fibras:
        y, z, A, mat = f
        A = float(A)
        s = max(4.0, min(60.0, 1100.0 * np.sqrt(A)))
        lbl = {1: "Nucleo confinado", 2: "Recubrimiento",
               3: "Refuerzo 8$\\phi$25"}[int(mat)]
        key = "lin" if mat in (1, 2) else "bar"
        ax.scatter(y, z, s=s, c=col[int(mat)], linewidths=0,
                   label=lbl if key not in seen else None)
        seen.add(key)
    ax.add_patch(plt.Rectangle((-0.35, -0.35), 0.70, 0.70,
                               fill=False, ec="k", lw=1.4))
    ax.set_xlim(-0.45, 0.45)
    ax.set_ylim(-0.45, 0.45)
    ax.set_xlabel("y [m] (eje fuerte)")
    ax.set_ylabel("z [m]")
    ax.set_aspect("equal")
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close(fig)
    return {"n_fibras_concreto": sum(1 for f in fibras if f[3] < 3),
            "n_barras": sum(1 for f in fibras if f[3] == 3),
            "As_bar_m2": A_br}


def graficar_phi_mu(Ps, curvas, out_file, demandas=None):
    """Figura con las curvas M-phi para los niveles de carga axial."""
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for P, (phis, Ms) in curvas.items():
        ax.plot(phis * 1e3, Ms, lw=1.6, label="P = %+6.0f kN" % P)
    if demandas:
        ax.axvline(0, c="k", lw=0.6, alpha=0.4)
    ax.set_xlabel("Curvatura $\\phi$ $\\times 10^{-3}$ [1/m]")
    ax.set_ylabel("Momento M [kN-m]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close(fig)


EC_U = 0.003                 # deformacion ultima del concreto no confinado
EPS_T_COM_CONTROL = 0.002    # limite control de compresion (ACI 318, estribada)


def _esfuerzos_seccion(cfg, c):
    """Esfuerzos P-M por compatibilidad de deformaciones (bloque de Whitney).

    Convencion: P positivo = compresion; M sobre el eje centroidal. Recibe
    la profundidad `c` del eje neutro desde el borde comprimido y devuelve
    (P, M) en kN y kN-m para una deformacion de borde eps_c = 0.003.

    Acero distribuido en 3 capas (eje fuerte de los 8 phi25: 3 + 2 + 3).
    """
    p = propiedades_seccion(cfg)
    b, h = p["b_m"], p["h_m"]
    fc = p["fc_MPa"] * 1e3
    fy = p["fy_MPa"] * 1e3
    Es = p["Es_kN_m2"]
    As_bar = p["As_m2"] / p["n_barras"]
    cover = p["recubrimiento_m"]
    capas = [(cover, 3), (h / 2.0, 2), (h - cover, 3)]  # (dist. desde borde, n)

    a = min(0.85 * c, h)                       # bloque rectangular de Whitney
    C = 0.85 * fc * b * a
    P = C
    M = C * (h / 2.0 - a / 2.0)
    for depth, n in capas:
        eps = EC_U * (1.0 - depth / c)
        fs = min(max(Es * eps, -fy), fy)
        F = n * As_bar * fs
        P += F
        M += F * (h / 2.0 - depth)
    return P, M


def _c_por_eps_t(cfg, eps_t):
    """Profundidad del eje neutro para una deformacion dada del acero
    extremo en traccion (traccion + compresion en borde = 0.003)."""
    p = propiedades_seccion(cfg)
    d = p["d_m"]
    return d * EC_U / (EC_U + eps_t)


def puntos_interaccion(cfg):
    """Los 5 puntos caracteristicos del diagrama de interaccion P-M.

        P_o  = compresion axial pura: 0.85 fc'(Ag-As) + fy As    (M = 0)
        cc   = control de compresion: eps_t = 0.002 en acero extremo
        bal  = condicion balanceada:  eps_t = eps_y = fy/Es
        flex = flexion pura:          P = 0
        ten  = tension axial pura:    todo el acero fluye a traccion (M = 0)

    Convencion: P positivo = compresion (los puntos axiales caen sobre el
    eje Y, con M = 0). Devuelve lista de dicts {nombre, P, M}.
    """
    p = propiedades_seccion(cfg)
    fc = p["fc_MPa"] * 1e3
    fy = p["fy_MPa"] * 1e3
    Es = p["Es_kN_m2"]
    ey = fy / Es

    # 1) Compresion axial pura (M = 0, sobre el eje Y)
    P_axial = 0.85 * fc * (p["Ag_m2"] - p["As_m2"]) + fy * p["As_m2"]

    # 2) Control de compresion (borde de la zona controlada por compresion)
    P_cc, M_cc = _esfuerzos_seccion(cfg, _c_por_eps_t(cfg, EPS_T_COM_CONTROL))

    # 3) Condicion balanceada (acero extremo en fluencia simultanea)
    P_bal, M_bal = _esfuerzos_seccion(cfg, _c_por_eps_t(cfg, ey))

    # 4) Flexion pura: P(c) decrece con c (P<0 en traccion profunda).
    #    mantener lo con P<0 y hi con P>0.
    lo, hi = 1e-4, _c_por_eps_t(cfg, ey)
    for _ in range(300):
        c = 0.5 * (lo + hi)
        if _esfuerzos_seccion(cfg, c)[0] > 0:
            hi = c
        else:
            lo = c
    P_flex, M_flex = _esfuerzos_seccion(cfg, 0.5 * (lo + hi))

    # 5) Tension axial pura (M = 0, sobre el eje Y)
    P_tension = -fy * p["As_m2"]

    return [
        {"nombre": "Compresion axial pura", "P": P_axial, "M": 0.0},
        {"nombre": "Control de compresion", "P": P_cc, "M": M_cc},
        {"nombre": "Condicion balanceada", "P": P_bal, "M": M_bal},
        {"nombre": "Flexion pura", "P": P_flex, "M": M_flex},
        {"nombre": "Tension axial pura", "P": P_tension, "M": 0.0},
    ]


def graficar_interaccion(cfg, puntos, out_file, demanda=None, fibras=None):
    """Diagrama de interaccion P-M (5 puntos caracteristicos).

    P (compresion positiva) en el eje Y, M en el eje X. El primer y el
    ultimo punto (compresion y tension axial pura) caen sobre el eje Y
    (M = 0) con el quiebre de pendiente alrededor del punto balanceado.

    `fibras` = {"Ps_ops": [...], "Mu": [...]} para superponer la validacion
    con el maximo M de las curvas M-phi (OpenSees, P ops negativo).
    """
    Ps = [pt["P"] for pt in puntos]
    Ms = [pt["M"] for pt in puntos]
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    ax.axvline(0, c="k", lw=0.8, alpha=0.4, ls="--", zorder=0)
    ax.plot(Ms, Ps, "o-", lw=1.8, ms=6, color="#0b5394", zorder=3,
            label="Curva de interaccion P-M")
    if fibras is not None:
        Ps_f = [-float(P) for P in fibras["Ps_ops"]]
        ax.plot(fibras["Mu"], Ps_f, "o", ms=6, mfc="none", mec="#e06666",
                zorder=2, label="Max M-phi (fiber section)")
    if demanda is not None:
        ax.plot(demanda["M"], -demanda["P_ops"], "r*", ms=20, zorder=4,
                label="Demanda max (G,Q,EX,EY)")
    for pt in puntos:
        ax.annotate(pt["nombre"], (pt["M"], pt["P"]),
                    textcoords="offset points", xytext=(9, 7),
                    fontsize=7.5, zorder=5)
    ax.set_xlabel("Momento flector M [kN-m]")
    ax.set_ylabel("Carga axial P [kN]   (arriba = compresion)")
    ax.set_ylim(1.25 * min(Ps), 1.03 * max(Ps))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close(fig)


def run_capacidad(cfg, outdir, P_grid=None, demanda=None):
    """Ejecuta toda la parte D y guarda figuras + curvas.

    1. Diagrama de interaccion P-M con los 5 puntos caracteristicos
       (compresion axial, control de compresion, balanceada, flexion pura,
       tension axial) por compatibilidad de deformaciones + Whitney.
    2. Validacion del overlay: maximo M de las curvas M-phi de la seccion
       de fibras (OpenSees) para los niveles de carga axial de P_grid.

    `demanda` (opcional): {"P_ops": float, "M": float} con P_ops negativo
    = compresion (convencion OpenSees). Devuelve dict con los puntos.
    """
    if P_grid is None:
        P_grid = np.array([-5000.0, -3000.0, -1000.0, 0.0])
    curvas = {}
    for P in P_grid:
        curvas[P] = moment_curvature(cfg, P)
    Mu = np.array([np.max(Ms) for _, Ms in curvas.values()])
    Ps_ops = np.array([float(p) for p in curvas.keys()])

    puntos = puntos_interaccion(cfg)
    P_axial = puntos[0]["P"]

    graficar_seccion(cfg, outdir / "fig_seccion_fibras.png")
    graficar_phi_mu(Ps_ops, curvas, outdir / "fig_M_phi_columna.png")
    graficar_interaccion(cfg, puntos, outdir / "fig_PM_columna.png",
                         demanda=demanda, fibras={"Ps_ops": Ps_ops, "Mu": Mu})
    return {"puntos": puntos, "P_axial": P_axial,
            "Ps_ops": Ps_ops, "Mu": Mu, "curvas": curvas}