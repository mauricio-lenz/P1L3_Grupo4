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


def capacidad_axial_pura(cfg):
    """Capacidad axial pura de la columna (compresion, kN, signo negativo):

        P_o = 0.85 fc' (Ag - As) + fy As     (ACI 318, columna estribada)

    Es el punto superior de la interaccion P-M (M = 0).
    """
    p = propiedades_seccion(cfg)
    fc = p["fc_MPa"] * 1e3              # MPa -> kN/m2
    fy = p["fy_MPa"] * 1e3
    return -(0.85 * fc * (p["Ag_m2"] - p["As_m2"]) + fy * p["As_m2"])


def graficar_punto_pm(cfg, Ps, Mu, out_file, demanda=None):
    """Primeros puntos de la curva de interaccion P-M (Mu max de M-phi).

    Incluye el punto de compresion axial pura (M = 0) para que la curva
    cierre arriba y muestre el cambio de pendiente alrededor del punto
    balanceado.
    """
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    idxy = np.argsort(Ps)
    ax.plot(Mu[idxy], Ps[idxy], "o-", lw=1.7, ms=5,
            label="Puntos P-M (M$_u$)")
    if demanda is not None:
        ax.plot(demanda["M"], demanda["P"], "r*", ms=16,
                label="Demanda max (G,Q,EX,EY)")
    ax.set_xlabel("Momento ultimo M$_u$ [kN-m]")
    ax.set_ylabel("Carga axial P [kN] (negativo = compresion)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close(fig)


def run_capacidad(cfg, outdir, P_grid=None, demanda=None):
    """Ejecuta toda la parte D y guarda figuras + curvas.

    Puntos de la curva P-M (5): compresion axial pura (M=0) + 4 niveles de
    carga axial de las curvas M-phi.

    Devuelve dict con curvas y Puntos P-M.
    """
    if P_grid is None:
        P_grid = np.array([-5000.0, -3000.0, -1000.0, 0.0])
    curvas = {}
    for P in P_grid:
        curvas[P] = moment_curvature(cfg, P)
    Mu = np.array([np.max(Ms) for _, Ms in curvas.values()])
    Ps = np.array([float(p) for p in curvas.keys()])

    # Punto de compresion axial pura (M = 0) al tope de la curva
    P_axial = capacidad_axial_pura(cfg)
    Ps_curve = np.concatenate(([P_axial], Ps))
    Mu_curve = np.concatenate(([0.0], Mu))

    graficar_seccion(cfg, outdir / "fig_seccion_fibras.png")
    graficar_phi_mu(Ps, curvas, outdir / "fig_M_phi_columna.png")
    graficar_punto_pm(cfg, Ps_curve, Mu_curve, outdir / "fig_PM_columna.png",
                      demanda=demanda)
    return {"Ps": Ps_curve, "Mu": Mu_curve, "curvas": curvas,
            "P_axial": P_axial}