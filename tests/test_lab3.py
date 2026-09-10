"""Pruebas P1L3 - Grupo 4.

Tests unitarios rapidos (config, seccion, superposicion algebraica) y uno de
integracion (construccion del edificio + caso G con conservacion exacta).
"""

import json
from pathlib import Path

import numpy as np
import pytest

from src import capacidad as cap
from src import casos as cs
from src import geometria as geo
from src import superposicion as sup
from src.edificio import cargar_config

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


@pytest.fixture(scope="session")
def cfg():
    c = cargar_config(DATA / "edificio_config.json")
    c["cap_hormigon"] = json.loads(
        (DATA / "casos_lab3.json").read_text(encoding="utf-8"))["cap_hormigon"]
    return c


@pytest.fixture(scope="session")
def casos():
    return json.loads((DATA / "casos_lab3.json").read_text(encoding="utf-8"))


def test_config_casos_lab(casos):
    assert set(casos["carga_viva_q_Q_kN_m2"]) >= {
        "Subterraneo", "Piso1", "Piso2", "Piso3", "Piso4"}
    assert casos["sismo_pseudoestatico"]["fraccion_g"] == 0.20
    lam = casos["superposicion"]["lambdas"]
    assert lam["G"] == 1.0 and lam["Q"] == 1.0
    assert lam["EX"] == 0.90 and lam["EY"] == 0.75
    assert casos["cap_hormigon"]["refuerzo_longitudinal"]["n_barras"] == 8


def test_propiedades_seccion(cfg):
    p = cap.propiedades_seccion(cfg)
    assert p["Ag_m2"] == pytest.approx(0.49, rel=1e-9)
    assert p["As_m2"] == pytest.approx(8 * np.pi * 0.025 ** 2 / 4, rel=1e-9)
    assert p["cuantia"] == pytest.approx(p["As_m2"] / p["Ag_m2"], rel=1e-9)
    assert p["fc_MPa"] == 25 and p["fy_MPa"] == 420


def test_discretizacion_fibras(cfg):
    fibras, A_br = cap.discretizacion_fibras(cfg)
    n_barras = sum(1 for f in fibras if f[3] == 3)
    assert n_barras == 8
    assert A_br == pytest.approx(np.pi * 0.025 ** 2 / 4, rel=1e-9)
    assert len(fibras) >= 200


def test_combinar_superposicion():
    r = sup.Resultado.__new__(sup.Resultado)
    r.desp = {"Piso4": {"G": np.zeros(6), "Q": np.ones(6) * 2.0,
                        "EX": np.ones(6), "EY": np.zeros(6)}}
    r.reac = {}
    r.fuerzas = {}
    r._reac_tot = np.zeros(6)
    lam = {"G": 1.0, "Q": 1.0, "EX": 0.9, "EY": 0.75}
    # construye 4 resultados ficticios
    res = {}
    base = {"desp": {"Piso4": np.zeros(6)}, "fuerzas": {},
            "reac": {}, "tot": np.zeros(6)}
    for c in lam:
        rr = sup.Resultado.__new__(sup.Resultado)
        rr.desp = {"Piso4": base["desp"]["Piso4"].copy()}
        rr.reac, rr.fuerzas, rr._reac_tot = {}, {}, np.zeros(6)
        res[c] = rr
    res["Q"].desp["Piso4"] = np.ones(6) * 2.0
    comb = sup.combinar({"orden_niveles": ["Piso4"]}, res, lam)
    assert np.allclose(comb.desp["Piso4"], np.ones(6) * 2.0)
    assert comb.reacciones_totales.shape == (6,)


def test_caso_G_conservacion(cfg, casos):
    import openseespy.opensees as ops
    ed = cs.construir_edificio(cfg)
    tags = cs.seleccionar_elementos(ed)
    assert tags["columna_critica"] is not None
    assert tags["viga_verifica"] is not None
    qG, _ = geo.carga_losa(cfg, casos)
    w_map, resumen, apisos, ax, ay, w_vol, total = \
        cs.correr_caso_gravedad(cfg, ed, qG, con_puntuales=True)
    res = cs.registrar(cfg, ed, tags)
    fz = sum(r[2] for r in res.reac.values())
    assert total == pytest.approx(sum(r["carga_total"]
                                      for r in resumen.values()), rel=1e-9)
    assert fz == pytest.approx(total, rel=1e-6)
    assert len(res.desp) == len(cfg["orden_niveles"])
    ops.wipe()