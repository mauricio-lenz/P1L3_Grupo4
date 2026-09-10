# P1L3 Grupo 4 — Lab: Carga viva, sismo, superposición y capacidad HA

Laboratorio P1L3 (curso Estructuras, Edificio U. de los Andes 2017_67).
Modelo estructural 3D OpenSeesPy del edificio (geometría de Semana 2) con:

- **Parte A** — Caso de carga viva Q (NCh1537, 3.0/1.5 kN/m²) sobre la misma
  geometría tributaria que el caso G; verificación de conservación
  `sum(Q transferida) = q_Q·A` y contra reacciones basales.
- **Parte B** — Sismo pseudoestático EX y EY: `F_i = 0.20·W_i` aplicada en el
  centro de masa de cada piso (con su torsión de piso). Verificación de corte
  basal, sentido de la deformada y rotación Rz por piso.
- **Parte C** — Superposición `R = λG·G + λQ·Q + λEX·EX + λEY·EY` sumando las
  respuestas unitarias y comparándola con una corrida OpenSees explícita con
  todas las cargas combinadas (desplazamiento, reacción y fuerza interna).
- **Parte D** — Capacidad de la columna P-70x70 con `section Fiber` de
  OpenSees (Concrete01 + Steel01): curvas momento-curvatura M-φ para varios
  niveles de carga axial y primeros puntos de la curva de interacción P-M.

## Estructura

```
data/
  edificio_config.json    geometría/materiales (Semana 2, no modificar)
  casos_lab3.json         parámetros del lab (q_Q, sismo, lambdas, HA)
src/
  edificio.py             modelo 3D paramétrico (heredado de Semana 2)
  geometria.py            áreas, CM, masas sísmicas y fuerzas laterales
  casos.py                casos base G/Q/EX/EY y verificaciones
  superposicion.py        Parte C: combinar y comparar
  capacidad.py            Parte D: Fiber Section, M-φ y P-M
  run_all.py              orquestador
tests/
  test_lab3.py            pruebas (pytest)
results/                  figuras y CSV generados
Informe_P1L3_Grupo4.md    informe de entrega
```

## Uso

```powershell
.venv\Scripts\Activate.ps1
python -m src.run_all          # ejecuta Partes A-D y genera figuras/CSV
pytest                        # pruebas
```

Todos los parámetros del lab están en `data/casos_lab3.json` (el profesor
puede cambiarlos sin tocar código): intensidad `q_Q`, `fraccion_g` del sismo
(o fuerzas manuales), `lambdas` de superposición y la sección HA de la
Parte D.

## Convenciones

- Unidades SI coherentes: m, kN, kN/m², kN·m. El JSON `casos_lab3.json`
  declara las unidades de cada bloque.
- OpenSees: compresión/curvatura negativas. En Parte D `P < 0 = compresión`
  (convención de OpenSees, igual que el ejemplo oficial de M-φ).
- Cada caso base (G, Q, EX, EY) se construye como modelo independiente
  (determinístico, tags idénticos entre corridas).

## Reproducibilidad

Repositorio de Grupo 4. Commit/tag de entrega indicado en el informe.
`.gitignore` excluye `.venv/`, `results/` y cachés; el código y el informe
son la fuente de verdad.