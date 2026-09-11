# P1L3 — Laboratorio: Carga viva, sismo pseudoestático, superposición y capacidad HA

**Edificio U. de los Andes — Proyecto 2017_67 · Grupo 4**
Repositorio: `github.com/mauricio-lenz/P1L3_Grupo4` — tag de entrega **`P1L3-entrega`**

- **Unidades:** SI coherente (m, kN, kN/m², kN·m).
- **Modelo:** misma geometría y áreas tributarias de Semana 2 (`data/edificio_config.json`).
- **Parámetros del lab:** `data/casos_lab3.json` (editables por el profesor sin tocar código).
- **Ejecución:** `python -m src.run_all` (figuras en `results/figures/`, CSV en `results/csv/`).

---

## Parte A — Caso de carga viva Q

La misma geometría tributaria del caso G se reutiliza con la intensidad de
carga viva de uso `q_Q` (NCh1537, oficina/educacional):

| Nivel | `q_Q` [kN/m²] | `q_Q·A` [kN] |
|------:|:---:|------------:|
| Subterráneo–Piso3 | 3.0 | 3.0 × 726.8 = 2180.4 |
| Piso4 (cubierta) | 1.5 | 1.5 × 823.7 = 1235.6 |

Verificación de conservación (transferida por áreas tributarias ↔ reacción basal):

| Cantidad | Valor |
|---------:|------:|
| Q total por áreas (`sum q_Q·A`) | **10247.8 kN** |
| Q transferida a vigas (suma resumen) | **10247.8 kN** |
| Suma de reacciones Fz (base) | **10247.8 kN** |
| Diferencia relativa (áreas) | 0.0 (ok) |
| Diferencia relativa (reacciones) | 0.0 (ok) |

La carga viva se agrega al modelo con el mismo algoritmo de áreas tributarias
de Semana 2 (`src/edificio.areas_tributarias(cfg, q_map)`); el caso Q no
incluye cargas puntuales de equipos (`tag_patron=1`).

## Parte B — Sismo pseudoestático EX y EY

Para cada piso se aplica en el centro de masa (CM) una fuerza lateral
`F_i = 0.20·W_i` donde `W_i = peso propio (losa + estructura) + 0.5·Q_i`
(`fraccion_sobrecarga_masa = 0.50`).

Pesos sísmicos por nivel (`results/csv/pesos_sismicos.csv`):

| Nivel | Área [m²] | PP [kN] | Q [kN] | W = PP + 0.5·Q [kN] | F = 0.20·W [kN] |
|------:|---------:|--------:|-------:|--------------------:|----------------:|
| Subterráneo | 726.8 | 7989.5 | 2180.4 | 9079.7 | 1815.9 |
| Piso1 | 726.8 | 7989.5 | 2180.4 | 9079.7 | 1815.9 |
| Piso2 | 726.8 | 7989.5 | 2180.4 | 9079.7 | 1815.9 |
| Piso3 | 823.7 | 8498.3 | 2471.1 | 9733.8 | 1946.8 |
| Piso4 | 823.7 | 7674.6 | 1235.6 | 8292.3 | 1658.5 |
| **Σ** | | | | **45265.3** | **9053.1** |

Como el CM de piso no coincide con el nodo master del diafragma rígido, cada
fuerza se aplica con su momento torsor (wrench) equivalente en el master
(`src/casos.aplicar_sismo`):

- EX: `(F, 0, 0)` + `M_z = −F·(y_CM − y_master)`
- EY: `(0, F, 0)` + `M_z = F·(x_CM − x_master)`

Verificaciones (`results/csv/cortes_basales.csv`):

| Caso | F lateral [kN] | Corte basal [kN] | δ techo [m] | Sentido | Torción techo Rz [rad] |
|-----:|--------------:|----------------:|------------:|:-------:|-----------------------:|
| EX | 9053.1 | −9053.1 | 0.00393 | ok | −3.6e-5 |
| EY | 9053.1 | −9053.1 | 0.04530 | ok | 2.7e-4 |

El corte basal es exactamente la reacción de equilibrio de la fuerza aplicada
y la deformada adopta el sentido de la carga; existe rotación Rz por piso
(torsión), mayor en EY por la excentricidad del CM en X.

## Parte C — Superposición de casos

Combinación evaluada (λ editables en `casos_lab3.json`):

> **R = 1.0·G + 1.0·Q + 0.9·EX + 0.75·EY**

Método de verificación: se suman linealmente las **respuestas unitarias** de
los cuatro casos (desplazamientos de master, reacciones y fuerzas internas) y
se comparan contra una **corrida explícita** de OpenSees con todas las cargas
combinadas aplicadas de una vez (patrón 1 = gravedad `λG·G + λQ·Q`, patrón 2 =
sismo `λEX·EX + λEY·EY`).

| Variable | Superposición | Explícita | Error relativo |
|----------|---------------|-----------|---------------:|
| Desplazamiento techo (master) | [0.004166, 0.033988, −0.003176, …] | idéntico | 2.2e-15 |
| Reacciones basales totales (global) | [−8147.8, −6789.8, 29520.0, …] kN | idéntico | 2.2e-15 |
| Reacción nodo base 1 (local) | [−6.57, −333.8, −90.4, …] kN | idéntico | 1.4e-15 |
| Fuerza interna columna crítica (tag 87) | [−1315.5, 443.1, 480.7, …] kN/kN·m | idéntico | 7.3e-16 |
| Fuerza interna viga verificación (tag 118) | [−334.8, 136.9, 136.9, …] kN/kN·m | idéntico | 4.4e-15 |

**Error relativo máximo = 4.4e-15** — la superposición reproduce la corrida
explícita con precisión de máquina (el sistema es lineal). Las verificaciones
mínimas pedidas (desplazamiento, reacción y fuerza interna) están cumplidas.

## Parte D — Capacidad de la columna HA con Fiber Section

Columna crítica (I′, Eje 1, Piso1), sección P-70x70 (`casos_lab3.json`):

| Propiedad | Valor |
|----------:|------:|
| b × h [m] | 0.70 × 0.70 (Ag = 0.49 m²) |
| Hormigón | H30, `fc = 25 MPa` (Concrete01) |
| Acero | G420, `fy = 420 MPa`, Es = 200 GPa (Steel01) |
| Recubrimiento libre | 0.05 m |
| Refuerzo | 8φ25 → A_s = 39.27 cm², ρ = 0.80 % |
| EI inicial (nominal, ACI) | 470 196 kN·m² |

**Fiber Section** construida con el comando `section Fiber` de OpenSees
(`src/capacidad.construir_seccion_fibra`): núcleo → `patch rect` 12×12, 4
cintas de recubrimiento `patch`, refuerzo → `layer straight` 3/3/1/1 barras.
Figuras: `results/figures/fig_seccion_fibras.png` (discretización: 192
fibras de concreto + 8 de acero), materiales y convenciones en la figura.

**Diagrama de interacción P–M** con los **5 puntos característicos**,
calculado por **compatibilidad de deformaciones** (deformación de borde
`εc = 0.003`, bloque rectangular de Whitney `a = β1·c`, `β1 = 0.85`) y
validado contra el máximo M de las curvas M-φ de la Fiber Section
(`results/csv/cap_pm.csv`, ![P][p-pos] positivo = compresión):

[p-pos]: # "convención del diagrama"

| Punto | P [kN] | M [kN·m] |
|------:|-------:|---------:|
| Compresión axial pura (M = 0, sobre el eje Y) | +11978 | 0 |
| Control de compresión (εt = 0.002) | +5021 | 1271 |
| Condición balanceada (εt = fy/Es) | +4884 | 1278 |
| Flexión pura (P = 0) | 0 | 512 |
| Tensión axial pura (M = 0, sobre el eje Y) | −1649 | 0 |

- **Compresión axial pura**: `P_o = 0.85·fc'·(Ag − As) + fy·As` (ACI 318).
- **Balanceada**: `c_b = d·εcu/(εcu + εy)`; acero extremo en fluencia
  simultánea con falla de concreto → aquí el diagrama cambia de pendiente
  bruscamente (controlado por compresión arriba, por tracción abajo).
- **Tensión axial pura**: todo el acero fluye a tracción → `P = −fy·As`.
- El primer y el último punto están sobre el eje Y (M = 0), como corresponde
  al diagrama de interacción.

En la curva P–M se marca la **demanda** de la columna crítica en el caso
combinado (Parte C): **P = −1315.5 kN, M = 26.5 kN·m**, holgadamente dentro
de la capacidad (factor de uso ≈ 3 %), coherente con el diseño lineal de
Semana 2. Figuras: `fig_M_phi_columna.png` (todas las curvas M-φ) y
`fig_PM_columna.png` (interacción P–M con la demanda).

## Resumen y conclusiones

1. El caso Q conserva exactamente la carga (áreas tributarias ↔ reacción
   basal), indistinguible del caso G en cuanto a topología de la transferencia.
2. El sismo pseudoestático aplica correctamente F = 0.20·W en los CM con su
   torsión de piso; corte basal, deformada y Rz verificados.
3. La superposición R = G + Q + 0.9EX + 0.75EY es idéntica a la corrida
   explícita (error 4.4e-15), validando el uso de casos unitarios lineales.
4. El diagrama de interacción P–M con los 5 puntos característicos
   (compresión axial, control de compresión, balanceada, flexión pura,
   tensión axial) reproduce la curva de capacidad HA esperada y la demanda
   queda muy por debajo de la capacidad (uso ≈ 3 %).

**Reproducibilidad:** `requirements.txt`, `pytest` (5 pruebas), ejecutar
`python -m src.run_all`. Entrega: tag `P1L3-entrega`.