# P1L3 - Grupo 4

Laboratorio P1L3 (carga viva Q, sismo pseudoestático EX/EY, superposición
y capacidad HA con Fiber Section) sobre el modelo 3D OpenSeesPy del edificio
U. de los Andes 2017_67.

## Comandos

- Ejecutar el pipeline completo (Partes A-D): `python -m src.run_all`
- Pruebas: `pytest`
- Ambiente: usar el venv `.venv` (openseespy 3.8.x). Nunca usar paquetes
  globales o venvs rotos (P1L1-lenz-mauricio, prueba-clon).

## Convenciones críticas

- Unidades SI coherentes: m, kN, kN/m², kN·m. No mezclar m/mm ni kN/N.
- En Parte D, OpenSees escribe compresión como negativo (`Concrete01`, carga
  axial y curvatura); `moment_curvature(P)` con `P < 0` = compresión.
- Los parámetros del lab van en `data/casos_lab3.json` (editables sin tocar
  código). `data/edificio_config.json` es geometría de Semana 2.
- Cada caso base reconstruye un modelo independiente (mismo tag layer);
  `constrain(... Transformation)` + solver estático lineal para el edificio.
- Los resultados y figuras se generan en `results/` (gitignoreado); commitear
  siempre el informe Markdown y el código.
- No editar `src/edificio.py` salvo necesidad real (heredado de Semana 2).
- Tareas heredadas: 10 puntos, entrega jue 10/09/2026 10:30.