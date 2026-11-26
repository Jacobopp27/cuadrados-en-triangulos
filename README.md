# Cuadrados en triángulos

Búsqueda numérica de mejores empaquetamientos de n cuadrados unitarios dentro
de un triángulo equilátero (tabla "Squares in Triangles" del Packing Center de
Erich Friedman). Para cada n, intenta meter los cuadrados en un triángulo
0.0003 más pequeño que el récord actual.

- `packer.c`: optimizador (energía de penetración + FIRE + basin hopping con lado fijo).
- `verify.py`: verificación exacta con aritmética racional.
- `init/`: configuraciones récord digitalizadas, usadas como punto de partida.
- `.github/workflows/buscar.yml`: lanza una máquina por récord (pestaña Actions → "Buscar récords" → Run workflow).

## Polígonos regulares (pentágonos, hexágonos, octágonos)

- `polypack.c`: generalización a un polígono regular de K lados (`./polypack K n segundos semilla s fixed`).
- `verify_poly.py`: superposiciones verificadas con aritmética racional exacta; contención con 80 dígitos y margen de 1e-30.
- `objetivos.txt`: los 120 récords objetivo (K, n, valor publicado, autor).
- `.github/workflows/poligonos.yml`: reparte los objetivos entre varias máquinas (Actions → "Buscar récords en polígonos").

## Pulir desde los récords

- `digitize_poly.py`: lee la imagen de un récord y estima centro y ángulo de cada cuadrado.
- `initpoly/`: las 100 configuraciones récord digitalizadas que la búsqueda a ciegas no alcanzó.
- `objetivos_pulir.txt` y `.github/workflows/pulir.yml`: compresión y basin hopping partiendo de cada récord (Actions → "Pulir récords en polígonos").

## Cuadrículas sin triángulos isósceles (problema 6.59 de AlphaEvolve)

El mayor subconjunto de la cuadrícula n×n sin tres puntos que formen un triángulo isósceles (incluido el caso degenerado del punto medio). Mejores valores de AlphaEvolve: 112 puntos (n=64) y 164 (n=100).

- `isos/isos.c`: búsqueda tabú con conteo incremental de violaciones (mediatrices y círculos), con simetría opcional.
- `isos/verify_isos.py`: verificación exacta con enteros.
- `isos/trabajos.txt` y `.github/workflows/isosceles.yml`: reparto de 80 búsquedas entre 20 máquinas.

## Cociente máx/mín de distancias en 3D

Tabla de Friedman: <https://erich-friedman.github.io/packing/maxmin3/>. Se buscan
n puntos en el espacio que minimicen r² = (distancia máxima / distancia mínima)².

- `mm3/mm3.c`: optimizador (umbral creciente + L-BFGS + saltos entre cuencas).
  Uso: `./mm3 n segundos semilla [inicial|-] [récord_r2]`.
- `mm3/verify_mm3.py`: calcula r² exactamente (aritmética racional) a partir de
  las coordenadas y lo compara con el récord publicado.
- Flujo `maxmin3.yml`: 80 búsquedas para n = 15–40 (n = 31–40 no tienen valor publicado).
