# Cuadrados en triángulos

Búsqueda numérica de mejores empaquetamientos de n cuadrados unitarios dentro
de un triángulo equilátero (tabla "Squares in Triangles" del Packing Center de
Erich Friedman). Para cada n, intenta meter los cuadrados en un triángulo
0.0003 más pequeño que el récord actual.

- `packer.c`: optimizador (energía de penetración + FIRE + basin hopping con lado fijo).
- `verify.py`: verificación exacta con aritmética racional.
- `init/`: configuraciones récord digitalizadas, usadas como punto de partida.
- `.github/workflows/buscar.yml`: lanza una máquina por récord (pestaña Actions → "Buscar récords" → Run workflow).
