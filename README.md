# Generar, expandir y visualizar redes

El flujo del proyecto tiene tres etapas:

1. **C:** generar las redes iniciales con `networks/networkBuilder.c`.
2. **Rust:** aplicar una expansión por pasos y relajar cada estado con FIRE.
3. **Python:** dibujar cada estado guardado para comparar el movimiento de las partículas.

## Requisitos y carpeta de trabajo

- Un compilador de C disponible como `cc` (Clang o GCC).
- Rust y Cargo con soporte para la edición 2024. El programa Rust no usa dependencias externas.
- Python 3 con NumPy y Matplotlib para los gráficos:

```bash
python3 -m pip install numpy matplotlib
```

Los ejemplos parten desde la raíz del proyecto:

```bash
cd /Users/cjvillarroel/Desktop/network
```

## 1. Crear las redes iniciales con C

Compila y ejecuta el generador **dentro de `networks/`**, porque guarda los `.dat`
en el directorio desde el que se ejecuta:

```bash
cd networks
cc -O3 -std=c11 networkBuilder.c -o networkBuilder -lm
./networkBuilder 0 4.3 4
cd ..
```

El ejecutable recibe tres argumentos:

```text
./networkBuilder <serial> <z> <cantidad_de_redes>
```

| Argumento del ejemplo | Significado |
| --- | --- |
| `0` | Identificador de la tanda; determina la semilla aleatoria como `1000 + serial`. |
| `4.3` | Coordinación media objetivo: vecinos por partícula. |
| `4` | Cantidad de redes que se generarán. |

Con la configuración actual de **4.096 partículas**, se guardan:

```text
networks/network_4096_4.300_0000.dat
networks/network_4096_4.300_0001.dat
networks/network_4096_4.300_0002.dat
networks/network_4096_4.300_0003.dat
```

La numeración es `serial * cantidad_de_redes + k`, con `k` desde cero. Para otra tanda
con cuatro redes, `./networkBuilder 1 4.3 4` usa otra semilla y genera los índices
`0004` a `0007`. Repetir los mismos argumentos reproduce la tanda y sobrescribe
sus archivos. La generación incluye simulación y minimización, por lo que puede tardar.

El tamaño se configura en `networks/networkBuilder.c` mediante **ambas** constantes:

```c
#define N        4096
#define DOUBLE_N 4096.0
```

Si lo cambias, `N` debe ser un cuadrado perfecto, como `1024`, `4096` o `10000`,
y debes recompilar. Los nombres de salida cambiarán según `N` y `z`.
La coordinación obtenida puede quedar ligeramente por debajo del objetivo porque
los enlaces se eliminan de uno en uno.

Los `.dat` iniciales contienen coordenadas normalizadas entre 0 y 1 y vecinos con
índices desde cero. `-1` es relleno. Su cabecera empieza en `0.0` como marcador;
Rust la omite y asume una caja inicial de **1 × 1**.

## 2. Ejecutar la expansión con Rust

Desde la raíz del proyecto, selecciona **una red**, el largo final y la cantidad de pasos:

```bash
cargo run --release -- network_4096_4.300_0000.dat 2.0 10
```

También puedes indicar la ruta de entrada:

```bash
cargo run --release -- networks/network_4096_4.300_0000.dat 2.0 10
```

La sintaxis es:

```text
cargo run --release -- <archivo.dat> <largo_final> <pasos>
```

| Argumento del ejemplo | Significado |
| --- | --- |
| `network_4096_4.300_0000.dat` | Nombre de una red en `networks/`, o ruta al archivo. |
| `2.0` | Largo final de cada lado: la caja pasa de `1 × 1` a `2 × 2`. Debe ser finito y mayor que 1. |
| `10` | Cantidad de pasos de expansión. Debe ser un entero mayor o igual que 1. |

Cargo compila el proyecto y ejecuta la simulación. En cada paso se expande la caja
y las posiciones de forma afín; después FIRE relaja las posiciones manteniendo
fija la caja de ese paso. Los contactos y sus longitudes de reposo se conservan.
El potencial actual es elástico, con rigidez `k = 1.0`. El archivo inicial no se modifica.

### Espaciado logarítmico y salidas

El strain de expansión es `ε = L - 1`, positivo. Como `log10(0)` no está definido,
se usa por defecto `ε_inicial = ε_final / 1000`: tres décadas espaciadas de forma
uniforme en `log10(ε)`. Esta convención está en `STRAIN_DECADES`, en `src/main.rs`.

```text
ε_final = L_final - 1
ε_step = ε_final * 10^(-3 + 3*(step-1)/(pasos-1))
L_step = 1 + ε_step
```

Con un solo paso se va directamente a `L_final`. Por ejemplo, `L_final = 1.1` y
cuatro pasos produce `1.0001`, `1.001`, `1.01`, `1.1`. Los largos objetivo se validan
para que sean distinguibles y estrictamente crecientes con precisión `f64`.
En cada paso se aplica el factor relativo `L_step / L_actual`, se relaja con FIRE
manteniendo fija esa caja y se guarda el estado convergido. Si FIRE falla o no
converge, se detiene el proceso y ese paso no se agrega a las salidas.

Las salidas siempre están en `output/` en la raíz del proyecto, junto a `Cargo.toml`,
independientemente del directorio desde el que se ejecute Cargo. Por ejemplo, ejecutar
desde `networks/` también guarda en el `output/` de la raíz, no en `networks/output/`:

```text
network_4096_4.300_0000_step_0001.dat
network_4096_4.300_0000_step_0002.dat
...
network_4096_4.300_0000_resumen.dat
```

Cada snapshot conserva el formato de texto separado por tabulaciones. La primera
columna de la cabecera es el largo actual `L`; las otras nueve contienen `-1.1`.
Luego hay una fila por partícula: `x`, `y` y los índices de sus vecinos desde cero,
con relleno `-1` hasta ocho vecinos (sin truncar si hay más). Las coordenadas son las
actuales en `[0, L)`, sin normalizarlas. Los contactos se conservan.
Estos archivos permiten visualizar el proceso; no guardan longitudes de reposo ni
el potencial, por lo que no sirven para reiniciar por sí solos la simulación.
`from_file` sigue reservado para las redes iniciales en una caja de lado 1.

El resumen tiene la cabecera `# box_length pressure` y dos columnas numéricas:
el largo y la presión media **después de FIRE**, una fila por paso completado.
Una ejecución reemplaza los archivos de los mismos nombres; no elimina snapshots
de otros índices que puedan existir de una ejecución anterior.

## 3. Crear los gráficos con Python

```bash
python3 plot_box.py network_4096_4.300_0000
```

Lee los archivos `output/<sistema>_step_*.dat` en orden numérico y guarda un PNG
por paso en `output/plots/`, con el mismo nombre del snapshot. Normaliza las posiciones
con el largo de la cabecera: `x_plot = x/L - 0.5`, `y_plot = y/L - 0.5`.
Todas las figuras muestran una caja de lado 1 centrada en el origen (`-0.5` a `0.5`).
Así se elimina la expansión afín del dibujo y se observan los movimientos relativos
de las partículas durante la relajación. La simulación y los `.dat` no se modifican.
Las figuras no tienen títulos, ejes ni leyenda: solo la caja, la red y sus cruces.

Las partículas y los contactos son azules. Cada contacto se dibuja una sola vez,
siguiendo su imagen periódica más cercana. Los enlaces que salen de la caja se
cortan y continúan en el lado opuesto. Se marcan en rojo las intersecciones
geométricas entre dos o más enlaces que no compartan una partícula. Un cruce
múltiple tiene un único punto; atravesar el borde de la caja por sí solo no genera
marcas. También se detectan cruces entre los tramos de enlaces periódicos.
Los segmentos colineales superpuestos no tienen una intersección puntual y no
se marcan. La búsqueda usa una cuadrícula para evitar comparar todos los pares.

Por ejemplo, después de la expansión de diez pasos anterior, puedes dibujar solo
el último paso o aumentar su resolución:

```bash
python3 plot_box.py network_4096_4.300_0000 --step 10 --dpi 240
```

`--output-dir` permite leer otra carpeta y `--plot-dir` elegir el destino de los PNG.

Para el ejemplo anterior se generan diez imágenes:

```text
output/plots/network_4096_4.300_0000_step_0001.png
...
output/plots/network_4096_4.300_0000_step_0010.png
```

El script toma por defecto `output/` junto a `plot_box.py`, la misma carpeta donde
Rust guarda los resultados. Repetir el comando de gráficos reemplaza los PNG del mismo nombre.
Si cambias la cantidad de pasos de una simulación anterior, revisa los snapshots
sobrantes: el gráfico lee todos los archivos de pasos presentes para ese sistema.

### Promediar resúmenes: −presión frente a strain en log-log

`plot_resumen.py` recibe listas de números de partículas y coordinaciones objetivo.
Busca todos los archivos `output/network_<N>_<z>_<id>_resumen.dat` para cada combinación
y dibuja una curva promedio por `(N, z)`:

```bash
python3 plot_resumen.py --sizes 1024 4096 --z 4.1 4.3
```

Para los resúmenes de `N = 4096`, `z = 4.3` que ya existen:

```bash
python3 plot_resumen.py --sizes 4096 --z 4.3
```

El strain es `ε = L - 1`. Se promedia la presión con igual peso para cada simulación,
sin mezclar tamaños ni coordinaciones. El gráfico muestra **−presión media frente a
strain**, con ambos ejes logarítmicos. En tensión, la presión es negativa y −presión
es positiva. No se usa el valor absoluto. `--linear-x` permite usar strain lineal,
manteniendo el eje de −presión logarítmico.
El gráfico incluye una banda de ±1 desviación estándar muestral
entre simulaciones; no es un intervalo de confianza ni el error estándar de la media.
Con una sola simulación se muestra su curva sin banda y se exporta desviación cero
como convención, porque no se puede estimar la dispersión entre corridas.

Los puntos con −presión <= 0 no pueden representarse en log y se omiten con un aviso;
la curva queda interrumpida en esos puntos. La banda se omite en los intervalos donde
sus límites no sean positivos y finitos, sin cambiar los valores para forzarlos a entrar
en la escala logarítmica. Si no queda ningún punto válido, se informa un error.

Los resúmenes de cada grupo deben contener los mismos valores de strain, admitiendo
solo pequeñas diferencias numéricas. Si difieren los pasos o el rango, se informa un
error: no se interpola ni se promedian presiones de strains diferentes. Las combinaciones
sin archivos se indican y se omiten; si ninguna tiene datos, el programa termina con error.

Se generan en `output/plots/`:

```text
presion_vs_strain.png
network_4096_4.300_promedio.dat
```

Se exporta una tabla por grupo con columnas `strain`, `pressure_mean`, `pressure_std`
y `simulations`. La tabla conserva la presión media con su signo original; el cambio
a −presión se aplica solo al gráfico. `--output-dir` cambia la carpeta de entrada y `--save ruta.png`
elige el gráfico de salida; las tablas quedan junto al gráfico. Las salidas de los
mismos nombres se reemplazan al repetir el comando.

Las pruebas de promedios y agrupamiento se ejecutan con:

```bash
python3 -m unittest test_plot_resumen.py
```

## Referencia del código

### Almacenamiento

- `Network::particles: Vec<Particle>`: cada partícula guarda `position: [f64; 2]`
  y `force: [f64; 2]`. La fuerza se inicializa en `[0.0, 0.0]` al cargar la red.
  Su identificador es su índice en el vector, comenzando en cero.
- `Network::contacts: Vec<Contact>`: cada contacto guarda `particles: [usize; 2]`
  y `initial_length: f64`, además de `contact_force: f64`, inicialmente cero.
  Cada enlace aparece una sola vez, con el índice menor primero.
  Los índices hacen referencia al vector de partículas, sin copiar posiciones.
- `Network::box_length: f64`: lado de la caja cuadrada periódica, inicialmente `1.0`.
- `Network::potential: P`: potencial que implementa el trait `Potential`.

Las estructuras y el lector están en `src/network.rs`. `src/main.rs` declara ese
módulo y se encarga de los argumentos y de mostrar los datos de la red.

El lector omite la cabecera, descarta el relleno `-1` y comprueba los índices,
la reciprocidad y la ausencia de vecinos duplicados o enlaces consigo mismo.
La longitud inicial de cada contacto se calcula en la caja de lado 1 usando la imagen
mínima periódica, incluidos los enlaces que cruzan un borde de la caja.

Con `ElasticPotential`, la longitud inicial es la longitud de reposo del resorte:
`energy = 0.5 * k * (length - initial_length)^2`.
Mientras los contactos usen estos índices, no se deben reordenar ni eliminar partículas
sin actualizar los contactos.

### Uso desde Rust

```rust
use crate::network::Network;
use crate::potential::ElasticPotential;
use crate::fire::Fire;

let potential = ElasticPotential::new(1.0);
let mut network = Network::from_file("networks/network_4096_4.300_0000.dat", potential)?;
network.affine_expand(1.01)?;
let result = Fire::default().minimize(&mut network)?;
if !result.converged {
    // Se agotó el límite de pasos; todavía no se alcanzó la tolerancia.
}
```

`from_file` construye la red a partir del archivo y devuelve errores de lectura o de
formato. `from_reader` permite cargar el mismo texto desde cualquier `BufRead`.
Los nombres de archivo no están restringidos por extensión; el formato esperado es
el del `.dat` del generador.

`affine_expand` escala las posiciones y la caja por el mismo factor. Conserva los
contactos y sus longitudes iniciales. Acepta factores finitos y positivos (menores que
1 para comprimir). `set_box_length` cambia solo la caja, sin escalar las posiciones.
Ambas operaciones devuelven `Result`, por lo que sus errores deben manejarse con `?`
o de forma explícita. Las fuerzas se inicializan en cero al cargar el archivo.
Después de cambiar posiciones o caja, se llama explícitamente a `update_forces()`.

`update_forces` reinicia las fuerzas y recorre cada contacto una sola vez.
Guarda la fuerza escalar `k * (length - initial_length)` en `contact_force`:
positiva en tensión y negativa en compresión. Multiplica por la dirección periódica
de i hacia j, suma la contribución a i y la resta de j. Así acumula las contribuciones
de todos los contactos en cada partícula.

Devuelve un error si un contacto tiene longitud cero o no finita, o si el potencial
produce una fuerza no finita. Si falla, no se deben usar las fuerzas hasta completar
una actualización válida. El potencial elástico supone una rigidez finita y positiva.

### Relajación con FIRE

`src/fire.rs` implementa [FIRE de Bitzek et al.](https://doi.org/10.1103/PhysRevLett.97.170201)
con velocity Verlet y masas unitarias. Las velocidades auxiliares se crean de cero en
cada llamada. Se mezclan con las fuerzas usando las normas globales; el paso de tiempo
crece después de cinco pasos con potencia positiva y se reduce al detener un movimiento
que sube la energía. Las posiciones se envuelven periódicamente en la caja.

Durante FIRE se mantienen fijos el lado de la caja deformada, los contactos y sus
longitudes de reposo. Por defecto usa `force_tolerance = 1e-8`, `max_steps = 100_000`,
`initial_dt = 0.01` y `max_dt = 0.1`, configurables con `Fire { .. }`.

`FireResult` informa `converged`, `steps`, `energy` y `max_force` (la máxima norma
de la fuerza total por partícula). Alcanzar el límite de pasos devuelve
`converged = false`; los errores numéricos devuelven `Err`.
La convergencia depende de fuerzas, no de que la energía sea cero: una red expandida
puede conservar tensiones y energía elástica en equilibrio. Es una búsqueda local;
no certifica el mínimo global ni la estabilidad mediante el Hessiano.

### Presión y fuerzas equivalentes sobre la caja

```rust
let [px, py] = network.normal_pressures()?;
let pressure = network.pressure()?;             // (px + py) / 2
let [fx, fy] = network.wall_forces()?;           // [px * L, py * L]
let mean_force = network.mean_wall_force()?;    // pressure * L
```

Se usa el [virial estático](https://docs.lammps.org/compute_pressure.html), con área
`A = L²`. Para cada contacto único, `d = [dx, dy]` es el desplazamiento periódico,
`l = |d|` y `f = potential.force(l, initial_length) = dU/dl`:

```text
Px = -sum(f * dx² / l) / A
Py = -sum(f * dy² / l) / A
P  = (Px + Py) / 2 = -sum(f * l) / (2*A)
```

La presión es positiva en compresión y negativa en tensión. Las velocidades de FIRE
son auxiliares del minimizador: no aportan presión cinética. Los métodos evalúan el
potencial con la geometría actual, sin depender de fuerzas almacenadas previamente.

En 2D, la presión tiene unidades de fuerza/longitud. Las fuerzas normales equivalentes
por lado son `Fx = Px*L` y `Fy = Py*L`. Su promedio con signo es `P*L`.
El signo positivo indica empuje hacia afuera; el negativo, tracción hacia adentro.
Los lados opuestos tienen direcciones opuestas: no se suman los cuatro lados.
Con periodicidad no hay paredes físicas; estos valores representan las fuerzas
normales medias equivalentes del esfuerzo global, no la fuerza local en un corte.
No incluyen componentes tangenciales de cizalla.

## Comprobación

```bash
cargo test
cargo clippy --all-targets -- -D warnings
```
