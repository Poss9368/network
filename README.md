# Expansión y análisis de redes elásticas desordenadas

Este proyecto genera redes bidimensionales periódicas, las expande por pasos y
relaja las posiciones con FIRE. Permite visualizar sus estados, crear GIF y
comparar la fuerza media frente al strain entre distintas simulaciones.

El flujo es **C → Rust → Python**:

1. `networks/networkBuilder.c` genera las redes iniciales.
2. Rust impone una expansión de la caja en ambas direcciones, relaja cada estado y
   guarda un archivo `summary`, con detalles de posiciones opcionales.
3. `plot_box.py` dibuja los estados y puede crear un GIF; `plot_resumen.py` promedia
   los resúmenes, exporta tablas `mean` y grafica fuerza media frente al strain.

## Requisitos

- Compilador de C como `cc` (Clang o GCC), si vas a generar redes nuevas.
- Rust y Cargo con soporte para la edición 2024. Rust no usa dependencias externas.
- Python 3 con NumPy, Matplotlib y Pillow:

```bash
python3 -m pip install numpy matplotlib Pillow
```

Los comandos siguientes se ejecutan desde la raíz del proyecto, salvo que se
indique otra carpeta:

```bash
cd /Users/carlos/Documents/codigos_universidad/network
```

## Guía rápida

Si ya tienes `networks/network_4096_4.100_0000.dat`, puedes ejecutar:

```bash
# Simular 10 puntos desde L = 1.001 hasta L = 1.1.
cargo run --release -- networks/network_4096_4.100_0000.dat 1.1 10 1.001

# Dibujar los estados y crear un GIF de 200 ms por paso.
python3 plot_box.py network_4096_4.100_0000 --gif

# Promediar todas las corridas disponibles con N = 4096 y z = 4.1.
python3 plot_resumen.py --sizes 4096 --z 4.1
```

El último comando requiere que las corridas de ese grupo tengan los mismos puntos
de strain. Con una sola corrida también funciona, sin banda de dispersión.

Si solo necesitas fuerza frente a strain, puedes omitir los detalles de posiciones:

```bash
cargo run --release -- networks/network_4096_4.100_0000.dat 1.1 10 1.001 --no-snapshots
python3 plot_resumen.py --sizes 4096 --z 4.1
```

## Carpetas y archivos de salida

Ejemplo para el sistema `network_4096_4.100_0000`:

```text
networks/
└── network_4096_4.100_0000.dat                  # Red inicial
output/
├── network_4096_4.100_0000_summary.dat          # Largo y presión por paso
├── steps_details/
│   ├── network_4096_4.100_0000_step_0001.dat    # Posiciones y contactos
│   └── ...
├── plots/
│   ├── network_4096_4.100_0000_step_0001.png
│   ├── ...
│   └── network_4096_4.100_0000.gif             # Con --gif
└── processed/
    ├── network_4096_4.100_mean.dat             # Promedios del grupo (N, z)
    └── force_vs_strain.png                    # Comparación de grupos
```

Rust guarda siempre en el `output/` de la raíz del proyecto, aunque se invoque
Cargo desde otra carpeta. Los scripts Python también resuelven sus rutas
predeterminadas respecto de su ubicación.

Las ejecuciones reemplazan archivos con los mismos nombres. No eliminan estados
sobrantes de corridas anteriores: si reduces la cantidad de pasos o usas
`--no-snapshots`, los estados antiguos pueden seguir presentes y `plot_box.py`
puede incluirlos. Usa una carpeta de entrada con los estados que quieras visualizar.
El nombre predeterminado `force_vs_strain.png` es común a todas las selecciones;
usa `--save` si quieres conservar varios gráficos.

## 1. Generar redes iniciales con C

Compila y ejecuta el generador **dentro de `networks/`**, porque escribe los archivos
en el directorio desde el que se ejecuta:

```bash
cd networks
cc -O3 -std=c11 networkBuilder.c -o networkBuilder -lm
./networkBuilder 0 4.1 4
cd ..
```

```text
./networkBuilder <serial> <z> <cantidad_de_redes>
```

| Argumento | Significado en el ejemplo |
| --- | --- |
| `serial = 0` | Identifica la tanda y determina la semilla aleatoria: `1000 + serial`. |
| `z = 4.1` | Coordinación media objetivo: vecinos por partícula. |
| `cantidad_de_redes = 4` | Genera cuatro redes. |

Con la configuración actual de 4096 partículas, el ejemplo genera
`network_4096_4.100_0000.dat` hasta `network_4096_4.100_0003.dat`.
El índice es `serial * cantidad_de_redes + k`, comenzando con `k = 0`.
Otra tanda con `./networkBuilder 1 4.1 4` genera los índices `0004` a `0007`.
Repetir los argumentos reproduce la tanda y sobrescribe sus archivos.
La generación incluye simulación y minimización, por lo que puede tardar.

Para cambiar el tamaño, modifica **ambas** constantes en `networks/networkBuilder.c`
y recompila:

```c
#define N        4096
#define DOUBLE_N 4096.0
```

`N` debe ser un cuadrado perfecto, como 1024, 4096 o 10000. La coordinación obtenida
puede quedar ligeramente por debajo del objetivo porque se eliminan enlaces de
uno en uno.

## 2. Expandir y relajar con Rust

```text
cargo run --release -- <archivo.dat> <largo_final> <pasos> [largo_inicial] [--no-snapshots]
```

`--release` activa las optimizaciones. El separador `--` indica que los argumentos
siguientes pertenecen al programa.

| Argumento u opción | Comportamiento |
| --- | --- |
| `archivo.dat` | Ruta de una red inicial. Un nombre sin ruta también se busca en `networks/`. |
| `largo_final` | Lado final de la caja, finito y mayor que 1. |
| `pasos` | Cantidad total de puntos simulados, entero mayor o igual que 1. |
| `largo_inicial` | Opcional: lado de la caja en el primer punto simulado; debe ser finito y mayor que 1. |
| `--no-snapshots` | Guarda solo el `summary`, sin generar detalles en `output/steps_details/`. Mantiene los mensajes de progreso. |

Por defecto se guardan tanto el resumen como los estados. El archivo de entrada
no se modifica. Los estados se numeran desde `0001`.

### Primer punto y espaciado

La red se carga en una caja de lado `L0 = 1`. El strain es:

```text
ε = (L - L0) / L0 = L - 1
```

`largo_inicial` es el lado del **primer estado expandido**, no el tamaño de la red
al cargarla. Con más de un paso debe ser menor que el largo final. Ambos extremos
se incluyen en el total de puntos.

Con `largo_inicial = 1.001`, `largo_final = 1.1` y tres pasos, los lados son:

```text
1.001 → 1.01 → 1.1
```

El espaciado es uniforme en `log10(ε)`, no en `L`. Si omites `largo_inicial`, se usa
`ε_inicial = ε_final / 1000`, un rango de tres décadas definido por
`STRAIN_DECADES` en `src/main.rs`. Para más de un paso:

```text
ε_final = L_final - 1
ε_inicial = L_inicial - 1          # O ε_final / 1000 si se omite
ε_step = ε_inicial * (ε_final / ε_inicial)^((step - 1) / (pasos - 1))
L_step = 1 + ε_step
```

Con un solo paso se simula directamente `L_final`; si indicas `largo_inicial`, debe
coincidir con él. Se rechazan rangos cuyos puntos no sean distinguibles y
estrictamente crecientes con precisión `f64`.

### Qué minimiza FIRE y cómo se obtiene la presión

En cada paso Rust:

1. Expande de forma afín la caja y las posiciones hasta el lado objetivo.
2. Mantiene fija la caja y relaja las posiciones con FIRE.
3. Calcula la presión del estado relajado y guarda las salidas habilitadas.

**Se impone la expansión; la presión es un resultado.** No se configura una presión
objetivo ni se modifica la caja durante FIRE. Los contactos y sus longitudes de
reposo se conservan. El potencial es elástico con `k = 1.0`.

FIRE converge cuando la mayor norma de fuerza neta por partícula es como máximo
`1e-8`, con un límite de 100000 iteraciones por relajación. Estos parámetros están
en `Fire::default()`, en `src/fire.rs`; no son argumentos de la terminal.
Fuerza neta casi cero por partícula puede coexistir con presión y energía no nulas.
La búsqueda es local y no certifica un mínimo global.

Si FIRE falla o no converge, el programa termina con error y no guarda ese paso.
El `summary` conserva los pasos ya completados.

## 3. Dibujar estados y crear un GIF

```bash
python3 plot_box.py network_4096_4.100_0000
python3 plot_box.py network_4096_4.100_0000 --gif --gif-duration 500
```

Lee los estados de `output/steps_details/` en orden numérico y genera un PNG por
paso en `output/plots/`. Para visualizar una corrida necesitas haber guardado sus
estados con Rust, sin `--no-snapshots`.

| Opción | Uso y valor predeterminado |
| --- | --- |
| `--step 10` | Dibuja únicamente el paso 10. Sin esta opción dibuja todos los encontrados. |
| `--dpi 240` | Resolución de los PNG; por defecto 180, debe ser positiva. |
| `--gif` | Crea también `<sistema>.gif` junto a los PNG. |
| `--gif-duration 500` | Milisegundos por cuadro; por defecto 200. Debe ser un múltiplo positivo de 10. |
| `--output-dir ruta` | Carpeta de entrada de estados; por defecto `output/steps_details/`. |
| `--plot-dir ruta` | Carpeta de salida de PNG y GIF; por defecto `output/plots/`. |

Si cambias `--output-dir` por otra carpeta y no indicas `--plot-dir`, las imágenes
se guardan en `plots/` dentro de esa carpeta de entrada.

El GIF se repite continuamente y usa solo los PNG generados en esa ejecución,
respetando `--step`. Todos los pasos tienen igual duración, independientemente del
strain; la animación no representa tiempo físico. Sin `--gif` se generan solo PNG.

### Cómo interpretar las imágenes

Las posiciones se normalizan con el lado de cada estado:
`x_plot = x/L - 0.5`, `y_plot = y/L - 0.5`. Todas las figuras muestran una caja
unitaria centrada en el origen, sin ejes, título ni leyenda. Así se elimina la
expansión afín y se observan los movimientos relativos durante la relajación.

Las partículas y los contactos son azules. Los enlaces siguen la imagen periódica
más cercana y se cortan y continúan al cruzar los bordes. Los puntos rojos marcan
cruces geométricos entre enlaces sin partículas compartidas; no crean contactos
nuevos en la simulación. Los cruces múltiples se marcan una sola vez y los
segmentos colineales superpuestos no se marcan como intersecciones puntuales.

## 4. Promediar corridas y graficar fuerza frente a strain

```bash
python3 plot_resumen.py --sizes 4096 --z 4.1
python3 plot_resumen.py --sizes 1024 4096 --z 4.1 4.3
```

El script busca `network_<N>_<z>_<id>_summary.dat` en `output/` y produce una curva
por combinación de tamaño y coordinación. Promedia las presiones con igual peso
por corrida, sin mezclar grupos. También acepta el sufijo antiguo `_resumen.dat`;
si existen ambos para la misma corrida, da prioridad a `_summary.dat`.

Los resúmenes de un grupo deben tener los mismos puntos de strain, admitiendo solo
pequeñas diferencias numéricas. Usa los mismos largos inicial y final y cantidad
de pasos. No se interpola. Las combinaciones sin archivos se omiten con un aviso;
si no hay datos para ninguna, el programa termina con error.

### Fuerza representada

El gráfico utiliza **`force_mean` frente a strain**, con ambos ejes logarítmicos:

```text
L = 1 + strain
force_mean = -pressure_mean * L
force_std = pressure_std * L
```

`force_mean` es la fuerza normal media equivalente por lado de la caja, positiva
en tracción y negativa en compresión. En 2D, la presión tiene unidades de
fuerza/longitud, por eso se multiplica por el lado actual. Con condiciones
periódicas no hay paredes físicas: es una fuerza equivalente del esfuerzo global.
La expansión es en ambas direcciones, por lo que no equivale a un ensayo uniaxial.

La banda representa ±1 desviación estándar muestral de la **fuerza**, no un intervalo
de confianza ni el error de la media. Con una sola corrida se omite la banda y se
exporta `pressure_std = 0` como convención.

Los puntos con fuerza media no positiva o no finita se omiten dejando huecos; no
se usa valor absoluto. La banda se omite donde sus límites no sean positivos y
finitos. Si ningún grupo tiene puntos representables, se informa un error.

| Opción | Uso y valor predeterminado |
| --- | --- |
| `--sizes 1024 4096` | Tamaños a comparar; obligatorio. Alias: `--tamanos`. |
| `--z 4.1 4.3` | Coordinaciones objetivo; obligatorio. Alias: `--coordinations`. |
| `--output-dir ruta` | Entrada de resúmenes; por defecto `output/`. |
| `--save ruta.png` | Ruta del gráfico; por defecto `output/processed/force_vs_strain.png`. |
| `--linear-x` | Usa strain lineal; la fuerza sigue en escala logarítmica. |
| `--dpi 240` | Resolución del gráfico; por defecto 180, debe ser positiva. |

Cambiar `--output-dir` solo cambia la entrada: la salida predeterminada sigue en
`output/processed/`. Con `--save`, las tablas `mean` quedan junto al gráfico.

## Formatos de datos

### Red inicial: `networks/<sistema>.dat`

Coordenadas normalizadas entre 0 y 1 y vecinos con índices desde cero; `-1` es
relleno. La cabecera comienza con `0.0` como marcador. Rust la omite y asume una
caja inicial de lado 1.

### Estado: `output/steps_details/<sistema>_step_XXXX.dat`

La primera columna de la cabecera contiene el lado actual `L`; las otras nueve
contienen `-1.1`. Después hay una fila por partícula: `x`, `y` y los índices de sus
vecinos, con relleno `-1` hasta ocho vecinos, sin truncar si hay más.
Las coordenadas son las actuales en `[0, L)`, sin normalizarlas.

Los estados conservan los contactos, pero no guardan longitudes de reposo ni el
potencial. Sirven para visualizar; no permiten reiniciar por sí solos la simulación.
`Network::from_file` está reservado para las redes iniciales de lado 1.

### Resumen: `output/<sistema>_summary.dat`

```text
# box_length pressure
```

Dos columnas numéricas separadas por tabulaciones, con una fila por paso completado:
el lado de la caja y la presión media calculada **después de FIRE**. La presión es
negativa en tracción y positiva en compresión.

### Promedios: `output/processed/network_<N>_<z>_mean.dat`

```text
# strain pressure_mean pressure_std simulations force_mean
```

| Columna | Contenido |
| --- | --- |
| `strain` | `L - 1`, porque el lado de referencia es 1. |
| `pressure_mean` | Presión promedio entre corridas, con su signo original. |
| `pressure_std` | Desviación estándar muestral de la presión; cero con una sola corrida. |
| `simulations` | Cantidad de corridas del grupo. |
| `force_mean` | Fuerza equivalente media por lado: `-pressure_mean * (1 + strain)`. |

La fuerza graficada y la quinta columna usan el mismo cálculo. La desviación
estándar de la fuerza se calcula para la banda y no se exporta como columna.

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
módulo y se encarga de los argumentos, la expansión por pasos y las salidas.

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
let mut network = Network::from_file("networks/network_4096_4.100_0000.dat", potential)?;
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
con potencia no positiva (`velocidad · fuerza <= 0`). Las posiciones se envuelven periódicamente en la caja.

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

El signo de `mean_wall_force()` es el de la presión. Para obtener la convención
positiva en tracción que usa Python:

```rust
let tensile_force = -network.mean_wall_force()?; // -pressure * L
```

## Comprobación

Desde la raíz del proyecto:

```bash
cargo fmt --check
cargo test
cargo clippy --all-targets -- -D warnings
python3 -m unittest test_plot_box.py test_plot_resumen.py
```

Las pruebas cubren la lectura de redes, las fuerzas y la presión, la relajación
FIRE, el rango de expansión, las carpetas de salida, `--no-snapshots`, los gráficos,
el orden de los cuadros del GIF y los promedios de presión y fuerza.
