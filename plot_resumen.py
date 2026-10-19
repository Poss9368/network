"""Promedia resúmenes por número de partículas y coordinación y grafica fuerza media vs strain en log-log.

Ejemplo: python3 plot_resumen.py --sizes 1024 4096 --z 4.1 4.3
Dependencias: numpy y matplotlib.
Lee resúmenes de output/ y guarda el gráfico y las tablas en output/processed/.
"""

import argparse
import math
from pathlib import Path
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def find_summaries(output_dir, size, coordination):
    """Selecciona por N y z del nombre, ignorando snapshots y otros resultados."""
    pattern = re.compile(r"network_(\d+)_(\d+(?:\.\d+)?)_(\d+)_(?:summary|resumen)\.dat")
    matches = {}
    # Los nombres nuevos tienen prioridad para no contar dos veces una corrida.
    candidates = list(output_dir.glob("*_resumen.dat")) + list(output_dir.glob("*_summary.dat"))
    for path in candidates:
        match = pattern.fullmatch(path.name)
        if match and int(match[1]) == size and math.isclose(
            float(match[2]), coordination, rel_tol=0, abs_tol=1e-9
        ):
            matches[int(match[3])] = path
    return [path for _, path in sorted(matches.items())]


def read_summary(path):
    """Convierte el largo a strain ε = L - 1 y conserva la presión con signo."""
    with path.open(encoding="utf-8") as stream:
        rows = [line for line in stream if line.strip() and not line.lstrip().startswith("#")]
    if not rows:
        raise ValueError(f"{path}: el resumen no tiene datos")
    try:
        data = np.loadtxt(rows, ndmin=2)
    except ValueError as error:
        raise ValueError(f"{path}: datos numéricos inválidos: {error}") from error
    if data.shape[1] != 2 or not np.isfinite(data).all():
        raise ValueError(f"{path}: se requieren dos columnas finitas: largo y presión")
    strain = data[:, 0] - 1.0
    if np.any(strain <= 0) or np.any(np.diff(strain) <= 0):
        raise ValueError(f"{path}: los largos deben ser mayores que 1 y estrictamente crecientes")
    return strain, data[:, 1]


def average_summaries(paths):
    """Media y desviación estándar muestral, dando igual peso a cada simulación."""
    if not paths:
        raise ValueError("No hay resúmenes para promediar")
    strain, first = read_summary(paths[0])
    pressures = [first]
    for path in paths[1:]:
        other_strain, pressure = read_summary(path)
        if strain.shape != other_strain.shape or not np.allclose(
            strain, other_strain, rtol=1e-10, atol=1e-14
        ):
            raise ValueError(
                f"{path}: sus valores de strain no coinciden con {paths[0]}. "
                "Usa corridas con los mismos pasos y largo final; no se interpola automáticamente."
            )
        pressures.append(pressure)
    values = np.stack(pressures)
    mean = values.mean(axis=0)
    std = values.std(axis=0, ddof=1) if len(paths) > 1 else np.zeros_like(mean)
    return strain, mean, std


def force_statistics(strain, pressure_mean, pressure_std):
    """Fuerza equivalente por lado y su dispersión: F = −P L, con L0 = 1."""
    length = 1.0 + strain
    return -pressure_mean * length, pressure_std * length


def plot_averages(groups, destination, linear_x=False, dpi=180):
    """Grafica fuerza media en escala logarítmica; la banda es ±1 desviación estándar."""
    fig, ax = plt.subplots(figsize=(8, 5.5), layout="constrained")
    plotted = False
    for size, coordination, paths, strain, mean, std in groups:
        force_mean, force_std = force_statistics(strain, mean, std)
        valid = np.isfinite(force_mean) & (force_mean > 0)
        if not valid.all():
            print(f"N={size}, z={coordination:g}: se omiten {(~valid).sum()} puntos con fuerza media no positiva o no finita en la escala logarítmica.", file=sys.stderr)
        if not valid.any():
            continue
        plotted = True
        line, = ax.plot(strain, np.ma.masked_where(~valid, force_mean), '-o', linewidth=1.8,
                        label=f"N={size}, z={coordination:g} (n={len(paths)})")
        if len(paths) > 1:
            lower, upper = force_mean - force_std, force_mean + force_std
            band = valid & np.isfinite(lower) & np.isfinite(upper) & (lower > 0)
            if np.any(valid & ~band):
                print(f"N={size}, z={coordination:g}: se omite la banda donde sus límites no son positivos y finitos.", file=sys.stderr)
            ax.fill_between(strain, lower, upper, where=band, color=line.get_color(), alpha=0.18)
    if not plotted:
        plt.close(fig)
        raise ValueError("No hay valores con fuerza media > 0 para generar el gráfico logarítmico")
    ax.set_xscale("linear" if linear_x else "log")
    ax.set_yscale("log")
    ax.set_xlabel("Strain ε = L − 1")
    ax.set_ylabel("Mean tensile force")
    ax.grid(alpha=0.2)
    ax.legend(title="Media; banda ±1 desviación estándar", fontsize=9, title_fontsize=9)
    fig.savefig(destination, dpi=dpi)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", "--tamanos", type=int, nargs="+", required=True,
                        help="Números de partículas que se quieren comparar")
    parser.add_argument("--z", "--coordinations", type=float, nargs="+", required=True,
                        help="Coordinaciones objetivo indicadas en los nombres de los archivos")
    project_dir = Path(__file__).resolve().parent
    parser.add_argument("--output-dir", type=Path, default=project_dir / "output",
                        help="Carpeta de entrada de los resúmenes; por defecto output/ junto al script")
    parser.add_argument("--save", type=Path, help="Gráfico de salida; por defecto output/processed/force_vs_strain.png junto al script")
    parser.add_argument("--linear-x", action="store_true", help="Usar strain lineal en lugar de logarítmico")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    if any(size <= 0 for size in args.sizes):
        parser.error("Los tamaños deben ser enteros positivos")
    if any(not math.isfinite(z) or z <= 0 for z in args.z):
        parser.error("Las coordinaciones deben ser positivas y finitas")
    if args.dpi <= 0:
        parser.error("--dpi debe ser positivo")
    try:
        groups = []
        for size in dict.fromkeys(args.sizes):
            for coordination in dict.fromkeys(args.z):
                paths = find_summaries(args.output_dir, size, coordination)
                if not paths:
                    print(f"Sin resúmenes para N={size}, z={coordination:g}; se omite esa combinación.", file=sys.stderr)
                    continue
                strain, mean, std = average_summaries(paths)
                groups.append((size, coordination, paths, strain, mean, std))
        if not groups:
            raise ValueError(f"No hay resúmenes para las combinaciones solicitadas en {args.output_dir}")

        destination = args.save or project_dir / "output" / "processed" / "force_vs_strain.png"
        destination.parent.mkdir(parents=True, exist_ok=True)
        for size, coordination, paths, strain, mean, std in groups:
            table = destination.parent / f"network_{size}_{coordination:.3f}_mean.dat"
            # L0 = 1: fuerza normal media equivalente por lado, positiva en tracción.
            force_mean, _ = force_statistics(strain, mean, std)
            np.savetxt(table, np.column_stack((strain, mean, std, np.full(len(strain), len(paths)), force_mean)),
                       header="strain pressure_mean pressure_std simulations force_mean",
                       fmt=["%.17e", "%.17e", "%.17e", "%d", "%.17e"], delimiter="\t")
            print(f"N={size}, z={coordination:g}: {len(paths)} simulaciones, {len(strain)} puntos → {table}")
        plot_averages(groups, destination, args.linear_x, args.dpi)
        print(f"Gráfico: {destination}")
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
