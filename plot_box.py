"""Dibuja cada red en una caja unitaria centrada y marca cruces entre enlaces.

Uso: python3 plot_box.py network_4096_4.300_0000
Dependencias: numpy y matplotlib.
"""

import argparse
import math
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np


def find_steps(system, output_dir):
    """Busca solo snapshots del sistema y los ordena por número de paso."""
    name = Path(system).name
    name = re.sub(r"_step_\d+\.dat$", "", name).removesuffix(".dat")
    pattern = re.compile(rf"{re.escape(name)}_step_(\d+)\.dat")
    steps = []
    for path in output_dir.glob("*_step_*.dat"):
        match = pattern.fullmatch(path.name)
        if match:
            steps.append((int(match[1]), path))
    if not steps:
        raise ValueError(f"No hay archivos de pasos para {name!r} en {output_dir}")
    return name, sorted(steps)


def read_step(path):
    """Lee L de la cabecera, coordenadas actuales y contactos únicos."""
    with path.open(encoding="utf-8") as stream:
        header = stream.readline().split()
        if not header:
            raise ValueError(f"{path}: falta la cabecera")
        length = float(header[0])
        if not math.isfinite(length) or length <= 0:
            raise ValueError(f"{path}: el largo de la caja debe ser positivo y finito")
        positions, neighbors = [], []
        for line_number, line in enumerate(stream, start=2):
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: faltan las coordenadas")
            position = tuple(map(float, fields[:2]))
            if any(not math.isfinite(x) or not 0 <= x <= length for x in position):
                raise ValueError(f"{path}:{line_number}: posición fuera de la caja")
            adjacent = [int(value) for value in fields[2:] if value != "-1"]
            if len(set(adjacent)) != len(adjacent):
                raise ValueError(f"{path}:{line_number}: vecinos duplicados")
            positions.append(position)
            neighbors.append(set(adjacent))
    if not positions:
        raise ValueError(f"{path}: no hay partículas")
    contacts = []
    for i, adjacent in enumerate(neighbors):
        for j in sorted(adjacent):
            if j == i or not 0 <= j < len(positions):
                raise ValueError(f"{path}: contacto inválido {i}–{j}")
            if i not in neighbors[j]:
                raise ValueError(f"{path}: contacto no recíproco {i}–{j}")
            if i < j:
                contacts.append((i, j))
    return length, np.asarray(positions), contacts


def split_periodic_contact(start, end, length):
    """Devuelve segmentos dentro de [0,L]² siguiendo la imagen mínima.

    Sigue la imagen mínima. Si cruza dos bordes en momentos diferentes, produce
    tres segmentos; si pasa por una esquina, produce dos.
    """
    delta = []
    for a, b in zip(start, end):
        d = b - a
        # Redondeo de empates como f64::round en Rust.
        d -= length * math.copysign(math.floor(abs(d / length) + 0.5), d)
        delta.append(d)

    times = [0.0, 1.0]
    for axis in range(2):
        if delta[axis] != 0:
            for border in (0.0, length):
                t = (border - start[axis]) / delta[axis]
                if 0 < t < 1:
                    times.append(t)
    times.sort()
    cuts = [times[0]]
    for t in times[1:]:
        if not math.isclose(t, cuts[-1], rel_tol=0, abs_tol=1e-12):
            cuts.append(t)

    segments = []
    for lo, hi in zip(cuts, cuts[1:]):
        middle = (lo + hi) / 2
        shift = [length * math.floor((start[a] + middle * delta[a]) / length) for a in range(2)]
        segment = []
        for t in (lo, hi):
            point = tuple(min(length, max(0.0, start[a] + t * delta[a] - shift[a])) for a in range(2))
            segment.append(point)
        segments.append(tuple(segment))
    return segments


def segment_intersection(first, second):
    """Intersección puntual; segmentos colineales superpuestos no definen un punto."""
    a, b = first
    c, d = second
    r = (b[0] - a[0], b[1] - a[1])
    s = (d[0] - c[0], d[1] - c[1])
    denominator = r[0] * s[1] - r[1] * s[0]
    if abs(denominator) <= 1e-12 * math.hypot(*r) * math.hypot(*s):
        return None
    q = (c[0] - a[0], c[1] - a[1])
    t = (q[0] * s[1] - q[1] * s[0]) / denominator
    u = (q[0] * r[1] - q[1] * r[0]) / denominator
    if -1e-12 <= t <= 1 + 1e-12 and -1e-12 <= u <= 1 + 1e-12:
        t = min(1.0, max(0.0, t))
        return (a[0] + t * r[0], a[1] + t * r[1])
    return None


def find_link_intersections(segments, owners, contacts, length):
    """Busca cruces entre contactos distintos que no compartan una partícula.

    Una cuadrícula reduce las comparaciones a segmentos de celdas cercanas.
    owners identifica el contacto original de cada tramo periódico.
    """
    if not segments:
        return []
    grid_size = max(1, int(math.sqrt(len(segments))))
    grid, points = {}, {}
    padding = length * 1e-12

    def cell(value):
        return min(grid_size - 1, max(0, math.floor(value / length * grid_size)))

    for index, segment in enumerate(segments):
        a, b = segment
        cells = [(x, y)
                 for x in range(cell(min(a[0], b[0]) - padding), cell(max(a[0], b[0]) + padding) + 1)
                 for y in range(cell(min(a[1], b[1]) - padding), cell(max(a[1], b[1]) + padding) + 1)]
        candidates = set()
        for key in cells:
            candidates.update(grid.get(key, ()))
        i, j = contacts[owners[index]]
        for other in sorted(candidates):
            if owners[index] == owners[other]:
                continue
            if i in contacts[owners[other]] or j in contacts[owners[other]]:
                continue
            point = segment_intersection(segment, segments[other])
            if point is not None:
                # Un cruce de tres o más enlaces se dibuja una sola vez.
                key = tuple(round(coordinate / length, 12) for coordinate in point)
                points.setdefault(key, point)
        for key in cells:
            grid.setdefault(key, []).append(index)
    return list(points.values())


def network_geometry(positions, contacts, length):
    segments, owners = [], []
    for index, (i, j) in enumerate(contacts):
        parts = split_periodic_contact(positions[i], positions[j], length)
        segments.extend(parts)
        owners.extend([index] * len(parts))
    return segments, find_link_intersections(segments, owners, contacts, length)


def plot_step(path, destination, dpi):
    length, positions, contacts = read_step(path)
    # Coordenadas reducidas: la expansión afín pura queda eliminada del dibujo.
    positions = positions / length
    segments, crossings = network_geometry(positions, contacts, 1.0)
    half = 0.5
    positions = positions - half
    segments = np.asarray(segments) - half

    blue, red = "#2878B5", "#E34A33"
    fig, ax = plt.subplots(figsize=(8, 8), layout="constrained")
    ax.add_collection(LineCollection(segments, colors=blue, linewidths=0.45, alpha=0.75))
    ax.scatter(positions[:, 0], positions[:, 1], s=2, color=blue, linewidths=0, zorder=3)
    ax.add_patch(Rectangle((-half, -half), 1.0, 1.0, fill=False, edgecolor="#222222", linewidth=1.2))
    if crossings:
        points = np.asarray(crossings) - half
        ax.scatter(points[:, 0], points[:, 1], s=16, color=red, linewidths=0, zorder=4)
    limit = 0.54
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit))
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(destination, dpi=dpi)
    plt.close(fig)
    return len(contacts), len(crossings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("system", help="Nombre del sistema, con o sin .dat")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument("--plot-dir", type=Path, help="Destino de los PNG; por defecto output/plots")
    parser.add_argument("--step", type=int, help="Dibujar solo este paso en la misma caja unitaria")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi debe ser positivo")
    try:
        system, steps = find_steps(args.system, args.output_dir)
        if args.step is not None:
            steps = [(step, path) for step, path in steps if step == args.step]
            if not steps:
                raise ValueError(f"No existe el paso {args.step} para {system}")
        plot_dir = args.plot_dir or args.output_dir / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)
        for step, path in steps:
            destination = plot_dir / f"{system}_step_{step:04d}.png"
            contacts, crossings = plot_step(path, destination, args.dpi)
            print(f"{destination} ({contacts} contactos, {crossings} cruces entre enlaces)")
    except (OSError, ValueError, IndexError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
