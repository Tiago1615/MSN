from __future__ import annotations
import math
from pathlib import Path
import matplotlib.pyplot as plt

# -------------------------
# Datos del problema
# -------------------------
TA = 50.0
K = 200.0
H_CONV = 60.0
Q_BASE = -26087.0
LE = 0.005

W = 0.025
BASE_H = 0.005
H = 0.025
X_MID = 0.0125

OUTDIR = Path(__file__).resolve().parent
MEF_FILE = OUTDIR / "Puntos_internos_MEF.txt"

# ---------------------------------------------------
# IDs:
# 1 -> laterales (q=0)
# 2 -> convección
# 3 -> base caliente
# ---------------------------------------------------
BOUNDARIES = [
    (1, 0.0),
    (2, H_CONV),
    (1, Q_BASE),
]


def almost(a: float, b: float, tol: float = 1e-12) -> bool:
    return abs(a - b) <= tol


def segment_points(
    p0: tuple[float, float],
    p1: tuple[float, float],
    le: float,
) -> list[tuple[float, float]]:
    x0, y0 = p0
    x1, y1 = p1

    length = math.hypot(x1 - x0, y1 - y0)
    n = max(1, round(length / le))

    return [
        (x0 + (x1 - x0) * i / n, y0 + (y1 - y0) * i / n)
        for i in range(n + 1)
    ]


def boundary_id(p0: tuple[float, float], p1: tuple[float, float]) -> int:
    x0, y0 = p0
    x1, y1 = p1

    # base inferior
    if almost(y0, 0.0) and almost(y1, 0.0):
        return 3

    # laterales exteriores inferiores
    if almost(x0, 0.0) and almost(x1, 0.0):
        return 1

    if almost(x0, W) and almost(x1, W):
        return 1

    # resto: convección
    return 2


def build_boundary_mesh() -> tuple[
    list[tuple[float, float]],
    list[tuple[int, int, int]],
]:
    """
    Mantiene la lógica original, pero:
    - duplica nodos entre segmentos
    - conserva generación automática
    - IDs de frontera corregidos
    """

    vertices = [
        (0.000, 0.000),
        (0.025, 0.000),
        (0.025, 0.025),
        (0.020, 0.025),
        (0.020, 0.005),
        (0.015, 0.005),
        (0.015, 0.025),
        (0.010, 0.025),
        (0.010, 0.005),
        (0.005, 0.005),
        (0.005, 0.025),
        (0.000, 0.025),
        (0.000, 0.000),
    ]

    nodes: list[tuple[float, float]] = []
    elems: list[tuple[int, int, int]] = []

    for a, b in zip(vertices[:-1], vertices[1:]):

        pts = segment_points(a, b, LE)

        start = len(nodes) + 1

        # DUPLICA nodos en cambios de tramo
        for p in pts:
            nodes.append(p)

        bid = boundary_id(a, b)

        for i in range(len(pts) - 1):
            elems.append((start + i, start + i + 1, bid))

    return nodes, elems


def read_mef_internal_points() -> list[tuple[float, float]]:
    """
    Usa alturas del MEF y elimina contorno.
    """
    if not MEF_FILE.exists():
        return [(X_MID, i * 0.0005) for i in range(1, 50)]

    ys: list[float] = []

    with MEF_FILE.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("%"):
                continue

            parts = line.split()

            if len(parts) >= 2:
                ys.append(float(parts[0]))

    ys = [y for y in ys if 1e-12 < y < H - 1e-12]

    return [(X_MID, y) for y in ys]


def write_model(
    filename: Path,
    nodes: list[tuple[float, float]],
    elems: list[tuple[int, int, int]],
    internal: list[tuple[float, float]],
) -> None:
    with filename.open("w", encoding="utf-8") as f:

        f.write(f"{TA:.12g}\n")
        f.write(f"{K:.12g}\n")

        f.write(f"{len(BOUNDARIES)}\n")
        for tcc, vcc in BOUNDARIES:
            f.write(f"{tcc:d} {vcc:.12g}\n")

        f.write(f"{len(nodes)}\n")
        for x, y in nodes:
            f.write(f"{x:.12g} {y:.12g}\n")

        f.write(f"{len(elems)}\n")
        for n1, n2, bid in elems:
            f.write(f"{n1:d} {n2:d} {bid:d}\n")

        f.write(f"{len(internal)}\n")
        for x, y in internal:
            f.write(f"{x:.12g} {y:.12g}\n")


def plot_preview(
    nodes: list[tuple[float, float]],
    elems: list[tuple[int, int, int]],
    internal: list[tuple[float, float]],
) -> None:
    plt.figure(figsize=(6, 6))

    colors = {1: "blue", 2: "cyan", 3: "yellow"}

    for n1, n2, bid in elems:
        x0, y0 = nodes[n1 - 1]
        x1, y1 = nodes[n2 - 1]

        plt.plot(
            [x0, x1],
            [y0, y1],
            color=colors[bid],
            linewidth=2,
        )

    plt.scatter(
        [p[0] for p in internal],
        [p[1] for p in internal],
        c="red",
        s=8,
    )

    plt.axis("equal")
    plt.grid(True)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.tight_layout()
    plt.savefig(OUTDIR / "malla_disipador_preview.png", dpi=200)


def main() -> None:
    nodes, elems = build_boundary_mesh()
    internal = read_mef_internal_points()

    write_model(OUTDIR / "malla_disipador.txt", nodes, elems, internal)
    plot_preview(nodes, elems, internal)

    print(f"Nodos de contorno: {len(nodes)}")
    print(f"Elementos de contorno: {len(elems)}")
    print(f"Puntos internos: {len(internal)}")
    print("Generado: malla_disipador.txt")
    print("Generado: malla_disipador_preview.png")


if __name__ == "__main__":
    main()