
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Optional
import numpy as np


class MeshAlignmentError(ValueError):
    pass


class ThermalMode(Enum):
    AUTO = "auto"
    CONDUCTION = "conduction"
    CONVECTION = "convection"
    FLUID_TRANSPORT = "fluid_transport"
    RADIATION = "radiation"


def normalize_thermal_mode(mode: Optional[str | ThermalMode]) -> Optional[ThermalMode]:
    if mode is None:
        return None

    if isinstance(mode, ThermalMode):
        return None if mode == ThermalMode.AUTO else mode

    text = str(mode).strip().lower()

    if text in ("", "none", "auto", "default", "padrao", "padrão"):
        return None

    aliases = {
        "conduction": ThermalMode.CONDUCTION,
        "conducao": ThermalMode.CONDUCTION,
        "condução": ThermalMode.CONDUCTION,
        "cond": ThermalMode.CONDUCTION,
        "convection": ThermalMode.CONVECTION,
        "conveccao": ThermalMode.CONVECTION,
        "convecção": ThermalMode.CONVECTION,
        "conv": ThermalMode.CONVECTION,
        "fluid_transport": ThermalMode.FLUID_TRANSPORT,
        "transport": ThermalMode.FLUID_TRANSPORT,
        "transporte": ThermalMode.FLUID_TRANSPORT,
        "enthalpy": ThermalMode.FLUID_TRANSPORT,
        "entalpico": ThermalMode.FLUID_TRANSPORT,
        "entálpico": ThermalMode.FLUID_TRANSPORT,
        "radiation": ThermalMode.RADIATION,
        "radiacao": ThermalMode.RADIATION,
        "radiação": ThermalMode.RADIATION,
        "rad": ThermalMode.RADIATION,
    }

    if text not in aliases:
        raise ValueError(
            f"thermal_mode inválido: {mode!r}. "
            "Use None, 'conduction', 'convection', 'fluid_transport' ou 'radiation'."
        )

    return aliases[text]


def material_name(material: Any) -> str:
    return getattr(material, "name", str(material))


@dataclass
class MaterialRegion:
    name: str
    material: Any
    x0: float
    x1: float
    y0: float
    y1: float
    tag: Optional[str] = None
    thermal_mode: Optional[ThermalMode | str] = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.thermal_mode = normalize_thermal_mode(self.thermal_mode)

        if self.x1 <= self.x0:
            raise ValueError(f"Região '{self.name}': x1 deve ser maior que x0.")
        if self.y1 <= self.y0:
            raise ValueError(f"Região '{self.name}': y1 deve ser maior que y0.")

    @property
    def width(self) -> float:
        return self.x1 - self.x0

    @property
    def height(self) -> float:
        return self.y1 - self.y0

    @property
    def area(self) -> float:
        return self.width * self.height

    def contains_center(self, x: float, y: float) -> bool:
        return self.x0 <= x <= self.x1 and self.y0 <= y <= self.y1

    def overlaps(self, other: "MaterialRegion", tol: float = 1e-14) -> bool:
        ox = (self.x0 < other.x1 - tol) and (self.x1 > other.x0 + tol)
        oy = (self.y0 < other.y1 - tol) and (self.y1 > other.y0 + tol)
        return ox and oy

    def shifted_to_grid(self, x_edges: np.ndarray, y_edges: np.ndarray) -> "MaterialRegion":
        def nearest(value: float, edges: np.ndarray) -> float:
            return float(edges[np.argmin(np.abs(edges - value))])

        return MaterialRegion(
            name=self.name,
            material=self.material,
            x0=nearest(self.x0, x_edges),
            x1=nearest(self.x1, x_edges),
            y0=nearest(self.y0, y_edges),
            y1=nearest(self.y1, y_edges),
            tag=self.tag,
            thermal_mode=self.thermal_mode,
            metadata=dict(self.metadata),
        )


@dataclass
class BoundaryCondition:
    side: str
    kind: str
    data: dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None

    def __post_init__(self) -> None:
        allowed = {"left", "right", "bottom", "top"}
        if self.side not in allowed:
            raise ValueError(f"side deve ser um de {sorted(allowed)}.")
        self.kind = str(self.kind).strip().lower()


HeatFunction = Callable[[float, float], float]


@dataclass
class HeatSource:
    region: str
    kind: str
    value: Optional[float] = None
    function: Optional[HeatFunction] = None
    data: dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None

    def __post_init__(self) -> None:
        self.kind = str(self.kind).strip().lower()

        allowed = {"volumetric", "surface"}
        if self.kind not in allowed:
            raise ValueError(f"kind de fonte deve ser um de {sorted(allowed)}.")

        if self.value is None and self.function is None:
            raise ValueError("HeatSource precisa de value ou function.")

    def evaluate(self, x: float, y: float) -> float:
        if self.function is not None:
            return float(self.function(x, y))
        return float(self.value)


@dataclass
class InterfaceRule:
    region_a: str
    region_b: str
    mode: ThermalMode | str
    data: dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None

    def __post_init__(self) -> None:
        mode = normalize_thermal_mode(self.mode)
        if mode is None:
            raise ValueError("InterfaceRule precisa de modo explícito.")
        self.mode = mode

    def matches(self, a: str, b: str) -> bool:
        return {self.region_a, self.region_b} == {a, b}


@dataclass
class Geometry2D:
    width: float
    height: float

    regions: list[MaterialRegion] = field(default_factory=list)
    bcs: list[dict[str, Any]] = field(default_factory=list)
    boundary_conditions: list[BoundaryCondition] = field(default_factory=list)
    heat_sources: list[HeatSource] = field(default_factory=list)
    interface_rules: list[InterfaceRule] = field(default_factory=list)

    nx: Optional[int] = None
    ny: Optional[int] = None
    dx: Optional[float] = None
    dy: Optional[float] = None
    x_edges: Optional[np.ndarray] = None
    y_edges: Optional[np.ndarray] = None
    x_centers: Optional[np.ndarray] = None
    y_centers: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.width <= 0:
            raise ValueError("width deve ser positivo.")
        if self.height <= 0:
            raise ValueError("height deve ser positivo.")

    def material(
        self,
        name: str,
        material: Any,
        x0: float,
        x1: float,
        y0: float,
        y1: float,
        tag: Optional[str] = None,
        thermal_mode: Optional[str | ThermalMode] = None,
        **metadata: Any,
    ) -> MaterialRegion:
        region = MaterialRegion(
            name=name,
            material=material,
            x0=float(x0),
            x1=float(x1),
            y0=float(y0),
            y1=float(y1),
            tag=tag,
            thermal_mode=thermal_mode,
            metadata=dict(metadata),
        )

        self._check_inside_domain(region)

        for old in self.regions:
            if region.overlaps(old):
                raise ValueError(
                    f"Região '{region.name}' sobrepõe a região '{old.name}'."
                )

        self.regions.append(region)
        return region

    def add_region(
        self,
        name: str,
        material: Any,
        x_min: float,
        x_max: float,
        y_min: float,
        y_max: float,
        tag: Optional[str] = None,
        validate_alignment: bool = False,
        thermal_mode: Optional[str | ThermalMode] = None,
        **metadata: Any,
    ) -> MaterialRegion:
        region = self.material(
            name=name,
            material=material,
            x0=x_min,
            x1=x_max,
            y0=y_min,
            y1=y_max,
            tag=tag,
            thermal_mode=thermal_mode,
            **metadata,
        )

        if validate_alignment and self.nx is not None:
            problems = self.alignment_problems()
            if problems:
                raise MeshAlignmentError(self._format_alignment_error(problems))

        return region

    def bc(self, side: str, kind: str, **data: Any) -> BoundaryCondition:
        bc = BoundaryCondition(side=side, kind=kind, data=dict(data))
        self.boundary_conditions.append(bc)
        self.bcs.append({"side": side, "kind": kind, **data})
        return bc

    def add_boundary_condition(
        self,
        side: str,
        bc_type: str,
        value: Optional[float] = None,
        data: Optional[dict[str, Any]] = None,
        name: Optional[str] = None,
    ) -> BoundaryCondition:
        data_dict = dict(data or {})
        kind = str(bc_type).lower()

        if value is not None:
            if kind in ("temperature", "dirichlet"):
                data_dict.setdefault("T", value)
            elif kind in ("heat_flux", "flux", "neumann"):
                data_dict.setdefault("q", value)
            elif kind in ("convection", "convective"):
                data_dict.setdefault("h", value)
            else:
                data_dict.setdefault("value", value)

        bc = BoundaryCondition(side=side, kind=bc_type, data=data_dict, name=name)
        self.boundary_conditions.append(bc)
        self.bcs.append({"side": side, "kind": bc_type, **data_dict})
        return bc

    def source(
        self,
        region: str,
        kind: str,
        value: Optional[float] = None,
        function: Optional[HeatFunction] = None,
        name: Optional[str] = None,
        **data: Any,
    ) -> HeatSource:
        self.get_region(region)

        src = HeatSource(
            region=region,
            kind=kind,
            value=value,
            function=function,
            data=dict(data),
            name=name,
        )

        self.heat_sources.append(src)
        return src

    def interface(
        self,
        region_a: str,
        region_b: str,
        mode: str | ThermalMode,
        name: Optional[str] = None,
        **data: Any,
    ) -> InterfaceRule:
        self.get_region(region_a)
        self.get_region(region_b)

        rule = InterfaceRule(
            region_a=region_a,
            region_b=region_b,
            mode=mode,
            data=dict(data),
            name=name,
        )

        self.interface_rules.append(rule)
        return rule

    def get_region(self, name: str) -> MaterialRegion:
        for region in self.regions:
            if region.name == name:
                return region
        raise KeyError(f"Região '{name}' não encontrada.")

    def source_for_region(self, region_name: str) -> list[HeatSource]:
        return [source for source in self.heat_sources if source.region == region_name]

    def interface_rule_between(
        self,
        region_a: str,
        region_b: str,
    ) -> Optional[InterfaceRule]:
        for rule in self.interface_rules:
            if rule.matches(region_a, region_b):
                return rule
        return None

    def transfer_mode_between_regions(
        self,
        region_a: str | MaterialRegion,
        region_b: str | MaterialRegion,
    ) -> str:
        A = self.get_region(region_a) if isinstance(region_a, str) else region_a
        B = self.get_region(region_b) if isinstance(region_b, str) else region_b

        rule = self.interface_rule_between(A.name, B.name)
        if rule is not None:
            return rule.mode.value

        if A.thermal_mode is not None and B.thermal_mode is not None:
            if A.thermal_mode != B.thermal_mode:
                raise ValueError(
                    f"Conflito de thermal_mode entre regiões '{A.name}' "
                    f"({A.thermal_mode.value}) e '{B.name}' ({B.thermal_mode.value}). "
                    "Use geom.interface(...) para resolver explicitamente."
                )
            return A.thermal_mode.value

        if A.thermal_mode is not None:
            return A.thermal_mode.value

        if B.thermal_mode is not None:
            return B.thermal_mode.value

        if hasattr(A.material, "default_interaction_with"):
            return A.material.default_interaction_with(B.material)

        if material_name(A.material) == material_name(B.material):
            return "conduction"

        return "conduction"

    def mesh(
        self,
        nx: int,
        ny: int,
        *,
        fix: str = "raise",
        max_factor: int = 20,
    ) -> None:
        if nx <= 0 or ny <= 0:
            raise ValueError("nx e ny devem ser positivos.")

        fix = fix.lower().strip()
        allowed = {"raise", "auto", "snap", "ask"}
        if fix not in allowed:
            raise ValueError(f"fix deve ser um de {sorted(allowed)}.")

        self._set_mesh(nx, ny)
        problems = self.alignment_problems()

        if not problems:
            return

        if fix == "raise":
            raise MeshAlignmentError(self._format_alignment_error(problems))

        if fix == "ask":
            choice = self._ask_fix_choice(problems)

            if choice == "auto":
                fix = "auto"
            elif choice == "snap":
                fix = "snap"
            else:
                raise MeshAlignmentError(self._format_alignment_error(problems))

        if fix == "snap":
            self.snap_to_grid()
            return

        if fix == "auto":
            new_nx, new_ny = self.suggest_aligned_mesh(
                target_nx=nx,
                target_ny=ny,
                max_factor=max_factor,
            )
            self._set_mesh(new_nx, new_ny)

            problems = self.alignment_problems()
            if problems:
                raise MeshAlignmentError(
                    "Não foi possível ajustar automaticamente a malha.\n"
                    + self._format_alignment_error(problems)
                )

            return

    def mesh_by_points(
        self,
        npx: int,
        npy: int,
        *,
        fix: str = "raise",
        max_factor: int = 20,
    ) -> None:
        if npx < 2 or npy < 2:
            raise ValueError("npx e npy devem ser pelo menos 2.")

        self.mesh(nx=npx - 1, ny=npy - 1, fix=fix, max_factor=max_factor)

    def _set_mesh(self, nx: int, ny: int) -> None:
        self.nx = int(nx)
        self.ny = int(ny)

        self.dx = self.width / self.nx
        self.dy = self.height / self.ny

        self.x_edges = np.linspace(0.0, self.width, self.nx + 1)
        self.y_edges = np.linspace(0.0, self.height, self.ny + 1)

        self.x_centers = 0.5 * (self.x_edges[:-1] + self.x_edges[1:])
        self.y_centers = 0.5 * (self.y_edges[:-1] + self.y_edges[1:])

    def alignment_problems(self, tol: float = 1e-12) -> list[str]:
        self._require_mesh()

        problems = []

        for region in self.regions:
            checks = [
                ("x0", region.x0, self.x_edges),
                ("x1", region.x1, self.x_edges),
                ("y0", region.y0, self.y_edges),
                ("y1", region.y1, self.y_edges),
            ]

            for label, value, edges in checks:
                dist = float(np.min(np.abs(edges - value)))
                if dist > tol:
                    nearest = float(edges[np.argmin(np.abs(edges - value))])
                    problems.append(
                        f"{region.name}.{label}={value:g} não está na malha; "
                        f"linha mais próxima={nearest:g}; erro={dist:g}"
                    )

        return problems

    def check(self) -> list[str]:
        self._require_mesh()

        problems = []
        problems.extend(self.alignment_problems())

        mat_grid = self.material_grid()
        void_count = int(np.sum(mat_grid == "VOID"))

        if void_count > 0:
            problems.append(f"{void_count} célula(s) sem material definido.")

        for source in self.heat_sources:
            try:
                self.get_region(source.region)
            except KeyError:
                problems.append(
                    f"Fonte '{source.name}' aponta para região inexistente "
                    f"'{source.region}'."
                )

        return problems

    def snap_to_grid(self) -> None:
        self._require_mesh()

        new_regions = []
        for region in self.regions:
            new_regions.append(region.shifted_to_grid(self.x_edges, self.y_edges))

        for i, r in enumerate(new_regions):
            if r.x1 <= r.x0 or r.y1 <= r.y0:
                raise MeshAlignmentError(f"Snap tornou a região '{r.name}' inválida.")

            for old in new_regions[:i]:
                if r.overlaps(old):
                    raise MeshAlignmentError(
                        f"Snap gerou sobreposição entre '{r.name}' e '{old.name}'."
                    )

        self.regions = new_regions

    def suggest_aligned_mesh(
        self,
        target_nx: int,
        target_ny: int,
        max_factor: int = 20,
        tol: float = 1e-10,
    ) -> tuple[int, int]:
        x_values, y_values = self._all_boundaries()

        nx_max = max(target_nx, 1) * max_factor
        ny_max = max(target_ny, 1) * max_factor

        possible_nx = []
        possible_ny = []

        for nx in range(target_nx, nx_max + 1):
            dx = self.width / nx
            if all(self._is_multiple(v, dx, tol=tol) for v in x_values):
                possible_nx.append(nx)

        for ny in range(target_ny, ny_max + 1):
            dy = self.height / ny
            if all(self._is_multiple(v, dy, tol=tol) for v in y_values):
                possible_ny.append(ny)

        if not possible_nx or not possible_ny:
            raise MeshAlignmentError(
                "Não foi encontrada malha alinhada no intervalo de busca. "
                "Aumente max_factor ou use fix='snap'."
            )

        return possible_nx[0], possible_ny[0]

    def _all_boundaries(self) -> tuple[list[float], list[float]]:
        x_values = [0.0, self.width]
        y_values = [0.0, self.height]

        for r in self.regions:
            x_values.extend([r.x0, r.x1])
            y_values.extend([r.y0, r.y1])

        return x_values, y_values

    @staticmethod
    def _is_multiple(value: float, step: float, tol: float = 1e-10) -> bool:
        q = value / step
        return abs(q - round(q)) <= tol

    def _format_alignment_error(self, problems: list[str]) -> str:
        preview = "\n".join(f"- {p}" for p in problems[:20])

        if len(problems) > 20:
            preview += f"\n... mais {len(problems) - 20} problema(s)."

        return (
            "A malha corta fronteiras de materiais.\n"
            "Uma célula não pode conter dois materiais simultaneamente.\n\n"
            f"{preview}\n\n"
            "Escolha uma das opções:\n"
            "1) usar fix='auto' para ajustar nx/ny automaticamente;\n"
            "2) usar fix='snap' para mover as fronteiras para a malha;\n"
            "3) alterar manualmente nx, ny ou as coordenadas das regiões."
        )

    def _ask_fix_choice(self, problems: list[str]) -> str:
        print(self._format_alignment_error(problems))
        print()
        print("Digite:")
        print("  1 para fix='auto'")
        print("  2 para fix='snap'")
        print("  qualquer outra coisa para cancelar")

        answer = input("Escolha: ").strip()

        if answer == "1":
            return "auto"
        if answer == "2":
            return "snap"
        return "cancel"

    def material_grid(self) -> np.ndarray:
        self._require_mesh()

        grid = np.full((self.ny, self.nx), "VOID", dtype=object)

        for j, y in enumerate(self.y_centers):
            for i, x in enumerate(self.x_centers):
                for region in self.regions:
                    if region.contains_center(x, y):
                        grid[j, i] = region.material
                        break

        return grid

    def region_grid(self) -> np.ndarray:
        self._require_mesh()

        grid = np.full((self.ny, self.nx), "VOID", dtype=object)

        for j, y in enumerate(self.y_centers):
            for i, x in enumerate(self.x_centers):
                for region in self.regions:
                    if region.contains_center(x, y):
                        grid[j, i] = region.name
                        break

        return grid

    def region_object_grid(self) -> np.ndarray:
        self._require_mesh()

        grid = np.full((self.ny, self.nx), None, dtype=object)

        for j, y in enumerate(self.y_centers):
            for i, x in enumerate(self.x_centers):
                for region in self.regions:
                    if region.contains_center(x, y):
                        grid[j, i] = region
                        break

        return grid

    def cell_bounds(self, i: int, j: int) -> tuple[float, float, float, float]:
        self._require_mesh()

        if i < 0 or i >= self.nx:
            raise IndexError("i fora da malha.")
        if j < 0 or j >= self.ny:
            raise IndexError("j fora da malha.")

        return (
            float(self.x_edges[i]),
            float(self.x_edges[i + 1]),
            float(self.y_edges[j]),
            float(self.y_edges[j + 1]),
        )

    def cell_center(self, i: int, j: int) -> tuple[float, float]:
        self._require_mesh()
        return float(self.x_centers[i]), float(self.y_centers[j])

    def show(
        self,
        grid: bool = True,
        labels: bool = False,
        ax: Optional[Any] = None,
        *,
        mode: str = "schematic",
        legend: bool = True,
        title: str = "Geometria e discretização",
    ) -> Any:
        self._require_mesh()

        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        from matplotlib.ticker import ScalarFormatter

        if ax is None:
            figsize = (7.5, 7.5) if mode == "physical" else (12.0, 6.5)
            _, ax = plt.subplots(figsize=figsize)

        colors = plt.rcParams["axes.prop_cycle"].by_key().get(
            "color", ["C0", "C1", "C2", "C3", "C4", "C5"]
        )

        mat_order: list[str] = []
        mat_to_color: dict[str, str] = {}

        for region in self.regions:
            key = material_name(region.material)
            if key not in mat_to_color:
                mat_to_color[key] = colors[len(mat_order) % len(colors)]
                mat_order.append(key)

        mat_to_regions: dict[str, list[str]] = {key: [] for key in mat_order}
        handles = []

        domain = patches.Rectangle(
            (0, 0),
            self.width,
            self.height,
            fill=False,
            linewidth=2,
            edgecolor="black",
        )
        ax.add_patch(domain)

        for region in self.regions:
            key = material_name(region.material)
            color = mat_to_color[key]
            mat_to_regions[key].append(region.name)

            rect = patches.Rectangle(
                (region.x0, region.y0),
                region.width,
                region.height,
                facecolor=color,
                alpha=0.35,
                edgecolor="black",
                linewidth=1.2,
            )
            ax.add_patch(rect)

            if labels:
                label = region.name
                if region.thermal_mode is not None:
                    label += f"\n{region.thermal_mode.value}"

                ax.text(
                    region.x0 + region.width / 2,
                    region.y0 + region.height / 2,
                    label,
                    ha="center",
                    va="center",
                    fontsize=8,
                    bbox=dict(
                        boxstyle="round,pad=0.2",
                        facecolor="white",
                        alpha=0.75,
                        edgecolor="none",
                    ),
                )

        for mat in mat_order:
            regions = ", ".join(mat_to_regions[mat])
            handles.append(
                patches.Patch(
                    facecolor=mat_to_color[mat],
                    edgecolor="black",
                    alpha=0.35,
                    label=f"{mat} ({regions})",
                )
            )

        if grid:
            for x in self.x_edges:
                ax.axvline(x, linewidth=0.2, alpha=0.35)
            for y in self.y_edges:
                ax.axhline(y, linewidth=0.2, alpha=0.35)

        if mode == "physical":
            ax.set_aspect("equal", adjustable="box")
        else:
            ax.set_aspect("auto")

        ax.set_xlim(0, self.width)
        ax.set_ylim(0, self.height)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(title)
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))

        if legend and handles:
            ax.legend(
                handles=handles,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                fontsize=8,
                frameon=True,
            )

        return ax

    def preview(
        self,
        *,
        grid: bool = False,
        mode: str = "schematic",
        title: str = "Preview da geometria",
    ) -> None:
        import matplotlib.pyplot as plt

        plt.close("all")
        self.show(grid=grid, labels=False, mode=mode, legend=True, title=title)
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.1)

    def summary(self) -> str:
        mesh_text = "malha ainda não definida"

        if self.nx is not None:
            mesh_text = (
                f"nx={self.nx}, ny={self.ny}, "
                f"dx={self.dx:g}, dy={self.dy:g}, "
                f"pontos=({self.nx + 1}, {self.ny + 1})"
            )

        lines = [
            "=== Geometry2D ===",
            f"Domínio: width={self.width:g} m, height={self.height:g} m",
            "Origem: (0,0) no canto inferior esquerdo",
            f"Malha: {mesh_text}",
            "",
            "Regiões:",
        ]

        if not self.regions:
            lines.append("  nenhuma")
        else:
            for region in self.regions:
                mat = material_name(region.material)
                mode = "padrão" if region.thermal_mode is None else region.thermal_mode.value
                lines.append(
                    f"  - {region.name}: material={mat}, "
                    f"x=[{region.x0:g}, {region.x1:g}], "
                    f"y=[{region.y0:g}, {region.y1:g}], "
                    f"thermal_mode={mode}"
                )

        lines.append("")
        lines.append("Condições de contorno:")

        if not self.boundary_conditions:
            lines.append("  nenhuma")
        else:
            for bc in self.boundary_conditions:
                lines.append(f"  - {bc.side}: {bc.kind}, data={bc.data}")

        lines.append("")
        lines.append("Fontes de calor:")

        if not self.heat_sources:
            lines.append("  nenhuma")
        else:
            for source in self.heat_sources:
                src_name = source.name or "-"
                value_text = "function" if source.function is not None else f"value={source.value:g}"
                lines.append(
                    f"  - {src_name}: region={source.region}, "
                    f"kind={source.kind}, {value_text}"
                )

        lines.append("")
        lines.append("Regras explícitas de interface:")

        if not self.interface_rules:
            lines.append("  nenhuma")
        else:
            for rule in self.interface_rules:
                rule_name = rule.name or "-"
                lines.append(
                    f"  - {rule_name}: {rule.region_a} <-> {rule.region_b}, "
                    f"mode={rule.mode.value}, data={rule.data}"
                )

        return "\n".join(lines)

    def _check_inside_domain(self, region: MaterialRegion) -> None:
        tol = 1e-14

        if region.x0 < -tol or region.x1 > self.width + tol:
            raise ValueError(f"Região '{region.name}' fora do domínio em x.")

        if region.y0 < -tol or region.y1 > self.height + tol:
            raise ValueError(f"Região '{region.name}' fora do domínio em y.")

    def _require_mesh(self) -> None:
        if self.nx is None:
            raise RuntimeError("Defina a malha antes: geom.mesh(nx, ny).")
