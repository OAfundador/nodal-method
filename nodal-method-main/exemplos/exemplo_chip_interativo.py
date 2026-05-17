"""
exemplo_chip_interativo.py  (+ comandos reator_*)

Modo script:   python exemplos/exemplo_chip_interativo.py --script exemplos/comandos_reator.txt
Modo interativo: python exemplos/exemplo_chip_interativo.py --interativo
Modo demo:     python exemplos/exemplo_chip_interativo.py
"""

from __future__ import annotations
import shlex, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from condutancias import conduction_G, convection_G
from geometria import Geometry2D
from materiais import Material, MaterialPhase
from nos import NodalNetwork, NodeKind, TransferKind, build_network_from_geometry
from solver import solve_steady_state

try:
    from reator_placa.geometria_reator import GeometriaReator
    from reator_placa.modelo_nodal_reator import (
        ConfigCaso, MODO_COS, MODO_CONST,
        construir_rede, extrair_resultados, gerar_chute_inicial,
    )
    from reator_placa.exemplo_reator_placa import (
        texto_resumo, salvar_csv, plot_temperaturas, plot_fluxo,
    )
    _REATOR_DISPONIVEL = True
except Exception:
    _REATOR_DISPONIVEL = False

NODE_KIND_MAP = {
    "diffusion": NodeKind.DIFFUSION,
    "arithmetic": NodeKind.ARITHMETIC,
    "fluid": NodeKind.FLUID,
    "boundary": NodeKind.BOUNDARY,
}
LINK_KIND_MAP = {
    "cond": TransferKind.CONDUCTION,
    "conv": TransferKind.CONVECTION,
    "fluid": TransferKind.FLUID_TRANSPORT,
    "rad": TransferKind.RADIATION,
}

HELP = """
============================================================
  METODO NODAL - REPL INTERATIVO
  Comandos disponiveis
============================================================

--- MATERIAIS ---
  material <nome> [phase=solid|fluid] [k=v] [rho=v] [cp=v] [mu=v]
       Define um material reutilizavel. Exemplos:
         material cobre k=380 rho=8900 cp=385
         material agua  phase=fluid k=0.6 rho=997 cp=4182 mu=8.9e-4
  materiais        Lista os materiais disponiveis (builtin + sessao)

--- GEOMETRIA 2D ---
  domain W=v H=v
       Cria dominio fisico largura W x altura H [m]
  region <nome> [mat=<material> | k=v] x0=v x1=v y0=v y1=v
       Define regiao retangular com material
  source <regiao> q=v
       Fonte volumetrica q_triplo [W/m3] na regiao
  bc <face> <tipo> [T=v] [h=v] [T_inf=v]
       face: bottom | top | left | right
       tipo: temperature | convection | adiabatic
  mesh nx=v ny=v [fix=auto|none]
  show_geom        Resume dominio, regioes, BCs e malha
  build_from_geom  [tz=v]   Gera rede nodal automaticamente

--- REDE NODAL DIRETA ---
  node <nome> <tipo> [Q=v] [T_fixed=v] [T_init=v] [V=v]
       tipos: diffusion | fluid | arithmetic | boundary
  link <a> <b> <tipo> G=v [dir=fwd|bwd|undirected]
       tipos: cond | conv | fluid | rad | equiv
  fluid_chain n1 n2 [n3 ...] mdot=v cp=v [T_in=v]
       Cadeia FLUID_TRANSPORT direcional n1->n2->...
       T_in cria no BOUNDARY de entrada automaticamente.
  show             Resumo da rede (nos e ligacoes)
  solve            Resolve estado permanente (Newton-Raphson)
  reset            Limpa rede, geometria e materiais

--- CONDUTANCIAS ---
  g_cond k=v A=v L=v   -> G = k*A/L  [W/K]
  g_conv h=v A=v        -> G = h*A    [W/K]

--- VISUALIZACAO ---
  viz  [png=arquivo.png] [titulo=Texto]
       PNG da rede: nos coloridos por T, links com setas
  gif  [gif=arquivo.gif] [fps=10] [steps=40] [dpi=100] [titulo=Texto]
       GIF pseudo-transiente: frio -> regime permanente

--- REATOR TIPO PLACA (atalhos) ---
  reator_config [n_axial=40] [vazao=3e-4] [potencia=5000]
                [modo=cos|constante] [T_in=30] [P_in=150000] [dP=10000]
  reator_geom   [Lx=0.5] [Ly=0.06] [df=0.001] [dcl=0.0005]
                [dch=0.003] [Lcanal=0.065]
  reator_solve / reator_show / reator_plot [dir=.]
  reator_csv [dir=.] / reator_gif [dir=.] [fps=12] [dpi=100]

--- CONTROLE ---
  help / quit / exit
============================================================
"""

class EstadoREPL:
    def __init__(self):
        self.net = NodalNetwork()
        self.ids = {}
        self.x_counter = 0.0
        self.geom = None
        self.materials = {}
        self.reator_cfg = None
        self.reator_geom = None
        self.reator_result = None

    def reset(self):
        self.__init__()

def parse_kwargs(tokens, lower_keys=True):
    out = {}
    for t in tokens:
        if "=" not in t:
            raise ValueError(f"opcao invalida: {t!r}")
        k, v = t.split("=", 1)
        out[k.strip().lower() if lower_keys else k.strip()] = v.strip()
    return out

def fnum(s): return float(str(s).replace(",", "."))

def cmd_node(state, args):
    nome, tipo = args[0], args[1].lower()
    if tipo not in NODE_KIND_MAP: raise ValueError(f"tipo invalido: {tipo!r}")
    if nome in state.ids: raise ValueError(f"no {nome!r} ja existe.")
    kw = parse_kwargs(args[2:])
    kind = NODE_KIND_MAP[tipo]
    if kind is NodeKind.BOUNDARY:
        Tf = fnum(kw["t_fixed"])
        nid = state.net.add_node(name=nome, kind=kind, x=state.x_counter, y=0.0, fixed_temperature=Tf, temperature=Tf)
    else:
        nid = state.net.add_node(name=nome, kind=kind, x=state.x_counter, y=0.0,
                                  volume=fnum(kw.get("v","1")), source=fnum(kw.get("q","0")), temperature=fnum(kw.get("t_init","30")))
    state.ids[nome] = nid; state.x_counter += 1.0
    return f"  -> no {nome} criado (id={nid}, kind={kind.value})"

def cmd_link(state, args):
    from nos import LinkDirection
    a, b, tipo = args[0], args[1], args[2].lower()
    kw = parse_kwargs(args[3:])
    dir_map = {
        "fwd": LinkDirection.I_TO_J, "i_to_j": LinkDirection.I_TO_J,
        "bwd": LinkDirection.J_TO_I, "j_to_i": LinkDirection.J_TO_I,
        "undirected": LinkDirection.UNDIRECTED,
    }
    dir_kw = kw.get("dir", "undirected").lower()
    direction = dir_map.get(dir_kw, LinkDirection.UNDIRECTED)
    state.net.add_link(state.ids[a], state.ids[b], LINK_KIND_MAP[tipo],
                       conductance=fnum(kw["g"]), direction=direction,
                       name=a + "_" + b)
    G_val = fnum(kw["g"])
    return "  -> ligacao " + a + "-" + b + " (" + tipo + ") G=" + f"{G_val:.6g}" + " W/K  dir=" + dir_kw

def cmd_show(state): return state.net.summary(max_nodes=20, max_links=30)

def cmd_solve(state):
    if len(state.net.nodes) == 0: return "  rede vazia."
    r = solve_steady_state(state.net, tol=1e-9, max_iter=500)
    out = [f"  convergiu={r.success}, |R|={r.residual_norm:.2e}, iter={r.iterations}"]
    if state.ids:
        out += ["", "  Temperaturas:"] + [f"    {n:8s} = {state.net.nodes[i].temperature:9.4f} C" for n,i in state.ids.items()]
    return "\n".join(out)

def _req(state):
    if not state.geom: raise ValueError("Use domain W= H= primeiro.")
    return state.geom


# ---------------------------------------------------------------------------
# Materiais
# ---------------------------------------------------------------------------
_MATERIAIS_BUILTIN = {
    "cobre":    dict(phase="solid",  k=380.0,  rho=8900.0, cp=385.0),
    "aluminio": dict(phase="solid",  k=205.0,  rho=2700.0, cp=900.0),
    "aco":      dict(phase="solid",  k=50.0,   rho=7800.0, cp=500.0),
    "silicio":  dict(phase="solid",  k=150.0,  rho=2330.0, cp=712.0),
    "uo2":      dict(phase="solid",  k=3.6,    rho=10970.0,cp=247.0),
    "agua":     dict(phase="fluid",  k=0.600,  rho=997.0,  cp=4182.0, mu=8.9e-4),
    "ar":       dict(phase="fluid",  k=0.026,  rho=1.18,   cp=1007.0, mu=1.85e-5),
}

def cmd_material(state, args):
    if not args:
        raise ValueError("material precisa de nome. Ex: material cobre k=380")
    nome = args[0].lower()
    kw   = parse_kwargs(args[1:]) if len(args) > 1 else {}
    props = dict(_MATERIAIS_BUILTIN.get(nome, {}))
    for campo in ("k","rho","cp","mu"):
        if campo in kw:
            props[campo] = fnum(kw[campo])
    if "phase" in kw:
        props["phase"] = kw["phase"].lower()
    if not props.get("phase"):
        props["phase"] = "solid"
    state.materials[nome] = props
    partes = []
    for campo, unid in [("k","W/(mK)"),("rho","kg/m3"),("cp","J/(kgK)"),("mu","Pa.s")]:
        if campo in props:
            partes.append(campo + "=" + str(props[campo]) + " " + unid)
    return (
        "  -> material " + repr(nome)
        + " [" + props["phase"] + "]  "
        + "  ".join(partes)
    )

def cmd_materiais(state, args):
    linhas = ['  --- Builtin ---']
    for n, d in _MATERIAIS_BUILTIN.items():
        k_s = ("k=" + str(d["k"])) if "k" in d else ""
        linhas.append("    " + n.ljust(10) + " [" + d["phase"] + "]  " + k_s)
    if state.materials:
        linhas.append('  --- Definidos na sessao ---')
        for n, d in state.materials.items():
            k_s = "k=" + str(d.get("k","-"))
            linhas.append("    " + n.ljust(10) + " [" + d.get("phase","solid") + "]  " + k_s)
    return "\n".join(linhas)

def cmd_domain(state, args):
    kw=parse_kwargs(args); state.geom=Geometry2D(width=fnum(kw["w"]),height=fnum(kw["h"]))
    return f"  -> Geometry2D {kw['w']} x {kw['h']} m"
def cmd_region(state, args):
    geom = _req(state); kw = parse_kwargs(args[1:]); nome = args[0]
    if "mat" in kw:
        mat_nome = kw["mat"].lower()
        props = dict(_MATERIAIS_BUILTIN.get(mat_nome, {}))
        props.update(state.materials.get(mat_nome, {}))
        if not props:
            raise ValueError("material " + repr(mat_nome) + " nao definido.")
        ph = MaterialPhase.FLUID if props.get("phase","solid") == "fluid" else MaterialPhase.SOLID
        k_val = props.get("k")
        if k_val is None:
            raise ValueError("material " + repr(mat_nome) + " sem k definido.")
        mat = Material(name=mat_nome, phase=ph, k=k_val,
                       rho=props.get("rho"), cp=props.get("cp"), mu=props.get("mu"))
    else:
        k_val = fnum(kw["k"])
        mat   = Material(name="mat_" + nome, phase=MaterialPhase.SOLID, k=k_val)
    geom.material(nome, mat,
                  x0=fnum(kw["x0"]), x1=fnum(kw["x1"]),
                  y0=fnum(kw["y0"]), y1=fnum(kw["y1"]))
    return "  -> regiao " + repr(nome) + " k=" + str(round(mat.k, 6))
def cmd_source(state, args):
    geom=_req(state); kw=parse_kwargs(args[1:])
    geom.source(region=args[0],kind="volumetric",value=fnum(kw["q"]))
    return f"  -> fonte {kw['q']} W/m3 em {args[0]!r}"
def cmd_bc(state, args):
    geom=_req(state); kw=parse_kwargs(args[2:],lower_keys=False)
    geom.bc(args[0].lower(),args[1].lower(),**{k:fnum(v) for k,v in kw.items()})
    return f"  -> BC {args[1]!r} em {args[0]!r}"
def cmd_mesh(state, args):
    geom=_req(state); kw=parse_kwargs(args)
    geom.mesh(nx=int(fnum(kw["nx"])),ny=int(fnum(kw["ny"])),fix=kw.get("fix","auto"))
    return f"  -> malha nx={geom.nx} ny={geom.ny}"
def cmd_build_from_geom(state, args):
    geom=_req(state); kw=parse_kwargs(args) if args else {}
    state.net=build_network_from_geometry(geom,thickness_z=fnum(kw.get("tz","1"))); state.ids.clear()
    return f"  -> NodalNetwork: {len(state.net.nodes)} nos, {len(state.net.links)} ligacoes"

def _chk(): 
    if not _REATOR_DISPONIVEL: raise RuntimeError("reator_placa nao encontrado. Execute da raiz.")

def cmd_reator_config(state, args):
    _chk(); kw=parse_kwargs(args); b=state.reator_cfg or ConfigCaso()
    n = int(fnum(kw["n_axial"])) if "n_axial" in kw else b.n_axial
    v = fnum(kw["vazao"])         if "vazao"   in kw else b.vazao_canal_m3_s
    p = fnum(kw["potencia"])      if "potencia"in kw else b.P_placa_W
    Ti= fnum(kw["t_in"])          if "t_in"    in kw else b.T_in_C
    Pi= fnum(kw["p_in"])          if "p_in"    in kw else b.P_in_Pa
    dP= fnum(kw["dp"])            if "dp"      in kw else b.dP_canal_Pa
    tl= fnum(kw["tol"])           if "tol"     in kw else b.tol
    it= int(fnum(kw["iter"]))     if "iter"    in kw else b.max_iter
    mr= kw.get("modo", b.modo_fluxo).lower()
    modo = MODO_COS if mr in ("cos","cossenoidal") else MODO_CONST if mr in ("constante","const") else (_ for _ in ()).throw(ValueError(f"modo invalido: {mr!r}"))
    state.reator_cfg = ConfigCaso(n_axial=n,vazao_canal_m3_s=v,P_placa_W=p,
        modo_fluxo=modo,T_in_C=Ti,P_in_Pa=Pi,dP_canal_Pa=dP,tol=tl,max_iter=it)
    state.reator_result = None
    return f"  -> ConfigCaso: n_axial={n}, vazao={v:.3e} m3/s, potencia={p:.1f} W, modo={modo}, T_in={Ti} C"

def cmd_reator_geom(state, args):
    _chk(); kw=parse_kwargs(args); b=state.reator_geom or GeometriaReator()
    state.reator_geom = GeometriaReator(
        Lx    =fnum(kw["lx"])     if "lx"     in kw else b.Lx,
        Ly    =fnum(kw["ly"])     if "ly"     in kw else b.Ly,
        df    =fnum(kw["df"])     if "df"     in kw else b.df,
        dcl   =fnum(kw["dcl"])    if "dcl"    in kw else b.dcl,
        dch   =fnum(kw["dch"])    if "dch"    in kw else b.dch,
        Lcanal=fnum(kw["lcanal"]) if "lcanal" in kw else b.Lcanal)
    g=state.reator_geom; state.reator_result=None
    return f"  -> GeometriaReator: Lx={g.Lx*1e3:.1f}mm df={g.df*1e3:.2f}mm dcl={g.dcl*1e3:.2f}mm dch={g.dch*1e3:.2f}mm Dh={g.Dh*1e3:.3f}mm"

def cmd_reator_solve(state):
    import numpy as np
    _chk()
    cfg = state.reator_cfg or ConfigCaso()
    g   = state.reator_geom or GeometriaReator()
    rede = construir_rede(cfg, g)
    z0   = gerar_chute_inicial(rede)
    sol  = solve_steady_state(rede.net, z0=z0, tol=cfg.tol, max_iter=cfg.max_iter,
                              update_network=True, prefer_scipy=True)
    res  = extrair_resultados(rede)
    res["solver_success"]       = bool(sol.success)
    res["solver_residual_norm"] = float(sol.residual_norm)
    res["solver_iterations"]    = sol.iterations
    res["n_nos"]                = len(rede.net.nodes)
    res["n_links"]              = len(rede.net.links)
    state.reator_result = res
    Tch_out = float(res["Tch"][-1])
    Tf_max  = float(np.max(res["Tf"]))
    return (f"  -> convergiu={sol.success}, |R|={sol.residual_norm:.2e}, iter={sol.iterations}\n"
            f"  -> Tch_saida={Tch_out:.2f} C, Tf_max={Tf_max:.2f} C")

def cmd_reator_show(state):
    _chk()
    if state.reator_result is None: return "  sem resultado. Execute reator_solve."
    return texto_resumo(state.reator_result, state.reator_result["cfg"].modo_fluxo)

def cmd_reator_plot(state, args):
    _chk()
    if state.reator_result is None: return "  sem resultado. Execute reator_solve."
    kw = parse_kwargs(args) if args else {}
    outdir = Path(kw.get("dir", "saidas_reator_placa")); outdir.mkdir(parents=True, exist_ok=True)
    res = state.reator_result; modo = res["cfg"].modo_fluxo
    plot_temperaturas(res, outdir/f"temperaturas_{modo}.png", f"EC placa q'' {modo}")
    plot_fluxo(res, outdir/f"fluxo_{modo}.png", f"Fluxo axial q'' {modo}")
    return f"  -> salvo em {outdir}: temperaturas_{modo}.png, fluxo_{modo}.png"

def cmd_reator_csv(state, args):
    _chk()
    if state.reator_result is None: return "  sem resultado. Execute reator_solve."
    kw = parse_kwargs(args) if args else {}
    outdir = Path(kw.get("dir", "saidas_reator_placa")); outdir.mkdir(parents=True, exist_ok=True)
    res = state.reator_result; modo = res["cfg"].modo_fluxo
    path = outdir/f"resultado_{modo}.csv"; salvar_csv(res, path)
    return f"  -> CSV salvo em {path}"


def cmd_reator_gif(state, args):
    """Gera GIF animado da rede nodal: campo 2D de temperatura com varredura axial."""
    _chk()
    if state.reator_result is None:
        return "  sem resultado. Execute reator_solve primeiro."
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import Normalize

    kw  = parse_kwargs(args) if args else {}
    outdir = Path(kw.get("dir", "saidas_reator_placa"))
    outdir.mkdir(parents=True, exist_ok=True)
    fps = int(fnum(kw.get("fps", "12")))
    dpi = int(fnum(kw.get("dpi", "100")))

    res  = state.reator_result
    cfg  = res["cfg"]
    modo = cfg.modo_fluxo
    N    = cfg.n_axial
    z    = res["z"] * 100.0  # cm (entrada=0)

    # --- Monta campo 2D: shape (9, N), de baixo para cima transversalmente ---
    # Linha 0 = centro do combustivel inferior (simetria -> igual ao superior)
    y_labels = ["comb.(inf)", "interf.(inf)", "clad(inf)", "sup.(inf)",
                 "fluido", "sup.(sup)", "clad(sup)", "interf.(sup)", "comb.(sup)"]
    T2d = np.array([res["Tf"], res["Ti"], res["Tcl"], res["Ts"], res["Tch"],
                    res["Ts"], res["Tcl"], res["Ti"],  res["Tf"]])  # (9, N)

    T_lo, T_hi = T2d.min(), T2d.max()
    norm = Normalize(vmin=T_lo, vmax=T_hi)
    cmap = plt.cm.hot

    # --- Layout ---
    fig = plt.figure(figsize=(13, 5))
    ax_map  = fig.add_axes([0.06, 0.12, 0.58, 0.78])   # mapa 2D
    ax_prof = fig.add_axes([0.70, 0.12, 0.18, 0.78])   # perfil transversal
    ax_cb   = fig.add_axes([0.90, 0.12, 0.025, 0.78])  # colorbar

    # Mapa 2D completo (exibido do inicio — animamos o scanner)
    img = ax_map.imshow(T2d, aspect="auto", origin="lower", cmap=cmap, norm=norm,
                        extent=[z[0], z[-1], -0.5, 8.5], interpolation="nearest")
    ax_map.set_xlabel("Posicao axial desde a entrada [cm]", fontsize=9)
    ax_map.set_yticks(range(9))
    ax_map.set_yticklabels(y_labels, fontsize=7)
    ax_map.set_title(f"Campo de temperatura nodal  |  q\'\' {modo}  |  {N} camadas axiais", fontsize=9)

    plt.colorbar(img, cax=ax_cb, label="T [°C]")

    # Linha vertical de varredura
    vline = ax_map.axvline(x=z[0], color="cyan", lw=1.8, ls="--", alpha=0.85)

    # Painel de perfil transversal
    ax_prof.set_xlim(T_lo - 0.5, T_hi + 0.5)
    ax_prof.set_ylim(-0.5, 8.5)
    ax_prof.set_yticks(range(9))
    ax_prof.set_yticklabels(y_labels, fontsize=7)
    ax_prof.set_xlabel("T [°C]", fontsize=8)
    ax_prof.set_title("Perfil\nz = --", fontsize=8)
    ax_prof.grid(True, alpha=0.3)
    pts_prof = ax_prof.scatter(T2d[:, 0], range(9), c=T2d[:, 0],
                               cmap=cmap, norm=norm, s=60, zorder=5)
    line_prof, = ax_prof.plot(T2d[:, 0], range(9), color="gray", lw=1, alpha=0.6)

    # Texto de temperatura do fluido
    txt_tch = ax_map.text(0.02, 0.97, "", transform=ax_map.transAxes,
                          fontsize=8, va="top", color="cyan",
                          bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.5))

    # --- Frames: subsample para GIF leve ---
    n_frames = min(N, 60)
    frame_idx = np.linspace(0, N - 1, n_frames, dtype=int)
    # adiciona ultimo frame duplicado para pausa visual
    frame_idx = np.concatenate([frame_idx, np.full(6, N - 1, dtype=int)])

    def update(fi):
        j = frame_idx[fi]
        vline.set_xdata([z[j]])
        col = T2d[:, j]
        pts_prof.set_offsets(np.column_stack([col, range(9)]))
        pts_prof.set_array(col)
        line_prof.set_xdata(col)
        ax_prof.set_title("Perfil\nz = " + f"{z[j]:.1f}" + " cm", fontsize=8)
        tch_j = float(res["Tch"][j])
        txt_tch.set_text("Tch = " + f"{tch_j:.2f}" + " degC")
        return vline, pts_prof, line_prof, txt_tch

    ani = animation.FuncAnimation(fig, update, frames=len(frame_idx),
                                   interval=1000 // fps, blit=True)
    gif_path = outdir / f"animacao_{modo}.gif"
    ani.save(str(gif_path), writer="pillow", dpi=dpi)
    plt.close(fig)
    return f"  -> GIF salvo em {gif_path}  ({len(frame_idx)} frames, {fps} fps)"


# ===========================================================================
# Transporte entalpico generico
# ===========================================================================

def cmd_fluid_chain(state, args):
    """fluid_chain n1 n2 ... mdot=v cp=v [T_in=v]
    Cria links FLUID_TRANSPORT direcionais n1->n2->n3->...
    Se T_in for dado, adiciona no BOUNDARY de entrada ligado a n1.
    """
    from nos import LinkDirection, TransferKind
    nos_names = [a for a in args if "=" not in a]
    kw = parse_kwargs([a for a in args if "=" in a])
    if len(nos_names) < 2:
        raise ValueError("fluid_chain precisa de pelo menos 2 nos: fluid_chain n1 n2 ... mdot=v cp=v")
    if "mdot" not in kw or "cp" not in kw:
        raise ValueError("fluid_chain precisa de mdot=v e cp=v")
    mdot = fnum(kw["mdot"])
    cp   = fnum(kw["cp"])
    G    = mdot * cp

    T_in_val = fnum(kw["t_in"]) if "t_in" in kw else None

    def resolve(name):
        if name in state.ids:
            return state.ids[name]
        try:
            nid = int(name)
            if nid in state.net.nodes:
                return nid
        except ValueError:
            pass
        raise ValueError("No nao encontrado: " + repr(name))

    chain = [resolve(n) for n in nos_names]
    msgs = []

    if T_in_val is not None:
        first_node = state.net.nodes[chain[0]]
        inlet_id = state.net.add_node(
            name="inlet_" + nos_names[0], kind=NodeKind.BOUNDARY,
            x=first_node.x - 1.0, y=first_node.y,
            fixed_temperature=T_in_val, temperature=T_in_val)
        state.net.add_link(inlet_id, chain[0], TransferKind.FLUID_TRANSPORT,
                           direction=LinkDirection.I_TO_J, conductance=G,
                           name="fluid_inlet_" + nos_names[0])
        msgs.append("  -> inlet BOUNDARY T=" + str(T_in_val) + " C -> " + nos_names[0])

    for k in range(len(chain) - 1):
        na, nb = chain[k], chain[k+1]
        lname = "fluid_" + nos_names[k] + "_" + nos_names[k+1]
        state.net.add_link(na, nb, TransferKind.FLUID_TRANSPORT,
                           direction=LinkDirection.I_TO_J, conductance=G,
                           name=lname)
        msgs.append("  -> " + nos_names[k] + " -> " + nos_names[k+1] + "  G=" + str(round(G,4)) + " W/K")

    return "\n".join(msgs)


# ===========================================================================
# Visualizacao generica
# ===========================================================================

def _render_network(net, title="Rede Nodal", fig=None, ax=None, alpha_links=0.55):
    """Desenha a NodalNetwork como grafo 2D colorido por temperatura."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from nos import LinkDirection, TransferKind, NodeKind

    ids   = list(net.nodes.keys())
    xs    = np.array([net.nodes[i].x for i in ids])
    ys    = np.array([net.nodes[i].y for i in ids])
    Ts    = np.array([net.nodes[i].temperature for i in ids])
    kinds = [net.nodes[i].kind for i in ids]

    T_lo = float(Ts.min())
    T_hi = float(Ts.max())
    if abs(T_hi - T_lo) < 0.01:
        T_hi = T_lo + 1.0
    norm = Normalize(vmin=T_lo, vmax=T_hi)
    cmap = plt.cm.hot

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(11, 6))

    link_style = {
        TransferKind.CONDUCTION:      ("#888888", "-",  0.9),
        TransferKind.CONVECTION:      ("#4499ff", "--", 1.0),
        TransferKind.FLUID_TRANSPORT: ("#ff4444", "-",  1.4),
        TransferKind.RADIATION:       ("#cc44ff", ":",  1.0),
        TransferKind.EQUIVALENT:      ("#aaaaaa", "-.", 0.8),
    }
    id_pos = {i: (net.nodes[i].x, net.nodes[i].y) for i in ids}
    for lk in net.links:
        xi, yi = id_pos[lk.node_i]
        xj, yj = id_pos[lk.node_j]
        col, ls, lw = link_style.get(lk.kind, ("#888888", "-", 0.9))
        ax.plot([xi, xj], [yi, yj], color=col, ls=ls, lw=lw,
                alpha=alpha_links, zorder=1)
        if lk.direction != LinkDirection.UNDIRECTED:
            mx = 0.5*(xi+xj); my = 0.5*(yi+yj)
            dx = xj-xi; dy = yj-yi
            if lk.direction == LinkDirection.J_TO_I:
                dx = -dx; dy = -dy
            norm_d = (dx**2+dy**2)**0.5
            if norm_d > 1e-9:
                ax.annotate("",
                    xy=(mx + dx/norm_d*0.12, my + dy/norm_d*0.12),
                    xytext=(mx, my),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.3),
                    zorder=2)

    mk_map = {
        "diffusion":  ("o", 90),
        "fluid":      ("s", 90),
        "arithmetic": ("^", 55),
        "boundary":   ("D", 65),
    }
    for kv in set(k.value for k in kinds):
        mask = [k.value == kv for k in kinds]
        idx  = [i for i,m in enumerate(mask) if m]
        xm = xs[idx]; ym = ys[idx]; Tm = Ts[idx]
        mk, sz = mk_map.get(kv, ("o", 80))
        ax.scatter(xm, ym, c=Tm, cmap=cmap, norm=norm,
                   marker=mk, s=sz, edgecolors="black", lw=0.6,
                   zorder=4, label=kv)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="T [degC]", shrink=0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(title + "  |  " + str(len(ids)) + " nos  " + str(len(net.links)) + " links")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.grid(True, alpha=0.18)
    return fig, ax, norm, cmap


def cmd_viz(state, args):
    """Visualiza a rede nodal atual: nos coloridos por temperatura, links tipados."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if len(state.net.nodes) == 0:
        return "  rede vazia. Adicione nos primeiro."
    kw      = parse_kwargs(args) if args else {}
    outfile = kw.get("png", "rede_nodal.png")
    titulo  = kw.get("titulo", "Rede Nodal")

    fig, ax, _, _ = _render_network(state.net, title=titulo)
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=130, bbox_inches="tight")
    plt.close(fig)
    return "  -> viz salvo em " + str(outpath)


def cmd_gif_generico(state, args):
    """GIF animado: pseudo-transiente do estado frio ate o regime permanente."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from nos import LinkDirection, TransferKind

    if len(state.net.nodes) == 0:
        return "  rede vazia."
    kw      = parse_kwargs(args) if args else {}
    outfile = kw.get("gif", "rede_nodal.gif")
    fps     = int(fnum(kw.get("fps", "10")))
    dpi     = int(fnum(kw.get("dpi", "100")))
    n_steps = int(fnum(kw.get("steps", "40")))
    titulo  = kw.get("titulo", "Rede Nodal")

    net  = state.net
    ids  = list(net.nodes.keys())
    xs   = np.array([net.nodes[i].x for i in ids])
    ys   = np.array([net.nodes[i].y for i in ids])
    T_final = np.array([net.nodes[i].temperature for i in ids])
    kinds = [net.nodes[i].kind for i in ids]

    T_lo = float(T_final.min())
    T_hi = float(T_final.max())
    if abs(T_hi - T_lo) < 0.01:
        T_hi = T_lo + 1.0
    norm = Normalize(vmin=T_lo, vmax=T_hi)
    cmap = plt.cm.hot

    is_boundary = np.array([net.nodes[i].is_boundary() for i in ids])
    T_start = np.where(is_boundary, T_final, T_lo)

    alphas = np.linspace(0.0, 1.0, n_steps)
    alphas = np.concatenate([alphas, np.ones(6)])

    fig, ax = plt.subplots(figsize=(11, 6))

    link_style = {
        TransferKind.CONDUCTION:      ("#888888", "-",  0.8),
        TransferKind.CONVECTION:      ("#4499ff", "--", 0.9),
        TransferKind.FLUID_TRANSPORT: ("#ff4444", "-",  1.3),
        TransferKind.RADIATION:       ("#cc44ff", ":",  0.9),
        TransferKind.EQUIVALENT:      ("#aaaaaa", "-.", 0.7),
    }
    id_pos = {i: (net.nodes[i].x, net.nodes[i].y) for i in ids}
    from nos import LinkDirection as LD
    for lk in net.links:
        xi, yi = id_pos[lk.node_i]
        xj, yj = id_pos[lk.node_j]
        col, ls, lw = link_style.get(lk.kind, ("#888888", "-", 0.8))
        ax.plot([xi, xj], [yi, yj], color=col, ls=ls, lw=lw, alpha=0.45, zorder=1)
        if lk.direction != LD.UNDIRECTED:
            mx = 0.5*(xi+xj); my = 0.5*(yi+yj)
            dx = xj-xi; dy = yj-yi
            if lk.direction == LD.J_TO_I:
                dx=-dx; dy=-dy
            nd = (dx**2+dy**2)**0.5
            if nd > 1e-9:
                ax.annotate("", xy=(mx+dx/nd*0.12, my+dy/nd*0.12),
                    xytext=(mx, my),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.2), zorder=2)

    mk_map  = {"diffusion":"o","fluid":"s","arithmetic":"^","boundary":"D"}
    sz_map  = {"diffusion":90,"fluid":90,"arithmetic":55,"boundary":65}

    scs = {}
    for kv in set(k.value for k in kinds):
        mask = np.array([k.value == kv for k in kinds])
        T0   = T_start[mask]
        scs[kv] = ax.scatter(xs[mask], ys[mask], c=T0, cmap=cmap, norm=norm,
                              marker=mk_map.get(kv,"o"), s=sz_map.get(kv,80),
                              edgecolors="black", lw=0.6, zorder=4, label=kv)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="T [degC]", shrink=0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(titulo)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.grid(True, alpha=0.18)
    pct_txt = ax.text(0.02, 0.97, "0%", transform=ax.transAxes,
                      fontsize=9, va="top",
                      bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    def update(fi):
        a = float(alphas[fi])
        T_cur = T_start + a * (T_final - T_start)
        for kv in scs:
            mask = np.array([k.value == kv for k in kinds])
            scs[kv].set_array(T_cur[mask])
        pct = int(round(a * 100))
        T_max_cur = float(T_cur.max())
        pct_txt.set_text(str(pct) + "%   Tmax=" + f"{T_max_cur:.1f}" + " degC")
        return list(scs.values()) + [pct_txt]

    ani = animation.FuncAnimation(fig, update, frames=len(alphas),
                                   interval=1000//fps, blit=True)
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(outpath), writer="pillow", dpi=dpi)
    plt.close(fig)
    n_frames = len(alphas)
    return "  -> GIF salvo em " + str(outpath) + "  (" + str(n_frames) + " frames, " + str(fps) + " fps)"

def executar(state, linha):
    raw = linha.strip()
    if not raw or raw.startswith("#"): return "", False
    tokens = shlex.split(raw); cmd = tokens[0].lower(); args = tokens[1:]
    if cmd in ("quit","exit","q"): return "Saindo.", True
    cmds = {
        "help":            lambda: HELP,
        "reset":           lambda: state.reset() or "  rede zerada.",
        "show":            lambda: cmd_show(state),
        "solve":           lambda: cmd_solve(state),
        "node":            lambda: cmd_node(state, args),
        "link":            lambda: cmd_link(state, args),
        "g_cond":          lambda: cmd_gcond(args),
        "g_conv":          lambda: cmd_gconv(args),
        "material":        lambda: cmd_material(state, args),
        "materiais":       lambda: cmd_materiais(state, args),
        "domain":          lambda: cmd_domain(state, args),
        "region":          lambda: cmd_region(state, args),
        "source":          lambda: cmd_source(state, args),
        "bc":              lambda: cmd_bc(state, args),
        "mesh":            lambda: cmd_mesh(state, args),
        "show_geom":       lambda: state.geom.summary() if state.geom else "sem geometria.",
        "build_from_geom": lambda: cmd_build_from_geom(state, args),
        "reator_config":   lambda: cmd_reator_config(state, args),
        "reator_geom":     lambda: cmd_reator_geom(state, args),
        "reator_solve":    lambda: cmd_reator_solve(state),
        "reator_show":     lambda: cmd_reator_show(state),
        "reator_plot":     lambda: cmd_reator_plot(state, args),
        "reator_csv":      lambda: cmd_reator_csv(state, args),
        "reator_gif":      lambda: cmd_reator_gif(state, args),
        "fluid_chain":     lambda: cmd_fluid_chain(state, args),
        "viz":             lambda: cmd_viz(state, args),
        "gif":             lambda: cmd_gif_generico(state, args),
    }
    if cmd not in cmds: raise ValueError(f"comando desconhecido: {cmd!r}. Digite help.")
    return cmds[cmd](), False

def cmd_gcond(args):
    kw=parse_kwargs(args); G=conduction_G(fnum(kw["k"]),fnum(kw["a"]),fnum(kw["l"]))
    return f"  G_cond = {G:.6g} W/K"
def cmd_gconv(args):
    kw=parse_kwargs(args); G=convection_G(fnum(kw["h"]),fnum(kw["a"]))
    return f"  G_conv = {G:.6g} W/K"

def rodar_demo():
    state = EstadoREPL()
    demo = ["# Exemplo chip - item (a)","g_cond k=380 A=7.854e-7 L=0.02",
            "node Tc diffusion Q=4 T_init=30","node Tp boundary T_fixed=44",
            "link Tc Tp cond G=0.17908","solve",
            "# item (b): conveccao","g_conv h=30 A=9e-4",
            "node Tar boundary T_fixed=20","link Tc Tar conv G=0.027","solve"]
    print("="*72, "\nExemplo chip - demo\n" + "="*72)
    for l in demo:
        if l.strip(): print(f"\n>>> {l}")
        try:
            r,_=executar(state,l)
            if r: print(r)
        except Exception as e: print(f"  ERRO: {e}")

def rodar_repl():
    state=EstadoREPL()
    print("="*72+"\nMetodo Nodal - interativo\nDigite help ou quit\n"+"="*72)
    while True:
        try: l=input(">>> ")
        except (EOFError,KeyboardInterrupt): print(); break
        try:
            r,enc=executar(state,l)
            if r: print(r)
            if enc: break
        except Exception as e: print(f"  ERRO: {e}")

def rodar_script(path):
    state=EstadoREPL()
    print("="*72+f"\nMetodo Nodal - script ({Path(path).name})\n"+"="*72)
    for l in Path(path).read_text(encoding="utf-8").splitlines():
        if l.strip(): print(f"\n>>> {l}")
        try:
            r,_=executar(state,l)
            if r: print(r)
        except Exception as e: print(f"  ERRO: {e}")

def main(argv):
    if len(argv)>=2 and argv[1]=="--interativo": rodar_repl()
    elif len(argv)>=3 and argv[1]=="--script": rodar_script(Path(argv[2]))
    else: rodar_demo()

if __name__=="__main__":
    main(sys.argv)
