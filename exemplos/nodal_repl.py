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
                  [k_poly=a0,a1,...] [k_expr="expressao em T"]
       k pode ser constante, polinomio ou expressao de T. Exemplos:
         material cobre  k=380 rho=8900 cp=385
         material agua   phase=fluid k=0.6 rho=997 cp=4182 mu=8.9e-4
         material uo2    k_poly=10.5,-0.012,4e-6 rho=10970 cp=330
         material u3si2  k_expr="15.0+0.002*T" rho=4300 cp=836
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
  nucleo_init [rows=5] [cols=5]
       Inicializa matriz vazia do nucleo para entrada manual
  nucleo_linha row=N vals=v1,v2,... [tipos=EC,CR,...]
       Define uma linha da matriz (row=0 e linha de topo)
  nucleo_mapa [png=mapa.png] [titulo=Texto]
       Gera PNG do mapa de potencia radial (usa matriz inserida ou padrao 5x5)
  nucleo_solve [P_total=5e6] [n_canais=19] [png=mapa_T.png] [plots=no] [dir=saidas_nucleo]
       Resolve todos os ECs e gera mapa 5x5 de temperatura maxima do combustivel
  ec_associacao [n_int=17] [n_ext=2] [dch_int=0.00289] [dch_ext=0.00452]
                [potencia=P] [vazao=V] [png=ec_assoc.png]
       Distribui vazao entre canais internos e externos (mesmo dP, Blasius)

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
        self.nucleo_data = None   # dict com rows, cols, vals, tipos (inserido via nucleo_init/nucleo_linha)

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
    for campo in ("k_poly","rho_poly","cp_poly","mu_poly"):
        if campo in kw:
            props[campo] = [fnum(v) for v in kw[campo].split(",")]
    for campo in ("k_expr","rho_expr","cp_expr","mu_expr"):
        if campo in kw:
            props[campo] = kw[campo]
    state.materials[nome] = props
    partes = []
    if "k_expr" in props:
        partes.append("k=f(T)[" + props["k_expr"] + "]")
    elif "k_poly" in props:
        partes.append("k=poly" + str(props["k_poly"]))
    elif "k" in props:
        partes.append("k=" + str(props["k"]) + " W/(mK)")
    for campo, unid in [("rho","kg/m3"),("cp","J/(kgK)"),("mu","Pa.s")]:
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
        import math as _math
        k_coeffs = props.get("k_poly")
        k_expr   = props.get("k_expr")
        if k_expr:
            _expr = k_expr
            k_func = lambda T, P=None, _e=_expr: float(
                eval(_e, {"T": T, "P": P, "math": _math, **vars(_math)}))
            mat = Material(name=mat_nome, phase=ph, k_func=k_func,
                           rho=props.get("rho"), cp=props.get("cp"), mu=props.get("mu"))
        elif k_coeffs:
            mat = Material(name=mat_nome, phase=ph, k_coeffs=k_coeffs,
                           rho=props.get("rho"), cp=props.get("cp"), mu=props.get("mu"))
        else:
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


# ---------------------------------------------------------------------------
# Mapa de potencia do nucleo
# ---------------------------------------------------------------------------
# Mapa de fatores de potencia radial (5x5) — Projeto Final TNR5703
# Linhas 0-4 de cima para baixo, colunas 0-4 da esquerda para direita
_MAPA_PADRAO = [
    [1.321, 1.563, 0.981, 1.628, 1.030],
    [0.857, 0.515, 1.129, 0.402, 0.826],
    [1.050, 1.877, 0.000, 1.914, 0.979],
    [0.860, 0.411, 1.151, 0.519, 0.822],
    [0.906, 1.028, 0.878, 1.044, 0.867],
]
# Tipos de elemento: "EC" = combustivel, "CR" = controle/refletor
_TIPOS_PADRAO = [
    ["EC","EC","EC","EC","EC"],
    ["EC","CR","EC","CR","EC"],
    ["EC","EC","CR","EC","EC"],
    ["EC","CR","EC","CR","EC"],
    ["EC","EC","EC","EC","EC"],
]


# ---------------------------------------------------------------------------
# Mapa de temperaturas do nucleo e associacao de canais paralelos
# ---------------------------------------------------------------------------

def _solve_one_canal(geom, cfg_template, P_canal):
    import numpy as np
    cfg = ConfigCaso(
        n_axial          = cfg_template.n_axial,
        vazao_canal_m3_s = cfg_template.vazao_canal_m3_s,
        P_placa_W        = P_canal,
        modo_fluxo       = cfg_template.modo_fluxo,
        T_in_C           = cfg_template.T_in_C,
        P_in_Pa          = cfg_template.P_in_Pa,
        dP_canal_Pa      = cfg_template.dP_canal_Pa,
        tol              = cfg_template.tol,
        max_iter         = cfg_template.max_iter,
    )
    rede = construir_rede(cfg, geom)
    z0   = gerar_chute_inicial(rede)
    sol  = solve_steady_state(rede.net, z0=z0, tol=cfg.tol, max_iter=cfg.max_iter,
                              update_network=True, prefer_scipy=True)
    return extrair_resultados(rede)

def _gerar_mapa_temperaturas(mapa_fat, tipos, T_fuel, nrows, ncols, outfile, titulo):
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    vals = [T_fuel[i][j] for i in range(nrows) for j in range(ncols)
            if T_fuel[i][j] == T_fuel[i][j]]
    T_min_v = min(vals) if vals else 30.0
    T_max_v = max(vals) if vals else 100.0
    pos_max = (0, 0); T_abs_max = 0.0
    for i in range(nrows):
        for j in range(ncols):
            v = T_fuel[i][j]
            if v == v and v > T_abs_max:
                T_abs_max = v; pos_max = (i, j)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_xlim(0, ncols); ax.set_ylim(0, nrows)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(titulo, fontsize=12, fontweight="bold", pad=14)
    cmap = plt.cm.plasma
    for i in range(nrows):
        for j in range(ncols):
            row_plot = nrows - 1 - i
            tipo = tipos[i][j]; fat = mapa_fat[i][j]; T = T_fuel[i][j]
            if fat == 0.0:
                face="#c0c0c0"; txt="---"; tcol="#333"
            elif tipo == "CR":
                face="#8ab4d4"; txt="CR"; tcol="white"
            elif T != T:
                face="#eeeeee"; txt="ERR"; tcol="#c00"
            else:
                norm = max(0.0, min(1.0, (T - T_min_v)/(T_max_v - T_min_v + 1e-12)))
                rgba = cmap(norm)
                face = "#{:02x}{:02x}{:02x}".format(
                    int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
                txt = f"{T:.1f}"; tcol = "white" if norm > 0.5 else "#222"
            rect = Rectangle((j, row_plot), 1, 1,
                              facecolor=face, edgecolor="white", linewidth=2.5)
            ax.add_patch(rect)
            if (i, j) == pos_max:
                ax.add_patch(Rectangle((j+0.04, row_plot+0.04), 0.92, 0.92,
                    fill=False, edgecolor="#cc0000", linewidth=3.0))
            ax.text(j+0.5, row_plot+0.5, txt, ha="center", va="center",
                    fontsize=11, fontweight="bold", color=tcol)
    sm = plt.cm.ScalarMappable(cmap=cmap,
         norm=plt.Normalize(vmin=T_min_v, vmax=T_max_v))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04)
    cbar.set_label("T combustivel max [C]", fontsize=9)
    fig.tight_layout()
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=150, bbox_inches="tight")
    plt.close(fig)

def cmd_nucleo_solve(state, args):
    if not _REATOR_DISPONIVEL:
        raise RuntimeError("modulo reator_placa nao disponivel.")
    if state.nucleo_data is None:
        raise ValueError("use nucleo_init + nucleo_linha para definir o mapa.")
    if state.reator_geom is None:
        raise ValueError("use reator_geom para definir a geometria.")
    if state.reator_cfg is None:
        raise ValueError("use reator_config para definir a configuracao.")
    import numpy as np
    kw       = parse_kwargs(args) if args else {}
    P_total  = fnum(kw.get("p_total",  "5e6"))
    n_canais = int(fnum(kw.get("n_canais", "19")))
    outfile  = kw.get("png", "saidas/mapa_temperaturas.png")
    do_plots = kw.get("plots", "no").lower() in ("yes", "sim", "1", "true")
    plotdir  = Path(kw.get("dir", "saidas_nucleo"))
    mapa  = state.nucleo_data["vals"]
    tipos = state.nucleo_data["tipos"]
    nrows = state.nucleo_data["rows"]
    ncols = state.nucleo_data["cols"]
    n_ec = sum(1 for i in range(nrows) for j in range(ncols)
               if tipos[i][j] == "EC" and mapa[i][j] > 0.0)
    if n_ec == 0:
        raise ValueError("nenhum EC ativo no mapa.")
    P_base_canal = P_total / n_ec / n_canais
    geom = state.reator_geom; cfg = state.reator_cfg
    T_fuel = [[float("nan")]*ncols for _ in range(nrows)]
    if do_plots:
        plotdir.mkdir(parents=True, exist_ok=True)
        import matplotlib; matplotlib.use("Agg")
        from reator_placa.exemplo_reator_placa import plot_temperaturas as _pt
    log = []
    for i in range(nrows):
        for j in range(ncols):
            fator = mapa[i][j]; tipo = tipos[i][j]
            if tipo != "EC" or fator <= 0.0:
                continue
            try:
                res = _solve_one_canal(geom, cfg, fator * P_base_canal)
                Tf  = float(np.max(res["Tf"]))
                Ts  = float(np.max(res["Ts"]))
                Tch = float(res["Tch"][-1])
                T_fuel[i][j] = Tf
                log.append(f"  EC({i+1},{j+1}) f={fator:.3f}  Tf_max={Tf:.1f}C  Ts_max={Ts:.1f}C  Tch_out={Tch:.1f}C")
                if do_plots:
                    _pt(res, plotdir/f"ec_{i+1}_{j+1}.png",
                        f"EC ({i+1},{j+1}) fator={fator:.3f}")
            except Exception as e:
                log.append(f"  EC({i+1},{j+1}) ERRO: {e}")
    _gerar_mapa_temperaturas(mapa, tipos, T_fuel, nrows, ncols, outfile,
                             titulo="T max combustivel por EC [C]")
    log.append(f"  -> mapa salvo em {outfile}")
    if do_plots:
        log.append(f"  -> plots em {plotdir}/")
    # GIF: scan EC a EC revelando o mapa conforme vai sendo resolvido
    gif_file = kw.get("gif", "")
    if gif_file:
        r = _gif_nucleo_scan(mapa, tipos, T_fuel, nrows, ncols, gif_file,
                             fps=int(fnum(kw.get("fps","6"))))
        log.append(r)
    return "\n".join(log)

def _plot_ec_associacao(res_int, res_ext, n_int, n_ext, outfile, titulo):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(titulo, fontsize=12, fontweight="bold")
    for ax, res, label in [
        (axes[0], res_int, f"Canal interno  (x{n_int})"),
        (axes[1], res_ext, f"Canal externo  (x{n_ext})"),
    ]:
        z = res["z"] * 100
        ax.plot(z, res["Tf"],  "r-",  lw=2,   label="T combustivel")
        ax.plot(z, res["Ti"],  "m--", lw=1.5, label="T interface")
        ax.plot(z, res["Tcl"], "b--", lw=1.5, label="T revestimento")
        ax.plot(z, res["Ts"],  "g-",  lw=1.5, label="T superficie")
        ax.plot(z, res["Tch"], "c-",  lw=2,   label="T fluido")
        ax.set_xlabel("Posicao axial z [cm]")
        ax.set_ylabel("Temperatura [C]")
        ax.set_title(label); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=130, bbox_inches="tight")
    plt.close(fig)

def cmd_ec_associacao(state, args):
    if not _REATOR_DISPONIVEL:
        raise RuntimeError("modulo reator_placa nao disponivel.")
    if state.reator_geom is None:
        raise ValueError("use reator_geom primeiro.")
    if state.reator_cfg is None:
        raise ValueError("use reator_config primeiro.")
    import numpy as np
    from reator_placa.geometria_reator import GeometriaReator as _GR
    kw      = parse_kwargs(args) if args else {}
    n_int   = int(fnum(kw.get("n_int",   "17")))
    n_ext   = int(fnum(kw.get("n_ext",    "2")))
    dch_int = fnum(kw.get("dch_int", str(state.reator_geom.dch)))
    dch_ext = fnum(kw.get("dch_ext", "0.00452"))
    Lcanal  = fnum(kw.get("lcanal",  str(state.reator_geom.Lcanal)))
    P_EC    = fnum(kw.get("potencia", str(state.reator_cfg.P_placa_W)))
    V_EC    = fnum(kw.get("vazao",    str(state.reator_cfg.vazao_canal_m3_s)))
    outfile = kw.get("png", "saidas/ec_associacao.png")
    geom    = state.reator_geom; cfg = state.reator_cfg
    def dh_canal(dch, Lc):
        A = dch * Lc; P = 2*(dch + Lc); return 4*A/P, A
    Dh_int, A_int = dh_canal(dch_int, Lcanal)
    Dh_ext, A_ext = dh_canal(dch_ext, Lcanal)
    w_int = A_int * Dh_int**(5/7)
    w_ext = A_ext * Dh_ext**(5/7)
    W_tot = n_int * w_int + n_ext * w_ext
    V_int = V_EC * (w_int / W_tot)
    V_ext = V_EC * (w_ext / W_tot)
    n_pl  = n_int + 1
    P_int = P_EC / n_pl
    P_ext = P_EC / (2 * n_pl)
    linhas = [
        f"  Canais internos: {n_int}  dch={dch_int*1e3:.2f}mm  Dh={Dh_int*1e3:.3f}mm",
        f"  Canais externos: {n_ext}  dch={dch_ext*1e3:.2f}mm  Dh={Dh_ext*1e3:.3f}mm",
        f"  Vazao  V_int={V_int*1e6:.2f}cm3/s ({V_int/V_EC*100:.1f}%)  V_ext={V_ext*1e6:.2f}cm3/s ({V_ext/V_EC*100:.1f}%)",
        f"  Potencia  P_int={P_int:.1f}W/canal  P_ext={P_ext:.1f}W/canal",
        f"  Resolvendo...",
    ]
    geom_int = _GR(Lx=geom.Lx, Ly=geom.Ly, df=geom.df, dcl=geom.dcl, dch=dch_int, Lcanal=Lcanal)
    geom_ext = _GR(Lx=geom.Lx, Ly=geom.Ly, df=geom.df, dcl=geom.dcl, dch=dch_ext, Lcanal=Lcanal)
    cfg_int = ConfigCaso(n_axial=cfg.n_axial, vazao_canal_m3_s=V_int, P_placa_W=P_int,
        modo_fluxo=cfg.modo_fluxo, T_in_C=cfg.T_in_C, P_in_Pa=cfg.P_in_Pa,
        dP_canal_Pa=cfg.dP_canal_Pa, tol=cfg.tol, max_iter=cfg.max_iter)
    cfg_ext = ConfigCaso(n_axial=cfg.n_axial, vazao_canal_m3_s=V_ext, P_placa_W=P_ext,
        modo_fluxo=cfg.modo_fluxo, T_in_C=cfg.T_in_C, P_in_Pa=cfg.P_in_Pa,
        dP_canal_Pa=cfg.dP_canal_Pa, tol=cfg.tol, max_iter=cfg.max_iter)
    res_int = _solve_one_canal(geom_int, cfg_int, P_int)
    res_ext = _solve_one_canal(geom_ext, cfg_ext, P_ext)
    Tf_i  = float(np.max(res_int["Tf"]));  Ts_i = float(np.max(res_int["Ts"]))
    Tch_i = float(res_int["Tch"][-1])
    Tf_e  = float(np.max(res_ext["Tf"]));  Ts_e = float(np.max(res_ext["Ts"]))
    Tch_e = float(res_ext["Tch"][-1])
    linhas += [
        f"  Canal INTERNO: Tf_max={Tf_i:.2f}C   Ts_max={Ts_i:.2f}C   Tch_saida={Tch_i:.2f}C",
        f"  Canal EXTERNO: Tf_max={Tf_e:.2f}C   Ts_max={Ts_e:.2f}C   Tch_saida={Tch_e:.2f}C",
    ]
    _plot_ec_associacao(res_int, res_ext, n_int, n_ext, outfile,
                        f"Associacao de canais  P={P_EC:.0f}W  V={V_EC*1e6:.0f}cm3/s")
    linhas.append(f"  -> figura salva em {outfile}")
    gif_file = kw.get("gif", "")
    if gif_file:
        r = _gif_ec_associacao(res_int, res_ext, n_int, n_ext, gif_file,
                               fps=int(fnum(kw.get("fps","10"))),
                               steps=int(fnum(kw.get("steps","40"))))
        linhas.append(r)
    state.reator_result = res_int
    return "\n".join(linhas)


def _gif_nucleo_scan(mapa_fat, tipos, T_fuel, nrows, ncols, outfile, fps=6):
    """GIF: revela o mapa de temperaturas EC por EC, na ordem em que foram resolvidos."""
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.patches import Rectangle

    vals = [T_fuel[i][j] for i in range(nrows) for j in range(ncols)
            if T_fuel[i][j] == T_fuel[i][j]]
    T_min_v = min(vals) if vals else 30.0
    T_max_v = max(vals) if vals else 100.0
    pos_max = (0,0); T_abs = 0.0
    for i in range(nrows):
        for j in range(ncols):
            v = T_fuel[i][j]
            if v == v and v > T_abs:
                T_abs = v; pos_max = (i,j)

    # Ordem de resolucao: ECs ativos linha por linha
    ordem = [(i,j) for i in range(nrows) for j in range(ncols)
             if tipos[i][j] == "EC" and mapa_fat[i][j] > 0.0]

    cmap = plt.cm.plasma

    def draw_frame(ax, revelados):
        ax.cla()
        ax.set_xlim(0, ncols); ax.set_ylim(0, nrows)
        ax.set_aspect("equal"); ax.axis("off")
        for i in range(nrows):
            for j in range(ncols):
                row_plot = nrows - 1 - i
                fat = mapa_fat[i][j]; tipo = tipos[i][j]; T = T_fuel[i][j]
                if fat == 0.0:
                    face="#c0c0c0"; txt="---"; tcol="#333"
                elif tipo == "CR":
                    face="#8ab4d4"; txt="CR"; tcol="white"
                elif (i,j) not in revelados:
                    face="#f0f0f0"; txt=f"{fat:.3f}"; tcol="#aaa"
                elif T != T:
                    face="#eeeeee"; txt="ERR"; tcol="#c00"
                else:
                    norm = max(0.0, min(1.0,(T-T_min_v)/(T_max_v-T_min_v+1e-12)))
                    rgba = cmap(norm)
                    face = "#{:02x}{:02x}{:02x}".format(
                        int(rgba[0]*255),int(rgba[1]*255),int(rgba[2]*255))
                    txt = f"{T:.1f}"; tcol = "white" if norm>0.5 else "#222"
                ax.add_patch(Rectangle((j,row_plot),1,1,
                    facecolor=face,edgecolor="white",linewidth=2))
                if (i,j)==pos_max and (i,j) in revelados:
                    ax.add_patch(Rectangle((j+0.04,row_plot+0.04),0.92,0.92,
                        fill=False,edgecolor="#cc0000",linewidth=3))
                ax.text(j+0.5,row_plot+0.5,txt,ha="center",va="center",
                        fontsize=10,fontweight="bold",color=tcol)

    fig, ax = plt.subplots(figsize=(6,6))
    sm = plt.cm.ScalarMappable(cmap=cmap,
         norm=plt.Normalize(vmin=T_min_v,vmax=T_max_v))
    sm.set_array([]); fig.colorbar(sm,ax=ax,fraction=0.03,pad=0.03).set_label("Tf max [C]",fontsize=8)

    frames = []
    # Frame inicial: tudo cinza
    frames.append(set())
    # Um frame por EC revelado
    rev = set()
    for pos in ordem:
        rev = rev | {pos}
        frames.append(frozenset(rev))
    # 3 frames finais com tudo revelado
    for _ in range(3):
        frames.append(frozenset(rev))

    def animate(k):
        ax.set_title(f"Mapa de T combustivel — EC {k}/{len(ordem)}", fontsize=11, fontweight="bold")
        draw_frame(ax, frames[k])
        return []

    ani = animation.FuncAnimation(fig, animate, frames=len(frames),
                                  interval=1000//fps, blit=False)
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(outpath), writer="pillow", fps=fps)
    plt.close(fig)
    return f"  -> GIF nucleo salvo em {outpath}  ({len(frames)} frames, {fps} fps)"


def _gif_ec_associacao(res_int, res_ext, n_int, n_ext, outfile, fps=10, steps=40):
    """GIF: pseudo-transiente dos dois tipos de canal do EC (interno e externo)."""
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    T_in = res_int["Tch"][0]  # temperatura de entrada (cold)

    # Interpola de frio (T_in) ate regime permanente
    alphas = np.linspace(0.0, 1.0, steps)
    # Adiciona alguns frames extras no final
    alphas = np.concatenate([alphas, np.ones(6)])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Associacao de canais — evolucao termica", fontsize=12, fontweight="bold")

    curves = [
        ("Tf",  "r-",  2.0,  "T combustivel"),
        ("Ti",  "m--", 1.5,  "T interface"),
        ("Tcl", "b--", 1.5,  "T revestimento"),
        ("Ts",  "g-",  1.5,  "T superficie"),
        ("Tch", "c-",  2.0,  "T fluido"),
    ]

    # Pre-calcula limites de T para eixo fixo
    all_T = []
    for res in (res_int, res_ext):
        for k,_,_,_ in curves:
            all_T.extend(res[k].tolist())
    T_lo = min(all_T) - 1; T_hi = max(all_T) + 1

    lines_int = {}; lines_ext = {}
    for ax, res, lns, label in [
        (axes[0], res_int, lines_int, f"Canal interno  (x{n_int})"),
        (axes[1], res_ext, lines_ext, f"Canal externo  (x{n_ext})"),
    ]:
        z = res["z"] * 100
        ax.set_xlim(z[0], z[-1]); ax.set_ylim(T_lo, T_hi)
        ax.set_xlabel("z [cm]"); ax.set_ylabel("T [C]")
        ax.set_title(label); ax.grid(True, alpha=0.3)
        for k, ls, lw, lab in curves:
            ln, = ax.plot([], [], ls, lw=lw, label=lab)
            lns[k] = ln
        ax.legend(fontsize=8)

    def animate(frame_idx):
        alpha = alphas[frame_idx]
        for res, lns in [(res_int, lines_int), (res_ext, lines_ext)]:
            z = res["z"] * 100
            for k, _, _, _ in curves:
                T_ss = res[k]
                T_now = T_in + alpha * (T_ss - T_in)
                lns[k].set_data(z, T_now)
        pct = int(alpha * 100)
        fig.suptitle(f"Associacao de canais — evolucao termica  [{pct}%]",
                     fontsize=12, fontweight="bold")
        return list(lines_int.values()) + list(lines_ext.values())

    ani = animation.FuncAnimation(fig, animate, frames=len(alphas),
                                  interval=1000//fps, blit=True)
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(outpath), writer="pillow", fps=fps)
    plt.close(fig)
    return f"  -> GIF associacao salvo em {outpath}  ({len(alphas)} frames, {fps} fps)"


def cmd_nucleo_init(state, args):
    """Inicializa matriz de mapa de nucleo vazia. Uso: nucleo_init [rows=5] [cols=5]"""
    kw   = parse_kwargs(args) if args else {}
    rows = int(fnum(kw.get("rows", "5")))
    cols = int(fnum(kw.get("cols", "5")))
    state.nucleo_data = {
        "rows":  rows,
        "cols":  cols,
        "vals":  [[0.0] * cols for _ in range(rows)],
        "tipos": [["EC"] * cols for _ in range(rows)],
    }
    return f"  nucleo {rows}x{cols} inicializado — use nucleo_linha para preencher."

def cmd_nucleo_linha(state, args):
    """Define uma linha do mapa. Uso: nucleo_linha row=N vals=v1,v2,... [tipos=EC,CR,...]"""
    if state.nucleo_data is None:
        raise ValueError("use nucleo_init antes de nucleo_linha.")
    kw   = parse_kwargs(args)
    row  = int(fnum(kw["row"]))
    vals = [fnum(v) for v in kw["vals"].split(",")]
    cols = state.nucleo_data["cols"]
    if len(vals) != cols:
        raise ValueError(f"esperados {cols} valores, recebidos {len(vals)}.")
    tipos = kw["tipos"].split(",") if "tipos" in kw else ["EC"] * cols
    if len(tipos) != cols:
        raise ValueError(f"esperados {cols} tipos, recebidos {len(tipos)}.")
    state.nucleo_data["vals"][row]  = vals
    state.nucleo_data["tipos"][row] = tipos
    vals_str = "  ".join(f"{v:.3f}" for v in vals)
    return f"  linha {row}: {vals_str}"

def cmd_nucleo_mapa(state, args):
    """Gera mapa de distribuicao de potencia radial do nucleo (5x5)."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    kw      = parse_kwargs(args) if args else {}
    outfile = kw.get("png", "saidas/mapa_nucleo.png")
    titulo  = kw.get("titulo", "Distribuicao de Potencia Radial do Nucleo")

    # Usa dados inseridos manualmente, se existirem; senao usa padrao
    if state.nucleo_data is not None:
        mapa  = state.nucleo_data["vals"]
        tipos = state.nucleo_data["tipos"]
    else:
        mapa  = _MAPA_PADRAO
        tipos = _TIPOS_PADRAO
    nrows, ncols = len(mapa), len(mapa[0])

    # Encontra maximo e sua posicao
    val_max = 0.0; pos_max = (0,0)
    for i in range(nrows):
        for j in range(ncols):
            if mapa[i][j] > val_max:
                val_max = mapa[i][j]; pos_max = (i, j)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_xlim(0, ncols); ax.set_ylim(0, nrows)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(titulo, fontsize=13, fontweight="bold", pad=14)

    for i in range(nrows):
        for j in range(ncols):
            val  = mapa[i][j]
            tipo = tipos[i][j]
            # linha 0 do mapa = topo da figura -> invertemos y
            row_plot = nrows - 1 - i

            if val == 0.0:
                face = "#c0c0c0"; txt_col = "#333333"  # barra de controle
            elif tipo == "CR":
                face = "#8ab4d4"; txt_col = "white"     # refletor/controle
            else:
                # gradiente salmon proporcional ao fator
                intensity = 0.55 + 0.45 * (val / val_max)
                r = int(255 * intensity)
                g = int(160 * (1 - 0.4*(val/val_max)))
                b = int(130 * (1 - 0.5*(val/val_max)))
                face = "#{:02x}{:02x}{:02x}".format(min(r,255), max(g,0), max(b,0))
                txt_col = "#222222"

            rect = Rectangle((j, row_plot), 1, 1,
                              facecolor=face, edgecolor="white", linewidth=2.5)
            ax.add_patch(rect)

            # Destaque no EC mais quente
            if (i, j) == pos_max:
                rect2 = Rectangle((j+0.04, row_plot+0.04), 0.92, 0.92,
                                   fill=False, edgecolor="#cc0000", linewidth=3.0)
                ax.add_patch(rect2)

            label = f"{val:.3f}" if val > 0 else "---"
            ax.text(j + 0.5, row_plot + 0.5, label,
                    ha="center", va="center", fontsize=13,
                    fontweight="bold", color=txt_col)

    # Legenda
    from matplotlib.patches import Patch
    legenda = [
        Patch(facecolor="#d97b5a", edgecolor="white", label="EC combustivel"),
        Patch(facecolor="#8ab4d4", edgecolor="white", label="Controle/Refletor"),
        Patch(facecolor="#c0c0c0", edgecolor="white", label="Barra inserida (q=0)"),
        Patch(facecolor="white",   edgecolor="#cc0000", linewidth=2, label="EC mais quente"),
    ]
    ax.legend(handles=legenda, loc="lower center",
              bbox_to_anchor=(0.5, -0.06), ncol=2, fontsize=9, framealpha=0.8)

    fig.tight_layout()
    outpath = Path(outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=150, bbox_inches="tight")
    plt.close(fig)
    pos_str = "(" + str(pos_max[0]+1) + "," + str(pos_max[1]+1) + ")"
    return "  -> mapa salvo em " + str(outpath) + "  | EC mais quente: posicao " + pos_str + " fator=" + str(val_max)

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
        "show_geom":       lambda: (
            state.geom.summary() if state.geom
            else (
                "GeometriaReator:\n"
                + "  Lx="     + str(round(state.reator_geom.Lx*1e3,2))    + " mm"
                + "  Ly="     + str(round(state.reator_geom.Ly*1e3,2))    + " mm"
                + "  df="     + str(round(state.reator_geom.df*1e3,3))    + " mm"
                + "  dcl="    + str(round(state.reator_geom.dcl*1e3,3))   + " mm"
                + "  dch="    + str(round(state.reator_geom.dch*1e3,3))   + " mm"
                + "  Lcanal=" + str(round(state.reator_geom.Lcanal*1e3,2))+ " mm"
                + "  Dh="     + str(round(state.reator_geom.Dh*1e3,4))   + " mm"
                + "  A_flow=" + str(round(state.reator_geom.area_flow*1e6,4))+ " mm2"
            ) if state.reator_geom
            else "sem geometria."
        ),
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
        "nucleo_init":     lambda: cmd_nucleo_init(state, args),
        "nucleo_linha":    lambda: cmd_nucleo_linha(state, args),
        "nucleo_mapa":     lambda: cmd_nucleo_mapa(state, args),
        "nucleo_solve":    lambda: cmd_nucleo_solve(state, args),
        "ec_associacao":   lambda: cmd_ec_associacao(state, args),
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

def rodar_arquivo(caminho):
    state = EstadoREPL()
    print("="*72+f"\nExecutando: {caminho}\n"+"="*72)
    with open(caminho) as f:
        linhas = f.readlines()
    for linha in linhas:
        raw = linha.rstrip()
        if not raw or raw.lstrip().startswith("#"):
            if raw: print(f"  {raw}")
            continue
        print(f"\n>>> {raw}")
        try:
            r, stop = executar(state, raw)
            if r: print(r)
            if stop: break
        except Exception as e:
            print(f"  ERRO: {e}")

def rodar_repl():
    state=EstadoREPL()
    print("="*72+"\nMetodo Nodal - interativo\nDigite help ou quit\n"+"="*72)
    while True:
        try:
            l=input(">>> ")
        except (EOFError, KeyboardInterrupt):
            print("\nSaindo."); break
        try:
            r, stop = executar(state, l)
            if r: print(r)
            if stop: break
        except Exception as e:
            print(f"  ERRO: {e}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "--demo":
            rodar_demo()
        else:
            rodar_arquivo(arg)
    else:
        rodar_repl()
