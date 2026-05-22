"""
nodal_repl.py

REPL / runner interativo do Metodo Nodal.

Runner 100% generico — nao conhece nenhum dominio fisico especifico.
Usa apenas os 5 modulos core: condutancias, geometria, materiais, nos, solver.
O que resolver e como resolver e definido exclusivamente pelo arquivo .txt.

Uso:
  python exemplos/nodal_repl.py                    # modo interativo
  python exemplos/nodal_repl.py --demo             # demo rapido
  python exemplos/nodal_repl.py exemplos/meu_caso.txt
"""

from __future__ import annotations
import csv as _csv
import math as _math
import shlex
import sys
from dataclasses import dataclass as _dataclass, field as _field
from pathlib import Path

import numpy as _np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from condutancias import conduction_G, convection_G
from geometria import Geometry2D
from materiais import Material, MaterialPhase
from nos import (
    NodalNetwork, NodeKind, TransferKind, LinkDirection,
    build_network_from_geometry, conduction_func_between_nodes,
)
from solver import solve_steady_state
















def _make_expr_func(expr_str):
    """Compila uma expressao em T (e opcionalmente P) num callable f(T,P=None)."""
    _ns = {"math": _math, **vars(_math)}
    return lambda T, P=None, _e=expr_str, _n=_ns: float(
        eval(_e, {"T": T, "P": P, **_n}))


def _make_t_sat_func(expr_str):
    """Compila uma expressao em P (Pa) num callable f(P_Pa) -> T_sat (C)."""
    _ns = {"math": _math, **vars(_math)}
    return lambda P, _e=expr_str, _n=_ns: float(
        eval(_e, {"P": P, **_n}))


def _material_from_props(nome, props):
    """Cria Material a partir do dict de propriedades armazenado no REPL.

    Para cada propriedade (k, rho, cp, mu), aceita:
      - <prop>_expr : string de expressao em T (e opcionalmente P) -> *_func
      - <prop>_poly : lista de coeficientes [a0, a1, ...]           -> *_coeffs
      - <prop>      : constante                                      -> escalar

    Se "t_sat_expr" estiver presente (expressao em P em Pa), e atachado
    como mat.t_sat_func — usado em extrair_resultados para a margem ate
    a saturacao.
    """
    ph = (MaterialPhase.FLUID if props.get("phase", "solid") == "fluid"
          else MaterialPhase.SOLID)

    kwargs = {"name": nome, "phase": ph}
    for prop in ("k", "rho", "cp", "mu"):
        expr = props.get(prop + "_expr")
        poly = props.get(prop + "_poly")
        const = props.get(prop)
        if expr:
            kwargs[prop + "_func"] = _make_expr_func(expr)
        elif poly:
            kwargs[prop + "_coeffs"] = poly
        elif const is not None:
            kwargs[prop] = const

    mat = Material(**kwargs)

    t_sat_expr = props.get("t_sat_expr")
    if t_sat_expr:
        mat.t_sat_func = _make_t_sat_func(t_sat_expr)
    else:
        mat.t_sat_func = None

    return mat


# ===========================================================================
# REPL
# ===========================================================================

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
         material u3si2  k_expr="1.73073*(3978.1/(724.61+1.8*T)+6.02366e-12*(1.8*T+492)**3)" rho=4300 cp=836
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

--- VARIAVEIS E LOOPS ---
  let <nome> = <expressao>      Define variavel (escalar ou lista).
                                Pode usar math.*, materiais, helpers, vars.
  for <var>=<lo>..<hi>          Inicia loop. Termina com 'end'.
    ...                         Linhas podem usar $var e $(expr).
  end
  reset_net                     Zera so a rede (preserva materiais/vars/geom).

--- EXPRESSOES (G_expr, Q_expr, ...) ---
  Em qualquer comando voce pode usar X_expr="..." em vez de X=v.
  Disponivel no namespace: math.*, materiais (ex.: agua.prop('k', T, P)),
  conduction_G(k,A,L), convection_G(h,A),
  dittus_h(material, T, P, mdot, A_flow, Dh),
  sum_over("var", lo, hi, "expr"),
  e em G_expr: T['nome_no'] -> temperatura atual durante a iteracao.

--- VISUALIZACAO ---
  viz  [png=arquivo.png] [titulo=Texto]
  gif  [gif=arquivo.gif] [fps=10] [steps=40] [dpi=100] [titulo=Texto]

--- CONTROLE ---
  help / quit / exit
============================================================
"""


# ===========================================================================
# Primitivas genericas: let / for / $var / $(expr) / G_expr / helpers de eval
# ===========================================================================

import re as _re_prim

def _dittus_h(material, T, P, mdot, A_flow, Dh):
    """h Dittus-Boelter Nu=0.023 Re^0.8 Pr^0.4. Retorna [W/(m2 K)]."""
    rho = material.prop("rho", T, P); mu = material.prop("mu", T, P)
    cp  = material.prop("cp",  T, P); k  = material.prop("k",  T, P)
    v   = mdot / (rho * A_flow)
    Re  = rho * v * Dh / mu; Pr = cp * mu / k
    Nu  = 0.023 * Re**0.8 * Pr**0.4
    return Nu * k / Dh

def _expr_namespace(state, extra=None):
    """Namespace para eval (let / Q_expr / G_expr / interpolacao)."""
    ns = {"math": _math, **vars(_math)}
    ns["dittus_h"]     = _dittus_h
    ns["conduction_G"] = conduction_G
    ns["convection_G"] = convection_G
    for nm, props in _MATERIAIS_BUILTIN.items():
        ns[nm] = _material_from_props(nm, props)
    for nm, props in state.materials.items():
        ns[nm] = _material_from_props(nm, props)
    ns.update(state.vars)
    if extra: ns.update(extra)
    def _sum_over(var, lo, hi, expr_str, _ns=ns):
        total = 0.0
        for v in range(int(lo), int(hi) + 1):
            ns2 = dict(_ns); ns2[var] = v
            total += float(eval(expr_str, ns2))
        return total
    ns["sum_over"] = _sum_over
    return ns

_VAR_RE_PRIM = _re_prim.compile(r'\$([A-Za-z_]\w*)')

def _fmt_val(v):
    if isinstance(v, bool):  return str(int(v))
    if isinstance(v, int):   return str(v)
    if isinstance(v, float):
        if v.is_integer() and abs(v) < 1e15: return str(int(v))
        return repr(v)
    return str(v)

def _interpolate(line, state):
    out = []; i = 0; n = len(line)
    while i < n:
        c = line[i]
        if c == '$' and i+1 < n and line[i+1] == '(':
            depth = 0; j = i+1
            while j < n:
                if line[j] == '(': depth += 1
                elif line[j] == ')':
                    depth -= 1
                    if depth == 0: break
                j += 1
            if depth != 0:
                raise ValueError(f"$(...) nao balanceado em: {line!r}")
            out.append(_fmt_val(eval(line[i+2:j], _expr_namespace(state))))
            i = j + 1
        elif c == '$':
            m = _VAR_RE_PRIM.match(line, i)
            if m:
                nm = m.group(1)
                if nm not in state.vars:
                    raise ValueError(f"variavel ${nm} indefinida")
                out.append(_fmt_val(state.vars[nm])); i = m.end()
            else:
                out.append('$'); i += 1
        else:
            out.append(c); i += 1
    return "".join(out)

def cmd_let(state, raw):
    """let <nome> = <expr> (escalar ou lista)."""
    if "=" not in raw: raise ValueError("uso: let <nome> = <expr>")
    name, expr = raw.split("=", 1); name = name.strip()
    val = eval(expr.strip(), _expr_namespace(state))
    state.vars[name] = val
    if isinstance(val, (list, tuple)):
        return f"  -> {name} = [{len(val)} elementos]"
    return f"  -> {name} = {_fmt_val(val)}"

def cmd_reset_net(state):
    """Zera apenas a rede (nos, links, ids). Preserva materiais/vars/geom."""
    state.net = NodalNetwork(); state.ids = {}; state.x_counter = 0.0
    return "  rede zerada (materiais/vars/geom preservados)."

def _parse_for_header(line):
    rest = line.strip()[4:].strip()
    if "=" not in rest or ".." not in rest:
        raise ValueError("uso: for <var>=<lo>..<hi>")
    np_, rp = rest.split("=", 1); lo_s, hi_s = rp.split("..", 1)
    return {"var": np_.strip(), "lo": lo_s.strip(), "hi": hi_s.strip(),
            "lines": [], "depth": 0}


class EstadoREPL:
    def __init__(self):
        self.net = NodalNetwork()
        self.ids = {}
        self.x_counter = 0.0
        self.geom = None
        self.materials = {}
        self.vars = {}          # let var = expr
        self.for_stack = []     # bloco for ... end

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
    ns = _expr_namespace(state)
    def _val(key, default):
        ek = key + "_expr"
        if ek in kw:  return float(eval(kw[ek], ns))
        if key in kw: return fnum(kw[key])
        return default
    if kind is NodeKind.BOUNDARY:
        Tf = _val("t_fixed", None)
        if Tf is None: raise ValueError("node boundary precisa de T_fixed=...")
        nid = state.net.add_node(name=nome, kind=kind, x=state.x_counter, y=0.0,
                                  fixed_temperature=Tf, temperature=Tf)
    else:
        V_def = 0.0 if kind is NodeKind.ARITHMETIC else 1.0
        nid = state.net.add_node(name=nome, kind=kind, x=state.x_counter, y=0.0,
                                  volume=_val("v", V_def),
                                  source=_val("q", 0.0),
                                  temperature=_val("t_init", 30.0))
    state.ids[nome] = nid; state.x_counter += 1.0
    return f"  -> no {nome} criado (id={nid}, kind={kind.value})"

def cmd_link(state, args):
    a, b, tipo = args[0], args[1], args[2].lower()
    kw = parse_kwargs(args[3:])
    dir_map = {
        "fwd": LinkDirection.I_TO_J, "i_to_j": LinkDirection.I_TO_J,
        "bwd": LinkDirection.J_TO_I, "j_to_i": LinkDirection.J_TO_I,
        "undirected": LinkDirection.UNDIRECTED,
    }
    dir_kw    = kw.get("dir", "undirected").lower()
    direction = dir_map.get(dir_kw, LinkDirection.UNDIRECTED)
    if "g_expr" in kw:
        expr_str       = kw["g_expr"]
        materials_objs = {nm: _material_from_props(nm, p)
                          for nm, p in _MATERIAIS_BUILTIN.items()}
        for nm, p in state.materials.items():
            materials_objs[nm] = _material_from_props(nm, p)
        names_to_ids = dict(state.ids)
        base_ns = {"math": _math, **vars(_math)}
        base_ns["dittus_h"]     = _dittus_h
        base_ns["conduction_G"] = conduction_G
        base_ns["convection_G"] = convection_G
        base_ns.update(materials_objs); base_ns.update(state.vars)
        def gf(T_map, _e=expr_str, _ns=base_ns, _ids=names_to_ids):
            T_named = {nm: T_map[nid] for nm, nid in _ids.items()
                       if nid in T_map}
            local = dict(_ns); local["T"] = T_named
            return float(eval(_e, local))
        state.net.add_link(state.ids[a], state.ids[b], LINK_KIND_MAP[tipo],
                           conductance_func=gf, direction=direction,
                           name=a + "_" + b)
        return f"  -> ligacao {a}-{b} ({tipo}) G_expr  dir={dir_kw}"
    G_val = float(eval(kw["g"], _expr_namespace(state)))
    state.net.add_link(state.ids[a], state.ids[b], LINK_KIND_MAP[tipo],
                       conductance=G_val, direction=direction,
                       name=a + "_" + b)
    return "  -> ligacao " + a + "-" + b + f" ({tipo}) G={G_val:.6g} W/K  dir=" + dir_kw

def cmd_show(state): return state.net.summary(max_nodes=20, max_links=30)

def cmd_solve(state):
    if len(state.net.nodes) == 0: return "  rede vazia."
    r = solve_steady_state(state.net, tol=1e-9, max_iter=500)
    out = [f"  convergiu={r.success}, |R|={r.residual_norm:.2e}, iter={r.iterations}"]
    if state.ids:
        out += ["", "  Temperaturas:"] + [
            f"    {n:8s} = {state.net.nodes[i].temperature:9.4f} C"
            for n, i in state.ids.items()
        ]
    return "\n".join(out)

def _req(state):
    if not state.geom: raise ValueError("Use domain W= H= primeiro.")
    return state.geom


# ---------------------------------------------------------------------------
# Materiais
# ---------------------------------------------------------------------------
# Builtin "agua" carrega as correlacoes de Incropera (Tab. A.6) e a equacao
# de Antoine como expressoes em T (C) e P (Pa). Quando o usuario faz
#   material agua phase=fluid k=0.620 rho=994 cp=4182 mu=7.7e-4
# cmd_material derruba os respectivos *_expr e usa as constantes do usuario,
# mas mantem t_sat_expr (saturacao continua disponivel para a margem termica).
_MATERIAIS_BUILTIN = {
    "cobre":    dict(phase="solid",  k=380.0,   rho=8900.0,  cp=385.0),
    "aluminio": dict(phase="solid",  k=205.0,   rho=2700.0,  cp=900.0),
    "aco":      dict(phase="solid",  k=50.0,    rho=7800.0,  cp=500.0),
    "silicio":  dict(phase="solid",  k=150.0,   rho=2330.0,  cp=712.0),
    "uo2":      dict(phase="solid",  k=3.6,     rho=10970.0, cp=247.0),
    "agua":     dict(
        phase="fluid",
        # Incropera Tab. A.6 (T em C, k em W/(m K))
        k_expr="0.5615 + 1.939e-3*T - 7.51e-6*T**2",
        # Polinomio de Kell para densidade da agua liquida (T em C, rho em kg/m3)
        rho_expr=("999.842594 + 6.793952e-2*T - 9.095290e-3*T**2 "
                  "+ 1.001685e-4*T**3 - 1.120083e-6*T**4 + 6.536332e-9*T**5"),
        # Calor especifico (T em C, cp em J/(kg K))
        cp_expr="4217.0 - 3.358*T + 0.04148*T**2 - 1.6e-4*T**3",
        # Viscosidade dinamica (Reichardt) (T em C, mu em Pa.s)
        mu_expr="2.414e-5 * 10.0**(247.8 / (T + 273.15 - 140.0))",
        # Saturacao por Antoine (P em Pa, T_sat em C)
        t_sat_expr="1730.63 / (8.07131 - log10(P/133.322)) - 233.426",
    ),
    "ar":       dict(phase="fluid",  k=0.026,   rho=1.18,    cp=1007.0, mu=1.85e-5),
    # Combustivel U3Si2-Al (PDF TNR5703). k(T)[BTU/h.ft.F] convertido pra SI:
    "u3si2_al": dict(phase="solid",
                     k_expr=("1.73073*(3978.1/(724.61 + 1.8*T) "
                             "+ 6.02366e-12*(1.8*T + 492)**3)"),
                     rho=4300.0, cp=836.0),
}

def cmd_material(state, args):
    if not args:
        raise ValueError("material precisa de nome. Ex: material cobre k=380")
    nome = args[0].lower()
    kw   = parse_kwargs(args[1:]) if len(args) > 1 else {}
    # Merge sobre o builtin: usuario sobrepoe campos individuais.
    props = dict(_MATERIAIS_BUILTIN.get(nome, {}))

    # Constantes explicitas derrubam expr/poly herdados para a mesma propriedade.
    for campo in ("k", "rho", "cp", "mu"):
        if campo in kw:
            props[campo] = fnum(kw[campo])
            props.pop(campo + "_expr", None)
            props.pop(campo + "_poly", None)

    # Polinomios derrubam constante e expr da mesma propriedade.
    for campo in ("k_poly", "rho_poly", "cp_poly", "mu_poly"):
        if campo in kw:
            props[campo] = [fnum(v) for v in kw[campo].split(",")]
            base = campo[:-5]
            props.pop(base, None)
            props.pop(base + "_expr", None)

    # Expressoes derrubam constante e poly da mesma propriedade.
    for campo in ("k_expr", "rho_expr", "cp_expr", "mu_expr"):
        if campo in kw:
            props[campo] = kw[campo]
            base = campo[:-5]
            props.pop(base, None)
            props.pop(base + "_poly", None)

    # Saturacao (fluidos): expressao em P (Pa).
    if "t_sat_expr" in kw:
        props["t_sat_expr"] = kw["t_sat_expr"]

    if "phase" in kw:
        props["phase"] = kw["phase"].lower()
    if not props.get("phase"):
        props["phase"] = "solid"

    state.materials[nome] = props

    partes = []
    if "k_expr" in props:
        partes.append("k=f(T)[" + props["k_expr"][:40] + "]")
    elif "k_poly" in props:
        partes.append("k=poly" + str(props["k_poly"]))
    elif "k" in props:
        partes.append(f"k={props['k']} W/(mK)")
    for campo, unid in [("rho", "kg/m3"), ("cp", "J/(kgK)"), ("mu", "Pa.s")]:
        if campo + "_expr" in props:
            partes.append(campo + "=f(T)")
        elif campo + "_poly" in props:
            partes.append(campo + "=poly")
        elif campo in props:
            partes.append(campo + "=" + str(props[campo]) + " " + unid)
    if "t_sat_expr" in props:
        partes.append("T_sat=f(P)")
    return "  -> material " + repr(nome) + " [" + props["phase"] + "]  " + "  ".join(partes)

def _resumo_propriedade_k(d):
    if "k_expr" in d:
        return "k=f(T)"
    if "k_poly" in d:
        return "k=poly"
    if "k" in d:
        return "k=" + str(d["k"])
    return "k=-"

def cmd_materiais(state, args):
    linhas = ["  --- Builtin ---"]
    for n, d in _MATERIAIS_BUILTIN.items():
        linhas.append("    " + n.ljust(10) + " [" + d["phase"] + "]  "
                      + _resumo_propriedade_k(d))
    if state.materials:
        linhas.append("  --- Definidos na sessao ---")
        for n, d in state.materials.items():
            linhas.append("    " + n.ljust(10) + " [" + d.get("phase", "solid") + "]  "
                          + _resumo_propriedade_k(d))
    return "\n".join(linhas)

def cmd_domain(state, args):
    kw = parse_kwargs(args)
    state.geom = Geometry2D(width=fnum(kw["w"]), height=fnum(kw["h"]))
    return f"  -> Geometry2D {kw['w']} x {kw['h']} m"

def cmd_region(state, args):
    geom = _req(state); kw = parse_kwargs(args[1:]); nome = args[0]
    if "mat" in kw:
        mat_nome = kw["mat"].lower()
        # cmd_material ja merge builtin+usuario; preferir o dict da sessao
        # quando existe (preserva a derrubada de *_expr feita por overrides).
        if mat_nome in state.materials:
            props = dict(state.materials[mat_nome])
        else:
            props = dict(_MATERIAIS_BUILTIN.get(mat_nome, {}))
        if not props:
            raise ValueError("material " + repr(mat_nome) + " nao definido.")
        mat = _material_from_props(mat_nome, props)
    else:
        k_val = fnum(kw["k"])
        mat   = Material(name="mat_" + nome, phase=MaterialPhase.SOLID, k=k_val)
    geom.material(nome, mat,
                  x0=fnum(kw["x0"]), x1=fnum(kw["x1"]),
                  y0=fnum(kw["y0"]), y1=fnum(kw["y1"]))
    try:
        k_show = mat.prop("k", 30.0)
        k_str  = "{:.6g}".format(k_show)
    except Exception:
        k_str = "f(T)" if mat.k_func is not None else "?"
    return "  -> regiao " + repr(nome) + " k=" + k_str

def cmd_source(state, args):
    geom = _req(state); kw = parse_kwargs(args[1:])
    geom.source(region=args[0], kind="volumetric", value=fnum(kw["q"]))
    return f"  -> fonte {kw['q']} W/m3 em {args[0]!r}"

def cmd_bc(state, args):
    geom = _req(state); kw = parse_kwargs(args[2:], lower_keys=False)
    geom.bc(args[0].lower(), args[1].lower(), **{k: fnum(v) for k, v in kw.items()})
    return f"  -> BC {args[1]!r} em {args[0]!r}"

def cmd_mesh(state, args):
    geom = _req(state); kw = parse_kwargs(args)
    geom.mesh(nx=int(fnum(kw["nx"])), ny=int(fnum(kw["ny"])), fix=kw.get("fix", "auto"))
    return f"  -> malha nx={geom.nx} ny={geom.ny}"

def cmd_build_from_geom(state, args):
    geom = _req(state); kw = parse_kwargs(args) if args else {}
    state.net = build_network_from_geometry(geom, thickness_z=fnum(kw.get("tz", "1")))
    state.ids.clear()
    return f"  -> NodalNetwork: {len(state.net.nodes)} nos, {len(state.net.links)} ligacoes"


# ---------------------------------------------------------------------------
# Comandos do reator tipo placa
# ---------------------------------------------------------------------------













# ===========================================================================
# Transporte entalpico generico
# ===========================================================================

def cmd_fluid_chain(state, args):
    """fluid_chain n1 n2 ... mdot=v cp=v [T_in=v]"""
    nos_names = [a for a in args if "=" not in a]
    kw        = parse_kwargs([a for a in args if "=" in a])
    if len(nos_names) < 2:
        raise ValueError("fluid_chain precisa de pelo menos 2 nos")
    if "mdot" not in kw or "cp" not in kw:
        raise ValueError("fluid_chain precisa de mdot=v e cp=v")
    mdot = fnum(kw["mdot"]); cp = fnum(kw["cp"]); G = mdot * cp
    T_in_val = fnum(kw["t_in"]) if "t_in" in kw else None

    def resolve(name):
        if name in state.ids: return state.ids[name]
        try:
            nid = int(name)
            if nid in state.net.nodes: return nid
        except ValueError:
            pass
        raise ValueError("No nao encontrado: " + repr(name))

    chain = [resolve(n) for n in nos_names]; msgs = []
    if T_in_val is not None:
        fn = state.net.nodes[chain[0]]
        inlet_id = state.net.add_node(
            name="inlet_" + nos_names[0], kind=NodeKind.BOUNDARY,
            x=fn.x - 1.0, y=fn.y,
            fixed_temperature=T_in_val, temperature=T_in_val)
        state.net.add_link(inlet_id, chain[0], TransferKind.FLUID_TRANSPORT,
                           direction=LinkDirection.I_TO_J, conductance=G,
                           name="fluid_inlet_" + nos_names[0])
        msgs.append("  -> inlet BOUNDARY T=" + str(T_in_val) + " C -> " + nos_names[0])

    for k in range(len(chain) - 1):
        na, nb = chain[k], chain[k+1]
        state.net.add_link(na, nb, TransferKind.FLUID_TRANSPORT,
                           direction=LinkDirection.I_TO_J, conductance=G,
                           name="fluid_" + nos_names[k] + "_" + nos_names[k+1])
        msgs.append("  -> " + nos_names[k] + " -> " + nos_names[k+1]
                    + "  G=" + str(round(G, 4)) + " W/K")
    return "\n".join(msgs)


# ===========================================================================
# Visualizacao generica
# ===========================================================================

def _render_network(net, title="Rede Nodal", fig=None, ax=None, alpha_links=0.55):
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    ids   = list(net.nodes.keys())
    xs    = _np.array([net.nodes[i].x for i in ids])
    ys    = _np.array([net.nodes[i].y for i in ids])
    Ts    = _np.array([net.nodes[i].temperature for i in ids])
    kinds = [net.nodes[i].kind for i in ids]

    T_lo = float(Ts.min()); T_hi = float(Ts.max())
    if abs(T_hi - T_lo) < 0.01: T_hi = T_lo + 1.0
    norm = Normalize(vmin=T_lo, vmax=T_hi); cmap = plt.cm.hot

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
        xi, yi = id_pos[lk.node_i]; xj, yj = id_pos[lk.node_j]
        col, ls, lw = link_style.get(lk.kind, ("#888888", "-", 0.9))
        ax.plot([xi, xj], [yi, yj], color=col, ls=ls, lw=lw,
                alpha=alpha_links, zorder=1)
        if lk.direction != LinkDirection.UNDIRECTED:
            mx = 0.5*(xi+xj); my = 0.5*(yi+yj)
            dx = xj-xi; dy = yj-yi
            if lk.direction == LinkDirection.J_TO_I: dx=-dx; dy=-dy
            nd = (dx**2+dy**2)**0.5
            if nd > 1e-9:
                ax.annotate("",
                    xy=(mx+dx/nd*0.12, my+dy/nd*0.12), xytext=(mx, my),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.3), zorder=2)

    mk_map = {"diffusion":("o",90),"fluid":("s",90),"arithmetic":("^",55),"boundary":("D",65)}
    for kv in set(k.value for k in kinds):
        mask = [k.value == kv for k in kinds]; idx = [i for i,m in enumerate(mask) if m]
        xm = xs[idx]; ym = ys[idx]; Tm = Ts[idx]
        mk, sz = mk_map.get(kv, ("o", 80))
        ax.scatter(xm, ym, c=Tm, cmap=cmap, norm=norm,
                   marker=mk, s=sz, edgecolors="black", lw=0.6, zorder=4, label=kv)

    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    plt.colorbar(sm, ax=ax, label="T [degC]", shrink=0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(title + "  |  " + str(len(ids)) + " nos  " + str(len(net.links)) + " links")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.grid(True, alpha=0.18)
    return fig, ax, norm, cmap


def cmd_viz(state, args):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if len(state.net.nodes) == 0: return "  rede vazia."
    kw = parse_kwargs(args) if args else {}
    outfile = kw.get("png", "rede_nodal.png"); titulo = kw.get("titulo", "Rede Nodal")
    fig, ax, _, _ = _render_network(state.net, title=titulo)
    outpath = Path(outfile); outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=130, bbox_inches="tight"); plt.close(fig)
    return "  -> viz salvo em " + str(outpath)


def cmd_gif_generico(state, args):
    """GIF pseudo-transiente: frio -> regime permanente."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    if len(state.net.nodes) == 0: return "  rede vazia."
    kw      = parse_kwargs(args) if args else {}
    outfile = kw.get("gif", "rede_nodal.gif")
    fps     = int(fnum(kw.get("fps", "10")))
    dpi     = int(fnum(kw.get("dpi", "100")))
    n_steps = int(fnum(kw.get("steps", "40")))
    titulo  = kw.get("titulo", "Rede Nodal")

    net = state.net; ids = list(net.nodes.keys())
    xs  = _np.array([net.nodes[i].x for i in ids])
    ys  = _np.array([net.nodes[i].y for i in ids])
    T_final = _np.array([net.nodes[i].temperature for i in ids])
    kinds   = [net.nodes[i].kind for i in ids]
    T_lo = float(T_final.min()); T_hi = float(T_final.max())
    if abs(T_hi - T_lo) < 0.01: T_hi = T_lo + 1.0
    norm = Normalize(vmin=T_lo, vmax=T_hi); cmap = plt.cm.hot

    is_boundary = _np.array([net.nodes[i].is_boundary() for i in ids])
    T_start = _np.where(is_boundary, T_final, T_lo)
    alphas  = _np.concatenate([_np.linspace(0.0, 1.0, n_steps), _np.ones(6)])

    fig, ax = plt.subplots(figsize=(11, 6))
    link_style = {
        TransferKind.CONDUCTION:      ("#888888", "-",  0.8),
        TransferKind.CONVECTION:      ("#4499ff", "--", 0.9),
        TransferKind.FLUID_TRANSPORT: ("#ff4444", "-",  1.3),
        TransferKind.RADIATION:       ("#cc44ff", ":",  0.9),
        TransferKind.EQUIVALENT:      ("#aaaaaa", "-.", 0.7),
    }
    id_pos = {i: (net.nodes[i].x, net.nodes[i].y) for i in ids}
    for lk in net.links:
        xi, yi = id_pos[lk.node_i]; xj, yj = id_pos[lk.node_j]
        col, ls, lw = link_style.get(lk.kind, ("#888888", "-", 0.8))
        ax.plot([xi, xj], [yi, yj], color=col, ls=ls, lw=lw, alpha=0.45, zorder=1)
        if lk.direction != LinkDirection.UNDIRECTED:
            mx=0.5*(xi+xj); my=0.5*(yi+yj); dx=xj-xi; dy=yj-yi
            if lk.direction == LinkDirection.J_TO_I: dx=-dx; dy=-dy
            nd=(dx**2+dy**2)**0.5
            if nd > 1e-9:
                ax.annotate("", xy=(mx+dx/nd*0.12, my+dy/nd*0.12), xytext=(mx,my),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.2), zorder=2)

    mk_map = {"diffusion":"o","fluid":"s","arithmetic":"^","boundary":"D"}
    sz_map = {"diffusion":90,"fluid":90,"arithmetic":55,"boundary":65}
    scs    = {}
    for kv in set(k.value for k in kinds):
        mask = _np.array([k.value == kv for k in kinds])
        scs[kv] = ax.scatter(xs[mask], ys[mask], c=T_start[mask], cmap=cmap, norm=norm,
                              marker=mk_map.get(kv,"o"), s=sz_map.get(kv,80),
                              edgecolors="black", lw=0.6, zorder=4, label=kv)

    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    plt.colorbar(sm, ax=ax, label="T [degC]", shrink=0.85)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title(titulo)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.grid(True, alpha=0.18)
    pct_txt = ax.text(0.02, 0.97, "0%", transform=ax.transAxes, fontsize=9, va="top",
                      bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    def update(fi):
        a = float(alphas[fi]); T_cur = T_start + a*(T_final - T_start)
        for kv in scs:
            mask = _np.array([k.value == kv for k in kinds])
            scs[kv].set_array(T_cur[mask])
        pct_txt.set_text(str(int(round(a*100))) + "%   Tmax=" + f"{float(T_cur.max()):.1f}" + " C")
        return list(scs.values()) + [pct_txt]

    ani = animation.FuncAnimation(fig, update, frames=len(alphas),
                                  interval=1000//fps, blit=True)
    outpath = Path(outfile); outpath.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(outpath), writer="pillow", dpi=dpi); plt.close(fig)
    return "  -> GIF salvo em " + str(outpath) + "  (" + str(len(alphas)) + " frames, " + str(fps) + " fps)"


# ---------------------------------------------------------------------------
# Mapa de potencia do nucleo
# ---------------------------------------------------------------------------




















# ---------------------------------------------------------------------------
# Dispatcher principal
# ---------------------------------------------------------------------------

def executar(state, linha):
    raw = linha.strip()
    if not raw or raw.startswith("#"): return "", False

    # Bloco for ... end: bufferiza, rastreia depth, executa no end
    if state.for_stack:
        low = raw.lower(); top = state.for_stack[-1]
        if low == "end":
            if top.get("depth", 0) > 0:
                top["depth"] -= 1; top["lines"].append(raw)
                return "", False
            frame = state.for_stack.pop(); outs = []
            ns = _expr_namespace(state)
            lo = int(eval(frame["lo"], ns)); hi = int(eval(frame["hi"], ns))
            saved = state.vars.get(frame["var"])
            for v in range(lo, hi + 1):
                state.vars[frame["var"]] = v
                for bline in frame["lines"]:
                    r, stop = executar(state, bline)
                    if r: outs.append(r)
                    if stop: break
            if saved is None: state.vars.pop(frame["var"], None)
            else: state.vars[frame["var"]] = saved
            return "\n".join(outs), False
        if low.startswith("for "):
            top["depth"] = top.get("depth", 0) + 1
            top["lines"].append(raw); return "", False
        top["lines"].append(raw); return "", False
    if raw.lower().startswith("for "):
        state.for_stack.append(_parse_for_header(raw)); return "", False

    # Interpolacao $var e $(expr)
    raw = _interpolate(raw, state)

    # let antes do shlex (expr pode ter '=')
    if raw.lower().startswith("let "):
        return cmd_let(state, raw[4:]), False

    tokens = shlex.split(raw); cmd = tokens[0].lower(); args = tokens[1:]
    if cmd in ("quit","exit","q"): return "Saindo.", True
    cmds = {
        "help":            lambda: HELP,
        "reset":           lambda: state.reset() or "  rede zerada.",
        "reset_net":       lambda: cmd_reset_net(state),
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
        "show_geom":       lambda: (state.geom.summary() if state.geom
                                    else "sem geometria."),
        "build_from_geom": lambda: cmd_build_from_geom(state, args),
        "fluid_chain":     lambda: cmd_fluid_chain(state, args),
        "viz":             lambda: cmd_viz(state, args),
        "gif":             lambda: cmd_gif_generico(state, args),
    }
    if cmd not in cmds:
        raise ValueError(f"comando desconhecido: {cmd!r}. Digite help.")
    return cmds[cmd](), False


def cmd_gcond(args):
    kw=parse_kwargs(args); G=conduction_G(fnum(kw["k"]),fnum(kw["a"]),fnum(kw["l"]))
    return f"  G_cond = {G:.6g} W/K"

def cmd_gconv(args):
    kw=parse_kwargs(args); G=convection_G(fnum(kw["h"]),fnum(kw["a"]))
    return f"  G_conv = {G:.6g} W/K"


def rodar_demo():
    state = EstadoREPL()
    demo = [
        "# Exemplo chip - item (a)",
        "g_cond k=380 A=7.854e-7 L=0.02",
        "node Tc diffusion Q=4 T_init=30",
        "node Tp boundary T_fixed=44",
        "link Tc Tp cond G=0.17908",
        "solve",
        "# item (b): conveccao",
        "g_conv h=30 A=9e-4",
        "node Tar boundary T_fixed=20",
        "link Tc Tar conv G=0.027",
        "solve",
    ]
    print("="*72, "\nExemplo chip - demo\n" + "="*72)
    for l in demo:
        if l.strip(): print(f"\n>>> {l}")
        try:
            r, _ = executar(state, l)
            if r: print(r)
        except Exception as e:
            print(f"  ERRO: {e}")


def rodar_arquivo(caminho):
    state = EstadoREPL()
    print("="*72 + f"\nExecutando: {caminho}\n" + "="*72)
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
    state = EstadoREPL()
    print("="*72 + "\nMetodo Nodal - interativo\nDigite help ou quit\n" + "="*72)
    while True:
        try:
            l = input(">>> ")
        except (EOFError, KeyboardInterrupt):
            print("\nSaindo."); break
        try:
            r, stop = executar(state, l)
            if r: print(r)
            if stop: break
        except Exception as e:
            print(f"  ERRO: {e}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "--demo":
            rodar_demo()
        else:
            rodar_arquivo(arg)
    else:
        rodar_repl()
