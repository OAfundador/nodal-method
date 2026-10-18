"""Exemplo Python: reutiliza um componente em 25 canais independentes."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from associacao import ModelAssembly
from visual.relatorio_associacao import export_report


def build_example():
    model = json.loads((Path(__file__).parent / 'componentes/canal_aquecido.json').read_text(encoding='utf-8'))
    group = ModelAssembly(name='Associação didática de 25 canais')
    group.add_component('canal', model)
    for index in range(1, 26):
        group.associate(f'canal_{index}', 'canal', source_factor=(index - 1) / 12, name=f'Canal {index}')
    group.set_view({'axis': 'y', 'region': 'aquecedor', 'unit': 'C'})
    return group


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saida', default='saidas/associacao', help='pasta dos arquivos gerados')
    args = parser.parse_args()
    destination = Path(args.saida)
    destination.mkdir(parents=True, exist_ok=True)
    group = build_example()
    result = group.solve(on_progress=lambda i, n, name: print(f'{i}/{n}: {name}'))
    group.export_csv(destination / 'temperaturas.csv')
    (destination / 'conjunto.json').write_text(json.dumps(group.to_dict(), ensure_ascii=False, indent=2), encoding='utf-8')
    (destination / 'conjunto.txt').write_text(group.to_repl(), encoding='utf-8')
    export_report(group, destination / 'relatorio.html')
    print(f'{result["count"]} canais; máxima {result["maximum_C"]:.6f} °C; resíduo {result["max_residual_W"]:.3e} W')
    print(f'Arquivos em: {destination.resolve()}')


if __name__ == '__main__':
    main()
