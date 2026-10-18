"""Relatório HTML portátil da associação: abre localmente, sem servidor."""
import json
from pathlib import Path

STATIC = Path(__file__).with_name('static')


def export_report(assembly, path):
    if assembly.result is None:
        raise ValueError('Resolva a associação antes de salvar o relatório.')
    payload = json.dumps({'assembly': assembly.to_dict(), 'result': assembly.result}, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    html = (STATIC / 'assembly.html').read_text(encoding='utf-8')
    css = (STATIC / 'assembly.css').read_text(encoding='utf-8')
    js = (STATIC / 'assembly.js').read_text(encoding='utf-8')
    html = html.replace('<link rel="stylesheet" href="/assembly.css">', '<style>' + css + '</style>')
    html = html.replace('<script src="/assembly.js"></script>', '<script id="saved-assembly" type="application/json">' + payload + '</script><script>' + js + '</script>')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding='utf-8')
