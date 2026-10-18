"""Servidor de loopback: sem publicação, sem serviços externos."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import json
import mimetypes
import threading
from urllib.parse import urlsplit

from visual.modelo import VisualModel, demo_cases

ROOT = Path(__file__).resolve().parents[1]
STATIC = Path(__file__).with_name("static")
COMPUTE_LOCK = threading.Lock()


class NodalRequestHandler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

    def _reply(self, status, content, kind="application/json; charset=utf-8"):
        if not isinstance(content, bytes):
            content = json.dumps(content, ensure_ascii=False, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(content)

    def _local_request(self):
        port = self.server.server_port
        if self.headers.get("Host") not in (f"127.0.0.1:{port}", f"localhost:{port}"):
            self._reply(403, {"error": "Acesso permitido apenas pelo endereço local da ferramenta."})
            return False
        origin = self.headers.get("Origin")
        if origin and origin not in (f"http://127.0.0.1:{port}", f"http://localhost:{port}"):
            self._reply(403, {"error": "Origem não permitida."})
            return False
        return True

    def do_GET(self):
        if not self._local_request():
            return
        path = urlsplit(self.path).path
        if path == "/api/examples":
            self._reply(200, demo_cases())
            return
        files = {"/": STATIC / "index.html", "/app.js": STATIC / "app.js", "/style.css": STATIC / "style.css"}
        if path.startswith(("/benchmarks/", "/comparacao/")):
            name = path.split("/", 2)[2]
            if name in {"literatura.csv", "literatura.json", "benchmarks_online.json", "onur.html",
                        "onur_curvas.png", "onur_transiente.csv", "onur_permanente.csv", "onur_resumo.json"}:
                files[path] = ROOT / "saidas_comparacao" / name
        target = files.get(path)
        if target and target.is_file():
            mime = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
            self._reply(200, target.read_bytes(), mime + ("; charset=utf-8" if mime.startswith("text/") else ""))
        else:
            self._reply(404, {"error": "Arquivo não encontrado."})

    def do_POST(self):
        if not self._local_request():
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 1024 * 1024:
                raise ValueError("O modelo deve ter até 1 MB.")
            payload = json.loads(self.rfile.read(length))
            route = urlsplit(self.path).path
            if route not in ("/api/mesh", "/api/solve", "/api/literature", "/api/onur"):
                self._reply(404, {"error": "Operação não encontrada."})
                return
            if not COMPUTE_LOCK.acquire(blocking=False):
                self._reply(409, {"error": "Há um cálculo em andamento. Aguarde sua conclusão."})
                return
            try:
                if route == "/api/literature":
                    from validacao_livro import BookBenchmarks
                    from validacao_online import OnlineBenchmarks
                    from validacao_onur import combined_benchmark
                    response = {"results": [r.to_dict() for r in BookBenchmarks().run()],
                                "online": [r.to_dict() for r in OnlineBenchmarks().run()],
                                "onur": combined_benchmark()}
                elif route == "/api/onur":
                    from validacao_onur import reproduce_figures
                    destination = ROOT / "saidas_comparacao"
                    if not (destination / "onur.html").is_file():
                        reproduce_figures(destination)
                    response = {"url": "/benchmarks/onur.html"}
                else:
                    model = VisualModel(payload)
                    response = model.solve() if route == "/api/solve" else model.mesh_description()
            finally:
                COMPUTE_LOCK.release()
            self._reply(200, response)
        except (ValueError, TypeError, KeyError, RuntimeError, ImportError) as error:
            self._reply(400, {"error": str(error)})
        except Exception:
            self._reply(500, {"error": "O cálculo falhou. Verifique o modelo e o terminal do servidor."})
            import traceback
            traceback.print_exc()


def create_server(port=8765):
    return ThreadingHTTPServer(("127.0.0.1", port), NodalRequestHandler)
