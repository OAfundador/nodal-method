"""Abra o editor local com: python app_visual.py."""
import argparse
import webbrowser

from visual.servidor import create_server


def main(argv=None):
    parser = argparse.ArgumentParser(description="Editor visual local do Método Nodal.")
    parser.add_argument("--porta", type=int, default=8765)
    parser.add_argument("--sem-abrir", action="store_true")
    args = parser.parse_args(argv)
    try:
        server = create_server(args.porta)
    except OSError as error:
        parser.exit(1, f"Não foi possível abrir a porta {args.porta}: {error}. Use --porta OUTRA_PORTA.\n")
    address = f"http://127.0.0.1:{server.server_port}"
    print(f"Método Nodal visual: {address}\nMantenha este terminal aberto. Ctrl+C encerra a ferramenta.", flush=True)
    if not args.sem_abrir:
        webbrowser.open(address)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nFerramenta encerrada.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
