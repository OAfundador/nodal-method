from contextlib import redirect_stdout
from pathlib import Path
import csv
import io
import tempfile
import unittest

from exemplos.nodal_repl import CaseRunner, EstadoREPL, executar, main

ROOT = Path(__file__).resolve().parents[1]


class ReplTests(unittest.TestCase):
    def test_solve_options_are_used(self):
        for option in ("tol=0", "max_iter=0", "unknown=1"):
            runner = CaseRunner()
            with self.subTest(option=option), self.assertRaises(ValueError):
                runner.run_lines(["node a diffusion Q=1", "node b boundary T_fixed=0",
                                  "link a b cond G=1", "solve " + option])

    def test_invalid_direction_does_not_create_a_link(self):
        state = EstadoREPL()
        executar(state, "node a fluid")
        executar(state, "node b fluid")
        with self.assertRaisesRegex(ValueError, "dir"):
            executar(state, "link a b fluid G=1 dir=typo")
        self.assertEqual(state.net.links, [])
        with self.assertRaisesRegex(ValueError, "phase"):
            executar(state, "material agua phase=gas")

    def test_case_stops_on_error_and_reports_line(self):
        runner = CaseRunner()
        with self.assertRaisesRegex(ValueError, "caso:2"):
            runner.run_lines(["node a fluid", "nodo b fluid", "node c fluid"], source="caso")
        self.assertNotIn("c", runner.state.ids)

    def test_utf8_and_csv_units(self):
        (ROOT / "tmp").mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / "tmp") as directory:
            path = Path(directory) / "caso.txt"
            path.write_text("# condução\nnode água boundary T_fixed=30\n", encoding="utf-8-sig")
            runner = CaseRunner()
            runner.run_file(path)
            output = Path(directory) / "nos.csv"
            runner.export_csv(output)
            with output.open(encoding="utf-8-sig", newline="") as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row["nome"], "água")
            self.assertEqual(float(row["temperatura_C"]), 30)
            self.assertEqual(float(row["temperatura_K"]), 303.15)

    def test_unclosed_and_unmatched_loop_are_errors(self):
        for lines in (["for i=1..2", "node a_$i fluid"], ["end"]):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                CaseRunner().run_lines(lines)

    def test_loop_variable_is_restored_after_failure(self):
        runner = CaseRunner()
        with self.assertRaises(ValueError):
            runner.run_lines(["let i = 9", "for i=1..2", "invalid", "end"])
        self.assertEqual(runner.state.vars["i"], 9)

    def test_nested_loops_and_inline_comments(self):
        runner = CaseRunner()
        runner.run_lines(["for i=1..2", "for j=1..2", "node n_$(i)_$j fluid # comentário", "end", "end"])
        self.assertEqual(len(runner.state.net.nodes), 4)

    def test_original_examples_still_converge(self):
        for path in sorted((ROOT / "exemplos").glob("*.txt")):
            with self.subTest(path=path.name):
                lines = [line for line in path.read_text(encoding="utf-8").splitlines()
                         if line.strip().split(maxsplit=1)[:1] not in (["viz"], ["gif"])]
                runner = CaseRunner()
                runner.run_lines(lines, source=path.name)
                self.assertIsNotNone(runner.state.last_result)
                self.assertTrue(runner.state.last_result.success)

    def test_list_cli_does_not_enter_repl(self):
        with redirect_stdout(io.StringIO()) as stream:
            self.assertEqual(main(["--listar"]), 0)
        self.assertIn("interativo_chip.txt", stream.getvalue())
