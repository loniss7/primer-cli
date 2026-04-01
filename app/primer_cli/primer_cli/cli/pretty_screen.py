from __future__ import annotations

import argparse
import shlex
import sys
import termios
import tty
from dataclasses import dataclass
from typing import Callable

from primer_cli.core.exceptions import PrimerCliError

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.prompt import Confirm, Prompt
    from rich.table import Table
except Exception:
    Console = None
    Panel = None
    Prompt = None
    Confirm = None
    Table = None


Handler = Callable[[argparse.Namespace], int]
Runner = Callable[[Handler, argparse.Namespace], int]


def _ask(console: Console, prompt: str, default: str | None = None) -> str:
    return Prompt.ask(prompt, default=default, console=console).strip()


def _ask_required(console: Console, prompt: str) -> str:
    while True:
        value = _ask(console, prompt)
        if value:
            return value
        console.print("[bold red]Value is required.[/bold red]")


def _ask_extra_args(console: Console) -> list[str]:
    raw = _ask(console, "Extra args (exactly as in CLI, optional)", "")
    if not raw:
        return []
    try:
        return shlex.split(raw)
    except ValueError as e:
        console.print(f"[bold red]Cannot parse extra args:[/bold red] {e}")
        return []


def _build_fetch_argv(console: Console) -> list[str]:
    gene = _ask_required(console, "Gene")
    output = _ask_required(console, "Output FASTA path")
    argv = ["fetch", "--gene", gene, "--output", output]
    max_records = _ask(console, "Max sequences", "100")
    if max_records:
        argv += ["--max", max_records]
    email = _ask(console, "NCBI email (optional)", "")
    if email:
        argv += ["--email", email]
    query = _ask(console, "Custom query (optional)", "")
    if query:
        argv += ["--query", query]
    batch_size = _ask(console, "Batch size", "20")
    if batch_size:
        argv += ["--batch-size", batch_size]
    argv += _ask_extra_args(console)
    return argv


def _build_align_argv(console: Console) -> list[str]:
    inp = _ask_required(console, "Input FASTA path")
    out = _ask_required(console, "Output aligned FASTA path")
    mafft = _ask(console, "MAFFT binary", "mafft")
    mafft_args = _ask(console, "MAFFT args", "--auto")
    argv = ["align", "--in", inp, "--out", out]
    if mafft:
        argv += ["--mafft", mafft]
    if mafft_args:
        argv += ["--mafft-args", mafft_args]
    argv += _ask_extra_args(console)
    return argv


def _build_conserved_argv(console: Console) -> list[str]:
    inp = _ask_required(console, "Input aligned FASTA path")
    out = _ask_required(console, "Output regions JSON path")
    window = _ask_required(console, "Window size")
    quantile = _ask_required(console, "Top quantile (0..1)")
    argv = ["conserved", "--inp", inp, "--out", out, "--window", window, "--quantile", quantile]
    argv += _ask_extra_args(console)
    return argv


def _build_run_argv(console: Console) -> list[str]:
    gene_name = _ask_required(console, "Gene or comma-separated genes")
    workdir = _ask_required(console, "Workdir path")
    out = _ask_required(console, "Output directory")
    argv = ["run", "--gene_name", gene_name, "--workdir", workdir, "--out", out]
    max_records = _ask(console, "Max sequences", "100")
    if max_records:
        argv += ["--max", max_records]
    mafft = _ask(console, "MAFFT binary", "mafft")
    if mafft:
        argv += ["--mafft", mafft]
    mafft_args = _ask(console, "MAFFT args", "--auto")
    if mafft_args:
        argv += ["--mafft-args", mafft_args]
    query = _ask(console, "Custom query (optional)", "")
    if query:
        argv += ["--query", query]
    email = _ask(console, "NCBI email (optional)", "")
    if email:
        argv += ["--email", email]
    batch_size = _ask(console, "Batch size", "20")
    if batch_size:
        argv += ["--batch-size", batch_size]
    window = _ask(console, "Conserved window", "25")
    if window:
        argv += ["--window", window]
    quantile = _ask(console, "Conserved quantile", "0.8")
    if quantile:
        argv += ["--quantile", quantile]
    top_n = _ask(console, "Top N output", "20")
    if top_n:
        argv += ["--top-n", top_n]
    argv += _ask_extra_args(console)
    return argv


def _build_predict_argv(console: Console) -> list[str]:
    raw = _ask_required(console, "Raw FASTA path")
    alignment = _ask_required(console, "Alignment FASTA path")
    regions = _ask_required(console, "Conserved regions JSON path")
    out = _ask_required(console, "Output directory")
    argv = [
        "predict",
        "--raw",
        raw,
        "--alignment",
        alignment,
        "--regions",
        regions,
        "--out",
        out,
    ]
    top_n = _ask(console, "Top N output", "20")
    if top_n:
        argv += ["--top-n", top_n]
    argv += _ask_extra_args(console)
    return argv


@dataclass(frozen=True)
class _MenuItem:
    label: str
    argv_builder: Callable[[Console], list[str]] | None


_MENU: list[_MenuItem] = [
    _MenuItem("fetch", _build_fetch_argv),
    _MenuItem("align", _build_align_argv),
    _MenuItem("conserved", _build_conserved_argv),
    _MenuItem("run", _build_run_argv),
    _MenuItem("predict", _build_predict_argv),
    _MenuItem("quit", None),
]


def _read_key() -> str:
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(fd)
        ch = sys.stdin.read(1)
        if ch == "\x1b":
            ch2 = sys.stdin.read(1)
            ch3 = sys.stdin.read(1)
            if ch2 == "[" and ch3 == "A":
                return "up"
            if ch2 == "[" and ch3 == "B":
                return "down"
            return "esc"
        if ch in {"\r", "\n"}:
            return "enter"
        if ch in {"q", "Q"}:
            return "quit"
        if ch == "\x03":
            raise KeyboardInterrupt
        return ch
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)


def _select_menu_item(console: Console) -> _MenuItem:
    selected_idx = 0
    while True:
        console.clear()
        console.print(
            Panel.fit(
                "[bold cyan]Primer CLI Pretty Screen[/bold cyan]\n"
                "Interactive wrapper over standard CLI commands",
                border_style="cyan",
            )
        )
        console.print("[dim]Use ↑/↓ to navigate, Enter to choose, q to quit[/dim]")
        console.print()

        table = Table(title="Choose command", show_edge=True, border_style="cyan")
        table.add_column("", width=3, justify="center")
        table.add_column("Command", style="bold white")
        for i, item in enumerate(_MENU):
            pointer = "[bold cyan]▸[/bold cyan]" if i == selected_idx else " "
            label = f"[bold cyan]{item.label}[/bold cyan]" if i == selected_idx else item.label
            table.add_row(pointer, label)
        console.print(table)

        key = _read_key()
        if key == "up":
            selected_idx = (selected_idx - 1) % len(_MENU)
        elif key == "down":
            selected_idx = (selected_idx + 1) % len(_MENU)
        elif key == "enter":
            return _MENU[selected_idx]
        elif key == "quit":
            return _MENU[-1]


def run_pretty_screen(
    parser: argparse.ArgumentParser,
    run_handler: Runner,
    *,
    default_log_level: str = "INFO",
) -> int:
    if Console is None:
        raise PrimerCliError(
            "Pretty screen requires 'rich'. Install dependency and retry (pip install rich)."
        )
    if not sys.stdin.isatty():
        raise PrimerCliError("Pretty screen requires an interactive TTY terminal")

    console = Console()

    while True:
        selected = _select_menu_item(console)
        if selected.argv_builder is None:
            console.clear()
            console.print("[green]Bye.[/green]")
            return 0

        try:
            console.clear()
            argv = selected.argv_builder(console)
            log_level = Prompt.ask("Log level", default=default_log_level, console=console).upper()
            if log_level:
                argv += ["--log-level", log_level]

            console.print()
            console.print(
                f"[bold magenta]Running:[/bold magenta] "
                f"[white]primer-cli {' '.join(shlex.quote(x) for x in argv)}[/white]"
            )
            args = parser.parse_args(argv)
            handler: Handler | None = getattr(args, "func", None)
            if handler is None:
                console.print("[bold red]No handler found for selected command.[/bold red]")
                continue
            rc = run_handler(handler, args)
            if rc == 0:
                console.print("[bold green]Completed successfully.[/bold green]")
            else:
                console.print(f"[bold yellow]Command finished with code {rc}.[/bold yellow]")
        except KeyboardInterrupt:
            console.print()
            console.print("[bold yellow]Interrupted by user.[/bold yellow]")
            return 130
        except Exception as e:
            console.print(f"[bold red]Error:[/bold red] {e}")

        console.print()
        again = Confirm.ask("Run another command?", default=True, console=console)
        if not again:
            console.print("[green]Bye.[/green]")
            return 0
        console.print()


def has_pretty_flag(argv: list[str]) -> bool:
    return "--pretty-screen" in argv


def clean_pretty_flag(argv: list[str]) -> list[str]:
    return [x for x in argv if x != "--pretty-screen"]
