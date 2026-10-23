"""Command line interface for the CABOMB V1 model."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

from cabomb_model.analysis import run_all_analyses, run_baseline
from cabomb_model.model_modes import ModelModes, V1_1_MODES, V1_2_MODES, V1_3_MODES
from cabomb_model.parameters import default_registry_path, load_registry, validate_registry
from cabomb_model.resource_paths import parameter_registry
from cabomb_model.reporting import build_html_report, build_markdown_report
from cabomb_model.v13_analysis import run_v13_all_analyses, run_v13_baseline
from cabomb_model.v13_reporting import build_v13_markdown_report


def _validate_registry() -> int:
    validation = validate_registry(load_registry())
    if validation.warnings:
        print("Warnings:")
        for warning in validation.warnings:
            print(f"- {warning}")
    if not validation.ok:
        print("Errors:")
        for error in validation.errors:
            print(f"- {error}")
        return 1
    print("Parameter registry validation passed.")
    return 0


def _modes_from_version(version: str) -> ModelModes:
    """Resolve the public model-version CLI flag to its explicit modes."""

    if version == "v1.1":
        return V1_1_MODES
    if version == "v1.2":
        return V1_2_MODES
    if version == "v1.3":
        return V1_3_MODES
    raise ValueError(f"Unknown model version {version!r}.")


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""

    parser = argparse.ArgumentParser(prog="cabomb-model")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validation_parser = subparsers.add_parser("validate-registry")

    baseline_parser = subparsers.add_parser("run-baseline")
    baseline_parser.add_argument("--output-dir", default="outputs")
    baseline_parser.add_argument(
        "--model-version", choices=("v1.1", "v1.2", "v1.3"), default="v1.1"
    )

    analysis_parser = subparsers.add_parser("run-analysis")
    analysis_parser.add_argument("--output-dir", default="outputs")
    analysis_parser.add_argument(
        "--model-version", choices=("v1.1", "v1.2", "v1.3"), default="v1.1"
    )

    report_parser = subparsers.add_parser("build-report")
    report_parser.add_argument("--output-dir", default="outputs")
    report_parser.add_argument("--report-path", default="reports/model_report.md")
    report_parser.add_argument("--html-path", default="reports/model_report.html")
    report_parser.add_argument(
        "--model-version", choices=("v1.1", "v1.2", "v1.3"), default="v1.1"
    )

    for command_parser in (validation_parser, baseline_parser, analysis_parser, report_parser):
        command_parser.add_argument(
            "--registry", type=Path,
            help="Use this parameter CSV instead of the bundled default registry.",
        )

    args = parser.parse_args(argv)
    with parameter_registry(args.registry):
        result = _execute(args)
        if result == 0 and args.command != "validate-registry":
            output = Path(args.output_dir)
            output.mkdir(parents=True, exist_ok=True)
            registry = default_registry_path()
            snapshot = output / "registry_used.csv"
            if registry.resolve() != snapshot.resolve():
                shutil.copyfile(registry, snapshot)
            record = {
                "selection": "explicit --registry" if args.registry else "bundled default",
                "sha256": hashlib.sha256(registry.read_bytes()).hexdigest(),
                "snapshot": "registry_used.csv",
            }
            (output / "registry_source.json").write_text(
                json.dumps(record, indent=2) + "\n", encoding="utf-8"
            )
        return result


def _execute(args: argparse.Namespace) -> int:
    """Execute one parsed command under its selected parameter registry."""
    if args.command == "validate-registry":
        return _validate_registry()
    if args.command == "run-baseline":
        if args.model_version == "v1.3":
            outputs = run_v13_baseline(Path(args.output_dir))
        else:
            outputs = run_baseline(
                Path(args.output_dir), modes=_modes_from_version(args.model_version)
            )
        print("Baseline simulation completed.")
        for key, value in outputs.items():
            if key not in ("result", "params"):
                print(f"{key}: {value}")
        return 0
    if args.command == "run-analysis":
        if args.model_version == "v1.3":
            run_v13_all_analyses(Path(args.output_dir))
        else:
            run_all_analyses(
                Path(args.output_dir), modes=_modes_from_version(args.model_version)
            )
        print(f"Analysis outputs written under {args.output_dir}.")
        return 0
    if args.command == "build-report":
        if args.model_version == "v1.3":
            markdown = build_v13_markdown_report(
                args.output_dir, args.report_path
            )
            html = build_html_report(
                markdown, args.html_path, title="CABOMB V1.3 Model Report"
            )
        else:
            markdown = build_markdown_report(
                args.output_dir, args.report_path, modes=_modes_from_version(args.model_version)
            )
            html = build_html_report(markdown, args.html_path)
        print(f"Report written to {markdown} and {html}.")
        return 0
    raise AssertionError(f"Unhandled command {args.command}")

