"""
utils/citation.py
-----------------
How to cite CASCADE: shared by the --citation flag of all three CLI tools.
"""

import argparse

from cascade import __version__

DOI = "10.1093/bioinformatics/btag570"

CITATION = (
    "Avila F, Kim J, Park J-C, Kang R, Lee H, Park S-H. CASCADE: criticality "
    "avalanche spike cross-platform analysis detection engine, a "
    "multi-manufacturer MEA bash analysis pipeline. Bioinformatics. "
    f"2026;42(8):btag570. https://doi.org/{DOI}"
)

BIBTEX = f"""@article{{Avila2026CASCADE,
  author  = {{Avila, Forbes and Kim, Jin and Park, Jong-Chan and Kang, Rian and Lee, Hyunsu and Park, Sun-Hyun}},
  title   = {{{{CASCADE}}: criticality avalanche spike cross-platform analysis detection engine, a multi-manufacturer {{MEA}} bash analysis pipeline}},
  journal = {{Bioinformatics}},
  year    = {{2026}},
  volume  = {{42}},
  number  = {{8}},
  pages   = {{btag570}},
  doi     = {{{DOI}}}
}}"""


def citation_text() -> str:
    return (
        "If you use CASCADE in your research, please cite:\n\n"
        f"  {CITATION}\n\n"
        "BibTeX:\n\n"
        f"{BIBTEX}\n\n"
        f"Please also report the CASCADE version you used (this is v{__version__}),\n"
        "because criticality metrics differ between versions (see CHANGELOG.md).\n"
    )


class CitationAction(argparse.Action):
    """argparse action: print the citation and exit, like --version."""

    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help="Show how to cite CASCADE and exit."):
        super().__init__(option_strings=option_strings, dest=dest,
                         default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print(citation_text(), end="")
        parser.exit()


def add_citation_argument(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--citation", action=CitationAction)
