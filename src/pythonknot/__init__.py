import os
from pathlib import Path


_TABLE_ENV = "PYTHONKNOT_ALEXANDER_TABLE"
_TABLE_NAME = "table_knot_Alexander_polynomial.txt"


def _configure_default_alexander_table() -> None:
    if os.environ.get(_TABLE_ENV):
        return

    packaged_table = Path(__file__).resolve().parent / "data" / _TABLE_NAME
    if packaged_table.exists():
        os.environ[_TABLE_ENV] = str(packaged_table)


_configure_default_alexander_table()
