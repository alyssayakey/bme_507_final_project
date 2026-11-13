"""Active research scope; upstream multi-drug data remain archived unchanged."""
from .common import DATA

DRUG = "quinidine"
DRUG_TABLE = DATA / "ap/quinidine.csv"
ALLOWED_CONDITIONS = ("control", DRUG)


def validate_condition(drug):
    if drug not in ALLOWED_CONDITIONS:
        raise ValueError("This study supports quinidine and its no-drug control only")
