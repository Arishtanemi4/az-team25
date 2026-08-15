"""Response models for /genes endpoints."""

from pydantic import BaseModel


class GeneMatch(BaseModel):
    """One gene search result -- exactly the two fields a combobox needs to display a choice
    and, once picked, resolve it unambiguously (ensembl_id is what scoring_service actually
    scores on; symbol is what a researcher recognises)."""

    ensembl_id: str
    symbol: str
