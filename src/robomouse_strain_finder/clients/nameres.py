"""Client for the Translator NameResolution service.

https://github.com/NCATSTranslator/NameResolution
"""

from __future__ import annotations

import httpx

from ..models import Candidate, Species


class NameResolverError(RuntimeError):
    pass


class NameResolverClient:
    def __init__(self, base_url: str, client: httpx.AsyncClient) -> None:
        self._base_url = base_url.rstrip("/")
        self._client = client

    async def lookup(
        self,
        string: str,
        *,
        species: Species | None = None,
        limit: int = 10,
        autocomplete: bool = False,
        biolink_types: tuple[str, ...] = (
            "biolink:Disease",
            "biolink:PhenotypicFeature",
        ),
        restrict_prefixes: bool = True,
    ) -> list[Candidate]:
        """Resolve free text to ranked CURIE candidates.

        `/lookup` takes query parameters, not a JSON body. Prefix and taxon
        filters are pipe-separated and case-sensitive.

        Disease and phenotype terms are not taxon-tagged the way genes are, so
        `species` is applied as an ontology-prefix filter (mouse -> MP, human ->
        MONDO/HP) rather than via `only_taxa`.
        """
        params: list[tuple[str, str]] = [
            ("string", string),
            ("limit", str(limit)),
            ("autocomplete", str(autocomplete).lower()),
        ]
        for biolink_type in biolink_types:
            params.append(("biolink_type", biolink_type))
        if species is not None and restrict_prefixes:
            params.append(("only_prefixes", "|".join(species.prefixes)))

        try:
            response = await self._client.post(f"{self._base_url}/lookup", params=params)
            response.raise_for_status()
            payload = response.json()
        except httpx.HTTPError as exc:
            raise NameResolverError(f"Name Resolver lookup failed for {string!r}: {exc}") from exc

        return [
            Candidate(
                curie=hit["curie"],
                label=hit.get("label") or hit["curie"],
                score=hit.get("score", 0.0),
                types=hit.get("types", [])[:4],
                synonyms=hit.get("synonyms", [])[:8],
                taxa=hit.get("taxa", []),
            )
            for hit in payload
        ]