"""Client for an Automat graph deployment (ROBOKOP / RoboMouse KG).

The Cypher endpoint accepts only a raw query string -- it has no parameter
binding -- so every CURIE interpolated into a query is validated against
`CURIE_RE` first and then emitted as a JSON string literal. Anything carrying a
quote, backslash, or brace is rejected before it can reach the database.
"""

from __future__ import annotations

import json
import re
from typing import Any

import httpx

#: Deliberately excludes quotes, backslashes, braces, and whitespace.
CURIE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._+-]{0,31}:[A-Za-z0-9._~%+*#/()=,;:-]{1,128}$")


class AutomatError(RuntimeError):
    pass


class InvalidCurieError(ValueError):
    pass


def validate_curie(curie: str) -> str:
    curie = curie.strip()
    if not CURIE_RE.match(curie):
        raise InvalidCurieError(f"Not a well-formed CURIE: {curie!r}")
    return curie


def cypher_literal(value: str) -> str:
    """Render a validated string as a Cypher literal."""
    return json.dumps(value)


def curie_list_literal(curies: list[str]) -> str:
    return "[" + ", ".join(cypher_literal(validate_curie(c)) for c in curies) + "]"


class AutomatClient:
    def __init__(self, graph_url: str, client: httpx.AsyncClient) -> None:
        self._graph_url = graph_url.rstrip("/")
        self._client = client

    async def cypher(self, query: str) -> list[dict[str, Any]]:
        """POST /cypher and flatten the Neo4j-style response to dicts."""
        try:
            response = await self._client.post(f"{self._graph_url}/cypher", json={"query": query})
            response.raise_for_status()
            payload = response.json()
        except httpx.HTTPError as exc:
            raise AutomatError(f"Cypher request failed: {exc}") from exc

        if payload.get("errors"):
            raise AutomatError(f"Cypher error: {payload['errors']}")

        results = payload.get("results") or []
        if not results:
            return []
        block = results[0]
        columns = block.get("columns", [])
        return [dict(zip(columns, row["row"])) for row in block.get("data", [])]

    async def node(self, curie: str) -> dict[str, Any] | None:
        curie = validate_curie(curie)
        try:
            response = await self._client.get(f"{self._graph_url}/node/{curie}")
        except httpx.HTTPError as exc:
            raise AutomatError(f"Node lookup failed for {curie}: {exc}") from exc
        if not response.is_success:
            return None
        payload = response.json()
        return payload or None

    async def graph_version(self) -> str | None:
        try:
            response = await self._client.get(f"{self._graph_url}/graph-metadata")
            response.raise_for_status()
        except httpx.HTTPError:
            return None
        payload = response.json()
        version = payload.get("version")
        name = payload.get("name")
        if name and version:
            return f"{name} {version}"
        return version or name