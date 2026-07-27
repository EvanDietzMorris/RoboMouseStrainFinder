"""The Cypher endpoint has no parameter binding, so CURIE validation is the
only thing standing between a query string and injected Cypher."""

from __future__ import annotations

import httpx
import pytest

from robomouse_strain_finder.clients.automat import (
    AutomatClient,
    AutomatError,
    InvalidCurieError,
    curie_list_literal,
    validate_curie,
)


@pytest.mark.parametrize(
    "curie",
    ["MONDO:0005027", "MP:0001261", "NCBIGene:57249", "MGI:1888498", "HP:0000002"],
)
def test_real_curies_validate(curie: str) -> None:
    assert validate_curie(curie) == curie


@pytest.mark.parametrize(
    "curie",
    [
        "MONDO:0005027' RETURN 1 //",
        'MONDO:0005027" OR true',
        "MONDO:0005027\\",
        "no-colon",
        "",
        "  ",
        ":0005027",
        "MONDO:0005027 DETACH DELETE n",
        "MONDO:{id}",
        "MONDO:0005027\nMATCH (n) DELETE n",
    ],
)
def test_injection_shaped_input_is_rejected(curie: str) -> None:
    with pytest.raises(InvalidCurieError):
        validate_curie(curie)


def test_list_literal_renders_quoted_curies() -> None:
    assert curie_list_literal(["MP:0001261", "HP:0000002"]) == '["MP:0001261", "HP:0000002"]'


def test_list_literal_rejects_a_bad_member() -> None:
    with pytest.raises(InvalidCurieError):
        curie_list_literal(["MP:0001261", "x' OR '1'='1"])


@pytest.mark.asyncio
async def test_cypher_flattens_the_neo4j_response() -> None:
    payload = {
        "results": [
            {
                "columns": ["curie", "symbol"],
                "data": [
                    {"row": ["NCBIGene:57249", "Gabrq"]},
                    {"row": ["NCBIGene:12367", "Casp3"]},
                ],
            }
        ],
        "errors": [],
    }
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json=payload))
    async with httpx.AsyncClient(transport=transport) as http:
        client = AutomatClient("https://example.org/robomousekg", http)
        rows = await client.cypher("MATCH (n) RETURN n")
    assert rows == [
        {"curie": "NCBIGene:57249", "symbol": "Gabrq"},
        {"curie": "NCBIGene:12367", "symbol": "Casp3"},
    ]


@pytest.mark.asyncio
async def test_cypher_surfaces_server_errors() -> None:
    payload = {"results": [], "errors": [{"code": "SyntaxError", "message": "bad"}]}
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json=payload))
    async with httpx.AsyncClient(transport=transport) as http:
        client = AutomatClient("https://example.org/robomousekg", http)
        with pytest.raises(AutomatError, match="SyntaxError"):
            await client.cypher("MATCH (n) RETURN n")


@pytest.mark.asyncio
async def test_missing_node_returns_none_not_an_empty_dict() -> None:
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json={}))
    async with httpx.AsyncClient(transport=transport) as http:
        client = AutomatClient("https://example.org/robomousekg", http)
        assert await client.node("MGI:1888498") is None