"""Route-level tests. The app state is populated directly so no network or
147 MB catalog load is involved."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from robomouse_strain_finder.catalog import MmrrcCatalog
from robomouse_strain_finder.clients.automat import AutomatError
from robomouse_strain_finder.config import Settings
from robomouse_strain_finder.main import app
from robomouse_strain_finder.models import Candidate, Species
from robomouse_strain_finder.pipeline import StrainFinder

from .test_pipeline import HEADER, FakeAutomat, FakeNameResolver


@pytest.fixture
def client(tmp_path: Path) -> TestClient:
    path = tmp_path / "catalog.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        row = dict.fromkeys(HEADER, "")
        row["STRAIN/STOCK_ID"] = "MMRRC:000012-UNC"
        row["MGI_GENE_ACCESSION_ID"] = "MGI:106185"
        row["GENE_SYMBOL"] = "Ccr2"
        writer.writerow([row[column] for column in HEADER])
    catalog = MmrrcCatalog(path)
    catalog.load()

    settings = Settings(mmrrc_catalog_path=path, auto_download_catalog=False)
    app.state.settings = settings
    app.state.catalog = catalog
    app.state.finder = StrainFinder(
        settings=settings,
        nameres=FakeNameResolver(
            [Candidate(curie="MONDO:0005027", label="epilepsy", score=1.0)]
        ),
        automat=FakeAutomat({}),
        catalog=catalog,
    )
    # The lifespan would overwrite app.state, so it is intentionally not run.
    return TestClient(app)


def test_healthz_reports_catalog_state(client: TestClient) -> None:
    body = client.get("/healthz").json()
    assert body["status"] == "ok"
    assert body["catalog_stocks"] == 1


def test_resolve_returns_candidates(client: TestClient) -> None:
    body = client.get("/api/resolve", params={"term": "epilepsy", "species": "human"}).json()
    assert body["species"] == "human"
    assert body["candidates"][0]["curie"] == "MONDO:0005027"


def test_search_returns_a_full_envelope(client: TestClient) -> None:
    body = client.get("/api/search", params={"term": "epilepsy", "species": "human"}).json()
    assert body["selected"]["curie"] == "MONDO:0005027"
    assert "provenance" in body
    assert body["provenance"]["catalog_rows"] == 1


def test_a_malformed_curie_is_rejected_before_reaching_the_graph(client: TestClient) -> None:
    response = client.get(
        "/api/search", params={"term": "x", "curie": "MONDO:0005027' RETURN 1 //"}
    )
    assert response.status_code == 422


def test_an_unknown_species_is_rejected(client: TestClient) -> None:
    assert client.get("/api/search", params={"term": "x", "species": "rat"}).status_code == 422


def test_a_blank_term_is_rejected(client: TestClient) -> None:
    assert client.get("/api/search", params={"term": ""}).status_code == 422


def test_graph_failures_surface_as_bad_gateway(client: TestClient) -> None:
    class BrokenAutomat(FakeAutomat):
        async def cypher(self, query: str):
            raise AutomatError("upstream graph is down")

    app.state.finder = StrainFinder(
        settings=app.state.settings,
        nameres=FakeNameResolver(
            [Candidate(curie="MONDO:0005027", label="epilepsy", score=1.0)]
        ),
        automat=BrokenAutomat({}),
        catalog=app.state.catalog,
    )
    response = client.get("/api/search", params={"term": "epilepsy", "species": "human"})
    assert response.status_code == 502
    assert "upstream graph is down" in response.json()["detail"]


def test_index_serves_the_spa_shell(client: TestClient) -> None:
    response = client.get("/app")
    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/html")
    assert "RoboMouse Strain Finder" in response.text
    # The shell is static: results arrive from the API in the browser.
    assert "/static/app.js" in response.text
    assert "MONDO:0005027" not in response.text


def test_root_serves_the_landing_page(client: TestClient) -> None:
    response = client.get("/")
    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/html")
    assert "RoboMouse Strain Finder" in response.text
    # The landing page is inert: no client script, and it points at the tool.
    assert "/static/app.js" not in response.text
    assert 'href="/app"' in response.text


@pytest.mark.parametrize("asset", ["/static/app.js", "/static/app.css", "/static/landing.css"])
def test_static_assets_are_served(client: TestClient, asset: str) -> None:
    response = client.get(asset)
    assert response.status_code == 200
    assert response.content


def test_the_ui_never_builds_markup_from_api_data(client: TestClient) -> None:
    """Catalog designations contain HTML, so the client must not use innerHTML."""
    source = client.get("/static/app.js").text
    assert ".innerHTML" not in source
    assert "insertAdjacentHTML" not in source


def test_the_ui_keeps_shareable_urls_on_the_app_route(client: TestClient) -> None:
    """The tool lives at /app; syncUrl must not rewrite deep links back to /."""
    source = client.get("/static/app.js").text
    assert "`/app?${keep}`" in source


def test_search_accepts_a_species_enum_value(client: TestClient) -> None:
    for species in (Species.MOUSE, Species.HUMAN):
        response = client.get("/api/search", params={"term": "epilepsy", "species": species.value})
        assert response.status_code == 200