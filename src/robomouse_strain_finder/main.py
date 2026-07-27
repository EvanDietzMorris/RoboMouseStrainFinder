"""FastAPI application wiring.

The MMRRC catalog is loaded once into memory during startup (it takes a few
seconds and roughly 150 MB on disk), then shared by every request.
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Annotated

import httpx
from fastapi import Depends, FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .catalog import MmrrcCatalog, load_catalog
from .clients import AutomatClient, NameResolverClient
from .clients.automat import AutomatError, InvalidCurieError
from .clients.nameres import NameResolverError
from .config import Settings, get_settings
from .models import Candidate, ResolveResponse, SearchResponse, Species
from .pipeline import StrainFinder

logger = logging.getLogger(__name__)

STATIC_DIR = Path(__file__).parent / "static"


@asynccontextmanager
async def lifespan(app: FastAPI):
    settings = get_settings()

    if settings.use_system_trust_store:
        # Networks that terminate TLS with their own CA break Python's bundled
        # trust store while the OS store still validates. Opt out with
        # RMSF_USE_SYSTEM_TRUST_STORE=false.
        import truststore

        truststore.inject_into_ssl()

    logger.info("Loading MMRRC catalog from %s", settings.mmrrc_catalog_path)
    catalog = load_catalog(
        settings.mmrrc_catalog_path,
        settings.mmrrc_catalog_url,
        settings.auto_download_catalog,
    )
    logger.info("Loaded %d catalog rows for %d stocks", catalog.rows, len(catalog.stocks))

    async with httpx.AsyncClient(timeout=settings.http_timeout, follow_redirects=True) as client:
        app.state.settings = settings
        app.state.catalog = catalog
        app.state.finder = StrainFinder(
            settings=settings,
            nameres=NameResolverClient(settings.name_resolver_url, client),
            automat=AutomatClient(settings.graph_url, client),
            catalog=catalog,
        )
        yield


app = FastAPI(
    title="RoboMouse Strain Finder",
    version="0.1.0",
    description=(
        "Resolve a phenotype or disease to a CURIE, traverse RoboMouse KG for "
        "associated genes, orthologs and pathways, and surface the MMRRC mouse "
        "strains annotated with those genes."
    ),
    lifespan=lifespan,
)


def get_finder(request: Request) -> StrainFinder:
    return request.app.state.finder


def get_catalog(request: Request) -> MmrrcCatalog:
    return request.app.state.catalog


def app_settings(request: Request) -> Settings:
    return request.app.state.settings


FinderDep = Annotated[StrainFinder, Depends(get_finder)]


@app.exception_handler(InvalidCurieError)
async def _invalid_curie_handler(request: Request, exc: InvalidCurieError) -> JSONResponse:
    # A safety net for any path that does not already translate this itself.
    return JSONResponse(status_code=422, content={"detail": str(exc)})


@app.get("/healthz")
async def healthz(
    catalog: Annotated[MmrrcCatalog, Depends(get_catalog)],
    settings: Annotated[Settings, Depends(app_settings)],
) -> dict[str, object]:
    return {
        "status": "ok",
        "graph_url": settings.graph_url,
        "name_resolver_url": settings.name_resolver_url,
        "catalog_path": str(catalog.path),
        "catalog_rows": catalog.rows,
        "catalog_stocks": len(catalog.stocks),
        "catalog_modified_at": catalog.modified_at,
    }


@app.get("/api/resolve", response_model=ResolveResponse)
async def resolve(
    finder: FinderDep,
    term: Annotated[str, Query(min_length=1, description="Phenotype or disease text")],
    species: Species = Species.MOUSE,
    limit: Annotated[int, Query(ge=1, le=50)] = 10,
    restrict_prefixes: Annotated[
        bool,
        Query(
            description=(
                "Restrict to the species' phenotype vocabulary "
                "(mouse -> MP, human -> MONDO/HP)."
            )
        ),
    ] = True,
) -> ResolveResponse:
    """Resolve free text to ranked disease/phenotype CURIE candidates."""
    try:
        candidates: list[Candidate] = await finder.resolve(
            term, species, limit=limit, restrict_prefixes=restrict_prefixes
        )
    except NameResolverError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
    return ResolveResponse(query=term, species=species, candidates=candidates)


@app.get("/api/search", response_model=SearchResponse)
async def search(
    finder: FinderDep,
    term: Annotated[str, Query(min_length=1, description="Phenotype or disease text")],
    species: Species = Species.MOUSE,
    curie: Annotated[
        str | None,
        Query(description="Skip name resolution and traverse from this CURIE."),
    ] = None,
    max_genes: Annotated[int, Query(ge=1, le=500)] = 100,
    max_pathways: Annotated[int, Query(ge=0, le=200)] = 40,
    max_strains: Annotated[int, Query(ge=1, le=2000)] = 200,
    include_orthologs: bool = True,
    bridge_species: Annotated[
        bool,
        Query(description="Expand seed terms through UPheno cross-species phenotype homology."),
    ] = True,
    exclude_mutation_types: Annotated[
        list[str] | None,
        Query(
            description=(
                "Raw MMRRC mutation-type codes to drop, e.g. `CI` for the chemically "
                "induced ENU mutagenesis archive. A stock is dropped only when all of "
                "its mutation types are excluded."
            )
        ),
    ] = None,
    exclude_tool_lines: Annotated[
        bool,
        Query(
            description=(
                "Drop reporter and cre-driver transgenes (GENSAT `Tg(gene-EGFP)` "
                "lines and similar). They leave the gene intact, so they are research "
                "tools rather than models of its disease."
            )
        ),
    ] = False,
) -> SearchResponse:
    """Run the full phenotype -> genes -> pathways -> MMRRC strains traversal."""
    try:
        return await finder.search(
            term=term,
            species=species,
            curie=curie,
            max_genes=max_genes,
            max_pathways=max_pathways,
            max_strains=max_strains,
            include_orthologs=include_orthologs,
            bridge_species=bridge_species,
            exclude_mutation_types=frozenset(exclude_mutation_types or ()),
            exclude_tool_lines=exclude_tool_lines,
        )
    except InvalidCurieError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except (AutomatError, NameResolverError) as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/", include_in_schema=False)
async def landing() -> FileResponse:
    """Serve the marketing landing page. Static, no JS, links through to /app."""
    return FileResponse(STATIC_DIR / "landing.html", media_type="text/html")


@app.get("/app", include_in_schema=False)
async def index() -> FileResponse:
    """Serve the single-page UI shell; it fetches the API from the browser."""
    return FileResponse(STATIC_DIR / "index.html", media_type="text/html")


def main() -> None:
    """Entry point for `robomouse-strain-finder`."""
    import uvicorn

    logging.basicConfig(level=logging.INFO)
    uvicorn.run("robomouse_strain_finder.main:app", host="127.0.0.1", port=8000)