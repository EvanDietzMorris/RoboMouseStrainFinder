"""Application settings.

Every remote endpoint is overridable so the app can be pointed at a different
Automat deployment or a local NameRes instance without code changes.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

DEFAULT_CATALOG_PATH = (
    Path.home() / ".cache" / "query-mmrrc-catalog" / "mmrrc_catalog_data.csv"
)


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="RMSF_", env_file=".env", extra="ignore")

    name_resolver_url: str = "https://name-resolution-exp.apps.renci.org"
    automat_url: str = "https://robokop-automat.apps.renci.org"
    graph: str = "robomousekg"

    mmrrc_catalog_path: Path = DEFAULT_CATALOG_PATH
    mmrrc_catalog_url: str = "https://www.mmrrc.org/about/mmrrc_catalog_data.csv"
    mmrrc_api_url: str = "https://api.mmrrc.org"

    # Download the catalog on startup if it is not already on disk. The file is
    # ~147 MB, so this is opt-in rather than silent.
    auto_download_catalog: bool = True

    http_timeout: float = 120.0

    # Corporate TLS interception breaks Python's bundled trust store while the
    # macOS keychain still validates fine. truststore delegates to the OS store.
    use_system_trust_store: bool = True

    @property
    def graph_url(self) -> str:
        return f"{self.automat_url.rstrip('/')}/{self.graph}"


@lru_cache
def get_settings() -> Settings:
    return Settings()