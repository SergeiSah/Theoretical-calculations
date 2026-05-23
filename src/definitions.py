from functools import lru_cache
from pathlib import Path

from pydantic import Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ENV_FILE = PROJECT_ROOT / ".env"


class ProjectPaths(BaseSettings):
    """Project filesystem paths."""

    model_config = SettingsConfigDict(
        env_file=ENV_FILE,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    root_dir: Path = Field(default=PROJECT_ROOT, validation_alias="ROOT_DIR")
    data_dir: Path | None = Field(default=None, validation_alias="DATA_DIR")
    compounds_db_path: Path | None = Field(default=None, validation_alias="COMPOUNDS_DB_PATH")
    elements_db_path: Path | None = Field(default=None, validation_alias="ELEMENTS_DB_PATH")
    pics_db_path: Path | None = Field(default=None, validation_alias="PICS_DB_PATH")
    atomic_factors_db_path: Path | None = Field(default=None, validation_alias="ATOMIC_FACTORS_DB_PATH")
    f1f2_dir: Path | None = Field(default=None, validation_alias="F1F2_DIR")

    @model_validator(mode="after")
    def set_default_paths(self):
        self.data_dir = self.data_dir or self.root_dir / "data"
        self.compounds_db_path = self.compounds_db_path or self.data_dir / "Compounds.db"
        self.elements_db_path = self.elements_db_path or self.data_dir / "Elements.db"
        self.pics_db_path = self.pics_db_path or self.data_dir / "PICS.sqlite"
        self.atomic_factors_db_path = self.atomic_factors_db_path or self.data_dir / "Atomic_factors.db"
        self.f1f2_dir = self.f1f2_dir or self.data_dir / "f1f2"
        return self


class MaterialsProjectCredentials(BaseSettings):
    """Materials Project credentials loaded from the project root .env file."""

    model_config = SettingsConfigDict(
        env_file=ENV_FILE,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    api_key: str | None = Field(default=None, validation_alias="MAPI_KEY")


@lru_cache
def get_project_paths() -> ProjectPaths:
    return ProjectPaths()


@lru_cache
def get_materials_project_credentials() -> MaterialsProjectCredentials:
    return MaterialsProjectCredentials()


PROJECT_PATHS = get_project_paths()
MATERIALS_PROJECT_CREDENTIALS = get_materials_project_credentials()

ROOT_DIR = str(PROJECT_PATHS.root_dir)
DATA_DIR = str(PROJECT_PATHS.data_dir)
COMPOUNDS_DB_PATH = str(PROJECT_PATHS.compounds_db_path)
ELEMENTS_DB_PATH = str(PROJECT_PATHS.elements_db_path)
PICS_DB_PATH = str(PROJECT_PATHS.pics_db_path)
ATOMIC_FACTORS_DB_PATH = str(PROJECT_PATHS.atomic_factors_db_path)
F1F2_DIR = str(PROJECT_PATHS.f1f2_dir)
