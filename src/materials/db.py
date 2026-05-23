"""Database helper utilities for materials modules."""

from sqlalchemy import MetaData, Table, create_engine


def connect_sqlite_table(db_path: str, table_name: str):
    """Open SQLite connection and reflect a table."""
    metadata = MetaData()
    engine = create_engine(f"sqlite+pysqlite:///{db_path}")
    connection = engine.connect()
    table = Table(table_name, metadata, autoload_with=engine)
    return engine, connection, table


def disconnect(engine, connection):
    """Close SQLAlchemy connection and dispose engine."""
    connection.close()
    engine.dispose()
