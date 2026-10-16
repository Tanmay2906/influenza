"""Shared helpers: config loading and NCBI access."""
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent.parent


def load_config(path=None):
    with open(path or ROOT / "config.yaml") as fh:
        return yaml.safe_load(fh)


def use_system_trust_store():
    """Use the OS certificate store for HTTPS (needed behind TLS-inspecting proxies)."""
    try:
        import truststore
        truststore.inject_into_ssl()
    except ImportError:
        pass
