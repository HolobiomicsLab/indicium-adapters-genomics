"""Tests that the genomics adapter is registered and discoverable."""
from __future__ import annotations

from importlib.metadata import entry_points

from indicium_adapters import discover_adapters


def test_entry_point_registered():
    names = {ep.name for ep in entry_points(group="indicium.adapters")
             if ep.module.startswith("indicium_adapters_genomics")}
    assert names == {"genomics"}


def test_discover_adapters_finds_genomics():
    found = discover_adapters()
    assert "genomics" in found, f"discover_adapters() missing 'genomics'; found={sorted(found)}"


def test_discovered_adapter_has_expected_api_version():
    assert discover_adapters()["genomics"].api_version == 1
