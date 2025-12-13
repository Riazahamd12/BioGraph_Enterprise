from fastapi.testclient import TestClient
from main import app
import os
import pytest

client = TestClient(app)

def test_home():
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {"status": "BioGraph AI Engine Ready 🟢"}

def test_analyze_manual_valid():
    payload = {
        "mode": "manual",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O", # Aspirin
        "target_id": "6LU7"
    }
    response = client.post("/analyze", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert data["type"] == "single"
    assert "score" in data
    assert "confidence" in data
    assert "graph" in data

def test_analyze_manual_invalid_structure():
    payload = {
        "mode": "manual",
        "smiles": "INVALID_SMILES",
        "target_id": "6LU7"
    }
    response = client.post("/analyze", json=payload)
    assert response.status_code == 200
    assert response.json() == {"error": "Invalid Structure"}

def test_get_image_valid():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    response = client.get(f"/get_image?smiles={smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "image/png"
    # Check if we got some bytes
    assert len(response.content) > 0

def test_get_image_invalid():
    smiles = "INVALID_SMILES_STRING"
    response = client.get(f"/get_image?smiles={smiles}")
    assert response.status_code == 404

def test_analyze_auto_no_file(monkeypatch):
    # Mock os.path.exists to return False
    monkeypatch.setattr(os.path, "exists", lambda x: False)
    # Also ensure DRUGS_DF is None (it might be loaded from real run)
    import main
    monkeypatch.setattr(main, "DRUGS_DF", None)

    payload = {"mode": "auto"}
    response = client.post("/analyze", json=payload)
    assert response.status_code == 200
    assert response.json() == {"error": "drugs.csv not found"}
