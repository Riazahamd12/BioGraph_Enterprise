from fastapi import FastAPI, UploadFile, File, Form 
from fastapi.responses import Response 
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
import random
import time
import pandas as pd
import io
import os
import requests 

# --- RDKit IMPORTS ---
from rdkit import Chem
from rdkit.Chem import Descriptors, QED
from rdkit.Chem.Draw import rdMolDraw2D

app = FastAPI()

app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)

class AnalysisRequest(BaseModel):
    target_id: str = None
    smiles: str = None
    mode: str

# --- GLOBAL DATA CACHE ---
DRUGS_DF = None

def load_drugs_data():
    global DRUGS_DF
    if os.path.exists("drugs.csv"):
        try:
            DRUGS_DF = pd.read_csv("drugs.csv")
            print("Loaded drugs.csv into memory.")
        except Exception as e:
            print(f"Error loading drugs.csv: {e}")

# Load data on startup
load_drugs_data()

# --- HELPER: RDKit PROPERTIES ---
def calculate_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        return {
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "qed": round(QED.qed(mol), 2)
        }
    except:
        return None

def process_drug_row(row):
    """Helper to process a drug row from DataFrame"""
    score = round(random.uniform(5.0, 9.9), 2)
    return {
        "name": row['name'],
        "smiles": row['smiles'],
        "score": score,
        "status": "ACTIVE" if score > 7.5 else "INACTIVE"
    }

@app.get("/")
def home():
    return {"status": "BioGraph AI Engine Ready 🟢"}

@app.get("/get_image")
def get_image(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return Response(status_code=404)

        # Use MolDraw2DCairo for transparent PNG support
        drawer = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
        opts = drawer.drawOptions()
        opts.clearBackground = False # Transparent background

        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()

        img_byte_arr = drawer.GetDrawingText()

        return Response(content=img_byte_arr, media_type="image/png")
    except Exception as e:
        print(f"Error generating image: {e}")
        return Response(status_code=500)

@app.post("/analyze")
def analyze_drug(data: AnalysisRequest):
    # 1. MANUAL MODE
    if data.mode == "manual":
        # Removed artificial delay time.sleep(1)
        props = calculate_properties(data.smiles)
        if not props: return {"error": "Invalid Structure"}
        
        score = round(random.uniform(6.0, 9.8), 2)
        
        # --- NEW: CONFIDENCE CALCULATION ---
        # Score jitna high, Confidence utna high (Simulation)
        confidence_val = int((score / 10) * 100) + random.randint(-2, 2)
        confidence_val = min(max(confidence_val, 85), 99) # 85% se 99% ke beech rakho

        return {
            "type": "single",
            "score": score,
            "status": "ACTIVE" if score > 7.0 else "INACTIVE",
            "confidence": f"{confidence_val}%", # Ab ye Percentage bheje ga
            "color": "#00f3ff" if score > 7.0 else "#ff0055",
            "graph": [
                {"subject": "Binding", "A": int(score * 10), "fullMark": 100},
                {"subject": "LogP", "A": min(int(((props['logp'] + 2)/8)*100), 100), "fullMark": 100},
                {"subject": "Mol.Weight", "A": min(int((props['mw']/800)*100), 100), "fullMark": 100},
                {"subject": "Toxicity", "A": random.randint(10, 40), "fullMark": 100},
                {"subject": "QED", "A": int(props['qed'] * 100), "fullMark": 100},
            ]
        }

    # 2. AUTO SCAN MODE
    elif data.mode == "auto":
        # Removed artificial delay time.sleep(1.5)
        results = []
        df = None

        if DRUGS_DF is not None:
            df = DRUGS_DF
        elif os.path.exists("drugs.csv"):
             # Fallback if global load failed
            df = pd.read_csv("drugs.csv")

        if df is not None:
            for _, row in df.iterrows():
                results.append(process_drug_row(row))
            results.sort(key=lambda x: x["score"], reverse=True)
            return {"type": "batch", "results": results[:20]}
        else:
            return {"error": "drugs.csv not found"}

@app.post("/upload")
async def upload_file(target_id: str = Form(...), file: UploadFile = File(...)):
    try:
        contents = await file.read()
        if file.filename.endswith('.csv'):
            df = pd.read_csv(io.BytesIO(contents))
        elif file.filename.endswith(('.xls', '.xlsx')):
            df = pd.read_excel(io.BytesIO(contents))
        else:
            return {"error": "Invalid file format"}

        results = []
        for _, row in df.head(50).iterrows():
            name = str(row.iloc[0])
            smiles = str(row.iloc[1]) 
            props = calculate_properties(smiles)
            if props:
                score = round(random.uniform(4.0, 9.5), 2)
                results.append({
                    "name": name,
                    "score": score,
                    "status": "ACTIVE" if score > 7.5 else "INACTIVE"
                })
        
        results.sort(key=lambda x: x["score"], reverse=True)
        return {"type": "batch", "results": results}

    except Exception as e:
        return {"error": str(e)}
