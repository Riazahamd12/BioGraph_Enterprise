from fastapi import FastAPI, UploadFile, File, HTTPException
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
import random
import time
import pandas as pd
import io
import os

# --- RDKit IMPORTS ---
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

app = FastAPI()

app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)

class AnalysisRequest(BaseModel):
    target_id: str = None
    smiles: str = None
    mode: str

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

@app.get("/")
def home():
    return {"status": "BioGraph AI Engine (Pandas + RDKit) 🟢"}

# --- MAIN ANALYSIS ENDPOINT ---
@app.post("/analyze")
def analyze_drug(data: AnalysisRequest):
    
    # 1. MANUAL MODE (Single Drug)
    if data.mode == "manual":
        time.sleep(1) # Simulation
        props = calculate_properties(data.smiles)
        if not props: return {"error": "Invalid Structure"}
        
        # Fake AI Score (Phase 4 main Real Model aayega)
        score = round(random.uniform(6.0, 9.8), 2)
        
        return {
            "type": "single",
            "score": score,
            "status": "ACTIVE" if score > 7.0 else "INACTIVE",
            "confidence": "RDKit Verified",
            "color": "#00f3ff" if score > 7.0 else "#ff0055",
            "graph": [
                {"subject": "Binding", "A": int(score * 10), "fullMark": 100},
                {"subject": "LogP", "A": min(int(((props['logp'] + 2)/8)*100), 100), "fullMark": 100},
                {"subject": "Mol.Weight", "A": min(int((props['mw']/800)*100), 100), "fullMark": 100},
                {"subject": "Toxicity", "A": random.randint(10, 40), "fullMark": 100},
                {"subject": "QED", "A": int(props['qed'] * 100), "fullMark": 100},
            ]
        }

    # 2. AUTO SCAN MODE (Read form CSV)
    elif data.mode == "auto":
        time.sleep(2) # Feel good delay
        results = []
        
        # CSV Read karo
        if os.path.exists("drugs.csv"):
            df = pd.read_csv("drugs.csv")
            
            # Har drug ko process karo
            for _, row in df.iterrows():
                props = calculate_properties(row['smiles'])
                if props:
                    # Fake Binding Score Logic
                    score = round(random.uniform(5.0, 9.9), 2)
                    results.append({
                        "name": row['name'],
                        "smiles": row['smiles'],
                        "score": score,
                        "status": "ACTIVE" if score > 7.5 else "INACTIVE"
                    })
            
            # Sort by Best Score
            results.sort(key=lambda x: x["score"], reverse=True)
            # Top 20 bhejo
            return {"type": "batch", "results": results[:20]}
        else:
            return {"error": "Database not found"}

# --- REAL FILE UPLOAD ENDPOINT ---
@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    try:
        contents = await file.read()
        
        # Excel ya CSV detect karo
        if file.filename.endswith('.csv'):
            df = pd.read_csv(io.BytesIO(contents))
        elif file.filename.endswith(('.xls', '.xlsx')):
            df = pd.read_excel(io.BytesIO(contents))
        else:
            return {"error": "Invalid file format. Use CSV or Excel."}

        results = []
        # Pehle 2 columns uthao (Name, SMILES)
        for _, row in df.head(50).iterrows(): # Limit to 50 for speed
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