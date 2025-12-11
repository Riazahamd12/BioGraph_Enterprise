import React from 'react';
import { Database, Search, Zap, Upload } from 'lucide-react';

/**
 * Sidebar: left input panel. Presentational wrapper so App.jsx is shorter.
 * All handlers are passed as props (same signatures as before).
 */
export default function Sidebar({
  activeTab,
  setActiveTab,
  target,
  setTarget,
  smiles,
  setSmiles,
  selectedFile,
  fileInputRef,
  handleFileSelect,
  handleScan,
  loading
}) {
  return (
    <div className="glass-panel panel-left" style={{ zIndex: 50 }}>
      <div className="panel-header"><Database size={20} color="#00f3ff" /><h3 className="panel-title">INPUT CONFIGURATION</h3></div>
      <div className="tab-group" style={{ position: 'relative', zIndex: 51 }}>
        <button className={`tab-btn ${activeTab === 'manual' ? 'active' : ''}`} onClick={() => setActiveTab('manual')}>MANUAL</button>
        <button className={`tab-btn ${activeTab === 'auto' ? 'active' : ''}`} onClick={() => setActiveTab('auto')}>AUTO DB</button>
        <button className={`tab-btn ${activeTab === 'upload' ? 'active' : ''}`} onClick={() => setActiveTab('upload')}>UPLOAD</button>
      </div>

      <div className="input-group">
        <label className="input-label">TARGET PROTEIN (PDB ID)</label>
        <div className="input-wrapper"><Search size={16} color="#888" className="input-icon" /><input className="cyber-input" placeholder="Ex: 6LU7" value={target} onChange={(e) => setTarget(e.target.value)} /></div>
        <div className="suggestions-box" style={{ pointerEvents: 'auto' }}>
          <span>Try:</span><span className="suggestion-text" onClick={() => setTarget('6LU7')}>Covid-19</span><span className="suggestion-text" onClick={() => setTarget('3PP0')}>Cancer</span><span className="suggestion-text" onClick={() => setTarget('1G0H')}>Diabetes</span>
        </div>
      </div>

      {activeTab === 'manual' && (
        <div className="input-group">
          <label className="input-label">LIGAND STRUCTURE (SMILES)</label>
          <textarea className="cyber-input textarea" rows="4" placeholder="Enter Chemical SMILES..." value={smiles} onChange={(e) => setSmiles(e.target.value)} />
          <div className="suggestions-box" style={{ pointerEvents: 'auto' }}>
            <span>Try:</span><span className="suggestion-text" onClick={() => setSmiles('CC(=O)Nc1ccc(O)cc1')}>Panadol</span><span className="suggestion-text" onClick={() => setSmiles('CC(=O)Oc1ccccc1C(=O)O')}>Aspirin</span><span className="suggestion-text" onClick={() => setSmiles('CC(C)Cc1ccc(C(C)C(=O)O)cc1')}>Ibuprofen</span>
          </div>
        </div>
      )}

      {activeTab === 'auto' && (
        <div style={{ padding: '15px', background: 'rgba(255,255,255,0.03)', borderRadius: '8px', border: '1px solid rgba(255,255,255,0.1)', marginBottom: '15px' }}>
          <div style={{ color: '#00f3ff', fontSize: '12px', fontWeight: 'bold', marginBottom: '10px' }}>DATABASE FILTERS</div>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '8px' }}>
            <div><label className="input-label" style={{ fontSize: '9px' }}>MIN MOL WEIGHT</label><input className="cyber-input" style={{ height: '30px', fontSize: '11px' }} placeholder="100" /></div>
            <div><label className="input-label" style={{ fontSize: '9px' }}>MAX MOL WEIGHT</label><input className="cyber-input" style={{ height: '30px', fontSize: '11px' }} placeholder="600" /></div>
            <div><label className="input-label" style={{ fontSize: '9px' }}>MAX LOGP</label><input className="cyber-input" style={{ height: '30px', fontSize: '11px' }} placeholder="5.0" /></div>
            <div><label className="input-label" style={{ fontSize: '9px' }}>MIN H-BOND</label><input className="cyber-input" style={{ height: '30px', fontSize: '11px' }} placeholder="0" /></div>
          </div>
          <div style={{ fontSize: '10px', color: '#666', marginTop: '10px' }}>Searching ChEMBL & PubChem libraries...</div>
        </div>
      )}

      {activeTab === 'upload' && (
        <div className="input-group" onClick={() => fileInputRef.current.click()} style={{ textAlign: 'center', border: '1px dashed rgba(255,255,255,0.2)', padding: '30px', borderRadius: '10px', color: '#888', marginBottom: '20px', cursor: 'pointer', transition: '0.3s', pointerEvents: 'auto' }} onMouseOver={(e) => e.currentTarget.style.borderColor = '#00f3ff'} onMouseOut={(e) => e.currentTarget.style.borderColor = 'rgba(255,255,255,0.2)'}>
          <input type="file" ref={fileInputRef} style={{ display: 'none' }} accept=".csv, .xlsx" onChange={handleFileSelect} />
          <Upload size={30} style={{ marginBottom: '10px', opacity: 0.5, color: selectedFile ? '#00f3ff' : 'inherit' }} />
          <div style={{ fontSize: '12px', fontWeight: 'bold', color: selectedFile ? '#fff' : '#888' }}>{selectedFile ? selectedFile.name : "Click to Upload File"}</div>
          <div style={{ fontSize: '10px', marginTop: '5px', color: '#555' }}>{selectedFile ? `${(selectedFile.size / 1024).toFixed(1)} KB` : "Supports CSV & Excel"}</div>
        </div>
      )}

      <button className="cyber-btn" onClick={handleScan} disabled={loading} style={{ pointerEvents: 'auto', zIndex: 52 }}>
        <div className="btn-content">{loading ? <span className="animate-spin" style={{ display: 'inline-block' }}><Zap size={20} /></span> : <Zap size={20} />}{loading ? "PROCESSING..." : (activeTab === 'manual' ? "INITIATE ANALYSIS" : (activeTab === 'auto' ? "SEARCH DATABASE" : "PROCESS BATCH"))}</div>
      </button>
    </div>
  );
}
