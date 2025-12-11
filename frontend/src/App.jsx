import React, { useState, useRef, useEffect, useCallback } from 'react';
import axios from 'axios';
import { Hexagon, Database, Search, Zap, Atom, Dna, Activity, ShieldCheck, Brain, Share2, Globe, Home, Info, Settings, Upload, List, Download, ArrowUpDown, Maximize, X, RefreshCw } from 'lucide-react';
import { ResponsiveContainer, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar as RechartsRadar, Tooltip } from 'recharts';
import jsPDF from 'jspdf';
import { toPng } from 'html-to-image';
import * as $3Dmol from '3dmol/build/3Dmol.js';

function App() {
  const [activeTab, setActiveTab] = useState('manual');
  const [target, setTarget] = useState('');
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [batchResults, setBatchResults] = useState([]);
  const [selectedId, setSelectedId] = useState(null);
  const [showAbout, setShowAbout] = useState(false);

  // VIEW & VISUALIZATION STATES
  const [resultView, setResultView] = useState('visualization');
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [vizMode, setVizMode] = useState('3D');
  const [style3D, setStyle3D] = useState('cartoon');
  const [viewLoading, setViewLoading] = useState(false); // New state for view transitions

  // REFS
  const cardRef = useRef(null);
  const fileInputRef = useRef(null);
  const mainViewerRef = useRef(null);
  const fullViewerRef = useRef(null);
  const viewerRef = useRef(null);  // Store the viewer instance for style updates
  const [selectedFile, setSelectedFile] = useState(null);

  // --- 3D VIEWER LOGIC ---
  const render3D = useCallback((element, bgColor) => {
    if (!element) return;

    setViewLoading(true); // Start loading

    // Clear previous viewer if it exists to ensure clean state
    element.innerHTML = '';

    const config = { backgroundColor: bgColor };
    const viewer = $3Dmol.createViewer(element, config);
    const pdbId = target || '6LU7';
    const ligandSmiles = smiles;

    // Store viewer reference for style updates
    if (element === mainViewerRef.current) {
      viewerRef.current = viewer;
    }

    // Safety timeout: stop loading if it takes too long (e.g. network error)
    const safetyTimeout = setTimeout(() => setViewLoading(false), 5000);

    $3Dmol.download(`pdb:${pdbId}`, viewer, { doAssembly: true, noSecondaryStructure: false }, function () {
      clearTimeout(safetyTimeout); // Clear timeout on success

      try {
        // 1. CLEAR EXISTING STYLES
        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();

        // 2. APPLY SELECTED STYLE TO PROTEIN (model 0)
        if (style3D === 'stick') {
          viewer.setStyle({ model: 0 }, { stick: { radius: 0.15, colorscheme: 'Jmol' } });
        } else if (style3D === 'surface') {
          viewer.setStyle({ model: 0 }, { line: { hidden: true } });
          viewer.addSurface($3Dmol.VDW, { opacity: 0.85, color: 'spectrum' }, { model: 0 });
        } else {
          // Default: CARTOON
          viewer.setStyle({ model: 0 }, { cartoon: { color: 'spectrum', thickness: 1.2, opacity: 1.0 } });
        }

        // 3. ADD LIGAND (Drug) as model 1
        if (ligandSmiles) {
          const mol = $3Dmol.read(ligandSmiles, { format: 'smi' });
          viewer.addModel(mol, "mol");
          // Always show Ligand as distinct Green Sticks
          viewer.setStyle({ model: 1 }, { stick: { colorscheme: 'greenCarbon', radius: 0.5 } });
        }

        viewer.zoomTo();
        viewer.spin('y', 0.5);
        viewer.render();
      } catch (e) {
        console.error("3D render error:", e);
      } finally {
        setViewLoading(false); // Stop loading properly
      }
    });
  }, [target, smiles, style3D]);

  // Function to update 3D style without reloading protein
  const updateStyle3D = useCallback((viewer) => {
    if (!viewer) return;

    setViewLoading(true);

    // Timeout allows UI to render loader
    const styleTimeout = setTimeout(() => {
      try {
        // Check if viewer is still valid in DOM
        if (!viewerRef.current) return;

        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();

        if (style3D === 'stick') {
          viewer.setStyle({ model: 0 }, { stick: { radius: 0.15, colorscheme: 'Jmol' } });
        } else if (style3D === 'surface') {
          viewer.setStyle({ model: 0 }, { line: { hidden: true } });
          viewer.addSurface($3Dmol.VDW, { opacity: 0.85, color: 'spectrum' }, { model: 0 });
        } else {
          viewer.setStyle({ model: 0 }, { cartoon: { color: 'spectrum', thickness: 1.2, opacity: 1.0 } });
        }

        if (smiles) {
          viewer.setStyle({ model: 1 }, { stick: { colorscheme: 'greenCarbon', radius: 0.5 } });
        }

        viewer.render();
      } catch (e) {
        console.error("Style error", e);
      } finally {
        setViewLoading(false);
      }
    }, 50);

    return () => clearTimeout(styleTimeout);
  }, [style3D, smiles]);

  // --- EFFECTS (UPDATED) ---
  useEffect(() => {
    let renderTimeout;

    // 3D Rendering Logic
    if (resultView === 'visualization' && vizMode === '3D' && result && !isFullscreen) {
      viewerRef.current = null;
      if (mainViewerRef.current) mainViewerRef.current.innerHTML = '';

      setViewLoading(true);

      renderTimeout = setTimeout(() => {
        if (mainViewerRef.current) {
          render3D(mainViewerRef.current, 'transparent');
        }
      }, 100);
    }

    // 2D Rendering Logic (New Fix: Ensure loading state is set)
    if (resultView === 'visualization' && vizMode === '2D') {
        setViewLoading(true);
    }

    return () => {
      if (renderTimeout) clearTimeout(renderTimeout);
      setViewLoading(false); // Ensure loading stops if unmounted
    };
  }, [resultView, vizMode, result, isFullscreen, render3D]);

  // Separate effect for styles
  useEffect(() => {
    let cleanupFn;
    if (viewerRef.current && vizMode === '3D' && !isFullscreen) {
      cleanupFn = updateStyle3D(viewerRef.current);
    }
    return () => {
      if (cleanupFn && typeof cleanupFn === 'function') cleanupFn();
    };
  }, [style3D, updateStyle3D, vizMode, isFullscreen]);

  // Added vizMode to deps to avoid missing dependency warning
  useEffect(() => {
    if (isFullscreen && vizMode === '3D') {
      setTimeout(() => render3D(fullViewerRef.current, 'black'), 100);
    }
  }, [isFullscreen, vizMode, render3D]);


  // --- DATA GENERATORS ---
  const generateGraphData = (score) => {
    return [
      { subject: 'Binding Affinity', A: Math.floor(score * 10), fullMark: 100, val: (Math.random() * 2 + 7).toFixed(1) },
      { subject: 'LogP', A: Math.floor(Math.random() * 60) + 40, fullMark: 100, val: (Math.random() * 4 + 1).toFixed(2) },
      { subject: 'Mol.Weight', A: Math.floor(Math.random() * 50) + 40, fullMark: 100, val: Math.floor(Math.random() * 300 + 200) },
      { subject: 'H-Bond Donors', A: Math.floor(Math.random() * 40) + 60, fullMark: 100, val: Math.floor(Math.random() * 4) + 1 },
      { subject: 'TPSA', A: Math.floor(Math.random() * 70) + 30, fullMark: 100, val: Math.floor(Math.random() * 100 + 40) },
    ];
  };

  const getMetricsTableData = () => {
    if (!result || !result.graph) return [];
    const details = {
      'Binding Affinity': { unit: 'k/mol', range: '> 7.0', desc: 'Target Interaction Strength' },
      'LogP': { unit: '', range: '1.0 - 5.0', desc: 'Lipophilicity (Cell Absorption)' },
      'Mol.Weight': { unit: 'Da', range: '< 500', desc: 'Molecular Size (Cell Permeability)' },
      'H-Bond Donors': { unit: 'Count', range: '< 5', desc: 'Lipinski Rule of 5 (Safety)' },
      'TPSA': { unit: 'Å²', range: '< 140', desc: 'Polar Surface Area (Absorption)' },
    };
    return result.graph.map(item => {
      let status = 'OPTIMAL';
      let color = '#00f3ff';
      let displayValue = item.val || item.A;
      // Note: item.val may be string if toFixed used earlier, convert for numeric checks
      const numericVal = Number(item.val) || Number(item.A) || 0;
      if (item.subject === 'Mol.Weight' && numericVal > 500) { status = 'HEAVY'; color = '#ff0055'; }
      if (item.subject === 'LogP' && numericVal > 5) { status = 'HIGH'; color = '#ff0055'; }
      if (item.subject === 'H-Bond Donors' && numericVal > 5) { status = 'FAIL'; color = '#ff0055'; }
      if (item.subject === 'TPSA' && numericVal > 140) { status = 'POOR'; color = '#ffaa00'; }
      return {
        name: item.subject,
        value: displayValue,
        unit: details[item.subject]?.unit || '',
        range: details[item.subject]?.range || '-',
        desc: details[item.subject]?.desc || '',
        status: status,
        color: color,
        barWidth: item.A
      };
    });
  };

  const handleDownload = async () => {
    if (!cardRef.current) return;
    const btn = cardRef.current.querySelector('.download-btn');
    if (btn) btn.style.display = 'none';
    try {
      const dataUrl = await toPng(cardRef.current, { quality: 1.0, backgroundColor: '#050508', style: { transform: 'none' }, filter: (node) => node.tagName !== 'BUTTON' });
      const pdf = new jsPDF('l', 'mm', 'a4');
      const pdfWidth = pdf.internal.pageSize.getWidth();
      const printHeight = (pdfWidth * 0.6);
      pdf.addImage(dataUrl, 'PNG', 0, 10, pdfWidth, printHeight);
      pdf.save(`BioGraph_Report.pdf`);
    } catch (error) {
      // use error so linter doesn't complain about unused variable
      console.error("Error generating report:", error);
      alert("Error generating report. Check console for details.");
    } finally {
      if (btn) btn.style.display = 'flex';
    }
  };

  const handleDrugClick = (drug) => {
    setSelectedId(drug.name);
    const calculatedConfidence = Math.min(Math.floor(drug.score * 10) + 2, 99);
    setResult({
      score: drug.score,
      status: drug.status,
      confidence: `${calculatedConfidence}%`,
      color: drug.status === 'ACTIVE' ? '#00f3ff' : '#ff0055',
      graph: generateGraphData(drug.score),
      name: drug.name,
      smiles: drug.smiles || smiles  // Use drug's SMILES or fallback to current SMILES
    });

    // Update the smiles state for 2D/3D viewers
    if (drug.smiles) {
      setSmiles(drug.smiles);
    }

    setResultView('visualization');
  };

  const handleTabChange = (tabName) => {
    setActiveTab(tabName); setResult(null); setBatchResults([]); setSelectedId(null); setResultView('visualization');
  };

  const handleFileSelect = (event) => {
    if (event.target.files && event.target.files[0]) { setSelectedFile(event.target.files[0]); }
  };

  const handleScan = async () => {
    // Validation
    if (!target) {
      alert("⚠️ Please enter a Target Protein ID (PDB)!");
      return;
    }
    if (activeTab === 'manual' && !smiles) {
      alert("⚠️ Please enter a valid SMILES structure!");
      return;
    }
    if (activeTab === 'upload' && !selectedFile) {
      alert("⚠️ Please select a file to upload!");
      return;
    }

    // Reset state
    setLoading(true);
    setResult(null);
    setBatchResults([]);
    setSelectedId(null);
    setResultView('visualization');

    const API_BASE = 'http://127.0.0.1:8000';

    try {
      if (activeTab === 'manual') {
        // MANUAL MODE - Single Analysis
        const response = await axios.post(`${API_BASE}/analyze`, {
          target_id: target,
          smiles: smiles,
          mode: 'manual'
        });

        if (response.data.error) {
          alert(`❌ Error: ${response.data.error}`);
          setLoading(false);
          return;
        }

        // Set result with all data from backend
        setResult({
          score: response.data.score,
          status: response.data.status,
          confidence: response.data.confidence,
          color: response.data.color,
          graph: response.data.graph,
          name: 'Manual Ligand',
          smiles: smiles  // Store SMILES for 3D/2D viewers
        });

      } else if (activeTab === 'auto') {
        // AUTO MODE - Database Search
        const response = await axios.post(`${API_BASE}/analyze`, {
          target_id: target,
          mode: 'auto'
        });

        if (response.data.error) {
          alert(`❌ Error: ${response.data.error}`);
          setLoading(false);
          return;
        }

        // Store batch results with SMILES data
        setBatchResults(response.data.results);

      } else if (activeTab === 'upload') {
        // UPLOAD MODE - File Processing
        const formData = new FormData();
        formData.append('file', selectedFile);
        formData.append('target_id', target);

        const response = await axios.post(`${API_BASE}/upload`, formData, {
          headers: { 'Content-Type': 'multipart/form-data' }
        });

        if (response.data.error) {
          alert(`❌ Error: ${response.data.error}`);
          setLoading(false);
          return;
        }

        // Store batch results
        setBatchResults(response.data.results);
      }

    } catch (error) {
      console.error('API Error:', error);

      // Better error messages based on error type
      if (error.code === 'ERR_NETWORK' || error.message.includes('Network Error')) {
        alert('🔌 Backend server is not running!\n\nPlease start the backend:\ncd backend\nuvicorn main:app --reload --port 8000');
      } else if (error.response?.status === 404) {
        alert('❌ API endpoint not found. Please check backend configuration.');
      } else if (error.response?.data?.error) {
        alert(`❌ ${error.response.data.error}`);
      } else {
        alert(`❌ Unexpected error: ${error.message}`);
      }
    } finally {
      setLoading(false);
    }
  };

  // prepare metrics for infotable view using getMetricsTableData
  const metricsForTable = getMetricsTableData();

  return (
    <div className="app-container">
      <div className="bg-grid" style={{ pointerEvents: 'none' }}></div>

      {isFullscreen && (
        <div className="viewer-container-fullscreen">
          <button className="fullscreen-close-btn" onClick={() => setIsFullscreen(false)}><X size={20} style={{ marginRight: '5px' }} /> EXIT STUDIO</button>
          <div ref={fullViewerRef} style={{ width: '100%', height: '100%' }}></div>
        </div>
      )}

      <div className="floating-elements" style={{ pointerEvents: 'none' }}>
        <div className="float-icon" style={{ top: '15%', left: '10%', animationDuration: '25s', color: '#00f3ff' }}><Atom size={180} strokeWidth={1} /></div>
        <div className="float-icon" style={{ top: '60%', right: '15%', animationDuration: '35s', color: '#bc13fe' }}><Dna size={200} strokeWidth={1} /></div>
      </div>

      <nav className="glass-nav">
        <div className="brand-identity" onClick={() => setShowAbout(false)}>
          <div className="glass-logo-container"><Hexagon size={28} color="#00f3ff" fill="rgba(0, 243, 255, 0.3)" strokeWidth={2} /></div>
          <div className="brand-text">BioGraph <span style={{ color: '#00f3ff' }}>AI</span></div>
        </div>
        <div className="nav-right">
          <div className={`nav-link ${!showAbout ? 'active-btn' : ''}`} onClick={() => setShowAbout(false)}><Home size={18} /><span>Home</span></div>
          <div className={`nav-link ${showAbout ? 'active-btn' : ''}`} onClick={() => setShowAbout(!showAbout)}><Info size={18} /><span>About</span></div>
          <div className="nav-link"><Settings size={18} /><span>Config</span></div>
        </div>
      </nav>

      <div className={`slider-container ${showAbout ? 'slide-active' : ''}`}>
        <div className="page-section" style={{ position: 'relative' }}>
          <div className="main-layout">

            <div className="glass-panel panel-left" style={{ zIndex: 50 }}>
              <div className="panel-header"><Database size={20} color="#00f3ff" /><h3 className="panel-title">INPUT CONFIGURATION</h3></div>
              <div className="tab-group" style={{ position: 'relative', zIndex: 51 }}>
                <button className={`tab-btn ${activeTab === 'manual' ? 'active' : ''}`} onClick={() => handleTabChange('manual')}>MANUAL</button>
                <button className={`tab-btn ${activeTab === 'auto' ? 'active' : ''}`} onClick={() => handleTabChange('auto')}>AUTO DB</button>
                <button className={`tab-btn ${activeTab === 'upload' ? 'active' : ''}`} onClick={() => handleTabChange('upload')}>UPLOAD</button>
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
                <div className="btn-content">{loading ? <Activity className="animate-spin" size={20} /> : (activeTab === 'manual' ? <Zap size={20} /> : <Database size={20} />)}{loading ? "PROCESSING..." : (activeTab === 'manual' ? "INITIATE ANALYSIS" : (activeTab === 'auto' ? "SEARCH DATABASE" : "PROCESS BATCH"))}</div>
              </button>
            </div>

            <div className="glass-panel panel-right" style={{ zIndex: 50 }}>
              {/* COMBINED TOP HEADER - System Status + Navigation Tabs */}
              <div style={{
                height: '60px',
                background: 'rgba(0,0,0,0.3)',
                borderBottom: '1px solid rgba(255,255,255,0.1)',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                padding: '0 20px',
                gap: '20px'
              }}>
                {/* LEFT: System Status Badge */}
                <div className="header-badge" style={{ margin: 0 }}>
                  <Activity size={14} color="#00f3ff" className={loading ? "animate-pulse" : ""} />
                  <span className="badge-text">SYSTEM: <span style={{ color: loading ? '#00f3ff' : '#fff' }}>{loading ? "BUSY" : (result || batchResults.length > 0 ? "ONLINE" : "IDLE")}</span></span>
                </div>

                {/* RIGHT: Tabs (only show when result exists) */}
                {result && (
                  <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                    {activeTab !== 'manual' && (
                      <button onClick={() => setResult(null)} className="nav-link" style={{ marginRight: '10px', display: 'flex', alignItems: 'center', gap: '5px' }}>
                        <List size={16} /> Back to List
                      </button>
                    )}
                    {[
                      { label: 'Studio', id: 'visualization' },
                      { label: 'Radar', id: 'radar' },
                      { label: 'InfoTable', id: 'infotable' }
                    ].map(tab => (
                      <button
                        key={tab.id}
                        onClick={() => setResultView(tab.id)}
                        className={`nav-link ${resultView === tab.id ? 'active-btn' : ''}`}
                      >
                        {tab.label}
                      </button>
                    ))}
                  </div>
                )}
              </div>

              <div className="hologram-wrapper">
                {(loading || (!result && batchResults.length === 0)) ? (
                  <div className="hologram-inner" style={{ pointerEvents: 'none' }}>
                    <div className="dna-spinner"><Dna size={120} color="#00f3ff" strokeWidth={1.5} /></div>
                    <div className="ring-1"></div><div className="ring-2"></div><div className="core-glow"></div>
                    {loading && <div className="loading-text">{activeTab === 'auto' ? "SCANNING DATABASE..." : (activeTab === 'upload' ? "READING FILE..." : "LOADING ANALYTICS...")}</div>}
                  </div>
                ) : (
                  result ? (
                    <div style={{ width: '100%', height: '100%', position: 'relative', display: 'flex', flexDirection: 'column' }}>

                      {/* RIGHT PANEL - TABBED INTERFACE */}
                      <div className="pro-studio-container" style={{ display: 'flex', flexDirection: 'column', width: '100%', height: '100%', overflow: 'hidden', background: 'rgba(5, 5, 8, 0.6)', borderRadius: '12px', border: '1px solid rgba(255,255,255,0.1)' }}>

                        {/* MAIN CONTENT AREA */}
                        <div style={{ flexGrow: 1, position: 'relative', overflow: 'hidden', width: '100%', height: '100%' }}>

                          {/* --- STUDIO VIEW --- */}
                          {resultView === 'visualization' && (
                            <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
                              {/* Studio Controls */}
                              <div style={{ padding: '10px 20px', borderBottom: '1px solid rgba(255,255,255,0.05)', display: 'flex', justifyContent: 'space-between', alignItems: 'center', background: 'rgba(255,255,255,0.02)' }}>
                                <div style={{ display: 'flex', gap: '5px' }}>
                                  <button onClick={() => setVizMode('2D')} className={`nav-link ${vizMode === '2D' ? 'active-btn' : ''}`}>2D View</button>
                                  <button onClick={() => setVizMode('3D')} className={`nav-link ${vizMode === '3D' ? 'active-btn' : ''}`}>3D View</button>
                                </div>
                                {vizMode === '3D' && (
                                  <div style={{ display: 'flex', gap: '5px' }}>
                                    {['Cartoon', 'Stick', 'Surface'].map(s => (
                                      <button key={s} onClick={() => setStyle3D(s.toLowerCase())} className={`nav-link ${style3D === s.toLowerCase() ? 'active-btn' : ''}`}>{s}</button>
                                    ))}
                                    <button onClick={() => render3D(mainViewerRef.current, 'transparent')} className="nav-link"><RefreshCw size={14} /></button>
                                  </div>
                                )}
                              </div>

                              {/* Viewport (UPDATED) */}
                              <div key={vizMode} style={{ flexGrow: 1, position: 'relative', background: vizMode === '3D' ? 'radial-gradient(circle at center, rgba(30,30,40,0.5) 0%, rgba(5,5,8,0.8) 100%)' : '#fff', overflow: 'hidden' }}>

                                {/* LOADING OVERLAY */}
                                {viewLoading && (
                                  <div style={{ position: 'absolute', top: 0, left: 0, width: '100%', height: '100%', background: 'rgba(0,0,0,0.7)', zIndex: 100, display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center', backdropFilter: 'blur(2px)' }}>
                                    <div className="animate-spin" style={{ marginBottom: '10px' }}><Atom size={40} color="#00f3ff" /></div>
                                    <div style={{ color: '#00f3ff', fontSize: '12px', letterSpacing: '2px', fontWeight: 'bold' }}>LOADING STRUCTURE...</div>
                                  </div>
                                )}

                                {vizMode === '2D' ? (
                                  <div style={{ width: '100%', height: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', padding: '20px' }}>
                                    <img
                                      key={smiles}
                                      src={`http://127.0.0.1:8000/get_image?smiles=${encodeURIComponent(smiles)}`}
                                      alt="Molecular Structure"
                                      onLoadStart={() => setViewLoading(true)}
                                      onLoad={() => setViewLoading(false)}
                                      onError={() => setViewLoading(false)}
                                      style={{
                                        maxWidth: '90%',
                                        maxHeight: '90%',
                                        width: 'auto',
                                        height: 'auto',
                                        objectFit: 'contain',
                                        display: 'block',
                                        filter: 'drop-shadow(0 0 10px rgba(0,0,0,0.2))'
                                      }}
                                    />
                                  </div>
                                ) : (
                                  <div ref={mainViewerRef} style={{ width: '100%', height: '100%' }}></div>
                                )}
                              </div>
                            </div>
                          )}

                          {/* --- RADAR VIEW --- */}
                          {(resultView === 'radar') && (
                            <div style={{ width: '100%', height: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', padding: '20px' }}>
                              <ResponsiveContainer width="100%" height="100%">
                                <RadarChart cx="50%" cy="50%" outerRadius="80%" data={result.graph || [
                                  { subject: 'Mol. Weight', A: 120, fullMark: 150 },
                                  { subject: 'LogP', A: 98, fullMark: 150 },
                                  { subject: 'H-Bond Donors', A: 86, fullMark: 150 },
                                  { subject: 'H-Bond Acceptors', A: 99, fullMark: 150 },
                                  { subject: 'TPSA', A: 85, fullMark: 150 },
                                  { subject: 'Rotatable Bonds', A: 65, fullMark: 150 },
                                ]}>
                                  <PolarGrid stroke="rgba(255,255,255,0.2)" />
                                  <PolarAngleAxis dataKey="subject" tick={{ fill: 'white', fontSize: 12 }} />
                                  <PolarRadiusAxis angle={30} domain={[0, 150]} tick={false} axisLine={false} />
                                  <RechartsRadar name="Molecule" dataKey="A" stroke="#00f3ff" strokeWidth={3} fill="#00f3ff" fillOpacity={0.3} />
                                  <Tooltip contentStyle={{ backgroundColor: '#000', border: '1px solid #333', borderRadius: '8px' }} />
                                </RadarChart>
                              </ResponsiveContainer>
                            </div>
                          )}

                          {/* --- INFOTABLE VIEW --- */}
                          {(resultView === 'infotable') && (
                            <div style={{ width: '100%', height: '100%', overflowY: 'auto', padding: '20px' }}>
                              <table style={{ width: '100%', borderCollapse: 'collapse', color: '#fff' }}>
                                <thead>
                                  <tr style={{ borderBottom: '1px solid #444', textAlign: 'left' }}>
                                    <th style={{ padding: '10px', color: '#888' }}>PROPERTY</th>
                                    <th style={{ padding: '10px', color: '#888' }}>VALUE</th>
                                    <th style={{ padding: '10px', color: '#888' }}>STATUS</th>
                                  </tr>
                                </thead>
                                <tbody>
                                  {metricsForTable.length > 0 ? metricsForTable.map((row, i) => (
                                    <tr key={i} style={{ borderBottom: '1px solid rgba(255,255,255,0.05)' }}>
                                      <td style={{ padding: '12px', fontWeight: '600' }}>{row.name}</td>
                                      <td style={{ padding: '12px', fontFamily: 'monospace', color: '#00f3ff' }}>{`${row.value} ${row.unit}`}</td>
                                      <td style={{ padding: '12px' }}><span style={{ background: row.status === 'OPTIMAL' ? 'rgba(0,255,100,0.06)' : 'rgba(255,0,85,0.06)', color: row.status === 'OPTIMAL' ? '#0f0' : '#ff6b9a', padding: '2px 8px', borderRadius: '12px', fontSize: '10px' }}>{row.status}</span></td>
                                    </tr>
                                  )) : (
                                    // fallback static rows if no metrics
                                    [
                                      { k: 'Molecular Weight', v: '428.5 g/mol', s: 'Optimal' },
                                      { k: 'Generic Name', v: result.name || 'N/A', s: '-' },
                                      { k: 'Target ID', v: target || 'N/A', s: 'Active' },
                                      { k: 'LogP', v: '2.4', s: 'Lipophilic' },
                                      { k: 'H-Bond Donors', v: '2', s: 'Good' },
                                      { k: 'H-Bond Acceptors', v: '4', s: 'Good' },
                                      { k: 'Ro5 Violations', v: '0', s: 'Pass' },
                                      { k: 'TPSA', v: '85.4', s: 'Moderate' },
                                      { k: 'Rotatable Bonds', v: '5', s: 'Flexible' },
                                    ].map((row, i) => (
                                      <tr key={i} style={{ borderBottom: '1px solid rgba(255,255,255,0.05)' }}>
                                        <td style={{ padding: '12px', fontWeight: '600' }}>{row.k}</td>
                                        <td style={{ padding: '12px', fontFamily: 'monospace', color: '#00f3ff' }}>{row.v}</td>
                                        <td style={{ padding: '12px' }}><span style={{ background: 'rgba(0,255,100,0.2)', color: '#0f0', padding: '2px 8px', borderRadius: '12px', fontSize: '10px' }}>{row.s}</span></td>
                                      </tr>
                                    ))
                                  )}
                                </tbody>
                              </table>
                            </div>
                          )}
                        </div>
                      </div>
                    </div>
                  ) : (
                    // BATCH LIST VIEW
                    <div className="scan-results-list" style={{ position: 'relative', zIndex: 20, width: '100%', marginTop: '20px', height: '350px', overflowY: 'auto', overflowX: 'hidden', paddingRight: '5px', background: 'rgba(0,0,0,0.2)', borderRadius: '15px', pointerEvents: 'auto' }}>
                      <div className="list-header" style={{ position: 'sticky', top: 0, zIndex: 30, background: 'rgba(5, 5, 10, 0.95)', backdropFilter: 'blur(10px)', borderBottom: '1px solid #00f3ff', padding: '15px 20px', display: 'flex', justifyContent: 'space-between', alignItems: 'center', boxShadow: '0 5px 20px rgba(0,0,0,0.5)', color: '#00f3ff', fontSize: '13px', fontWeight: 'bold', letterSpacing: '1px' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}><List size={14} /> DRUG NAME</div>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>SCORE <ArrowUpDown size={12} /></div>
                      </div>
                      <div style={{ padding: '10px' }}>
                        {batchResults.map((item, index) => (
                          <div key={index} onClick={() => handleDrugClick(item)} style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', padding: '15px 20px', marginBottom: '8px', background: selectedId === item.name ? 'rgba(0,243,255,0.1)' : 'rgba(255,255,255,0.03)', border: selectedId === item.name ? '1px solid #00f3ff' : '1px solid transparent', borderRadius: '10px', cursor: 'pointer', transition: 'all 0.2s' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                              <div style={{ width: '8px', height: '8px', borderRadius: '50%', background: item.status === 'ACTIVE' ? '#00f3ff' : '#ff0055', boxShadow: item.status === 'ACTIVE' ? '0 0 10px #00f3ff' : 'none' }}></div>
                              <div><div style={{ fontWeight: 'bold', color: '#fff' }}>{item.name}</div><div style={{ fontSize: '10px', color: '#666' }}>Auto-generated</div></div>
                            </div>
                            <div style={{ color: item.status === 'ACTIVE' ? '#00f3ff' : '#ff0055', fontWeight: 'bold' }}>{item.score}</div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )
                )}
                <div className="grid-overlay" style={{ zIndex: 0, pointerEvents: 'none' }}></div>
              </div>
            </div>
          </div >

          {result && (
            <div className={`result-card ${result.status === 'ACTIVE' ? 'active' : 'inactive'}`} style={{ borderColor: result.color, zIndex: 60, pointerEvents: 'auto' }} ref={cardRef}>
              <div className="score-section">
                <div className="score-box"><div className="input-label">BINDING AFFINITY</div><div className="score-val fade-in-text" style={{ textShadow: `0 0 20px ${result.color}` }}>{result.score}</div></div>
                <div className="vertical-div"></div>
                <div><div className="input-label">PREDICTION</div><div className="status-val fade-in-text" style={{ color: result.color }}>{result.status}</div></div>
              </div>
              <div className="metrics-container fade-in-text">
                <div className="metric-box"><div className="input-label">TOXICITY</div><div className="metric-row" style={{ color: '#00f3ff' }}><ShieldCheck size={18} /> Low Risk</div></div>
                <div className="metric-box"><div className="input-label">CONFIDENCE</div><div className="metric-row" style={{ color: '#bc13fe' }}><Dna size={18} /> {result.confidence}</div></div>
              </div>
              <button className="download-btn" onClick={handleDownload} style={{ background: 'rgba(255,255,255,0.1)', border: '1px solid rgba(255,255,255,0.2)', color: '#fff', padding: '10px', borderRadius: '50%', cursor: 'pointer', display: 'flex', alignItems: 'center', justifyContent: 'center', transition: '0.3s' }}>
                <Download size={18} />
              </button>
            </div>
          )
          }
        </div >

        <div className="page-section" style={{ pointerEvents: 'auto' }}>
          <div className="about-page-content">
            <div className="hologram-inner" style={{ marginBottom: '30px', height: '150px' }}><div className="dna-spinner"><Hexagon size={80} color="#00f3ff" fill="rgba(0, 243, 255, 0.3)" strokeWidth={2} /></div></div>
            <div className="brand-text">BioGraph <span style={{ color: '#00f3ff' }}>AI</span></div> <div><br /></div>
            <div className="about-hero-subtitle">Next-Generation Drug Discovery</div>
            <div className="features-grid">
              <div className="feature-card"><div className="f-icon"><Brain size={30} /></div><div className="f-title">Graph Neural Networks</div><div className="f-desc">Utilizing advanced GNNs to predict molecular interactions with high accuracy.</div></div>
              <div className="feature-card"><div className="f-icon"><Share2 size={30} /></div><div className="f-title">Binding Affinity</div><div className="f-desc">Calculates the strength of interactions between protein targets and ligands.</div></div>
              <div className="feature-card"><div className="f-icon"><Globe size={30} /></div><div className="f-title">ADMET Profiling</div><div className="f-desc">Instant analysis of Toxicity, Solubility, and Drug-likeness properties.</div></div>
            </div>
          </div>
        </div>

      </div >
    </div >
  );
}

export default App;