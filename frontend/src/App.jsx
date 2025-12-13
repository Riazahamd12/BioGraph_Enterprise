/* --- IMPORTS --- */
import React, { useState, useRef, useEffect, useCallback } from 'react';
import axios from 'axios';

// 1. ICONS (Yahan 'Radar' icon hai)
import { 
  Hexagon, Home, Info, Settings, List, ArrowUpDown, X, RefreshCw, 
  Dna, Atom, Upload, Database, Search, Zap, Activity, ShieldCheck, 
  Dna as DnaIcon, Download as DownloadIcon, Brain, Magnet, FlaskConical, 
  Box, Radar, Table // <--- Ye Icon hai
} from 'lucide-react';

// 2. CHARTS (Yahan 'Radar' ko rename karke 'RechartsRadar' banaya hai)
import { 
  ResponsiveContainer, RadarChart, PolarGrid, PolarAngleAxis, 
  PolarRadiusAxis, Radar as RechartsRadar, Tooltip // <--- YEH ZAROORI HAI
} from 'recharts';

import jsPDF from 'jspdf';
import { toPng } from 'html-to-image';
import InteractiveImage from './components/InteractiveImage';
import Viewer3D from './components/Viewer3D';
import ResultCard from './components/ResultCard';
import Sidebar from './components/Sidebar';
import "./styles/index.css";


/**
 * App.jsx
 * - This file wires components together.
 * - I kept your original state names & API calls so behaviour is identical.
 * - Minor reformatting + comments only.
 */

function App() {

  // --- core app state (same as your original) ---
  const [activeTab, setActiveTab] = useState('manual');
  const [target, setTarget] = useState('');
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [batchResults, setBatchResults] = useState([]);
  const [selectedId, setSelectedId] = useState(null);
  const [showAbout, setShowAbout] = useState(false);

  // view & viz
  const [resultView, setResultView] = useState('visualization');
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [vizMode, setVizMode] = useState('3D');
  const [style3D, setStyle3D] = useState('cartoon');
  const [viewLoading, setViewLoading] = useState(false);

  // refs
  const cardRef = useRef(null);
  const fileInputRef = useRef(null);
  const mainViewerRef = useRef(null);
  const fullViewerRef = useRef(null);
  const viewerRef = useRef(null);
  const [selectedFile, setSelectedFile] = useState(null);

  // --- 3D render callback (used by Viewer3D too) ---
  const render3D = useCallback((element, bgColor) => {
    // We keep same behavior as before; Viewer3D also handles rendering
    if (!element) return;
    setViewLoading(true);
    element.innerHTML = '';
    const config = { backgroundColor: bgColor };
    const viewer = window.$3Dmol?.createViewer ? window.$3Dmol.createViewer(element, config) : null;
    if (!viewer) {
      // Fallback: if 3dmol not global, Viewer3D component will manage it
      setViewLoading(false);
      return;
    }
    if (element === mainViewerRef.current) viewerRef.current = viewer;

    const pdbId = target || '6LU7';
    const ligandSmiles = smiles;

    const safetyTimeout = setTimeout(() => setViewLoading(false), 5000);

    window.$3Dmol.download(`pdb:${pdbId}`, viewer, { doAssembly: true, noSecondaryStructure: false }, function () {
      clearTimeout(safetyTimeout);
      try {
        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();
        if (style3D === 'stick') {
          viewer.setStyle({}, { stick: { radius: 0.15, colorscheme: 'Jmol' } });
        } else if (style3D === 'surface') {
          viewer.setStyle({}, { line: { hidden: true } });
          viewer.addSurface(window.$3Dmol.VDW, { opacity: 0.85, color: 'spectrum' }, { hetflag: false });
        } else {
          viewer.setStyle({}, { cartoon: { color: 'spectrum', thickness: 1.2, opacity: 1.0 } });
        }

        if (ligandSmiles) {
          const mol = window.$3Dmol.read(ligandSmiles, { format: 'smi' });
          viewer.addModel(mol, 'mol');
          viewer.setStyle({ model: -1 }, { stick: { colorscheme: 'greenCarbon', radius: 0.5 } });
        }

        viewer.zoomTo();
        viewer.spin('y', 0.5);
        viewer.render();
      } catch (e) {
        console.error('3D error (App):', e);
      } finally {
        setViewLoading(false);
      }
    });
  }, [target, smiles, style3D]);

  // update style without reload
  const updateStyle3D = useCallback((viewer) => {
    if (!viewer) return;
    setViewLoading(true);
    const t = setTimeout(() => {
      try {
        if (!viewerRef.current) return;
        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();
        if (style3D === 'stick') {
          viewer.setStyle({}, { stick: { radius: 0.15, colorscheme: 'Jmol' } });
        } else if (style3D === 'surface') {
          viewer.setStyle({}, { line: { hidden: true } });
          viewer.addSurface(window.$3Dmol.VDW, { opacity: 0.85, color: 'spectrum' }, { hetflag: false });
        } else {
          viewer.setStyle({}, { cartoon: { color: 'spectrum', thickness: 1.2, opacity: 1.0 } });
        }
        if (smiles) viewer.setStyle({ model: -1 }, { stick: { colorscheme: 'greenCarbon', radius: 0.5 } });
        viewer.render();
      } catch (e) {
        console.error('Style error (App):', e);
      } finally {
        setViewLoading(false);
      }
    }, 50);
    return () => clearTimeout(t);
  }, [style3D, smiles]);

  // effect: render 3D when result available (same logic)
  useEffect(() => {
    let renderTimeout;
    if (resultView === 'visualization' && vizMode === '3D' && result && !isFullscreen) {
      viewerRef.current = null;
      if (mainViewerRef.current) mainViewerRef.current.innerHTML = '';
      setViewLoading(true);
      renderTimeout = setTimeout(() => {
        if (mainViewerRef.current) render3D(mainViewerRef.current, 'transparent');
      }, 100);
    }
    if (resultView === 'visualization' && vizMode === '2D') {
      setViewLoading(true);
    }
    return () => {
      if (renderTimeout) clearTimeout(renderTimeout);
      setViewLoading(false);
    };
  }, [resultView, vizMode, result, isFullscreen, render3D, smiles]);

  useEffect(() => {
    let cleanupFn;
    if (viewerRef.current && vizMode === '3D' && !isFullscreen) {
      cleanupFn = updateStyle3D(viewerRef.current);
    }
    return () => {
      if (cleanupFn && typeof cleanupFn === 'function') cleanupFn();
    };
  }, [style3D, updateStyle3D, vizMode, isFullscreen]);

  useEffect(() => {
    if (isFullscreen && vizMode === '3D') {
      setTimeout(() => render3D(fullViewerRef.current, 'black'), 100);
    }
  }, [isFullscreen, vizMode, render3D]);

  // ------------------------
  // Helper data functions - same as original
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

  // ------------------------
  // API / user handlers (kept same)
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
      smiles: drug.smiles || smiles
    });

    if (drug.smiles) setSmiles(drug.smiles);
    setResultView('visualization');
  };

  const handleTabChange = (tabName) => {
    setActiveTab(tabName); setResult(null); setBatchResults([]); setSelectedId(null); setResultView('visualization');
  };

  const handleFileSelect = (event) => {
    if (event.target.files && event.target.files[0]) { setSelectedFile(event.target.files[0]); }
  };

  const handleScan = async () => {
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

    setLoading(true);
    setResult(null);
    setBatchResults([]);
    setSelectedId(null);
    setResultView('visualization');

    const API_BASE = 'http://127.0.0.1:8000';

    try {
      if (activeTab === 'manual') {
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

        setResult({
          score: response.data.score,
          status: response.data.status,
          confidence: response.data.confidence,
          color: response.data.color,
          graph: response.data.graph,
          name: 'Manual Ligand',
          smiles: smiles
        });

      } else if (activeTab === 'auto') {
        const response = await axios.post(`${API_BASE}/analyze`, {
          target_id: target,
          mode: 'auto'
        });

        if (response.data.error) {
          alert(`❌ Error: ${response.data.error}`);
          setLoading(false);
          return;
        }

        setBatchResults(response.data.results);

      } else if (activeTab === 'upload') {
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

        setBatchResults(response.data.results);
      }

    } catch (error) {
      console.error('API Error:', error);
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

  const metricsForTable = getMetricsTableData();

  // ------------------------
  // Render (unchanged visuals, but using our small components)
  return (
    <div className="app-container">
      <div className="bg-grid" style={{ pointerEvents: 'none' }}></div>

      {isFullscreen && (
        <div className="viewer-container-fullscreen">
          <button className="fullscreen-close-btn" onClick={() => setIsFullscreen(false)}><X size={20} style={{ marginRight: '5px' }} /> EXIT STUDIO</button>
          <div ref={fullViewerRef} style={{ width: '100%', height: '100%' }} />
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
        {/* --- NEW: CENTER TEXT ADDED HERE --- */}
        <div className="nav-center-text">Next-Generation Drug Discovery</div>
        <div className="nav-right">
          <div className={`nav-link ${!showAbout ? 'active-btn' : ''}`} onClick={() => setShowAbout(false)}><Home size={18} /><span>Home</span></div>
          <div className={`nav-link ${showAbout ? 'active-btn' : ''}`} onClick={() => setShowAbout(!showAbout)}><Info size={18} /><span>About</span></div>
          <div className="nav-link"><Settings size={18} /><span>Settings</span></div>
        </div>
      </nav>

      <div className={`slider-container ${showAbout ? 'slide-active' : ''}`}>
        <div className="page-section" style={{ position: 'relative' }}>
          <div className="main-layout">
            {/* Left sidebar (extracted component) */}
            <Sidebar
              activeTab={activeTab}
              setActiveTab={handleTabChange}
              target={target}
              setTarget={setTarget}
              smiles={smiles}
              setSmiles={setSmiles}
              selectedFile={selectedFile}
              fileInputRef={fileInputRef}
              handleFileSelect={handleFileSelect}
              handleScan={handleScan}
              loading={loading}
            />

            {/* Right panel (main visualization + lists) */}
            <div className="glass-panel panel-right" style={{ zIndex: 50 }}>
              <div style={{
                height: '60px',
                background: 'rgba(0,0,0,0)',
                borderBottom: '1px solid rgba(255,255,255,0.1)',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                padding: '0 20px',
                gap: '20px'
              }}>
                <div className="header-badge" style={{ margin: 0, padding: '12px 15px' }}>
                  
                  {/* --- 1. DESKTOP VIEW (Full Details) --- */}
                  <div className="status-desktop">
                    <Activity size={14} color="#00f3ff" className={loading ? "animate-pulse" : ""} />
                    <span className="badge-text" style={{ marginLeft: '8px' }}>
                      SYSTEM: <span style={{ color: loading ? '#00f3ff' : '#fff' }}>
                        {loading ? "BUSY" : (result || batchResults.length > 0 ? "OPERATION COMPLETED SUCCESSFULLY" : "IDLE")}
                      </span>
                    </span>
                  </div>

                  {/* --- 2. MOBILE VIEW (Only "SYSTEM" + Lights) --- */}
                  <div className="status-mobile">
                    {/* Yahan sirf 'SYSTEM' likha hai, agay koi status variable nahi lagaya */}
                    <span className="badge-text" style={{ fontSize: '10px' }}>SYSTEM
                    </span>
                    
                    <div className="status-lights">
                      <div className={`light red ${!loading && !result && batchResults.length === 0 ? 'active' : ''}`}></div>
                      <div className={`light blue ${loading ? 'active animate-pulse' : ''}`}></div>
                      <div className={`light green ${!loading && (result || batchResults.length > 0) ? 'active' : ''}`}></div>
                    </div>
                  </div>

                </div>

                {result && (
                  <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                    {activeTab !== 'manual' && (
                      <button onClick={() => setResult(null)} className="nav-link" style={{ marginRight: '10px', display: 'flex', alignItems: 'center', gap: '5px' }}>
                        <List size={16} /> <span className="desktop-text">Back</span>
                      </button>
                    )}
                    {[
                      { label: 'Studio', id: 'visualization', icon: <Box size={16} /> },
                      { label: 'Radar', id: 'radar', icon: <Radar size={16} /> },
                      { label: 'InfoTable', id: 'infotable', icon: <Table size={16} /> }
                    ].map(tab => (
                      <button
                        key={tab.id}
                        onClick={() => setResultView(tab.id)}
                        className={`nav-link ${resultView === tab.id ? 'active-btn' : ''}`}
                        title={tab.label} // Mobile par press hold karne par naam dikhega
                        style={{ display: 'flex', alignItems: 'center', gap: '8px', justifyContent: 'center' }}
                      >
                        {tab.icon}
                        <span className="mobile-hide-text">{tab.label}</span>
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
                      <div className="pro-studio-container" style={{ display: 'flex', flexDirection: 'column', width: '100%', height: '100%', overflow: 'hidden', background: 'rgba(5, 5, 8, 0.6)', borderRadius: '12px', border: '1px solid rgba(255,255,255,0.1)' }}>
                        <div style={{ flexGrow: 1, position: 'relative', overflow: 'hidden', width: '100%', height: '100%' }}>
                          {resultView === 'visualization' && (
                            <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
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

                              <div key={vizMode} style={{ flexGrow: 1, position: 'relative', background: vizMode === '3D' ? 'radial-gradient(circle at center, rgba(30,30,40,0.5) 0%, rgba(5,5,8,0.8) 100%)' : '#fff', overflow: 'hidden' }}>
                                {viewLoading && (
                                  <div style={{ position: 'absolute', top: 0, left: 0, width: '100%', height: '100%', background: 'rgba(0,0,0,0.7)', zIndex: 100, display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center', backdropFilter: 'blur(2px)' }}>
                                    <div className="animate-spin" style={{ marginBottom: '10px' }}><Atom size={40} color="#00f3ff" /></div>
                                    <div style={{ color: '#00f3ff', fontSize: '12px', letterSpacing: '2px', fontWeight: 'bold' }}>LOADING STRUCTURE...</div>
                                  </div>
                                )}

                                {vizMode === '2D' ? (
                                  <div style={{ width: '100%', height: '100%', padding: '20px', boxSizing: 'border-box' }}>
                                    <InteractiveImage
                                      key={smiles}
                                      src={`http://127.0.0.1:8000/get_image?smiles=${encodeURIComponent(smiles)}`}
                                      alt="Molecular Structure"
                                      onLoadStart={() => setViewLoading(true)}
                                      onLoad={() => setViewLoading(false)}
                                      onError={() => setViewLoading(false)}
                                    />
                                  </div>
                                ) : (
                                  <Viewer3D
                                    target={target || '6LU7'}
                                    smiles={smiles}
                                    style3D={style3D}
                                    bgColor={'transparent'}
                                    onLoadingChange={setViewLoading}
                                    viewerRefProp={viewerRef}
                                  />
                                )}
                              </div>
                            </div>
                          )}

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
          </div>

          {result && <ResultCard result={result} cardRef={cardRef} onDownload={handleDownload} />}
        </div>

        {/* About page */}
        <div className="page-section" style={{ pointerEvents: 'auto' }}>
          <div className="about-page-content">
            <div className="hologram-inner" style={{ marginBottom: '30px', height: '150px' }}><div className="dna-spinner"><Hexagon size={80} color="#00f3ff" fill="rgba(0, 243, 255, 0.3)" strokeWidth={2} /></div></div>
            <div className="brand-text">BioGraph <span style={{ color: '#00f3ff' }}>AI</span></div> <div><br /></div>
            <div className="about-hero-subtitle">Next-Generation Drug Discovery</div>
           <div className="features-grid">
              
              {/* Card 1: GNN (Brain Icon) */}
              <div className="feature-card">
                <div className="f-icon">
                  <Brain size={30} />
                </div>
                <div className="f-title">AI-Driven GNNs</div>
                <div className="f-desc">
                  Advanced Graph Neural Networks (GNNs) analyze molecular graphs to predict drug-target interactions with state-of-the-art accuracy.
                </div>
              </div>

              {/* Card 2: Binding (Magnet Icon) */}
              <div className="feature-card">
                <div className="f-icon">
                  <Magnet size={30} style={{ transform: 'rotate(45deg)' }} />
                </div>
                <div className="f-title">Binding Affinity</div>
                <div className="f-desc">
                  Quantifies the binding strength (Kd/Ki) between ligands and protein pockets, simulating molecular docking to find potent inhibitors.
                </div>
              </div>

              {/* Card 3: ADMET (Shield Icon) */}
              <div className="feature-card">
                <div className="f-icon">
                  <ShieldCheck size={30} />
                </div>
                <div className="f-title">ADMET & Safety</div>
                <div className="f-desc">
                  Comprehensive profiling of Absorption, Toxicity, and Drug-Likeness (QED, Lipinski Rule) to ensure clinical safety.
                </div>
              </div>

              {/* Card 4: Visualization (Atom Icon) */}
              <div className="feature-card">
                <div className="f-icon">
                  <Atom size={30} />
                </div>
                <div className="f-title">3D Visualization</div>
                <div className="f-desc">
                  Interactive 3D molecular rendering engine to inspect chemical structures, binding pockets, and confirmations in real-time.
                </div>
              </div>

              {/* Card 5: Repurposing (Refresh Icon) */}
              <div className="feature-card">
                <div className="f-icon">
                  <RefreshCw size={30} />
                </div>
                <div className="f-title">Drug Repurposing</div>
                <div className="f-desc">
                  Screening libraries of existing FDA-approved drugs to identify novel therapeutic uses, accelerating the discovery timeline.
                </div>
              </div>

              {/* Card 6: Analytics (Activity Icon) */}
              <div className="feature-card">
                <div className="f-icon">
                  <Activity size={30} />
                </div>
                <div className="f-title">Smart Analytics</div>
                <div className="f-desc">
                  Instant graphical analysis using Radar Charts to visualize multi-parameter optimization scores and molecular properties.
                </div>
              </div>

            </div>
          </div>
        </div>

      </div>
    </div>
  );
}

export default App;
