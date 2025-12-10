import React, { useState, useRef } from 'react'; 
import { Hexagon, Settings, Info, Database, Layers, Search, Zap, Atom, Dna, FlaskConical, Activity, ShieldCheck, ArrowLeft, Brain, Share2, Globe, Home, Upload, FileText, List, Microscope, Download, ArrowUpDown, ChevronRight } from 'lucide-react';
import { ResponsiveContainer, Radar, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis } from 'recharts';

function App() {
  // --- STATE ---
  const [activeTab, setActiveTab] = useState('manual'); 
  const [target, setTarget] = useState(''); 
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  
  const [result, setResult] = useState(null); 
  const [batchResults, setBatchResults] = useState([]); 
  const [selectedId, setSelectedId] = useState(null); 
  const [showAbout, setShowAbout] = useState(false); 

  // File Upload State
  const [selectedFile, setSelectedFile] = useState(null);
  const fileInputRef = useRef(null); 

  // --- HELPER: GRAPH DATA GENERATOR ---
  const generateGraphData = (score) => {
    return [
      { subject: 'Binding', A: Math.floor(score * 10), fullMark: 100 },
      { subject: 'LogP', A: Math.floor(Math.random() * 80) + 20, fullMark: 100 },
      { subject: 'Mol.Weight', A: Math.floor(Math.random() * 90), fullMark: 100 },
      { subject: 'Toxicity', A: Math.floor(Math.random() * 30), fullMark: 100 },
      { subject: 'QED', A: Math.floor(Math.random() * 100), fullMark: 100 },
    ];
  };

  // --- INTERACTION: CLICK ON LIST ITEM ---
  const handleDrugClick = (drug) => {
    setSelectedId(drug.name);
    setResult({
      score: drug.score,
      status: drug.status,
      confidence: "94.2%",
      color: drug.status === 'ACTIVE' ? '#00f3ff' : '#ff0055',
      graph: generateGraphData(drug.score)
    });
  };

  // --- INTERACTION: FILE SELECTION ---
  const handleFileSelect = (event) => {
    if (event.target.files && event.target.files[0]) {
      setSelectedFile(event.target.files[0]);
    }
  };

  // --- MAIN LOGIC ---
  const handleScan = async () => {
    if (activeTab === 'manual') {
        if (!target) { alert("Please enter Target ID!"); return; }
        if (!smiles) { alert("Please enter SMILES!"); return; }
    }
    if (activeTab === 'auto' && !target) { alert("Please enter Target ID!"); return; }
    if (activeTab === 'upload' && !selectedFile) { alert("Please select a file first!"); return; }

    setLoading(true);
    setResult(null);
    setBatchResults([]);
    setSelectedId(null);

    try {
      let response;
      let data;

      // 1. MANUAL & AUTO SCAN
      if (activeTab === 'manual' || activeTab === 'auto') {
        response = await fetch('http://127.0.0.1:8000/analyze', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ target_id: target, smiles: smiles, mode: activeTab }),
        });
        data = await response.json();
      } 
      
      // 2. UPLOAD MODE
      else if (activeTab === 'upload') {
         const formData = new FormData();
         formData.append('file', selectedFile);

         response = await fetch('http://127.0.0.1:8000/upload', {
            method: 'POST',
            body: formData,
         });
         data = await response.json();
      }

      if (data.error) {
        alert(data.error);
      } else {
        if (data.type === 'single') {
          setResult(data);
        } else if (data.type === 'batch') {
          setBatchResults(data.results);
        }
      }

    } catch (error) {
      console.error("Backend Error:", error);
      alert("⚠️ Backend Error! Make sure server is running.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="app-container">
      <div className="bg-grid"></div>
      <div className="floating-elements">
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
          <div className={`nav-link ${showAbout ? 'active-btn' : ''}`} onClick={() => setShowAbout(true)}><Info size={18} /><span>About</span></div>
          <div className="nav-link"><Settings size={18} /><span>Config</span></div>
        </div>
      </nav>

      <div className={`slider-container ${showAbout ? 'slide-active' : ''}`}>
        
        {/* === PAGE 1: HOME === */}
        <div className="page-section" style={{ position: 'relative' }}>
          <div className="main-layout">
            
            {/* LEFT PANEL */}
            <div className="glass-panel panel-left">
              <div className="panel-header"><Database size={20} color="#00f3ff" /><h3 className="panel-title">INPUT CONFIGURATION</h3></div>
              
              <div className="tab-group">
                <button className={`tab-btn ${activeTab === 'manual' ? 'active' : ''}`} onClick={() => setActiveTab('manual')}>MANUAL</button>
                <button className={`tab-btn ${activeTab === 'auto' ? 'active' : ''}`} onClick={() => setActiveTab('auto')}>AUTO DB</button>
                <button className={`tab-btn ${activeTab === 'upload' ? 'active' : ''}`} onClick={() => setActiveTab('upload')}>UPLOAD</button>
              </div>

              {activeTab !== 'upload' && (
                <div className="input-group">
                  <label className="input-label">TARGET PROTEIN (PDB ID)</label>
                  <div className="input-wrapper"><Search size={16} color="#888" className="input-icon"/><input className="cyber-input" placeholder="Ex: 6LU7" value={target} onChange={(e) => setTarget(e.target.value)} /></div>
                  
                  {/* --- DISEASE SUGGESTIONS (FIXED: 3 ITEMS) --- */}
                  <div className="suggestions-box">
                    <span>Try:</span>
                    <span className="suggestion-text" onClick={() => setTarget('6LU7')}>Covid-19</span>
                    <span className="suggestion-text" onClick={() => setTarget('3PP0')}>Cancer</span>
                    <span className="suggestion-text" onClick={() => setTarget('1G0H')}>Diabetes</span>
                  </div>
                </div>
              )}

              {activeTab === 'manual' && (
                <div className="input-group">
                  <label className="input-label">LIGAND STRUCTURE (SMILES)</label>
                  <textarea className="cyber-input textarea" rows="4" placeholder="Enter Chemical SMILES..." value={smiles} onChange={(e) => setSmiles(e.target.value)} />
                  
                  {/* --- DRUG SUGGESTIONS (FIXED: 3 ITEMS) --- */}
                  <div className="suggestions-box">
                    <span>Try:</span>
                    <span className="suggestion-text" onClick={() => setSmiles('CC(=O)Nc1ccc(O)cc1')}>Panadol</span>
                    <span className="suggestion-text" onClick={() => setSmiles('CC(=O)Oc1ccccc1C(=O)O')}>Aspirin</span>
                    <span className="suggestion-text" onClick={() => setSmiles('CC(C)Cc1ccc(C(C)C(=O)O)cc1')}>Ibuprofen</span>
                  </div>
                </div>
              )}

              {activeTab === 'upload' && (
                <div 
                    className="input-group" 
                    onClick={() => fileInputRef.current.click()} 
                    style={{
                        textAlign:'center', border:'1px dashed rgba(255,255,255,0.2)', padding:'30px', 
                        borderRadius:'10px', color:'#888', marginBottom:'20px', cursor:'pointer', transition: '0.3s'
                    }}
                    onMouseOver={(e) => e.currentTarget.style.borderColor = '#00f3ff'}
                    onMouseOut={(e) => e.currentTarget.style.borderColor = 'rgba(255,255,255,0.2)'}
                >
                   <input 
                       type="file" 
                       ref={fileInputRef} 
                       style={{display: 'none'}} 
                       accept=".csv, .xlsx" 
                       onChange={handleFileSelect}
                   />
                   <Upload size={30} style={{marginBottom:'10px', opacity:0.5, color: selectedFile ? '#00f3ff' : 'inherit'}}/>
                   <div style={{fontSize:'12px', fontWeight: 'bold', color: selectedFile ? '#fff' : '#888'}}>
                       {selectedFile ? selectedFile.name : "Click to Upload File"}
                   </div>
                   <div style={{fontSize:'10px', marginTop:'5px', color:'#555'}}>
                       {selectedFile ? `${(selectedFile.size / 1024).toFixed(1)} KB` : "Supports CSV & Excel"}
                   </div>
                </div>
              )}

              <button className="cyber-btn" onClick={handleScan} disabled={loading}>
                <div className="btn-content">
                  {loading ? <Activity className="animate-spin" size={20}/> : (activeTab === 'manual' ? <Zap size={20} /> : <Microscope size={20} />)}
                  {loading ? "PROCESSING..." : (activeTab === 'upload' ? "PROCESS BATCH" : "INITIATE ANALYSIS")}
                </div>
              </button>
            </div>

            {/* RIGHT PANEL */}
            <div className="glass-panel panel-right">
               <div className="header-badge-container">
                  <div className="header-badge">
                    <Activity size={14} color="#00f3ff" className={loading ? "animate-pulse" : ""}/> 
                    <span className="badge-text">SYSTEM STATUS: <span style={{ color: loading ? '#00f3ff' : (result || batchResults.length > 0 ? '#00f3ff' : '#888') }}>{loading ? "PROCESSING..." : (result || batchResults.length > 0 ? "ANALYSIS COMPLETE" : "STANDBY")}</span></span>
                  </div>
               </div>

               <div className="hologram-wrapper">
                 {(loading || (!result && batchResults.length === 0)) ? (
                   <div className="hologram-inner">
                      <div className="dna-spinner"><Dna size={120} color="#00f3ff" strokeWidth={1.5} /></div>
                      <div className="ring-1"></div><div className="ring-2"></div><div className="core-glow"></div>
                      {loading && <div className="loading-text">{activeTab === 'auto' ? "SCANNING DATABASE..." : (activeTab === 'upload' ? "READING FILE..." : "LOADING ANALYTICS...")}</div>}
                   </div>
                 ) : (
                   activeTab === 'manual' && result ? (
                     // SINGLE RESULT GRAPH
                     <div className="result-container-flex">
                        <div className="result-content-row">
                            <div className="structure-box" style={{ borderColor: `${result.color}40`, padding: '10px', width: '45%' }}>
                                <img src={`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(smiles)}/image?width=300&height=300&bgcolor=transparent`} alt="Structure" className="structure-img" style={{ height: '180px' }} />
                                <div className="box-label">CHEMICAL STRUCTURE</div>
                            </div>
                            <div className="structure-box" style={{ borderColor: `${result.color}40`, padding: '0px', width: '50%', height: '220px', background: 'rgba(0,0,0,0.3)' }}>
                                <ResponsiveContainer width="100%" height="100%">
                                  <RadarChart cx="50%" cy="50%" outerRadius="70%" data={result.graph}>
                                    <PolarGrid stroke="#333" />
                                    <PolarAngleAxis dataKey="subject" tick={{ fill: '#888', fontSize: 10 }} />
                                    <PolarRadiusAxis angle={30} domain={[0, 100]} tick={false} axisLine={false} />
                                    <Radar name="Stats" dataKey="A" stroke={result.color} strokeWidth={2} fill={result.color} fillOpacity={0.4} />
                                  </RadarChart>
                                </ResponsiveContainer>
                                <div className="box-label" style={{marginTop:'-10px', paddingBottom:'10px'}}>ADMET PROFILE</div>
                            </div>
                        </div>
                     </div>
                   ) : (
                     // BATCH RESULTS LIST (INTERACTIVE & STICKY HEADER)
                     <div className="scan-results-list" style={{
                        position: 'relative', zIndex: 20, width: '100%', marginTop: '20px',
                        height: '350px', overflowY: 'auto', overflowX: 'hidden', paddingRight: '5px',
                        background: 'rgba(0,0,0,0.2)', borderRadius: '15px'
                     }}>
                        {/* STICKY HEADER */}
                        <div className="list-header" style={{
                            position: 'sticky', top: 0, zIndex: 30, 
                            background: 'rgba(5, 5, 10, 0.95)', backdropFilter: 'blur(10px)',
                            borderBottom: '1px solid #00f3ff', padding: '15px 15px',
                            display: 'flex', justifyContent: 'space-between', alignItems: 'center',
                            boxShadow: '0 5px 20px rgba(0,0,0,0.5)', color: '#00f3ff',
                            fontSize: '12px', fontWeight: 'bold', letterSpacing: '1px'
                        }}>
                          <div style={{display:'flex', alignItems:'center', gap:'5px'}}><List size={14}/> DRUG NAME</div>
                          <div style={{display:'flex', alignItems:'center', gap:'5px'}}>SCORE <ArrowUpDown size={12}/></div>
                        </div>
                        
                        {/* LIST ITEMS */}
                        <div style={{padding: '10px'}}>
                          {batchResults.map((item, index) => (
                            <div 
                              className={`scan-item ${selectedId === item.name ? 'active-row' : ''}`}
                              key={index} 
                              onClick={() => handleDrugClick(item)}
                              style={{
                                cursor:'pointer', display: 'flex', justifyContent: 'space-between', alignItems: 'center',
                                padding: '15px', marginBottom: '8px', 
                                background: selectedId === item.name ? 'rgba(0, 243, 255, 0.15)' : 'rgba(255,255,255,0.03)',
                                border: selectedId === item.name ? '1px solid #00f3ff' : '1px solid transparent',
                                borderRadius: '10px', transition: 'all 0.2s ease'
                              }}
                            >
                                <div style={{display:'flex', alignItems:'center', gap:'10px'}}>
                                  <div style={{
                                    width:'8px', height:'8px', borderRadius:'50%', 
                                    background: item.status === 'ACTIVE' ? '#00f3ff' : '#ff0055',
                                    boxShadow: item.status === 'ACTIVE' ? '0 0 10px #00f3ff' : 'none'
                                  }}></div>
                                  <div>
                                    <div className="drug-name" style={{fontWeight:'700', color:'#fff', fontSize:'14px'}}>{item.name}</div>
                                    <div className="drug-smiles" style={{fontSize:'10px', color:'#666', marginTop:'2px'}}>Tap for details</div>
                                  </div>
                                </div>
                                <div style={{display:'flex', alignItems:'center', gap:'10px'}}>
                                  <div className="match-score" style={{color: item.status === 'ACTIVE' ? '#00f3ff' : '#ff0055', fontWeight:'800', fontSize:'16px'}}>{item.score}</div>
                                  <ChevronRight size={14} color="#444" />
                                </div>
                            </div>
                          ))}
                        </div>
                     </div>
                   )
                 )}
                 <div className="grid-overlay" style={{zIndex: 0}}></div>
               </div>
            </div>
          </div>

          {/* === RESULT CARD === */}
          {result && (
            <div className={`result-card ${result.status === 'ACTIVE' ? 'active' : 'inactive'}`} style={{ borderColor: result.color }}>
              <div className="score-section">
                <div className="score-box"><div className="input-label">BINDING AFFINITY</div><div className="score-val fade-in-text" style={{ textShadow: `0 0 20px ${result.color}` }}>{result.score}</div></div>
                <div className="vertical-div"></div>
                <div><div className="input-label">PREDICTION</div><div className="status-val fade-in-text" style={{ color: result.color }}>{result.status}</div></div>
              </div>
              <div className="metrics-container fade-in-text">
                 <div className="metric-box"><div className="input-label">TOXICITY</div><div className="metric-row" style={{ color: '#00f3ff' }}><ShieldCheck size={18} /> Low Risk</div></div>
                 <div className="metric-box"><div className="input-label">CONFIDENCE</div><div className="metric-row" style={{ color: '#bc13fe' }}><Dna size={18} /> {result.confidence}</div></div>
              </div>
              <button style={{background:'rgba(255,255,255,0.1)', border:'1px solid rgba(255,255,255,0.2)', color:'#fff', padding:'10px', borderRadius:'50%', cursor:'pointer', display:'flex', alignItems:'center', justifyContent:'center', transition:'0.3s'}}>
                 <Download size={18} />
              </button>
            </div>
          )}

        </div>

        {/* === PAGE 2: ABOUT === */}
        <div className="page-section">
          <div className="about-page-content">
             <div className="hologram-inner" style={{marginBottom: '30px', height: '150px'}}>
                <div className="dna-spinner"><Hexagon size={80} color="#00f3ff" fill="rgba(0, 243, 255, 0.3)" strokeWidth={2} /></div>
             </div>
             <div className="brand-text">BioGraph <span style={{ color: '#00f3ff' }}>AI</span></div> <div><br /></div> 
             <div className="about-hero-subtitle">Next-Generation Drug Discovery</div>
             <div className="features-grid">
                <div className="feature-card"><div className="f-icon"><Brain size={30} /></div><div className="f-title">Graph Neural Networks</div><div className="f-desc">Utilizing advanced GNNs to predict molecular interactions with high accuracy.</div></div>
                <div className="feature-card"><div className="f-icon"><Share2 size={30} /></div><div className="f-title">Binding Affinity</div><div className="f-desc">Calculates the strength of interactions between protein targets and ligands.</div></div>
                <div className="feature-card"><div className="f-icon"><Globe size={30} /></div><div className="f-title">ADMET Profiling</div><div className="f-desc">Instant analysis of Toxicity, Solubility, and Drug-likeness properties.</div></div>
             </div>
          </div>
        </div>

      </div>
    </div>
  );
}

export default App;