import React from 'react';
import { ShieldCheck, Dna, Download } from 'lucide-react';

/**
 * Small presentational component for the floating result card.
 * Keeps markup centralized.
 */
export default function ResultCard({ result, cardRef, onDownload }) {
  if (!result) return null;

  return (
    <div className={`result-card ${result.status === 'ACTIVE' ? 'active' : 'inactive'}`} style={{ borderColor: result.color, zIndex: 60, pointerEvents: 'auto' }} ref={cardRef}>
      <div className="score-section">
        <div className="score-box">
          <div className="input-label">BINDING AFFINITY</div>
          <div className="score-val fade-in-text" style={{ textShadow: `0 0 20px ${result.color}` }}>{result.score}</div>
        </div>
        <div className="vertical-div"></div>
        <div>
          <div className="input-label">PREDICTION</div>
          <div className="status-val fade-in-text" style={{ color: result.color }}>{result.status}</div>
        </div>
      </div>

      <div className="metrics-container fade-in-text">
        <div className="metric-box">
          <div className="input-label">TOXICITY</div>
          <div className="metric-row" style={{ color: '#00f3ff' }}><ShieldCheck size={18} /> Low Risk</div>
        </div>
        <div className="metric-box">
          <div className="input-label">CONFIDENCE</div>
          <div className="metric-row" style={{ color: '#bc13fe' }}><Dna size={18} /> {result.confidence}</div>
        </div>
      </div>

      <button className="download-btn" onClick={onDownload} style={{ background: 'rgba(255,255,255,0.1)', border: '1px solid rgba(255,255,255,0.2)', color: '#fff', padding: '10px', borderRadius: '50%', cursor: 'pointer', display: 'flex', alignItems: 'center', justifyContent: 'center', transition: '0.3s' }}>
        <Download size={18} />
      </button>
    </div>
  );
}
