import React, { useState } from 'react';
import { Plus, Minus, RefreshCw, Move } from 'lucide-react';

/**
 * InteractiveImage
 * - Zoom (wheel + buttons)
 * - Drag to pan
 * - Minimal inline styles so main CSS stays same
 * - Keyed by parent (key={smiles}) to reset on SMILES change
 */
const controlBtnStyle = {
  background: 'transparent',
  border: 'none',
  color: '#fff',
  cursor: 'pointer',
  padding: '5px',
  borderRadius: '4px',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  transition: '0.2s'
};

export default function InteractiveImage({ src, alt, onLoadStart, onLoad, onError }) {
  const [scale, setScale] = useState(1);
  const [position, setPosition] = useState({ x: 0, y: 0 });
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState({ x: 0, y: 0 });

  const handleWheel = (e) => {
    e.preventDefault();
    e.stopPropagation();
    const delta = e.deltaY > 0 ? -0.1 : 0.1;
    const newScale = Math.min(Math.max(0.5, scale + delta), 5);
    setScale(newScale);
  };

  const handleMouseDown = (e) => {
    e.preventDefault();
    setIsDragging(true);
    setDragStart({ x: e.clientX - position.x, y: e.clientY - position.y });
  };

  const handleMouseMove = (e) => {
    if (isDragging) {
      e.preventDefault();
      setPosition({
        x: e.clientX - dragStart.x,
        y: e.clientY - dragStart.y
      });
    }
  };

  const handleMouseUp = () => {
    setIsDragging(false);
  };

  return (
    <div
      className="interactive-wrapper"
      style={{
        width: '100%',
        height: '100%',
        position: 'relative',
        overflow: 'hidden',
        cursor: isDragging ? 'grabbing' : 'grab',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        background: 'rgba(0,0,0,0.2)'
      }}
      onWheel={handleWheel}
      onMouseDown={handleMouseDown}
      onMouseMove={handleMouseMove}
      onMouseUp={handleMouseUp}
      onMouseLeave={handleMouseUp}
    >
      {/* Floating Controls */}
      <div style={{
        position: 'absolute',
        bottom: 20,
        right: 20,
        zIndex: 20,
        display: 'flex',
        flexDirection: 'column',
        gap: '8px',
        background: 'rgba(0,0,0,0.6)',
        padding: '8px',
        borderRadius: '12px',
        backdropFilter: 'blur(5px)',
        border: '1px solid rgba(255,255,255,0.1)'
      }}>
        <button onClick={(e) => { e.stopPropagation(); setScale(s => Math.min(5, s + 0.5)); }} style={controlBtnStyle} title="Zoom In"><Plus size={18} /></button>
        <button onClick={(e) => { e.stopPropagation(); setScale(s => Math.max(0.5, s - 0.5)); }} style={controlBtnStyle} title="Zoom Out"><Minus size={18} /></button>
        <button onClick={(e) => { e.stopPropagation(); setScale(1); setPosition({ x: 0, y: 0 }); }} style={controlBtnStyle} title="Reset View"><RefreshCw size={18} /></button>
      </div>

      {/* Indicator overlay */}
      <div style={{ position: 'absolute', top: 15, left: 15, color: 'rgba(255,255,255,0.3)', pointerEvents: 'none', display: 'flex', gap: 5, fontSize: '10px', fontWeight: 'bold' }}>
        <Move size={12} /> INTERACTIVE 2D
      </div>

      <img
        src={src}
        alt={alt}
        onLoadStart={onLoadStart}
        onLoad={onLoad}
        onError={onError}
        draggable={false}
        style={{
          transform: `translate(${position.x}px, ${position.y}px) scale(${scale})`,
          transition: isDragging ? 'none' : 'transform 0.1s ease-out',
          maxWidth: '85%',
          maxHeight: '85%',
          width: 'auto',
          height: 'auto',
          objectFit: 'contain',
          pointerEvents: 'none',
          userSelect: 'none',
          filter: 'drop-shadow(0 0 15px rgba(0,243,255,0.15))'
        }}
      />
    </div>
  );
}
