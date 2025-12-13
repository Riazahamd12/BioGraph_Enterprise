import React, { useEffect } from 'react';
import { X, CheckCircle, AlertCircle, Info } from 'lucide-react';

const Toast = ({ message, type = 'info', onClose, duration = 3000 }) => {
  useEffect(() => {
    if (duration) {
      const timer = setTimeout(() => {
        onClose();
      }, duration);
      return () => clearTimeout(timer);
    }
  }, [duration, onClose]);

  const styles = {
    info: { bg: 'rgba(0, 243, 255, 0.1)', border: '#00f3ff', icon: <Info size={18} color="#00f3ff" /> },
    success: { bg: 'rgba(0, 255, 100, 0.1)', border: '#00ff64', icon: <CheckCircle size={18} color="#00ff64" /> },
    error: { bg: 'rgba(255, 0, 85, 0.1)', border: '#ff0055', icon: <AlertCircle size={18} color="#ff0055" /> },
  };

  const style = styles[type] || styles.info;

  return (
    <div style={{
      position: 'fixed',
      bottom: '20px',
      right: '20px',
      background: '#050508',
      border: `1px solid ${style.border}`,
      borderRadius: '8px',
      padding: '12px 20px',
      display: 'flex',
      alignItems: 'center',
      gap: '12px',
      boxShadow: '0 5px 20px rgba(0,0,0,0.5)',
      color: '#fff',
      zIndex: 9999,
      animation: 'slideIn 0.3s ease-out',
      backdropFilter: 'blur(10px)'
    }}>
      {style.icon}
      <span style={{ fontSize: '14px' }}>{message}</span>
      <button onClick={onClose} style={{ background: 'transparent', border: 'none', color: '#666', cursor: 'pointer', padding: 0, marginLeft: '10px' }}>
        <X size={14} />
      </button>
      <style>{`
        @keyframes slideIn {
          from { transform: translateY(100%); opacity: 0; }
          to { transform: translateY(0); opacity: 1; }
        }
      `}</style>
    </div>
  );
};

export default Toast;
