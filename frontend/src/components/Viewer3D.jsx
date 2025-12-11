import React, { useEffect, useRef } from 'react';
import * as $3Dmol from '3dmol/build/3Dmol.js';

/**
 * Viewer3D
 * - lightweight wrapper around $3Dmol
 * - props: target (pdb id), smiles, style3D ('cartoon'|'stick'|'surface'), bgColor
 * - keeps viewer inside local ref and calls render on prop changes
 *
 * Note: This is a pure visual wrapper. Main logic lives in App (for shared viewerRef).
 */
export default function Viewer3D({ target = '6LU7', smiles = '', style3D = 'cartoon', bgColor = 'transparent', onLoadingChange = () => { }, viewerRefProp = null }) {
  const containerRef = useRef(null);
  const localViewerRef = useRef(null);

  useEffect(() => {
    if (!containerRef.current) return;
    // clear
    containerRef.current.innerHTML = '';
    onLoadingChange(true);

    const viewer = $3Dmol.createViewer(containerRef.current, { backgroundColor: bgColor });
    localViewerRef.current = viewer;
    if (viewerRefProp) viewerRefProp.current = viewer;

    const safety = setTimeout(() => onLoadingChange(false), 6000);

    $3Dmol.download(`pdb:${target}`, viewer, { doAssembly: true, noSecondaryStructure: false }, function () {
      clearTimeout(safety);
      try {
        viewer.setStyle({}, {});
        viewer.removeAllSurfaces();

        if (style3D === 'stick') {
          viewer.setStyle({}, { stick: { radius: 0.15, colorscheme: 'Jmol' } });
        } else if (style3D === 'surface') {
          viewer.setStyle({}, { line: { hidden: true } });
          viewer.addSurface($3Dmol.VDW, { opacity: 0.85, color: 'spectrum' }, { hetflag: false });
        } else {
          viewer.setStyle({}, { cartoon: { color: 'spectrum', thickness: 1.2, opacity: 1.0 } });
        }

        // Add ligand if smiles provided
        if (smiles) {
          const mol = $3Dmol.read(smiles, { format: 'smi' });
          viewer.addModel(mol, 'mol');
          viewer.setStyle({ model: -1 }, { stick: { colorscheme: 'greenCarbon', radius: 0.5 } });
        }

        viewer.zoomTo();
        viewer.spin('y', 0.5);
        viewer.render();
      } catch (e) {
        console.error('3D render error (Viewer3D):', e);
      } finally {
        onLoadingChange(false);
      }
    });

    // cleanup
    return () => {
      if (viewerRefProp) viewerRefProp.current = null;
      try {
        if (localViewerRef.current) {
          localViewerRef.current.stopAnimation && localViewerRef.current.stopAnimation();
          localViewerRef.current = null;
        }
      // eslint-disable-next-line no-unused-vars
      } catch (e) { /* ignore cleanup errors */ }
    };
  // intentionally include deps we care about
  }, [target, smiles, style3D, bgColor, onLoadingChange, viewerRefProp]);

  return <div ref={containerRef} style={{ width: '100%', height: '100%' }} />;
}
