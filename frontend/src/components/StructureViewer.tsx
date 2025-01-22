import React, { useEffect, useRef, useState } from 'react';
import { Maximize2, Minimize2, RotateCw } from 'lucide-react';

const MoleculeViewer = ({ pdbData }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [isSpinning, setIsSpinning] = useState(true);
  const viewerRef = useRef(null);
  const viewerInstanceRef = useRef(null);
  const containerRef = useRef(null);

  useEffect(() => {
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = '/3Dmol-min.js';
      script.async = true;

      script.onload = () => {
        initViewer();
      };

      document.head.appendChild(script);

      return () => {
        document.head.removeChild(script);
      };
    } else {
      initViewer();
    }
  }, []);

  useEffect(() => {
    if (!pdbData || !viewerInstanceRef.current) return;
    updateViewer();
  }, [pdbData, isSpinning]);

  const initViewer = () => {
    if (!viewerRef.current || !window.$3Dmol) return;

    try {
      const viewer = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'white',
        antialias: true
      });

      viewerInstanceRef.current = viewer;
      updateViewer();
    } catch (error) {
      console.error('Error initializing 3D viewer:', error);
    }
  };

  const updateViewer = () => {
    if (!viewerInstanceRef.current || !pdbData) return;

    try {
      const viewer = viewerInstanceRef.current;
      viewer.clear();
      viewer.addModel(pdbData, "pdb");

      // Set ball and stick style with adjusted sizes
      viewer.setStyle({}, {
        stick: { radius: 0.15, colorscheme: 'rasmol' },
        sphere: { radius: 0.4, colorscheme: 'rasmol' }
      });

      viewer.zoomTo();
      viewer.spin(isSpinning);
      viewer.render();
    } catch (error) {
      console.error('Error updating viewer:', error);
    }
  };

  const toggleFullscreen = () => {
    if (!containerRef.current) return;

    if (!isFullscreen) {
      if (containerRef.current.requestFullscreen) {
        containerRef.current.requestFullscreen();
      }
    } else {
      if (document.exitFullscreen) {
        document.exitFullscreen();
      }
    }
    setIsFullscreen(!isFullscreen);

    // Resize viewer after fullscreen change
    setTimeout(() => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.resize();
        viewerInstanceRef.current.render();
      }
    }, 100);
  };

  const toggleSpin = () => {
    setIsSpinning(!isSpinning);
    if (viewerInstanceRef.current) {
      viewerInstanceRef.current.spin(!isSpinning);
    }
  };

  // Handle window resize
  useEffect(() => {
    const handleResize = () => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.resize();
        viewerInstanceRef.current.render();
      }
    };

    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, []);

  return (
    <div className="flex flex-col h-full" ref={containerRef}>
      {/* Controls bar with shadow for better visual separation */}
      <div className="flex justify-end space-x-2 p-2 bg-white border-b">
        <button
          onClick={toggleSpin}
          className={`p-2 rounded-lg transition-colors ${
            isSpinning ? 'text-amber-500 bg-amber-50' : 'text-gray-500 hover:bg-gray-100'
          }`}
          title={isSpinning ? 'Stop rotation' : 'Start rotation'}
        >
          <RotateCw className="w-5 h-5" />
        </button>
        <button
          onClick={toggleFullscreen}
          className="p-2 text-gray-500 hover:bg-gray-100 rounded-lg transition-colors"
          title={isFullscreen ? 'Exit fullscreen' : 'Enter fullscreen'}
        >
          {isFullscreen ? (
            <Minimize2 className="w-5 h-5" />
          ) : (
            <Maximize2 className="w-5 h-5" />
          )}
        </button>
      </div>

      {/* Viewer container */}
      <div
        className={`relative flex-grow ${
          isFullscreen ? 'fixed inset-0 z-50 bg-white' : ''
        }`}
      >
        <div
          ref={viewerRef}
          className={`w-full ${isFullscreen ? 'h-screen' : 'h-[400px]'}`}
          style={{ position: 'relative' }}
        />
      </div>

      {/* Help text */}
      <div className="p-2 text-xs text-gray-500 bg-gray-50 rounded-b-lg">
        <p>Click and drag to rotate • Scroll to zoom • Shift + Click to center</p>
      </div>
    </div>
  );
};

export default MoleculeViewer;