import React, { useEffect, useRef, useState, useCallback } from 'react';
import { Maximize2, Minimize2, RotateCw, ZoomIn, ZoomOut } from 'lucide-react';
import { debounce } from 'lodash';

const OptimizedMoleculeViewer = ({ pdbData }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [isSpinning, setIsSpinning] = useState(true); // Default to true for continuous rotation
  const [quality, setQuality] = useState('medium');
  const viewerRef = useRef(null);
  const viewerInstanceRef = useRef(null);
  const containerRef = useRef(null);

  // Debounced render function to prevent too frequent updates
  const debouncedRender = useCallback(
    debounce(() => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.render();
      }
    }, 16), // ~60fps
    []
  );

  const initViewer = useCallback(() => {
    if (!viewerRef.current || !window.$3Dmol) return;

    try {
      const viewer = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'white',
        antialias: true,
        disableFog: true, // Disable fog for better performance
        preserveDrawingBuffer: false // Better performance
      });

      viewerInstanceRef.current = viewer;
      updateViewer();
    } catch (error) {
      console.error('Error initializing 3D viewer:', error);
    }
  }, []);

  const updateViewer = useCallback(() => {
    if (!viewerInstanceRef.current || !pdbData) return;

    try {
      const viewer = viewerInstanceRef.current;
      viewer.clear();
      viewer.removeAllModels();
      viewer.addModel(pdbData, "pdb");

      // Adjust style based on quality setting
      const styles = {
        low: {
          stick: { radius: 0.1 },
          sphere: { radius: 0.3 }
        },
        medium: {
          stick: { radius: 0.15 },
          sphere: { radius: 0.4 }
        },
        high: {
          stick: { radius: 0.2 },
          sphere: { radius: 0.5 }
        }
      };

      viewer.setStyle({}, styles[quality]);
      viewer.zoomTo();
      viewer.spin(isSpinning);
      debouncedRender();
    } catch (error) {
      console.error('Error updating viewer:', error);
    }
  }, [pdbData, quality, isSpinning, debouncedRender]);

  // Load 3Dmol.js script
  useEffect(() => {
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = '/3Dmol-min.js';
      script.async = true;
      script.onload = initViewer;
      document.head.appendChild(script);
      return () => {
        document.head.removeChild(script);
      };
    } else {
      initViewer();
    }
  }, [initViewer]);

  // Update viewer when dependencies change
  useEffect(() => {
    updateViewer();
  }, [pdbData, quality, isSpinning, updateViewer]);

  // Optimized resize handler
  useEffect(() => {
    const handleResize = debounce(() => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.resize();
        debouncedRender();
      }
    }, 100);

    window.addEventListener('resize', handleResize);
    return () => {
      window.removeEventListener('resize', handleResize);
      handleResize.cancel();
    };
  }, [debouncedRender]);

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

    // Delayed resize for fullscreen transition
    setTimeout(() => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.resize();
        debouncedRender();
      }
    }, 100);
  };

  const toggleSpin = () => {
    setIsSpinning(!isSpinning);
  };

  const handleZoom = (factor) => {
    if (viewerInstanceRef.current) {
      viewerInstanceRef.current.zoom(factor);
      debouncedRender();
    }
  };

  return (
    <div className="flex flex-col h-full" ref={containerRef}>
      <div className="flex justify-between items-center p-2 bg-white border-b">
        <div className="flex items-center space-x-2">
          <button
            onClick={() => handleZoom(1.2)}
            className="p-2 text-gray-500 hover:bg-gray-100 rounded-lg"
            title="Zoom in"
          >
            <ZoomIn className="w-5 h-5" />
          </button>
          <button
            onClick={() => handleZoom(0.8)}
            className="p-2 text-gray-500 hover:bg-gray-100 rounded-lg"
            title="Zoom out"
          >
            <ZoomOut className="w-5 h-5" />
          </button>
        </div>

        <div className="flex items-center space-x-2">
          <select
            value={quality}
            onChange={(e) => setQuality(e.target.value)}
            className="text-sm border rounded p-1"
          >
            <option value="low">Low Quality</option>
            <option value="medium">Medium Quality</option>
            <option value="high">High Quality</option>
          </select>
          
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
            className="p-2 text-gray-500 hover:bg-gray-100 rounded-lg"
            title={isFullscreen ? 'Exit fullscreen' : 'Enter fullscreen'}
          >
            {isFullscreen ? (
              <Minimize2 className="w-5 h-5" />
            ) : (
              <Maximize2 className="w-5 h-5" />
            )}
          </button>
        </div>
      </div>

      <div className={`relative flex-grow ${isFullscreen ? 'fixed inset-0 z-50 bg-white' : ''}`}>
        <div
          ref={viewerRef}
          className={`w-full ${isFullscreen ? 'h-screen' : 'h-96'}`}
          style={{ position: 'relative' }}
        />
      </div>

      <div className="p-2 text-xs text-gray-500 bg-gray-50">
        <p>Click and drag to rotate • Scroll to zoom • Shift + Click to center</p>
      </div>
    </div>
  );
};

export default OptimizedMoleculeViewer;