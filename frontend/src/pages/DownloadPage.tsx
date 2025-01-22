import React, { useState } from 'react';
import { Download, AlertCircle, Loader2, FileDown, FileJson, Table, Filter } from 'lucide-react';
import { downloadPeptides } from '../lib/api';

const DownloadPage = () => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [format, setFormat] = useState<'csv' | 'json'>('csv');

  const handleDownload = async () => {
    setLoading(true);
    setError(null);

    try {
      const blob = await downloadPeptides(format);
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      const timestamp = new Date().toISOString().split('T')[0];
      a.href = url;
      a.download = `peptides_${timestamp}.${format}`;
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
      document.body.removeChild(a);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Download failed');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-amber-50">
      <div className="max-w-4xl mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-8">Download Data</h1>

        <div className="bg-white rounded-2xl shadow-lg p-6 mb-8">
          <div className="space-y-6">
            {/* Format Selection */}
            <div>
              <h2 className="text-xl font-bold text-gray-900 mb-4">Select Format</h2>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <button
                  onClick={() => setFormat('csv')}
                  className={`flex items-center justify-center p-4 rounded-lg border-2 transition-colors ${
                    format === 'csv'
                      ? 'border-amber-500 bg-amber-50 text-amber-700'
                      : 'border-gray-200 hover:border-amber-200'
                  }`}
                >
                  <Table className="w-6 h-6 mr-3" />
                  <div className="text-left">
                    <div className="font-semibold">CSV Format</div>
                    <div className="text-sm text-gray-500">Excel compatible spreadsheet</div>
                  </div>
                </button>
                <button
                  onClick={() => setFormat('json')}
                  className={`flex items-center justify-center p-4 rounded-lg border-2 transition-colors ${
                    format === 'json'
                      ? 'border-amber-500 bg-amber-50 text-amber-700'
                      : 'border-gray-200 hover:border-amber-200'
                  }`}
                >
                  <FileJson className="w-6 h-6 mr-3" />
                  <div className="text-left">
                    <div className="font-semibold">JSON Format</div>
                    <div className="text-sm text-gray-500">Programmatic data access</div>
                  </div>
                </button>
              </div>
            </div>

            {/* Download Button */}
            <button
              onClick={handleDownload}
              disabled={loading}
              className="w-full flex items-center justify-center px-6 py-3 bg-amber-500 text-white rounded-lg hover:bg-amber-600 transition-colors disabled:opacity-50"
            >
              {loading ? (
                <>
                  <Loader2 className="w-5 h-5 mr-2 animate-spin" />
                  Preparing Download...
                </>
              ) : (
                <>
                  <FileDown className="w-5 h-5 mr-2" />
                  Download Dataset
                </>
              )}
            </button>

            {/* Error Message */}
            {error && (
              <div className="p-4 bg-red-50 border border-red-200 rounded-lg text-red-700 flex items-center">
                <AlertCircle className="w-5 h-5 mr-2" />
                {error}
              </div>
            )}
          </div>

          {/* Format Details */}
          <div className="mt-8 border-t pt-6">
            <h3 className="text-lg font-semibold mb-4">Format Details</h3>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div className="space-y-2">
                <h4 className="font-medium text-gray-900">CSV Format Includes:</h4>
                <ul className="text-sm text-gray-600 space-y-1">
                  <li>• Peptide sequences and species</li>
                  <li>• Physical properties (MW, pI, etc.)</li>
                  <li>• Protein Source and Function Description</li>
                  <li>• Activity data and references</li>
                </ul>
              </div>
              <div className="space-y-2">
                <h4 className="font-medium text-gray-900">JSON Format Includes:</h4>
                <ul className="text-sm text-gray-600 space-y-1">
                  <li>• Complete peptide data structure</li>
                  <li>• Nested property relationships</li>
                  <li>• External database references</li>
                  <li>• Full experimental details</li>
                </ul>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default DownloadPage;