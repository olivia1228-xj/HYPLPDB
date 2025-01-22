import React, { useState, useEffect } from 'react';
import { useParams, Link } from 'react-router-dom';
import {ArrowLeft, AlertCircle, ExternalLink, Database, FileText, Dna, Loader2, Beaker, Scale, Activity, ChevronLeft, ChevronRight} from 'lucide-react';
import { getPeptideData } from '../lib/api';
import type { PeptideData } from '../lib/db';
import StructureSection from '../components/StructureSection';

const PeptideInfoPage: React.FC = () => {
  const { id } = useParams<{ id: string }>();
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<PeptideData | null>(null);
  const [activeTab, setActiveTab] = useState<'overview' | 'structure' | 'references'>('overview');
  const [currentPage, setCurrentPage] = useState(1);
  const itemsPerPage = 5;

  useEffect(() => {
    const fetchData = async () => {
      try {
        if (!id) throw new Error('No peptide ID provided');
        setLoading(true);
        const peptideData = await getPeptideData(id);
        setData(peptideData);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to load peptide data');
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, [id]);

  // 分页逻辑
  const paginateItems = <T extends any>(items: T[]): T[] => {
    const startIndex = (currentPage - 1) * itemsPerPage;
    return items.slice(startIndex, startIndex + itemsPerPage);
  };
  const handleExternalLink = (url: string | null) => {
    if (!url) return;

    // 处理 PubMed URL
    if (url.includes('pubmed.ncbi.nlm.nih.gov')) {
      // 提取 PubMed ID
      const pmidMatch = url.match(/\d+/);
      if (pmidMatch) {
        const pmid = pmidMatch[0];
        window.open(`https://pubmed.ncbi.nlm.nih.gov/${pmid}/`, '_blank', 'noopener,noreferrer');
        return;
      }
    }
    // 处理其他 URL
    window.open(url, '_blank', 'noopener,noreferrer');
  };

  if (loading) {
    return (
      <div className="min-h-screen bg-amber-50 flex items-center justify-center">
        <div className="text-center">
          <Loader2 className="w-12 h-12 animate-spin text-amber-500 mx-auto" />
          <p className="mt-4 text-gray-600">Loading peptide data...</p>
        </div>
      </div>
    );
  }

  if (error || !data) {
    return (
      <div className="min-h-screen bg-amber-50 flex items-center justify-center">
        <div className="bg-white rounded-lg shadow-lg p-8 max-w-md w-full text-center">
          <AlertCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
          <h2 className="text-2xl font-bold text-gray-900 mb-2">Error Loading Data</h2>
          <p className="text-gray-600 mb-6">{error}</p>
          <Link
            to="/search"
            className="inline-flex items-center text-amber-500 hover:text-amber-600"
          >
            <ArrowLeft className="w-4 h-4 mr-2" />
            Back to Search
          </Link>
        </div>
      </div>
    );
  }

  const { basicInfo, physicalProperties, structure, externalReferences } = data;

  // 计算总页数
  const totalPages = Math.ceil(externalReferences.literature.length / itemsPerPage);

  return (
    <div className="min-h-screen bg-amber-50">
      <div className="max-w-7xl mx-auto px-4 py-8">
        {/* Back Button */}
        <Link
          to="/search"
          className="inline-flex items-center text-amber-500 hover:text-amber-600 mb-6"
        >
          <ArrowLeft className="w-4 h-4 mr-2" />
          Back to Search
        </Link>

        {/* Header */}
        <div className="bg-white rounded-2xl shadow-lg p-6 mb-8">
          <div className="flex justify-between items-start">
            <div>
              <h1 className="text-3xl font-bold text-gray-900 mb-2">{basicInfo.sequence}</h1>
              <p className="text-lg text-gray-700 mb-4">ID: {basicInfo.id}</p>
            </div>

          </div>
          <div className="bg-gray-50 p-4 rounded-lg mb-4">
          <p className="text-xl font-mono">Species:{basicInfo.species}</p>
            <p className="mt-2 text-sm text-gray-700">Sequence Length: {basicInfo.length} amino acids</p>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
            <div className="bg-gray-50 p-4 rounded-lg">
              <p className="text-sm text-gray-500">Protein Source</p>
              <p className="text-lg font-semibold">{basicInfo.source}</p>
            </div>
            <div className="bg-gray-50 p-4 rounded-lg">
              <p className="text-sm text-gray-500">Extraction Method</p>
              <p className="text-lg font-semibold">{basicInfo.extractionMethod}</p>
            </div>
            <div className="bg-gray-50 p-4 rounded-lg">
              <p className="text-sm text-gray-500">IC50</p>
              <p className="text-lg font-semibold">{basicInfo.activity.ic50}</p>
            </div>
            <div className="bg-gray-50 p-4 rounded-lg">
              <p className="text-sm text-gray-500">Validation Model</p>
              <p className="text-lg font-semibold">{basicInfo.activity.validationModel}</p>
            </div>
          </div>
        </div>

        {/* Navigation Tabs */}
        <div className="flex flex-wrap gap-4 mb-8">
          <button
            onClick={() => setActiveTab('overview')}
            className={`flex items-center px-6 py-2 rounded-full transition-colors ${
              activeTab === 'overview'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-gray-700 hover:bg-amber-100'
            }`}
          >
            <Activity className="w-4 h-4 mr-2" />
            Overview
          </button>
          <button
            onClick={() => setActiveTab('structure')}
            className={`flex items-center px-6 py-2 rounded-full transition-colors ${
              activeTab === 'structure'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-gray-700 hover:bg-amber-100'
            }`}
          >
            <Beaker className="w-4 h-4 mr-2" />
            Structure
          </button>
          <button
            onClick={() => setActiveTab('references')}
            className={`flex items-center px-6 py-2 rounded-full transition-colors ${
              activeTab === 'references'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-gray-700 hover:bg-amber-100'
            }`}
          >
            <FileText className="w-4 h-4 mr-2" />
            References
          </button>
        </div>

        {/* Content Area */}
        <div className="bg-white rounded-2xl shadow-lg p-6">
          {activeTab === 'overview' && (
            <div className="space-y-8">
              {/* Activity Information */}
              <div>
                <h2 className="text-xl font-bold mb-4">Activity Information</h2>
                <div className="bg-gray-50 rounded-lg p-6">
                  <div className="space-y-6">
                    <div>
                      <h3 className="font-medium text-gray-700 mb-2">Functions</h3>
                      <div className="flex flex-wrap gap-2">
                        {basicInfo.activity.functions.map((func, index) => (
                          <span
                            key={index}
                            className="px-3 py-1 bg-blue-100 text-blue-800 rounded-full text-sm"
                          >
                            {func}
                          </span>
                        ))}
                      </div>
                    </div>
                    <div>
                      <h3 className="font-medium text-gray-700 mb-2">Description</h3>
                      <p className="text-gray-600 leading-relaxed">
                        {basicInfo.activity.functionDescription}
                      </p>
                    </div>
                  </div>
                </div>
              </div>

              {/* Physical Properties */}
              <div>
                <h2 className="text-xl font-bold mb-4">Physical Properties</h2>
                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                  <div className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-500">Molecular Weight</p>
                    <p className="text-lg font-semibold">{physicalProperties.molecularWeight.toFixed(2)} Da</p>
                  </div>
                  <div className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-500">Isoelectric Point</p>
                    <p className="text-lg font-semibold">{physicalProperties.isoelectricPoint.toFixed(2)}</p>
                  </div>
                  <div className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-500">Charge at pH 7</p>
                    <p className="text-lg font-semibold">{physicalProperties.chargeAtPH7.toFixed(2)}</p>
                  </div>
                  <div className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-500">Hydrophobicity</p>
                    <p className="text-lg font-semibold">{physicalProperties.hydrophobicity.toFixed(2)}</p>
                  </div>
                  <div className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-500">Instability Index</p>
                    <p className="text-lg font-semibold">{physicalProperties.instabilityIndex.toFixed(2)}</p>
                  </div>
                  <div className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-500">Aromaticity</p>
                    <p className="text-lg font-semibold">{physicalProperties.aromaticity.toFixed(2)}</p>
                  </div>
                </div>
              </div>
            </div>
          )}

          {activeTab === 'structure' && data && (
            <StructureSection structure={data.structure} />
          )}

          {activeTab === 'references' && (
            <div className="space-y-8">
              {/* Literature References */}
              <div>
                <h3 className="text-lg font-semibold mb-4">Literature References</h3>
                <div className="space-y-3">
                  {paginateItems(externalReferences.literature).map((ref, index) => (
                    ref.display && (
                      <button
                        key={index}
                        onClick={() => handleExternalLink(ref.link)}
                        className="w-full text-left bg-gray-50 p-4 rounded-lg hover:bg-gray-100 transition-colors"
                      >
                        <div className="flex items-center justify-between">
                          <div>
                            <span className="font-medium">{ref.display}</span>
                            {ref.id && <p className="text-sm text-gray-500 mt-1">{ref.id}</p>}
                          </div>
                          <ExternalLink className="w-4 h-4 text-gray-400" />
                        </div>
                      </button>
                    )
                  ))}
                </div>

                {/* Pagination */}
                {totalPages > 1 && (
                  <div className="flex justify-center items-center space-x-2 mt-6">
                    <button
                      onClick={() => setCurrentPage(prev => Math.max(prev - 1, 1))}
                      disabled={currentPage === 1}
                      className="p-2 rounded-lg border text-gray-500 disabled:opacity-50"
                    >
                      <ChevronLeft className="w-5 h-5" />
                    </button>
                    <span className="text-gray-600">
                      Page {currentPage} of {totalPages}
                    </span>
                    <button
                      onClick={() => setCurrentPage(prev => Math.min(prev + 1, totalPages))}
                      disabled={currentPage === totalPages}
                      className="p-2 rounded-lg border text-gray-500 disabled:opacity-50"
                    >
                      <ChevronRight className="w-5 h-5" />
                    </button>
                  </div>
                )}
              </div>

              {/* Database References */}
              <div>
                <h3 className="text-lg font-semibold mb-4">Database References</h3>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  {Object.entries(externalReferences.databases).map(([key, value]) => (
                    value.display && (
                      <button
                        key={key}
                        onClick={() => handleExternalLink(value.link)}
                        className="w-full text-left bg-gray-50 p-4 rounded-lg hover:bg-gray-100 transition-colors"
                      >
                        <div className="flex items-center justify-between">
                          <div>
                            <p className="text-sm text-gray-500 capitalize">{key.replace(/([A-Z])/g, ' $1').trim()}</p>
                            <p className="font-medium">{value.display}</p>
                            {value.id !== value.display && <p className="text-sm text-gray-500 mt-1">{value.id}</p>}
                          </div>
                          <ExternalLink className="w-4 h-4 text-gray-400" />
                        </div>
                      </button>
                    )
                  ))}
                </div>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default PeptideInfoPage;