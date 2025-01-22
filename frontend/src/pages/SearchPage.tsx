import React, { useState, useEffect } from 'react';
import { useNavigate, useSearchParams, Link } from 'react-router-dom';
import { Search, AlertCircle, X, ChevronLeft, ChevronRight, Filter, Upload, Loader2, SlidersHorizontal, Atom, Dna } from 'lucide-react';
import {
  searchPeptides,
  performBlastComparison,
  performStructureSearch,
  type SearchParams,
  type SearchResult,
  type BlastComparisonResult,
  type StructureSearchResult
} from '../lib/api';

const SearchResultCard = ({ result, type }: {
  result: SearchResult | BlastComparisonResult | StructureSearchResult;
  type: 'simple' | 'blast' | 'structure';
}) => {
  return (
    <Link
      to={`/peptide/${result.id}`}
      className="block border border-gray-100 rounded-lg p-4 hover:bg-gray-50 transition-colors"
    >
      <div className="flex justify-between items-start mb-2">
        <h3 className="font-semibold text-gray-900">{result.id}</h3>
        {type === 'simple' && 'class' in result && (
          <span className="text-sm bg-amber-100 text-amber-800 px-2 py-1 rounded-full">
            {result.species}
          </span>
        )}
        {type === 'blast' && 'score' in result && (
          <div className="flex items-center space-x-2">
            <span className="text-sm bg-blue-100 text-blue-800 px-2 py-1 rounded-full">
              Score: {result.score.toFixed(1)}
            </span>
            <span className="text-sm bg-green-100 text-green-800 px-2 py-1 rounded-full">
              Identity: {(result.identity * 100).toFixed(1)}%
            </span>
          </div>
        )}
        {type === 'structure' && 'similarity' in result && (
          <span className="text-sm bg-purple-100 text-purple-800 px-2 py-1 rounded-full">
            Similarity: {(result.similarity * 100).toFixed(1)}%
          </span>
        )}
      </div>

      <p className="text-sm font-mono bg-gray-50 px-2 py-1 rounded mb-2">
        {result.sequence}
      </p>

      {type === 'simple' && 'source' in result && (
        <div className="text-sm text-gray-600">
          <p><span className="font-medium">Species:</span> {result.species}</p>
          <p><span className="font-medium">MW:</span> {result.molecularWeight.toFixed(1)} Da</p>
          <p><span className="font-medium">Length:</span> {result.sequenceLength} aa</p>
        </div>
      )}

      {type === 'blast' && 'score' in result && (
        <div className="text-sm text-gray-600 space-y-1">
          {result.sequenceLength && (
            <p><span className="font-medium">Length:</span> {result.sequenceLength} aa</p>
          )}
          {result.function && (
            <p><span className="font-medium">Function:</span> {result.function}</p>
          )}
        </div>
      )}

      {type === 'structure' && 'similarity' in result && (
        <div className="text-sm text-gray-600 space-y-1">
          {result.molecular_weight && (
            <p><span className="font-medium">MW:</span> {result.molecular_weight.toFixed(1)} Da</p>
          )}
          {result.smiles && (
            <p><span className="font-medium">SMILES:</span> {result.smiles}</p>
          )}
        </div>
      )}
    </Link>
  );
};

const SearchPage = () => {
  const [searchParams] = useSearchParams();
  const [searchMode, setSearchMode] = useState<'simple' | 'structure' | 'blast'>('simple');
  const [searchQuery, setSearchQuery] = useState(searchParams.get('q') || '');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [showAdvancedFilters, setShowAdvancedFilters] = useState(false);
  const [searchResults, setSearchResults] = useState<SearchResult[]>([]);
  const [blastResults, setBlastResults] = useState<BlastComparisonResult[]>([]);
  const [structureResults, setStructureResults] = useState<StructureSearchResult[]>([]);
  const [currentPage, setCurrentPage] = useState(1);
  const [pagination, setPagination] = useState<{
    total: number;
    page: number;
    per_page: number;
    total_pages: number;
  } | null>(null);
  const [filters, setFilters] = useState<SearchParams>({
    page: 1,
    per_page: 10,
    exact_match: false,
    id: '',
    species: '',
    sequence: '',
    functionDescription: '',
    mw_min: undefined,
    mw_max: undefined,
    pi_min: undefined,
    pi_max: undefined,
    charge_min: undefined,
    charge_max: undefined,
    ic50: '',
    pmid: '',
    uniprot: '',
    biopep: ''
  });

  const navigate = useNavigate();
  const validateSequence = (sequence: string): boolean => {
    const validChars = new Set('ACDEFGHIKLMNPQRSTVWY');
    // 去除 FASTA 头和空白字符
    const cleanSequence = sequence
      .split('\n')
      .filter(line => !line.startsWith('>'))
      .join('')
      .replace(/\s/g, '')
      .toUpperCase();

    return cleanSequence.split('').every(char => validChars.has(char));
  };

  const handleModeChange = (mode: 'simple' | 'structure' | 'blast') => {
    setSearchMode(mode);
    setSearchQuery('');
    setSelectedFile(null);
    setError(null);
    setSearchResults([]);
    setBlastResults([]);
    setStructureResults([]);
  };

  const handleSearch = async (e?: React.FormEvent) => {
    if (e) {
      e.preventDefault();
    }

    if (!searchQuery.trim() && !selectedFile) {
      setError('Please enter the search content or upload the file');
      return;
    }

    setLoading(true);
    setError(null);

    try {
      switch (searchMode) {
        case 'simple':
          const response = await searchPeptides({
            ...filters,
            q: searchQuery,
            page: currentPage,
            per_page: 10,
            exact_match: false // 默认使用模糊匹配
          });
          setSearchResults(response.items);
          setPagination(response.pagination);
          break;
        case 'blast':
          try {
            setLoading(true);
            setError('Sequence alignment in progress, please wait...');

            if (!searchQuery?.trim() && !selectedFile) {
              throw new Error('Please enter a sequence or upload a FASTA file');
            }

            if (searchQuery?.trim()) {
              const sequence = searchQuery.trim();
              // 检查是否是 FASTA 格式
              const lines = sequence.split('\n').filter(line => line.trim());
              const isFasta = lines[0].startsWith('>');

              // 获取实际序列部分
              const sequencePart = isFasta
                ? lines.slice(1).join('')
                : sequence;

              if (!validateSequence(sequencePart)) {
                throw new Error('Invalid amino acid sequence, only standard amino acid letters are allowed (eg:ACDEFGHIKLMNPQRSTVWY)');
              }

              // 如果不是 FASTA 格式，添加标准头
              const formattedSequence = isFasta
                ? sequence
                : `>Query\n${sequence}`;

              // 更新搜索框中的序列为格式化后的序列
              setSearchQuery(formattedSequence);
            }

            const blastResponse = await performBlastComparison(searchQuery, selectedFile);

            if (blastResponse.results && Array.isArray(blastResponse.results)) {
              setBlastResults(blastResponse.results);
              setPagination(blastResponse.pagination);
              setError(null);
            } else {
              throw new Error('未找到匹配的序列');
            }
          } catch (err) {
            console.error('BLAST search error:', err);
            setError(err instanceof Error ? err.message : 'BLAST比对失败，请重试');
            setBlastResults([]);
            setPagination(null);
          } finally {
            setLoading(false);
          }
          break;
        case 'structure':
          try {
            setLoading(true);
            const inputType = searchQuery.trim().startsWith('InChI=') ? 'inchi' : 'smiles';

            // 添加提示信息
            setError('Searching, please wait patiently...');

            const structureResponse = await performStructureSearch(
              searchQuery,
              selectedFile,
              inputType
            );

            if (structureResponse.results && Array.isArray(structureResponse.results)) {
              setStructureResults(structureResponse.results);
              setPagination(structureResponse.pagination);
              setError(null); // 清除提示信息
            } else {
              throw new Error('No matching structure found');
            }
          } catch (err) {
            console.error('Structure search error:', err);
            setError(err instanceof Error ? err.message : 'Search failed, please try again');
            setStructureResults([]);  // 清空结果
            setPagination(null);      // 清空分页
          } finally {
            setLoading(false);
          }
          break;
      }
    } catch (err) {
      console.error('Search error:', err);
      setError(err instanceof Error ? err.message : 'Search failed');
    } finally {
      setLoading(false);
    }
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      // Validate file type based on search mode
      const validTypes = {
        structure: ['.mol2', '.sdf', '.pdb'],
        blast: ['.fasta', '.fa', '.txt']
      };

      const fileExt = file.name.toLowerCase().substring(file.name.lastIndexOf('.'));
      const isValidType = searchMode === 'structure'
        ? validTypes.structure.includes(fileExt)
        : validTypes.blast.includes(fileExt);

      if (!isValidType) {
        setError(`Invalid file type. Allowed types for ${searchMode} search: ${
          searchMode === 'structure' 
            ? validTypes.structure.join(', ')
            : validTypes.blast.join(', ')
        }`);
        return;
      }

      setSelectedFile(file);
      setError(null);
    }
  };

  const handleClear = () => {
    setSearchQuery('');
    setSelectedFile(null);
    setError(null);
    setSearchResults([]);
    setBlastResults([]);
    setStructureResults([]);
  };

  useEffect(() => {
    if (searchParams.get('q')) {
      handleSearch();
    }
  }, [currentPage]);

  return (
    <div className="min-h-screen bg-amber-50">
      <div className="max-w-7xl mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-8">Search Database</h1>

        {/* Search Mode Tabs */}
        <div className="flex flex-wrap gap-4 mb-8">
          <button
            onClick={() => handleModeChange('simple')}
            className={`flex items-center px-6 py-2 rounded-full transition-colors ${
              searchMode === 'simple'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-gray-700 hover:bg-amber-100'
            }`}
          >
            <Search className="w-4 h-4 mr-2" />
            Simple Search
          </button>
          <button
            onClick={() => handleModeChange('structure')}
            className={`flex items-center px-6 py-2 rounded-full transition-colors ${
              searchMode === 'structure'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-gray-700 hover:bg-amber-100'
            }`}
          >
            <Atom className="w-4 h-4 mr-2" />
            Structure Search
          </button>
          <button
            onClick={() => handleModeChange('blast')}
            className={`flex items-center px-6 py-2 rounded-full transition-colors ${
              searchMode === 'blast'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-gray-700 hover:bg-amber-100'
            }`}
          >
            <Dna className="w-4 h-4 mr-2" />
            BLAST Search
          </button>
        </div>

        {/* Search Forms */}
        <div className="bg-white rounded-2xl shadow-lg p-6 mb-8">
          {searchMode === 'simple' && (
            <form onSubmit={handleSearch} className="space-y-4">
              <div>
                <label htmlFor="simple-search" className="block text-base font-medium text-gray-700 mb-2">
                  Search by ID, species, sequence, or function
                </label>
                <div className="relative">
                  <input
                    type="text"
                    id="simple-search"
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    placeholder="Enter ID (e.g., HYPLPDB0001), species, sequence,Inhibits cholesterol synthesis..."
                    className="w-full px-6 py-3 rounded-full border-2 border-amber-200 focus:border-amber-500 focus:ring-0 pr-24"
                  />
                  {searchQuery && (
                    <button
                      type="button"
                      onClick={handleClear}
                      className="absolute right-20 top-1/2 -translate-y-1/2 text-gray-400 hover:text-gray-600"
                    >
                      <X className="w-5 h-5" />
                    </button>
                  )}
                  <button
                    type="submit"
                    disabled={loading}
                    className="absolute right-2 top-1/2 -translate-y-1/2 px-4 py-1.5 bg-amber-500 text-white rounded-full hover:bg-amber-600 transition-colors disabled:opacity-50"
                  >
                    {loading ? (
                      <Loader2 className="w-4 h-4 animate-spin" />
                    ) : (
                      <Search className="w-4 h-4" />
                    )}
                  </button>
                </div>
              </div>

              {/* Advanced Filters Toggle */}
              <button
                type="button"
                onClick={() => setShowAdvancedFilters(!showAdvancedFilters)}
                className="flex items-center text-amber-600 hover:text-amber-700"
              >
                <SlidersHorizontal className="w-4 h-4 mr-2" />
                {showAdvancedFilters ? 'Hide Advanced Filters' : 'Show Advanced Filters'}
              </button>

              {/* Updated Advanced Filters */}
              {showAdvancedFilters && (
                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-4 bg-amber-50 rounded-lg">
                  {/* Basic Information */}
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      ID
                    </label>
                    <input
                      type="text"
                      placeholder="e.g., HYPLPDB0001"
                      value={filters.id || ''}
                      onChange={(e) => setFilters(prev => ({
                        ...prev,
                        id: e.target.value
                      }))}
                      className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    />
                  </div>

                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Species
                    </label>
                    <input
                      type="text"
                      placeholder="e.g., Homo sapiens"
                      value={filters.species || ''}
                      onChange={(e) => setFilters(prev => ({
                        ...prev,
                        species: e.target.value
                      }))}
                      className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    />
                  </div>

                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Sequence
                    </label>
                    <input
                      type="text"
                      placeholder="e.g., AAVAAP"
                      value={filters.sequence || ''}
                      onChange={(e) => setFilters(prev => ({
                        ...prev,
                        sequence: e.target.value
                      }))}
                      className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    />
                  </div>

                  {/* Physical Properties */}
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Molecular Weight Range (Da)
                    </label>
                    <div className="flex gap-2">
                      <input
                        type="number"
                        step="0.1"
                        placeholder="Min"
                        value={filters.mw_min || ''}
                        onChange={(e) => setFilters(prev => ({
                          ...prev,
                          mw_min: e.target.value ? Number(e.target.value) : undefined
                        }))}
                        className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                      />
                      <input
                        type="number"
                        step="0.1"
                        placeholder="Max"
                        value={filters.mw_max || ''}
                        onChange={(e) => setFilters(prev => ({
                          ...prev,
                          mw_max: e.target.value ? Number(e.target.value) : undefined
                        }))}
                        className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                      />
                    </div>
                  </div>

                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Function Description
                    </label>
                    <input
                      type="text"
                      placeholder="e.g., Antimicrobial"
                      value={filters.functionDescription || ''}
                      onChange={(e) => setFilters(prev => ({
                        ...prev,
                        functionDescription: e.target.value
                      }))}
                      className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    />
                  </div>

                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      IC50
                    </label>
                    <input
                      type="text"
                      placeholder="e.g., 10 µM"
                      value={filters.ic50 || ''}
                      onChange={(e) => setFilters(prev => ({
                        ...prev,
                        ic50: e.target.value
                      }))}
                      className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    />
                  </div>

                  {/* External References */}
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      External Database IDs
                    </label>
                    <div className="space-y-2">
                      <input
                        type="text"
                        placeholder="PMID"
                        value={filters.pmid || ''}
                        onChange={(e) => setFilters(prev => ({
                          ...prev,
                          pmid: e.target.value
                        }))}
                        className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                      />
                      <input
                        type="text"
                        placeholder="UniProt ID"
                        value={filters.uniprot || ''}
                        onChange={(e) => setFilters(prev => ({
                          ...prev,
                          uniprot: e.target.value
                        }))}
                        className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                      />
                      <input
                        type="text"
                        placeholder="BioPep ID"
                        value={filters.biopep || ''}
                        onChange={(e) => setFilters(prev => ({
                          ...prev,
                          biopep: e.target.value
                        }))}
                        className="w-full px-3 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                      />
                    </div>
                  </div>
                </div>
              )}
            </form>
          )}

          {searchMode === 'structure' && (
            <form onSubmit={handleSearch} className="space-y-4">
              <div>
                <label htmlFor="structure-search" className="block text-base font-medium text-gray-700 mb-2">
                  Enter chemical structure format
                </label>
                <div className="space-y-4">
                  <div className="relative">  {/* 添加relative定位 */}
                    <textarea
                      id="structure-search"
                      value={searchQuery}
                      onChange={(e) => setSearchQuery(e.target.value)}
                      placeholder="Enter SMILES or InChI format..."
                      rows={4}
                      className="w-full px-4 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    />
                    {/* 添加清空按钮 */}
                    {searchQuery && (
                      <button
                        type="button"
                        onClick={handleClear}
                        className="absolute right-2 top-2 text-gray-400 hover:text-gray-600"
                      >
                        <X className="w-5 h-5" />
                      </button>
                    )}
                  </div>
                  <div className="flex items-center gap-4">
                    <div className="relative">
                      <input
                        type="file"
                        id="structure-file"
                        onChange={handleFileUpload}
                        accept=".mol2,.sdf,.pdb"
                        className="hidden"
                      />
                      <label
                        htmlFor="structure-file"
                        className="flex items-center px-6 py-2 bg-amber-100 text-amber-700 rounded-full hover:bg-amber-200 transition-colors cursor-pointer"
                      >
                        <Upload className="w-4 h-4 mr-2" />
                        Upload Structure File (MOL2, SDF, PDB)
                      </label>
                    </div>
                  </div>
                  {selectedFile && (
                    <div className="flex items-center justify-between bg-gray-50 p-3 rounded-lg">
                      <span className="text-sm text-gray-600">{selectedFile.name}</span>
                      <button
                        type="button"
                        onClick={() => setSelectedFile(null)}
                        className="text-gray-400 hover:text-gray-600"
                      >
                        <X className="w-4 h-4" />
                      </button>
                    </div>
                  )}
                </div>
              </div>
              <button
                type="submit"
                disabled={loading || (!searchQuery && !selectedFile)}
                className="w-full px-6 py-3 bg-amber-500 text-white rounded-full hover:bg-amber-600 transition-colors disabled:opacity-50 flex items-center justify-center"
              >
                {loading ? (
                  <>
                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                    Searching...
                  </>
                ) : (
                  <>
                    <Search className="w-4 h-4 mr-2" />
                    Search Similar Structures
                  </>
                )}
              </button>
            </form>
          )}

          {searchMode === 'blast' && (
            <form onSubmit={handleSearch} className="space-y-4">
              <div>
                <label htmlFor="blast-search" className="block text-base font-medium text-gray-700 mb-2">
                  Enter protein sequence or upload FASTA file
                </label>
                <div className="space-y-4">
                  <div className="relative">  {/* 添加relative定位 */}
                    <textarea
                      id="blast-search"
                      value={searchQuery}
                      onChange={(e) => setSearchQuery(e.target.value)}
                      placeholder="Enter protein sequence in FASTA format..."
                      rows={6}
                      className="w-full px-4 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0 font-mono"
                    />
                    {/* 添加清空按钮 */}
                    {searchQuery && (
                      <button
                        type="button"
                        onClick={handleClear}
                        className="absolute right-2 top-2 text-gray-400 hover:text-gray-600"
                      >
                        <X className="w-5 h-5" />
                      </button>
                    )}
                  </div>
                  <div className="relative">
                    <input
                      type="file"
                      id="blast-file"
                      onChange={handleFileUpload}
                      accept=".fasta,.fa,.txt"
                      className="hidden"
                    />
                    <label
                      htmlFor="blast-file"
                      className="flex items-center px-6 py-2 bg-amber-100 text-amber-700 rounded-full hover:bg-amber-200 transition-colors cursor-pointer inline-block"
                    >
                      <Upload className="w-4 h-4 mr-2" />
                      Upload FASTA File (.fasta, .fa, .txt)
                    </label>
                  </div>
                  {selectedFile && (
                    <div className="flex items-center justify-between bg-gray-50 p-3 rounded-lg">
                      <span className="text-sm text-gray-600">{selectedFile.name}</span>
                      <button
                        type="button"
                        onClick={() => setSelectedFile(null)}
                        className="text-gray-400 hover:text-gray-600"
                      >
                        <X className="w-4 h-4" />
                      </button>
                    </div>
                  )}
                </div>
              </div>
              <button
                type="submit"
                disabled={loading || (!searchQuery && !selectedFile)}
                className="w-full px-6 py-3 bg-amber-500 text-white rounded-full hover:bg-amber-600 transition-colors disabled:opacity-50 flex items-center justify-center"
              >
                {loading ? (
                  <>
                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                    Running BLAST...
                  </>
                ) : (
                  <>
                    <Search className="w-4 h-4 mr-2" />
                    Run BLAST Comparison
                  </>
                )}
              </button>
            </form>
          )}
        </div>

        {/* Search Tips */}
        <div className="bg-white rounded-2xl shadow-lg p-6 mb-8">
          <h2 className="text-lg font-semibold mb-4">Search Tips</h2>
          <div className="space-y-2 text-gray-600">
            {searchMode === 'simple' && (
                <>
                  <p>• Search by ID (e.g., "HYPLPDB0001" or "1")</p>
                  <p>• Search by sequence (e.g., "AAVAAP")</p>
                  <p>• Search by species (e.g., "Tenebrio")</p>
                  <p>• Search by function (e.g., "Inhibits cholesterol synthesis")</p>
                  <p>• Use advanced filters for more specific results</p>
                </>
            )}
            {searchMode === 'structure' && (
              <>
                <p>• Enter SMILES notation (e.g., "CC(=O)NC(C)C(=O)OH")</p>
                <p>• Enter InChI string(e.g., "InChI=1S/C6H9N3O2.C3H7NO2/c7-5(6(10)11)1-4-2;1-2(4)3(5)6/h2-3,5H,1,7H2,(H,8,9)(H,10,11)")</p>
                <p>• Upload structure files in MOL2, SDF, or PDB format</p>
                <p>• Results are ranked by structural similarity</p>
              </>
            )}
            {searchMode === 'blast' && (
              <>
                <p>• Enter protein sequence (e.g., "AAVLPVESSVVK")</p>
                <p>• Upload FASTA files (.fasta, .fa, .txt)</p>
                <p>• Use standard amino acid codes (A, F, D, L, etc.)</p>
                <p>• Results show sequence identity and alignment scores</p>
              </>
            )}
          </div>
        </div>

        {/* Error Message */}
        {error && (
          <div className="mb-8 p-4 bg-red-50 border border-red-200 rounded-lg text-red-700 flex items-center">
            <AlertCircle className="w-5 h-5 mr-2" />
            {error}
          </div>
        )}

        {/* Results Display */}
        {/* No Results Message */}
        {!loading && ((searchQuery || selectedFile)) &&
         !error && (
          (searchMode === 'simple' && searchResults.length === 0) ||
          (searchMode === 'blast' && blastResults.length === 0) ||
          (searchMode === 'structure' && structureResults.length === 0)
        ) && (
          <div className="bg-white rounded-2xl shadow-lg p-6 mb-8">
            <div className="text-center py-8">
              <h3 className="text-xl font-semibold text-gray-900 mb-4">No matches found</h3>
              <p className="text-gray-600 mb-6">Try adjusting your search with these suggestions:</p>

              {searchMode === 'simple' && (
                <ul className="text-left max-w-md mx-auto space-y-3">
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Check for typos in your search term
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Try using fewer or different keywords
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Use the advanced filters for more specific results
                  </li>
                </ul>
              )}

              {searchMode === 'blast' && (
                <ul className="text-left max-w-md mx-auto space-y-3">
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Verify that your sequence uses valid amino acid codes
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Try a shorter sequence fragment
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Check the FASTA format if you're using a file
                  </li>
                </ul>
              )}

              {searchMode === 'structure' && (
                <ul className="text-left max-w-md mx-auto space-y-3">
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Verify your SMILES or InChI notation
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Try searching with a simpler substructure
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Check that your structure file is in MOL2, SDF, or PDB format
                  </li>
                  <li className="flex items-center text-gray-600">
                    <span className="w-2 h-2 bg-amber-400 rounded-full mr-3"></span>
                    Ensure your chemical notation follows standard format
                  </li>
                </ul>
              )}

              <button
                onClick={handleClear}
                className="mt-8 px-6 py-2 bg-amber-100 text-amber-700 rounded-full hover:bg-amber-200 transition-colors inline-flex items-center"
              >
                <X className="w-4 h-4 mr-2" />
                Clear Search
              </button>
            </div>
          </div>
        )}
        {/* Simple Search Results */}
        {searchResults.length > 0 && (
          <div className="bg-white rounded-2xl shadow-lg p-6">
            <div className="flex justify-between items-center mb-6">
              <h2 className="text-xl font-bold">搜索结果</h2>
              <span className="text-sm text-gray-500">
                找到 {pagination?.total || searchResults.length} 个结果
              </span>
            </div>
            <div className="space-y-4">
              {searchResults.map((result) => (
                <SearchResultCard
                  key={`simple-${result.id}-${result.sequence}`} // 确保唯一的 key
                  result={result}
                  type="simple"
                />
              ))}
            </div>
          </div>
        )}

        {/* Structure Results */}
        {structureResults.length > 0 && (
          <div className="bg-white rounded-2xl shadow-lg p-6">
            <div className="flex justify-between items-center mb-6">
              <h2 className="text-xl font-bold">结构搜索结果</h2>
              <span className="text-sm text-gray-500">
                找到 {structureResults.length} 个匹配
              </span>
            </div>
            <div className="space-y-4">
              {structureResults.map((result) => (
                <SearchResultCard
                  key={`structure-${result.id}-${result.sequence}`} // 确保唯一的 key
                  result={result}
                  type="structure"
                />
              ))}
            </div>
          </div>
        )}

        {/* BLAST Results */}
        {blastResults.length > 0 && (
          <div className="bg-white rounded-2xl shadow-lg p-6">
            <div className="flex justify-between items-center mb-6">
              <h2 className="text-xl font-bold">BLAST 比对结果</h2>
              <span className="text-sm text-gray-500">
                找到 {blastResults.length} 个匹配序列
              </span>
            </div>
            <div className="space-y-4">
              {blastResults.map((result, index) => {
                // 生成唯一的 key
                const uniqueKey = result.id
                  ? `blast-${result.id}`
                  : `blast-${index}-${result.sequence?.substring(0, 10)}`;

                return (
                  <SearchResultCard
                    key={uniqueKey}
                    result={{
                      ...result,
                      // 确保必要的字段存在
                      id: result.id || `result-${index}`,
                      name: result.name || 'Unknown Sequence',
                      sequence: result.sequence || '',
                      score: result.score || 0,
                      identity: result.identity || 0,
                      function: result.function || ''
                    }}
                    type="blast"
                  />
                );
              })}
            </div>

            {/* Pagination */}
            {pagination && pagination.total_pages > 1 && (
              <div className="mt-6 flex justify-center">
                <div className="flex items-center space-x-2">
                  <button
                    onClick={() => setCurrentPage(currentPage - 1)}
                    disabled={currentPage === 1}
                    className="p-2 rounded-lg border text-gray-500 disabled:opacity-50"
                  >
                    <ChevronLeft className="w-5 h-5" />
                  </button>

                  {Array.from({ length: pagination.total_pages }, (_, i) => i + 1)
                    .filter(page => {
                      const distance = Math.abs(page - currentPage);
                      return distance === 0 || distance === 1 || page === 1 || page === pagination.total_pages;
                    })
                    .map((page, index, array) => (
                      <React.Fragment key={page}>
                        {index > 0 && array[index - 1] !== page - 1 && (
                          <span className="text-gray-400">...</span>
                        )}
                        <button
                          onClick={() => setCurrentPage(page)}
                          className={`w-8 h-8 rounded-lg ${
                            currentPage === page
                              ? 'bg-amber-500 text-white'
                              : 'text-gray-500 hover:bg-gray-100'
                          }`}
                        >
                          {page}
                        </button>
                      </React.Fragment>
                    ))}

                  <button
                    onClick={() => setCurrentPage(currentPage + 1)}
                    disabled={currentPage === pagination.total_pages}
                    className="p-2 rounded-lg border text-gray-500 disabled:opacity-50"
                  >
                    <ChevronRight className="w-5 h-5" />
                  </button>
                </div>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
};

export default SearchPage;