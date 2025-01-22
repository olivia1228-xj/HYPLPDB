import axios from 'axios';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:5000';

// 创建 axios 实例，增加超时时间
const api = axios.create({
  baseURL: API_URL,
  timeout: 60000, // 增加到 30 秒
  headers: {
    'Content-Type': 'application/json'
  }
});

// 添加请求重试配置
api.interceptors.response.use(undefined, async (err) => {
  const { config, message } = err;
  if (!config || !config.retry) {
    return Promise.reject(err);
  }

  config.retryCount = config.retryCount || 0;

  if (config.retryCount >= config.retry) {
    return Promise.reject(err);
  }

  config.retryCount += 1;

  // 延迟重试
  const delay = new Promise(resolve => setTimeout(resolve, 1000));
  await delay;

  return api(config);
});

// FormData 拦截器保持不变
api.interceptors.request.use((config) => {
  if (config.data instanceof FormData) {
    delete config.headers['Content-Type'];
  }
  return config;
});

// Types
interface PaginatedResponse<T> {
  items: T[];
  total: number;
  page: number;
  per_page: number;
  total_pages: number;
}

export interface SearchParams {
  q?: string;
  page?: number;
  per_page?: number;
  exact_match?: boolean;
  id?: string;
  species?: string;
  sequence?: string;
  functionDescription?: string;
  mw_min?: number;
  mw_max?: number;
  pi_min?: number;
  pi_max?: number;
  charge_min?: number;
  charge_max?: number;
  ic50?: string;
  pmid?: string;
  uniprot?: string;
  biopep?: string;
}

export interface SearchResult {
  id: string;
  species: string;
  sequence: string;
  source: string;
  sequenceLength: number;
  molecularWeight: number;
  isoelectricPoint: number;
  chargeAtPH7: number;
  ic50: string;
  functionDescription: string;
}


export interface BlastMatch {
  id: string;
  species: string;
  sequence: string;
  score: number;
  identity: number;
  source: string;  // 改为必需字段
  sequenceLength: number;  // 改为必需字段
  function?: string;
}


export interface BlastComparisonResult {
  query_header: string;
  query_sequence: string;
  matches: BlastMatch[];
}

export interface StructureSearchResult {
  id: string;
  species: string;
  sequence: string;
  smiles?: string;
  inchi?: string;
  molecular_weight?: number;
  function?: string;
  source?: string;
  similarity: number;
}

export interface StructureSearchResponse {
  results: StructureSearchResult[];
  pagination: {
    total: number;
    page: number;
    per_page: number;
    total_pages: number;
  };
}



// API Functions
export const getPeptideData = async (id: string) => {
  try {
    const response = await api.get(`/api/peptides/${id}`);
    return response.data;
  } catch (error) {
    throw error;
  }
};

export const searchPeptides = async (params: SearchParams): Promise<PaginatedResponse<SearchResult>> => {
  try {
    const response = await api.get('/api/peptides/search', { params });
    return response.data;
  } catch (error) {
    throw error;
  }
};

export const performStructureSearch = async (input: string | null, file: File | null, type: 'smiles' | 'inchi' = 'smiles'): Promise<any> => {
  try {
    const formData = new FormData();

    if (file) {
      formData.append('file', file);
    } else if (input?.trim()) {
      formData.append(type, input.trim());
    } else {
      throw new Error('请输入结构信息或上传文件');
    }

    // 添加重试配置
    const response = await api.post('/api/peptides/structure-search', formData, {
      retry: 3,
      timeout: 30000,
      headers: {
        'Content-Type': 'multipart/form-data',
      }
    });

    if (!response.data || !response.data.results) {
      console.error('Invalid response structure:', response.data);
      throw new Error('The server returns data in an incorrect format');
    }

    return {
      results: response.data.results,
      pagination: response.data.pagination || {
        total: response.data.results.length,
        page: 1,
        per_page: response.data.results.length,
        total_pages: 1
      }
    };
  } catch (error) {
    console.error('Structure search error:', error);
    if (axios.isAxiosError(error)) {
      if (error.code === 'ECONNABORTED') {
        throw new Error('The search request timed out, please try again later');
      }
      throw new Error(error.response?.data?.error || 'Structure search failed, please try again');
    }
    throw error;
  }
};

export const performBlastComparison = async (sequence: string | null, file: File | null): Promise<any> => {
  try {
    const formData = new FormData();

    if (file) {
      formData.append('file', file);
    } else if (sequence?.trim()) {
      // 确保序列格式正确
      let formattedSequence = sequence.trim();
      // 如果输入的不是FASTA格式，添加FASTA头
      if (!formattedSequence.startsWith('>')) {
        formattedSequence = `>Query\n${formattedSequence}`;
      }
      formData.append('sequence', formattedSequence);
    } else {
      throw new Error('Please enter a sequence or upload a FASTA file');
    }

    const response = await api.post('/api/peptides/blast-compare', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      timeout: 30000, // 增加超时时间到30秒
      retry: 3
    });

    if (!response.data || (!response.data.results && !Array.isArray(response.data))) {
      console.error('Invalid response format:', response.data);
      throw new Error('服务器返回数据格式错误');
    }

    // 处理返回的数据
    const results = Array.isArray(response.data) ? response.data : response.data.results;

    return {
      results: results.map((result: any, index: number) => ({
        id: result.id || `result-${index}`,
        name: result.name || 'Unknown',
        sequence: result.sequence || '',
        score: result.score,
        identity: result.identity,
        function: result.function
      })),
      pagination: response.data.pagination || {
        total: results.length,
        page: 1,
        per_page: results.length,
        total_pages: 1
      }
    };
  } catch (error) {
    console.error('BLAST comparison error:', error);
    if (axios.isAxiosError(error)) {
      if (error.code === 'ECONNABORTED') {
        throw new Error('请求超时，请稍后重试');
      }
      if (error.response?.status === 500) {
        throw new Error('服务器处理错误，请检查输入格式是否正确');
      }
      throw new Error(error.response?.data?.error || '比对失败，请重试');
    }
    throw error;
  }
};
// Analytics API
export const getAnalyticsData = async () => {
  try {
    const response = await api.get('/api/sequence-analysis');
    return response.data;
  } catch (error) {
    throw error;
  }
};

// Download API
export const downloadPeptides = async (format: 'csv' | 'json' = 'csv'): Promise<Blob> => {
  try {
    const response = await api.get('/api/peptides/download', {
      params: { format },
      responseType: 'blob'
    });
    return response.data;
  } catch (error) {
    throw error;
  }
};

// Submit API
export interface SubmissionData {
  sequence: string;
  link: string;
  submittedBy: string;
  name?: string;
  source?: string;
  extractionMethod?: string;
  functionDescription?: string;
  ic50?: string;
  validation?: string;
}

export const submitPeptideData = async (data: SubmissionData): Promise<{ id: string }> => {
  try {
    const response = await api.post('/api/peptides/submit', data);
    return response.data;
  } catch (error) {
    throw error;
  }
};