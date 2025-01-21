import React, { useState } from 'react';
import { AlertCircle, Check, Loader2, Link as LinkIcon, Mail, Info } from 'lucide-react';
import { submitPeptideData } from '../lib/api';

interface FormData {
  sequence: string;
  link: string;
  submittedBy: string;
  name: string;
  source: string;
  extractionMethod: string;
  functionDescription: string;
  ic50: string;
  validation: string;
}

const SubmitPage = () => {
  const [formData, setFormData] = useState<FormData>({
    sequence: '',
    link: '',
    submittedBy: '',
    name: '',
    source: '',
    extractionMethod: '',
    functionDescription: '',
    ic50: '',
    validation: ''
  });

  const [errors, setErrors] = useState<Partial<FormData>>({});
  const [submitStatus, setSubmitStatus] = useState<'idle' | 'submitting' | 'success' | 'error'>('idle');
  const [errorMessage, setErrorMessage] = useState<string>('');
  const [submissionId, setSubmissionId] = useState<string | null>(null);

  const validateForm = () => {
    const newErrors: Partial<FormData> = {};
    let isValid = true;

    // 验证序列
    if (!formData.sequence.trim()) {
      newErrors.sequence = 'Peptide sequence is required';
      isValid = false;
    } else {
      // 验证序列格式（只允许有效的氨基酸代码）
      const validAminoAcids = /^[ACDEFGHIKLMNPQRSTVWY]+$/i;
      if (!validAminoAcids.test(formData.sequence.toUpperCase())) {
        newErrors.sequence = 'Invalid sequence. Use only standard amino acid codes (A-Y)';
        isValid = false;
      }
    }

    // 验证链接
    if (!formData.link.trim()) {
      newErrors.link = 'Reference link is required';
      isValid = false;
    } else {
      try {
        new URL(formData.link);
      } catch {
        newErrors.link = 'Please enter a valid URL';
        isValid = false;
      }
    }

    // 验证邮箱
    if (!formData.submittedBy.trim()) {
      newErrors.submittedBy = 'Email is required';
      isValid = false;
    } else {
      const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
      if (!emailRegex.test(formData.submittedBy)) {
        newErrors.submittedBy = 'Please enter a valid email address';
        isValid = false;
      }
    }

    setErrors(newErrors);
    return isValid;
  };

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name]: value }));
    // 清除该字段的错误
    if (errors[name as keyof FormData]) {
      setErrors(prev => ({ ...prev, [name]: undefined }));
    }
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!validateForm()) {
      return;
    }

    setSubmitStatus('submitting');
    setErrorMessage('');

    try {
      const response = await submitPeptideData(formData);
      setSubmissionId(response.id);
      setSubmitStatus('success');
      
      // 重置表单
      setFormData({
        sequence: '',
        link: '',
        submittedBy: '',
        name: '',
        source: '',
        extractionMethod: '',
        functionDescription: '',
        ic50: '',
        validation: ''
      });
    } catch (error) {
      setSubmitStatus('error');
      setErrorMessage(error instanceof Error ? error.message : 'Submission failed. Please try again.');
    }
  };

  return (
    <div className="min-h-screen bg-amber-50">
      <div className="max-w-4xl mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-8">Submit New Peptide Data</h1>

        <div className="bg-white rounded-2xl shadow-lg p-6 mb-8">
          <form onSubmit={handleSubmit} className="space-y-6">
            {/* Required Information */}
            <div>
              <h2 className="text-xl font-bold text-gray-900 mb-4">Required Information</h2>
              <div className="space-y-4">
                <div>
                  <label htmlFor="sequence" className="block text-sm font-medium text-gray-700 mb-1">
                    Peptide Sequence <span className="text-red-500">*</span>
                  </label>
                  <input
                    type="text"
                    id="sequence"
                    name="sequence"
                    value={formData.sequence}
                    onChange={handleInputChange}
                    className={`w-full px-4 py-2 rounded-lg border-2 ${
                      errors.sequence ? 'border-red-300' : 'border-amber-200'
                    } focus:border-amber-500 focus:ring-0 font-mono`}
                    placeholder="Enter amino acid sequence (e.g., AKVPLG)"
                    disabled={submitStatus === 'submitting'}
                  />
                  {errors.sequence && (
                    <p className="mt-1 text-sm text-red-600 flex items-center">
                      <AlertCircle className="w-4 h-4 mr-1" />
                      {errors.sequence}
                    </p>
                  )}
                  <p className="mt-1 text-sm text-gray-500 flex items-center">
                    <Info className="w-4 h-4 mr-1" />
                    Use standard one-letter amino acid codes (A-Y)
                  </p>
                </div>

                <div>
                  <label htmlFor="link" className="block text-sm font-medium text-gray-700 mb-1">
                    Reference Link <span className="text-red-500">*</span>
                  </label>
                  <div className="relative">
                    <LinkIcon className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400 w-5 h-5" />
                    <input
                      type="url"
                      id="link"
                      name="link"
                      value={formData.link}
                      onChange={handleInputChange}
                      className={`w-full pl-10 pr-4 py-2 rounded-lg border-2 ${
                        errors.link ? 'border-red-300' : 'border-amber-200'
                      } focus:border-amber-500 focus:ring-0`}
                      placeholder="https://doi.org/..."
                      disabled={submitStatus === 'submitting'}
                    />
                  </div>
                  {errors.link && (
                    <p className="mt-1 text-sm text-red-600 flex items-center">
                      <AlertCircle className="w-4 h-4 mr-1" />
                      {errors.link}
                    </p>
                  )}
                  <p className="mt-1 text-sm text-gray-500 flex items-center">
                    <Info className="w-4 h-4 mr-1" />
                    Provide DOI or publication URL
                  </p>
                </div>

                <div>
                  <label htmlFor="submittedBy" className="block text-sm font-medium text-gray-700 mb-1">
                    Submitter Email <span className="text-red-500">*</span>
                  </label>
                  <div className="relative">
                    <Mail className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400 w-5 h-5" />
                    <input
                      type="email"
                      id="submittedBy"
                      name="submittedBy"
                      value={formData.submittedBy}
                      onChange={handleInputChange}
                      className={`w-full pl-10 pr-4 py-2 rounded-lg border-2 ${
                        errors.submittedBy ? 'border-red-300' : 'border-amber-200'
                      } focus:border-amber-500 focus:ring-0`}
                      placeholder="your.email@example.com"
                      disabled={submitStatus === 'submitting'}
                    />
                  </div>
                  {errors.submittedBy && (
                    <p className="mt-1 text-sm text-red-600 flex items-center">
                      <AlertCircle className="w-4 h-4 mr-1" />
                      {errors.submittedBy}
                    </p>
                  )}
                  <p className="mt-1 text-sm text-gray-500 flex items-center">
                    <Info className="w-4 h-4 mr-1" />
                    We'll use this to contact you about your submission
                  </p>
                </div>
              </div>
            </div>

            {/* Optional Information */}
            <div>
              <h2 className="text-xl font-bold text-gray-900 mb-4">Additional Information (Optional)</h2>
              <div className="space-y-4">
                <div>
                  <label htmlFor="name" className="block text-sm font-medium text-gray-700 mb-1">
                    Peptide Name
                  </label>
                  <input
                    type="text"
                    id="name"
                    name="name"
                    value={formData.name}
                    onChange={handleInputChange}
                    className="w-full px-4 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    disabled={submitStatus === 'submitting'}
                  />
                </div>

                <div>
                  <label htmlFor="source" className="block text-sm font-medium text-gray-700 mb-1">
                    Source
                  </label>
                  <input
                    type="text"
                    id="source"
                    name="source"
                    value={formData.source}
                    onChange={handleInputChange}
                    className="w-full px-4 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    placeholder="e.g., Milk protein, Soybean"
                    disabled={submitStatus === 'submitting'}
                  />
                </div>

                <div>
                  <label htmlFor="functionDescription" className="block text-sm font-medium text-gray-700 mb-1">
                    Function Description
                  </label>
                  <textarea
                    id="functionDescription"
                    name="functionDescription"
                    value={formData.functionDescription}
                    onChange={handleInputChange}
                    rows={3}
                    className="w-full px-4 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    disabled={submitStatus === 'submitting'}
                  />
                </div>

                <div>
                  <label htmlFor="ic50" className="block text-sm font-medium text-gray-700 mb-1">
                    IC50
                  </label>
                  <input
                    type="text"
                    id="ic50"
                    name="ic50"
                    value={formData.ic50}
                    onChange={handleInputChange}
                    className="w-full px-4 py-2 rounded-lg border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
                    placeholder="e.g., 50 µM"
                    disabled={submitStatus === 'submitting'}
                  />
                </div>
              </div>
            </div>

            {/* Submit Button */}
            <div className="flex items-center justify-between pt-6 border-t">
              <div className="text-sm text-gray-500 flex items-center">
                <span className="text-red-500 mr-1">*</span> Required fields
              </div>
              <button
                type="submit"
                disabled={submitStatus === 'submitting'}
                className="flex items-center px-6 py-3 bg-amber-500 text-white rounded-full hover:bg-amber-600 transition-colors disabled:opacity-50"
              >
                {submitStatus === 'submitting' ? (
                  <>
                    <Loader2 className="w-5 h-5 mr-2 animate-spin" />
                    Submitting...
                  </>
                ) : (
                  'Submit Data'
                )}
              </button>
            </div>
          </form>
        </div>

        {/* Status Messages */}
        {submitStatus === 'success' && (
          <div className="fixed bottom-4 right-4 flex items-center bg-green-100 border border-green-400 text-green-700 px-4 py-3 rounded shadow-lg">
            <Check className="w-5 h-5 mr-2" />
            <div>
              <p className="font-medium">Submission successful!</p>
              <p className="text-sm">Submission ID: {submissionId}</p>
            </div>
          </div>
        )}
        
        {submitStatus === 'error' && (
          <div className="fixed bottom-4 right-4 flex items-center bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded shadow-lg">
            <AlertCircle className="w-5 h-5 mr-2" />
            <div>
              <p className="font-medium">Submission failed</p>
              <p className="text-sm">{errorMessage}</p>
            </div>
          </div>
        )}
        {/* Footer */}
      <footer className="bg-white py-6 mt-12">
        <div className="max-w-7xl mx-auto px-4">
          <div className="border-t border-gray-200 pt-6">
            <div className="text-center text-gray-600 text-sm">
              <p className="font-medium mb-2">If you need any help, please contact us</p>
              <p>wudongchimian@163.com</p>
            </div>
          </div>
        </div>
      </footer>
      </div>
    </div>
  );
};

export default SubmitPage;